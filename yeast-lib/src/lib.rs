#![deny(unused_crate_dependencies)]
use itertools::Itertools;
use lazy_static::lazy_static;
use petgraph::{
    dot::{Config, Dot},
    graphmap::DiGraphMap,
    visit::NodeIndexable,
    Direction,
};
use plotters::{
    drawing::IntoDrawingArea,
    prelude::{
        ChartBuilder, LabelAreaPosition, RGBAColor, Rectangle, SVGBackend, ShapeStyle, WHITE,
    },
};
use rayon::prelude::*;
use regex::Regex;
use std::fmt::{Debug, Display, Formatter, Result as FmtResult};
use std::{io::Write, path::Path};
use time::{format_description::well_known::Iso8601, Date, OffsetDateTime};

pub mod averager_subjective_values;
pub mod calc;
pub mod common;
pub mod density;
pub mod kinetics;
pub mod raw_data_processing;
pub mod strain_summary;
pub mod uncertain;

use crate::averager_subjective_values::Averager;
use crate::calc::{
    precision_probability, t_doubling_probability, PRECISION_SCALE, PRECISION_SHAPE,
    T_DOUBLING_MEAN_PRIOR, T_DOUBLING_VARIANCE_PRIOR,
};
use crate::common::Density;
use crate::kinetics::{
    ConcentrationFitLogistics, FitParams, FitPoint, N_POINTS_BEFORE_MEAN_PRECISION,
    N_POINTS_T_DOUBLING,
};
use crate::raw_data_processing::{
    liquid::{Clustering, KineticsLiquid, Liquid},
    organoleptic::{
        LiquidAppearance, Organoleptic, OrganolepticScoreSets, PassedOrganoleptic, YeastAppearance,
    },
    plate::{Plate, ThermalOutcome},
    slant::Slant,
    stock::Stock,
    strain::Strain,
};
use crate::strain_summary::StrainSummary;
use crate::uncertain::Uncertain;

#[derive(Debug)]
pub struct DataSet {
    pub liquids: Vec<Liquid>,
    pub organoleptics: Vec<Organoleptic>,
    pub plates: Vec<Plate>,
    pub slants: Vec<Slant>,
    pub stocks: Vec<Stock>,
    pub strains: Vec<Strain>,
}

impl DataSet {
    pub fn init(enable_log: bool) -> Self {
        Self {
            liquids: Liquid::from_dir(enable_log),
            organoleptics: Organoleptic::from_dir(enable_log),
            plates: Plate::from_dir(enable_log),
            slants: Slant::from_dir(enable_log),
            stocks: Stock::from_dir(enable_log),
            strains: Strain::from_file(enable_log),
        }
    }

    pub fn passed_organoleptic_by_id(&self, id: &str) -> Option<&PassedOrganoleptic> {
        for organoleptic in self.organoleptics.iter() {
            if let Organoleptic::Passed(passed_organoleptic) = organoleptic {
                if passed_organoleptic.id == id {
                    return Some(passed_organoleptic);
                }
            }
        }
        None
    }

    pub fn kinetics_liquid_by_id(&self, id: &str) -> Option<&KineticsLiquid> {
        for liquid in self.liquids.iter() {
            match liquid {
                Liquid::QA(kinetics_liquid) => {
                    if kinetics_liquid.id == id {
                        return Some(kinetics_liquid);
                    }
                }
                Liquid::Uniformity(kinetics_liquid) => {
                    if kinetics_liquid.id == id {
                        return Some(kinetics_liquid);
                    }
                }
                _ => {}
            }
        }
        None
    }

    pub fn uniformity_kinetics_liquid_by_id(&self, id: &str) -> Option<&KineticsLiquid> {
        for liquid in self.liquids.iter() {
            if let Liquid::Uniformity(kinetics_liquid) = liquid {
                if kinetics_liquid.id == id {
                    return Some(kinetics_liquid);
                }
            }
        }
        None
    }

    pub fn qa_kinetics_liquid_by_id(&self, id: &str) -> Option<&KineticsLiquid> {
        for liquid in self.liquids.iter() {
            if let Liquid::QA(kinetics_liquid) = liquid {
                if kinetics_liquid.id == id {
                    return Some(kinetics_liquid);
                }
            }
        }
        None
    }

    pub fn slant_by_id(&self, id: &str) -> Option<&Slant> {
        for slant in self.slants.iter() {
            match slant {
                Slant::Origin(origin_slant) => {
                    if origin_slant.id == id {
                        return Some(slant);
                    }
                }
                Slant::Lab(regular_slant) => {
                    if regular_slant.id == id {
                        return Some(slant);
                    }
                }
                Slant::Prod(regular_slant) => {
                    if regular_slant.id == id {
                        return Some(slant);
                    }
                }
            }
        }
        None
    }

    /// Find strain code corresponding to provided sample id
    pub fn strain_code_for_id(&self, id: &str) -> Option<String> {
        for strain in self.strains.iter() {
            let di_graph_map_whole = self.di_graph_map_whole(Some(&strain.code));
            for node in di_graph_map_whole.nodes() {
                if node.id == id {
                    return Some(strain.code.to_owned());
                }
            }
        }
        None
    }

    pub fn orphans(&self) {
        let di_graph_map_whole = self.di_graph_map_whole(None);

        macro_rules! orphan {
            ($data_part: ident, $text: literal) => {
                for element in self.$data_part.iter() {
                    let mut is_orphan = true;
                    let id = element.id();
                    for node in di_graph_map_whole.nodes() {
                        if node.id == id {
                            is_orphan = false;
                            break;
                        }
                    }
                    if is_orphan {
                        println!("ORPHAN: {} {id} does not belong to any strain.", $text)
                    }
                }
            };
        }

        orphan!(liquids, "Liquid");
        orphan!(organoleptics, "Organoleptic");
        orphan!(plates, "Plate");
        orphan!(slants, "Slant");
    }

    pub fn stock_density_by_id(&self, id: &str) -> Option<&Density> {
        for stock in self.stocks.iter() {
            if stock.id == id {
                return stock.density.as_ref();
            }
        }
        None
    }

    pub fn progeny_plaques_slant_details(&self, parent_id: &str) -> Vec<Plaque> {
        let mut out: Vec<Plaque> = Vec::new();
        for liquid in self.liquids.iter() {
            match liquid {
                Liquid::Package(package_liquid) => {
                    if package_liquid.parent == parent_id {
                        out.push(Plaque {
                            id: package_liquid.id.as_ref(),
                            kind: PlaqueKind::LiquidPackage {
                                count: package_liquid.batch_count,
                            },
                        })
                    }
                }
                Liquid::Pitch(pitch_liquid) => {
                    if pitch_liquid.parent == parent_id {
                        if pitch_liquid.fermentation_temperature.is_lagered() {
                            out.push(Plaque {
                                id: pitch_liquid.id.as_ref(),
                                kind: PlaqueKind::LiquidPitchLagered {
                                    style: pitch_liquid.style.as_ref(),
                                },
                            })
                        } else {
                            out.push(Plaque {
                                id: pitch_liquid.id.as_ref(),
                                kind: PlaqueKind::LiquidPitch {
                                    style: pitch_liquid.style.as_ref(),
                                },
                            })
                        }
                    }
                }
                Liquid::Propagation(regular_liquid) => {
                    if regular_liquid.parent == parent_id {
                        out.push(Plaque {
                            id: regular_liquid.id.as_ref(),
                            kind: PlaqueKind::LiquidRegular,
                        })
                    }
                }
                Liquid::QA(kinetics_liquid) => {
                    if kinetics_liquid.parent == parent_id {
                        if kinetics_liquid.temperature.is_lagered() {
                            out.push(Plaque {
                                id: kinetics_liquid.id.as_ref(),
                                kind: PlaqueKind::LiquidRegularLagered,
                            })
                        } else {
                            out.push(Plaque {
                                id: kinetics_liquid.id.as_ref(),
                                kind: PlaqueKind::LiquidRegular,
                            })
                        }
                    }
                }
                Liquid::Revive(regular_liquid) => {
                    if regular_liquid.parent == parent_id {
                        out.push(Plaque {
                            id: regular_liquid.id.as_ref(),
                            kind: PlaqueKind::LiquidRegular,
                        })
                    }
                }
                Liquid::Uniformity(kinetics_liquid) => {
                    if kinetics_liquid.parent == parent_id {
                        if kinetics_liquid.temperature.is_lagered() {
                            out.push(Plaque {
                                id: kinetics_liquid.id.as_ref(),
                                kind: PlaqueKind::LiquidRegularLagered,
                            })
                        } else {
                            out.push(Plaque {
                                id: kinetics_liquid.id.as_ref(),
                                kind: PlaqueKind::LiquidRegular,
                            })
                        }
                    }
                }
            }
        }
        for organoleptic in self.organoleptics.iter() {
            match organoleptic {
                Organoleptic::Failed(failed_organoleptic) => {
                    if failed_organoleptic
                        .participants
                        .set
                        .contains(&parent_id.to_owned())
                    {
                        out.push(Plaque {
                            id: failed_organoleptic.id.as_ref(),
                            kind: PlaqueKind::OrganolepticFailed,
                        })
                    }
                }
                Organoleptic::Passed(passed_organoleptic) => {
                    if passed_organoleptic
                        .participants
                        .set
                        .contains(&parent_id.to_owned())
                    {
                        out.push(Plaque {
                            id: passed_organoleptic.id.as_ref(),
                            kind: PlaqueKind::OrganolepticPassed,
                        })
                    }
                }
            }
        }
        for plate in self.plates.iter() {
            match plate {
                Plate::Regular(regular_plate) => {
                    if regular_plate.parent == parent_id {
                        out.push(Plaque {
                            id: regular_plate.id.as_ref(),
                            kind: PlaqueKind::PlateRegular,
                        })
                    }
                }
                Plate::Thermal(thermal_plate) => {
                    if thermal_plate.parent == parent_id {
                        match thermal_plate.thermal_outcome {
                            ThermalOutcome::Ale => out.push(Plaque {
                                id: thermal_plate.id.as_ref(),
                                kind: PlaqueKind::PlateThermalAle,
                            }),
                            ThermalOutcome::Lager => out.push(Plaque {
                                id: thermal_plate.id.as_ref(),
                                kind: PlaqueKind::PlateThermalLager,
                            }),
                        }
                    }
                }
            }
        }
        out
    }

    pub fn progeny_plaques(&self, parent_id: &str) -> Vec<Plaque> {
        let mut out = self.progeny_plaques_slant_details(parent_id);
        for slant in self.slants.iter() {
            match slant {
                Slant::Lab(regular_slant) => {
                    if regular_slant.parent == parent_id {
                        out.push(Plaque {
                            id: regular_slant.id.as_ref(),
                            kind: PlaqueKind::SlantRegular,
                        })
                    }
                }
                Slant::Prod(regular_slant) => {
                    if regular_slant.parent == parent_id {
                        out.push(Plaque {
                            id: regular_slant.id.as_ref(),
                            kind: PlaqueKind::SlantProd,
                        })
                    }
                }
                _ => {}
            }
        }
        out
    }

    pub fn organoleptics_for_slant<'a>(&'a self, slant_node: Plaque<'a>) -> OrganolepticsDone {
        let mut active_ends = vec![slant_node];
        let mut edges: Vec<(Plaque, Plaque)> = Vec::new();
        while let Some(current_active_end) = active_ends.pop() {
            for progeny in self
                .progeny_plaques_slant_details(current_active_end.id)
                .into_iter()
            {
                edges.push((current_active_end, progeny));
                active_ends.push(progeny);
            }
        }
        let di_graph_map_slant_details = DiGraphMap::<Plaque, Empty>::from_edges(edges);
        for external in di_graph_map_slant_details
            .clone()
            .into_graph::<usize>()
            .externals(Direction::Outgoing)
        {
            match di_graph_map_slant_details.from_index(external.index()).kind {
                PlaqueKind::OrganolepticFailed => {
                    return OrganolepticsDone::Failed;
                }
                PlaqueKind::OrganolepticPassed => {
                    return OrganolepticsDone::Passed;
                }
                _ => {}
            }
        }
        OrganolepticsDone::None
    }

    pub fn di_graph_map_whole<'a>(
        &'a self,
        filter_strain: Option<&'a str>,
    ) -> DiGraphMap<Plaque, Empty> {
        let mut active_ends: Vec<Plaque> = Vec::new();
        if let Some(strain_code) = filter_strain {
            for slant in self.slants.iter() {
                if let Slant::Origin(origin_slant) = slant {
                    if origin_slant.strain_code == strain_code {
                        active_ends.push(Plaque {
                            id: origin_slant.id.as_ref(),
                            kind: PlaqueKind::SlantOrigin { strain_code },
                        })
                    }
                }
            }
        } else {
            for strain in self.strains.iter() {
                for slant in self.slants.iter() {
                    if let Slant::Origin(origin_slant) = slant {
                        if origin_slant.strain_code == strain.code {
                            active_ends.push(Plaque {
                                id: origin_slant.id.as_ref(),
                                kind: PlaqueKind::SlantOrigin {
                                    strain_code: origin_slant.strain_code.as_ref(),
                                },
                            })
                        }
                    }
                }
            }
        }
        let mut edges: Vec<(Plaque, Plaque)> = Vec::new();
        while let Some(current_active_end) = active_ends.pop() {
            for progeny in self.progeny_plaques(current_active_end.id).into_iter() {
                edges.push((current_active_end, progeny));
                active_ends.push(progeny);
            }
        }
        DiGraphMap::<Plaque, Empty>::from_edges(edges)
    }

    pub fn plot_genealogy(&self, filter_strain: Option<&str>) {
        let graph = self.di_graph_map_whole(filter_strain);
        let graph_content = Dot::with_attr_getters(
            &graph,
            &[Config::EdgeNoLabel, Config::GraphContentOnly],
            &|_, _| String::new(),
            &|_, plaque| {
                format!(
                    "shape = \"{}\" fillcolor = \"{}\"",
                    plaque.1.shape(),
                    plaque.1.color()
                )
            },
        );
        let directory = format!("{OUTPUT}/{GENEALOGY_DIR}");
        if !Path::new(&directory).exists() {
            if let Err(e) = std::fs::create_dir_all(&directory) {
                println!("Error making genealogy data directory: {e}")
            }
        }
        let file_path = match filter_strain {
            Some(strain_code) => format!("{directory}/{strain_code}.dot"),
            None => format!("{directory}/{ALL_PLOTS}.dot"),
        };
        let mut file = std::fs::File::create(file_path).unwrap();
        write!(file, "{}", full_graph(graph_content.to_string()))
            .expect("Error while writing into {file_path}");
    }

    pub fn slant_status<'a>(&'a self, strain_code: &'a str) -> StrainStatus {
        let di_graph_map_whole = self.di_graph_map_whole(Some(strain_code));
        let mut slant_nodes_set: Vec<Plaque> = Vec::new();
        for node in di_graph_map_whole.nodes() {
            match node.kind {
                PlaqueKind::SlantOrigin { strain_code: _ } => slant_nodes_set.push(node),
                PlaqueKind::SlantRegular => slant_nodes_set.push(node),
                _ => {}
            }
        }
        slant_nodes_set.sort();

        let most_recent_slant_node = slant_nodes_set
            .pop()
            .expect("There must be at least one slant for each strain.");
        let id_most_recent = most_recent_slant_node.id;

        let is_backed_up_most_recent =
            is_slant_backed_up(&di_graph_map_whole, most_recent_slant_node);
        let organoleptics_done_most_recent = self.organoleptics_for_slant(most_recent_slant_node);

        let age_most_recent = {
            if let Some(ref captured_date) = REG_SAMPLE_NAME.captures(most_recent_slant_node.id) {
                if let Ok(date) = Date::parse(&captured_date["date"], &Iso8601::DEFAULT) {
                    let current_date = OffsetDateTime::now_local().unwrap().date();
                    let age_whole_weeks = (current_date - date).whole_weeks();
                    if age_whole_weeks < AGE_WEEKS_ATTN {
                        AgeChecked::Recent
                    } else if age_whole_weeks < AGE_WEEKS_URGENT {
                        AgeChecked::RenewSoon
                    } else {
                        AgeChecked::RenewImmediately
                    }
                } else {
                    AgeChecked::RenewImmediately
                }
            } else {
                AgeChecked::RenewImmediately
            }
        };

        let mut previous: Vec<SlantOld> = Vec::new();
        while let Some(slant_node) = slant_nodes_set.pop() {
            let slant_old = SlantOld {
                id: slant_node.id,
                is_backed_up: is_slant_backed_up(&di_graph_map_whole, slant_node),
                organoleptics_done: self.organoleptics_for_slant(slant_node),
            };
            previous.push(slant_old);
        }
        StrainStatus {
            strain_code,
            id_most_recent,
            is_backed_up_most_recent,
            organoleptics_done_most_recent,
            age_most_recent,
            previous,
        }
    }

    pub fn strain_exists(&self, strain_code: &str) -> Option<&Strain> {
        self.strains
            .iter()
            .find(|&strain| strain.code == strain_code)
    }

    pub fn update_caches_by_strain_code(&self, strain_code: &str) {
        let di_graph_map_whole = self.di_graph_map_whole(Some(strain_code));
        let mut kinetics_liquid_set: Vec<&KineticsLiquid> = Vec::new();
        for node in di_graph_map_whole.nodes() {
            // At this point only regular liquids here (e.g. not lagered ones)!
            if let PlaqueKind::LiquidRegular = node.kind {
                if let Some(kinetics_liquid) = self.kinetics_liquid_by_id(node.id) {
                    kinetics_liquid_set.push(kinetics_liquid)
                }
            }
        }
        let write_log = true;
        collect_concentration_fits(&kinetics_liquid_set, strain_code, write_log);
    }

    pub fn report_uniformity(&self, strain_code: &str) {
        let di_graph_map_whole = self.di_graph_map_whole(Some(strain_code));

        let write_log = true;

        let concentration_fit_groups_uniformity =
            self.concentration_fit_groups_uniformity(&di_graph_map_whole, strain_code, write_log);

        if let Some(mixed_uniformity_points) = mixed_points_if_strain_is_uniform(
            &concentration_fit_groups_uniformity,
            strain_code,
            write_log,
        ) {
            plot_heatmap(
                &mixed_uniformity_points,
                &format!("{OUTPUT}/{HEATMAPS}/{UNIFORMITY}"),
                &format!("{strain_code}_{UNIFORMITY}_heatmap.svg"),
            )
        }
    }

    pub fn report_variability(&self, strain_code: &str) {
        let di_graph_map_whole = self.di_graph_map_whole(Some(strain_code));

        let write_log = true;

        let concentration_fit_groups_uniformity =
            self.concentration_fit_groups_uniformity(&di_graph_map_whole, strain_code, write_log);

        let concentration_fit_groups_qa =
            self.concentration_fit_groups_qa(&di_graph_map_whole, strain_code, write_log);

        if let Some(mixed_uniformity_points) = mixed_points_if_strain_is_uniform(
            &concentration_fit_groups_uniformity,
            strain_code,
            write_log,
        ) {
            mixed_points_if_strain_is_stable_over_time(
                mixed_uniformity_points,
                &concentration_fit_groups_qa,
                strain_code,
                write_log,
            );
        } else {
            println!("Strain {strain_code} failed kinetics uniformity.")
        }
    }

    pub fn concentration_fit_groups_uniformity(
        &self,
        di_graph_map_whole: &DiGraphMap<Plaque, Empty>,
        strain_code: &str,
        write_log: bool,
    ) -> Vec<Vec<ConcentrationFitLogistics>> {
        let mut uniformity_grouped_kinetics_liquid_sets: Vec<Vec<&KineticsLiquid>> = Vec::new();
        for node in di_graph_map_whole.nodes() {
            // At this point only regular liquids here (e.g. not lagered ones)!
            if let PlaqueKind::LiquidRegular = node.kind {
                let parent_nodes: Vec<Plaque> = di_graph_map_whole
                    .neighbors_directed(node, Direction::Incoming)
                    .collect();
                if parent_nodes.len() == 1 && parent_nodes[0].kind == PlaqueKind::PlateRegular {
                    let mut current_group: Vec<&KineticsLiquid> = Vec::new();
                    let mut current_node = node;
                    let mut progeny_nodes: Vec<Plaque> = di_graph_map_whole
                        .neighbors_directed(current_node, Direction::Outgoing)
                        .collect();

                    while let Some(kinetics_liquid) =
                        self.uniformity_kinetics_liquid_by_id(current_node.id)
                    {
                        current_group.push(kinetics_liquid);
                        if progeny_nodes.len() == 1
                            && progeny_nodes[0].kind == PlaqueKind::LiquidRegular
                        {
                            current_node = progeny_nodes[0];
                            progeny_nodes = di_graph_map_whole
                                .neighbors_directed(current_node, Direction::Outgoing)
                                .collect();
                        } else {
                            break;
                        }
                    }

                    if !current_group.is_empty() {
                        uniformity_grouped_kinetics_liquid_sets.push(current_group)
                    }
                }
            }
        }
        uniformity_grouped_kinetics_liquid_sets
            .par_iter()
            .filter_map(|uniformity_kinetics_liquid_set| {
                let fit_group = collect_concentration_fits(
                    uniformity_kinetics_liquid_set,
                    strain_code,
                    write_log,
                );
                if !fit_group.is_empty() {
                    Some(fit_group)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn concentration_fit_groups_qa(
        &self,
        di_graph_map_whole: &DiGraphMap<Plaque, Empty>,
        strain_code: &str,
        write_log: bool,
    ) -> Vec<Vec<ConcentrationFitLogistics>> {
        let mut qa_grouped_kinetics_liquid_sets: Vec<Vec<&KineticsLiquid>> = Vec::new();
        for node in di_graph_map_whole.nodes() {
            // At this point only regular liquids here (e.g. not lagered ones)!
            if let PlaqueKind::LiquidRegular = node.kind {
                let parent_nodes: Vec<Plaque> = di_graph_map_whole
                    .neighbors_directed(node, Direction::Incoming)
                    .collect();
                if parent_nodes.len() == 1 && parent_nodes[0].kind == PlaqueKind::PlateRegular {
                    let mut current_group: Vec<&KineticsLiquid> = Vec::new();
                    let mut current_node = node;
                    let mut progeny_nodes: Vec<Plaque> = di_graph_map_whole
                        .neighbors_directed(current_node, Direction::Outgoing)
                        .collect();

                    while let Some(kinetics_liquid) = self.qa_kinetics_liquid_by_id(current_node.id)
                    {
                        current_group.push(kinetics_liquid);
                        if progeny_nodes.len() == 1
                            && progeny_nodes[0].kind == PlaqueKind::LiquidRegular
                        {
                            current_node = progeny_nodes[0];
                            progeny_nodes = di_graph_map_whole
                                .neighbors_directed(current_node, Direction::Outgoing)
                                .collect();
                        } else {
                            break;
                        }
                    }

                    if !current_group.is_empty() {
                        current_group.sort_by(|a, b| a.id.cmp(&b.id));
                        qa_grouped_kinetics_liquid_sets.push(current_group)
                    }
                }
            }
        }

        qa_grouped_kinetics_liquid_sets.sort_by(|a, b| {
            a.first()
                .expect("added only non-empty groups")
                .id
                .cmp(&b.first().expect("added only non-empty groups").id)
        });

        qa_grouped_kinetics_liquid_sets
            .par_iter()
            .filter_map(|qa_kinetics_liquid_set| {
                let fit_group =
                    collect_concentration_fits(qa_kinetics_liquid_set, strain_code, write_log);
                if !fit_group.is_empty() {
                    Some(fit_group)
                } else {
                    None
                }
            })
            .collect()
    }

    pub fn kinetics_heatmap_by_strain(&self, strain_code: &str) {
        let di_graph_map_whole = self.di_graph_map_whole(Some(strain_code));
        let mut kinetics_liquid_set: Vec<&KineticsLiquid> = Vec::new();
        for node in di_graph_map_whole.nodes() {
            // At this point only regular liquids here (e.g. not lagered ones)!
            if let PlaqueKind::LiquidRegular = node.kind {
                if let Some(kinetics_liquid) = self.kinetics_liquid_by_id(node.id) {
                    kinetics_liquid_set.push(kinetics_liquid)
                }
            }
        }

        let write_log = true;
        let concentration_fit_logistics_set =
            collect_concentration_fits(&kinetics_liquid_set, strain_code, write_log);

        let mixed_points = mix_points(MixerInput::RealData(
            concentration_fit_logistics_set
                .iter()
                .collect::<Vec<&ConcentrationFitLogistics>>()
                .as_slice(),
        ));

        plot_heatmap(
            &mixed_points,
            &format!("{OUTPUT}/{HEATMAPS}/{DATA_FIT}"),
            &format!("{strain_code}_{DATA_FIT}_heatmap.svg"),
        )
    }

    pub fn kinetics_plots_by_strain(&self, strain_code: &str) {
        let di_graph_map_whole = self.di_graph_map_whole(Some(strain_code));
        for node in di_graph_map_whole.nodes() {
            if let PlaqueKind::LiquidRegular = node.kind {
                if let Some(kinetics_liquid) = self.kinetics_liquid_by_id(node.id) {
                    if let Err(e) = kinetics_liquid.concentration_plot(strain_code) {
                        println!("{e}")
                    }
                    if let Err(e) = kinetics_liquid.density_plot(strain_code) {
                        println!("{e}")
                    }
                }
            }
            if let PlaqueKind::LiquidRegularLagered = node.kind {
                if let Some(kinetics_liquid) = self.kinetics_liquid_by_id(node.id) {
                    if let Err(e) = kinetics_liquid.concentration_plot(strain_code) {
                        println!("{e}")
                    }
                    if let Err(e) = kinetics_liquid.density_plot(strain_code) {
                        println!("{e}")
                    }
                }
            }
        }
    }

    pub fn kinetics_plots_by_id(&self, id: &str) {
        if let Some(strain_code) = self.strain_code_for_id(id) {
            if let Some(kinetics_liquid) = self.kinetics_liquid_by_id(id) {
                if let Err(e) = kinetics_liquid.concentration_plot(&strain_code) {
                    println!("{e}")
                }
                if let Err(e) = kinetics_liquid.density_plot(&strain_code) {
                    println!("{e}")
                }
            } else {
                println!("Id {id} does not correspond to known kinetics liquid record");
            }
        } else {
            println!("No node with {id}")
        }
    }

    pub fn strain_summary_data(&self, strain: &Strain) {
        let di_graph_map_whole = self.di_graph_map_whole(Some(&strain.code));
        let mut appearance_liquid_set: Vec<LiquidAppearance> = Vec::new();
        let mut scores: Vec<OrganolepticScoreSets> = Vec::new();
        let mut appearance_yeast_set: Vec<YeastAppearance> = Vec::new();
        let mut clustering_set: Vec<Clustering> = Vec::new();
        let mut real_attenuation_set: Vec<Uncertain> = Vec::new();
        let mut thermal_test = None;

        for node in di_graph_map_whole.nodes() {
            match node.kind {
                PlaqueKind::OrganolepticPassed => {
                    let passed_organoleptic = self.passed_organoleptic_by_id(node.id).unwrap();
                    if let Some(ref appearance) = passed_organoleptic.appearance {
                        appearance_liquid_set.push(appearance.liquid.clone());
                        appearance_yeast_set.push(appearance.yeast.clone());
                    }
                    if let Some(ref organoleptic_score_sets) =
                        passed_organoleptic.organoleptic_score_sets
                    {
                        scores.push(organoleptic_score_sets.clone());
                    }
                    let mut stock = None;
                    for neighbor in di_graph_map_whole.neighbors_directed(node, Direction::Incoming)
                    {
                        let kinetics_liquid = self.kinetics_liquid_by_id(neighbor.id).expect(
                            "150 mL vial, kinetics liquid, a standard predecessor for organoleptics",
                        );
                        stock = kinetics_liquid.stock.as_ref();
                        for neighbor2 in
                            di_graph_map_whole.neighbors_directed(neighbor, Direction::Incoming)
                        {
                            let kinetics_liquid2 = self.kinetics_liquid_by_id(neighbor2.id).expect(
                                "40 mL vial, kinetics liquid, a standard predecessor for 150 mL vial",
                            );
                            if stock == kinetics_liquid2.stock.as_ref() {
                                for neighbor3 in di_graph_map_whole
                                    .neighbors_directed(neighbor2, Direction::Incoming)
                                {
                                    let kinetics_liquid3 =
                                        self.kinetics_liquid_by_id(neighbor3.id).expect("10 mL vial, kinetics liquid, a standard predecessor for 40 mL vial");
                                    if stock != kinetics_liquid3.stock.as_ref() {
                                        stock = None
                                    }
                                }
                            } else {
                                stock = None
                            }
                        }
                    }
                    if let Some(stock_id) = stock {
                        if let Some(density) = self.stock_density_by_id(stock_id) {
                            let og = density.value / 1000.0;
                            let extract_og =
                                -668.962 + 1262.45 * og - 776.43 * og.powi(2) + 182.94 * og.powi(3);
                            if let Some(density_data) = &passed_organoleptic.density_data {
                                for density_entry in density_data.density_entries.iter() {
                                    let fg = (density_entry.density_value) / 1000.0;
                                    let extract_fg = -668.962 + 1262.45 * fg - 776.43 * fg.powi(2)
                                        + 182.94 * fg.powi(3);
                                    let q = 0.22 + 0.001 * extract_og;
                                    let real_extract = (q * extract_og + extract_fg) / (1.0 + q);
                                    let real_attenuation =
                                        100.0 * (extract_og - real_extract) / extract_og;
                                    real_attenuation_set.push(Uncertain {
                                        value: real_attenuation,
                                        error: 0.0,
                                    })
                                }
                            }
                        }
                    }
                }
                PlaqueKind::LiquidRegular | PlaqueKind::LiquidRegularLagered => {
                    if let Some(kinetics_liquid) = self.kinetics_liquid_by_id(node.id) {
                        for measurement in kinetics_liquid.measurement_set.set.iter() {
                            if let Some(ref clustering) = measurement.clustering {
                                if measurement.dilution.is_none()
                                    || clustering != &Clustering::Low
                                    || clustering_set.is_empty()
                                {
                                    clustering_set.push(clustering.clone())
                                }
                            }
                        }
                    }
                }
                PlaqueKind::PlateThermalAle => match thermal_test {
                    Some(ref previous_thermal_test) => {
                        if previous_thermal_test != &ThermalOutcome::Ale {
                            thermal_test = None
                        }
                    }
                    None => thermal_test = Some(ThermalOutcome::Ale),
                },
                PlaqueKind::PlateThermalLager => match thermal_test {
                    Some(ref previous_thermal_test) => {
                        if previous_thermal_test != &ThermalOutcome::Lager {
                            thermal_test = None
                        }
                    }
                    None => thermal_test = Some(ThermalOutcome::Lager),
                },
                _ => {}
            }
        }

        let write_log = false;

        let concentration_fit_groups_uniformity =
            self.concentration_fit_groups_uniformity(&di_graph_map_whole, &strain.code, write_log);

        let concentration_fit_groups_qa =
            self.concentration_fit_groups_qa(&di_graph_map_whole, &strain.code, write_log);

        let average_peak_doubling_time_minutes = {
            if let Some(mixed_uniformity_points) = mixed_points_if_strain_is_uniform(
                &concentration_fit_groups_uniformity,
                &strain.code,
                write_log,
            ) {
                match mixed_points_if_strain_is_stable_over_time(
                    mixed_uniformity_points,
                    &concentration_fit_groups_qa,
                    &strain.code,
                    write_log,
                ) {
                    Some(mixed_points) => mixed_points
                        .iter()
                        .max_by(|point1, point2| point1.probability.total_cmp(&point2.probability))
                        .map(|point| point.t_doubling.round()),
                    None => None,
                }
            } else {
                None
            }
        };

        let strain_summary = StrainSummary {
            code: strain.code.to_owned(),
            name: strain.name.to_owned(),
            description: strain.description.to_owned(),
            styles: strain.styles.to_owned(),
            liquid_appearance: LiquidAppearance::average(appearance_liquid_set),
            yeast_appearance: YeastAppearance::average(appearance_yeast_set),
            organoleptic_score_sets: OrganolepticScoreSets::average(scores),
            clustering: Clustering::average(clustering_set),
            real_attenuation: Uncertain::mean_with_standard_deviation(&real_attenuation_set),
            thermal_test,
            average_peak_doubling_time_minutes,
        };

        match toml::to_string(&strain_summary) {
            Ok(tomled) => {
                let directory = format!("{OUTPUT}/{STRAINS_DIR}");
                if !Path::new(&directory).exists() {
                    if let Err(e) = std::fs::create_dir_all(&directory) {
                        println!("Error making strain summary data directory: {e}")
                    }
                }
                if let Err(e) =
                    std::fs::write(format!("{directory}/{}_summary.toml", strain.code), tomled)
                {
                    println!("Error writing summary for {}: {e}", strain.code)
                }
            }
            Err(e) => println!(
                "Error converting summary for {} into toml format: {e}",
                strain.code
            ),
        }
    }
}

pub fn collect_concentration_fits(
    kinetics_liquid_set: &[&KineticsLiquid],
    strain_code: &str,
    write_log: bool,
) -> Vec<ConcentrationFitLogistics> {
    kinetics_liquid_set
            .par_iter()
            .filter_map(|kinetics_liquid| {
                let mut good_cached_data = None;
                if let Ok(tomled_data) = std::fs::read_to_string(format!("{OUTPUT}/{CACHES_DIR}/{strain_code}/{}.toml", kinetics_liquid.id)) {
                    if let Ok(concentration_fit) = toml::from_str::<ConcentrationFitLogistics>(&tomled_data) {
                        if concentration_fit.fit_params == FitParams::default() && concentration_fit.id == kinetics_liquid.id {
                            if write_log {println!("Extracted good cache for {}", kinetics_liquid.id)}
                            good_cached_data = Some(concentration_fit);
                        }
                    }
                }
                good_cached_data.or_else( || {
                    if let Some(concentration_fit) = kinetics_liquid.concentration_fit_logistics() {
                        if !concentration_fit.points.iter().any(|element| element.probability.is_nan()) {
                            println!("Caching data for {}", kinetics_liquid.id);
                            match toml::to_string(&concentration_fit) {
                                Ok(tomled) => {
                                    let directory = format!("{OUTPUT}/{CACHES_DIR}/{strain_code}");
                                    if !Path::new(&directory).exists() {
                                        if let Err(e) = std::fs::create_dir_all(directory) {
                                            println!("Error making cached data directory for {strain_code}: {e}")
                                        }
                                    }
                                    if let Err(e) = std::fs::write(
                                        format!("{OUTPUT}/{CACHES_DIR}/{strain_code}/{}.toml", kinetics_liquid.id),
                                        tomled,
                                    ) {
                                        println!("Error writing cached fit data for {}: {e}", kinetics_liquid.id)
                                    }
                                }
                                Err(e) => println!(
                                    "Error converting fit data for {} into toml format: {e}", kinetics_liquid.id
                                ),
                            }
                            Some(concentration_fit)
                        } else {
                            if write_log {println!("SKIPPING FILE {}, NaN encountered", kinetics_liquid.id)}
                            None
                        }
                    }
                    else {
                        if write_log {println!("No valid concentration fit for {}", kinetics_liquid.id)}
                        None
                    }
                })
            })
            .collect()
}

pub fn plot_heatmap(fit_points_mixed: &[FitPoint], directory: &str, filename: &str) {
    let precision_max = fit_points_mixed
        .iter()
        .max_by(|fit_point_1, fit_point_2| fit_point_1.precision.total_cmp(&fit_point_2.precision))
        .expect("combined fit points set is not empty")
        .precision;
    let probability_max = fit_points_mixed
        .iter()
        .max_by(|fit_point_1, fit_point_2| {
            fit_point_1.probability.total_cmp(&fit_point_2.probability)
        })
        .expect("combined fit points set is not empty")
        .probability;

    if !Path::new(&directory).exists() {
        if let Err(e) = std::fs::create_dir_all(directory) {
            println!("Error making directory {directory}: {e}")
        }
    }
    let output_name = format!("{directory}/{filename}");
    let root_drawing_area = SVGBackend::new(
        &output_name,
        (
            SQUARE_SIDE * N_POINTS_T_DOUBLING as u32 + 2 * MARGIN,
            SQUARE_SIDE * fit_points_mixed.len() as u32 / N_POINTS_T_DOUBLING as u32 + 2 * MARGIN,
        ),
    )
    .into_drawing_area();
    root_drawing_area.fill(&WHITE).unwrap();

    let t_doubling_half_square_side =
        3.0 * T_DOUBLING_VARIANCE_PRIOR / (N_POINTS_T_DOUBLING - 1) as f64;
    let precision_half_square_side =
        0.5 * PRECISION_SHAPE * PRECISION_SCALE / (N_POINTS_BEFORE_MEAN_PRECISION - 1) as f64;

    let mut ctx = ChartBuilder::on(&root_drawing_area)
        .set_label_area_size(LabelAreaPosition::Left, MARGIN)
        .set_label_area_size(LabelAreaPosition::Right, MARGIN)
        .set_label_area_size(LabelAreaPosition::Bottom, MARGIN)
        .set_label_area_size(LabelAreaPosition::Top, MARGIN)
        .build_cartesian_2d(
            T_DOUBLING_MEAN_PRIOR - 3.0 * T_DOUBLING_VARIANCE_PRIOR - t_doubling_half_square_side
                ..T_DOUBLING_MEAN_PRIOR
                    + 3.0 * T_DOUBLING_VARIANCE_PRIOR
                    + t_doubling_half_square_side,
            precision_half_square_side..precision_max + precision_half_square_side,
        )
        .unwrap();

    ctx.configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .x_desc("doubling time, min")
        .axis_desc_style(("sans-serif", 40))
        .x_label_style(("sans-serif", 20))
        .y_desc("precision")
        .y_label_style(("sans-serif", 20))
        .draw()
        .unwrap();

    let plotting_area = ctx.plotting_area();

    for fit_point in fit_points_mixed.iter() {
        let color_value = (fit_point.probability * (u8::MAX as f64) / probability_max) as u8;
        plotting_area
            .draw(&Rectangle::new(
                [
                    (
                        fit_point.t_doubling - t_doubling_half_square_side,
                        fit_point.precision + precision_half_square_side,
                    ),
                    (
                        fit_point.t_doubling + t_doubling_half_square_side,
                        fit_point.precision - precision_half_square_side,
                    ),
                ],
                ShapeStyle {
                    color: RGBAColor(color_value, color_value, color_value, 1.0),
                    filled: true,
                    stroke_width: 0,
                },
            ))
            .unwrap();
    }
}

pub enum MixerInput<'a> {
    NormalizedFitAndRealData {
        normalized_fit: &'a [FitPoint],
        real_data: &'a [&'a ConcentrationFitLogistics],
    },
    RealData(&'a [&'a ConcentrationFitLogistics]),
    SisterSetDummy(&'a [FitPoint]),
}

/// Mixes fit points and afterwards applies fit parameters probability
/// correction for each point.
pub fn mix_points(mixer_input: MixerInput<'_>) -> Vec<FitPoint> {
    match mixer_input {
        MixerInput::NormalizedFitAndRealData {
            normalized_fit,
            real_data,
        } => {
            if !real_data.is_empty() {
                let len = real_data[0].points.len();
                let mut fit_points_mixed = normalized_fit.to_owned();
                for set_element in real_data.iter() {
                    assert_eq!(set_element.points.len(), len);
                    for (i, fit_points_mixed_element) in fit_points_mixed.iter_mut().enumerate() {
                        fit_points_mixed_element.probability *= set_element.points[i].probability;
                    }
                }
                fit_points_mixed
            } else {
                normalized_fit.to_owned()
            }
        }
        MixerInput::RealData(concentration_fit_logistics_set) => {
            assert!(
                !concentration_fit_logistics_set.is_empty(),
                "expect concentration fit logistics set to be not empty at this point"
            );
            let len = concentration_fit_logistics_set[0].points.len();
            let mut fit_points_mixed: Vec<FitPoint> = concentration_fit_logistics_set[0]
                .points
                .iter()
                .map(|fit_point| FitPoint {
                    t_doubling: fit_point.t_doubling,
                    precision: fit_point.precision,
                    probability: t_doubling_probability(fit_point.t_doubling)
                        * precision_probability(1.0 / (fit_point.precision).sqrt()),
                })
                .collect();
            for set_element in concentration_fit_logistics_set.iter() {
                assert_eq!(set_element.points.len(), len);
                for (i, fit_points_mixed_element) in fit_points_mixed.iter_mut().enumerate() {
                    fit_points_mixed_element.probability *= set_element.points[i].probability;
                }
            }
            fit_points_mixed
        }
        MixerInput::SisterSetDummy(sister_set) => sister_set
            .iter()
            .map(|fit_point| FitPoint {
                t_doubling: fit_point.t_doubling,
                precision: fit_point.precision,
                probability: t_doubling_probability(fit_point.t_doubling)
                    * precision_probability(1.0 / (fit_point.precision).sqrt()),
            })
            .collect(),
    }
}

pub fn mixed_points_if_strain_is_uniform(
    concentration_fit_groups_uniformity: &[Vec<ConcentrationFitLogistics>],
    strain_code: &str,
    write_log: bool,
) -> Option<Vec<FitPoint>> {
    if !concentration_fit_groups_uniformity.is_empty() {
        let mut output_fit_points = None;
        let first_available_set_points = &concentration_fit_groups_uniformity
                .first()
                .expect("concentration fit groups set is not empty - just checked - first indexing always works")
                .first()
                .expect("empty groups not allowed - second indexing always works")
                .points;
        let points_len = first_available_set_points.len();
        let set_len = concentration_fit_groups_uniformity.iter().flatten().count();

        let mut ids_and_probability: Vec<(Vec<&str>, f64)> = Vec::new();
        for group_set_1_indices in (1..concentration_fit_groups_uniformity.len()).powerset() {
            let mut ids_subset_1: Vec<&str> = Vec::new();
            let mut subset1: Vec<&ConcentrationFitLogistics> = Vec::new();
            let mut subset2: Vec<&ConcentrationFitLogistics> = Vec::new();
            for (group_index, group) in concentration_fit_groups_uniformity.iter().enumerate() {
                for concentration_fit in group.iter() {
                    assert_eq!(concentration_fit.points.len(), points_len);
                    if group_index == 0 || group_set_1_indices.contains(&group_index) {
                        subset1.push(concentration_fit);
                        ids_subset_1.push(&concentration_fit.id);
                    } else {
                        subset2.push(concentration_fit);
                    }
                }
            }
            let mixed_points_1 = mix_points(MixerInput::RealData(&subset1));
            let mixed_points_2 = {
                if subset2.is_empty() {
                    mix_points(MixerInput::SisterSetDummy(&mixed_points_1))
                } else {
                    mix_points(MixerInput::RealData(&subset2))
                }
            };
            let sum1 = mixed_points_1
                .iter()
                .fold(0f64, |acc, x| acc + x.probability);
            let sum2 = mixed_points_2
                .iter()
                .fold(0f64, |acc, x| acc + x.probability);
            let sum_mult = sum1 * sum2;
            ids_and_probability.push((ids_subset_1, sum_mult));
            if subset2.is_empty() {
                output_fit_points = Some(mixed_points_1)
            }
        }

        let output_fit_points = output_fit_points.expect(
            "must have encountered a variant with all elements in subset1 and empty subset2",
        );

        let most_probable = ids_and_probability
            .iter()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .unwrap();
        if most_probable.0.len() == set_len {
            if most_probable.1 != 0.0 {
                if write_log {
                    println!(
                        "{strain_code} uniformity is good! highest probability is {:e}",
                        most_probable.1
                    )
                }
                Some(output_fit_points)
            } else {
                if write_log {
                    println!(
                        "ERROR in {strain_code} uniformity! highest probability is {:e}",
                        most_probable.1
                    )
                }
                None
            }
        } else {
            if write_log {
                println!("{strain_code} uniformity is suspicious; most probable: {most_probable:?}\nwhole set: {ids_and_probability:?}")
            }
            None
        }
    } else {
        if write_log {
            println!("Insufficient data for uniformity conclusion")
        }
        None
    }
}

pub fn mixed_points_if_strain_is_stable_over_time(
    mixed_normalized_points_uniformity_set: Vec<FitPoint>,
    concentration_fit_groups_qa: &[Vec<ConcentrationFitLogistics>],
    strain_code: &str,
    write_log: bool,
) -> Option<Vec<FitPoint>> {
    if !concentration_fit_groups_qa.is_empty() {
        let mut output_fit_points = None;
        let mut last_index_into_subset1_and_probability: Vec<(usize, f64)> = Vec::new();
        let full_length = concentration_fit_groups_qa.len();

        for last_index_into_subset1 in 0..=full_length {
            let subset1: Vec<&ConcentrationFitLogistics> = concentration_fit_groups_qa
                [..last_index_into_subset1]
                .iter()
                .flatten()
                .collect();
            let subset2: Vec<&ConcentrationFitLogistics> = concentration_fit_groups_qa
                [last_index_into_subset1..]
                .iter()
                .flatten()
                .collect();

            let mixed_points_1 = mix_points(MixerInput::NormalizedFitAndRealData {
                normalized_fit: &mixed_normalized_points_uniformity_set,
                real_data: &subset1,
            });
            let mixed_points_2 = {
                if subset2.is_empty() {
                    mix_points(MixerInput::SisterSetDummy(&mixed_points_1))
                } else {
                    mix_points(MixerInput::RealData(&subset2))
                }
            };
            let sum1 = mixed_points_1
                .iter()
                .fold(0f64, |acc, x| acc + x.probability);
            let sum2 = mixed_points_2
                .iter()
                .fold(0f64, |acc, x| acc + x.probability);
            let sum_mult = sum1 * sum2;
            last_index_into_subset1_and_probability.push((last_index_into_subset1, sum_mult));
            if subset2.is_empty() {
                output_fit_points = Some(mixed_points_1)
            }
        }

        let output_fit_points = output_fit_points.expect(
            "must have encountered a variant with all elements in subset1 and empty subset2",
        );

        let most_probable = last_index_into_subset1_and_probability
            .iter()
            .max_by(|a, b| a.1.total_cmp(&b.1))
            .expect("always has at least one element");
        let last_index_into_subset1 = most_probable.0;
        if last_index_into_subset1 == full_length {
            if most_probable.1 != 0.0 {
                if write_log {
                    println!(
                    "{strain_code} remains good for all tested slant renewals! highest probability is {:e}",
                    most_probable.1
                )
                }
                Some(output_fit_points)
            } else {
                if write_log {
                    println!(
                        "ERROR in {strain_code} variability! highest probability is {:e}",
                        most_probable.1
                    )
                }
                None
            }
        } else {
            let edge_set: Vec<&str> = concentration_fit_groups_qa
                .get(last_index_into_subset1)
                .expect("at this point last index is below full length and thus accessible")
                .iter()
                .map(|fit_group_element| fit_group_element.id.as_ref())
                .collect();
            if write_log {
                println!("{strain_code} variability is suspicious; change in kinetics fit parameters suspected for measurement set {edge_set:?} and afterwards\nwhole set (last_index_into_subset1, probability): {last_index_into_subset1_and_probability:?}")
            }
            None
        }
    } else {
        if write_log {
            println!("No kinetics data outside uniformity testing is available for {strain_code}")
        }
        Some(mixed_normalized_points_uniformity_set)
    }
}

pub const SQUARE_SIDE: u32 = 70;
pub const MARGIN: u32 = 3 * SQUARE_SIDE;

pub const OUTPUT: &str = "output";

pub const CACHES_DIR: &str = "caches";
pub const HEATMAPS: &str = "heatmaps";
pub const UNIFORMITY: &str = "uniformity";
pub const DATA_FIT: &str = "data_fit";
pub const GENEALOGY_DIR: &str = "genealogy_plots";
pub const ALL_PLOTS: &str = "all_plots";

pub const STRAINS_DIR: &str = "strain_summaries";

fn is_slant_backed_up(di_graph_map_whole: &DiGraphMap<Plaque, Empty>, slant_node: Plaque) -> bool {
    for neighbor in di_graph_map_whole.neighbors_directed(slant_node, Direction::Outgoing) {
        if let PlaqueKind::SlantProd = neighbor.kind {
            return true;
        }
    }
    false
}

lazy_static! {
    static ref REG_SAMPLE_NAME: Regex = Regex::new(r#"^(?P<date>\d{4}-\d{2}-\d{2})-.+$"#)
        .expect("constructed from checked static value");
}

pub const AGE_WEEKS_ATTN: i64 = 16;
pub const AGE_WEEKS_URGENT: i64 = 24;

#[derive(Debug, Eq, PartialEq)]
pub enum OrganolepticsDone {
    Failed,
    None,
    Passed,
}

#[derive(Debug, Eq, PartialEq)]
pub enum AgeChecked {
    Recent,
    RenewSoon,
    RenewImmediately,
}

#[derive(Debug)]
pub struct StrainStatus<'a> {
    pub strain_code: &'a str,
    pub id_most_recent: &'a str,
    pub is_backed_up_most_recent: bool,
    pub organoleptics_done_most_recent: OrganolepticsDone,
    pub age_most_recent: AgeChecked,
    pub previous: Vec<SlantOld<'a>>,
}

impl<'a> StrainStatus<'a> {
    pub fn last_good_organoleptics_in_previous(&self) -> Option<&'a str> {
        for slant_old in self.previous.iter() {
            if let OrganolepticsDone::Passed = slant_old.organoleptics_done {
                return Some(slant_old.id);
            }
        }
        None
    }
}

impl<'a> Display for StrainStatus<'a> {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        match self.age_most_recent {
            AgeChecked::Recent => match self.organoleptics_done_most_recent {
                OrganolepticsDone::Failed => match self.last_good_organoleptics_in_previous() {
                    Some(last_good_id) => write!(f, "{}: latest slant {}, ATTN: organoleptics FAILED, redo using slant {}", self.strain_code, self.id_most_recent, last_good_id),
                    None => write!(f, "{}: latest slant {}, ATTN: organoleptics FAILED, strain was never purified before", self.strain_code, self.id_most_recent),
                },
                OrganolepticsDone::None => match self.last_good_organoleptics_in_previous() {
                    Some(last_good_id) => write!(f, "{}: latest slant {}, ATTN: organoleptics missing, last slant verified organoleptically {}", self.strain_code, self.id_most_recent, last_good_id),
                    None => write!(f, "{}: latest slant {}, ATTN: organoleptics missing, strain was never purified before", self.strain_code, self.id_most_recent),
                },
                OrganolepticsDone::Passed => {
                    if self.is_backed_up_most_recent {
                        write!(f, "{}: latest slant {}, no action needed", self.strain_code, self.id_most_recent)
                    } else {
                        write!(f, "{}: latest slant {}, ATTN: backup", self.strain_code, self.id_most_recent)
                    }
                }
            },
            AgeChecked::RenewSoon => match self.organoleptics_done_most_recent {
                OrganolepticsDone::Failed => match self.last_good_organoleptics_in_previous() {
                    Some(last_good_id) => write!(f, "{}: latest slant {}, ATTN: renew soon, use slant {}", self.strain_code, self.id_most_recent, last_good_id),
                    None => write!(f, "{}: latest slant {}, ATTN: renew soon, strain was never purified before", self.strain_code, self.id_most_recent),
                },
                OrganolepticsDone::None => match self.last_good_organoleptics_in_previous() {
                    Some(last_good_id) => write!(f, "{}: latest slant {}, ATTN: renew soon, last slant verified organoleptically {}", self.strain_code, self.id_most_recent, last_good_id),
                    None => write!(f, "{}: latest slant {}, ATTN: renew soon, strain was never purified before", self.strain_code, self.id_most_recent),
                },
                OrganolepticsDone::Passed => {
                    if self.is_backed_up_most_recent {
                        write!(f, "{}: latest slant {}, ATTN: renew soon", self.strain_code, self.id_most_recent)
                    } else {
                        write!(f, "{}: latest slant {}, ATTN: renew soon, no backup exists", self.strain_code, self.id_most_recent)
                    }
                }
            },
            AgeChecked::RenewImmediately => match self.organoleptics_done_most_recent {
                OrganolepticsDone::Failed => match self.last_good_organoleptics_in_previous() {
                    Some(last_good_id) => write!(f, "{}: latest slant {}, ATTN: renew ASAP, use slant {}", self.strain_code, self.id_most_recent, last_good_id),
                    None => write!(f, "{}: latest slant {}, ATTN: renew ASAP, strain was never purified before", self.strain_code, self.id_most_recent),
                },
                OrganolepticsDone::None => match self.last_good_organoleptics_in_previous() {
                    Some(last_good_id) => write!(f, "{}: latest slant {}, ATTN: renew ASAP, last slant verified organoleptically {}", self.strain_code, self.id_most_recent, last_good_id),
                    None => write!(f, "{}: latest slant {}, ATTN: renew ASAP, strain was never purified before", self.strain_code, self.id_most_recent),
                },
                OrganolepticsDone::Passed => {
                    if self.is_backed_up_most_recent {
                        write!(f, "{}: latest slant {}, ATTN: renew ASAP", self.strain_code, self.id_most_recent)
                    } else {
                        write!(f, "{}: latest slant {}, ATTN: renew ASAP, no backup exists", self.strain_code, self.id_most_recent)
                    }
                }
            },
        }
    }
}

#[derive(Debug)]
pub struct SlantOld<'a> {
    pub id: &'a str,
    pub is_backed_up: bool,
    pub organoleptics_done: OrganolepticsDone,
}

pub fn full_graph(graph_content: String) -> String {
    format!(
        "
digraph {{
    rankdir=LR
    node [style=filled]

    {}

{graph_content}
}}",
        legend()
    )
}

pub fn legend() -> String {
    format!("subgraph cluster_legend {{
        label=\"Legend\"
        color = black;
        fontsize = 20;
        penwidth = 3;
        legend1 [ label = \"Slant\", shape = \"{SLANT_SHAPE}\", style = filled, fillcolor = \"{SLANT_COLOR}\" ];
        legend2 [ label = \"Slant\nProduction\", shape = \"{SLANT_SHAPE}\", style = filled, fillcolor = \"{SLANT_PROD_COLOR}\" ];
        legend3 [ label = \"Plate\", shape = \"{PLATE_SHAPE}\", style = filled, fillcolor = \"{PLATE_REGULAR_COLOR}\" ];
        legend4 [ label = \"Ale 37C\", shape = \"{PLATE_SHAPE}\", style = filled, fillcolor = \"{PLATE_THERMAL_ALE}\" ];
        legend5 [ label = \"Lager 37C\", shape = \"{PLATE_SHAPE}\", style = filled, fillcolor = \"{PLATE_THERMAL_LAGER}\" ];
        legend6 [ label = \"Liquid\", shape = \"{LIQUID_SHAPE}\", style = filled, fillcolor = \"{LIQUID_COLOR}\" ];
        legend7 [ label = \"Liquid\ncold\", shape = \"{LIQUID_SHAPE}\", style = filled, fillcolor = \"{LIQUID_COLOR_LAGER}\" ];
        legend8 [ label = \"Packaged\", shape = \"{PACKAGE_SHAPE}\", style = filled, fillcolor = \"{PACKAGE_COLOR}\" ];
        legend9 [ label = \"Pitch\", shape = \"{PACKAGE_SHAPE}\", style = filled, fillcolor = \"{LIQUID_COLOR}\" ];
        legend10 [ label = \"Pitch\ncold\", shape = \"{PACKAGE_SHAPE}\", style = filled, fillcolor = \"{LIQUID_COLOR_LAGER}\" ];
        legend11 [ label = \"Organoleptic\nPassed\", shape = \"{ORGANOLEPTIC_SHAPE}\", style = filled, fillcolor = \"{ORGANOLEPTIC_PASSED_COLOR}\" ];
        legend12 [ label = \"Organoleptic\nFailed\", shape = \"{ORGANOLEPTIC_SHAPE}\", style = filled, fillcolor = \"{ORGANOLEPTIC_FAILED_COLOR}\" ];
        {{legend1 -> legend2 -> legend3 -> legend4 -> legend5 -> legend6 -> legend7 -> legend8 -> legend9 -> legend10 -> legend11 -> legend12 [style=invis];}}
    }}")
}

#[derive(Clone, Default)]
pub struct Empty;

impl Display for Empty {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "")
    }
}

impl Debug for Empty {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        write!(f, "")
    }
}

#[derive(Copy, Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum PlaqueKind<'a> {
    LiquidPackage { count: u8 },
    LiquidPitch { style: &'a str },
    LiquidPitchLagered { style: &'a str },
    LiquidRegular,
    LiquidRegularLagered,
    OrganolepticFailed,
    OrganolepticPassed,
    PlateRegular,
    PlateThermalAle,
    PlateThermalLager,
    SlantOrigin { strain_code: &'a str },
    SlantProd,
    SlantRegular,
}

#[derive(Copy, Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct Plaque<'a> {
    pub id: &'a str,
    pub kind: PlaqueKind<'a>,
}

pub const LIQUID_SHAPE: &str = "oval";
pub const ORGANOLEPTIC_SHAPE: &str = "house";
pub const PACKAGE_SHAPE: &str = "note";
pub const PLATE_SHAPE: &str = "hexagon";
pub const SLANT_SHAPE: &str = "rectangle";

pub const LIQUID_COLOR: &str = "#F9E076";
pub const LIQUID_COLOR_LAGER: &str = PLATE_THERMAL_LAGER;
pub const ORGANOLEPTIC_FAILED_COLOR: &str = "#CC3300";
pub const ORGANOLEPTIC_PASSED_COLOR: &str = "#9ACD32";
pub const PACKAGE_COLOR: &str = "#DB7093";
pub const PLATE_REGULAR_COLOR: &str = "#FAF0E6";
pub const PLATE_THERMAL_ALE: &str = "#FF7F50";
pub const PLATE_THERMAL_LAGER: &str = "#91B6D4";
pub const SLANT_COLOR: &str = "#CE954B";
pub const SLANT_PROD_COLOR: &str = "#52B2BF";

impl<'a> Plaque<'a> {
    pub fn shape(&self) -> &str {
        match &self.kind {
            PlaqueKind::LiquidPackage { .. }
            | PlaqueKind::LiquidPitch { .. }
            | PlaqueKind::LiquidPitchLagered { .. } => PACKAGE_SHAPE,
            PlaqueKind::LiquidRegular | PlaqueKind::LiquidRegularLagered => LIQUID_SHAPE,
            PlaqueKind::OrganolepticFailed | PlaqueKind::OrganolepticPassed => ORGANOLEPTIC_SHAPE,
            PlaqueKind::PlateRegular
            | PlaqueKind::PlateThermalAle
            | PlaqueKind::PlateThermalLager => PLATE_SHAPE,
            PlaqueKind::SlantOrigin { .. } | PlaqueKind::SlantProd | PlaqueKind::SlantRegular => {
                SLANT_SHAPE
            }
        }
    }
    pub fn color(&self) -> &str {
        match &self.kind {
            PlaqueKind::LiquidPackage { .. } => PACKAGE_COLOR,
            PlaqueKind::LiquidRegular | PlaqueKind::LiquidPitch { .. } => LIQUID_COLOR,
            PlaqueKind::LiquidRegularLagered | PlaqueKind::LiquidPitchLagered { .. } => {
                LIQUID_COLOR_LAGER
            }
            PlaqueKind::OrganolepticFailed => ORGANOLEPTIC_FAILED_COLOR,
            PlaqueKind::OrganolepticPassed => ORGANOLEPTIC_PASSED_COLOR,
            PlaqueKind::PlateRegular => PLATE_REGULAR_COLOR,
            PlaqueKind::PlateThermalAle => PLATE_THERMAL_ALE,
            PlaqueKind::PlateThermalLager => PLATE_THERMAL_LAGER,
            PlaqueKind::SlantOrigin { .. } | PlaqueKind::SlantRegular => SLANT_COLOR,
            PlaqueKind::SlantProd => SLANT_PROD_COLOR,
        }
    }
}

impl<'a> Display for Plaque<'a> {
    fn fmt(&self, f: &mut Formatter) -> FmtResult {
        match self.kind {
            PlaqueKind::LiquidPackage { count } => write!(f, "{}\n{} vial(s)", self.id, count),
            PlaqueKind::LiquidPitch { style } => write!(f, "{}\n{}", self.id, style),
            PlaqueKind::LiquidPitchLagered { style } => write!(f, "{}\n{}", self.id, style),
            PlaqueKind::SlantOrigin { strain_code } => write!(f, "{}\n({})", self.id, strain_code),
            _ => write!(f, "{}", self.id),
        }
    }
}
