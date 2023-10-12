// #![allow(unused_must_use)]
// #![allow(dead_code)]
// TODO: get reference as pitch rate from experiment description file

mod dim_analysis;
use dim_analysis::*;

use core::panic;
use std::fs::{self, OpenOptions,File};
use std::io::{BufWriter, Write};
use plotters::prelude::*;
use toml::{Table, Value};
use time::PrimitiveDateTime;
use std::path::{Path, PathBuf};
use walkdir::WalkDir;
use std::collections::{HashMap, BTreeMap};
use regex::Regex;
use petgraph::{graphmap::DiGraphMap, dot::{Dot, Config}, visit::{Dfs, Reversed, Bfs}, prelude::GraphMap};
use toml::map::Map;
use argmin::core::{State, Error, Executor, CostFunction};
use argmin::solver::neldermead::NelderMead;
use ndarray::Array1;
use std::sync::Mutex; 
use qrcode_generator::QrCodeEcc;
use itertools::Itertools;
use std::process::Command;

lazy_static::lazy_static! { static ref WARNINGS: Mutex<Vec<String>> = Mutex::new(Vec::new()); }

const BASE_URL_FOR_QR_CODES: &str = "https://feature-main.alzymologist-github-io.pages.dev/info/slants/";
const INPUT_DIR: &str = "data/"; 
const OUTPUT_DIR: &str = "output/";
const OUTPUT_COUNT_DIR: &str = "output/count/";
const OUTPUT_DENSITY_DIR: &str = "output/density/";
const OUTPUT_FITTING_DIR: &str = "output/fitting/";
const OUTPUT_GENEALOGY_DIR: &str = "output/genealogy/";
const OUTPUT_QRCODES_DIR: &str = "output/qrcodes/";

const GENEALOGY_NAME: &str = "genealogy";
const YEAST_PAGE_PATH: &str = "../content/info/yeasts.md"; 
const DOT_REPLACEMENT: &str =r##"
digraph {
compound=true;
rankdir=LR;
node [style=filled, colorscheme=prgn6]
    subgraph cluster_legend {
        label="Legend"
        color=black;
        penwidth=3;
        subgraph media {
            label="Media"
            ranksep=0.2;
            rankdir=TB;
            legendA0 [ label = "Medium:", style=filled, fillcolor="#ffffff", shape="plaintext"];
            legendA1 [ label = "missing", style=filled];
            legendA2 [ label = "stock", style=filled, fillcolor="#EE6AA7" ];
            legendA3 [ label = "slant", style=filled, fillcolor=2 ];
            legendA4 [ label = "plate", style=filled, fillcolor=3 ];
            legendA5 [ label = "liquid", style=filled, fillcolor=4 ];
            legendA6 [ label = "organoleptic", style=filled, fillcolor=5 ];
            {legendA0 -> legendA1 -> legendA2 -> legendA3 -> legendA4 -> legendA5 -> legendA6 [style=invis];}
        }
        subgraph protocol {
            label="Protocol"
            ranksep=0.2;
            rankdir=TB;
            legendB0 [ label = "Protocol:", style=filled, fillcolor="#ffffff", shape="plaintext"];
            legendB1 [ label = "revive", style=filled, fillcolor="#ffffff", shape="pentagon"];
            legendB2 [ label = "scale-up", style=filled, fillcolor="#ffffff", shape="octagon"];
            legendB3 [ label = "propagation", style=filled, fillcolor="#ffffff", shape="box"];
            legendB4 [ label = "package", style=filled, fillcolor="#ffffff", shape="box3d"];
            legendB5  [ label = "isolation", style=filled, fillcolor="#ffffff", shape="cylinder"];
            legendB6  [ label = "slant-qa", style=filled, fillcolor="#ffffff", shape="folder"];
            legendB7  [ label = "yeast-uniformity", style=filled, fillcolor="#ffffff", shape="egg"];
            {legendB0 -> legendB1 -> legendB2 -> legendB3 -> legendB4 -> legendB5 -> legendB6 -> legendB7  [style=invis];}
        }
    }
spacer [style=invis, height=3]  // spacer node"##;

// Constants for Nelder Mead problems
const NELDER_MEAD_BOUNDS_C0: (f64, f64) = (0.0, 1E20f64);
const NELDER_MEAD_BOUNDS_CMAX: (f64, f64) = (0.0, 1E20f64);
const NELDER_MEAD_BOUNDS_GS: (f64, f64) = (0.0, 10.0);
const NELDER_MEAD_BOUNDS_3D: [(f64, f64); 3] = [NELDER_MEAD_BOUNDS_C0, NELDER_MEAD_BOUNDS_CMAX, NELDER_MEAD_BOUNDS_GS];
const C0_INITIAL: f64 = 1E9f64;
const C0_DELTA:f64 = 5E9f64;
const CMAX_INITIAL: f64 = 1E14f64;
const CMAX_DELTA: f64 = 5E14f64; 
const GS_INITIAL: f64 = 0.5f64;
const GS_DELTA: f64 = 0.5f64;
const MAX_ITERATIONS: u64 = 10000;

const MIN_POINTS_TO_PLOT_CONC: usize = 3;

type Nodes = HashMap<String, Map<String, Value>>;

fn log_warnings(s: &str) -> () {
    let mut warnings = WARNINGS.lock().unwrap();
    warnings.push(s.to_string());
}

fn tomls_into_nodes_and_links (input_dir: &str) -> (Nodes, HashMap<String, String>) {
    let toml_files: Vec<PathBuf> = WalkDir::new(input_dir)
        .into_iter()
        .filter_map(|entry| entry.ok())
        .filter(|entry| 
            entry.file_type().is_file() &&
            entry.file_name().to_string_lossy().ends_with(".toml"))
        .map(|entry| entry.into_path() )
        .collect();
        println!("TOML files found: {:?}", toml_files.len());

    let mut checked_nodes: Nodes = HashMap::new();
    let mut local_links_to_nodes: HashMap<String, String> = HashMap::new();

    for file in toml_files {
        let contents = fs::read_to_string(&file.as_path()).unwrap();
        if let Ok(mut toml_map) = contents.parse::<Table>() {
            let id = try_to_read_field_as_string(&toml_map, "id");
            let time = try_to_read_reference_time(&toml_map);
            let medium = try_to_read_field_as_string(&toml_map, "medium");
            let parent = try_to_read_field_as_vec(&toml_map, "parent");
            let participants = try_to_read_field_as_vec(&toml_map, "participants");
            match (&id, medium, time, parent, participants) {
                (None, _, _, _, _) => { log_warnings( &format!("{:<25} {:<25} {:?}", "Id (at least)", "missing from", file) );},
                (_, None, _, _, _) => { log_warnings( &format!("{:<25} {:<25} {:?}", "Medium (at least)", "missing from", file) );}
                (_, _, None, _, _) => { log_warnings( &format!("{:<25} {:<25} {:?}", "Time (at least)", "missing from", file) );}
                (_, _, _, Some(_), Some(_)) => { log_warnings( &format!("{:<25} {:<25} {:?}", "Parent and participants", "can't be both present", file))}
                (Some(_), Some(_), Some(_), Some(_), None) |
                (Some(_), Some(_), Some(_), None, Some(_)) => {
                    checked_nodes.insert(id.clone().unwrap(), toml_map);
                    local_links_to_nodes.insert(id.unwrap(), String::from(file.canonicalize().unwrap().to_str().unwrap()));
                },
                (Some(_), Some(m), Some(_), None, None)  => {
                    if m == "stock" {
                        let new = Value::Array(vec![Value::String("new".into())]); 
                        toml_map.insert(String::from("parent"),new); // Modifing toml to mark impilictly that sample is new
                        checked_nodes.insert(id.clone().unwrap(), toml_map);
                        local_links_to_nodes.insert(id.unwrap(), String::from(file.canonicalize().unwrap().to_str().unwrap()));
                    } else { log_warnings( &format!("{:<25} {:<25} {:?}", "Parent or participants", "missing from", file));}
                },
            }
        }
    } 
    println!("TOML files converted into checked nodes: {:?}", checked_nodes.len());
    (checked_nodes, local_links_to_nodes)
}

fn nodes_into_components(nodes: Nodes, dust_mode: bool) -> HashMap<String, Nodes> {
    // Splits nodes into connectivity components (dust_mode == false) or just dust (dust_mode == true).
    // With (dust_mode == true) each component recieves exactly one Node.
    // Thus (dust_mode == true) allows to fit each Node independently when we later fit components.
    let mut components: HashMap<String, Nodes> = HashMap::new();
    let graph = build_digraphmap(nodes.clone());
    let reversed_graph = Reversed(&graph);

    for (id, toml_map) in nodes.clone().into_iter() {
        let mut farthest_ancestor_id = id.clone(); // Trivial ancestors id;
        
        if !dust_mode { 
            let mut bfs = Bfs::new(&reversed_graph, &id);
            while let Some(next_node_id) = bfs.next(&reversed_graph) {
                farthest_ancestor_id = next_node_id.to_owned(); 
            }
        }
        components.entry(farthest_ancestor_id)
            .or_insert_with(HashMap::new)
            .insert(id, toml_map);
    }
    println!("(Connectivity) components created from checked nodes: {}", components.len());
    log_component_separation_results(&components);
    components
}

fn components_into_problems(components: HashMap<String, Nodes>) -> HashMap<String, NelderMeadComponentProblem> {
    // Prepares data for the Nelder-Mead Algorithm. Does not run the optimiziation itself. 
    let mut problems_for_components: HashMap<String, NelderMeadComponentProblem> = HashMap::new();
    for (component_ancestor_id, nodes) in components {

        let mut problems_for_this_component: HashMap<String, NelderMeadSingleProblem> = HashMap::new();

        for (id, toml_map) in nodes.into_iter() {
            let experiment_res = toml_map.clone().try_into::<UniformityExperiment>();
            let reference_time = try_to_read_reference_time(&toml_map).unwrap();
            
            if let Ok(res) = experiment_res { 
                if let Some(measurements) = res.measurement {  // TOML can contain no measurements
                    let points: Vec<NelderMeadPoint> = measurements
                        .into_iter()
                        .filter_map(|measurement| read_uniformity_point(measurement).ok())
                        .filter(|point| point.concentration.is_some() && point.density.is_some())
                        .map (|point| NelderMeadPoint {
                            conc: point.concentration.unwrap(),
                            density: point.density.unwrap(),
                            hours: (point.timestamp - reference_time).as_seconds_f64()/(60.0*60.0)})
                        .collect();
    
                    if points.len() >= MIN_POINTS_TO_PLOT_CONC {
                        let nm_problem_raw = NelderMeadSingleProblem {
                            points,
                            reference_time,
                        }; 
                        problems_for_this_component.insert(id, nm_problem_raw);
                    }
                }
            } 
        }
        
        if problems_for_this_component.len() != 0 {
            let dimensions = problems_for_this_component.len()*2 + 1;
            let initial_simplex = initial_simplex_nd(dimensions); 

            problems_for_components.insert(component_ancestor_id,
            NelderMeadComponentProblem {
                problems: problems_for_this_component,
                initial_simplex,
            });
        }
    }
    println!("(Connectivity) components converted into optimization problems: {:?}", problems_for_components.len());
    problems_for_components
} 

fn log_component_separation_results(components: &HashMap<String, Nodes>) -> () {
    fs::create_dir_all(OUTPUT_DIR).expect("Failed to create directory.");
    if let Ok(file) = File::create(format!("{}/components.txt", OUTPUT_DIR)) {
        let mut buffer = BufWriter::new(file);
        for (ancestor, component) in &components.clone() {
            writeln!(buffer, "Component ancestor ID: {}", ancestor).unwrap();
            let mut ids: Vec<_> = component.keys().cloned().collect::<Vec<_>>();
            ids.sort();
            writeln!(buffer, "{},\n", ids.join(", ")).unwrap();
        }
        buffer.flush().unwrap();
    }
}

fn build_digraphmap(nodes: Nodes) -> DiGraphMap<&'static str, ()> {
    // DiGraphMap<&'static str, ()>  has () because edges do not carry additional information.
    let edges: Vec<(&'static str, &'static str)> = nodes.values().flat_map(|mapping| {
        let id_string = try_to_read_field_as_string(&mapping, "id").unwrap();
        let id: &'static str = Box::leak(id_string.into_boxed_str());
        let parent_or_participants: Vec<String> = { 
            let parent = try_to_read_field_as_vec(mapping, "parent");
            let participant = try_to_read_field_as_vec(mapping, "participants");
            parent.or(participant).unwrap_or_else(|| {
                panic!("Node with id {:?} does not have parent or participants. It should have been checked before.", id);
            })
        };
    
        parent_or_participants
            .into_iter()
            .map(|par_id_string| Box::leak(par_id_string.into_boxed_str()) as &'static str)
            .map(|par_id| (par_id, id))
            .filter(|(par_id, _)| *par_id != "new" ) // Filter out edges with "new" parent or participant
            .collect::<Vec<_>>()
    }).collect();

    let mut graph = DiGraphMap::<&str, ()>::from_edges(edges);
    for mapping in nodes.values() { // Add nodes with no edges
        let id = Box::leak(try_to_read_field_as_string(&mapping, "id").unwrap().into_boxed_str()) as &'static str;
        if !graph.contains_node(id) { graph.add_node(id); }
    }
    graph
}

fn initial_simplex_nd(number_of_dimensions: usize) -> Vec<Vec<f64>> {
    // Creates initial simplex for the Nelder Mead Problem.
    // Space in which simplex is created must be 2N+1 dimensional, where N is a positive integer >= 1, e.g. 3, 5, 7 ...
    // Simplex will have 2N+2 vertices.
    assert!(number_of_dimensions >= 3 && number_of_dimensions % 2 == 1);
    let mut simplex: Vec<Vec<f64>> = vec![];
    let number_of_vertices = number_of_dimensions + 1;
    for vertex_number in 0..number_of_vertices {

        let mut base_vertex: Vec<f64> = (0..number_of_dimensions-1) 
            .map(|dim| if dim % 2 == 0 {C0_INITIAL} else {CMAX_INITIAL})
            .chain(std::iter::once(GS_INITIAL))
            .collect();
        
        // Simplex must be non-degenerate — all vertices must be unique.
        // So we introduce different perturbations to each copy of the base_vector.
        // Exact values of initial parameters and their perturbations do not matter.
        match vertex_number {
            0 => {},
            v if v == number_of_vertices - 1 => {base_vertex[v - 1] += GS_DELTA},
            v if v % 2 == 1 => {base_vertex[v - 1] += C0_DELTA},
            v if v % 2 == 0 => {base_vertex[v - 1] += CMAX_DELTA},
            _ => {panic!("Should be unreachable.")}
        }
        simplex.push(base_vertex);
    }
    simplex
}

#[derive(Debug, Clone)]
struct NelderMeadPoint {
    conc: O32<UnitDensity>,
    density: O32<MassDensity>,
    hours: f64,
}
#[derive(Debug, Clone)]
struct NelderMeadSingleProblem {
    points: Vec<NelderMeadPoint>,
    reference_time: PrimitiveDateTime,
}

#[derive(Debug, Clone)]
struct NelderMeadComponentProblem {
    problems: HashMap<String, NelderMeadSingleProblem>,
    initial_simplex: Vec<Vec<f64>>, 
}

type NelderMeadComponentSolution = argmin::core::OptimizationResult<NelderMeadComponentProblem, NelderMead<Vec<f64>, f64>, argmin::core::IterState<Vec<f64>, (), (), (), f64>>;

impl CostFunction for NelderMeadSingleProblem {
    type Param = Vec<f64>;
    type Output = f64;
    // params: conc_0, conc_max, growth_speed_max
    fn cost(&self, params: &Self::Param) -> Result<Self::Output, Error> {
        let times = &self.points.clone().into_iter().map(|p|p.hours).collect();
        let modeled_concentrations = model_cell_concentrations(params, times);

        // Cost using Mean Squared Logarithmic Error (MSLE)
        let cost = modeled_concentrations.into_iter().zip(self.points.clone().into_iter())
            .map(|(modeled_conc, point)| (modeled_conc, point.conc.v as f64, point.conc.e as f64))
            .filter(|&(modeled, actual, error)| modeled > 0.0 && actual > 0.0 && error > 0.0) 
            .map(|(modeled, actual, _error)|(modeled.log10() - actual.log10()).powi(2))
            .sum::<f64>() / self.points.len() as f64;

        // If parameters are out of bounds, the penalty becomes non-zero:
        let mut penalty: f64 = 0.0;
        let penalty_factor = 1E30f64;
        let bounds = &NELDER_MEAD_BOUNDS_3D[..];
        for (param, &(lower_bound, upper_bound)) in params.iter().zip(bounds) {
            if *param < lower_bound {
                penalty += penalty_factor * (*param - lower_bound).powi(2);
            }
            if *param > upper_bound {
                penalty += penalty_factor * (*param - upper_bound).powi(2);
            }
        }
        Ok(cost + penalty)
    }
}

// Cost for the NelderMeadComponentProblem is defined using cost for the NelderMeadSingleProblem 
impl CostFunction for NelderMeadComponentProblem {
    type Param = Vec<f64>;
    type Output = f64;

    fn cost(&self, params: &Self::Param) -> Result<Self::Output, Error> {
        assert!((params.len() - 1) % 2 == 0);
        let mut cost_over_component = 0.0;

        let unique_chunks = params[..(params.len() - 1)].chunks(2); //?
        let common_param = params.last().unwrap();

        for (problem, unique_chunk) in self.problems.values().zip(unique_chunks) {
            let mut combined_params = unique_chunk.to_vec();
            combined_params.push(*common_param);

            let problem_cost = problem.cost(&combined_params)?; 
            cost_over_component += problem_cost;
        }

        Ok(cost_over_component)
    }
}

fn run_nelder_mead_for_component(cp: NelderMeadComponentProblem) -> Result<NelderMeadComponentSolution, Error> { 
    let initial_simplex: NelderMead<Vec<f64>, f64> =  NelderMead::new(cp.clone().initial_simplex);
    let result = Executor::new(cp, initial_simplex)
        .configure(|state| state
            .max_iters(MAX_ITERATIONS)
            .target_cost(0.0))
        .run();
    
    result
}

fn log_nelder_mead_solutions(component_id: &String, nm: &NelderMeadComponentProblem, nm_solution: &NelderMeadComponentSolution) -> () {
    fs::create_dir_all(OUTPUT_DIR).expect("Failed to create directory.");
    fs::create_dir_all(OUTPUT_COUNT_DIR).expect("Failed to create directory.");
    fs::create_dir_all(OUTPUT_DENSITY_DIR).expect("Failed to create directory.");
    fs::create_dir_all(OUTPUT_FITTING_DIR).expect("Failed to create directory."); 

    let component_cost = nm_solution.state.cost;
    let component_params = nm_solution.state().get_best_param().unwrap();
    let component_params_for_printing: Vec<String> = component_params.iter().map(|&param| format!("{:.3e}", param)).collect();
    let mut unique_param_chunks = component_params[..(component_params.len() - 1)].chunks(2); // Split vector of component parameters to find parameters for this primitive problem

    if let Ok(file) = File::create(format!("{}/fitting-{}.txt", OUTPUT_FITTING_DIR, component_id)) {
        let mut buffer = BufWriter::new(file);
        writeln!(buffer, "Component ancestor ID: {}", component_id).unwrap();
        writeln!(buffer, "Component optimized parameters: {:?}", component_params_for_printing).unwrap();
        writeln!(buffer, "Component optimization cost: {}", component_cost).unwrap();
        writeln!(buffer, "Dimensionality of the optimization problem for this component: {}", nm.initial_simplex[0].len()).unwrap(); 
        writeln!(buffer, "Initial simplex used for the optimization problem for this component: ").unwrap();
        for vertex in &nm.initial_simplex {
            let formatted_row: Vec<String> = vertex.iter().map(|&value| format!("{:.3e}", value)).collect();
            writeln!(buffer, "[{}]", formatted_row.join(", ")).unwrap();
        }
        writeln!(buffer, "").unwrap();
        
        writeln!(buffer, "Parameters for Single Nelder-Mead problems are: [C0, CMax, µmax]").unwrap();
        writeln!(buffer, "C0 — Cell concentration at the reference time (cells/m^3)\nCMax — Maximal cell concentration (cells/m^3)\nµmax - Maximal specific cell growth rate (1/h)\n").unwrap();

            for (single_problem_id, single_problem) in &nm.problems {
                let single_problem_params_unique = unique_param_chunks.next().unwrap_or(&[]);
                let single_problem_params_all = vec![single_problem_params_unique[0], single_problem_params_unique[1], *component_params.last().unwrap()];
                let single_problem_params_for_printing: Vec<String> = single_problem_params_all.iter().map(|&param| format!("{:.3e}", param)).collect();
                let cell_conc_for_printing: Vec<String> = single_problem.clone().points.into_iter().map(|p| format!("{:.3e}", p.conc.v)).collect();
                let times_for_printing: Vec<String> = single_problem.clone().points.into_iter().map(|p| format!("{:.3}", p.hours)).collect();
                
                writeln!(buffer, "Single problem ID: {}", single_problem_id).unwrap();
                writeln!(buffer, "Optimized parameters: {:?}", single_problem_params_for_printing).unwrap();
                writeln!(buffer, "Used times (hours since reference):   {:?}", times_for_printing).unwrap();
                writeln!(buffer, "Used cell concentrations (cells/m^3): {:?}", cell_conc_for_printing).unwrap();
                writeln!(buffer, "Used reference time: {:?}\n", single_problem.reference_time).unwrap();

                let plot_name_count: String = OUTPUT_COUNT_DIR.to_owned() + "count-" + single_problem_id + ".svg";
                plot_count(&single_problem_id, &plot_name_count, &single_problem.clone(), component_cost, single_problem_params_all);

                let plot_name_density: String = OUTPUT_DENSITY_DIR.to_owned() + "density-" + single_problem_id + ".svg";
                plot_density(&single_problem_id, &plot_name_density, &single_problem.clone());
            }
    }
}

fn model_cell_concentrations(params: &[f64], times: &Vec<f64>) -> Vec<f64> {
    let (conc_0, conc_max, growth_speed_max) = (params[0], params[1], params[2]);
    let concentrations = times
        .iter()
        .map(|time| {
            let exp = (growth_speed_max * time).exp();
            (conc_0 * conc_max * exp) / (conc_max - conc_0 + conc_0 * exp)})
        .collect();
    concentrations
}

fn generate_points_with_logistic_model(params: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let number_of_generated_points = 100;
    let min_time = 0f64;
    let max_time = 75f64;
    let created_time_data = Array1::linspace(min_time, max_time, number_of_generated_points).to_vec();
    let created_conc_data = model_cell_concentrations(&params, &created_time_data);
    (created_time_data, created_conc_data)
 }

fn plot_count(id:&str, plot_name: &str, nm: &NelderMeadSingleProblem,  cost_per_datapoint: f64, optimized_params: Vec<f64>) ->  Result<(), DrawingAreaErrorKind<std::io::Error>> {
    let conc_0 = optimized_params[0];
    let conc_max = optimized_params[1];
    let growth_speed_max = optimized_params[2];
    let root_drawing_area = SVGBackend::new(&plot_name, (1024, 768)).into_drawing_area();
    root_drawing_area.fill(&WHITE).unwrap();

    let time_min_shown = 0f64; // hours
    let time_max_shown = 75f64; // hours
    let conc_min_shown = 1E10f64;
    let conc_max_shown = 5E14f64;

    let mut ctx = ChartBuilder::on(&root_drawing_area)
        .set_label_area_size(LabelAreaPosition::Left, 80)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption(format!("{}", id), ("sans-serif", 40))
        .build_cartesian_2d(time_min_shown..time_max_shown, (conc_min_shown..conc_max_shown).log_scale())?;

        ctx.configure_mesh()
            .x_desc("relative time, hours")
            .axis_desc_style(("sans-serif", 40))
            .x_label_formatter(&|x| format!("{}", x))
            .x_label_style(("sans-serif", 20))
            .y_desc("Cell concentration, cells/m^3")
            .y_label_formatter(&|x| format!("{:e}", x))
            .y_label_style(("sans-serif", 20))
            .draw()?;


        root_drawing_area.draw(&Text::new(
            format!("Component cost = {:.3e} ", cost_per_datapoint ),
            (500, 570), 
            ("sans-serif", 30).into_font(),
        ))?;
        root_drawing_area.draw(&Text::new(
            format!("Maximal specific cell growth rate = {:.3} 1/h", growth_speed_max),
            (500, 600), 
            ("sans-serif", 30).into_font(),
        ))?;

        root_drawing_area.draw(&Text::new(
            format!("Initial cell concentration = {:.3e} cells/m^3", conc_0),
            (500, 630), 
            ("sans-serif", 30).into_font(),
        ))?;

        root_drawing_area.draw(&Text::new(
            format!("Max cell concentration = {:.3e} cells/m^3", conc_max),
            (500, 660), 
            ("sans-serif", 30).into_font(),
        ))?;

        ctx.draw_series(
            nm.points.clone().into_iter().filter_map(|p| {
                Some(ErrorBar::new_vertical(
                    p.hours,
                    (p.conc.v - p.conc.e) as f64,
                    p.conc.v as f64,
                    (p.conc.v + p.conc.e) as f64,
                    BLUE.filled(), 10))
            }
        ))?;

    let (created_time_data, created_conc_data) = generate_points_with_logistic_model(&optimized_params);
    ctx.draw_series(LineSeries::new(
        created_time_data.iter().zip(created_conc_data.iter()).map(|(time, conc)| (*time, *conc)),
        &RED,
    ))?;
    // println!("Created plot: {:?}", &plot_name);
    Ok(())
}

fn plot_density(id:&str, plot_name: &str, nm: &NelderMeadSingleProblem) -> Result<(), DrawingAreaErrorKind<std::io::Error>> {
    let root_drawing_area = SVGBackend::new(&plot_name, (1024, 768)).into_drawing_area();
    root_drawing_area.fill(&WHITE).unwrap();

    let time_min_shown = 0f64; // hours
    let time_max_shown = 75f64; // hours
    let density_min_shown = 1000f64;
    let density_max_shown = 1100f64;

    let mut ctx = ChartBuilder::on(&root_drawing_area)
        .set_label_area_size(LabelAreaPosition::Left, 100)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption(id, ("sans-serif", 40))
        .build_cartesian_2d(time_min_shown..time_max_shown, density_min_shown..density_max_shown)?;

    ctx.configure_mesh()
        .x_desc("relative time, hours")
        .axis_desc_style(("sans-serif", 40))
        .x_label_formatter(&|x| format!("{}", x))
        .x_label_style(("sans-serif", 20))
        .y_desc("density, kg/m^3")
        .y_label_style(("sans-serif", 20))
        .draw()?;
    

    ctx.draw_series(
        nm.points.clone().into_iter().filter_map(|p| {
            Some(ErrorBar::new_vertical(
                p.hours,
                (p.density.v - p.density.e) as f64,
                p.density.v as f64,
                (p.density.v + p.density.e) as f64,
                BLUE.filled(), 10))}
    ))?;
    Ok(())
}

fn plot_genealogy(pathname: String, nodes: Nodes, file_links: HashMap<String, String>) {
    //// Uses the graphwiz program, which is called via shell.
    fs::create_dir_all(OUTPUT_GENEALOGY_DIR).expect("Failed to create directory.");
    let graph = build_digraphmap(nodes.clone());
    let dot = Dot::with_config(&graph,&[Config::EdgeNoLabel]);
    let mut content = format!("{:?}", dot);
    let mut dotfile = File::create(pathname.to_owned() + ".dot").unwrap();
    
    let colours_for_media: HashMap<&str, &str> = HashMap::from([
        ("slant", "2"),
        ("plate", "3"),
        ("liquid", "4"),
        ("organoleptic", "5"),
        ("stock", r##" "#EE6AA7" "##),
        ]);

    let shapes_for_protocols: HashMap<&str, &str> = HashMap::from([
        ("yeast-uniformity", "egg"),
        ("isolation", "cylinder"),
        ("revive", "pentagon"),
        ("scale-up", "octagon"),
        ("package", "box3d"),
        ("propagation", "box"),
        ("slant-qa", "folder"),
        ]);

    //// Editing of .dot after generation 
    content = content.replacen("digraph {", DOT_REPLACEMENT, 1);
    let node_pattern = Regex::new(r#"(\d+) \[ label = "\\"(.+?)\\"" \]"#).unwrap();
    for node_captures in node_pattern.captures_iter(&content.clone()) {
        let captured_nodenumber = node_captures.get(1).unwrap().as_str();
        let captured_id = node_captures.get(2).unwrap().as_str();

        if let Some(toml_map) = nodes.get(captured_id)  {
            let medium = toml_map.get("medium").unwrap().as_str().unwrap(); 
            if let Some(colour) = colours_for_media.get(medium) {
                let link = format!("file://{}", file_links.get(captured_id).unwrap());

                let maybe_protocol = try_to_read_field_as_map(&toml_map, "protocol");
                let shape_name = maybe_protocol
                .and_then(|map| map.get("name"))
                .and_then(|v| v.as_str())
                .unwrap_or("ellipse");
            
                let shape = shapes_for_protocols.get(shape_name).unwrap_or(&"ellipse");
                 let pattern_with_colour = format!(r#"{} [label="{}", URL="{}", fillcolor={}, shape={}]"#, captured_nodenumber, captured_id, link, colour, shape);
                content = content.replace(&node_captures[0], &pattern_with_colour);
            }
        } else { //// If captured_id was not found in the graph
                log_warnings( &format!("{:<25} {:<25}", captured_id, "id does not correspond to any converted TOML") );
                let simple_pattern = format!(r#"{} [ label = "{}" ]"#, captured_nodenumber, captured_id);
                content = content.replace(&node_captures[0], &simple_pattern); 
        }
    }
    write!(dotfile, "{}", content).expect("Error while writing into {dotfile}");
    let shell_command = format!("cat {}.dot | dot -Tsvg > {}.svg", pathname, pathname); // Shell command
    Command::new("sh").arg("-c").arg(shell_command).status().expect("Failed to execute command"); 
}

// fn plot_and_log_qrcode(id: &str, character: String) -> () {
//     fs::create_dir_all(OUTPUT_QRCODES_DIR).expect("Failed to create directory.");
//     let qrcode_pathname = OUTPUT_QRCODES_DIR.to_owned() + id + ".svg";
//     let full_weblink = BASE_URL_FOR_QR_CODES.to_owned() + id + ".svg";
//     qrcode_generator::to_svg_to_file_from_str(&full_weblink, QrCodeEcc::Low, 512, None::<&str>,&qrcode_pathname).unwrap();
    
//     // let md_link = format!("{}", character);  
//     // write!(buffer, "{}   \n", md_link).expect("unable to write");
//     // let md_link = format!("* [{}](@/info/slants/{}.md)\n", character, id);  

//     // let md_pathname = format!("../content/info/slants/{}.md", id);
//     // let mut file = File::create(md_pathname).unwrap();
//     // let page_text = format!("+++\ntitle = \"Slant {}\"\ndate = 2023-06-16\n+++\n\n![QR Code](/data/yeast/{}.svg)\n\n[Slant {} Data](/data/yeast/{}.toml)\n\n[All slants](@/info/yeast.md)\n\n ## Propagations\n", id, id, id, id);
//     // file.write_all(page_text.as_bytes()).unwrap();
// }

    // let toml_insides = toml_map.to_string();
    // let toml_pathname = format!("{}{}.toml", OUTPUT_DIR, id);
    // let mut file = File::create(toml_pathname).expect("Could not create sample toml file");
    // file.write_all(toml_insides.as_bytes()).expect("Could not write data to sample toml file");

fn plot_and_log_qrcode_OLD(id: &str, toml_map:Map<String, Value>, buffer: &mut BufWriter<File>) -> () {
    fs::create_dir_all(OUTPUT_DIR.to_owned() + "/qrcodes/").expect("Failed to create directory.");
    let slant_qrcode_pathname = OUTPUT_DIR.to_owned() + "/qrcodes/" + id + ".svg";
    let slant_full_weblink = BASE_URL_FOR_QR_CODES.to_owned() + id;
    qrcode_generator::to_svg_to_file_from_str(&slant_full_weblink, QrCodeEcc::Low, 512, None::<&str>,&slant_qrcode_pathname).unwrap();
    // println!("Created QR code `{}` linking to the `{}`", &slant_qrcode_image_pathname, slant_full_weblink);

    let slant_md_link = format!("* [{}](@/info/slants/{}.md)\n", id, id);  
    write!(buffer, "{}", slant_md_link).expect("unable to write");

    let slant_md_pathname = format!("../content/info/slants/{}.md", id);
    let mut slant_file = File::create(slant_md_pathname).unwrap();
    let slant_page_text = format!("+++\ntitle = \"Slant {}\"\ndate = 2023-06-16\n+++\n\n![QR Code](/data/yeast/{}.svg)\n\n[Slant {} Data](/data/yeast/{}.toml)\n\n[All slants](@/info/yeast.md)\n\n ## Propagations\n", id, id, id, id);
    slant_file.write_all(slant_page_text.as_bytes()).unwrap();

    let slant_toml_insides = toml_map.to_string();
    let slant_toml_pathname = format!("{}{}.toml", OUTPUT_DIR, id);
    let mut file = File::create(slant_toml_pathname).expect("Could not create sample toml file");
    file.write_all(slant_toml_insides.as_bytes()).expect("Could not write data to sample toml file");
}

fn find_organoleptic_node_in_component(component_nodes: &Nodes) -> Option<Map<String, Value>> {
    for (_, toml_map) in component_nodes.iter() {
        if let Some(medium) = try_to_read_field_as_string(&toml_map, "medium") {
            if medium == "organoleptic" {
                return Some(toml_map.clone());
            }
        }
    }
    None
}

fn populate_site_pages(nodes: Nodes) {
    fs::create_dir_all(OUTPUT_QRCODES_DIR).expect("Failed to create directory.");
    fs::create_dir_all("../content/info/yeasts/").expect("Failed to create directory.");
    fs::create_dir_all("../static/yeast-component-output/").expect("Failed to create directory.");

    let yeast_md = OpenOptions::new().write(true).append(true).open(YEAST_PAGE_PATH).expect("Unable to open yeast page.");
    let mut yeast_buffer = BufWriter::new(yeast_md);

    let components = nodes_into_components(nodes.clone(), false);

    //// Ordering the hashmap
    let mut ordered_nodes: Vec<(&str, &Map<String, Value>)> = nodes.iter().map(|(key, value)| (key.as_str(), value)).collect();
    ordered_nodes.sort_by_key(|&(key, _)| key);

    //// Filter for items for yeast page with depth 1:
    let mut depth1_items = HashMap::new();
    for (id, toml_map) in ordered_nodes.clone().into_iter(){
        let maybe_protocol = try_to_read_field_as_map(&toml_map, "protocol");
        if let Some(protocol) = maybe_protocol {
            if protocol.get("name").unwrap().as_str().unwrap() == "package" {
                let ancestors_id = components.iter()
                    .find(|&(_, nodes_map)| nodes_map.contains_key(id))
                    .map(|(component_id, nodes_map)| {component_id.clone()})
                    .unwrap();

                let ancestors_toml_map = nodes.get(ancestors_id.as_str()).unwrap(); 
                let character = try_to_read_field_as_string(ancestors_toml_map, "character").unwrap();
                depth1_items.insert(ancestors_id, character);
            }
        }
    }
    let ordered_depth1_items: BTreeMap<_, _> = depth1_items.into_iter()
        .collect::<Vec<_>>()
        .into_iter()
        .sorted_by(|a, b| a.1.cmp(&b.1))
        .collect();
    
    // let organoleptic_results: 
    //// Populate yeast page with depth 1: 
    for (component_id, character) in ordered_depth1_items {
        let qrcode_pathname = OUTPUT_QRCODES_DIR.to_owned() + &component_id + ".svg";
        let qrcode_weblink = BASE_URL_FOR_QR_CODES.to_owned() + &component_id + ".svg";
        qrcode_generator::to_svg_to_file_from_str(&qrcode_weblink, QrCodeEcc::Low, 512, None::<&str>,&qrcode_pathname).unwrap();
        
        let package_page_link = format!("* [{}](@/info/yeasts/{}.md)\n", character, component_id);
        write!(yeast_buffer, "{}", package_page_link).expect("unable to write");
        
        let pakage_page_pathname = format!("../content/info/yeasts/{}.md", component_id);
        let mut slant_file = File::create(pakage_page_pathname).unwrap();
    
        let maybe_organoleptic_toml = find_organoleptic_node_in_component(components.get(&component_id).unwrap());
    
        let mut character_placeholder = String::from(""); 
        let mut appearence_yeast_placeholder = String::from(""); 
        let mut appearence_liquid_placeholder = String::from("");  

        if let Some(toml_map) = maybe_organoleptic_toml {
            if let Some(character) = try_to_read_field_as_string(&toml_map, "character") {
                character_placeholder = format!("Yeast character: {}   \n", character);
            }
            if let Some(appearance_map) = try_to_read_field_as_map(&toml_map, "appearance") {
                if let Some(liquid) = try_to_read_field_as_string(&appearance_map, "liquid") {
                    appearence_liquid_placeholder = format!("Liquid appearance: {}   \n", liquid);
                }
                if let Some(yeast) = try_to_read_field_as_string(&appearance_map, "yeast") {
                    appearence_yeast_placeholder = format!("Yeast appearance: {}   \n", yeast);
                }
            }
        }
    
        let slant_page_text = format!(
            "+++\n\
            title = \"{}\"\n\
            date = 0001-01-01\n\
            +++\n\
            {}{}{}\n\
            Yeast genealogy:\n\
            ![Genealogy](/yeast-component-output/genealogy/genealogy-{}.svg)\n\n\
            [All yeasts](@/info/yeasts.md)",
            character, character_placeholder, appearence_yeast_placeholder, appearence_liquid_placeholder, component_id
        );
        // [Slant {} Data](/yeast-component-output/{}.toml)\n\n\
        slant_file.write_all(slant_page_text.as_bytes()).unwrap();
        
    }
}


            
            // if let Some(ancestor_id) = maybe_ancestor {
            //     let ancestor_slant_md_pathname = format!("../content/info/slants/{}.md", ancestor_id);

            //     let mut ancestor_slant_file = OpenOptions::new()
            //         .write(true)
            //         .append(true)
            //         .open(ancestor_slant_md_pathname)
            //         .expect("Unable to open slant page.");

            //     let expected_count_plot_pathname = OUTPUT_DIR.to_owned() + &id + "-count.svg";
            //     let expected_density_plot_pathname = OUTPUT_DIR.to_owned() + &id + "-density.svg"; 

            //     if Path::new(&expected_count_plot_pathname).exists() || Path::new(&expected_density_plot_pathname).exists() {
            //         let mut sample_section_text = format!("Sample {} data in graphic format:\n", id); 
            //         if Path::new(&expected_count_plot_pathname).exists() {
            //             sample_section_text += &String::from(format!("![Sample {} count plot](/data/yeast/{}-count.svg)\n", id, id));}
            //         if Path::new(&expected_density_plot_pathname).exists() {
            //             sample_section_text += &String::from(format!("![Sample {} density plot](/data/yeast/{}-density.svg)\n", id, id));
            //         }
            //         ancestor_slant_file.write_all(sample_section_text.as_bytes()).unwrap();
            //     }
            // }
        


fn main() {
    let (nodes, local_file_links) = tomls_into_nodes_and_links(INPUT_DIR);
    let components = nodes_into_components(nodes.clone(), false);
    let problems = components_into_problems(components.clone());
    
    let mut total_cost: f64 = 0.0;
    for (component_id, component_problem) in problems {
        let component_solution =  run_nelder_mead_for_component(component_problem.clone()).unwrap();
        log_nelder_mead_solutions(&component_id, &component_problem, &component_solution);
        total_cost += component_solution.state.cost;
    }
    println!("Total cost for all optimization problems: {:.4}", total_cost);

    for (component_id, component) in components {
        let genealogy_pathname = OUTPUT_GENEALOGY_DIR.to_owned() + "genealogy-" + &component_id;
        plot_genealogy(genealogy_pathname, component, local_file_links.clone());
    }
    let main_genealogy_pathname = OUTPUT_DIR.to_owned() + GENEALOGY_NAME;
    plot_genealogy(main_genealogy_pathname, nodes.clone(), local_file_links.clone());
    println!("Genealogy graph plotted using graphwiz.");

    if Path::new(&YEAST_PAGE_PATH).exists() {
        println!("Yeast page is found at '{}'. Populating it with data.", YEAST_PAGE_PATH);
        populate_site_pages(nodes.clone());
    }
    
    let warnings = WARNINGS.lock().unwrap();
    if !warnings.is_empty() {
        println!("\nWarnings:", ); 
        for w in  warnings.iter() {
            println!("{}", w);
        }
    };
}