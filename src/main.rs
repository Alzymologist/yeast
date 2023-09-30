#![allow(unused_must_use)]

// TODO: get reference as pitch rate from experiment description file

use core::panic;
use std::{ops, fs::{self, OpenOptions,File}};
use std::io::{BufWriter, Write};
use plotters::prelude::*;
use serde::Deserialize;
use toml::{Table, Value};
use time::{format_description::well_known::Iso8601, PrimitiveDateTime};
use std::path::{Path, PathBuf};
use walkdir::WalkDir;
use std::collections::HashMap;
use regex::Regex;
use petgraph::{graphmap::DiGraphMap, dot::{Dot, Config}, visit::{Dfs, Reversed}, prelude::GraphMap};
use toml::map::Map;
use argmin::core::{State, Error, Executor, CostFunction};
use argmin::solver::neldermead::NelderMead;
use ndarray::Array1;
use std::sync::Mutex; 
use qrcode_generator::QrCodeEcc;
use itertools::{join, multiunzip, multizip, Itertools};

lazy_static::lazy_static! { static ref WARNINGS: Mutex<Vec<String>> = Mutex::new(Vec::new()); }
const BASE_URL_FOR_QR_CODES: &str = "https://feature-main.alzymologist-github-io.pages.dev/info/slants/";
const INPUT_DIR: &str = "data/"; 
const OUTPUT_DIR: &str = "output/";
const DOTFILE_NAME: &str = "genealogy";
const YEAST_PAGE_PATH: &str = "../content/info/yeast.md"; 
const NELDER_MEAD_BOUNDS: [(f64, f64); 3] = [(0.0, 1E20f64), (0.0, 1E20f64), (0.0, 10.0)];
const DOT_REPLACEMENT: &str =r##"
digraph {
rankdir=LR
node [style=filled, colorscheme=prgn6]

subgraph cluster_legend {
    label="Legend"
    color=black;
    fontsize=20;
    penwidth=3;
    ranksep=0.2;
    rankdir=TB;
    legend0 [ label = "Missing", style=filled];
    legend1 [ label = "Stock", style=filled, fillcolor="#EE6AA7" ];
    legend2 [ label = "Slant", style=filled, fillcolor=2 ];
    legend3 [ label = "Plate", style=filled, fillcolor=3 ];
    legend4 [ label = "Liquid", style=filled, fillcolor=4 ];
    legend5 [ label = "Organoleptic", style=filled, fillcolor=5 ];
    {legend0 -> legend1 -> legend2 -> legend3 -> legend4 -> legend5 [style=invis];}
}
    "##;

/// All dimensional types
trait Dimension: Copy + Clone {}

/// All dimensional types except unitless
trait PDimension: Dimension {}

#[derive(Debug, PartialEq, Copy, Clone)]
/// Unitless dimension for counts and other unitless values (logarithms, activities, etc)
struct E {}

#[derive(Debug, PartialEq, Copy, Clone)]
/// Density in kg/m3 SI
struct MassDensity {}

#[derive(Debug, PartialEq, Copy, Clone)]
/// Density in counts per m3 in SI
struct UnitDensity {}

#[derive(Debug, PartialEq, Copy, Clone)]
/// Volume in m3 in SI
struct Volume {}

#[derive(Debug, PartialEq, Copy, Clone)]
/// Mass in kg in SI
struct Mass {}

impl Dimension for E {}
impl Dimension for MassDensity {}
impl Dimension for UnitDensity {}
impl Dimension for Volume {}
impl Dimension for Mass {}
impl PDimension for MassDensity {}
impl PDimension for UnitDensity {}
impl PDimension for Volume {}
impl PDimension for Mass {}

#[derive(Debug, Copy, Clone)]
/// Real-world evidence observation data, with uncertainty and units
struct O32<U: Dimension> {
    v: f32,
    e: f32,
    u: U,
}

impl O32<E> {
    pub fn new(value: u32) -> Self {
        O32::<E> {
            v: value as f32,
            e: (value as f32).sqrt(),
            u: E{},
        }
    }

    pub fn exact(value: u32) -> Self {
        O32::<E> {
            v: value as f32,
            e: 0.0,
            u: E{},
        }
    }

    pub fn log2(self) -> Self {
        O32::<E> {
            v: self.v.log2(),
            e: self.e/(self.v*(2f32.ln())),
            u: E{},
        }
    }

    pub fn ln(self) -> Self {
        O32::<E> {
            v: self.v.ln(),
            e: self.e/self.v,
            u: E{},
        }
    }
}

impl <U: PDimension> O32<U> {
    /// Get unitless value of something, by dividing by SI standard. Should be used to get log.
    fn unitless(self) -> O32<E> {
        O32::<E> {
            v: self.v,
            e: self.e,
            u: E{}
        }
    }

    /// Convert value to its unitless log value by normalizing with SI unit.
    pub fn log2(self) -> O32<E> {
        self.unitless().log2()
    }
}

impl O32<Volume> {
    pub fn parse(input: Measurement) -> Result<Self, ReadError> {
        let v = input.value;
        let e = match input.error {
            Some(a) => a,
            None => return Err(ReadError::UncertaintyMissing),
        };
        match input.unit {
            Some(a) => match a.as_str() {
                "ml" => Ok(Self{
                    v: v*1E-6,
                    e: e*1E-6,
                    u: Volume{},
                }),
                _ => Err(ReadError::UnitNotImplemented),
            },
            None => Ok(Self{
                v: v,
                e: e,
                u: Volume{},
            }),
        }
    }
}

impl O32<Mass> {
    pub fn parse(input: Measurement) -> Result<Self, ReadError> {
        let v = input.value;
        let e = match input.error {
            Some(a) => a,
            None => return Err(ReadError::UncertaintyMissing),
        };
        match input.unit {
            Some(a) => match a.as_str() {
                "g" => Ok(Self{
                    v: v*1E-3,
                    e: e*1E-3,
                    u: Mass{},
                }),
                _ => Err(ReadError::UnitNotImplemented),
            },
            None => Ok(Self{
                v: v,
                e: e,
                u: Mass{},
            }),
        }
    }
}

impl O32<MassDensity> {
    pub fn parse(input: Measurement) -> Result<Self, ReadError> {
        let v = input.value;
        let e = match input.error {
            Some(a) => a,
            None => return Err(ReadError::UncertaintyMissing),
        };
        match input.unit {
            Some(a) => match a.as_str() {
                "g/ml" => Ok(Self{
                    v: v*1E3,
                    e: e*1E3,
                    u: MassDensity{},
                }),
                _ => Err(ReadError::UnitNotImplemented),
            },
            None => Ok(Self{
                v: v,
                e: e,
                u: MassDensity{},
            }),
        }
    }

}

impl <U: Dimension> ops::Add for O32<U> {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self {
            v: self.v + other.v,
            e: (self.e.powi(2) + other.e.powi(2)).sqrt(),
            u: self.u,
        }
    }
}

impl <U: Dimension> ops::Sub for O32<U> {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        Self {
            v: self.v - other.v,
            e: (self.e.powi(2) + other.e.powi(2)).sqrt(),
            u: self.u,
        }
    }
}

impl <U: Dimension> ops::Mul<f32> for O32<U> {
    type Output = Self;

    fn mul(self, rhs: f32) -> Self {
        Self {
            v: self.v*rhs,
            e: self.e*rhs,
            u: self.u,
        }
    }
}

fn mul_error(lv: f32, le: f32, rv: f32, re: f32) -> f32 {
    ((le*rv).powi(2) + (re*lv).powi(2)).sqrt()
}

impl <U: Dimension> ops::Mul<O32<E>> for O32<U> {
    type Output = Self;

    fn mul(self, rhs: O32<E>) -> Self {
        Self {
            v: self.v*rhs.v,
            e: mul_error(self.v, self.e, rhs.v, rhs.e),
            u: self.u,
        }
    }
}

impl <U: Dimension> ops::Div<f32> for O32<U> {
    type Output = Self;

    fn div(self, rhs: f32) -> Self {
        Self {
            v: self.v/rhs,
            e: self.e/rhs,
            u: self.u,
        }
    }
}

fn div_error(nv: f32, ne: f32, dv: f32, de: f32) -> f32 {
    ((ne/dv).powi(2) + ((nv*de)/(dv.powi(2))).powi(2)).sqrt()
}

impl <U: PDimension> ops::Div<O32<E>> for O32<U> {
    type Output = Self;

    fn div(self, rhs: O32<E>) -> Self {
        Self {
            v: self.v/rhs.v,
            e: div_error(self.v, self.e, rhs.v, rhs.e),
            u: self.u,
        }
    }
}

impl <U: Dimension> ops::Div<O32<U>> for O32<U> {
    type Output = O32<E>;

    fn div(self, rhs: O32<U>) -> O32<E> {
        O32 {
            v: self.v/rhs.v,
            e: div_error(self.v, self.e, rhs.v, rhs.e),
            u: E{},
        }
    }
}

impl ops::Div<O32<Volume>> for O32<E> {
    type Output = O32<UnitDensity>;

    fn div(self, rhs: O32<Volume>) -> O32<UnitDensity> {
        O32 {
            v: self.v/rhs.v,
            e: div_error(self.v, self.e, rhs.v, rhs.e),
            u: UnitDensity{},
        }
    }
}

impl ops::Div<O32<MassDensity>> for O32<Mass> {
    type Output = O32<Volume>;

    fn div(self, rhs: O32<MassDensity>) -> O32<Volume> {
        O32 {
            v: self.v/rhs.v,
            e: div_error(self.v, self.e, rhs.v, rhs.e),
            u: Volume{},
        }
    }
}

fn mean_weighted<U: Dimension>(data: &[O32<U>]) -> O32<U> {
    let mut numerator = 0f32;
    let mut denominator = 0f32;
    for i in data.iter() {
        numerator += (i.v)/(i.e.powi(2));
        denominator += i.e.powi(-2);
    }
    O32::<U> {
        v: numerator/denominator,
        e: denominator.sqrt().recip(),
        u: data[0].u,
    }
}

fn average<U: Dimension>(data: &[O32<U>]) -> O32<U> {
    mean_weighted(data)
}

#[derive(Debug, Deserialize)]
struct Measurement {
    value: f32,
    error: Option<f32>,
    unit: Option<String>,
}

#[derive(Debug, Deserialize)]
struct UniformityExperiment {
    measurement: Option<Vec<UniformityMeasurement>>,
}

#[derive(Debug, Deserialize)]
struct UniformityMeasurement {
    time: Value,
    density: Option<Measurement>,
    count: Option<Count>,
    dilution: Option<Dilution>,
}

#[derive(Debug, Deserialize)]
struct Dilution {
    sample: Measurement,
    diluent: Measurement,
}

#[derive(Debug, Deserialize)]
struct Count {
    top: Vec<u32>,
    bottom: Vec<u32>,
}

/// Read single measurement file for uniformity test experiment
fn read_uniformity_point(data: UniformityMeasurement) -> Result<UniformityPoint, ReadError> {
    let volume = O32::<Volume> {
        v: 1e-10,
        e: 1e-12,
        u: Volume{},
    };

    let diluent_density = O32::<MassDensity> {
        v: 997.0,
        e: 3.0,
        u: MassDensity{},
    };
    
    let timestamp: PrimitiveDateTime = match data.time {
        Value::Datetime(a) => PrimitiveDateTime::parse(&a.to_string(), &Iso8601::DEFAULT).or(Err(ReadError::TimeFormat(a.to_string())))?,
        Value::String(ref a) => PrimitiveDateTime::parse(&a.to_string(), &Iso8601::DEFAULT).or(Err(ReadError::TimeFormat(a.to_string())))?,
        _ => return Err(ReadError::TimeFormat(data.time.to_string())),
    };
   
    
    let concentration = match data.count {
        Some(count) => {
            let count_top = O32::<E>::new(count.top.iter().sum());
            let count_bottom = O32::<E>::new(count.bottom.iter().sum());

            let count = mean_weighted(&[count_top, count_bottom]);

            if count.v.is_nan() {
                None
            } else {
                match data.dilution {
                    Some(a) => {
                        let sample = O32::<Volume>::parse(a.sample)?;
                        let diluent = O32::<Mass>::parse(a.diluent)?;
                        let dilution = sample/(diluent/diluent_density);
                        Some(count/volume/dilution)
                    },
                    None => Some(count/volume),
                }
            }
        },
        None => None,
    };

    let density = match data.density {
        Some(density) => Some(O32::<MassDensity>::parse(density)?),
        None => None,
    };

    Ok(UniformityPoint { timestamp, concentration, density })
}

#[derive(Debug, PartialEq)]
enum ReadError {
    TimeFormat(String),
    UncertaintyMissing,
    UnitNotImplemented,
    ValueMissing(String),
    ValueType(String),
}

fn plot_count(id:&str, plot_name: &str, nm: &NelderMeadProblemRaw,  cost_per_datapoint: f64, optimized_params: Vec<f64>) ->  Result<(), DrawingAreaErrorKind<std::io::Error>>  {
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
            format!("Cost per datapoint = {:.3e} ", cost_per_datapoint ),
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

fn plot_density(id:&str, plot_name: &str, data: &[UniformityPoint], reference_time: PrimitiveDateTime) -> Result<(), DrawingAreaErrorKind<std::io::Error>> {
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

    ctx
        .draw_series(
            data.iter().filter_map(|point| {
                match point.density {
                    Some(ref density) => {
                        let relative_time_in_hours = (point.timestamp-reference_time).as_seconds_f32() as f64/(60.0*60.0);
                        Some(ErrorBar::new_vertical(relative_time_in_hours, (density.v - density.e) as f64, density.v as f64, (density.v + density.e) as f64, BLUE.filled(), 10))
                    }
                    None => None,
                }
            }
        ))?;
    // println!("Created plot: {:?}", &plot_name);
    Ok(())
}

fn try_to_read_field_as_string (map: &Map<String, Value>, key: &str) -> Option<String>{
    match map.get(key) {
        Some(Value::String(s)) => {Some(s.clone())},
        _ => {None}
    }
}

// We expect TOML field to be a list of strings or just a single string.
// Vec of one string is returned to make handling similar for the parent and participants cases;
fn try_to_read_field_as_vec (map: &Map<String, Value>, key: &str) -> Option<Vec<String>>{
    match map.get(key) {
        Some(Value::String(s)) => Some(vec![s.clone()]),
        Some(Value::Array(arr)) => {
            let strings_vec: Vec<String> = arr.iter()
            .filter_map(|v| v.as_str().map(String::from))
            .collect();
            Some(strings_vec)
        },
        _ => {None}
    }
}

fn try_to_read_reference_time(read_map: &Map<String, Value>) -> Option<PrimitiveDateTime> {
    match read_map.get("time") {
        Some(Value::Datetime(a)) => PrimitiveDateTime::parse(&a.to_string(), &Iso8601::DEFAULT).ok(),
        Some(Value::String(a)) => {
            PrimitiveDateTime::parse(a, &Iso8601::DEFAULT)
                .or_else(|_| PrimitiveDateTime::parse(&format!("{}T00:00", a), &Iso8601::DEFAULT))
                .ok()
        },
        _ => None,
    }
}

fn push_into_warnings(s: &str) -> () {
    let mut warnings = WARNINGS.lock().unwrap();
    warnings.push(s.to_string());
}

// Creates initial simplex for Nelder Mead Problem. Simplex must be non-degenerate. Exact values of parameters and shifts do not matter.
fn create_initial_simplex() -> Vec<Vec<f64>> {
    // const PARAMETER_NAMES: [&str; 3] = ["conc_0", "conc_max", "growth_speed_max"];
    let vertex0: Vec<f64> = vec![1E9f64         , 1E14f64          , 0.5f64];
    let vertex1: Vec<f64> = vec![1E9f64 + 5E9f64, 1E14f64          , 0.5f64];
    let vertex2: Vec<f64> = vec![1E9f64         , 1E14f64 + 5E14f64, 0.5f64];
    let vertex3: Vec<f64> = vec![1E9f64         , 1E14f64          , 0.5f64 + 0.5f64];
    
    vec![vertex0, vertex1, vertex2, vertex3]
}

#[derive(Debug, Clone, Copy)]
/// One cell counting experiment data point
struct UniformityPoint {
    timestamp: PrimitiveDateTime,
    concentration: Option<O32<UnitDensity>>,
    density: Option<O32<MassDensity>>,
}
#[derive(Debug, Clone)]
struct RawNelderMeadPoint {
    conc: O32<UnitDensity>,
    hours: f64,
}
#[derive(Debug, Clone)]
struct NelderMeadProblemRaw {
    points: Vec<RawNelderMeadPoint>,
    reference_time: PrimitiveDateTime,
    bounds: [(f64, f64); 3], 
}

impl CostFunction for NelderMeadProblemRaw {
    type Param = Vec<f64>;
    type Output = f64;
    // params: conc_0, conc_max, growth_speed_max
    fn cost(&self, params: &Self::Param) -> Result<Self::Output, Error> {
        let times = &self.points.clone().into_iter().map(|p|p.hours).collect();
        let modeled_concentrations = calculate_cell_concentrations_using_logistic_model(params, times);

        // Cost using Mean Squared Logarithmic Error (MSLE)
        let cost = modeled_concentrations.into_iter().zip(self.points.clone().into_iter())
        .map(|(modeled_conc, point)| (modeled_conc, point.conc.v as f64, point.conc.e as f64))
        .filter(|&(modeled, actual, error)| modeled > 0.0 && actual > 0.0 && error > 0.0) 
        .map(|(modeled, actual, _error)|(modeled.log10() - actual.log10()).powi(2))
        .sum::<f64>() / self.points.len() as f64;

        // If parameters are out of bounds, the penalty becomes non-zero:
        let mut penalty: f64 = 0.0;
        let penalty_factor = 1E30f64;
        for (param, &(lower_bound, upper_bound)) in params.iter().zip(&self.bounds) {
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

type R = argmin::core::OptimizationResult<NelderMeadProblemRaw, NelderMead<Vec<f64>, f64>, argmin::core::IterState<Vec<f64>, (), (), (), f64>>;

fn run_nelder_mead(nm_problem_raw: NelderMeadProblemRaw) -> Result<R, Error> { 
    let nm_initial_simplex: NelderMead<Vec<f64>, f64> = NelderMead::new(create_initial_simplex());
    let result = Executor::new(nm_problem_raw, nm_initial_simplex)
        .configure(|state| state
            .max_iters(1000)
            .target_cost(0.0)
        )
        .run();

    result
}

fn calculate_cell_concentrations_using_logistic_model(params: &[f64], times: &Vec<f64>) -> Vec<f64> {
    let (conc_0, conc_max, growth_speed_max) = (params[0], params[1], params[2]);
    let concentrations = times
        .iter()
        .map(|time| {
            let exp = (growth_speed_max * time).exp();
            (conc_0 * conc_max * exp) / (conc_max - conc_0 + conc_0 * exp)
        })
        .collect();
    concentrations
}

fn generate_points_with_logistic_model(params: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let number_of_generated_points = 100;
    let min_time = 0f64;
    let max_time = 75f64;
    let created_time_data = Array1::linspace(min_time, max_time, number_of_generated_points).to_vec();
    let created_conc_data = calculate_cell_concentrations_using_logistic_model(&params, &created_time_data);
    (created_time_data, created_conc_data)
 }

fn build_digraphmap(raw_nodes: HashMap<String, Map<String, Value>>) -> DiGraphMap<&'static str, ()> {
    let edges: Vec<(&'static str, &'static str)> = raw_nodes.values().flat_map(|mapping| {
        let id = Box::leak(try_to_read_field_as_string(&mapping, "id").unwrap().into_boxed_str());
        let id_immutable: &'static str = &*id; // Convert to an immutable reference

        let parent_or_participants: Vec<String> = { 
            let parent = try_to_read_field_as_vec(mapping, "parent");
            let participant = try_to_read_field_as_vec(mapping, "participants");
            if let Some(p) = parent {p}
            else if let Some(p) = participant {p}
            else {panic!("Node with id {:?} does not have parent or participants. It should have been checked before.", id);}
        };

        parent_or_participants.into_iter().map(|par_id| {
            let par_id = Box::leak(par_id.into_boxed_str());
            let par_id_immutable: &'static str = &*par_id; 
            (par_id_immutable, id_immutable)
        })
        .filter(|(par_id, _)| *par_id != "new" ) // Filter out edges with "new" parent or participant
        .collect::<Vec<_>>()
        
    }).collect();

    let mut graph = DiGraphMap::<&str, ()>::from_edges(edges);
    for mapping in raw_nodes.values() { // Add graph components with no edges (lone nodes)
        let id = Box::leak(try_to_read_field_as_string(&mapping, "id").unwrap().into_boxed_str());
        let id_immutable: &'static str = &*id; // Convert to an immutable reference
        if !graph.contains_node(id_immutable) {
            graph.add_node(id_immutable);
        }
    }
    graph
}

fn prepare_dot_file(graph: &GraphMap<&str, (), petgraph::Directed>, checked_nodes: HashMap<String, Map<String, Value>>) {
    //// Creating initial .dot string from the graph
    let dot = Dot::with_config(&graph,&[Config::EdgeNoLabel]);
    let mut content = format!("{:?}", dot);
    let dotfile_path = OUTPUT_DIR.to_owned() + DOTFILE_NAME + ".dot";
    let mut dotfile = File::create(&dotfile_path).unwrap();
    
    let colours_for_marked_dirs: HashMap<&str, &str> = HashMap::from([
        ("slant", "2"),
        ("plate", "3"),
        ("liquid", "4"),
        ("organoleptic", "5"),
        ("stock", r##" "#EE6AA7" "##),
        ]);

    //// Editing of .dot after generation 
    content = content.replacen("digraph {", DOT_REPLACEMENT, 1);

    //// Searching for lines describing nodes in .dot string
    let node_pattern = Regex::new(r#"(\d+) \[ label = "\\"(.+?)\\"" \]"#).unwrap();
    for node_captures in node_pattern.captures_iter(&content.clone()) {
        let captured_nodenumber = node_captures.get(1).unwrap().as_str();
        let captured_label = node_captures.get(2).unwrap().as_str();

        if let Some(tomlmap) = checked_nodes.get(captured_label)  { 
            let medium = tomlmap.get("medium").unwrap().as_str().unwrap(); 
            if let Some(colour) = colours_for_marked_dirs.get(medium) {  //// If there is corresponding colour for this medium
                let pattern_with_colour = format!(r#"{} [ label = "{}" fillcolor={} ]"#, captured_nodenumber, captured_label, colour);
                content = content.replace(&node_captures[0], &pattern_with_colour);
            }
        } else { //// If Captured_label was not found in the graph
                push_into_warnings( &format!("{:<25} {:<25}", captured_label, "id does not correspond to any converted TOML") );
                let simple_pattern = format!(r#"{} [ label = "{}" ]"#, captured_nodenumber, captured_label);
                content = content.replace(&node_captures[0], &simple_pattern); 
        }
    }
    write!(dotfile, "{}", content).expect("Error while writing into {dotfile}");
    println!("\nCreated {:?} for Graphviz. You can create image with this shell command for example:\ncat {} | dot -Tpng > {}\n", dotfile_path, dotfile_path, OUTPUT_DIR.to_owned() + DOTFILE_NAME + ".png");
}

fn populate_site_pages(graph: &GraphMap<&str, (), petgraph::Directed>, checked_nodes: HashMap<String, Map<String, Value>>) {
    let yeast_md = OpenOptions::new()
    .write(true)
    .append(true)
    .open(YEAST_PAGE_PATH)
    .expect("Unable to open yeast page.");

    fs::create_dir_all("../content/info/slants/").expect("Failed to create directory.");
    fs::create_dir_all("../static/data/yeast/").expect("Failed to create directory.");

    let mut yeast_buffer = BufWriter::new(yeast_md);

    //// Ordering the hashmap:
    let mut ordered_nodes: Vec<(&str, &Map<String, Value>)> = checked_nodes.iter().map(|(key, value)| (key.as_str(), value)).collect();
    ordered_nodes.sort_by_key(|&(key, _)| key);

    //// First pass over nodes (to write slant data).
    for (id, tomlmap) in ordered_nodes.iter() {
        if tomlmap.get("medium").unwrap().as_str().unwrap() == "slant" {
            let slant_full_weblink = BASE_URL_FOR_QR_CODES.to_owned() + id;
            let slant_qrcode_image_pathname = OUTPUT_DIR.to_owned() + id + ".svg";
            qrcode_generator::to_svg_to_file_from_str(&slant_full_weblink, QrCodeEcc::Low, 512, None::<&str>,&slant_qrcode_image_pathname).unwrap();
            println!("Created QR code `{}` linking to the `{}`", &slant_qrcode_image_pathname, slant_full_weblink);

            let slant_md_link = format!("* [{}](@/info/slants/{}.md)\n", id, id);  
            write!(yeast_buffer, "{}", slant_md_link).expect("unable to write");
            let slant_md_pathname = format!("../content/info/slants/{}.md", id);
            let mut slant_file = File::create(slant_md_pathname).unwrap();
            let slant_page_text = format!("+++\ntitle = \"Slant {}\"\ndate = 2023-06-16\n+++\n\n![QR Code](/data/yeast/{}.svg)\n\n[Slant {} Data](/data/yeast/{}.toml)\n\n[All slants](@/info/yeast.md)\n\n ## Propagations\n", id, id, id, id);
            slant_file.write_all(slant_page_text.as_bytes()).unwrap();

            let slant_toml_insides = tomlmap.to_string();
            let slant_toml_pathname = format!("{}{}.toml", OUTPUT_DIR, id);
            let mut file = File::create(slant_toml_pathname).expect("Could not create sample toml file");
            file.write_all(slant_toml_insides.as_bytes()).expect("Could not write data to sample toml file");
        }
    }
    //// Second pass over nodes (to write other data).
    for (id, tomlmap) in ordered_nodes.iter() {
        if tomlmap.get("medium").unwrap().as_str().unwrap() != "slant" { 

            //// Deep first search for the ancestor slant on the reversed graph:
            let mut maybe_ancestor = None; 
            let reversed_graph = Reversed(&graph);
            let mut dfs = Dfs::new(&reversed_graph, id);
 'dfs_loop: while let Some(nx) = dfs.next(&reversed_graph) {
                if let Some(tomlmap) = checked_nodes.get(nx) {
                    let medium = tomlmap.get("medium").unwrap().as_str().unwrap();
                    if medium == "slant" {
                        maybe_ancestor = Some(nx); 
                        break 'dfs_loop; 
                    }
                }
            }
            if let Some(ancestor_id) = maybe_ancestor {
                let ancestor_slant_md_pathname = format!("../content/info/slants/{}.md", ancestor_id);

                let mut ancestor_slant_file = OpenOptions::new()
                .write(true)
                .append(true)
                .open(ancestor_slant_md_pathname)
                .expect("Unable to open slant page.");

                let expected_count_plot_pathname = OUTPUT_DIR.to_owned() + &id + "-count.svg";
                let expected_density_plot_pathname = OUTPUT_DIR.to_owned() + &id + "-density.svg"; 

                if Path::new(&expected_count_plot_pathname).exists() || Path::new(&expected_density_plot_pathname).exists() {
                    let mut sample_section_text = format!("Sample {} data in graphic format:\n", id); 
                    if Path::new(&expected_count_plot_pathname).exists() {
                        sample_section_text += &String::from(format!("![Sample {} count plot](/data/yeast/{}-count.svg)\n", id, id));}
                    if Path::new(&expected_density_plot_pathname).exists() {
                        sample_section_text += &String::from(format!("![Sample {} density plot](/data/yeast/{}-density.svg)\n", id, id));
                    }
                    ancestor_slant_file.write_all(sample_section_text.as_bytes()).unwrap();
                }
            }
        }
     }
}

fn tomls_into_checked_nodes (input_dir: &str) -> HashMap<String, Map<String, Value>> {
    let toml_files: Vec<PathBuf> = WalkDir::new(input_dir)
    .into_iter()
    .filter_map(|entry| entry.ok())
    .filter(|entry| 
        entry.file_type().is_file() &&
        entry.file_name().to_string_lossy().ends_with(".toml"))
    .map(|entry| entry.into_path() )
    .collect();
    println!("TOML files found: {:?}", toml_files.len());

    let mut checked_nodes: HashMap<String, Map<String, Value>> = HashMap::new();
    for file in toml_files {
        let contents = fs::read_to_string(&file.as_path()).unwrap();
        let mut toml_map: Map<String, Value> = contents.parse::<Table>().unwrap();
        let id = try_to_read_field_as_string(&toml_map, "id");
        let time = try_to_read_reference_time(&toml_map);
        let medium = try_to_read_field_as_string(&toml_map, "medium");
        let parent = try_to_read_field_as_vec(&toml_map, "parent");
        let participants = try_to_read_field_as_vec(&toml_map, "participants");
        match (&id, medium, time, parent, participants) {
            (None, _, _, _, _) => { push_into_warnings( &format!("{:<25} {:<25} {:?}", "Id (at least)", "missing from", file) );},
            (_, None, _, _, _) => { push_into_warnings( &format!("{:<25} {:<25} {:?}", "Medium (at least)", "missing from", file) );}
            (_, _, None, _, _) => { push_into_warnings( &format!("{:<25} {:<25} {:?}", "Time (at least)", "missing from", file) );}
            (_, _, _, Some(_), Some(_)) => { push_into_warnings( &format!("{:<25} {:<25} {:?}", "Parent and participants", "can't be both present", file))}
            (Some(_), Some(_), Some(_), Some(_), None)  => {checked_nodes.insert(id.unwrap(), toml_map);},
            (Some(_), Some(_), Some(_), None, Some(_))  => {checked_nodes.insert(id.unwrap(), toml_map);},
            (Some(_), Some(m), Some(_), None, None)  => {
                if m == "stock" {
                    let new = Value::Array(vec![Value::String("new".into())]); 
                    toml_map.insert(String::from("parent"),new); // Modifing toml to mark impilictly that sample is new
                    checked_nodes.insert(id.unwrap(), toml_map); 
                } else { push_into_warnings( &format!("{:<25} {:<25} {:?}", "Parent or participants", "missing from", file));}
            },
        }
    } 
    println!("TOML files converted into checked nodes: {:?}", checked_nodes.len());
    checked_nodes
}

fn checked_nodes_into_problems(nodes: HashMap<String, Map<String, Value>>) -> HashMap<String, NelderMeadProblemRaw> {
    let mut nm_problems: HashMap<String, NelderMeadProblemRaw> = HashMap::new();

    for (id, toml_map) in nodes.into_iter() {
        let experiment_res = toml_map.clone().try_into::<UniformityExperiment>();
        let reference_time = try_to_read_reference_time(&toml_map).unwrap();

        if let Some(measurements) = experiment_res.unwrap().measurement {  // TOML can contain no measurements
            let points: Vec<RawNelderMeadPoint> = measurements
            .into_iter()
            .filter_map(|measurement| read_uniformity_point(measurement).ok())
            .filter(|point| point.concentration.is_some())
            .map (|point| RawNelderMeadPoint {
                conc: point.concentration.unwrap(),
                hours: (point.timestamp - reference_time).as_seconds_f64()/(60.0*60.0)
            })
            .collect();

            if points.len() > 2 {
                let nm_problem_raw = NelderMeadProblemRaw {
                    points,
                    reference_time,
                    bounds: NELDER_MEAD_BOUNDS,
                }; 
                nm_problems.insert(id, nm_problem_raw);
            }
        }
    }
    println!("Checked nodes converted into optimization problems: {:?}", nm_problems.len());
    nm_problems
} 

fn output_nm_results(id: &String, nm: &NelderMeadProblemRaw, nm_result: &R, cost_per_datapoint: f64) -> () {
    fs::create_dir_all(OUTPUT_DIR).expect("Failed to create directory.");
    if let Ok(file) = File::create(format!("{}/{}-count-fitting.txt", OUTPUT_DIR, id)) {
        let mut buffer = BufWriter::new(file);
        let checked_cell_conc_for_printing = join(nm.points.iter().map(|point| format!("{:.3e}", point.conc.v)),", ");
        let checked_times_for_printing = join(nm.points.iter().map(|point| format!("{:.5}", point.hours)),", ");
        let fitting_results = format!("{}", &nm_result);
        
        writeln!(buffer, "Fitting results for id {:?}", id).unwrap();
        writeln!(buffer, "Used cell concentrations:   [{}]", checked_cell_conc_for_printing).unwrap();
        writeln!(buffer, "Used times (hours since reference): [{}]", checked_times_for_printing,).unwrap();
        writeln!(buffer, "Reference time: {:?}\n", &nm.reference_time).unwrap();
        writeln!(buffer, "{}", fitting_results.trim()).unwrap();
        writeln!(buffer, "Parameters for Nelder-Mead problem are:\nCell concentration at the reference time (cells/m^3),\nMaximal cell concentration (cells/m^3),\nMaximal specific cell growth rate (1/h).").unwrap();
        buffer.flush().unwrap();
    }

    let plot_name_count: String = OUTPUT_DIR.to_owned() + id + "-count.svg";
    let plot_name_density: String = OUTPUT_DIR.to_owned() + id + "-density.svg";
    let optimized_params = nm_result.state().get_best_param().unwrap().clone();
    plot_count(id, &plot_name_count, &nm.clone(), cost_per_datapoint, optimized_params);
    // plot_density(&id, &plot_name_density, nm.checked_points, nm.reference_time);
}

fn main() {
    let mut total_cost: f64 = 0.0;
    let checked_nodes = tomls_into_checked_nodes(INPUT_DIR);
    let problems = checked_nodes_into_problems(checked_nodes.clone());

    for (id, problem) in problems {
        let nm_result =  run_nelder_mead(problem.clone()).unwrap();
        output_nm_results(&id, &problem, &nm_result, nm_result.state.cost );
        total_cost += nm_result.state.cost ;
    } 

    println!("Sum of costs for all datasets: {:.4}", total_cost);
    let graph = build_digraphmap(checked_nodes.clone());
    prepare_dot_file(&graph, checked_nodes.clone());

    if Path::new(&YEAST_PAGE_PATH).exists() {
        println!("Yeast page is found at '{}'. Populating it with data.", YEAST_PAGE_PATH);
        populate_site_pages(&graph, checked_nodes.clone());
    } else {
        println!("Yeast page is missing at '{}'. Can't populate it with data.", YEAST_PAGE_PATH);
    };
    
    let warnings = WARNINGS.lock().unwrap();
    if !warnings.is_empty() {
        println!("\nWarnings:", ); 
        for w in  warnings.iter() {
            println!("{}", w);
        }
    };
}