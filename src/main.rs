
use std::{marker::PhantomData, ops, fs::{self, OpenOptions,File, read}};
use std::ffi::OsString;
use std::io::{BufWriter, Write};
use nalgebra::{DVector, DMatrix, linalg::SVD};
use plotters::prelude::*;
use serde::{Deserialize, Serialize};
use toml::{Table, Value};
use time::{format_description::well_known::Iso8601, PrimitiveDateTime};
use std::path::Path;
use std::collections::HashMap;
use regex::Regex;
use petgraph::{graphmap::DiGraphMap, dot::{Dot, Config}, visit::{Dfs, Reversed}};
use toml::map::Map;

const OUTPUT_DIR: &str = "output/";

/// All dimensional types
trait Dimension: Copy + Clone {}

/// All dimensional types except unitless
trait PDimension: Dimension {
/*    fn log2(&self) -> LogOf<Self> {
        O32::<LogOf<Self>> {
            v: self.v.log2(),
            e: self.e/self.v,
            u: LogOf::<Self> {_p: PhantomData::<Self>},
        }
    }*/
}

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

/// Struct to fit data with line; determinictic up to computational error
///
/// TODO: add observation typesafety and error propagtion
#[derive(Debug, Deserialize)]
struct LinearFit {
    pub tan: f32,
    pub intercept: f32,
}

impl LinearFit {
    /// free_coefficient coefficient to multiply by free member; should be 1.0, but might be causing computational
    /// error if x data magnitude is too far. Adjust to get sane results; think of it as initial
    /// fit parameter.
    pub fn solve(x: Vec<f32>, y: Vec<f32>, w: Vec<f32>, free_coefficient: f32) -> Option<Self> {
        let length = x.len();
        if length < 2 {return None};
        let x_vector = DVector::<f32>::from_vec(x);
        let ones = DVector::<f32>::repeat(length, free_coefficient);
        let indep_t = DMatrix::from_columns(&[x_vector, ones]).transpose();

        let w = DVector::<f32>::from_vec(w);

        let w = DMatrix::from_diagonal(&w);

        let y_vector = DVector::from_vec(y);
        let dep = DMatrix::from_columns(&[y_vector]);

        let left_m = indep_t.clone()*w.clone()*indep_t.clone().transpose();

        let left = SVD::new(left_m, true, true);
                    
        let right = indep_t * w * dep;
        let fit = left.solve(&right, 1e-32);
                    
        if let Ok(fit) = fit {
            if fit[0] != 0f32 || fit[1] != 0f32 {
                return Some(Self{
                    tan: fit[0],
                    intercept: fit[1]*free_coefficient,
                });
            }
        }
        return None;
    }
}

#[derive(Debug, Deserialize)]
struct Measurement {
    value: f32,
    error: Option<f32>,
    unit: Option<String>,
}

#[derive(Debug)]
/// One cell counting experiment data point
struct UniformityPoint {
    timestamp: PrimitiveDateTime,
    concentration: Option<O32<UnitDensity>>,
    density: Option<O32<MassDensity>>,
}

#[derive(Debug, Deserialize)]
struct UniformityExperiment {
    id: Option<String>,
    time: Option<Value>,
    measurement: Option<Vec<UniformityMeasurement>>,
    parent: Option<Value>,
}

// #[derive(Debug, Deserialize, Serialize)]
// #[serde(untagged)]
// enum ParentEnum {
//     VecVar(Vec<String>),
//     NoneVar,
// }

// impl ParentEnum {
//     fn is_some(&self) -> bool {
//         match self {
//             ParentEnum::NoneVar => false,
//             _ => true,
//         }
//     }
// }

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

    Ok(UniformityPoint {
        timestamp: timestamp,
        concentration: concentration,
        density: density,
    })
}

#[derive(Debug, PartialEq)]
enum ReadError {
    TimeFormat(String),
    UncertaintyMissing,
    UnitNotImplemented,
    ValueMissing(String),
    ValueType(String),
}

fn plot_counts(name: &str, data: &[UniformityPoint], fit: Option<LinearFit>, reference_time: PrimitiveDateTime) -> Result<(), DrawingAreaErrorKind<std::io::Error>> {
    let output_name = OUTPUT_DIR.to_owned() + name + "-count.svg";
    // println!("writing {output_name}");
    let root_drawing_area = SVGBackend::new(&output_name, (1024, 768)).into_drawing_area();
    root_drawing_area.fill(&WHITE).unwrap();

    let x_min = 0f32;
    let x_max = 800000f32;
    let y_min = 0f32;//1E10f32;
    let y_max = 5E14f32;

    let mut ctx = ChartBuilder::on(&root_drawing_area)
        .set_label_area_size(LabelAreaPosition::Left, 80)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption(name, ("sans-serif", 40))
        .build_cartesian_2d(x_min..x_max, (y_min..y_max).log_scale())?;

    ctx.configure_mesh()
        .x_desc("relative time, s")
        .axis_desc_style(("sans-serif", 40))
        .x_label_formatter(&|x| format!("{:e}", x))
        .x_label_style(("sans-serif", 20))
        .y_desc("cell concentration, CFU/m^3")
        .y_label_formatter(&|x| format!("{:e}", x))
        .y_label_style(("sans-serif", 20))
        .draw()?;

    ctx
        .draw_series(
            data.iter().filter_map(|x| {
                match x.concentration {
                    Some(ref concentration) => 
                        Some(ErrorBar::new_vertical((x.timestamp-reference_time).as_seconds_f32(), concentration.v - concentration.e, concentration.v, concentration.v+concentration.e, BLUE.filled(), 10)),
                    None => None,
                }
            }
        ))?;
    if let Some(fit) = fit {
        let doubling_rate = fit.tan.recip()/3600f32;
        if doubling_rate>0f32 {
            // println!("Doubling rate {doubling_rate} h");
            ctx.draw_series(LineSeries::new(
                    ((x_min as u32)..(x_max as u32))
                        .step_by(1000)
                        .map(|x| ((x) as f32, (fit.tan*(x as f32)+fit.intercept)
                                .exp2()
                                .clamp(1f32, f32::MAX))), RED))?
                .label(format!("Doubling rate {doubling_rate} h"))
                .legend(|(x, y)| PathElement::new(vec![(x, y), (x+20, y)], RED));
            ctx.configure_series_labels()
                .border_style(&BLACK)
                .background_style(&WHITE.mix(0.8))
                .draw()?;
        }
    }
    println!("{}", name);
    let mut seconds: Vec<i64> = vec!();
    let mut concentrations: Vec<f32> = vec!();
    for point in data.iter(){
        let s = point.timestamp.assume_utc().unix_timestamp();
        if let Some(c) = point.concentration
        {
            concentrations.push(c.v);
        };
        seconds.push(s);
    }
    println!("{:?}\n{:?}\n", seconds, concentrations); 

    Ok(())
}

fn plot_density(name: &str, data: &[UniformityPoint], reference_time: PrimitiveDateTime) -> Result<(), DrawingAreaErrorKind<std::io::Error>> {
    let output_name = OUTPUT_DIR.to_owned() + name + "-density.svg";
    let root_drawing_area = SVGBackend::new(&output_name, (1024, 768)).into_drawing_area();
    root_drawing_area.fill(&WHITE).unwrap();

    let mut ctx = ChartBuilder::on(&root_drawing_area)
        .set_label_area_size(LabelAreaPosition::Left, 100)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .caption(name, ("sans-serif", 40))
        .build_cartesian_2d(0f32..800000f32, 1000f32..1100f32)?;

    ctx.configure_mesh()
        .x_desc("relative time, s")
        .axis_desc_style(("sans-serif", 40))
        .x_label_formatter(&|x| format!("{:e}", x))
        .x_label_style(("sans-serif", 20))
        .y_desc("density, kg/m^3")
        .y_label_style(("sans-serif", 20))
        .draw()?;

    ctx
        .draw_series(
            data.iter().filter_map(|x| {
                match x.density {
                    Some(ref density) => 
                        Some(ErrorBar::new_vertical((x.timestamp-reference_time).as_seconds_f32(), density.v - density.e, density.v, density.v+density.e, BLUE.filled(), 10)),
                    None => None,
                }
            }
        ))?;
    Ok(())
}

fn try_to_read_field (map: &Map<String, Value>, key: &str, filename_for_printing: Option<OsString>) -> Option<String>{
    match map.get(key) {
        Some(Value::String(s)) => {
            Some(s.clone()) 
            },
        _ => {
            if let Some(f) = filename_for_printing {
            println!("{:<15} is missing in {:?}", key, f);
        }
            None
        }
    }
}

fn try_to_read_parent (map: &Map<String, Value>, filename_for_printing: Option<OsString>) -> Option<Vec<String>>{
    match map.get("parent") {
        Some(Value::String(s)) => Some(vec![s.clone()]),
        Some(Value::Array(arr)) => {
            let strings_vec: Vec<String> = arr.iter()
            .filter_map(|v| v.as_str().map(String::from))
            .collect();
            Some(strings_vec)
        },
        _ => {
            if let Some(f) = filename_for_printing {
                println!("{:<15} is missing in {:?}", "parent", f);
            }
            None
        }
    }
}

fn main() {
    //let reference_time = PrimitiveDateTime::parse("2023-05-22T00:00", &Iso8601::DEFAULT).expect("Time was parsed by toml, stricter than ISO");
    // TODO: get reference as pitch rate from experiment description file

    let data_dirs = Vec::from([
       String::from("data/slants/2023/"), 
       String::from("data/liquid/2023/"),
       String::from("data/plates/2023/"), 
       String::from("data/stock/2023/"),
       String::from("data/organoleptic/2023"),
    ]);

    fs::create_dir_all(OUTPUT_DIR).expect("Failed to create directory.");

    let mut raw_nodes: HashMap<String, Map<String, Value>> = HashMap::new();
    for data_dir in data_dirs {
        for file in fs::read_dir(&data_dir).unwrap() {
            let mut points = Vec::new();

        //// Check if TOML data is okay:
        if let Ok(file) = file {
            if file.file_name().into_string().unwrap().split(".").last() == Some("toml") && file.file_type().expect("file should have type in this platform").is_file() {
                let contents = fs::read_to_string(&file.path()).unwrap();
                let read_map = contents.parse::<Table>().unwrap();
                
                // let parent = try_to_read_field(&read_map, "parent", file.file_name());
                let id = try_to_read_field(&read_map, "id", Some(file.file_name()));
                let medium = try_to_read_field(&read_map, "medium", Some(file.file_name()));
                let parent = try_to_read_parent(&read_map, Some(file.file_name())); 

                let time = match read_map.get("time") {
                    Some(x) => match x {
                        Value::Datetime(a) => Some(PrimitiveDateTime::parse(&a.to_string(), &Iso8601::DEFAULT).or(Err(ReadError::TimeFormat(a.to_string()))).unwrap()),
                        Value::String(ref a) => {
                            match PrimitiveDateTime::parse(&a.to_string(), &Iso8601::DEFAULT) {
                                Ok(a) => Some(a),
                                Err(_) => Some(PrimitiveDateTime::parse(&(a.to_string() + "T00:00"), &Iso8601::DEFAULT).or(Err(ReadError::TimeFormat(a.to_string()))).unwrap()),
                            }
                        },
                        _ => { 
                            println!("{:<15} is incorrect in {:?}", "time", file.file_name());
                            None
                        }
                    },
                    _ => {
                        println!("{:<15} is missing in {:?}", "time", file.file_name());
                        None
                    },
                };
                
                //// Section where existence of id, parent and medium is guaranteed:
                if id.is_some() && parent.is_some() && medium.is_some() {
                    raw_nodes.insert(id.clone().unwrap(), read_map.clone());
                    
                    let data: UniformityExperiment = toml::from_str(&contents).unwrap();
                    if let Some(measurements) = data.measurement {
                        for measurement in measurements {
                            points.push(read_uniformity_point(measurement).unwrap());
                        }
                    }

                    //// Section where existence of correct time is additionally guaranteed: 
                    if let Some(time_guaranteed) = time {
                        let (x, (y, w)): (Vec<f32>, (Vec<f32>, Vec<f32>)) = points.iter().filter_map(|x| {
                            match x.concentration {
                                Some(a) => {
                                    let log_a = a.log2();
                                    let timestamp_diff: f32 =  (x.timestamp - time.unwrap()).as_seconds_f32();
                                    Some((timestamp_diff, (log_a.v, log_a.e.powi(-2))))
                                },
                                None => None,
                            }
                        }).unzip();
                                    
                        let fit = LinearFit::solve(x, y, w, 1e5);
                        plot_counts(&id.clone().unwrap(), &points, fit, time_guaranteed);
                        plot_density(&id.clone().unwrap(), &points, time_guaranteed);
                    }
                }
            }
        } 
    }
}
    //// Creating initial .dot string from graph edges
    let edges: Vec<(&'static str, &'static str)> = raw_nodes.values().flat_map(|mapping| {

        let id = Box::leak(try_to_read_field(&mapping, "id", None).unwrap().into_boxed_str());
        let id_immutable: &'static str = &*id; // Convert to immutable reference
        let parent = try_to_read_parent(&mapping, None).unwrap(); 

        parent.into_iter().map(|parent_id| {
            let parent_id = Box::leak(parent_id.into_boxed_str());
            let parent_id_immutable: &'static str = &*parent_id; 
            (parent_id_immutable, id_immutable)
        }).collect::<Vec<_>>()
    }).collect();

    let graph = DiGraphMap::<_, ()>::from_edges(edges);
    let dot = Dot::with_config(&graph,&[Config::EdgeNoLabel]);
    let dotfile = "genealogy.dot";
    let mut file = File::create(dotfile).unwrap();
    println!("\nCreating {:?}", dotfile);
    let mut content = format!("{:?}", dot);
    ////
    
    //// Colors for media
    let colours_for_media: HashMap<&str, &str> = HashMap::from([
        ("slant", "2"),
        ("plate", "3"),
        ("liquid", "4"),
        ("organoleptic", "5"),
        ]);


    //// Editing of .dot after generation 
let replacement_str =r#"
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
        legend1 [ label = "Slant", style=filled, fillcolor=2 ];
        legend2 [ label = "Plate", style=filled, fillcolor=3 ];
        legend3 [ label = "Liquid", style=filled, fillcolor=4 ];
        legend4 [ label = "Organoleptic", style=filled, fillcolor=5 ];
        {legend1 -> legend2 -> legend3 -> legend4 [style=invis];}
    }
        "#;

    content = content.replacen("digraph {", replacement_str, 1); //// Configuration for vertical layout and colouring

    let re = Regex::new(r#"\[ label = "\\"(.+?)\\"" \]"#).unwrap();

    for captures in re.captures_iter(&content.clone()) {
        let scanned_id = captures.get(1).unwrap().as_str();

        if let Some(tomlmap) = raw_nodes.get(scanned_id) {
            let medium = tomlmap.get("medium").unwrap().as_str().unwrap(); 
            if let Some(colour) = colours_for_media.get(medium) {
                let complex_pattern = format!(r#"[ label = "{}" fillcolor={} ]"#, scanned_id, colour);
                content = content.replace(&captures[0], &complex_pattern);
            } 
        } else {
            print!("Id '{}' was not found in node keys. \n", scanned_id);
            let simple_pattern = format!(r#"[ label = "{}" ]"#, scanned_id);
            content = content.replace(&captures[0], &simple_pattern); 
        }
    }
    write!(file, "{}", content).expect("Error while writing into {dotfile}");
    // Use this shell command to create image from .dot file:
    // cat genealogy.dot | dot -Tpng > genealogy.png 
    ////
    
    //// Will we populate site pages?
    let yeast_page_path = "../content/info/yeast.md"; 
    let populating_site_pages = match Path::new(&yeast_page_path).exists() {
        true => { 
            println!("File {} is found.", yeast_page_path);
            true
        },
        false => {
            println!("File {} does not exist.", yeast_page_path);
            false
        },
    };
    
    if populating_site_pages {
        let yeast_md = OpenOptions::new()
        .write(true)
        .append(true)
        .open(yeast_page_path)
        .expect("Unable to open yeast page.");

        fs::create_dir_all("../content/info/slants/").expect("Failed to create directory.");
        fs::create_dir_all("../static/data/yeast/").expect("Failed to create directory.");

        let mut yeast_buffer = BufWriter::new(yeast_md);

        //// Ordering the Hashmap:
        let mut ordered_nodes: Vec<(&str, &Map<String, Value>)> = raw_nodes.iter().map(|(key, value)| (key.as_str(), value)).collect();
        ordered_nodes.sort_by_key(|&(key, _)| key);

        //// First pass over nodes. Used to write slant data.
        for (id, tomlmap) in ordered_nodes.iter() {
            if tomlmap.get("medium").unwrap().as_str().unwrap() == "slant" {
                let slant_link = format!("* [{}](@/info/slants/{}.md)\n", id, id);  
                write!(yeast_buffer, "{}", slant_link).expect("unable to write");
                let slant_filename = format!("../content/info/slants/{}.md", id);
                let mut slant_file = File::create(slant_filename).unwrap();
                let slant_page_text = format!("+++\ntitle = \"Slant {}\"\ndate = 2023-06-16\n+++\n\n[Slant {} Data](/data/yeast/{}.toml)\n\n[All slants](@/info/yeast.md)\n\nPropagations:\n", id, id, id);
                slant_file.write_all(slant_page_text.as_bytes()).unwrap();

                let slant_toml_string = tomlmap.to_string();
                let slant_toml_filname = format!("../static/data/yeast/{}.toml", id);
                let mut file = File::create(slant_toml_filname).expect("Could not create sample toml file");
                file.write_all(slant_toml_string.as_bytes()).expect("Could not write data to sample toml file");
            }
        }
        //// Second pass over nodes. Used to write other data.
        for (id, tomlmap) in ordered_nodes.iter() {
            if tomlmap.get("medium").unwrap().as_str().unwrap() != "slant" { 

                //// Deep first search of ancestor slant on the reversed graph:
                let mut ancestor_option = None; 
                let reversed_graph = Reversed(&graph);
                let mut dfs = Dfs::new(&reversed_graph, id);
    'dfs_loop: while let Some(nx) = dfs.next(&reversed_graph) {
                    if let Some(tomlmap) = raw_nodes.get(nx) {
                        let medium = tomlmap.get("medium").unwrap().as_str().unwrap();
                        if medium == "slant" {
                            // println!("Ancestor of {} is {}", id, nx);
                            ancestor_option = Some(nx); 
                            break 'dfs_loop; 
                        }
                    }
                }
                if let Some(ancestor_id) = ancestor_option {
                    let slant_page_filename = format!("../content/info/slants/{}.md", ancestor_id);

                    let mut slant_file = OpenOptions::new()
                    .write(true)
                    .append(true)
                    .open(slant_page_filename)
                    .expect("Unable to open slant page.");

                    let slant_page_text = format!("* [Sample {} Data](/data/yeast/{}.toml)\n", id, id);
                    slant_file.write_all(slant_page_text.as_bytes()).unwrap();
                    
                    let sample_toml_string = tomlmap.to_string();

                    //// Writing TOML with sample data:
                    let sample_toml_filname = format!("../static/data/yeast/{}.toml", id);
                    let mut file = File::create(sample_toml_filname).expect("Could not create sample toml file");
                    file.write_all(sample_toml_string.as_bytes()).expect("Could not write data to sample toml file");
                }
            }
        }
    
}

}
