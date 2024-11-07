use semver::Version;
use time::{format_description::well_known::Iso8601, PrimitiveDateTime};
use toml::{Table, Value};

use super::keywords::{
    BATCH, BOTTOM, CLUSTERING, COUNT, DEG_C, DENSITY, DILUENT, DILUTION, ERROR, G, ID, LIQUID, LOW,
    MEASUREMENT, METHOD, ML, MODERATE, NAME, PACKAGE, PARENT, PITCH, PROPAGATION, PROTOCOL, QA,
    REFRACTOMETRY_METHOD, REVIVE, SAMPLE, STOCK, STRONG, STYLE, TEMPERATURE, TIME, TOP, UNIFORMITY,
    UNIT, VALUE, VERSION,
};
use super::utils::{try_table_and_id_from_file, Extractable};
use super::valid_methods::{
    VALID_VERSIONS_PACKAGE, VALID_VERSIONS_PITCH, VALID_VERSIONS_PROPAGATION, VALID_VERSIONS_QA,
    VALID_VERSIONS_REVIVE, VALID_VERSIONS_UNIFORMITY,
    VALID_VERSIONS_UNIFORMITY_MISSING_TEMPERATURE,
};
use crate::common::Density;

#[derive(Clone, Debug, PartialEq)]
pub enum Liquid {
    Package(PackagedLiquid),
    Pitch(PitchLiquid),
    Propagation(RegularLiquid),
    QA(KineticsLiquid),
    Revive(RegularLiquid),
    Uniformity(KineticsLiquid),
}

pub const LIQUIDS_DIR: &str = "../yeast/data/liquid";

impl Liquid {
    pub fn from_dir(enable_log: bool) -> Vec<Self> {
        let mut liquids: Vec<Liquid> = Vec::new();
        for year_dir in std::fs::read_dir(LIQUIDS_DIR).unwrap().flatten() {
            if year_dir.metadata().unwrap().is_dir() {
                for year_dir_entry in std::fs::read_dir(year_dir.path()).unwrap().flatten() {
                    let year_dir_entry_metadata = year_dir_entry.metadata().unwrap();
                    if year_dir_entry_metadata.is_dir() {
                        for file in std::fs::read_dir(year_dir_entry.path()).unwrap().flatten() {
                            if let Some(a) = Liquid::kinetics_liquid_from_file(file, enable_log) {
                                liquids.push(a)
                            }
                        }
                    } else if let Ok(contents_toml) = std::fs::read_to_string(year_dir_entry.path())
                    {
                        if let Ok(all_auxiliary_liquids) = contents_toml.parse::<Table>() {
                            if let Some(auxiliary_liquid_table_set) =
                                <&[Value]>::try_single_key(&all_auxiliary_liquids, LIQUID)
                            {
                                for auxiliary_liquid_table_value in auxiliary_liquid_table_set {
                                    if let Value::Table(auxiliary_liquid_table) =
                                        auxiliary_liquid_table_value
                                    {
                                        if let Some(current_auxiliary_liquid) =
                                            Liquid::auxiliary_liquid_from_table(
                                                auxiliary_liquid_table,
                                                enable_log,
                                            )
                                        {
                                            liquids.push(current_auxiliary_liquid)
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if enable_log & liquids.is_empty() {
            println!("No valid liquid entries in directory {LIQUIDS_DIR}.");
        }

        liquids
    }

    pub fn auxiliary_liquid_from_table(
        auxiliary_liquid_table: &Table,
        enable_log: bool,
    ) -> Option<Self> {
        if let Some(id) = <&str>::try_single_key(auxiliary_liquid_table, ID) {
            let id = id.to_owned();
            let parent = match <&str>::try_single_key(auxiliary_liquid_table, PARENT) {
                Some(a) => a.to_owned(),
                _ => {
                    if enable_log {
                        println!(
                            "Warning: auxiliary liquid {id} has no valid parent. Data is skipped."
                        )
                    }
                    return None;
                }
            };
            let protocol_version = match <&str>::try_double_key(
                auxiliary_liquid_table,
                PROTOCOL,
                VERSION,
            ) {
                Some(version_str) => {
                    if let Ok(version) = Version::parse(version_str) {
                        version
                    } else {
                        if enable_log {
                            println!("Warning: auxiliary liquid {id} protocol version could not be parsed. Data is skipped.")
                        }
                        return None;
                    }
                }
                None => {
                    if enable_log {
                        println!(
                        "Warning: auxiliary liquid {id} has no protocol version recorded. Data is skipped."
                    )
                    }
                    return None;
                }
            };
            match <&str>::try_double_key(auxiliary_liquid_table, PROTOCOL, NAME) {
                Some(PACKAGE) => {
                    if enable_log & !VALID_VERSIONS_PACKAGE.contains(&protocol_version) {
                        println!("Warning: auxiliary liquid package run {id} has undocumented protocol version {protocol_version}. Data is processed anyways.");
                    }
                    let batch_count: u8 =
                        match <i64>::try_double_key(auxiliary_liquid_table, BATCH, COUNT) {
                            Some(a_i64) => match a_i64.try_into() {
                                Ok(a_u8) => a_u8,
                                Err(_) => return None,
                            },
                            None => return None,
                        };
                    Some(Self::Package(PackagedLiquid {
                        id,
                        parent,
                        batch_count,
                    }))
                }
                Some(PITCH) => {
                    if enable_log & !VALID_VERSIONS_PITCH.contains(&protocol_version) {
                        println!("Warning: auxiliary pitch liquid {id} has undocumented protocol version {protocol_version}. Data is processed anyways.");
                    }
                    let style = match <&str>::try_single_key(auxiliary_liquid_table, STYLE) {
                        Some(a) => a.to_owned(),
                        _ => {
                            if enable_log {
                                println!(
                                    "Warning: auxiliary pitch liquid {id} has no valid style. Data is skipped."
                                )
                            }
                            return None;
                        }
                    };
                    Some(Self::Pitch(PitchLiquid { id, parent, style }))
                }
                Some(PROPAGATION) => {
                    if enable_log & !VALID_VERSIONS_PROPAGATION.contains(&protocol_version) {
                        println!("Warning: liquid propagate run {id} has undocumented protocol version {protocol_version}. Data is processed anyways.");
                    }
                    Some(Self::Propagation(RegularLiquid { id, parent }))
                }
                Some(REVIVE) => {
                    if enable_log & !VALID_VERSIONS_REVIVE.contains(&protocol_version) {
                        println!("Warning: liquid revive run {id} has undocumented protocol version {protocol_version}. Data is processed anyways.");
                    }
                    Some(Self::Revive(RegularLiquid { id, parent }))
                }
                _ => {
                    if enable_log {
                        println!("Warning: unexpected protocol name for auxiliary liquid {id}. Data skipped.")
                    }
                    None
                }
            }
        } else {
            if enable_log {
                println!("Encountered auxiliary liquid with no id.")
            }
            None
        }
    }

    pub fn kinetics_liquid_from_file(file: std::fs::DirEntry, enable_log: bool) -> Option<Self> {
        if let Some((table_liquid, id)) = try_table_and_id_from_file(file) {
            let parent = match <&str>::try_single_key(&table_liquid, PARENT) {
                Some(a) => a.to_owned(),
                _ => {
                    if enable_log {
                        println!(
                            "Warning: kinetics liquid {id}.toml has no valid parent. File skipped."
                        )
                    }
                    return None;
                }
            };
            let protocol_version = match <&str>::try_double_key(&table_liquid, PROTOCOL, VERSION) {
                Some(version_str) => {
                    if let Ok(version) = Version::parse(version_str) {
                        version
                    } else {
                        if enable_log {
                            println!("Warning: kinetics liquid {id}.toml protocol version could not be parsed. File skipped.")
                        }
                        return None;
                    }
                }
                None => {
                    if enable_log {
                        println!(
                        "Warning: kinetics liquid {id}.toml has no protocol version recorded. File skipped."
                    )
                    }
                    return None;
                }
            };
            match <&str>::try_double_key(&table_liquid, PROTOCOL, NAME) {
                Some(QA) => {
                    if enable_log & !VALID_VERSIONS_QA.contains(&protocol_version) {
                        println!("Warning: liquid QA run {id}.toml has undocumented protocol version {protocol_version}. File is processed anyways.");
                    }
                    let measurement_set = MeasurementSet::try_from(&table_liquid, &id, enable_log);
                    let pitch_time = match <&str>::try_single_key(&table_liquid, TIME) {
                        Some(timestamp_str) => {
                            match PrimitiveDateTime::parse(timestamp_str, &Iso8601::DEFAULT) {
                                Ok(valid_time) => valid_time,
                                Err(_) => {
                                    if enable_log {
                                        println!("Pitch timestamp on measurement in kinetics liquid {id}.toml file is invalid.")
                                    }
                                    return None;
                                }
                            }
                        }
                        None => {
                            if enable_log {
                                println!("No pitch timestamp on measurement in kinetics liquid {id}.toml file.")
                            }
                            return None;
                        }
                    };
                    let temperature = match TemperatureSet::try_from(&table_liquid) {
                        Some(a) => Temperature::Measured(a),
                        None => {
                            if VALID_VERSIONS_UNIFORMITY_MISSING_TEMPERATURE
                                .contains(&protocol_version)
                            {
                                Temperature::Room
                            } else {
                                if enable_log {
                                    println!("Invalid or missing experiment temperature record in liquid {id}.toml for version {protocol_version}. File skipped.")
                                }
                                return None;
                            }
                        }
                    };
                    let stock = <&str>::try_single_key(&table_liquid, STOCK).map(|a| a.to_owned());
                    Some(Self::QA(KineticsLiquid {
                        id,
                        parent,
                        pitch_time,
                        measurement_set,
                        temperature,
                        stock,
                    }))
                }
                Some(UNIFORMITY) => {
                    if enable_log & !VALID_VERSIONS_UNIFORMITY.contains(&protocol_version) {
                        println!("Warning: liquid uniformity run {id}.toml has undocumented protocol version {protocol_version}. File is processed anyways.");
                    }
                    let measurement_set = MeasurementSet::try_from(&table_liquid, &id, enable_log);
                    let pitch_time = match <&str>::try_single_key(&table_liquid, TIME) {
                        Some(timestamp_str) => {
                            match PrimitiveDateTime::parse(timestamp_str, &Iso8601::DEFAULT) {
                                Ok(valid_time) => valid_time,
                                Err(_) => {
                                    if enable_log {
                                        println!("Pitch timestamp on measurement in kinetics liquid {id}.toml file is invalid.")
                                    }
                                    return None;
                                }
                            }
                        }
                        None => {
                            if enable_log {
                                println!("No pitch timestamp on measurement in kinetics liquid {id}.toml file.")
                            }
                            return None;
                        }
                    };
                    let temperature = match TemperatureSet::try_from(&table_liquid) {
                        Some(a) => Temperature::Measured(a),
                        None => {
                            if VALID_VERSIONS_UNIFORMITY_MISSING_TEMPERATURE
                                .contains(&protocol_version)
                            {
                                Temperature::Room
                            } else {
                                if enable_log {
                                    println!("Invalid or missing experiment temperature record in liquid {id}.toml for version {protocol_version}. File skipped.")
                                }
                                return None;
                            }
                        }
                    };
                    let stock = <&str>::try_single_key(&table_liquid, STOCK).map(|a| a.to_owned());
                    Some(Self::Uniformity(KineticsLiquid {
                        id,
                        parent,
                        pitch_time,
                        measurement_set,
                        temperature,
                        stock,
                    }))
                }
                _ => {
                    if enable_log {
                        println!(
                        "Warning: unexpected protocol in kinetics liquid {id}.toml. File skipped."
                    )
                    }
                    None
                }
            }
        } else {
            None
        }
    }

    pub fn id(&self) -> &str {
        match &self {
            Liquid::Package(packaged_liquid) => packaged_liquid.id.as_ref(),
            Liquid::Pitch(pitch_liquid) => pitch_liquid.id.as_ref(),
            Liquid::Propagation(regular_liquid) => regular_liquid.id.as_ref(),
            Liquid::QA(kinetics_liquid) => kinetics_liquid.id.as_ref(),
            Liquid::Revive(regular_liquid) => regular_liquid.id.as_ref(),
            Liquid::Uniformity(kinetics_liquid) => kinetics_liquid.id.as_ref(),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PackagedLiquid {
    pub id: String,
    pub parent: String,
    pub batch_count: u8,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RegularLiquid {
    pub id: String,
    pub parent: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PitchLiquid {
    pub id: String,
    pub parent: String,
    pub style: String,
}

#[derive(Clone, Debug, PartialEq)]
pub struct KineticsLiquid {
    pub id: String,
    pub parent: String,
    pub pitch_time: PrimitiveDateTime,
    pub measurement_set: MeasurementSet,
    pub temperature: Temperature,
    pub stock: Option<String>,
}

pub const COUNTER_SIZE: usize = 25;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Count {
    pub top: [u16; COUNTER_SIZE],
    pub bottom: [u16; COUNTER_SIZE],
}

pub fn try_counter_set(table_count: &Table, counter_where: &str) -> Option<[u16; COUNTER_SIZE]> {
    match <&[Value]>::try_single_key(table_count, counter_where) {
        Some(counter_array) => {
            if counter_array.len() != COUNTER_SIZE {
                return None;
            }
            let mut out: Vec<u16> = Vec::with_capacity(COUNTER_SIZE);
            for element in counter_array.iter() {
                if let Value::Integer(a) = element {
                    if let Ok(a_u16) = (*a).try_into() {
                        out.push(a_u16);
                    } else {
                        return None;
                    }
                } else {
                    return None;
                }
            }
            Some(
                out.try_into()
                    .expect("pre-checked length with valid elements"),
            )
        }
        None => None,
    }
}

impl Count {
    pub fn try_from(table_count: &Table) -> Option<Self> {
        let top = match try_counter_set(table_count, TOP) {
            Some(top_array) => top_array,
            None => return None,
        };
        let bottom = match try_counter_set(table_count, BOTTOM) {
            Some(bottom_array) => bottom_array,
            None => return None,
        };
        Some(Count { top, bottom })
    }
    pub fn total(&self) -> u16 {
        let mut total = 0u16;
        for value in self.top.iter() {
            total += value;
        }
        for value in self.bottom.iter() {
            total += value;
        }
        total
    }
}

impl Density {
    pub fn try_from(table_measurement: &Table, id: &str, enable_log: bool) -> Option<Self> {
        match <&Table>::try_single_key(table_measurement, DENSITY) {
            Some(table_density) => {
                if let Some(REFRACTOMETRY_METHOD) = <&str>::try_single_key(table_density, METHOD) {
                    let value = match <f64>::try_single_key(table_density, VALUE) {
                        Some(a) => a,
                        None => {
                            if enable_log {
                                println!("Missing or malformed density value in one of the measurements in kinetics liquid {id}.toml.")
                            }
                            return None;
                        }
                    };
                    if enable_log
                        & (value.total_cmp(&900.0) == std::cmp::Ordering::Less
                            || value.total_cmp(&1100.0) == std::cmp::Ordering::Greater)
                    {
                        println!("Check density units in kinetics liquid {id}.toml")
                    }
                    let error = match <f64>::try_single_key(table_density, ERROR) {
                        Some(a) => a,
                        None => {
                            if enable_log {
                                println!("Missing or malformed density error in one of the measurements in kinetics liquid {id}.toml.")
                            }
                            return None;
                        }
                    };
                    Some(Self { value, error })
                } else {
                    if enable_log {
                        println!("Unexpected or missing density method in one of the measurements in kinetics liquid {id}.toml.")
                    }
                    None
                }
            }
            None => None,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct DilutionInfo {
    pub sample_volume: f64,
    pub sample_volume_error: f64,
    pub diluent_mass: f64,
    pub diluent_mass_error: f64,
}

impl DilutionInfo {
    pub fn try_from(table_measurement: &Table, id: &str, enable_log: bool) -> Option<Self> {
        if let Some(table_dilution) = <&Table>::try_single_key(table_measurement, DILUTION) {
            if let Some(ML) = <&str>::try_double_key(table_dilution, SAMPLE, UNIT) {
                if let Some(G) = <&str>::try_double_key(table_dilution, DILUENT, UNIT) {
                    let sample_volume = match <f64>::try_double_key(table_dilution, SAMPLE, VALUE) {
                        Some(a) => a,
                        None => {
                            if enable_log {
                                println!("Invalid or missing sample volume in one of the dilutions in kinetics liquid {id}.toml.")
                            }
                            return None;
                        }
                    };
                    let sample_volume_error = match <f64>::try_double_key(
                        table_dilution,
                        SAMPLE,
                        ERROR,
                    ) {
                        Some(a) => a,
                        None => {
                            if enable_log {
                                println!("Invalid or missing sample volume error in one of the dilutions in kinetics liquid {id}.toml.")
                            }
                            return None;
                        }
                    };
                    let diluent_mass = match <f64>::try_double_key(table_dilution, DILUENT, VALUE) {
                        Some(a) => a,
                        None => {
                            if enable_log {
                                println!("Invalid or missing diluent mass in one of the dilutions in kinetics liquid {id}.toml.")
                            }
                            return None;
                        }
                    };
                    let diluent_mass_error = match <f64>::try_double_key(
                        table_dilution,
                        DILUENT,
                        ERROR,
                    ) {
                        Some(a) => a,
                        None => {
                            if enable_log {
                                println!("Invalid or missing diluent mass error in one of the dilutions in kinetics liquid {id}.toml.")
                            }
                            return None;
                        }
                    };
                    Some(Self {
                        sample_volume,
                        sample_volume_error,
                        diluent_mass,
                        diluent_mass_error,
                    })
                } else {
                    if enable_log {
                        println!("Invalid or missing diluent mass measurement unit in one of the dilutions in kinetics liquid {id}.toml.")
                    }
                    None
                }
            } else {
                if enable_log {
                    println!("Invalid or missing sample volume measurement unit in one of the dilutions in kinetics liquid {id}.toml.")
                }
                None
            }
        } else {
            None
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Clustering {
    Low,
    Moderate,
    Strong,
}

impl Clustering {
    pub fn try_from(table_measurement: &Table, id: &str, enable_log: bool) -> Option<Self> {
        match <&str>::try_single_key(table_measurement, CLUSTERING) {
            Some(LOW) => Some(Clustering::Low),
            Some(MODERATE) => Some(Clustering::Moderate),
            Some(STRONG) => Some(Clustering::Strong),
            Some(_) => {
                if enable_log {
                    println!("Invalid clustering record in kinetics liquid {id}.toml file.")
                }
                None
            }
            _ => None,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Measurement {
    pub timestamp: PrimitiveDateTime,
    pub clustering: Option<Clustering>,
    pub count: Option<Count>,
    pub density_refractometer: Option<Density>,
    pub dilution: Option<DilutionInfo>,
}

impl Measurement {
    pub fn try_from(table_measurement: &Table, id: &str, enable_log: bool) -> Option<Self> {
        let timestamp = match <&str>::try_single_key(table_measurement, TIME) {
            Some(timestamp_str) => {
                match PrimitiveDateTime::parse(timestamp_str, &Iso8601::DEFAULT) {
                    Ok(valid_time) => valid_time,
                    Err(_) => {
                        if enable_log {
                            println!("Timestamp on measurement in kinetics liquid {id}.toml file is invalid.")
                        }
                        return None;
                    }
                }
            }
            None => {
                if enable_log {
                    println!("No timestamp on measurement in kinetics liquid {id}.toml file.")
                }
                return None;
            }
        };
        let clustering = Clustering::try_from(table_measurement, id, enable_log);
        let count = match <&Table>::try_single_key(table_measurement, COUNT) {
            Some(table_count) => match Count::try_from(table_count) {
                Some(a) => Some(a),
                None => {
                    if enable_log {
                        println!("Encountered invalid count set in kinetics liquid {id}.toml file.")
                    }
                    None
                }
            },
            None => None,
        };
        let density_refractometer = Density::try_from(table_measurement, id, enable_log);
        let dilution = DilutionInfo::try_from(table_measurement, id, enable_log);
        Some(Measurement {
            timestamp,
            clustering,
            count,
            density_refractometer,
            dilution,
        })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct MeasurementSet {
    pub set: Vec<Measurement>,
}

impl MeasurementSet {
    pub fn try_from(table_liquid: &Table, id: &str, enable_log: bool) -> Self {
        let mut set: Vec<Measurement> = Vec::new();
        if let Some(array_measurement) = <&[Value]>::try_single_key(table_liquid, MEASUREMENT) {
            for entry_measurement in array_measurement.iter() {
                if let Value::Table(table_measurement) = entry_measurement {
                    if let Some(measurement) =
                        Measurement::try_from(table_measurement, id, enable_log)
                    {
                        set.push(measurement)
                    } else if enable_log {
                        println!("Measurement not processable in kinetics liquid {id}.toml file.")
                    }
                } else if enable_log {
                    println!("Unexpected measurement format encountered in kinetics liquid {id}.toml file.")
                }
            }
        }

        MeasurementSet { set }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum Temperature {
    Measured(TemperatureSet),
    Room,
}

#[derive(Clone, Debug, PartialEq)]
pub struct TemperatureSet {
    pub temperature_experiment: f64,
    pub temperature_error: f64,
}

impl TemperatureSet {
    pub fn try_from(table_liquid: &Table) -> Option<Self> {
        if let Some(temperature_experiment) =
            <f64>::try_double_key(table_liquid, TEMPERATURE, VALUE)
        {
            if let Some(temperature_error) = <f64>::try_double_key(table_liquid, TEMPERATURE, ERROR)
            {
                if let Some(a) = <&str>::try_double_key(table_liquid, TEMPERATURE, UNIT) {
                    if a == DEG_C {
                        Some(TemperatureSet {
                            temperature_experiment,
                            temperature_error,
                        })
                    } else {
                        None
                    }
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            None
        }
    }
}
