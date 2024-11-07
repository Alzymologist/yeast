use semver::Version;
use serde::Serialize;
use toml::{Table, Value};

use super::keywords::{
    ALE, ID, LAGER, NAME, OUTCOME, PARENT, PLATE, PLATING, PROTOCOL, THERMAL, VERSION,
};
use super::utils::Extractable;
use super::valid_methods::{VALID_VERSIONS_PLATING, VALID_VERSIONS_THERMAL};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Plate {
    Regular(RegularPlate),
    Thermal(ThermalPlate),
}

pub const PLATES_DIR: &str = "../yeast/data/plates";

impl Plate {
    pub fn from_dir(enable_log: bool) -> Vec<Self> {
        let mut plates: Vec<Plate> = Vec::new();
        for file in std::fs::read_dir(PLATES_DIR).unwrap().flatten() {
            if let Ok(contents_toml) = std::fs::read_to_string(file.path()) {
                if let Ok(all_plates) = contents_toml.parse::<Table>() {
                    if let Some(plate_table_set) = <&[Value]>::try_single_key(&all_plates, PLATE) {
                        for plate_table_value in plate_table_set {
                            if let Value::Table(plate_table) = plate_table_value {
                                if let Some(current_plate) =
                                    Plate::try_from_plate_table(plate_table, enable_log)
                                {
                                    plates.push(current_plate)
                                }
                            }
                        }
                    }
                }
            }
        }

        if enable_log & plates.is_empty() {
            println!("No valid plate entries in directory {PLATES_DIR}.");
        }

        plates
    }

    pub fn try_from_plate_table(plate_table: &Table, enable_log: bool) -> Option<Self> {
        if let Some(id) = <&str>::try_single_key(plate_table, ID) {
            let id = id.to_owned();
            let parent = match <&str>::try_single_key(plate_table, PARENT) {
                Some(a) => a.to_owned(),
                _ => {
                    if enable_log {
                        println!("Warning: plate {id} has no valid parent. Entry skipped.")
                    }
                    return None;
                }
            };
            let protocol_version = match <&str>::try_double_key(plate_table, PROTOCOL, VERSION) {
                Some(version_str) => {
                    if let Ok(version) = Version::parse(version_str) {
                        version
                    } else {
                        if enable_log {
                            println!("Warning: plate {id} protocol version could not be parsed. Entry skipped.")
                        }
                        return None;
                    }
                }
                None => {
                    if enable_log {
                        println!(
                            "Warning: plate {id} has no protocol version recorded. Entry skipped."
                        )
                    }
                    return None;
                }
            };
            match <&str>::try_double_key(plate_table, PROTOCOL, NAME) {
                Some(PLATING) => {
                    if enable_log & !VALID_VERSIONS_PLATING.contains(&protocol_version) {
                        println!("Warning: regular plate {id} has undocumented protocol version {protocol_version}. Entry is processed anyways.");
                    }
                    Some(Self::Regular(RegularPlate { id, parent }))
                }
                Some(THERMAL) => {
                    let thermal_outcome = match <&str>::try_single_key(plate_table, OUTCOME) {
                        Some(ALE) => ThermalOutcome::Ale,
                        Some(LAGER) => ThermalOutcome::Lager,
                        _ => {
                            if enable_log {
                                println!("Warning: thermal plate {id} has no valid thermal outcome. Entry skipped.")
                            }
                            return None;
                        }
                    };
                    if enable_log & !VALID_VERSIONS_THERMAL.contains(&protocol_version) {
                        println!("Warning: thermal plate {id} has undocumented protocol version {protocol_version}. Entry is processed anyways.");
                    }
                    Some(Self::Thermal(ThermalPlate {
                        id,
                        parent,
                        thermal_outcome,
                    }))
                }
                _ => {
                    if enable_log {
                        println!("Warning: plate {id} has no valid protocol name. Entry skipped.")
                    }
                    None
                }
            }
        } else {
            if enable_log {
                println!("Encountered plate with no id.")
            }
            None
        }
    }

    pub fn id(&self) -> &str {
        match &self {
            Plate::Regular(regular_plate) => regular_plate.id.as_ref(),
            Plate::Thermal(thermal_plate) => thermal_plate.id.as_ref(),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RegularPlate {
    pub id: String,
    pub parent: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ThermalPlate {
    pub id: String,
    pub parent: String,
    pub thermal_outcome: ThermalOutcome,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub enum ThermalOutcome {
    Ale,
    Lager,
}
