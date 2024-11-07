use toml::{Table, Value};

use super::keywords::{CLONED, ID, NEW, PARENT, SLANT, STRAIN};
use super::utils::Extractable;

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum Slant {
    Origin(OriginSlant),
    Lab(RegularSlant),
    Prod(RegularSlant),
}

pub const SLANTS_DIR: &str = "../yeast/data/slants";

impl Slant {
    pub fn from_dir(enable_log: bool) -> Vec<Self> {
        let mut slants: Vec<Slant> = Vec::new();
        for file in std::fs::read_dir(SLANTS_DIR).unwrap().flatten() {
            if let Ok(contents_toml) = std::fs::read_to_string(file.path()) {
                if let Ok(all_slants) = contents_toml.parse::<Table>() {
                    if let Some(slant_table_set) = <&[Value]>::try_single_key(&all_slants, SLANT) {
                        for slant_table_value in slant_table_set {
                            if let Value::Table(slant_table) = slant_table_value {
                                if let Some(current_slant) =
                                    Slant::try_from_slant_table(slant_table, enable_log)
                                {
                                    slants.push(current_slant)
                                }
                            }
                        }
                    }
                }
            }
        }

        if enable_log & slants.is_empty() {
            println!("No valid slant entries in directory {SLANTS_DIR}.");
        }

        slants
    }

    pub fn try_from_slant_table(slant_table: &Table, enable_log: bool) -> Option<Self> {
        if let Some(id) = <&str>::try_single_key(slant_table, ID) {
            let id = id.to_owned();
            match <&str>::try_single_key(slant_table, PARENT) {
                Some(NEW) => {
                    if let Some(strain_code) = <&str>::try_single_key(slant_table, STRAIN) {
                        Some(Slant::Origin(OriginSlant {
                            id,
                            strain_code: strain_code.to_owned(),
                        }))
                    } else {
                        if enable_log {
                            println!(
                                "Warning: Origin slant {id} has no strain code. Entry skipped."
                            )
                        }
                        None
                    }
                }
                Some(parent) => match <bool>::try_single_key(slant_table, CLONED) {
                    Some(true) => Some(Slant::Prod(RegularSlant {
                        id,
                        parent: parent.to_owned(),
                    })),
                    _ => Some(Slant::Lab(RegularSlant {
                        id,
                        parent: parent.to_owned(),
                    })),
                },
                _ => {
                    if enable_log {
                        println!("Warning: slant {id} has no valid parent or Origin slant identifier. Entry skipped.")
                    }
                    None
                }
            }
        } else {
            if enable_log {
                println!("Encountered slant with no id.")
            }
            None
        }
    }

    pub fn id(&self) -> &str {
        match &self {
            Slant::Origin(origin_slant) => origin_slant.id.as_ref(),
            Slant::Lab(regular_slant) => regular_slant.id.as_ref(),
            Slant::Prod(regular_slant) => regular_slant.id.as_ref(),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OriginSlant {
    pub id: String,
    pub strain_code: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RegularSlant {
    pub id: String,
    pub parent: String,
}
