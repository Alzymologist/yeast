use toml::{Table, Value};

use super::keywords::{CODE, DESCRIPTION, NAME, STRAIN, STYLES};
use super::utils::Extractable;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Strain {
    pub code: String,
    pub name: Option<String>,
    pub description: Option<String>,
    pub styles: Vec<String>,
}

pub const STRAINS_FILE: &str = "../yeast/data/strains.toml";

impl Strain {
    pub fn from_file(enable_log: bool) -> Vec<Self> {
        let mut strains: Vec<Strain> = Vec::new();
        if let Ok(contents_toml) = std::fs::read_to_string(STRAINS_FILE) {
            if let Ok(all_strains) = contents_toml.parse::<Table>() {
                if let Some(strain_table_set) = <&[Value]>::try_single_key(&all_strains, STRAIN) {
                    for strain_table_value in strain_table_set {
                        if let Value::Table(strain_table) = strain_table_value {
                            if let Some(current_strain) =
                                Strain::try_from_strain_table(strain_table, enable_log)
                            {
                                strains.push(current_strain)
                            }
                        }
                    }
                }
            }
        }

        if enable_log & strains.is_empty() {
            println!("No valid strains in file {STRAINS_FILE}.");
        }

        strains
    }

    pub fn try_from_strain_table(strain_table: &Table, enable_log: bool) -> Option<Self> {
        if let Some(code) = <&str>::try_single_key(strain_table, CODE) {
            let name = <&str>::try_single_key(strain_table, NAME).map(|name| name.to_owned());
            let description = <&str>::try_single_key(strain_table, DESCRIPTION)
                .map(|description| description.to_owned());
            let mut styles: Vec<String> = Vec::new();
            if let Some(style_set) = <&[Value]>::try_single_key(strain_table, STYLES) {
                for style_value in style_set {
                    if let Value::String(style) = style_value {
                        styles.push(style.to_string())
                    }
                }
            }
            Some(Strain {
                code: code.to_owned(),
                name,
                description,
                styles,
            })
        } else {
            if enable_log {
                println!("Strain with no code encountered.")
            }
            None
        }
    }
}
