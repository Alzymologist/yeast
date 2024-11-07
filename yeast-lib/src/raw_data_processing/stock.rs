use toml::{Table, Value};

use super::keywords::{
    DENSITY, ERROR, HYDROMETER_METHOD, ID, METHOD, REFRACTOMETRY_METHOD, STOCK, U_TUBE_METHOD,
    VALUE,
};
use super::utils::Extractable;
use crate::common::Density;

#[derive(Clone, Debug, PartialEq)]
pub struct Stock {
    pub id: String,
    pub density: Option<Density>,
    pub density_refractometer: Option<Density>,
}

pub const STOCKS_DIR: &str = "../yeast/data/stock";

impl Stock {
    pub fn from_dir(enable_log: bool) -> Vec<Self> {
        let mut stocks: Vec<Stock> = Vec::new();
        for file in std::fs::read_dir(STOCKS_DIR).unwrap().flatten() {
            if let Ok(contents_toml) = std::fs::read_to_string(file.path()) {
                if let Ok(all_stocks) = contents_toml.parse::<Table>() {
                    if let Some(stock_table_set) = <&[Value]>::try_single_key(&all_stocks, STOCK) {
                        for stock_table_value in stock_table_set {
                            if let Value::Table(stock_table) = stock_table_value {
                                if let Some(current_stock) =
                                    Stock::try_from_stock_table(stock_table, enable_log)
                                {
                                    stocks.push(current_stock)
                                }
                            }
                        }
                    }
                }
            }
        }

        if enable_log & stocks.is_empty() {
            println!("No valid stock entries in directory {STOCKS_DIR}.");
        }

        stocks
    }

    pub fn try_from_stock_table(stock_table: &Table, enable_log: bool) -> Option<Self> {
        if let Some(id) = <&str>::try_single_key(stock_table, ID) {
            let id = id.to_owned();
            if let Some(table_density) = <&Table>::try_single_key(stock_table, DENSITY) {
                if let Some(value) = <f64>::try_single_key(table_density, VALUE) {
                    let density_single_entry = {
                        if enable_log
                            & (value.total_cmp(&900.0) == std::cmp::Ordering::Less
                                || value.total_cmp(&1100.0) == std::cmp::Ordering::Greater)
                        {
                            println!("Check density units in stock {id}")
                        }

                        <f64>::try_single_key(table_density, ERROR)
                            .map(|error| Density { value, error })
                    };
                    if enable_log & density_single_entry.is_none() {
                        println!("Invalid or missing density error in stock {id}.")
                    }

                    if let Some(method_name) = <&str>::try_single_key(table_density, METHOD) {
                        match method_name {
                            HYDROMETER_METHOD | U_TUBE_METHOD => Some(Self {
                                id,
                                density: density_single_entry,
                                density_refractometer: None,
                            }),
                            REFRACTOMETRY_METHOD => {
                                if enable_log {
                                    println!("The only density measure available in stock {id} is refractometer density.")
                                }
                                Some(Self {
                                    id,
                                    density: None,
                                    density_refractometer: density_single_entry,
                                })
                            }
                            _ => {
                                if enable_log {
                                    println!(
                                        "Invalid method name in density measure for stock {id}."
                                    )
                                }
                                None
                            }
                        }
                    } else {
                        if enable_log {
                            println!("No method info in density measure for stock {id}.")
                        }
                        None
                    }
                } else if let Some(value_set) = <&Table>::try_single_key(table_density, VALUE) {
                    if let Some(error_set) = <&Table>::try_single_key(table_density, ERROR) {
                        if let Some(method_set) = <&Table>::try_single_key(table_density, METHOD) {
                            if value_set.len() == error_set.len()
                                && value_set.len() == method_set.len()
                            {
                                let mut density = None;
                                let mut density_refractometer = None;
                                for (key, method) in method_set.iter() {
                                    if let Value::String(method_name) = method {
                                        let value = match <f64>::try_single_key(value_set, key) {
                                            Some(a) => a,
                                            None => {
                                                if enable_log {
                                                    println!(
                                                        "Invalid density value in density measure for stock {id}."
                                                    )
                                                }
                                                return None;
                                            }
                                        };
                                        let error = match <f64>::try_single_key(error_set, key) {
                                            Some(a) => a,
                                            None => {
                                                if enable_log {
                                                    println!(
                                                        "Invalid density error in density measure for stock {id}."
                                                    )
                                                }
                                                return None;
                                            }
                                        };
                                        match method_name.as_str() {
                                            HYDROMETER_METHOD | U_TUBE_METHOD => {
                                                if density.is_some() {
                                                    if enable_log {
                                                        println!(
                                                            "Duplicate method in density measurements for stock {id}."
                                                        )
                                                    }
                                                    return None;
                                                } else {
                                                    density = Some(Density { value, error })
                                                }
                                            }
                                            REFRACTOMETRY_METHOD => {
                                                if density_refractometer.is_some() {
                                                    if enable_log {
                                                        println!(
                                                            "Duplicate method in density measurements for stock {id}."
                                                        )
                                                    }
                                                    return None;
                                                } else {
                                                    density_refractometer =
                                                        Some(Density { value, error })
                                                }
                                            }
                                            _ => {
                                                if enable_log {
                                                    println!(
                                                        "Invalid method name in density measure for stock {id}."
                                                    )
                                                }
                                                return None;
                                            }
                                        }
                                    } else {
                                        if enable_log {
                                            println!(
                                                "Invalid method name in density measure for stock {id}."
                                            )
                                        }
                                        return None;
                                    }
                                }
                                Some(Self {
                                    id,
                                    density,
                                    density_refractometer,
                                })
                            } else {
                                return None;
                            }
                        } else {
                            return None;
                        }
                    } else {
                        return None;
                    }
                } else {
                    if enable_log {
                        println!("Invalid or missing density data in stock {id}. Entry skipped.")
                    }
                    None
                }
            } else {
                if enable_log {
                    println!("Warning: no density info in stock {id}. Entry skipped.")
                }
                None
            }
        } else {
            if enable_log {
                println!("Encountered stock with no id.")
            }
            None
        }
    }
}
