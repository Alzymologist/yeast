use toml::{Table, Value};

pub fn try_table_and_id_from_file(file: std::fs::DirEntry) -> Option<(Table, String)> {
    if let Ok(filename) = file.file_name().into_string() {
        if filename.ends_with(".toml") {
            match std::fs::read_to_string(file.path()) {
                Ok(contents_toml) => {
                    match contents_toml.parse::<Table>() {
                        Ok(a) => Some((a, filename.trim_end_matches(".toml").to_string())),
                        Err(_) => {
                            println!("Warning: file {filename} content is not a valid toml. File skipped.");
                            None
                        }
                    }
                }
                Err(_) => {
                    println!("Warning: unable to read file {filename}. File skipped.");
                    None
                }
            }
        } else {
            None
        }
    } else {
        None
    }
}

pub trait Extractable<'a>: Sized {
    fn try_single_key(table: &'a Table, key: &'a str) -> Option<Self>;
    fn try_double_key(table: &'a Table, key_outer: &'a str, key_inner: &'a str) -> Option<Self> {
        match table.get(key_outer) {
            Some(Value::Table(inner_table)) => Self::try_single_key(inner_table, key_inner),
            _ => None,
        }
    }
}
/*
impl<'a> Extractable<'a> for f32 {
    fn try_single_key(table: &'a Table, key: &'a str) -> Option<Self> {
        match table.get(key) {
            Some(Value::Float(ref a)) => Some(*a as f32),
            Some(Value::Integer(ref a)) => match <i64 as TryInto<u16>>::try_into(*a) {
                Ok(b) => Some(b.into()),
                Err(_) => None,
            },
            _ => None,
        }
    }
}
*/

impl<'a> Extractable<'a> for f64 {
    fn try_single_key(table: &'a Table, key: &'a str) -> Option<Self> {
        match table.get(key) {
            Some(Value::Float(ref a)) => Some(*a),
            Some(Value::Integer(ref a)) => match <i64 as TryInto<u32>>::try_into(*a) {
                Ok(b) => Some(b.into()),
                Err(_) => None,
            },
            _ => None,
        }
    }
}

macro_rules! impl_extractable_deref {
    ($($ty: ty, $variant: ident), *) => {
        $(
            impl <'a> Extractable <'a> for $ty {
                fn try_single_key(table: &'a Table, key: &'a str) -> Option<Self> {
                    match table.get(key) {
                        Some(Value::$variant(ref a)) => Some(*a),
                        _ => None,
                    }
                }
            }
        )*
    }
}

impl_extractable_deref!(bool, Boolean);
impl_extractable_deref!(i64, Integer);

macro_rules! impl_extractable_ref {
    ($($ty: ty, $variant: ident), *) => {
        $(
            impl <'a> Extractable <'a> for $ty {
                fn try_single_key(table: &'a Table, key: &'a str) -> Option<Self> {
                    match table.get(key) {
                        Some(Value::$variant(ref a)) => Some(a),
                        _ => None,
                    }
                }
            }
        )*
    }
}

impl_extractable_ref!(&'a str, String);
impl_extractable_ref!(&'a Table, Table);
impl_extractable_ref!(&'a [Value], Array);
