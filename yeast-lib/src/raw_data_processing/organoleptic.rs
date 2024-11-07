use toml::{Table, Value};

use super::keywords::{
    APPEARANCE, CONVERGE, DEG_C, DENSE, DENSITY, ERROR, FLUFFY, HAZY, ID, LIQUID, MATCH, METHOD,
    ORGANOLEPTIC, PARTICIPANTS, SCORES, SLIGHTLY_FLUFFY, SLIGHTLY_HAZY, SMELL, TASTE, TEMPERATURE,
    TRANSPARENT, UNIT, U_TUBE_METHOD, VALUE, YEAST,
};
use super::utils::Extractable;

#[derive(Clone, Debug, PartialEq)]
pub enum Organoleptic {
    Failed(FailedOrganoleptic),
    Passed(PassedOrganoleptic),
}

pub const ORGANOLEPTICS_DIR: &str = "../yeast/data/organoleptic";

impl Organoleptic {
    pub fn from_dir(enable_log: bool) -> Vec<Self> {
        let mut organoleptics: Vec<Organoleptic> = Vec::new();
        for file in std::fs::read_dir(ORGANOLEPTICS_DIR).unwrap().flatten() {
            if let Ok(contents_toml) = std::fs::read_to_string(file.path()) {
                if let Ok(all_organoleptics) = contents_toml.parse::<Table>() {
                    if let Some(organoleptic_table_set) =
                        <&[Value]>::try_single_key(&all_organoleptics, ORGANOLEPTIC)
                    {
                        for organoleptic_table_value in organoleptic_table_set {
                            if let Value::Table(organoleptic_table) = organoleptic_table_value {
                                if let Some(current_organoleptic) =
                                    Organoleptic::try_from_organoleptic_table(
                                        organoleptic_table,
                                        enable_log,
                                    )
                                {
                                    organoleptics.push(current_organoleptic)
                                }
                            }
                        }
                    }
                }
            }
        }

        if enable_log & organoleptics.is_empty() {
            println!("No valid organoleptic entries in directory {ORGANOLEPTICS_DIR}.");
        }

        organoleptics
    }

    pub fn try_from_organoleptic_table(
        organoleptic_table: &Table,
        enable_log: bool,
    ) -> Option<Self> {
        if let Some(id) = <&str>::try_single_key(organoleptic_table, ID) {
            let id = id.to_owned();
            let participants = match Participants::try_from(organoleptic_table) {
                Some(a) => a,
                None => {
                    if enable_log {
                        println!(
                        "Warning: organoleptic {id} has no valid participants list. Entry skipped."
                    )
                    }
                    return None;
                }
            };
            let organoleptic_type_found = match OrganolepticTypeFound::try_from(organoleptic_table)
            {
                Some(a) => a,
                None => {
                    if enable_log {
                        println!(
                        "Warning: organoleptic {id} has no valid organoleptic type. Entry skipped."
                    )
                    }
                    return None;
                }
            };
            let organoleptic_type = organoleptic_type_found.ty;
            if enable_log {
                let participants_number = participants.set.len();
                match organoleptic_type {
                    OrganolepticType::QA { .. } => {
                        if participants_number > 1 {
                            println!("Warning: organoleptic {id} has unusual number of participants ({participants_number}) for QA test. Entry is processed anyways.");
                        }
                    }
                    OrganolepticType::Uniformity { .. } => {
                        if participants_number < 3 {
                            println!("Warning: organoleptic {id} has unusual number of participants ({participants_number}) for Uniformity test. Entry is processed anyways.");
                        }
                    }
                }
            }
            if organoleptic_type_found.passed {
                let density_data = DensityAtTemperature::try_from(organoleptic_table);
                if enable_log & density_data.is_none() {
                    println!("Warning: organoleptic {id} has no valid density data. Entry is processed anyways.");
                }

                let appearance = Appearance::try_from(organoleptic_table);
                if enable_log & appearance.is_none() {
                    println!("Warning: organoleptic {id} has no valid appearance data. Entry is processed anyways.");
                }

                let organoleptic_score_sets = OrganolepticScoreSets::try_from(organoleptic_table);

                Some(Self::Passed(PassedOrganoleptic {
                    id,
                    participants,
                    organoleptic_type,
                    density_data,
                    appearance,
                    organoleptic_score_sets,
                }))
            } else {
                Some(Self::Failed(FailedOrganoleptic {
                    id,
                    participants,
                    organoleptic_type,
                }))
            }
        } else {
            if enable_log {
                println!("Encountered organoleptic with no id.")
            }
            None
        }
    }

    pub fn id(&self) -> &str {
        match &self {
            Organoleptic::Failed(failed_organoleptic) => failed_organoleptic.id.as_ref(),
            Organoleptic::Passed(passed_organoleptic) => passed_organoleptic.id.as_ref(),
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct FailedOrganoleptic {
    pub id: String,
    pub participants: Participants,
    pub organoleptic_type: OrganolepticType,
}

#[derive(Clone, Debug, PartialEq)]
pub struct PassedOrganoleptic {
    pub id: String,
    pub participants: Participants,
    pub organoleptic_type: OrganolepticType,
    pub density_data: Option<DensityAtTemperature>,
    pub appearance: Option<Appearance>,
    pub organoleptic_score_sets: Option<OrganolepticScoreSets>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DensityEntry {
    pub density_value: f64,
    pub temperature_value: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DensityAtTemperature {
    pub density_entries: Vec<DensityEntry>,
    pub density_error: f64,
    pub temperature_error: f64,
}

impl DensityAtTemperature {
    pub fn try_from(table_organoleptic: &Table) -> Option<Self> {
        let mut density_entries: Vec<DensityEntry> = Vec::new();
        if let Some(density_value) = <f64>::try_double_key(table_organoleptic, DENSITY, VALUE) {
            if let Some(temperature_value) =
                <f64>::try_double_key(table_organoleptic, TEMPERATURE, VALUE)
            {
                density_entries.push(DensityEntry {
                    density_value,
                    temperature_value,
                })
            } else {
                return None;
            }
        } else if let Some(multiple_density_entries) =
            <&Table>::try_double_key(table_organoleptic, DENSITY, VALUE)
        {
            if let Some(multiple_temperature_entries) =
                <&Table>::try_double_key(table_organoleptic, TEMPERATURE, VALUE)
            {
                if multiple_density_entries.len() == multiple_temperature_entries.len() {
                    for (key, value) in multiple_density_entries.iter() {
                        if let Value::Float(density_value) = value {
                            match multiple_temperature_entries.get(key) {
                                Some(Value::Float(temperature_value)) => {
                                    density_entries.push(DensityEntry {
                                        density_value: *density_value,
                                        temperature_value: *temperature_value,
                                    })
                                }
                                _ => return None,
                            }
                        }
                    }
                } else {
                    return None;
                }
            } else {
                return None;
            }
        }
        let density_error = match <f64>::try_double_key(table_organoleptic, DENSITY, ERROR) {
            Some(a) => a,
            None => return None,
        };
        match <&str>::try_double_key(table_organoleptic, DENSITY, METHOD) {
            Some(a) => {
                if a != U_TUBE_METHOD {
                    return None;
                }
            }
            None => return None,
        };
        let temperature_error = match <f64>::try_double_key(table_organoleptic, TEMPERATURE, ERROR)
        {
            Some(a) => a,
            None => return None,
        };
        match <&str>::try_double_key(table_organoleptic, TEMPERATURE, UNIT) {
            Some(a) => {
                if a != DEG_C {
                    return None;
                }
            }
            None => return None,
        };

        if density_entries.is_empty() {
            None
        } else {
            Some(Self {
                density_entries,
                density_error,
                temperature_error,
            })
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OrganolepticTypeFound {
    pub passed: bool,
    pub ty: OrganolepticType,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum OrganolepticType {
    Uniformity,
    QA,
}

impl OrganolepticTypeFound {
    pub fn try_from(table_organoleptic: &Table) -> Option<Self> {
        if let Some(passed_smell) = <bool>::try_double_key(table_organoleptic, CONVERGE, SMELL) {
            <bool>::try_double_key(table_organoleptic, CONVERGE, TASTE).map(|passed_taste| Self {
                passed: passed_smell & passed_taste,
                ty: OrganolepticType::Uniformity,
            })
        } else if let Some(passed_smell) = <bool>::try_double_key(table_organoleptic, MATCH, SMELL)
        {
            <bool>::try_double_key(table_organoleptic, MATCH, TASTE).map(|passed_taste| Self {
                passed: passed_smell & passed_taste,
                ty: OrganolepticType::QA,
            })
        } else {
            None
        }
    }
}

pub const SCORES_LEN: usize = 8;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Score {
    None,
    Mild,
    Strong,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ScoreSet {
    pub clean: Score,
    pub alcohol: Score,
    pub bitter: Score,
    pub tart: Score,
    pub fruits_berries: Score,
    pub flowers: Score,
    pub spice: Score,
    pub nuts: Score,
}

impl ScoreSet {
    fn try_from(table_organoleptic: &Table, key: &str) -> Option<Self> {
        if let Some(array_scores) = <&[Value]>::try_double_key(table_organoleptic, SCORES, key) {
            if array_scores.len() == SCORES_LEN {
                let mut scores_collector = Vec::with_capacity(SCORES_LEN);
                for array_element in array_scores.iter() {
                    match array_element {
                        Value::Integer(0) => scores_collector.push(Score::None),
                        Value::Integer(1) => scores_collector.push(Score::Mild),
                        Value::Integer(2) => scores_collector.push(Score::Strong),
                        _ => return None,
                    }
                }
                Some(Self {
                    clean: scores_collector[0],
                    alcohol: scores_collector[1],
                    bitter: scores_collector[2],
                    tart: scores_collector[3],
                    fruits_berries: scores_collector[4],
                    flowers: scores_collector[5],
                    spice: scores_collector[6],
                    nuts: scores_collector[7],
                })
            } else {
                None
            }
        } else {
            None
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct OrganolepticScoreSets {
    pub smell: ScoreSet,
    pub taste: ScoreSet,
}

impl OrganolepticScoreSets {
    pub fn try_from(table_organoleptic: &Table) -> Option<Self> {
        if let Some(smell) = ScoreSet::try_from(table_organoleptic, SMELL) {
            ScoreSet::try_from(table_organoleptic, TASTE).map(|taste| Self { smell, taste })
        } else {
            None
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum LiquidAppearance {
    Transparent,
    SlightlyHazy,
    Hazy,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum YeastAppearance {
    Dense,
    SlightlyFluffy,
    Fluffy,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Appearance {
    pub liquid: LiquidAppearance,
    pub yeast: YeastAppearance,
}

impl Appearance {
    pub fn try_from(table_organoleptic: &Table) -> Option<Self> {
        let liquid = match <&str>::try_double_key(table_organoleptic, APPEARANCE, LIQUID) {
            Some(TRANSPARENT) => LiquidAppearance::Transparent,
            Some(SLIGHTLY_HAZY) => LiquidAppearance::SlightlyHazy,
            Some(HAZY) => LiquidAppearance::Hazy,
            _ => return None,
        };
        let yeast = match <&str>::try_double_key(table_organoleptic, APPEARANCE, YEAST) {
            Some(DENSE) => YeastAppearance::Dense,
            Some(SLIGHTLY_FLUFFY) => YeastAppearance::SlightlyFluffy,
            Some(FLUFFY) => YeastAppearance::Fluffy,
            _ => return None,
        };
        Some(Self { liquid, yeast })
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Participants {
    pub set: Vec<String>,
}

impl Participants {
    pub fn try_from(table_organoleptic: &Table) -> Option<Self> {
        if let Some(array_participants) =
            <&[Value]>::try_single_key(table_organoleptic, PARTICIPANTS)
        {
            let mut set: Vec<String> = Vec::with_capacity(array_participants.len());
            for element in array_participants.iter() {
                if let Value::String(participant) = element {
                    set.push(participant.to_owned())
                } else {
                    return None;
                }
            }
            Some(Self { set })
        } else {
            None
        }
    }
}
