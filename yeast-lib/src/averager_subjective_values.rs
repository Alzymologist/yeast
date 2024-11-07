use serde::Serialize;

use crate::raw_data_processing::{
    liquid::Clustering,
    organoleptic::{LiquidAppearance, OrganolepticScoreSets, Score, ScoreSet, YeastAppearance},
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum GeneralizedProperty {
    Alpha,
    Beta,
    Gamma,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum GeneralizedPropertyAveraged {
    Alpha,
    AlphaToBeta,
    Beta,
    BetaToGamma,
    Gamma,
    Undefined,
}

impl GeneralizedProperty {
    pub fn average(set: Vec<Self>) -> Option<GeneralizedPropertyAveraged> {
        if set.is_empty() {
            return None;
        }
        let mut count_alpha = 0usize;
        let mut count_beta = 0usize;
        let mut count_gamma = 0usize;
        for set_element in set.iter() {
            match set_element {
                GeneralizedProperty::Alpha => count_alpha += 1,
                GeneralizedProperty::Beta => count_beta += 1,
                GeneralizedProperty::Gamma => count_gamma += 1,
            }
        }
        if count_alpha > 0 {
            if count_gamma > 0 {
                Some(GeneralizedPropertyAveraged::Undefined)
            } else if count_beta > 0 {
                Some(GeneralizedPropertyAveraged::AlphaToBeta)
            } else {
                Some(GeneralizedPropertyAveraged::Alpha)
            }
        } else if count_beta > 0 {
            if count_gamma > 0 {
                Some(GeneralizedPropertyAveraged::BetaToGamma)
            } else {
                Some(GeneralizedPropertyAveraged::Beta)
            }
        } else {
            Some(GeneralizedPropertyAveraged::Gamma)
        }
    }
}

pub trait Averager: Sized + Into<GeneralizedProperty> {
    type AverageEnum: From<GeneralizedPropertyAveraged>;
    fn average(set: Vec<Self>) -> Option<Self::AverageEnum> {
        let modified_set: Vec<GeneralizedProperty> =
            set.into_iter().map(|element| element.into()).collect();
        GeneralizedProperty::average(modified_set).map(|a| a.into())
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub enum LiquidAppearanceAveraged {
    Transparent,
    TransparentToSlightlyHazy,
    SlightlyHazy,
    SlightlyHazyToHazy,
    Hazy,
    Undefined,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub enum YeastAppearanceAveraged {
    Dense,
    DenseToSlightlyFluffy,
    SlightlyFluffy,
    SlightlyFluffyToFluffy,
    Fluffy,
    Undefined,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub enum ScoreAveraged {
    None,
    NoneToMild,
    Mild,
    MildToStrong,
    Strong,
    Undefined,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub enum ClusteringAveraged {
    Low,
    LowToModerate,
    Moderate,
    ModerateToStrong,
    Strong,
    Undefined,
}

macro_rules! into_averager {
    ($($property: ty, $averaged_property: ty, $alpha: ident, $alpha_to_beta: ident, $beta: ident, $beta_to_gamma: ident, $gamma: ident), *) => {
        $(
            impl From<$property> for GeneralizedProperty {
                fn from(item: $property) -> Self {
                    match item {
                        <$property>::$alpha => GeneralizedProperty::Alpha,
                        <$property>::$beta => GeneralizedProperty::Beta,
                        <$property>::$gamma => GeneralizedProperty::Gamma,
                    }
                }
            }
            impl From<GeneralizedPropertyAveraged> for $averaged_property {
                fn from(item: GeneralizedPropertyAveraged) -> Self {
                    match item {
                        GeneralizedPropertyAveraged::Alpha => <$averaged_property>::$alpha,
                        GeneralizedPropertyAveraged::AlphaToBeta => <$averaged_property>::$alpha_to_beta,
                        GeneralizedPropertyAveraged::Beta => <$averaged_property>::$beta,
                        GeneralizedPropertyAveraged::BetaToGamma => <$averaged_property>::$beta_to_gamma,
                        GeneralizedPropertyAveraged::Gamma => <$averaged_property>::$gamma,
                        GeneralizedPropertyAveraged::Undefined => <$averaged_property>::Undefined,
                    }
                }
            }
            impl Averager for $property {
                type AverageEnum = $averaged_property;
            }
        )*
    }
}

into_averager!(
    LiquidAppearance,
    LiquidAppearanceAveraged,
    Transparent,
    TransparentToSlightlyHazy,
    SlightlyHazy,
    SlightlyHazyToHazy,
    Hazy
);
into_averager!(
    YeastAppearance,
    YeastAppearanceAveraged,
    Dense,
    DenseToSlightlyFluffy,
    SlightlyFluffy,
    SlightlyFluffyToFluffy,
    Fluffy
);
into_averager!(
    Score,
    ScoreAveraged,
    None,
    NoneToMild,
    Mild,
    MildToStrong,
    Strong
);
into_averager!(
    Clustering,
    ClusteringAveraged,
    Low,
    LowToModerate,
    Moderate,
    ModerateToStrong,
    Strong
);

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct ScoreSetAveraged {
    pub clean: ScoreAveraged,
    pub alcohol: ScoreAveraged,
    pub bitter: ScoreAveraged,
    pub tart: ScoreAveraged,
    pub fruits_berries: ScoreAveraged,
    pub flowers: ScoreAveraged,
    pub spice: ScoreAveraged,
    pub nuts: ScoreAveraged,
}

macro_rules! average_value {
    ($($set: tt, $value: tt), *) => {
        $(
            {
                let scores: Vec<Score> = $set.iter().map(|element| element.$value.clone()).collect();
                Score::average(scores)
            }
        )*
    }
}

impl ScoreSet {
    pub fn average_not_empty(set: Vec<Self>) -> ScoreSetAveraged {
        assert!(!set.is_empty(), "Set must not be empty here.");
        ScoreSetAveraged {
            clean: average_value!(set, clean).expect("set is not empty"),
            alcohol: average_value!(set, alcohol).expect("set is not empty"),
            bitter: average_value!(set, bitter).expect("set is not empty"),
            tart: average_value!(set, tart).expect("set is not empty"),
            fruits_berries: average_value!(set, fruits_berries).expect("set is not empty"),
            flowers: average_value!(set, flowers).expect("set is not empty"),
            spice: average_value!(set, spice).expect("set is not empty"),
            nuts: average_value!(set, nuts).expect("set is not empty"),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct OrganolepticScoreSetsAveraged {
    pub smell: ScoreSetAveraged,
    pub taste: ScoreSetAveraged,
}

impl OrganolepticScoreSets {
    pub fn average(set: Vec<Self>) -> Option<OrganolepticScoreSetsAveraged> {
        if set.is_empty() {
            None
        } else {
            Some(OrganolepticScoreSetsAveraged {
                smell: ScoreSet::average_not_empty(
                    set.iter().map(|element| element.smell.clone()).collect(),
                ),
                taste: ScoreSet::average_not_empty(
                    set.iter().map(|element| element.taste.clone()).collect(),
                ),
            })
        }
    }
}
