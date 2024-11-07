use serde::Serialize;

use crate::averager_subjective_values::{
    ClusteringAveraged, LiquidAppearanceAveraged, OrganolepticScoreSetsAveraged,
    YeastAppearanceAveraged,
};
use crate::raw_data_processing::plate::ThermalOutcome;
use crate::uncertain::MeanWithStandardDeviation;

#[derive(Serialize)]
pub struct StrainSummary {
    pub code: String,
    pub name: Option<String>,
    pub description: Option<String>,
    pub styles: Vec<String>,
    pub liquid_appearance: Option<LiquidAppearanceAveraged>,
    pub yeast_appearance: Option<YeastAppearanceAveraged>,
    pub organoleptic_score_sets: Option<OrganolepticScoreSetsAveraged>,
    pub clustering: Option<ClusteringAveraged>,
    pub real_attenuation: Option<MeanWithStandardDeviation>,
    pub thermal_test: Option<ThermalOutcome>,
    pub average_peak_doubling_time_minutes: Option<f64>,
}
