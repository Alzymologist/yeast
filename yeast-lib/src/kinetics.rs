use multicalc::numerical_integration::{
    integrator::IntegratorMultiVariable, iterative_integration::MultiVariableSolver,
};
use plotters::{
    drawing::IntoDrawingArea,
    prelude::{
        ChartBuilder, DrawingAreaErrorKind, ErrorBar, IntoLogRange, LabelAreaPosition, SVGBackend,
        BLUE, WHITE,
    },
    style::Color,
};
use probability::distribution::{Distribution, Gamma};
use serde::{Deserialize, Serialize};
use std::path::Path;
use time::PrimitiveDateTime;

use crate::{
    calc::{
        ln_n1_probability, t_bend_probability, LN_N1_MEAN_PRIOR, LN_N1_VARIANCE_PRIOR,
        PRECISION_SCALE, PRECISION_SHAPE, T_BEND_MEAN_PRIOR, T_BEND_VARIANCE_PRIOR,
        T_DOUBLING_MEAN_PRIOR, T_DOUBLING_VARIANCE_PRIOR,
    },
    common::Density,
    raw_data_processing::liquid::{KineticsLiquid, Measurement},
    OUTPUT,
};

/// CFU concentration, CFU/mL
#[derive(Clone, Debug)]
pub struct Concentration {
    pub value: f64,
    pub error: f64,
}

/// Concentration data point
#[derive(Clone, Debug)]
pub struct ConcentrationDataPoint {
    pub time_minutes: f64,
    pub concentration: Concentration,
}

/// Density data point
#[derive(Clone, Debug)]
pub struct DensityDataPoint<'a> {
    pub time_minutes: u32,
    pub density: &'a Density,
}

/// Hemocytometer utilized area volume, mL
pub const COUNTER_VOLUME: f64 = 1e-10;

/// Diluent (water) density
pub const DILUENT_DENSITY: f64 = 0.9982067;

impl Measurement {
    pub fn concentration(&self) -> Option<Concentration> {
        if let Some(ref count) = self.count {
            let total_count = count.total() as f64;
            let total_count_error = total_count.sqrt();
            if let Some(ref dilution) = self.dilution {
                let value = total_count * dilution.diluent_mass
                    / (DILUENT_DENSITY * dilution.sample_volume * COUNTER_VOLUME);
                let error = value
                    * ((total_count_error / total_count).powi(2)
                        + (dilution.diluent_mass_error / dilution.diluent_mass).powi(2)
                        + (dilution.sample_volume_error / dilution.sample_volume).powi(2))
                    .sqrt();
                Some(Concentration { value, error })
            } else {
                let value = total_count / COUNTER_VOLUME;
                let error = total_count_error / COUNTER_VOLUME;
                Some(Concentration { value, error })
            }
        } else {
            None
        }
    }

    pub fn concentration_data_point(
        &self,
        start_time: PrimitiveDateTime,
    ) -> Option<ConcentrationDataPoint> {
        if let Some(concentration) = self.concentration() {
            let time_minutes = (self.timestamp - start_time).whole_minutes() as f64;
            Some(ConcentrationDataPoint {
                time_minutes,
                concentration: Concentration {
                    value: concentration.value,
                    error: concentration.error,
                },
            })
        } else {
            None
        }
    }

    pub fn density_data_point(
        &self,
        start_time: PrimitiveDateTime,
    ) -> Option<DensityDataPoint<'_>> {
        if let Some(ref density) = self.density_refractometer {
            let time_minutes = (self.timestamp - start_time).whole_minutes() as u32;
            Some(DensityDataPoint {
                time_minutes,
                density,
            })
        } else {
            None
        }
    }
}

pub const N_POINTS_T_DOUBLING: usize = 30;
pub const N_POINTS_BEFORE_MEAN_PRECISION: usize = 30;
pub const PRECISION_INTEGRAL_FRACTION: f64 = 0.9999;

#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct FitParams {
    pub n_points_t_doubling: usize,
    pub n_points_before_mean_precision: usize,
    pub precision_integral_fraction: f64,
    pub ln_n1_mean_prior: f64,
    pub ln_n1_variance_prior: f64,
    pub t_bend_mean_prior: f64,
    pub t_bend_variance_prior: f64,
    pub t_doubling_mean_prior: f64,
    pub t_doubling_variance_prior: f64,
    pub precision_shape: f64,
    pub precision_scale: f64,
}

impl Default for FitParams {
    fn default() -> Self {
        FitParams {
            n_points_t_doubling: N_POINTS_T_DOUBLING,
            n_points_before_mean_precision: N_POINTS_BEFORE_MEAN_PRECISION,
            precision_integral_fraction: PRECISION_INTEGRAL_FRACTION,
            ln_n1_mean_prior: LN_N1_MEAN_PRIOR,
            ln_n1_variance_prior: LN_N1_VARIANCE_PRIOR,
            t_bend_mean_prior: T_BEND_MEAN_PRIOR,
            t_bend_variance_prior: T_BEND_VARIANCE_PRIOR,
            t_doubling_mean_prior: T_DOUBLING_MEAN_PRIOR,
            t_doubling_variance_prior: T_DOUBLING_VARIANCE_PRIOR,
            precision_shape: PRECISION_SHAPE,
            precision_scale: PRECISION_SCALE,
        }
    }
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct FitPoint {
    pub t_doubling: f64,
    pub precision: f64,
    pub probability: f64,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ConcentrationFitLogistics {
    pub id: String,
    pub fit_params: FitParams,
    pub points: Vec<FitPoint>,
}

impl KineticsLiquid {
    pub fn concentration_fit_logistics(&self) -> Option<ConcentrationFitLogistics> {
        let mut plot_set: Vec<ConcentrationDataPoint> = Vec::new();
        for measurement in self.measurement_set.set.iter() {
            if let Some(concentration_data_point) =
                measurement.concentration_data_point(self.pitch_time)
            {
                if concentration_data_point.concentration.value
                    <= 3.0 * concentration_data_point.concentration.error
                {
                    //TODO: something? this is still a valid data point
                    continue;
                }
                plot_set.push(concentration_data_point)
            }
        }

        if !plot_set.is_empty() {
            // find first valid point, set its time to be a new origin
            let new_origin = plot_set
                .iter()
                .min_by(|a, b| a.time_minutes.total_cmp(&b.time_minutes))
                .unwrap()
                .time_minutes;
            for concentration_data_point in plot_set.iter_mut() {
                concentration_data_point.time_minutes -= new_origin;
            }

            // args: [ln_n1, t_bend, t_doubling, sigma = precision^(-0.5)]
            let func = |args: &[f64; 4]| -> f64 {
                let mut probability = ln_n1_probability(args[0]) * t_bend_probability(args[1]);
                for concentration_data_point in plot_set.iter() {
                    let new_addition = crate::calc::pr_log_logistics(
                        concentration_data_point.time_minutes,
                        args[0],
                        args[1],
                        args[2],
                        args[3],
                        concentration_data_point.concentration.value,
                        concentration_data_point.concentration.error,
                    );
                    probability *= new_addition;
                }
                probability
            };

            let mut points: Vec<FitPoint> = Vec::new();

            let precision_step =
                PRECISION_SHAPE * PRECISION_SCALE / (N_POINTS_BEFORE_MEAN_PRECISION - 1) as f64;
            let mut precision = precision_step;
            //println!("precision step: {precision_step:e}");
            let precision_distribution = Gamma::new(PRECISION_SHAPE, PRECISION_SCALE);
            //println!("precision distribution: {}", precision_distribution.distribution(precision));

            while precision_distribution.distribution(precision) < PRECISION_INTEGRAL_FRACTION {
                //println!("precision: {precision}");
                //println!("precision probability: {}", precision_distribution.density(precision));
                //println!("precision distribution: {}", precision_distribution.distribution(precision));

                let sigma = 1.0 / precision.sqrt();
                //println!("sigma {sigma:?}");

                for point_index_t_doubling in 0..N_POINTS_T_DOUBLING {
                    let t_doubling = T_DOUBLING_MEAN_PRIOR - 3.0 * T_DOUBLING_VARIANCE_PRIOR
                        + (point_index_t_doubling as f64) * (6.0 * T_DOUBLING_VARIANCE_PRIOR)
                            / (N_POINTS_T_DOUBLING - 1) as f64;
                    let point = [
                        LN_N1_MEAN_PRIOR + 3.0 * LN_N1_VARIANCE_PRIOR,
                        T_BEND_MEAN_PRIOR + 3.0 * T_BEND_VARIANCE_PRIOR,
                        t_doubling,
                        sigma,
                    ];

                    let mut integrator = MultiVariableSolver::default();
                    integrator.set_total_iterations(100);
                    let integration_limit = [
                        [
                            LN_N1_MEAN_PRIOR - 3.0 * LN_N1_VARIANCE_PRIOR,
                            LN_N1_MEAN_PRIOR + 3.0 * LN_N1_VARIANCE_PRIOR,
                        ],
                        [
                            T_BEND_MEAN_PRIOR - 3.0 * T_BEND_VARIANCE_PRIOR,
                            T_BEND_MEAN_PRIOR + 3.0 * T_BEND_VARIANCE_PRIOR,
                        ],
                    ];
                    let probability = integrator
                        .get_double_partial(&func, [0, 1], &integration_limit, &point)
                        .unwrap();
                    points.push(FitPoint {
                        t_doubling,
                        precision,
                        probability,
                    })
                }

                precision += precision_step;
            }
            Some(ConcentrationFitLogistics {
                id: self.id.to_owned(),
                fit_params: FitParams::default(),
                points,
            })
        } else {
            None
        }
    }

    pub fn concentration_plot(
        &self,
        strain_code: &str,
    ) -> Result<(), DrawingAreaErrorKind<std::io::Error>> {
        let mut plot_set: Vec<ConcentrationDataPoint> = Vec::new();
        for measurement in self.measurement_set.set.iter() {
            if let Some(concentration_data_point) =
                measurement.concentration_data_point(self.pitch_time)
            {
                plot_set.push(concentration_data_point)
            }
        }

        if !plot_set.is_empty() {
            let directory = format!("{OUTPUT}/{KINETICS}/{COUNT}");
            if !Path::new(&directory).exists() {
                if let Err(e) = std::fs::create_dir_all(&directory) {
                    println!("Error making kinetics data output directory: {e}")
                }
            }
            let output_name = format!("{directory}/{strain_code}_{}_{COUNT}.svg", self.id);
            let root_drawing_area = SVGBackend::new(&output_name, (1024, 768)).into_drawing_area();
            root_drawing_area.fill(&WHITE)?;

            let mut ctx = ChartBuilder::on(&root_drawing_area)
                .set_label_area_size(LabelAreaPosition::Left, 100)
                .set_label_area_size(LabelAreaPosition::Bottom, 60)
                .caption(&self.id, ("sans-serif", 40))
                .build_cartesian_2d(
                    0f64..LENGTH_SHORT_EXPERIMENT as f64,
                    (1e10..5e14).log_scale(),
                )?;

            ctx.configure_mesh()
                .x_desc("relative time, min")
                .axis_desc_style(("sans-serif", 40))
                .x_label_formatter(&|x| format!("{:e}", x))
                .x_label_style(("sans-serif", 20))
                .y_desc("cell concentration, CFU/m^3")
                .y_label_formatter(&|x| format!("{:e}", x))
                .y_label_style(("sans-serif", 20))
                .draw()?;

            ctx.draw_series(plot_set.iter().map(|concentration_data_point| {
                ErrorBar::new_vertical(
                    concentration_data_point.time_minutes,
                    concentration_data_point.concentration.value
                        - concentration_data_point.concentration.error,
                    concentration_data_point.concentration.value,
                    concentration_data_point.concentration.value
                        + concentration_data_point.concentration.error,
                    BLUE.filled(),
                    10,
                )
            }))?;
        }
        Ok(())
    }

    pub fn density_plot(
        &self,
        strain_code: &str,
    ) -> Result<(), DrawingAreaErrorKind<std::io::Error>> {
        let mut plot_set: Vec<DensityDataPoint> = Vec::new();
        let mut largest_time_mark = 0;
        for measurement in self.measurement_set.set.iter() {
            if let Some(density_data_point) = measurement.density_data_point(self.pitch_time) {
                if density_data_point.time_minutes > largest_time_mark {
                    largest_time_mark = density_data_point.time_minutes
                }
                plot_set.push(density_data_point)
            }
        }
        let experiment_length = ExperimentLength::from(largest_time_mark);
        if !plot_set.is_empty() {
            let directory = format!("{OUTPUT}/{KINETICS}/{DENSITY}");
            if !Path::new(&directory).exists() {
                if let Err(e) = std::fs::create_dir_all(&directory) {
                    println!("Error making kinetics data output directory: {e}")
                }
            }
            let output_name = format!(
                "{directory}/{strain_code}_{}_{}_{DENSITY}.svg",
                experiment_length.name(),
                self.id
            );
            let root_drawing_area = SVGBackend::new(&output_name, (1024, 768)).into_drawing_area();
            root_drawing_area.fill(&WHITE)?;

            let mut ctx = ChartBuilder::on(&root_drawing_area)
                .set_label_area_size(LabelAreaPosition::Left, 100)
                .set_label_area_size(LabelAreaPosition::Bottom, 60)
                .caption(&self.id, ("sans-serif", 40))
                .build_cartesian_2d(0f64..experiment_length.x_limit(), 1000f64..1060f64)?;

            ctx.configure_mesh()
                .x_desc("relative time, min")
                .axis_desc_style(("sans-serif", 40))
                .x_label_formatter(&|x| format!("{:e}", x))
                .x_label_style(("sans-serif", 20))
                .y_desc("refractometer density, a.u.")
                .y_label_style(("sans-serif", 20))
                .draw()?;

            ctx.draw_series(plot_set.iter().map(|density_data_point| {
                ErrorBar::new_vertical(
                    density_data_point.time_minutes as f64,
                    density_data_point.density.value - density_data_point.density.error,
                    density_data_point.density.value,
                    density_data_point.density.value + density_data_point.density.error,
                    BLUE.filled(),
                    10,
                )
            }))?;
        }
        Ok(())
    }
}

pub const KINETICS: &str = "kinetics";
pub const COUNT: &str = "count";
pub const DENSITY: &str = "density";

enum ExperimentLength {
    Short,
    Regular,
    Forgotten,
}

impl From<u32> for ExperimentLength {
    fn from(item: u32) -> ExperimentLength {
        if item < LENGTH_SHORT_EXPERIMENT {
            ExperimentLength::Short
        } else if item < LENGTH_REGULAR_EXPERIMENT {
            ExperimentLength::Regular
        } else {
            ExperimentLength::Forgotten
        }
    }
}

impl ExperimentLength {
    fn x_limit(&self) -> f64 {
        match &self {
            ExperimentLength::Short => LENGTH_SHORT_EXPERIMENT as f64,
            ExperimentLength::Regular => LENGTH_REGULAR_EXPERIMENT as f64,
            ExperimentLength::Forgotten => LENGTH_FORGOTTEN_EXPERIMENT as f64,
        }
    }
    fn name(&self) -> &str {
        match &self {
            ExperimentLength::Short => "short",
            ExperimentLength::Regular => "regular",
            ExperimentLength::Forgotten => "forgotten",
        }
    }
}

pub const LENGTH_SHORT_EXPERIMENT: u32 = 4_500;
pub const LENGTH_REGULAR_EXPERIMENT: u32 = 45_000;
pub const LENGTH_FORGOTTEN_EXPERIMENT: u32 = 100_000;
