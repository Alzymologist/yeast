use nalgebra::{linalg::SVD, DMatrix, DVector};
use plotters::{
    drawing::IntoDrawingArea,
    prelude::{
        ChartBuilder, Circle, DiscreteRanged, DrawingAreaErrorKind, IntoLinspace,
        LabelAreaPosition, LineSeries, SVGBackend, BLACK, BLUE, RED, WHITE,
    },
    style::Color,
};

use crate::uncertain::Uncertain;
use crate::OUTPUT;

/// Density of unfermented or fermented wort at standard temperature 20 C.
///
/// Most empirical calculations are using density at standard 20 C, whereas most
/// measurements are performed at slightly different temperatures. U-tube
/// measurement provides both the density and the temperature automatically, for
/// hydrometer measurements temperature reading with thermometer could be added.
///
/// There is an empirical correction function for density of unfermented wort,
/// made by fitting the measured density as a function of temperature for pure
/// water, see HBD #963, 9/7/92, by Christopher Lyons
/// (<https://brewery.org/library/HydromCorr0992.html>).
/// This correction seems quite far-fetched even for the unfermented wort
/// containing dissolved sugars, and particularly so for fermented one. However,
/// in lieu of better alternatives (except, of course, making all the
/// measurements in well-thermostated at 20 C liquids), the correction is used.
/// The density error gets higher in case of corrected measurements.
pub type CorrectedDensity = Uncertain;

impl CorrectedDensity {
    pub fn attenuation(density_wort: Self, density_fermented: Self) -> Uncertain {
        let extract_og = extract(density_wort / 1000.0);
        let extract_fg = extract(density_fermented / 1000.0);
        let q = 0.22 + 0.001 * extract_og;
        let real_extract = (q * extract_og + extract_fg) / (1.0 + q);
        ((extract_og - real_extract) / extract_og) * 100.0
    }

    pub fn abw(density_wort: Self, density_fermented: Self) -> Uncertain {
        let og = density_wort / 1000.0;
        let fg = density_fermented / 1000.0;
        76.08 * (og - fg) / (1.775 - og)
    }

    pub fn abv(density_wort: Self, density_fermented: Self) -> Uncertain {
        Self::abw(density_wort, density_fermented) / 0.794
    }
}

pub fn extract(density: CorrectedDensity) -> Uncertain {
    GRAVITY_TO_EXTRACT_COEFFICIENTS[0] * density * density * density
        + GRAVITY_TO_EXTRACT_COEFFICIENTS[1] * density * density
        + GRAVITY_TO_EXTRACT_COEFFICIENTS[2] * density
        + GRAVITY_TO_EXTRACT_COEFFICIENTS[3]
}

pub const GRAVITY_TO_EXTRACT_COEFFICIENTS: [f64; 4] = [182.94, -776.43, 1262.45, -668.962];

#[derive(Clone, Debug, PartialEq)]
pub struct DensityAtTemperature {
    pub density_value: f64,
    pub density_error: f64,
    pub temperature_value: f64,
    pub temperature_error: f64,
}

impl DensityAtTemperature {
    /// Correcting density for wort (i.e. water with mostly dissolved sugars).
    /// Value gets adjusted for temperature proportionally to corresponding pure
    /// water densities.
    /// Density error increases: fit is not perfect, assumption that wort
    /// contains only dissolved sugars is imprecise.
    pub fn correct(&self) -> CorrectedDensity {
        let pure_water_at_measured = pure_water_density(self.temperature_value);
        let value = self.density_value * PURE_WATER_DENSITY_AT_20C / pure_water_at_measured;
        let error = (self.density_error.powi(2) + (CORRECTION_DISPERSION).powi(2)).sqrt();
        CorrectedDensity { value, error }
    }
}

/// Fit coefficients for pure water density data.
pub const PURE_WATER_DENSITY_COEFFICIENTS: [f64; 4] = [
    1.5564259379785908e-5,
    -0.005882602497378819,
    0.01718147181469476,
    1000.0133336822623,
];

/// Standard temperature. Density at this temperature is used for all
/// calculations.
pub const STANDARD_TEMPERATURE: f64 = 20.0;

pub const PURE_WATER_DENSITY_AT_20C: f64 = 998.128436194643;

/// Pure water density calculated at specified temperature using fit
/// coefficients.
pub fn pure_water_density(temperature: f64) -> f64 {
    PURE_WATER_DENSITY_COEFFICIENTS[0] * temperature.powi(3)
        + PURE_WATER_DENSITY_COEFFICIENTS[1] * temperature.powi(2)
        + PURE_WATER_DENSITY_COEFFICIENTS[2] * temperature
        + PURE_WATER_DENSITY_COEFFICIENTS[3]
}

/// Water density data, from CRC Handbook of Chemistry and Physics,
/// Lide D.R. (ed.), 90 ed., chapter 6-4.
pub const WATER_TEMPERATURE_DENSITY_CRC: &[(f64, f64)] = &[
    (0.1, 999.8493),
    (4.0, 999.9750),
    (5.0, 999.9668),
    (10.0, 999.7021),
    (15.0, 999.1016),
    (18.0, 998.5976),
    (20.0, 998.2063),
    (25.0, 997.0480),
    (30.0, 995.6511),
    (35.0, 994.0359),
    (38.0, 992.9695),
    (40.0, 992.2204),
    (45.0, 990.21),
    (50.0, 988.04),
    (55.0, 985.69),
    (60.0, 983.20),
    (65.0, 980.55),
    (70.0, 977.76),
    (75.0, 974.84),
    (80.0, 971.79),
    (85.0, 968.61),
    (90.0, 965.31),
    (95.0, 961.89),
    (99.974, 958.37),
];

pub const DENSITY_PLOT: &str = "water_density_fit";

pub fn fit_density() -> [f64; 4] {
    let x_3_vector = DVector::<f64>::from_vec(
        WATER_TEMPERATURE_DENSITY_CRC
            .iter()
            .map(|(a, _)| a.powi(3))
            .collect(),
    );
    let x_2_vector = DVector::<f64>::from_vec(
        WATER_TEMPERATURE_DENSITY_CRC
            .iter()
            .map(|(a, _)| a.powi(2))
            .collect(),
    );
    let x_1_vector = DVector::<f64>::from_vec(
        WATER_TEMPERATURE_DENSITY_CRC
            .iter()
            .map(|(a, _)| *a)
            .collect(),
    );
    let x_0_vector = DVector::<f64>::repeat(WATER_TEMPERATURE_DENSITY_CRC.len(), 1.0);
    let indep_t =
        DMatrix::from_columns(&[x_3_vector, x_2_vector, x_1_vector, x_0_vector]).transpose();

    let y_vector = DVector::from_vec(
        WATER_TEMPERATURE_DENSITY_CRC
            .iter()
            .map(|(_, b)| *b)
            .collect(),
    );
    let dep = DMatrix::from_columns(&[y_vector]);

    let left_m = indep_t.transpose();

    let left = SVD::new(left_m, true, true);

    let right = dep;
    let fit = left
        .solve(&right, 1e-50)
        .expect("static data, fit was checked");

    [fit[0], fit[1], fit[2], fit[3]]
}

pub fn plot_water_density() -> Result<(), DrawingAreaErrorKind<std::io::Error>> {
    let output_name = format!("../{OUTPUT}/{DENSITY_PLOT}.svg");

    let root_drawing_area = SVGBackend::new(&output_name, (1024, 768)).into_drawing_area();
    root_drawing_area.fill(&WHITE)?;

    let mut ctx = ChartBuilder::on(&root_drawing_area)
        .set_label_area_size(LabelAreaPosition::Left, 100)
        .set_label_area_size(LabelAreaPosition::Right, 100)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .set_label_area_size(LabelAreaPosition::Top, 60)
        .build_cartesian_2d(0f64..100f64, 955f64..1015f64)?;

    ctx.configure_mesh()
        .x_desc("temperature, C")
        .axis_desc_style(("sans-serif", 40))
        .x_label_style(("sans-serif", 20))
        .y_desc("water density, kg/m^3")
        .y_label_style(("sans-serif", 20))
        .draw()?;

    ctx.draw_series(
        WATER_TEMPERATURE_DENSITY_CRC
            .iter()
            .map(|(x, y)| Circle::new((*x, *y), 4, BLUE.filled())),
    )?;

    let x_axis = (0f64..100f64).step(0.01);

    let coefficients = fit_density();

    ctx.draw_series(LineSeries::new(
        x_axis.values().map(|x| {
            (
                x,
                coefficients[0] * x.powi(3)
                    + coefficients[1] * x.powi(2)
                    + coefficients[2] * x
                    + coefficients[3],
            )
        }),
        &BLACK,
    ))?;

    Ok(())
}

/// Measured values for wort density at different temperatures using U-tube.
///
/// At high density values and typical temperature deviations, this could be
/// considered a worst-case scenario.
pub const MEASURED_DENSITY_WORST_CASE: &[(f64, f64)] = &[
    (30.8, 1090.6),
    (28.5, 1091.6),
    (25.5, 1092.6),
    (24.5, 1092.7),
    (23.3, 1092.6),
    (20.0, 1092.5),
    (19.6, 1092.8),
    (19.1, 1093.0),
    (17.2, 1093.7),
    (16.0, 1094.1),
    (15.7, 1093.9),
];

pub const CORRECTION_DISPERSION: f64 = 0.95;

pub fn plot_experimental_density_calculate_correction_dispersion(
) -> Result<f64, DrawingAreaErrorKind<std::io::Error>> {
    let output_name = format!("../{OUTPUT}/measured_density.svg");

    let root_drawing_area = SVGBackend::new(&output_name, (1024, 768)).into_drawing_area();
    root_drawing_area.fill(&WHITE)?;

    let mut ctx = ChartBuilder::on(&root_drawing_area)
        .set_label_area_size(LabelAreaPosition::Left, 100)
        .set_label_area_size(LabelAreaPosition::Right, 100)
        .set_label_area_size(LabelAreaPosition::Bottom, 60)
        .set_label_area_size(LabelAreaPosition::Top, 60)
        .build_cartesian_2d(14f64..32f64, 1090f64..1095f64)?;

    ctx.configure_mesh()
        .x_desc("temperature, C")
        .axis_desc_style(("sans-serif", 40))
        .x_label_style(("sans-serif", 20))
        .y_desc("wort density, kg/m^3")
        .y_label_style(("sans-serif", 20))
        .draw()?;

    ctx.draw_series(
        MEASURED_DENSITY_WORST_CASE
            .iter()
            .map(|(x, y)| Circle::new((*x, *y), 4, BLUE.filled())),
    )?;

    let mut corrected_series: Vec<(f64, f64)> = Vec::new();
    let mut density_at_20c = None;
    let mut dispersion_collector = 0f64;
    for (temperature, measured_density) in MEASURED_DENSITY_WORST_CASE {
        if *temperature == STANDARD_TEMPERATURE {
            density_at_20c = Some(*measured_density)
        } else {
            let pure_water_at_measured = pure_water_density(*temperature);
            let corrected_density =
                *measured_density * PURE_WATER_DENSITY_AT_20C / pure_water_at_measured;
            corrected_series.push((*temperature, corrected_density));
        }
    }

    let density_at_20c =
        density_at_20c.expect("bad dataset, no measurement at {STANDARD_TEMPERATURE}");

    ctx.draw_series(
        corrected_series
            .iter()
            .chain(std::iter::once(&(STANDARD_TEMPERATURE, density_at_20c)))
            .map(|(x, y)| Circle::new((*x, *y), 4, RED.filled())),
    )?;

    for (_, corrected_density) in corrected_series.iter() {
        dispersion_collector += (corrected_density - density_at_20c).powi(2);
    }
    Ok((dispersion_collector / corrected_series.len() as f64).sqrt())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verify_water_density_fit() {
        plot_water_density().unwrap();
        let fit_coefficients_calculated = fit_density();
        assert_eq!(fit_coefficients_calculated, PURE_WATER_DENSITY_COEFFICIENTS);
    }

    #[test]
    fn verify_pure_water_density_at_20c() {
        let pure_water_density_at_20c_calculated = pure_water_density(STANDARD_TEMPERATURE);
        assert_eq!(
            pure_water_density_at_20c_calculated,
            PURE_WATER_DENSITY_AT_20C
        );
    }

    #[test]
    fn verify_correction_dispersion() {
        let correction_dispersion =
            plot_experimental_density_calculate_correction_dispersion().unwrap();
        assert_eq!(
            (correction_dispersion * 100.0).round() / 100.0,
            CORRECTION_DISPERSION
        );
    }

    #[test]
    fn correct_density() {
        let value = DensityAtTemperature {
            density_value: 1070.0,
            density_error: 1.0,
            temperature_value: 35.0,
            temperature_error: 1.0,
        };
        let value_correct = value.correct();
        assert_eq!(value_correct.value.round(), 1074.0);
        assert_eq!(value_correct.error.round(), 1.0);
    }

    #[test]
    fn attenuation() {
        let wort_density = DensityAtTemperature {
            density_value: 1086.0,
            density_error: 1.0,
            temperature_value: 37.0,
            temperature_error: 1.0,
        };
        let corrected_wort_density = wort_density.correct();
        let fermented_density = DensityAtTemperature {
            density_value: 1025.5,
            density_error: 0.1,
            temperature_value: 22.8,
            temperature_error: 0.1,
        };
        let corrected_fermented_density = fermented_density.correct();

        let attenuation =
            CorrectedDensity::attenuation(corrected_wort_density, corrected_fermented_density);
        assert_eq!(attenuation.value.round(), 56.0);
        assert_eq!(attenuation.error.round(), 14.0);

        let abv = CorrectedDensity::abv(corrected_wort_density, corrected_fermented_density);
        assert_eq!((abv.value * 10.0).round() / 10.0, 9.1);
        assert_eq!((abv.error * 10.0).round() / 10.0, 0.2);
    }
}
