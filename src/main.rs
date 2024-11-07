use clap::{Parser, Subcommand};
use yeast_lib::DataSet;

#[derive(Parser)]
#[command(name = "YeastDataProcessing")]
#[command(version = "0.1.0")]
#[command(
    about = "Yeast data analysis",
    long_about = "Processing the data available from directory `../yeast/data`"
)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
    /// Cache kinetics fit data for all available strains
    ///
    /// Add missing caches and update existing ones if they used outdated fit
    /// settings
    CacheAll,

    /// Cache kinetics fit data for listed strains
    ///
    /// Add missing caches and update existing ones if they used outdated fit
    /// settings
    CacheStrain {
        /// Set of strains to use
        strains: Vec<String>,
    },

    /// Check available yeast data for possible typos and missing parts
    CheckData,

    /// Find strain code for sample with provided id
    FindStrain {
        /// Sample id
        id: String,
    },

    /// Plot genealogy data for all available strains
    GenealogyAll {
        /// Plot on separate graphs
        #[arg(short, long)]
        separate: bool,
    },

    /// Plot genealogy data for listed strains
    GenealogyStrain {
        /// Set of strains to use
        strains: Vec<String>,
    },

    /// Plot kinetics fit heatmaps for all available strains
    HeatMapAll,

    /// Plot kinetics fit heatmaps for listed strains
    HeatMapStrain {
        /// Set of strains to use
        strains: Vec<String>,
    },

    /// Plot and fit all available kinetics data
    PlotAll,

    /// Plot and fit kinetics data from listed files
    PlotId {
        /// Set of files to use
        ids: Vec<String>,
    },

    /// Plot and fit kinetics data for listed strains
    PlotStrain {
        /// Set of strains to use
        strains: Vec<String>,
    },

    /// Check strain kinetics uniformity
    ReportUniformityAll,

    /// Check strain kinetics uniformity
    ///
    /// Kinetics uniformity in concurrent experiment, 3 or more data points
    ReportUniformity {
        /// Set of strains to use
        strains: Vec<String>,
    },

    /// Check strain kinetics variability through strain history
    ///
    /// Kinetics retention as strain slants get renewed
    ReportVariability {
        /// Set of strains to use
        strains: Vec<String>,
    },

    /// Check status of all active slants
    ///
    /// Displays most recent lab slant (checked or unchecked), and renewal status
    StatusAll,

    /// Check status of listed slants
    StatusStrain {
        /// Set of strains to use
        strains: Vec<String>,
    },

    /// Generate and export text summary for all available strains
    SummaryAll,

    /// Generate and export text summary for listed strains
    SummaryStrain {
        /// Set of strains to use
        strains: Vec<String>,
    },
}

fn main() {
    let cli = Cli::parse();
    match cli.command {
        Command::CacheAll => {
            let data_set = DataSet::init(false);
            for strain in data_set.strains.iter() {
                data_set.update_caches_by_strain_code(&strain.code)
            }
        }
        Command::CacheStrain { strains } => {
            let data_set = DataSet::init(false);
            for strain in strains.iter() {
                data_set.update_caches_by_strain_code(strain)
            }
        }
        Command::CheckData => {
            let data_set = DataSet::init(true);
            data_set.orphans();
        }
        Command::FindStrain { id } => {
            let data_set = DataSet::init(false);
            match data_set.strain_code_for_id(&id) {
                Some(strain_code) => println!("{id}: corresponding strain code {strain_code}"),
                None => println!("Sample with id {id} not found"),
            }
        }
        Command::GenealogyAll { separate } => {
            let data_set = DataSet::init(false);
            if separate {
                for strain in data_set.strains.iter() {
                    data_set.plot_genealogy(Some(&strain.code))
                }
            } else {
                data_set.plot_genealogy(None)
            }
        }
        Command::GenealogyStrain { strains } => {
            let data_set = DataSet::init(false);
            for strain_code in strains.iter() {
                if data_set.strain_exists(strain_code).is_some() {
                    data_set.plot_genealogy(Some(strain_code))
                } else {
                    println!("{strain_code}: name unknown")
                }
            }
        }
        Command::HeatMapAll => {
            let data_set = DataSet::init(false);
            for strain in data_set.strains.iter() {
                data_set.kinetics_heatmap_by_strain(&strain.code);
            }
        }
        Command::HeatMapStrain { strains } => {
            let data_set = DataSet::init(false);
            for strain_code in strains.iter() {
                if data_set.strain_exists(strain_code).is_some() {
                    data_set.kinetics_heatmap_by_strain(strain_code);
                } else {
                    println!("{strain_code}: name unknown")
                }
            }
        }
        Command::PlotAll => {
            let data_set = DataSet::init(false);
            for strain in data_set.strains.iter() {
                data_set.kinetics_plots_by_strain(&strain.code)
            }
        }
        Command::PlotId { ids } => {
            let data_set = DataSet::init(false);
            for id in ids.iter() {
                data_set.kinetics_plots_by_id(id)
            }
        }
        Command::PlotStrain { strains } => {
            let data_set = DataSet::init(false);
            for strain_code in strains.iter() {
                if data_set.strain_exists(strain_code).is_some() {
                    data_set.kinetics_plots_by_strain(strain_code)
                } else {
                    println!("{strain_code}: name unknown")
                }
            }
        }
        Command::ReportUniformityAll => {
            let data_set = DataSet::init(false);
            for strain in data_set.strains.iter() {
                data_set.report_uniformity(&strain.code);
            }
        }
        Command::ReportUniformity { strains } => {
            let data_set = DataSet::init(false);
            for strain_code in strains.iter() {
                if data_set.strain_exists(strain_code).is_some() {
                    data_set.report_uniformity(strain_code);
                } else {
                    println!("{strain_code}: name unknown")
                }
            }
        }
        Command::ReportVariability { strains } => {
            let data_set = DataSet::init(false);
            for strain_code in strains.iter() {
                if data_set.strain_exists(strain_code).is_some() {
                    data_set.report_variability(strain_code);
                } else {
                    println!("{strain_code}: name unknown")
                }
            }
        }
        Command::StatusAll => {
            let data_set = DataSet::init(false);
            for strain in data_set.strains.iter() {
                println!("{}", data_set.slant_status(&strain.code))
            }
        }
        Command::StatusStrain { strains } => {
            let data_set = DataSet::init(false);
            for strain_code in strains.iter() {
                if data_set.strain_exists(strain_code).is_some() {
                    println!("{}", data_set.slant_status(strain_code))
                } else {
                    println!("{strain_code}: name unknown")
                }
            }
        }
        Command::SummaryAll => {
            let data_set = DataSet::init(false);
            for strain in data_set.strains.iter() {
                data_set.strain_summary_data(strain)
            }
        }
        Command::SummaryStrain { strains } => {
            let data_set = DataSet::init(false);
            for strain_code in strains.iter() {
                if let Some(strain) = data_set.strain_exists(strain_code) {
                    data_set.strain_summary_data(strain)
                } else {
                    println!("{strain_code}: name unknown")
                }
            }
        }
    }
}
