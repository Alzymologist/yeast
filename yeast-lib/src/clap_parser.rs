use clap::Parser;

#[derive(Parser)]
#[command(name = "YeastDataProcessing")]
#[command(version = "0.1.0")]
#[command(about = "Does awesome things", long_about = None)]
struct Cli {
    #[arg(long)]
    two: String,
    #[arg(long)]
    one: String,
}

