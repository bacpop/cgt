use statrs::distribution::{Beta,ContinuousCDF};
use statrs::statistics::*;
use clap::Parser;
use std::io::{BufReader, BufRead, BufWriter, Write};
use std::fs::File;
use rand::prelude::{*, Distribution};
use rand::distributions::Uniform;
use indicatif::ProgressIterator;
use std::collections::HashMap;

fn main() {

    let args = Args::parse();
    let breaks = parse_breaks(&args.breaks);

    let completeness: Vec<f64> = load_completeness(&args.completeness, args.completeness_column);
    let n_genomes = completeness.len();

    let countvec = load_counts(&args.matrix);
    let counts: HashMap<String, usize> = countvec.into_iter().collect();

    let mut rng = thread_rng();
    let prior_sample_core = prior_sample(args.beta_param1, args.beta_param2, args.n_samples, breaks.1, true, &mut rng);
    let prior_sample_notcore = prior_sample(args.beta_param1, args.beta_param2, args.n_samples, breaks.1, false, &mut rng);
    let prior_sample_rare = prior_sample(args.beta_param1, args.beta_param2, args.n_samples, breaks.0, false, &mut rng);
    let prior_sample_notrare = prior_sample(args.beta_param1, args.beta_param2, args.n_samples, breaks.0, true, &mut rng);


    let mut obs_sample_core = Vec::with_capacity(args.n_samples);
    let mut obs_sample_notcore = Vec::with_capacity(args.n_samples);
    let mut obs_sample_rare = Vec::with_capacity(args.n_samples);
    let mut obs_sample_notrare = Vec::with_capacity(args.n_samples);

    for i in (0..args.n_samples).progress() {
        obs_sample_core.push(poibin_sample(prior_sample_core[i], &completeness, n_genomes, &mut rng));
        obs_sample_notcore.push(poibin_sample(prior_sample_notcore[i], &completeness, n_genomes, &mut rng));
        obs_sample_rare.push(poibin_sample(prior_sample_rare[i], &completeness, n_genomes, &mut rng));
        obs_sample_notrare.push(poibin_sample(prior_sample_notrare[i], &completeness, n_genomes, &mut rng));

    }

    let ns = args.n_samples as f64;

    let mut prob_core_as_notcore: Vec<f64> = Vec::new();
    let mut prob_notcore_as_core: Vec<f64> = Vec::new();
    let mut prob_rare_as_notrare: Vec<f64> = Vec::new();
    let mut prob_notrare_as_rare: Vec<f64> = Vec::new();

    for i in 0..=n_genomes {
        prob_core_as_notcore.push((obs_sample_core.iter().filter(|x| *x < &i).count() as f64) / ns);
        prob_notcore_as_core.push((obs_sample_notcore.iter().filter(|x| *x >= &i).count() as f64) / ns);
        prob_rare_as_notrare.push((obs_sample_rare.iter().filter(|x| *x > &i).count() as f64) / ns);
        prob_notrare_as_rare.push((obs_sample_notrare.iter().filter(|x| *x <= &i).count() as f64) / ns);
    }

    let rare_threshold = prob_notrare_as_rare.iter().position(|x| x > &args.error).unwrap();
    let core_threshold = prob_core_as_notcore.iter().position(|x| x > &args.error).unwrap();
    
    write_labels(&counts, &rare_threshold, &core_threshold, args.output_file);

    println!("Core threshold: >= {core_threshold} observations or >= {:.2}% frequency", 100.0 * (core_threshold as f64) / (n_genomes as f64));
    println!("Rare threshold: <= {rare_threshold} observations or <= {:.2}% frequency", 100.0 * (rare_threshold as f64) / (n_genomes as f64));


}

fn write_labels(counts: &HashMap<String, usize>, rare_threshold: &usize, core_threshold: &usize, output_path: String) -> () {
    let f = File::create(output_path).unwrap();
    let mut fs = BufWriter::new(f);
    writeln!(fs, "gene\tcount\tlabel");

    let mut label: String;
    for (key, count) in counts.iter() {

        if count >= core_threshold {
            label = String::from("core");
        } else if count <= rare_threshold {
            label = String::from("rare")
        } else {
            label = String::from("middle")
        }
        
        writeln!(fs, "{}\t{}\t{}", key, count, label);
    }
}

fn poibin_sample(p: f64, completeness: &Vec<f64>, n: usize, rng: &mut ThreadRng) -> usize {
    let c: Vec<f64> = completeness.iter().map(|x| x * p).collect();
    let unif_all = Uniform::new(0.0, 1.0);
    let u = unif_all.sample_iter(rng).take(n);
    c.iter().zip(u).filter(|(c, u)| u <= *c).count()
}

fn prior_sample(betap1: f64, betap2: f64, n: usize, brk: f64, lower: bool, rng: &mut ThreadRng) -> Vec<f64> {
    let b = Beta::new(betap1, betap2).unwrap();
    let cdf_val = b.cdf(brk);

    let uniform_samples = match lower {
        true => {uniform_sample(cdf_val, 1.0, n, rng)},
        false => {uniform_sample(0.0, cdf_val, n, rng)}
    };

    uniform_samples.iter().map(|x| b.inverse_cdf(*x)).collect()
}

fn uniform_sample(lb: f64, ub: f64, n: usize, rng: &mut ThreadRng) -> Vec<f64> {
    let unif = Uniform::new(lb, ub);
    let mut out: Vec<f64> = Vec::new();

    for _ in 0..n {
        out.push(rng.sample(unif));
    }

    out
}

fn load_completeness(filename: &str, col_index: usize) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    let file = BufReader::new(File::open(filename).expect("Completeness file not found"));
    for line in file.lines().skip(1) {
        if let Ok(x) = line {
            out.push(x.split_whitespace().skip(col_index - 1).next().unwrap().parse::<f64>().unwrap() / 100.0);
        }
    }
    out
}

fn load_counts<'a>(filename: &str) -> Vec<(String, usize)> {

    let mut out: Vec<(String, usize)> = Vec::new();
    let file = BufReader::new(File::open(filename).expect("Counts file not found"));

    for line in file.lines().skip(1) {
        if let Ok(x) = line {
            let mut spl = x.split_whitespace();
            let name: &str = spl.next().unwrap();
            let count: usize = spl.collect::<Vec<&str>>()
            .iter().map(|el| el.parse::<usize>().unwrap()).collect::<Vec<usize>>().iter().sum();

            out.push((String::from(name), count));
        }
    }
    out
}

fn parse_breaks(string_breaks: &str) -> (f64, f64) {
    let x: Vec<f64> = string_breaks.split(',').collect::<Vec<&str>>().iter().map(|el| el.parse::<f64>().unwrap()).collect();
    let mi = x.clone().min();
    let ma = x.max();
    (mi, ma)
}

static DEFAULT_BREAKS: &str = "0.05,0.95";

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct Args {
    pub completeness: String,
    pub matrix: String,
    #[arg(long, default_value_t = String::from("celebrimbor_output.tsv"))]
    pub output_file: String,
    #[arg(long, default_value_t = 1)]
    pub completeness_column: usize,
    #[arg(long, default_value_t = DEFAULT_BREAKS.to_string())]
    pub breaks: String,
    #[arg(long, default_value_t = 0.05)]
    pub error: f64,
    #[arg(long, default_value_t = 10000)]
    pub n_samples: usize,
    #[arg(long, default_value_t = 0.1)]
    pub beta_param1: f64,
    #[arg(long, default_value_t = 0.1)]
    pub beta_param2: f64,
}

