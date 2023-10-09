use statrs::distribution::{Beta,ContinuousCDF};
use statrs::statistics::*;
use statrs::prec;
use clap::{Parser};
use std::io::{BufReader, BufRead};
use std::fs::File;
use rand::prelude::{*, Distribution};
use rand::distributions::Uniform;
use std::time::Instant;

fn main() {
    let start = Instant::now();
    
    let args = Args::parse();
    let breaks = parse_breaks(&args.breaks);
    // println!("{:?}", breaks);

    let completeness: Vec<f64> = load_completeness(&args.completeness, args.completeness_column);
    let n_genomes = completeness.len();
    // println!("{:?}", completeness);

    let counts: Vec<usize> = load_counts(&args.matrix);
    let n_genes = counts.len();
    // println!("{:?}", counts);

    let prior_sample_core = prior_sample(args.beta_param1, args.beta_param2, args.n_samples, breaks.1, true);
    let prior_sample_notcore = prior_sample(args.beta_param1, args.beta_param2, args.n_samples, breaks.1, false);

    let prior_sample_rare = prior_sample(args.beta_param1, args.beta_param2, args.n_samples, breaks.0, false);
    let prior_sample_notrare = prior_sample(args.beta_param1, args.beta_param2, args.n_samples, breaks.0, true);


    // // println!("{:?}", prior_sample_core);
    let mut rng = thread_rng();
    let unif_all = Uniform::new(0.0, 1.0);
    
    let mut obs_sample_core = Vec::new();
    let mut obs_sample_notcore = Vec::new();
    let mut obs_sample_rare = Vec::new();
    let mut obs_sample_notrare = Vec::new();

    for i in 0..args.n_samples {
        obs_sample_core.push(poibin_sample(prior_sample_core[i], &completeness, n_genomes));
        obs_sample_notcore.push(poibin_sample(prior_sample_notcore[i], &completeness, n_genomes));
        obs_sample_rare.push(poibin_sample(prior_sample_rare[i], &completeness, n_genomes));
        obs_sample_notrare.push(poibin_sample(prior_sample_notrare[i], &completeness, n_genomes));
        
    }

    // println!("{:?}", prior_sample_core);
    // println!("{:?}", obs_sample_core);

    // println!("{:?}", prior_sample_notcore);
    // println!("{:?}", obs_sample_notcore);

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

    // println!("{:?}", prob_core_as_notcore);
    println!("{:?}", prob_notcore_as_core);

    // println!("{:?}", prob_rare_as_notrare);
    // println!("{:?}", prob_notrare_as_rare);

    let mut thresholds: Vec<usize> = Vec::new();
    let quant = vec![0.05, 0.1, 0.25];
    for q in quant {
        thresholds.push(prob_core_as_notcore.iter().position(|x| x > &q).unwrap());
        println!("Probability of assigning rare as not rare <= {}% at >= {:?} gene observations, corresponding probability of assigning not rare as rare: {:?}%",
         q * 100.0, &thresholds.last().unwrap(), prob_notcore_as_core.get(*thresholds.last().unwrap()).unwrap());
    }
    
    let end = Instant::now();

    eprintln!("Done in {}s", end.duration_since(start).as_secs());
    eprintln!("Done in {}ms", end.duration_since(start).as_millis());
    eprintln!("Done in {}ns", end.duration_since(start).as_nanos());


}

fn poibin_sample(p: f64, completeness: &Vec<f64>, n: usize) -> usize {
    let c: Vec<f64> = completeness.iter().map(|x| x * p).collect();
    let mut rng = thread_rng();
    let unif_all = Uniform::new(0.0, 1.0);
    let u = unif_all.sample_iter(&mut rng).take(n);
    c.iter().zip(u).filter(|(c, u)| u <= *c).count()
}

fn prior_sample(betap1: f64, betap2: f64, n: usize, brk: f64, lower: bool) -> Vec<f64> {
    let b = Beta::new(betap1, betap2).unwrap();
    let cdf_val = b.cdf(brk);

    let uniform_samples = match lower {
        true => {uniform_sample(cdf_val, 1.0, n)},
        false => {uniform_sample(0.0, cdf_val, n)}
    };

    uniform_samples.iter().map(|x| b.inverse_cdf(*x)).collect()
}

fn uniform_sample(lb: f64, ub: f64, n: usize) -> Vec<f64> {
    let unif = Uniform::new(lb, ub);
    let mut out: Vec<f64> = Vec::new();
    let mut rng = thread_rng();

    for _ in 0..n {
        out.push(rng.sample(unif));
    }

    out
}

// fn poibinom_sample(completeness: &Vec<f64>, prior_samples: &Vec<f64>) -> Vec<usize> {
//     let mut out: Vec<usize> = Vec::new();

//     for ps in prior_samples {
//        out.push(completeness.iter().fold(0,|acc, el| acc + sample_bernoulli(el * ps)));
//     }

//     out

// }

fn sample_bernoulli(p: f64) -> usize {
    let mut rng = thread_rng();
    if rng.sample(Uniform::new(0.0, 1.0)) <= p {
        1
    } else {
        0
    }
}

fn load_completeness(filename: &str, col_index: usize) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    let file = BufReader::new(File::open(filename).expect("File not found"));
    for line in file.lines().skip(1) {
        if let Ok(x) = line {
            out.push(x.split_whitespace().skip(col_index - 1).next().unwrap().parse::<f64>().unwrap() / 100.0);
        }
    }
    out
}

fn load_counts(filename: &str) -> Vec<usize> {
    let mut out: Vec<usize> = Vec::new();
    let file = BufReader::new(File::open(filename).expect("File not found"));
    
    for line in file.lines().skip(1) {
        if let Ok(x) = line {
            let count: usize = x.split_whitespace().skip(1).collect::<Vec<&str>>()
            .iter().map(|el| el.parse::<usize>().unwrap()).collect::<Vec<usize>>().iter().sum();

            out.push(count);
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

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
pub struct Args {
    pub completeness: String,
    pub completeness_column: usize,
    pub matrix: String,
    pub breaks: String,
    pub error: f64,
    #[arg(default_value_t = 100)]
    pub n_samples: usize,
    #[arg(default_value_t = 0.1)]
    pub beta_param1: f64,
    #[arg(default_value_t = 0.1)]
    pub beta_param2: f64,
}

