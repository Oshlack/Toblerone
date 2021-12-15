
// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
// Copyright (c) 2021 Andrew Lonsdale tinyt version

use log::{debug, error, info, trace, warn, LevelFilter, SetLoggerError};
use serde::Deserialize;
use env_logger::fmt::Target;
use std::io::Write;
use bio::io::{fasta, fastq};
use docopt::Docopt;
use failure::Error;
use std::{env, fs};
use std::{path::PathBuf, str};

use tinyt::{
    build_index::build_index,
    pseudoaligner,
    pseudoaligner::process_reads,
};
use tinyt::{config, utils};
use crate::config::{LEFT_EXTEND_FRACTION, READ_COVERAGE_THRESHOLD, DEFAULT_ALLOWED_MISMATCHES, TRIM_VAL};

const PKG_NAME: &'static str = env!("CARGO_PKG_NAME");
const PKG_VERSION: &'static str = env!("CARGO_PKG_VERSION");
const USAGE: &'static str = "
tinyt

Usage:
  tinyt index [--num-threads=<n>] -i <index> <ref-fasta>
  tinyt map [--num-threads=<n>] [--read-length=<r>][--trim-size=<t>] [--skip-trim] [--output=<file>] -i <index> <reads-fastq> [<reads-pair-fastq>]
  tinyt -h | --help | -v | --version

Options:
  -n --num-threads N  Number of worker threads [default: 2]
  -t --trim-size T    Size of base pairs to trim when checking unique read matches [default: 5] 
  -s --skip-trim      Skip the trim read check for unqiue read matches
  -m --mismatch M     Number of allowed mismatches for per read [default: 2]
  -r --read-length R  Provide read length for depth estimation
  -o --output FILE    Output results to file instead of stdout
  -h --help           Show this screen.
  -v --version        Show version.
";

#[derive(Clone, Debug, Deserialize)]
struct Args {
    arg_ref_fasta: String,
    arg_index: String,
    arg_reads_fastq: String,
    arg_reads_pair_fastq: String,
    flag_output: Option<String>,
    flag_num_threads: usize,

    cmd_index: bool,

    cmd_map: bool,
    flag_trim_size: usize,
    flag_mismatch: usize,
    flag_skip_trim: bool,
    flag_read_length: Option<usize>,

    flag_version: bool,
    flag_v: bool,
}

fn main() -> Result<(), Error> {
    let mut args: Args = Docopt::new(USAGE)
        .and_then(|d| d.deserialize())
        .unwrap_or_else(|e| e.exit());



   // let level = log::LevelFilter::Info;
//      pretty_env_logger::formatted_timed_builder().target(Target::Stdout).init();
 // initialize logger
     if env::var_os("RUST_LOG").is_none() {                                                                                                                                                                                                                                      
       env::set_var("RUST_LOG", "tinyt=info");                                                                                                                                                                                                                             
     }   
    pretty_env_logger::init_timed();

    

	if args.flag_version || args.flag_v {
        println! {"{} {}", PKG_NAME, PKG_VERSION};
        return Ok(());
    }

 if args.flag_skip_trim && args.flag_trim_size != 2 {
        warn! {"--trim-size has no effect when --skip-trim used"};
    }
 if args.flag_trim_size ==0 {
        warn! {"--trim-size of 0 implies --skip-trim"};
       args.flag_skip_trim = true; 
    }

if args.flag_trim_size <= args.flag_mismatch, {
        warn! {"--trim-size less than or equal to mismatch may be ineffective"};
    }

    debug!("Command line args:\n{:?}", args);



    if args.cmd_index {
        info!("Building index from fasta: {}",&args.arg_ref_fasta);
        let fasta = fasta::Reader::from_file(&args.arg_ref_fasta)?;
        let (seqs, tx_names, tx_gene_map, gene_length_map) = utils::read_transcripts(fasta)?;
        let index =
            build_index::<config::KmerType>(&seqs, &tx_names, &tx_gene_map,&gene_length_map,  args.flag_num_threads)?;
        info!("Finished building index!");

        info!("Writing index to disk");
        utils::write_obj(&index, &args.arg_index)?;
        info!("Finished writing index!");
        info!("Total equivalence classes: {}", index.dbg.len() );

        use debruijn::Mer;
        use std::fs::File;
        use std::io::Write;


        // output index statistics to "args.arg_index.index.ec.csv"
        let mut w = File::create(format!("{}.ec.csv", &args.arg_index)).unwrap();
        let mut unique_ec = 0;
            //println!("EC,SeqLength,TranscriptCount,TranscriptNames");
            writeln!(&mut w,"EC,SeqLength,TranscriptCount,TranscriptNames").unwrap();
        for e in index.dbg.iter_nodes() {
            let eqid = e.data();
            let eq = &index.eq_classes[*eqid as usize];
            if eq.len() == 1 {
		unique_ec += 1;
            writeln!(&mut w,"EC{},{},{},{:?}", e.node_id, e.sequence().len(), eq.len(), index.tx_names[eq[0] as usize]).unwrap();
	    } else {

            //let nameslist = eq.iter().map(|x| &index.tx_names[*x as usize]).collect().join(",");
            let mut nameslist = String::new();
            for n in eq {

               //println!("{}",index.tx_names[*n as usize]);
            if !nameslist.is_empty() {
            nameslist.push(',');
            }
             nameslist.push_str(&index.tx_names[*n as usize]);


           //    &nameslist.append(index.tx_names[*n as usize]);
            }

            writeln!(&mut w,"EC{},{},{},{:?}", e.node_id, e.sequence().len(), eq.len(),nameslist).unwrap();
            //debug!("EC{}\t{}\t{}\t{:?}", e.node_id, e.sequence().len(), eq.len(),nameslist);

            }

        
        }

            info!("Unique equivalence classes: {}", unique_ec);


    } else if args.cmd_map {
        info!("Reading index from disk");
        let index = utils::read_obj(args.arg_index)?;
        info!("Finished reading index!");

        info!("Mapping reads from fastq");
        let reads = fastq::Reader::from_file(args.arg_reads_fastq)?;
    if args.arg_reads_pair_fastq == ""  {
        info!("Single end reads provided");
        process_reads::<config::KmerType>(reads,None, &index, args.flag_output, args.flag_num_threads,!args.flag_skip_trim,args.flag_trim_size,args.flag_mismatch,args.flag_read_length)?;
    } else {
        info!("Paired end reads provided");
        let reads_pair = fastq::Reader::from_file(args.arg_reads_pair_fastq)?;
        process_reads::<config::KmerType>(reads,Some(reads_pair), &index, args.flag_output, args.flag_num_threads,!args.flag_skip_trim,args.flag_trim_size,args.flag_mismatch,args.flag_read_length)?;
    }


    }

    info!("Done!");
    Ok(())
}
