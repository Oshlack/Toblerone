// Copyright (c) 2018 10x Genomics, Inc. All rights reserved.
// Copyright (c) 2021 Andrew Lonsdale tinyt

use std::collections::HashMap;
use std::fmt::Debug;
use std::io::{self, Write};
use std::path::Path;
use std::sync::{mpsc, Arc, Mutex};
use std::{self, fs::File, str};
use std::fs::OpenOptions;

use bio::io::fastq;
use boomphf::hashmap::NoKeyBoomHashMap;
use crossbeam_utils::thread::scope;
use debruijn::dna_string::DnaString;

use debruijn::graph::DebruijnGraph;
use debruijn::{Dir, Kmer, Mer, Vmer};
use failure::Error;
use log::{debug,info,warn,error};
use serde::{Deserialize, Serialize};
use itertools::Itertools;  // itertools = "0.8"


use crate::config::{LEFT_EXTEND_FRACTION, READ_COVERAGE_THRESHOLD, DEFAULT_ALLOWED_MISMATCHES, TRIM_VAL};
use crate::equiv_classes::EqClassIdType;
use crate::utils;





#[derive(Serialize, Deserialize, Debug)]
pub struct Pseudoaligner<K: Kmer> {
    pub dbg: DebruijnGraph<K, EqClassIdType>,
    pub eq_classes: Vec<Vec<u32>>,
    pub dbg_index: NoKeyBoomHashMap<K, (u32, u32)>,
    pub tx_names: Vec<String>,
    pub tx_gene_mapping: HashMap<String, String>,
    pub gene_length_mapping: HashMap<String, usize>,
}

impl<K: Kmer + Sync + Send> Pseudoaligner<K> {
    pub fn new(
        dbg: DebruijnGraph<K, EqClassIdType>,
        eq_classes: Vec<Vec<u32>>,
        dbg_index: NoKeyBoomHashMap<K, (u32, u32)>,
        tx_names: Vec<String>,
        tx_gene_mapping: HashMap<String, String>,
        gene_length_mapping: HashMap<String, usize>,
    ) -> Pseudoaligner<K> {
        Pseudoaligner {
            dbg,
            eq_classes,
            dbg_index,
            tx_names,
            tx_gene_mapping,
            gene_length_mapping,
        }
    }

    /// Pseudo-align `read_seq` and return a list of nodes that the read was aligned to, with mismatch = 2
    /// Used in build_index where mismatches default can remain 2 for testing
    pub fn map_read_to_nodes(&self, read_seq: &DnaString, nodes: &mut Vec<usize>) -> Option<usize> {
      match self.map_read_to_nodes_with_mismatch(read_seq, nodes, DEFAULT_ALLOWED_MISMATCHES) {
        Some((read_coverage, _mismatches,_read_length)) => Some(read_coverage),
        None => None
      }
    }

    /// Pseudo-align `read_seq` and return a list of nodes that the read was aligned to, with configurable # of allowed mismatches
    pub fn map_read_to_nodes_with_mismatch(&self, read_seq: &DnaString, nodes: &mut Vec<usize>, allowed_mismatches: usize) -> Option<(usize, usize, usize)> {
        let read_length = read_seq.len();
        let mut read_coverage: usize = 0;
        let mut mismatch_count: usize = 0;

        // We're filling out nodes
        nodes.clear();

        let left_extend_threshold = (LEFT_EXTEND_FRACTION * read_length as f64) as usize;

        let mut kmer_pos: usize = 0;
        let kmer_length = K::k();

        if read_seq.len() < kmer_length {
            return None;
        }

        let last_kmer_pos = read_length - kmer_length;
        let mut kmer_lookups = 0;

        {
            // Scan the read for the first kmer that exists in the reference
            let mut find_kmer_match = |kmer_pos: &mut usize| -> Option<(usize, usize)> {
                while *kmer_pos <= last_kmer_pos {
                    let read_kmer = read_seq.get_kmer(*kmer_pos);

                    kmer_lookups += 1;
                    match self.dbg_index.get(&read_kmer) {
                        None => (),
                        Some((nid, offset)) => {
                            // Verify that the kmer actually matches -- the MPHF can have false
                            // positives.
                            let node = self.dbg.get_node(*nid as usize);
                            let ref_seq_slice = node.sequence();
                            let ref_kmer: K = ref_seq_slice.get_kmer(*offset as usize);

                            if read_kmer == ref_kmer {
                                return Some((*nid as usize, *offset as usize));
                            }
                        }
                    };
                    *kmer_pos += 1;
                }

                None
            };

            // extract the first exact matching position of a kmer
            // from the read in the DBG
            let (mut node_id, mut kmer_offset) = match find_kmer_match(&mut kmer_pos) {
                None => (None, None),
                Some((nid, offset)) => (Some(nid), Some(offset)),
            };

	    debug!(" kmer pos is {:?}",kmer_pos);



            // check if we can extend back if there were SNP in every kmer query
            if kmer_pos >= left_extend_threshold && node_id.is_some() {
                let mut last_pos = kmer_pos - 1;
                let mut prev_node_id = node_id.unwrap();
                let mut prev_kmer_offset = if kmer_offset.unwrap() > 0 {
                    kmer_offset.unwrap() - 1
                } else {
                    0
                };

                loop {
                    let node = self.dbg.get_node(prev_node_id);
                    debug!("{:?}, {:?}, {:?}, {:?}, {:?}",
                             node, node.sequence(),
                             self.eq_classes[ *node.data() as usize],
                             prev_kmer_offset, last_pos);

                    // length of remaining read before kmer match
                    let skipped_read = last_pos + 1;

                    // length of the skipped node sequence before kmer match
                    let skipped_ref = prev_kmer_offset + 1;

                    // find maximum extention possbile before fork or eof read
                    let max_matchable_pos = std::cmp::min(skipped_read, skipped_ref);

                    let ref_seq_slice = node.sequence();
                    let mut premature_break = false;
                    let mut matched_bases = 0;
                    let mut seen_snp = 0;
                    for idx in 0..max_matchable_pos {
                        let ref_pos = prev_kmer_offset - idx;
                        let read_offset = last_pos - idx;

                        // compare base by base
                        if ref_seq_slice.get(ref_pos) != read_seq.get(read_offset) {
                            // Record mismatch
                            mismatch_count += 1;

                            // Allowing num_mismatch-SNP
                            seen_snp += 1;
                            if seen_snp > allowed_mismatches {
                                premature_break = true;
                                break;
                            }
                        }

                        matched_bases += 1;
                        read_coverage += 1;
                    }

                    //break the loop if end of read reached or a premature mismatch
                    if last_pos + 1 - matched_bases == 0 || premature_break {
                        break;
                    }

                    // adjust last position
                    last_pos -= matched_bases;

                    // If reached here then a fork is found in the reference.
                    let exts = node.exts();
                    let next_base = read_seq.get(last_pos);
                    if exts.has_ext(Dir::Left, next_base) {
                        // found a left extention.
                        let index = exts
                            .get(Dir::Left)
                            .iter()
                            .position(|&x| x == next_base)
                            .unwrap();

                        let edge = node.l_edges()[index];

                        //update the previous node's id
                        prev_node_id = edge.0;
                        let prev_node = self.dbg.get_node(prev_node_id);
                        prev_kmer_offset = prev_node.sequence().len() - kmer_length;

                        // extract colors
                        nodes.push(prev_node.node_id);
                    } else {
                        break;
                    }
                } // end-loop

		} //end-if

		// could enforce coverage minimum but for varying read lengths, easier to require first match to be close to start of read 
                // could set to mismatch value
            // forward search
            if kmer_pos <= last_kmer_pos {
                loop {
                    let node = self.dbg.get_node(node_id.unwrap());
                    debug!("{:?}, {:?}, {:?}, {:?}",
                             node, node.sequence(),
                             self.eq_classes[ *node.data() as usize],
                             kmer_offset);
                    kmer_pos += kmer_length;
                    read_coverage += kmer_length;

                    // extract colors
                    nodes.push(node.node_id);

                    // length of remaining read after kmer match
                    let remaining_read = read_length - kmer_pos;

                    // length of the remaining node sequence after kmer match
                    let ref_seq_slice = node.sequence();
                    let ref_length = ref_seq_slice.len();
                    let ref_offset = kmer_offset.unwrap() + kmer_length;
                    let informative_ref = ref_length - ref_offset;

                    // find maximum extention possbile before fork or eof read
                    let max_matchable_pos = std::cmp::min(remaining_read, informative_ref);

                    let mut premature_break = false;
                    let mut matched_bases = 0;
                    let mut seen_snp = 0;
                    for idx in 0..max_matchable_pos {
                        let ref_pos = ref_offset + idx;
                        let read_offset = kmer_pos + idx;

                        // compare base by base
                        if ref_seq_slice.get(ref_pos) != read_seq.get(read_offset) {
                            // Record mismatch
                            mismatch_count += 1;

                            // Allowing num_mismatch-SNP
                            seen_snp += 1;
                            if seen_snp > allowed_mismatches {
                                premature_break = true;
                                break;
                            }
                        }

                        matched_bases += 1;
                        read_coverage += 1;
                    }

                    kmer_pos += matched_bases;
                    //break the loop if end of read reached or a premature mismatch
                    if kmer_pos >= read_length {
                        break;
                    }

                    // If reached here then a fork is found in the reference.
                    let exts = node.exts();
                    let next_base = read_seq.get(kmer_pos);

                    if !premature_break && exts.has_ext(Dir::Right, next_base) {
                        // found a right extention.
                        let index = exts
                            .get(Dir::Right)
                            .iter()
                            .position(|&x| x == next_base)
                            .unwrap();

                        let edge = node.r_edges()[index];

                        //update the next node's id
                        node_id = Some(edge.0);
                        kmer_offset = Some(0);

                        //adjust for kmer_position
                        kmer_pos -= kmer_length - 1;
                        read_coverage -= kmer_length - 1;
                    } else {
                        // can't extend node in dbg extract read using mphf
                        // TODO: might have to check some cases
                        if kmer_pos > last_kmer_pos {
                            // can't search in mphf if no full kmer can be made
                            break;
                        }

                        // get the match through mphf
                        match find_kmer_match(&mut kmer_pos) {
                            None => break,
                            Some((nid, offset)) => {
                                node_id = Some(nid);
                                kmer_offset = Some(offset);
                            }
                        };
                    }
                } // end-loop
            } //end-if

        }

        if nodes.len() == 0 {
            if read_coverage != 0 {
                panic!(
                    "Different read coverage {:?} than num of eqclasses {:?}",
                    nodes.len(),
                    read_coverage
                );
            }
            //println!("lookups: {} -- no hit", kmer_lookups);
            None
        } else {
            //println!("lookups: {} -- cov: {}", kmer_lookups, read_coverage);
            Some((read_coverage, mismatch_count,read_length))
        }
    }

    /// Convert a list of nodes contacted by a read into an equivalence class.
    /// Supply node list in `nodes`. Equivalence class will be written to `eq_class`.
    pub fn nodes_to_eq_class(&self, nodes: &mut Vec<usize>, eq_class: &mut Vec<u32>) {
        eq_class.clear();

        if nodes.len() == 0 {
            return;
        }

        // Sort nodes to get the shorter equivalence class first.
        nodes.sort_by_key(|n| {
            let eqclass_id = self.dbg.get_node(*n).data();
            self.eq_classes[*eqclass_id as usize].len()
        });

        let _lens: Vec<_> = nodes
            .iter()
            .map(|n| {
                let eqclass_id = self.dbg.get_node(*n).data();
                self.eq_classes[*eqclass_id as usize].len()
            })
            .collect();
        //println!("nodes: {:?}, lens: {:?}", nodes, lens);

        // Intersect the equivalence classes
        let first_node = nodes[0];

        //println!("node: {}, seq: {:?}", first_node,  self.dbg.get_node(first_node).sequence());
        let first_color = self.dbg.get_node(first_node).data();
        eq_class.extend(&self.eq_classes[*first_color as usize]);

        for node in nodes.iter().skip(1) {
            let color = self.dbg.get_node(*node).data();
            intersect(eq_class, &self.eq_classes[*color as usize]);
        }
    }

    /// Pseudoalign the `read_seq` to the graph. Returns a tuple of the
    /// eqivalence class, the number of bases aligned on success,
    /// and the number of mismatched bases, or None is no alignment could be found.
    pub fn map_read_with_mismatch(&self, read_seq: &DnaString, allowed_mismatches: usize) -> Option<(Vec<u32>, usize, usize,usize)> {
        let mut nodes = Vec::new();

        match self.map_read_to_nodes_with_mismatch(read_seq, &mut nodes, allowed_mismatches) {
            Some((read_coverage, mismatches,read_length)) => {
                let mut eq_class = Vec::new();
                self.nodes_to_eq_class(&mut nodes, &mut eq_class);
                Some((eq_class, read_coverage, mismatches,read_length))
            }
            None => None,
        }
    }

    /// Pseudoalign the `read_seq` to the graph with # mismatches = 2. Returns a tuple of the
    /// eqivalence class and the number of bases aligned on success, and number of mismatches
    /// or None is no alignment could be found.
    pub fn map_read(&self, read_seq: &DnaString, mismatch_size: usize) -> Option<(Vec<u32>, usize, usize, usize)> {
      match self.map_read_with_mismatch(read_seq, mismatch_size) {
        Some((eq_class, read_coverage, mismatches, read_length)) => Some((eq_class, read_coverage,mismatches, read_length)),
        None => None,
      }
    }
}

/// Compute the intersection of v1 and v2 inplace on top of v1
/// v1 and v2 must be sorted and deduplicated.
pub fn intersect<T: Eq + Ord>(v1: &mut Vec<T>, v2: &[T]) {
    if v1.is_empty() {
        return;
    }

    if v2.is_empty() {
        v1.clear();
    }

    let mut fill_idx1 = 0;
    let mut idx1 = 0;
    let mut idx2 = 0;

    while idx1 < v1.len() && idx2 < v2.len() {
        let rem_slice = &v2[idx2..];
        match rem_slice.binary_search(&v1[idx1]) {
            Ok(pos) => {
                v1.swap(fill_idx1, idx1);
                fill_idx1 += 1;
                idx1 += 1;
                idx2 = pos + 1;
            }
            Err(pos) => {
                idx1 += 1;
                idx2 = pos;
            }
        }
    }
    v1.truncate(fill_idx1);
}




// high level function to call match_read and check revcomp, select best match for unique ec if in doubtc
pub fn match_strands<K: Kmer + Sync + Send>(record: &fastq::Record,  trim: bool,trimsize: usize,mismatchsize: usize, index: &Pseudoaligner<K>  ) -> Option<(Option<(bool,bool,String,Vec<u32>, usize,usize,bool,usize)>,String)> {


 // make next steps a function so can be called for R1 and R2  in a paired end version
                            let dna_string = str::from_utf8(&record.seq()).unwrap();
                            let seq = DnaString::from_dna_string(dna_string);
                            let read_data = index.map_read(&seq,mismatchsize);
                            // Some format: Hit to ec, unique ec, read, eq_class, coverage,mistmatches,mistmatches, trimmed
                            let wrapped_read_data = match_read(read_data,&dna_string.to_owned(),&record.id().to_owned(),trim,trimsize,mismatchsize,seq.len(),index);


                            let dna_string_revcomp =  bio::alphabets::dna::revcomp(record.seq());
                            let dna_string_revcomp = str::from_utf8(&dna_string_revcomp).unwrap();
                            let seq_revcomp = DnaString::from_dna_string(dna_string_revcomp);
			//	info!("test slice {:?}",&dna_string_revcomp[1..40]);
                            let read_data_revcomp = index.map_read(&seq_revcomp,mismatchsize);
                            let wrapped_read_data_revcomp = match_read(read_data_revcomp,&dna_string_revcomp.to_owned(),&record.id().to_owned(),trim,trimsize,mismatchsize,seq.len(),index);

			    // return one of the two reads, with origin
   			let return_result = match (&wrapped_read_data, &wrapped_read_data_revcomp) {


                        (Some((_,_,_ , forward_eq_class, _,_,forward_trim,_)),Some((_,_,_, rev_eq_class,_,_,rev_trim,_)) ) => {

                                if forward_eq_class == rev_eq_class  {
                                        if forward_trim == rev_trim {
                                                Some((wrapped_read_data,"forward-equal".to_owned()))

                                        } else {

                                                if *forward_trim {
                                                        Some((wrapped_read_data_revcomp,"reverse-equal-notrim".to_owned()))

                                                } else {
                                                        Some((wrapped_read_data,"forward-equal-notrim".to_owned()))

                                                }
                                        }

                                } else {
                                        if forward_eq_class.len() > rev_eq_class.len() &&  rev_eq_class.len() > 0 {
                                        Some((wrapped_read_data_revcomp,"reverse-vsN".to_owned()))

                                        } else {
                                                if forward_eq_class.len() == 0 {
                                        Some((wrapped_read_data_revcomp,"reverse".to_owned()))
                                        //Some((wrapped_read_data_revcomp,"reverse-vs0"))

                                                } else {
                                                if rev_eq_class.len() == 0 {
                                        Some((wrapped_read_data,"forward".to_owned()))
                                        //Some((wrapped_read_data,"forward-vs0".to_owned()))
                                        }       else {
                                        Some((wrapped_read_data,"forward-vsN".to_owned()))

}
                                                }

                                        }


                                }

}

                        (_,None) => Some((wrapped_read_data,"forward-none".to_owned())),
                        (None,_) =>  Some((wrapped_read_data_revcomp,"reverse-none".to_owned())),
                        (None,None) => None
                        };



		return return_result;

}




// function to match read to result
pub fn match_read<K: Kmer + Sync + Send>(optional: Option<(Vec<u32>, usize, usize,usize)>, seq: &String,  record_id: &String,  trim: bool, trimsize: usize, mismatchsize: usize,seqlength: usize, index: &Pseudoaligner<K>  ) -> Option<(bool,bool,String,Vec<u32>, usize,usize,bool,usize)> {

        // Some format: Hit to ec, unique ec, read, eq_class, coverage, mismatches,trimmed 

	match optional {
        	Some((eq_class,coverage,mismatches,readlen)) => {
                                    if coverage >= seq.len() && mismatches == 0 &&  eq_class.len() == 1 {
                                    //if coverage >= READ_COVERAGE_THRESHOLD && eq_class.len() == 1 {
                                        debug!("{:?}",&seq);
                                      if trim { // trim check

                                        // check if same unique ec matches using a trimmed version of the read, to catch hanging 1 or 2 bp matches for deletions

                                        let trim_end = seqlength-trimsize;
                                        //println!("{:?}",seq);
                                        debug!("{:?}",&seq[trimsize..trim_end]);
                                        let trim_seq = DnaString::from_dna_string(&seq[trimsize..trim_end]);
                                        let trim_read_data = index.map_read(&trim_seq,mismatchsize);
					//debug!("{} length {} trimlength {}",record_id,seqlength, trim_seq.len()); 

                                        let wrapped_trim_read_data = match trim_read_data {

                                        Some((trim_eq_class, trim_coverage,trim_mismatches,trim_readlen)) => {
                                         //println!("{:?}",&trim_eq_class);

                                          if trim_eq_class == eq_class {
					debug!(" {} no diff trimed {:?} {:?}",record_id,trim_eq_class,eq_class); 
                                       Some((true,true, record_id.to_owned(), eq_class, coverage,mismatches,false,readlen))

                                        } else {
					debug!("{} diff trimed {:?} {:?}",record_id,trim_eq_class,eq_class); 
                                       Some((true,false, record_id.to_owned(), trim_eq_class, coverage,mismatches,true,readlen))
					// use original readength

                                        }

                                        }
                                        // trim these also? to be consistent, yes, indicates different
                                        None =>  {
					debug!("orig {:?}",index.map_read(&DnaString::from_dna_string(&seq),mismatchsize));
					debug!(" trimsize... {:?}",index.map_read(&DnaString::from_dna_string(&seq[trimsize..]),mismatchsize));
					debug!(" ...trim_end {:?}",index.map_read(&DnaString::from_dna_string(&seq[..trim_end]),mismatchsize));
					debug!("coded trim {:?}",trim_read_data);
					debug!(" {} none trimed {:?}  returns None", record_id.to_owned(),eq_class); 
					debug!("orig {:?} trimmed {:?}",&DnaString::from_dna_string(&seq),&DnaString::from_dna_string(&seq[trimsize..trim_end])); 
					//debug!("orig {:?} trimmed {:?}",seq[trimsize..trim_end]); 
					 Some((true,false, record_id.to_owned(), eq_class, coverage,mismatches,true,readlen))
					}
                                        };


                                        wrapped_trim_read_data


                                       } else { 

                                       Some((true,true, record_id.to_owned(), eq_class, coverage,mismatches,false,readlen))

                                       }
                                    
				} else {
                                        Some((true,false, record_id.to_owned(), eq_class, coverage,mismatches,false,readlen))
                                    }
                                }
                                

				None => Some((false,false, record_id.to_owned(), Vec::new(), 0,0,false,0)),
    	}












}



pub fn process_reads<K: Kmer + Sync + Send>(
    reader: fastq::Reader<File>,
    reader_pair: Option<fastq::Reader<File>>,
    index: &Pseudoaligner<K>,
    outfile: Option<String>,
    num_threads: usize,
    trim: bool,
    trimsize: usize,
    mismatchsize: usize,
    read_length: Option<usize>,
) -> Result<(), Error> {
    info!("Done Reading index");
    info!("Starting Multi-threaded Mapping");

  let mut output_file = match outfile{
        Some(filename) => {
    	info!("Output file: {}", &filename);

                 Box::new(OpenOptions::new().append(true)
                             .create_new(true)
                             .open(filename).unwrap()) as Box<dyn Write>

        },
        None => {
    	info!("Output file: STDOUT");
	Box::new(std::io::stdout()) as Box<dyn Write>
	},
    };




    let (tx, rx) = mpsc::sync_channel(num_threads);

    // if paired end, need a differnet atomix_reader
    let is_paired =  if let Some(_) = &reader_pair {
		    true 
	} else {
  		false 
 	};


   let read_length_arg  = match read_length {

	Some(length) => length, 
	None => 0
	};
	
    let check_read = if read_length_arg == 0 { true } else {false}; 



    info!("Spawning {} threads for Mapping.\n", num_threads);
    scope(|scope| {

		if is_paired {

			let readerR2 = match reader_pair {
                                Some(readerR2) => readerR2,
                                None => panic!("Error {:?} in reading fastq R2", reader_pair),
                            };



       let atomic_reader_pair =  Arc::new(Mutex::new(reader.records().zip(readerR2.records()).into_iter()));
//let atomic_reader_pair =         Arc::new(Mutex::new(reader.records()));
info!("Spawning {} threads for Mapping.\n", num_threads);
        for _ in 0..num_threads {
            let tx = tx.clone();
            let reader = Arc::clone(&atomic_reader_pair);

            scope.spawn(move |_| {
                loop {
                    // If work is available, do that work.
                    match utils::get_next_record_pair(&reader) {
                        Some(result_record) => {
				let (record,recordR2) = match result_record {
                                (Ok(record),Ok(recordR2)) => (record,recordR2),
                                (Err(err),_) => panic!("Error {:?} in reading fastq R1", err),
                                (_,Err(err)) => panic!("Error {:?} in reading fastq R2", err),
                            };
	

				    if trim && trimsize > &record.seq().len()/2 {

                                                //panic!("trim too long");std::process::exit(1)
                                                error!{"Trimsize {:?} too long, entire sequence trimmed for read {:?}",trimsize, &record.id()};
                                        	tx.send(None).expect("Could not send data!");
						break;
					}


                            let compared_read_data = match_strands(&record,trim,trimsize,mismatchsize,index);
                            let compared_read_data_R2 = match_strands(&recordR2,trim,trimsize,mismatchsize,index);
			
				// compare EC class result from each read pair, return the best one for unique EC
				//Option<(Option<(bool,bool,String,Vec<u32>, usize,usize,bool)>,String)>
				//println!("R1: {:?}",compared_read_data); 				
				//println!("R2: {:?}",compared_read_data_R2);

				// if none pass through none, if only 1 return that one
				// choose between : shortest EC unless EC is 0
				// if equal, return with either (or emebd?)
				// need to udate donstream to extpct another string 				
				let selected_read =  match (compared_read_data, compared_read_data_R2) {

				(None,None) => None,
				(None, Some((None, _))) => None, 
				(Some((None, _)), None) => None,
				(Some((None, _)), Some((None, _))) => None,
				(Some((Some(read),strandinfo)),Some((None,_))) => {Some((Some(read),strandinfo))},
				(Some((Some(read),strandinfo)),None) => {Some((Some(read),strandinfo))},
				(None,Some((Some(read2),strandinfo2))) => {Some((Some(read2),strandinfo2))},
				(Some((None,_)),Some((Some(read2),strandinfo2))) => {Some((Some(read2),strandinfo2))},
				
				//both reads give some information - preference for unique (or shortest) 
				(Some((Some((mapped,unique,transcript,eq_class,coverage,mismatches,is_trimmed,readlen)),strandinfo)),
				Some((Some((mapped2,unique2,transcript2,eq_class2,coverage2,mismatches2,is_trimmed2,readlen2)),strandinfo2))) => {
		
					 if eq_class == eq_class2 {
				//Some((Some((mapped2,unique2,transcript2,eq_class2,coverage2,mismatches2,is_trimmed2)),"from R1 and R2".to_string()))	
					// TODO mix values, choose best here? at moment defaults to R2 unless R1 trimmed
				// any read that was trimmed is by deficitonio shortest (or eqauka shortest) , and want that reported, otherwise go trough and selecyt shortst
				// TODO 
					if (is_trimmed) {
				Some((Some((mapped,unique,transcript,eq_class,coverage,mismatches,is_trimmed,readlen)),"from R1 and R2".to_string()))	
					} else {
				Some((Some((mapped2,unique2,transcript2,eq_class2,coverage2,mismatches2,is_trimmed2,readlen2)),"from R1 and R2".to_string()))	
					}



					} else {
					// not equal - return the trimmed version first, otherwise shortest

					if (is_trimmed && is_trimmed2) {
                                //println!(". choose R1 and R2 -134 trmmed {:?} vs {:?} ",eq_class2, eq_class);
                                Some((Some((mapped,unique,transcript,eq_class,coverage,mismatches,is_trimmed,readlen)),"from R1 and R2".to_string()))

					} else if (is_trimmed) {
                                //println!("a choose R1 -134 trmmed {:?} vs {:?} ",eq_class2, eq_class);
                                Some((Some((mapped,unique,transcript,eq_class,coverage,mismatches,is_trimmed,readlen)),"from R1".to_string()))
                                        } else if (is_trimmed2) { 
                                //println!("b choose R2 -134 trmmed {:?} vs {:?} ",eq_class2, eq_class);
                                Some((Some((mapped2,unique2,transcript2,eq_class2,coverage2,mismatches2,is_trimmed2,readlen2)),"from R2".to_string()))

					} else if eq_class.len() == 0 {
						// choose shortest but not 0
                                //println!("1 choose R2 -134 {:?} vs {:?} ",eq_class2, eq_class);
 
				Some((Some((mapped2,unique2,transcript2,eq_class2,coverage2,mismatches2,is_trimmed2,readlen2)),"from R2".to_string()))	
					} else if eq_class2.len() == 0 { 
                                //println!(" 2choose R1 -134 {:?} vs {:?}",eq_class,eq_class2);
				Some((Some((mapped,unique,transcript,eq_class,coverage,mismatches,is_trimmed,readlen)),"from R1".to_string()))	
					} else if eq_class.len() > eq_class2.len() {
                                //println!(" 3 choose R2 -134 {:?} vs {:?}",eq_class2,eq_class);
					// non zero and different, R2 shortest
				Some((Some((mapped2,unique2,transcript2,eq_class2,coverage2,mismatches2,is_trimmed2,readlen2)),"from R2".to_string()))	

					} else {
                                //println!("4 choose R1 -134 {:?} vs {:?}",eq_class,eq_class2);
					// non zero and different, R1 shortest
				Some((Some((mapped,unique,transcript,eq_class,coverage,mismatches,is_trimmed,readlen)),"from R1".to_string()))	

				}

} 
				}	




				};

                            tx.send(selected_read).expect("Could not send data!");
                        }
                        None => {
                            // send None to tell receiver that the queue ended
                            tx.send(None).expect("Could not send data!");
                            break;
                        }
                    }; //end-match
                } // end loop
            }); //end-scope
        } // end-for



	} else {
let atomic_reader =         Arc::new(Mutex::new(reader.records()));
info!("Spawning {} threads for Mapping.\n", num_threads);
        for _ in 0..num_threads {
            let tx = tx.clone();
            let reader = Arc::clone(&atomic_reader);

            scope.spawn(move |_| {
                loop {
                    // If work is available, do that work.
                    match utils::get_next_record(&reader) {
                        Some(result_record) => {
                            let record = match result_record {
                                Ok(record) => record,
                                Err(err) => panic!("Error {:?} in reading fastq", err),
                            };

				   if trim && trimsize > &record.seq().len()/2 {

                                                //panic!("trim too long");std::process::exit(1)
                                                error!{"Trimsize {:?} too long, entire sequence trimmed for read {:?}",trimsize, &record.id()};
                                                tx.send(None).expect("Could not send data!");
                                                break;
                                        }



                                //match_strands


                            let compared_read_data = match_strands(&record,trim,trimsize,mismatchsize,index);

                            tx.send(compared_read_data).expect("Could not send data!");
                        }
                        None => {
                            // send None to tell receiver that the queue ended
                            tx.send(None).expect("Could not send data!");
                            break;
                        }
                    }; //end-match
                } // end loop
            }); //end-scope
        } // end-for


};




        let mut read_counter: usize = 0;
        let mut mapped_read_counter: usize = 0;
        let mut trimmed_read_counter: usize = 0;
       
        let mut dead_thread_count = 0;


	let mut avg_read_length = read_length_arg;

         // collect frequency of reads that match a unique ec
        let mut frequency: HashMap<&str, u32> = HashMap::new();
        let mut strandfrequency: HashMap<String, u32> = HashMap::new();




        for eq_class in rx.iter() {
            match eq_class {
                None => {
                    dead_thread_count += 1;
                    if dead_thread_count == num_threads {
                        drop(tx);
                        // can't continue with a flag check
                        // weird Rusty way !
                        // Consume whatever is remaining
                        // Not worrying about counters; hunch is their
                        // should be less
                        for eq_class in rx.iter() {
                            eq_class.map_or((), |eq_class| eprintln!("{:?}", eq_class));
                        }
                        break;
                    }
                }
                Some((None,_)) => { },
                Some((Some(read_data),strand)) => {
                    //println!("{:?}", read_data);
                    //println!("Strand:{:?}", strand);
                   *strandfrequency.entry(strand.to_string()).or_insert(0) += 1;

                    if read_data.0 {
                        mapped_read_counter += 1;

                         // read mapped to unqiue ec
                         if read_data.1 {
                   *frequency.entry(&index.tx_names[read_data.3[0] as usize]).or_insert(0) += 1;




                    } else {
			if read_data.6 { // was trimmed
                        trimmed_read_counter += 1;

			}
			}
                     } 

                    read_counter += 1;

		// read length check
		if check_read {avg_read_length += read_data.7};


                    if read_counter % 1_000_000 == 0 {
                        let frac_mapped = mapped_read_counter as f32 * 100.0 / read_counter as f32;
                        eprint!(
                            "\rDone Mapping {} reads w/ Rate: {}",
                            read_counter, frac_mapped
                        );
                        io::stderr().flush().expect("Could not flush stdout");
                    }
                } // end-Some
            } // end-match
        } // end-for


    // calculates states
    let mut unique_counter: u32 = 0;
        for (key, value) in &frequency {
            unique_counter = unique_counter + value
        }
    info!("Processed reads: {}", read_counter);
    info!("Unique  reads: {}", unique_counter);
    if trim {info!("Unique rejected reads: {}", trimmed_read_counter)} else { info!("Unique rejected reads: not run") };
    info!("Shared reads: {}", mapped_read_counter-(unique_counter as usize)-trimmed_read_counter);
    info!("Mapped reads: {}", mapped_read_counter);
    info!("Unmapped reads: {}", read_counter-mapped_read_counter);
    // calculates read strand sates
        for (key, value) in &strandfrequency {
    info!("Mapped {} reads: {}", key,value );
        }

	let average_read_length = if check_read {avg_read_length / read_counter } else {read_length_arg} ;


		// relative per base pair value
		let read_mult = if is_paired { 2*average_read_length} else {average_read_length}; 



		//manipulate the hash - fill missing entries with 0 from index

			for trans in &index.tx_names {
			if trans.contains("_"){

			if !frequency.contains_key(&trans as &str) {
			//println!(" {} is missing", trans);
  			frequency.insert(&trans as &str, 0);
			}		

		}

}
     // can add in the tx_gene map entry, and gene length entry here
    writeln!(output_file,"Gene,Deletion,Count,Total,GeneLength,ReadLength,ScaleFactor,Proportion,ScaledProportion");
        //for (key, value) in &frequency.keys().sorted() {
        for key in &frequency.keys().sorted() {
	let value = &frequency[&key as &str];
	let gene_id = &index.tx_gene_mapping.get(&key.to_string()).unwrap();
	let gene_length = &index.gene_length_mapping.get(&gene_id.to_string()).unwrap();
	 	let prop = *value as f32 / mapped_read_counter as f32;	
		let scalefactor =  (read_mult as f32 / (**gene_length as f32));
		let scaleprop = prop / scalefactor;


            writeln!(output_file,"{} ,  {} , {},  {}, {}, {}, {}, {}, {}",gene_id, key, value,mapped_read_counter , gene_length,average_read_length,scalefactor,prop, scaleprop);
        }

        //println!("gene map {:?}",&index.tx_gene_mapping);
    })
    .unwrap(); //end crossbeam

    eprintln!();

    info!("Done Mapping Reads");
    Ok(())
}




#[cfg(test)]
mod test {
    use super::*;
    use proptest::collection::vec;
    use proptest::prelude::*;
    use proptest::proptest;
    use std::collections::HashSet;
    use std::hash::Hash;
    use std::iter::FromIterator;

    fn test_intersect<T: Hash + Eq + Clone + Ord + Debug>(v1: &Vec<T>, v2: &Vec<T>) {
        let mut c1 = v1.clone();
        let c2 = v2.clone();

        let s1: HashSet<T> = HashSet::from_iter(c1.iter().cloned());
        let s2: HashSet<T> = HashSet::from_iter(c2.iter().cloned());
        let intersection = s1.intersection(&s2);

        let mut int1: Vec<T> = intersection.cloned().collect();
        int1.sort();

        intersect(&mut c1, &c2);

        assert_eq!(c1, int1);
    }

    #[test]
    fn intersect_test() {
        let v1 = vec![1, 2, 3, 4, 5, 6, 7, 8, 9];
        let v2 = vec![1, 2, 3];
        let v3 = vec![1, 4, 5];
        let v4 = vec![7, 8, 9];
        let v5 = vec![9];
        let v6: Vec<usize> = vec![];
        let v7 = vec![1, 2, 3, 6, 7, 8, 9];
        let v8 = vec![1, 7, 8, 9, 10];
        let v9 = vec![10, 15, 20];
        let v10 = vec![21, 22, 23];
        let v11 = vec![0];
        let v12 = vec![0, 1000, 5000];
        let v13 = vec![0, 1000, 1000001];
        let v14 = vec![5];
        let v15 = vec![100000000];
        let v16 = vec![1, 23, 45, 1000001, 100000000];

        let vecs = vec![
            v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11, v12, v13, v14, v15, v16,
        ];

        for v1 in vecs.iter() {
            for v2 in vecs.iter() {
                test_intersect(v1, v2);
                test_intersect(v2, v1);
            }
        }
    }

    proptest! {
        #![proptest_config(ProptestConfig { cases: 1000, .. ProptestConfig::default()})]
        #[test]
        fn intersect_prop_test(
            mut v1 in vec(0..100usize, 0..5000usize),
            mut v2 in vec(0..100usize, 0..5000usize),
        ) {

            v1.sort(); v1.dedup();
            v2.sort(); v2.dedup();
            test_intersect(&v1, &v2);
            test_intersect(&v2, &v1);
        }
    }
}
