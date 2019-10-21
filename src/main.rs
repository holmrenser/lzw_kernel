use std::fs::File;

use bio::io::fasta;

use clap::{Arg, App};

// use ndarray::Array2;
// use ndarray::prelude::*;
use ndarray_linalg::*;

mod lzw_kernel;
use crate::lzw_kernel::lzw_kernel;

mod neighbor_joining;
use crate::neighbor_joining::{
    NeighborJoining,
    neighbor_joining
};

mod tree;
use crate::tree::{Tree};

fn main() {
    let arguments = App::new("lzw_kernel")
        .version("0.0.1")
        .about("All vs All LZW kernel calculation for DNA or protein sequences")
        .arg(Arg::with_name("file")
                .short("f")
                .long("fasta_file")
                .takes_value(true)
                .required(true)
                .help("FASTA formatted file with DNA or protein sequences"))
        .get_matches();
    
    let fasta_file = arguments.value_of("file").unwrap();
    println!("{}", fasta_file);
    let file_handle = File::open(fasta_file).unwrap();
    let reader = fasta::Reader::new(file_handle);
    
    let records: Vec<fasta::Record> = reader.records().map(|r| {
        r.expect("error parsing fasta record")
    }).collect();

    let names: Vec<String> = records.iter()
        .map(|r| r.id().to_string())
        .collect();

    // println!("{:?}", names);

    let kernel = lzw_kernel(records);

    // let dist = 1.0 - kernel;
    let dist = kernel.inv().unwrap();

    let nj_tree = Tree::from_neighbor_joining(dist, names);

    println!("{:?}", nj_tree.to_newick());
    /*
    let d = array![
        [0., 1., 2., 3., 4., 5.],
        [1., 0., 3., 4., 5., 6.],
        [2., 3., 0., 5., 6., 7.],
        [3., 4., 5., 0., 7., 8.],
        [4., 5., 6., 7., 0., 9.],
        [5., 6., 7., 8., 9., 0.]
    ];
    let tree = neighbor_joining(d);
    */
}
