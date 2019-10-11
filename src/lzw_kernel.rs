use std::collections::HashSet;
use std::str::FromStr;
use std::error;
use std::fmt;

use ndarray::Array2;
use ndarray::parallel::prelude::*;

use bio::alphabets::{
    Alphabet, 
    dna, 
    protein
};
use bio::io::fasta;

enum ALPHABET {
    Dna,
    DnaN,
    DnaIupac,
    Protein,
    Test
}

impl ALPHABET {
    fn get_alphabet(&self) -> Alphabet {
        match self {
            ALPHABET::Dna => dna::alphabet(),
            ALPHABET::DnaN => dna::n_alphabet(),
            ALPHABET::DnaIupac => dna::iupac_alphabet(),
            ALPHABET::Protein => protein::alphabet(),
            ALPHABET::Test => Alphabet::new(&b"ab"[..]),
        }
    }
}

trait RecordExt {
    fn lzw_compress(&self, a: ALPHABET) -> LzwResult;
}

impl RecordExt for fasta::Record {
    fn lzw_compress(&self, a: ALPHABET) -> LzwResult {
        let seq = self.seq()
            .into_iter()
            .map(|c| *c as char)
            .collect::<String>();
        lzw_compress(&seq, a)
    }
}

#[derive(Debug, Clone)]
struct CodeDictionaryParseError;

impl fmt::Display for CodeDictionaryParseError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Error initializing CodeDictionary")
    }
}

impl error::Error for CodeDictionaryParseError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        // Generic error, underlying cause isn't tracked.
        None
    }
}

#[derive(Debug)]
struct CodeDictionary {
    keys: HashSet<String>,
    pub code: HashSet<String>
}

impl FromStr for CodeDictionary {
    type Err = CodeDictionaryParseError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let keys: HashSet<String> = s.chars().map(|c| c.to_string()).collect();
        Ok(CodeDictionary {
            keys,
            code: HashSet::new(),
        })
    }
}

impl CodeDictionary {
    fn from_alphabet(alphabet: Alphabet) -> LzwResult {
        let keys = alphabet.symbols
            .iter()
            .map(|c| (c as u8 as char).to_string())
            .collect::<HashSet<String>>();
        Ok(CodeDictionary {
            keys: keys,
            code: HashSet::new()
        })
    }
}

type LzwResult = Result<CodeDictionary, CodeDictionaryParseError>;

fn lzw_compress(seq: &str, d: ALPHABET) -> LzwResult {
    let alphabet = d.get_alphabet();
    let mut dictionary = CodeDictionary::from_alphabet(alphabet)?;
    let mut code = String::new();
    for c in seq.chars().take(seq.len() - 1) {
        let codec = format!("{}{}", code, c);
        if dictionary.keys.contains(&codec) {
            code.push(c);
        } else {
            dictionary.keys.insert(codec);
            dictionary.code.insert(code);
            code = c.to_string();
        }
    }
    code.push(seq.chars().last().unwrap());
    dictionary.code.insert(code);
    Ok(dictionary)
}

fn _raw_lzw_score(
    code_a: &CodeDictionary, 
    code_b: &CodeDictionary,
    w: f64,
    gamma: f64
) -> f64 {
    (gamma * w * code_a.code.intersection(&code_b.code).count() as f64).exp()
}

fn _normalized_lzw_score_reference(
    code_a: &CodeDictionary,
    code_b: &CodeDictionary,
    w: f64,
    gamma: f64
) -> f64 {
    let score_aa = _raw_lzw_score(&code_a, &code_a, w, gamma);
    let score_bb = _raw_lzw_score(&code_b, &code_b, w, gamma);
    let score_ab = _raw_lzw_score(&code_a, &code_b, w, gamma);
    score_ab / (score_aa * score_bb).sqrt()
}

fn normalized_lzw_score(
    code_a: &CodeDictionary,
    code_b: &CodeDictionary,
    w: f64,
    gamma: f64
) -> f64 {
    let matches = w * code_a.code.intersection(&code_b.code).count() as f64;
    let a_score = code_a.code.len() as f64 * w;
    let b_score = code_b.code.len() as f64 * w;
    ((gamma * matches) - (0.5_f64 * gamma * (a_score + b_score))).exp()
}

pub fn lzw_kernel(fasta_records: Vec<fasta::Record>) -> Array2<f64> {
    let codes = fasta_records.into_iter()
        .map(|r| r.lzw_compress(ALPHABET::DnaIupac)
            .expect("error during lzw compression"))
        .collect::<Vec<_>>();

    let w = 1.0;
    let gamma = 0.05;
    let n = codes.len();

    let mut scores: Vec<f64> = vec![0.0; n.pow(2)];
    
    scores.par_iter_mut()
        .zip(0..n.pow(2))
        .for_each(|(score, i)| {
            let row = i / n;
            let col = i % n;
            if row == col { *score = 0.5 };
            if row < col {
                *score = normalized_lzw_score(
                    &codes[row], &codes[col], w, gamma
                )
            }
        });
    
    let mut kernel = Array2::from_shape_vec((n, n), scores).unwrap();
    kernel = &kernel + &kernel.t();
    return kernel;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let str_a = &"ACGTACGT";
        let str_b = &"ACGTACGTT";
        let code_a = lzw_compress(str_a, ALPHABET::Dna).unwrap();
        let code_b = lzw_compress(str_b, ALPHABET::Dna).unwrap();
        println!("{:?}", code_a);
        println!("{:?}", code_b);

        let w = 1.0_f64;
        let gamma = 1.0_f64;

        let raw_score = _raw_lzw_score(&code_a, &code_b, w, gamma);
        println!("raw score: {:?}", raw_score);
        
        let normalized_score_1 = _normalized_lzw_score_reference(
            &code_a, &code_b, w, gamma);
        println!("normalized score 1: {:?}", normalized_score_1);
        
        let normalized_score_2 = normalized_lzw_score(
            &code_a, &code_b, w, gamma);
        println!("normalized score 2: {:?}", normalized_score_2);

        let score = normalized_lzw_score(&code_a, &code_a, w, gamma);
        println!("Assert equal score: {:?}", score);
    }
}