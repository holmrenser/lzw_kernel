use std::collections::HashSet;
use std::str::FromStr;
use std::error;
use std::fmt;
use rand::{thread_rng, Rng};
use rand::distributions::Alphanumeric;

enum ALPHABET {
    FULL,
    DNA,
    PROTEIN,
    TEST
}

static FULL_ALPHABET: &str = "ABCDEFGHIJKLMNOPQRSTUVWabcdefghijklmnopqrstuvwxyz";
static PROTEIN_ALPHABET: &str = "QIYDGIHS*";
static DNA_ALPHABET: &str = "ACGTN";
static TEST_ALPHABET: &str = "ab";

#[derive(Debug, Clone)]
struct CodeDictionaryError;

impl fmt::Display for CodeDictionaryError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Cannot initialize CodeDictionary from string")
    }
}

impl error::Error for CodeDictionaryError {
    fn source(&self) -> Option<&(dyn error::Error + 'static)> {
        // Generic error, underlying cause isn't tracked.
        None
    }
}

#[derive(Debug)]
struct CodeDictionary {
    keys: HashSet<String>,
    code: Vec<String>
}

impl FromStr for CodeDictionary {
    type Err = CodeDictionaryError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let keys: HashSet<String> = s.chars()
            .map(|c| c.to_string())
            .collect::<HashSet<_>>();
        let code: Vec<String> = vec![];
        Ok( CodeDictionary { keys: keys, code: code })
    }
}

type LzwResult = std::result::Result<CodeDictionary, CodeDictionaryError>;

fn lzw_compress(s: &str, d: ALPHABET) -> LzwResult {
    let alphabet: &str;
    match d {
        ALPHABET::FULL => {
            alphabet = FULL_ALPHABET
        },
        ALPHABET::DNA => {
            alphabet = DNA_ALPHABET
        },
        ALPHABET::PROTEIN => {
            alphabet = PROTEIN_ALPHABET
        },
        ALPHABET::TEST => {
            alphabet = TEST_ALPHABET
        }
    }
    let mut dictionary = CodeDictionary::from_str(alphabet)?;
    let mut i: usize = 0;
    let max_i: usize = s.chars().count();
    let mut reached_end: bool = false;
    loop {
        let mut lookahead: usize = 0;
        loop {
            lookahead += 1;
            if i + lookahead >= max_i {
                reached_end = true;
                break
            };
            let found = dictionary.keys.contains(&s[i..i + lookahead]);
            if !found {
                break
            }
        }
        if reached_end {
            lookahead = max_i - i;
            dictionary.keys.insert(
                s[i..i + lookahead].to_string()
            );
        } else {
            lookahead -= 1;
            dictionary.keys.insert(
                s[i..i + lookahead + 1].to_string()
            );
        }
        dictionary.code.push(
            s[i..i + lookahead].to_string()
        );
        
        i += lookahead;
        if i + lookahead > max_i { break };
    }
    Ok(dictionary)
}

fn raw_lzw_kernel(
    code_a: &CodeDictionary, 
    code_b: &CodeDictionary,
    w: f64,
    gamma: f64
) -> f64 {
    let mut matches: f64 = 0.0;
    for a in &code_a.code {
        for b in &code_b.code {
            if a == b {
                matches += w
            }
        }
    }
    (gamma * matches).exp()
}

fn _normalized_lzw_kernel_reference(
    code_a: &CodeDictionary,
    code_b: &CodeDictionary,
    w: f64,
    gamma: f64
) -> f64 {
    let score_aa = raw_lzw_kernel(&code_a, &code_a, w, gamma);
    let score_bb = raw_lzw_kernel(&code_b, &code_b, w, gamma);
    let score_ab = raw_lzw_kernel(&code_a, &code_b, w, gamma);
    score_ab / (score_aa * score_bb).sqrt()
}

fn normalized_lzw_kernel(
    code_a: &CodeDictionary,
    code_b: &CodeDictionary,
    w: f64,
    gamma: f64
) -> f64 {
    let mut matches = 0.0_f64;
    for a in &code_a.code {
        for b in &code_b.code {
            if a == b {
                matches += w
            }
        }
    }
    let a_score = code_a.code.len() as f64 * w;
    let b_score = code_b.code.len() as f64 * w;
    ((gamma * matches) - (0.5_f64 * gamma * (a_score + b_score))).exp()
}

fn test() {
    let str_a = &"ACGTACGT";
    let str_b = &"ACGTACGTT";
    let code_a = lzw_compress(str_a, ALPHABET::DNA).unwrap();
    let code_b = lzw_compress(str_b, ALPHABET::DNA).unwrap();
    println!("{:?}", code_a);
    println!("{:?}", code_b);

    let w = 1.0_f64;
    let gamma = 1.0_f64;

    let raw_score = raw_lzw_kernel(&code_a, &code_b, w, gamma);
    println!("raw score: {:?}", raw_score);
    
    let normalized_score_1 = _normalized_lzw_kernel_reference(
        &code_a, &code_b, w, gamma);
    println!("normalized score 1: {:?}", normalized_score_1);
    
    let normalized_score_2 = normalized_lzw_kernel(
        &code_a, &code_b, w, gamma);
    println!("normalized score 2: {:?}", normalized_score_2);

    let score = normalized_lzw_kernel(&code_a, &code_a, w, gamma);
    println!("Assert equal score: {:?}", score);
}

fn random_string() -> String {
    let s = thread_rng()
        .sample_iter(&Alphanumeric)
        .take(100)
        .filter(|c| c.is_alphabetic())
        .collect::<String>();
    s
} 

fn benchmark(n: i64) {
    let stepsize = n / 10;
    for i in 1..n {
        if i % stepsize == 0 {
            println!("{}", i)
        }
        let s_i = random_string();
        let code_i = lzw_compress(&s_i, ALPHABET::DNA).unwrap();
        for j in 1..n {
            let s_j = random_string();
            let code_j = lzw_compress(&s_j, ALPHABET::DNA).unwrap();
            let score = normalized_lzw_kernel(&code_i, &code_j, 1.0, 1.0);
        }
    }
}

fn main() {
    benchmark(1_000_i64)
}
