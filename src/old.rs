use std::collections::HashSet;
use std::str::FromStr;
use std::error;
use std::fmt;
/*
#[derive(Debug,Clone)]
struct ParseDictionaryError;

impl fmt::Display for ParseDictionaryError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "invalid dictionary string")
    }
}

#[derive(Debug)]
struct Dictionary<'d>  {
    values: HashSet<&'d str>
}

impl FromStr for Dictionary<'_> {
    type Err = ParseDictionaryError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let values = HashSet(["aa","ab"]);
        Ok(Dictionary { values: values })
    }
}

fn lzw_compress(string: &String) -> i64 {
    let mut count: i64 = 0;
    let mut dictionary: HashSet<String> = HashSet::new();
    for i in 1..string.chars().count() {
        println!("{}", &string[i-1..i]);
        dictionary.insert(string[..i].to_string());
        count += 1;
    }
    println!("{:?}",dictionary);
    return count
}
*/

static ALPHABET: [u8; 4] = "abcd";
fn main() {
    // println!("{}", ALPHABET)
    for c in &ALPHABET {
        println!("{}", c)
    }
}
