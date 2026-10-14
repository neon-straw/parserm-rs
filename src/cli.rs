//! Command line parsing compatible with parseRM.pl's Getopt::Long setup:
//! single or double dash, unique-prefix abbreviations, case-insensitive
//! names, `-opt value` and `-opt=value`, unknown options warned and ignored.

#[derive(Debug, Default, Clone)]
pub struct Opts {
    pub input: Option<String>,
    pub dir: bool,
    pub parse: bool,
    pub age: Option<String>,
    pub land: Option<String>,
    pub my: Option<String>,
    pub fa: bool,
    pub nrem: bool,
    pub glen: Option<String>,
    pub rlib: Option<String>,
    pub kim: bool,
    pub simple: Option<String>,
    pub edit: Option<String>,
    pub what: Option<String>,
    pub contain: bool,
    pub updates: bool,
    pub help: bool,
    pub v: bool,
    /// Not in parseRM.pl: opt into fixes of the original's bugs (no parity).
    pub fix_quirks: bool,
}

const NAMES: &[(&str, bool)] = &[
    ("in", true),
    ("dir", false),
    ("parse", false),
    ("age", true),
    ("land", true),
    ("my", true),
    ("fa", false),
    ("nrem", false),
    ("glen", true),
    ("rlib", true),
    ("kim", false),
    ("simple", true),
    ("edit", true),
    ("what", true),
    ("contain", false),
    ("updates", false),
    ("help", false),
    ("v", false),
];

fn resolve(name: &str) -> Result<(&'static str, bool), String> {
    let lname = name.to_ascii_lowercase();
    if let Some(&(n, v)) = NAMES.iter().find(|(n, _)| *n == lname) {
        return Ok((n, v));
    }
    let cands: Vec<_> = NAMES.iter().filter(|(n, _)| n.starts_with(&lname)).collect();
    match cands.as_slice() {
        [(n, v)] => Ok((n, *v)),
        [] => Err(format!("Unknown option: {name}")),
        _ => Err(format!("Option {name} is ambiguous")),
    }
}

pub fn parse_args(args: &[String]) -> Opts {
    let mut o = Opts::default();
    let mut i = 0;
    while i < args.len() {
        let a = &args[i];
        i += 1;
        if a == "--" {
            break;
        }
        let body = match a.strip_prefix("--").or_else(|| a.strip_prefix('-')) {
            Some(b) if !b.is_empty() => b,
            _ => continue, // non-option argument: left alone, like Getopt::Long's permute
        };
        let (name, inline) = match body.split_once('=') {
            Some((n, v)) => (n, Some(v.to_string())),
            None => (body, None),
        };
        if name.eq_ignore_ascii_case("fix-quirks") || name.eq_ignore_ascii_case("fix_quirks") {
            o.fix_quirks = true;
            continue;
        }
        let (opt, takes) = match resolve(name) {
            Ok(r) => r,
            Err(e) => {
                eprintln!("{e}");
                continue;
            }
        };
        let value = if takes {
            match inline {
                Some(v) => Some(v),
                None if i < args.len() => {
                    i += 1;
                    Some(args[i - 1].clone())
                }
                None => {
                    eprintln!("Option {opt} requires an argument");
                    continue;
                }
            }
        } else {
            if inline.is_some() {
                eprintln!("Option {opt} does not take an argument");
                continue;
            }
            None
        };
        match opt {
            "in" => o.input = value,
            "dir" => o.dir = true,
            "parse" => o.parse = true,
            "age" => o.age = value,
            "land" => o.land = value,
            "my" => o.my = value,
            "fa" => o.fa = true,
            "nrem" => o.nrem = true,
            "glen" => o.glen = value,
            "rlib" => o.rlib = value,
            "kim" => o.kim = true,
            "simple" => o.simple = value,
            "edit" => o.edit = value,
            "what" => o.what = value,
            "contain" => o.contain = true,
            "updates" => o.updates = true,
            "help" => o.help = true,
            "v" => o.v = true,
            _ => unreachable!(),
        }
    }
    o
}

/// Perl truthiness of an optional string option.
pub fn t(v: &Option<String>) -> bool {
    v.as_deref().is_some_and(|s| crate::perl::truthy(s.as_bytes()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(s: &str) -> Opts {
        parse_args(&s.split_whitespace().map(String::from).collect::<Vec<_>>())
    }

    #[test]
    fn abbreviations() {
        let o = p("-i x.out -p -a 5,15 -l=50,1 --kim -s all -w name,L1 -c -v -m 0.1");
        assert_eq!(o.input.as_deref(), Some("x.out"));
        assert!(o.parse && o.kim && o.contain && o.v);
        assert_eq!(o.age.as_deref(), Some("5,15"));
        assert_eq!(o.land.as_deref(), Some("50,1"));
        assert_eq!(o.simple.as_deref(), Some("all"));
        assert_eq!(o.what.as_deref(), Some("name,L1"));
        assert_eq!(o.my.as_deref(), Some("0.1"));
        let o = p("-in a.align -parse -glen g.txt -rlib l.fa -edit e.txt --fix-quirks");
        assert_eq!(o.glen.as_deref(), Some("g.txt"));
        assert_eq!(o.rlib.as_deref(), Some("l.fa"));
        assert_eq!(o.edit.as_deref(), Some("e.txt"));
        assert!(o.fix_quirks);
    }
}
