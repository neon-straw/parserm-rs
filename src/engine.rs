//! Port of parseRM.pl's main flow (`-p`, `-a`, `-l`). Line references
//! (`:NNN`) point at parseRM.pl v5.8.2.

use crate::cli::{parse_args, t, Opts};
use crate::lexmin::lexmin;
use crate::perl::{filename, fmt, get_path, has_word, is_space, lines, num, split_byte, split_ws, truthy};
use crate::sweep::{sweep, Run, SweepHit};
use rustc_hash::FxHashMap;
use std::path::Path;

/// A Perl `die`/`confess`: message to STDERR, non-zero exit.
#[derive(Debug)]
pub struct Die(pub String);

type R<T> = Result<T, Die>;

macro_rules! die {
    ($($a:tt)*) => { return Err(Die(format!($($a)*))) };
}

const VERSION: &str = "5.8.2";
const USAGE: &str = "
   Usage [v5.8.2]:
   parserm-rs -i <genome.(out|align)>
            [-p] [-g <X>] [-r <repeat_library.fa>]
            [-a <file.txt> OR <X,X>] [-m <file.txt> OR <X>] [-l <max,bin>]
            [-d] [-k] [-s <type>] [-e <file>] [-w <type,name>] [-c]
            [-v] [-u] [-h] [--fix-quirks]

   Rust port of parseRM.pl (Kapusta, Suh & Feschotte 2017, PNAS).
   Same options and outputs as parseRM.pl; see `perl parseRM.pl -h`.
   --fix-quirks enables fixes for bugs of the original (no output parity).
";

/// Age categories; index into per-type counters.
const AGE_TYPES: [&str; 5] = ["LS", "A", "nd", "na", ""];
const T_LS: u8 = 0;
const T_A: u8 = 1;
const T_ND: u8 = 2;
const T_NA: u8 = 3;
const T_UNDEF: u8 = 4;

/// Fake address printed for the `Rage` landscape row, where parseRM.pl
/// stringifies a reference (`REF(0x...)`, :883-886).
const REF_KEY: &str = "REF(0x600000000000)";

struct Filter {
    ftype: Vec<u8>,
    name_lc: Vec<u8>,
    re: Option<regex::bytes::Regex>,
}

struct Land {
    max: f64,
    /// (coord, j, printable column of `coord`)
    bins: Vec<(f64, f64, Option<usize>)>,
    /// header columns `[i;i2[`
    header: Vec<(String, String)>,
}

impl Land {
    fn new(spec: &str) -> R<Land> {
        let f = split_byte(spec.as_bytes(), b',');
        let max = num(f.first().copied().unwrap_or(b""));
        let bin = num(f.get(1).copied().unwrap_or(b""));
        if bin.is_nan() || bin <= 0.0 || !max.is_finite() {
            die!("\n --land {spec}: bin length must be > 0 (parseRM.pl would loop forever)\n\n");
        }
        // :904 / :936 (print loops)
        let mut header = Vec::new();
        let mut col_of: FxHashMap<String, usize> = FxHashMap::default();
        let mut k = 0.0f64;
        while k < max {
            let s = fmt(k);
            col_of.entry(s.clone()).or_insert(header.len());
            header.push((s, fmt(k + bin)));
            k += bin;
        }
        // :880 (bin search loop)
        let mut bins = Vec::new();
        let mut j = bin;
        while j <= max {
            let coord = j - bin;
            bins.push((coord, j, col_of.get(&fmt(coord)).copied()));
            j += bin;
        }
        Ok(Land { max, bins, header })
    }

    /// `None`: no bin matched (no landscape entry). `Some(col)`: entry created,
    /// counted in printable column `col` if any.
    fn locate(&self, div: f64) -> Option<Option<usize>> {
        if div > self.max {
            return None;
        }
        self.bins
            .iter()
            .find(|&&(c, j, _)| div >= c && div < j)
            .map(|&(_, _, col)| col)
    }
}

pub struct Ctx {
    o: Opts,
    input: String,
    nonte: Vec<u8>,
    p: bool,
    a: bool,
    teage: bool,
    my: bool,
    edit_crash: bool,
    /// `$TE`: lc(name) -> (class, family) (columns 1 and 2)
    te: FxHashMap<Vec<u8>, (Vec<u8>, Vec<u8>)>,
    genlen: FxHashMap<Vec<u8>, Vec<u8>>,
    liblen: FxHashMap<Vec<u8>, Option<Vec<u8>>>,
    srates: FxHashMap<Vec<u8>, Vec<u8>>,
    age: FxHashMap<Vec<u8>, (Vec<u8>, Vec<u8>)>,
    land: Option<Land>,
    filter: Option<Filter>,
}

fn exists(p: &str) -> bool {
    !p.is_empty() && Path::new(p).exists()
}

fn read(p: &str, ctx: &str) -> R<Vec<u8>> {
    std::fs::read(p).map_err(|e| Die(format!("\nERROR ({ctx}): could not open to read {p} {e}\n")))
}

fn write(p: &str, data: &[u8], ctx: &str) -> R<()> {
    std::fs::write(p, data).map_err(|e| Die(format!("\nERROR ({ctx}): could not open to write {p} {e}\n")))
}

fn lc(s: &[u8]) -> Vec<u8> {
    s.to_ascii_lowercase()
}

/// `$x =~ /.fa(sta|$)/`
fn looks_fasta(s: &str) -> bool {
    let b = s.as_bytes();
    (1..b.len()).any(|i| b[i..].starts_with(b"fa") && (i + 2 == b.len() || b[i + 2..].starts_with(b"sta")))
}

/// `$1 if ($x =~ /(.*)\.out/ || $x =~ /(.*)\.align/)`
fn strip_rm_ext(s: &str) -> Option<&str> {
    s.rfind(".out").or_else(|| s.rfind(".align")).map(|i| &s[..i])
}

/// Entry point; returns the process exit code.
pub fn run(args: &[String]) -> i32 {
    match run_inner(args) {
        Ok(()) => 0,
        Err(Die(msg)) => {
            eprint!("{msg}");
            255
        }
    }
}

fn run_inner(args: &[String]) -> R<()> {
    let o = parse_args(args);
    let mut ctx = Ctx::new(o)?;
    let files: Vec<String> = if ctx.o.dir {
        ls(&ctx.input)
    } else {
        vec![ctx.input.clone()]
    };
    let path = get_path(&ctx.input);
    if ctx.p {
        ctx.prep_parse(&path)?;
    }
    if ctx.a {
        ctx.prep_age()?;
    }
    if let Some(spec) = ctx.o.land.clone().filter(|s| truthy(s.as_bytes())) {
        ctx.land = Some(Land::new(&spec)?);
    }
    if ctx.my {
        let my = ctx.o.my.clone().unwrap_or_default();
        ctx.load_val(&my, ValKind::Rates)?;
    }
    if t(&ctx.o.edit) {
        let edit = ctx.o.edit.clone().unwrap_or_default();
        ctx.load_te_info(&edit)?;
        // :162 `my $TE = load_TE_info($TES)` assigns the sub's return value (1)
        // to $TE, so the first `$TE->{...}` dies under strict refs.
        ctx.edit_crash = !ctx.o.fix_quirks;
    }
    for file in files {
        if !truthy(file.as_bytes()) {
            continue;
        }
        let f = if ctx.o.dir {
            format!("{}/{}", ctx.input, filename(&file))
        } else {
            format!("{}/{}", path, filename(&file))
        };
        if ctx.o.v {
            eprintln!("     -> {f}..");
        }
        if !exists(&f) || !(f.ends_with(".align") || f.ends_with(".out")) {
            continue;
        }
        ctx.process_file(&f)?;
    }
    Ok(())
}

/// `ls $IN` (:141): the order comes from ls itself (locale collation), which
/// decides the content of `$IN.splitage_all.tab` and the first-seen `$TE` entries.
fn ls(dir: &str) -> Vec<String> {
    match std::process::Command::new("sh").arg("-c").arg(format!("ls {dir}")).output() {
        Ok(out) => String::from_utf8_lossy(&out.stdout).lines().map(String::from).collect(),
        Err(_) => list_dir(dir),
    }
}

fn list_dir(dir: &str) -> Vec<String> {
    let mut v: Vec<String> = match std::fs::read_dir(dir) {
        Ok(rd) => rd
            .filter_map(|e| e.ok())
            .map(|e| e.file_name().to_string_lossy().into_owned())
            .filter(|n| !n.starts_with('.'))
            .collect(),
        Err(_) => Vec::new(),
    };
    v.sort();
    v
}

enum ValKind {
    Rates,
    Len,
}

impl Ctx {
    /// :1058 check_opt
    fn new(o: Opts) -> R<Ctx> {
        let inp = t(&o.input);
        let p = o.parse;
        let a = t(&o.age);
        let l = t(&o.land);
        if !inp && !p && !a && !l && !o.help && !o.updates && o.v {
            die!("\n Script parseRM.pl version {VERSION}\n\n");
        }
        if o.updates || o.help {
            die!("{USAGE}");
        }
        if !inp && !p && !a && !l {
            die!("{USAGE}");
        }
        if !p && !a && !l {
            die!("\n please chose one of --parse, --age or --land (use -h to see the full usage)\n\n");
        }
        let mut input = o.input.clone().unwrap_or_default();
        if inp && !exists(&input) {
            die!("\n {input} does not exist?\n\n");
        }
        if t(&o.rlib) && !exists(o.rlib.as_deref().unwrap_or("")) {
            die!("\n {} does not exist?\n\n", o.rlib.as_deref().unwrap_or(""));
        }
        let check_file = |v: &Option<String>| -> R<()> {
            let s = v.as_deref().unwrap_or("");
            if !exists(s) && !s.bytes().any(|b| b.is_ascii_digit() || b == b'.' || b == b',') {
                die!("\n {s} does not exist as a file but also not numerical?\n\n");
            }
            Ok(())
        };
        if t(&o.glen) {
            check_file(&o.glen)?;
        }
        if a {
            check_file(&o.age)?;
        }
        let teage = a && exists(o.age.as_deref().unwrap_or(""));
        let my = t(&o.my);
        if my {
            check_file(&o.my)?;
        }
        if o.dir && input.len() > 1 && input.ends_with('/') {
            input.pop();
        }
        let filter = if t(&o.what) {
            let w = o.what.clone().unwrap_or_default();
            let f = split_byte(w.as_bytes(), b',');
            let ftype = f.first().map(|x| x.to_vec()).unwrap_or_default();
            let name_lc = lc(f.get(1).copied().unwrap_or(b""));
            let re = if o.contain {
                let pat = String::from_utf8_lossy(&name_lc).into_owned();
                Some(
                    regex::bytes::Regex::new(&pat)
                        .map_err(|e| Die(format!("\nERROR: invalid --what regexp {pat}: {e}\n")))?,
                )
            } else {
                None
            };
            Some(Filter { ftype, name_lc, re })
        } else {
            None
        };
        Ok(Ctx {
            nonte: o.simple.clone().unwrap_or_else(|| "no".into()).into_bytes(),
            input,
            p,
            a,
            teage,
            my,
            edit_crash: false,
            te: FxHashMap::default(),
            genlen: FxHashMap::default(),
            liblen: FxHashMap::default(),
            srates: FxHashMap::default(),
            age: FxHashMap::default(),
            land: None,
            filter,
            o,
        })
    }

    /// :422 load_val for "rates" and "len"
    fn load_val(&mut self, data: &str, kind: ValKind) -> R<()> {
        let map = match kind {
            ValKind::Rates => &mut self.srates,
            ValKind::Len => &mut self.genlen,
        };
        if exists(data) {
            let buf = read(data, "sub load_val")?;
            for l in lines(&buf) {
                let f = split_ws(l);
                let name = f.first().copied().unwrap_or(b"").to_vec();
                let d = f.get(1).copied().unwrap_or(b"").to_vec();
                map.insert(name, d);
            }
        } else {
            map.insert(filename(&self.input).as_bytes().to_vec(), data.as_bytes().to_vec());
        }
        Ok(())
    }

    /// :571 prep_parse
    fn prep_parse(&mut self, path: &str) -> R<()> {
        if t(&self.o.glen) {
            let g = self.o.glen.clone().unwrap_or_default();
            self.load_val(&g, ValKind::Len)?;
        } else if self.o.fa {
            self.get_tot_length(path);
        }
        if t(&self.o.rlib) {
            self.get_all_lengths()?;
        }
        Ok(())
    }

    /// :590 get_tot_length (-fa). Not covered by parity; implements the
    /// intended behaviour (N/n removed unless -nrem).
    fn get_tot_length(&mut self, path: &str) {
        let list: Vec<String> = if !self.o.dir {
            let mut fa = strip_rm_ext(&self.input).unwrap_or(&self.input).to_string();
            if !looks_fasta(&fa) {
                fa.push_str(".fa");
            }
            let p = format!("{path}/{fa}");
            if !exists(&p) {
                eprintln!("       WARN: {fa} could not be found, % genome won't be determined");
                return;
            }
            vec![p]
        } else {
            list_dir(&self.input)
                .into_iter()
                .filter(|n| n.ends_with("fa"))
                .map(|n| format!("{}/{}", self.input, n))
                .collect()
        };
        for fa in list {
            if fa == self.input || !looks_fasta(&fa) {
                continue;
            }
            let Ok(buf) = std::fs::read(&fa) else { continue };
            if buf.first() != Some(&b'>') {
                continue;
            }
            let mut len: u64 = 0;
            for l in lines(&buf) {
                if l.contains(&b'>') {
                    continue;
                }
                len += if self.o.nrem {
                    l.len() as u64
                } else {
                    l.iter().filter(|&&b| b != b'n' && b != b'N').count() as u64
                };
            }
            self.genlen.insert(filename(&fa).as_bytes().to_vec(), len.to_string().into_bytes());
        }
    }

    /// :628 get_all_lengths (-rlib); writes `$LIB.lengths` like the original.
    fn get_all_lengths(&mut self) -> R<()> {
        let lib = self.o.rlib.clone().unwrap_or_default();
        let buf = std::fs::read(&lib).unwrap_or_default();
        if buf.first() != Some(&b'>') {
            if self.o.v {
                eprintln!("       WARN: {lib} not a fasta file? => skipped");
            }
            return Ok(());
        }
        let lfile = format!("{lib}.lengths");
        if exists(&lfile) {
            if let Ok(lbuf) = std::fs::read(&lfile) {
                for l in lines(&lbuf) {
                    let f = split_ws(l);
                    let id = f.first().copied().unwrap_or(b"");
                    let rn = split_byte(id, b'#').first().copied().unwrap_or(b"");
                    self.liblen.insert(lc(rn), f.get(1).map(|x| x.to_vec()));
                }
                return Ok(());
            }
        }
        let mut out = Vec::new();
        let mut id: Vec<u8> = Vec::new();
        let mut ln: u64 = 0;
        let mut started = false;
        let record = |id: &[u8], ln: u64, out: &mut Vec<u8>, liblen: &mut FxHashMap<_, _>| {
            out.extend_from_slice(id);
            out.extend_from_slice(format!("\t{ln}\n").as_bytes());
            let rn = split_byte(id, b'#').first().copied().unwrap_or(b"");
            liblen.insert(lc(rn), Some(ln.to_string().into_bytes()));
        };
        for l in lines(&buf) {
            if l.first() == Some(&b'>') {
                if started {
                    record(&id, ln, &mut out, &mut self.liblen);
                }
                started = true;
                let f = split_ws(l);
                id = f.first().copied().unwrap_or(b"").to_vec();
                if let Some(pos) = id.iter().position(|&b| b == b'>') {
                    id.remove(pos);
                }
                ln = 0;
            } else {
                ln += l.len() as u64;
            }
        }
        record(&id, ln, &mut out, &mut self.liblen);
        let _ = std::fs::write(&lfile, &out);
        Ok(())
    }

    /// :948 prep_age
    fn prep_age(&mut self) -> R<()> {
        let aa = self.o.age.clone().unwrap_or_default();
        if self.teage {
            // NB: loads into $TE, not $AGE (:951); age types end up undef.
            self.load_te_info(&aa)
        } else {
            // :451 load_val_hash
            let d = aa.as_bytes();
            let (m, n): (Vec<u8>, Vec<u8>) = if d.contains(&b',') {
                let f = split_byte(d, b',');
                (
                    f.first().copied().unwrap_or(b"").to_vec(),
                    f.get(1).copied().unwrap_or(b"").to_vec(),
                )
            } else {
                (d.to_vec(), d.to_vec())
            };
            let v = if num(&m) < num(&n) { (m, n) } else { (n, m) };
            self.age.insert(filename(&self.input).as_bytes().to_vec(), v);
            Ok(())
        }
    }

    /// :960 load_TE_info
    fn load_te_info(&mut self, file: &str) -> R<()> {
        let buf = read(file, "sub get_TEs_infos")?;
        for l in lines(&buf) {
            if !has_word(l) {
                continue;
            }
            let f = split_byte(l, b'\t');
            let get = |i: usize| f.get(i).copied().unwrap_or(b"").to_vec();
            self.te.insert(lc(&get(0)), (get(1), get(2)));
        }
        Ok(())
    }

    /// :463 get_Rclass_Rfam (returns corrected class and family).
    fn get_rclass_rfam(&mut self, rname: &[u8], classfam: &[u8]) -> R<(Vec<u8>, Vec<u8>)> {
        if self.edit_crash {
            die!("Can't use string (\"1\") as a HASH ref while \"strict refs\" in use at parseRM.pl line 478.\n");
        }
        let key = lc(rname);
        if let Some((c, f)) = self.te.get(&key) {
            if truthy(c) {
                return Ok((c.clone(), f.clone()));
            }
        }
        let (rc, rf): (Vec<u8>, Vec<u8>) = if classfam.contains(&b'/') {
            let f = split_byte(classfam, b'/');
            (
                f.first().copied().unwrap_or(b"").to_vec(),
                f.get(1).copied().unwrap_or(b"").to_vec(),
            )
        } else {
            match classfam.iter().rposition(|&b| b == b'.') {
                Some(i) => (classfam[i + 1..].to_vec(), classfam[..i].to_vec()),
                None => (classfam.to_vec(), classfam.to_vec()),
            }
        };
        self.te.insert(key, (rc.clone(), rf.clone()));
        Ok((rc, rf))
    }

    /// :500 get_RM_array_skip
    fn skip(&self, rname: &[u8], rc: &[u8], rf: &[u8]) -> bool {
        let nonte = self.nonte.as_slice();
        let mut skip = false;
        if !(nonte == b"nonTE" || nonte == b"all") && rc == b"nonTE" {
            skip = true;
        }
        if nonte != b"all"
            && (rc == b"Simple_repeat"
                || rc == b"Low_complexity"
                || rc == b"Satellite"
                || rc.ends_with(b"RNA")
                || rc.ends_with(b"omeric")
                || rc == b"ARTEFACT")
        {
            skip = true;
        }
        if let Some(fl) = &self.filter {
            let ft = fl.ftype.as_slice();
            let keep = match &fl.re {
                Some(re) => {
                    (ft == b"name" && re.is_match(&lc(rname)))
                        || (ft == b"class" && re.is_match(&lc(rc)))
                        || (ft == b"family" && re.is_match(&lc(rf)))
                }
                None => {
                    (ft == b"name" && fl.name_lc == lc(rname))
                        || (ft == b"class" && fl.name_lc == lc(rc))
                        || (ft == b"family" && fl.name_lc == lc(rf))
                }
            };
            if !keep {
                skip = true;
            }
        }
        skip
    }

    fn process_file(&mut self, f: &str) -> R<()> {
        if Path::new(f).is_dir() {
            // open() on a directory succeeds in Perl, readline gives undef: no repeats.
            return Ok(());
        }
        let file = std::fs::File::open(f)
            .map_err(|e| Die(format!("\nERROR (sub get_RMalign_array): can't open to read {f} {e}\n")))?;
        let mut data = self.load(f, std::io::BufReader::with_capacity(1 << 20, file))?;
        if data.hits.is_empty() {
            return Ok(());
        }
        self.apply_my(f, &mut data)?;
        let agg = self.aggregate(&data);
        self.print_outputs(f, &data, &agg)
    }
}

// ---------------------------------------------------------------------------
// Loading
// ---------------------------------------------------------------------------

#[derive(Clone, Copy)]
struct Hit {
    start: i64,
    end: i64,
    div: f64,
    div_truthy: bool,
    div_s: u32,
    del_s: u32,
    ins_s: u32,
    name: u32,
    gname: u32,
    rfull: u32,
    rid: u32,
}

#[derive(Default)]
struct Interner {
    map: FxHashMap<Vec<u8>, u32>,
    strs: Vec<Vec<u8>>,
}

impl Interner {
    fn id(&mut self, s: &[u8]) -> u32 {
        if let Some(&i) = self.map.get(s) {
            return i;
        }
        let i = self.strs.len() as u32;
        self.map.insert(s.to_vec(), i);
        self.strs.push(s.to_vec());
        i
    }
    fn len(&self) -> usize {
        self.strs.len()
    }
}

/// String keys mapped to dense ids; the strings are never read back.
#[derive(Default)]
struct KeyIds {
    map: FxHashMap<Box<[u8]>, u32>,
}

impl KeyIds {
    fn id(&mut self, s: &[u8]) -> u32 {
        if let Some(&i) = self.map.get(s) {
            return i;
        }
        let i = self.map.len() as u32;
        self.map.insert(s.into(), i);
        i
    }
    fn len(&self) -> usize {
        self.map.len()
    }
}

#[derive(Default)]
struct FileData {
    hits: Vec<Hit>,
    names: Interner,
    gnames: Interner,
    vals: Interner,
    rfull: KeyIds,
    rid: KeyIds,
    is_align: bool,
    path: String,
}

/// Coordinates beyond 2^53 cannot be represented exactly (Perl would also
/// loop practically forever over them).
const MAX_COORD: f64 = 9_007_199_254_740_992.0;

/// Repeat name as recovered at aggregation time: `split("-#-", $RFULLID)`
/// keeps the part before the first `-#-`, which may cut names containing it.
fn parse_name(rname: &[u8]) -> &[u8] {
    // In `$RNAME . "-#-" . ...` the first `-#-` either lies inside the name or
    // straddles it only when the name ends with "-#".
    match memchr::memmem::find(rname, b"-#-") {
        Some(i) => &rname[..i],
        None if rname.ends_with(b"-#") => &rname[..rname.len() - 2],
        None => rname,
    }
}

/// Div value as loaded: a string from the file, or an NV after `-kim`.
enum Div<'a> {
    Str(&'a [u8]),
    Num(f64),
}

impl Ctx {
    /// :254 load_RM_in_array
    fn load(&mut self, f: &str, mut rdr: impl std::io::BufRead) -> R<FileData> {
        let is_align = f.contains(".align");
        let mut d = FileData {
            is_align,
            path: f.to_string(),
            ..Default::default()
        };
        let filter_on = !(self.nonte == b"all" && self.filter.is_none());
        let mut key: Vec<u8> = Vec::with_capacity(128);
        let mut prevskip_no = false;
        let mut seen_header = false;
        let mut skipped_noblock = 0usize;

        let mut line: Vec<u8> = Vec::with_capacity(256);
        loop {
            line.clear();
            // A read error ends the file, like Perl's readline returning undef.
            let n = rdr.read_until(b'\n', &mut line).unwrap_or(0);
            if n == 0 {
                break;
            }
            if line.last() == Some(&b'\n') {
                line.pop();
            }
            let raw = line.as_slice();
            // :488 check_line
            if is_align {
                if !has_word(raw)
                    || is_space(raw[0])
                    || raw.starts_with(b"C ")
                    || raw.starts_with(b"Matrix")
                    || raw.starts_with(b"Transitions")
                    || raw.starts_with(b"Gap_init")
                {
                    continue;
                }
            } else if !has_word(raw)
                || memchr::memmem::find(raw, b"position").is_some()
                || memchr::memmem::find(raw, b"score").is_some()
                || memchr::memmem::find(raw, b"Score").is_some()
                || memchr::memmem::find(raw, b"Gname").is_some()
            {
                continue;
            }

            let (l, fields): (&[u8], Vec<&[u8]>) = if is_align {
                (raw, if raw[0].is_ascii_digit() { split_ws(raw) } else { Vec::new() })
            } else {
                let st = raw.iter().position(|&b| !is_space(b)).unwrap_or(raw.len());
                let l = &raw[st..];
                (l, split_ws(l))
            };
            let fld = |i: usize| fields.get(i).copied().unwrap_or(b"");

            let (rname, classfam, div, strand, block): (&[u8], &[u8], Div, &[u8], &[u8]);
            if is_align {
                // :304 get_RMalign_val
                if !l[0].is_ascii_digit() {
                    // Kimura and other lines re-run class/filter logic on stale
                    // values, which is idempotent once a header was seen. The
                    // Kimura branch itself is dead code (`substr($l,...)` on a ref).
                    if self.edit_crash || !seen_header {
                        self.get_rclass_rfam(b"", b"")?;
                    }
                    if self.o.fix_quirks && prevskip_no && l.starts_with(b"Kimura") {
                        if let (Some(last), Some(eq)) = (d.hits.last_mut(), l.iter().position(|&b| b == b'=')) {
                            let v: Vec<u8> = l[eq + 1..].iter().copied().filter(|&b| !is_space(b)).collect();
                            let v = match split_byte(&v, b'=').first() {
                                Some(x) => x.to_vec(),
                                None => Vec::new(),
                            };
                            last.div = num(&v);
                            last.div_truthy = truthy(&v);
                            last.div_s = d.vals.id(&v);
                        }
                    }
                    continue;
                }
                seen_header = true;
                let (st, full) = if fld(8) == b"C" || fld(8) == b"+" {
                    (fld(8), fld(9))
                } else {
                    (&b"+"[..], fld(8))
                };
                strand = st;
                if full.contains(&b'#') {
                    let parts = split_byte(full, b'#');
                    rname = parts.first().copied().unwrap_or(b"");
                    classfam = parts.get(1).copied().unwrap_or(b"");
                } else {
                    rname = full;
                    classfam = b"Unknown";
                }
                div = Div::Str(fld(1));
                block = b"";
            } else {
                // :340 get_RMout_val
                rname = fld(9);
                classfam = fld(10);
                strand = fld(8);
                block = fld(14);
                div = if self.o.kim {
                    let arg = 1.0 - num(fld(1)) * 4.0 / 300.0;
                    if arg <= 0.0 {
                        die!("Can't take log of {} at parseRM.pl line 350.\n", fmt(arg));
                    }
                    Div::Num(-300.0 / 4.0 * arg.ln())
                } else {
                    Div::Str(fld(1))
                };
                if !truthy(block) {
                    skipped_noblock += 1;
                    continue;
                }
            }

            let (rc, rf) = self.get_rclass_rfam(rname, classfam)?;
            let skip = filter_on && self.skip(rname, &rc, &rf);
            prevskip_no = !skip;
            if skip {
                continue;
            }
            if !l[0].is_ascii_digit() {
                continue;
            }

            let (gname, gst, gen) = (fld(4), fld(5), fld(6));
            key.clear();
            key.extend_from_slice(rname);
            key.extend_from_slice(b"-#-");
            if is_align {
                key.extend_from_slice(gname);
                key.push(b':');
                key.extend_from_slice(gst);
                key.push(b'_');
                key.extend_from_slice(gen);
                key.push(b'_');
                key.extend_from_slice(strand);
            } else {
                key.extend_from_slice(block);
                key.push(b':');
                key.extend_from_slice(gname);
            }
            let rfull = d.rfull.id(&key);
            key.clear();
            key.extend_from_slice(rname);
            key.extend_from_slice(b"-#-");
            key.extend_from_slice(gname);
            key.push(b':');
            key.extend_from_slice(gst);
            key.push(b'-');
            key.extend_from_slice(gen);
            let rid = d.rid.id(&key);

            let (dv, dt, ds) = match div {
                Div::Str(s) => (num(s), truthy(s), d.vals.id(s)),
                Div::Num(v) => (v, v != 0.0, d.vals.id(fmt(v).as_bytes())),
            };
            let (start, end) = (num(gst), num(gen));
            if start.abs() > MAX_COORD || end.abs() > MAX_COORD || start.is_nan() || end.is_nan() {
                die!("\nERROR: coordinates out of range in {f}: {} {}\n", fmt(start), fmt(end));
            }
            d.hits.push(Hit {
                start: start as i64,
                end: end as i64,
                div: dv,
                div_truthy: dt,
                div_s: ds,
                del_s: d.vals.id(fld(2)),
                ins_s: d.vals.id(fld(3)),
                name: d.names.id(parse_name(rname)),
                gname: d.gnames.id(gname),
                rfull,
                rid,
            });
        }
        if self.o.v && skipped_noblock > 0 {
            eprintln!("          WARN: {skipped_noblock} lines skipped because they had no block info");
        }
        Ok(d)
    }

    /// :364 `-my`: div -> My, per hit, before the position loop.
    fn apply_my(&self, f: &str, d: &mut FileData) -> R<()> {
        if !self.my {
            return Ok(());
        }
        let rate = self.srates.get(filename(f).as_bytes()).map(|v| num(v)).unwrap_or(0.0);
        let denom = rate * 2.0;
        if denom == 0.0 {
            die!("Illegal division by zero at parseRM.pl line 366.\n");
        }
        for h in &mut d.hits {
            let v = h.div / 100.0 / denom;
            h.div = v;
            h.div_truthy = v != 0.0;
            h.div_s = d.vals.id(fmt(v).as_bytes());
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Aggregation (parse_all_parse / parse_all_age / parse_all_land)
// ---------------------------------------------------------------------------

#[derive(Default)]
struct Agg {
    tot: u64,
    double: u64,
    // per name id
    class_of: Vec<Vec<u8>>,
    fam_of: Vec<Vec<u8>>,
    name_nr: Vec<u64>,
    name_dbl: Vec<u64>,
    name_cnt_nr: Vec<u64>,
    name_cnt_tot: Vec<u64>,
    /// (name, value string) -> nt, for div / del / ins
    hist: [FxHashMap<(u32, u32), u64>; 3],
    /// Winning intervals in scaffold-group order, used only for ambiguous median ties.
    median_runs: Vec<Vec<MedianRun>>,
    age_nr: [u64; 5],
    age_dbl: [u64; 5],
    age_cnt_nr: [u64; 5],
    age_cnt_tot: [u64; 5],
    land_name: Vec<Option<Vec<u64>>>,
    land_rage: Option<Vec<u64>>,
}

struct MedianRun {
    group: usize,
    a: i64,
    b: i64,
    values: [u32; 3],
}

/// Fragment-ID "first seen" bookkeeping for the age counts (:1004-1012):
/// the first position in iteration order (groups in file order, positions in
/// string order within a group) decides the category.
struct FirstSeen {
    seen: Vec<bool>,
    cands: Vec<(u32, u8, i64, i64)>,
}

impl FirstSeen {
    fn new(n: usize) -> Self {
        FirstSeen {
            seen: vec![false; n],
            cands: Vec::new(),
        }
    }

    fn push(&mut self, key: u32, ty: u8, a: i64, b: i64) {
        if !self.seen[key as usize] {
            self.cands.push((key, ty, a, b));
        }
    }

    /// Resolve the candidates of one group.
    fn flush(&mut self, counts: &mut [u64; 5]) {
        let mut cands = std::mem::take(&mut self.cands);
        cands.sort_by_key(|c| c.0);
        for chunk in cands.chunk_by(|x, y| x.0 == y.0) {
            let ty = if chunk.iter().all(|c| c.1 == chunk[0].1) {
                chunk[0].1
            } else {
                chunk
                    .iter()
                    .map(|c| (lexmin(c.2, c.3), c.1))
                    .min_by(|x, y| x.0.as_bytes().cmp(y.0.as_bytes()))
                    .map(|x| x.1)
                    .expect("non-empty chunk")
            };
            counts[ty as usize] += 1;
            self.seen[chunk[0].0 as usize] = true;
        }
        cands.clear();
        self.cands = cands;
    }
}

impl Ctx {
    fn age_type(&self, f: &str, div: f64) -> u8 {
        if self.teage {
            return T_UNDEF;
        }
        let name = if self.o.dir { filename(&self.input) } else { filename(f) };
        let (min, max) = match self.age.get(name.as_bytes()) {
            Some((a, b)) => (num(a), num(b)),
            None => (0.0, 0.0),
        };
        let mut ty = T_NA;
        if div <= min {
            ty = T_LS;
        }
        if div >= max {
            ty = T_A;
        }
        if div < max && div > min {
            ty = T_ND;
        }
        ty
    }

    fn aggregate(&self, d: &FileData) -> Agg {
        let nn = d.names.len();
        let mut g = Agg {
            median_runs: (0..nn).map(|_| Vec::new()).collect(),
            name_nr: vec![0; nn],
            name_dbl: vec![0; nn],
            name_cnt_nr: vec![0; nn],
            name_cnt_tot: vec![0; nn],
            land_name: vec![None; nn],
            ..Default::default()
        };
        for s in &d.names.strs {
            let (c, f) = self.te.get(&lc(s)).cloned().unwrap_or_default();
            g.class_of.push(c);
            g.fam_of.push(f);
        }
        let p = self.p;
        let a = self.a;
        let rage = a && self.teage;
        let ncols = self.land.as_ref().map_or(0, |l| l.header.len());

        // Per-hit derived values.
        let types: Vec<u8> = if a {
            d.hits.iter().map(|h| self.age_type(&d.path, h.div)).collect()
        } else {
            Vec::new()
        };
        let cols: Vec<Option<Option<usize>>> = match &self.land {
            Some(land) => {
                let mut cache: FxHashMap<u64, Option<Option<usize>>> = FxHashMap::default();
                d.hits
                    .iter()
                    .map(|h| *cache.entry(h.div.to_bits()).or_insert_with(|| land.locate(h.div)))
                    .collect()
            }
            None => Vec::new(),
        };

        let mut seen_full_p = vec![false; if p { d.rfull.len() } else { 0 }];
        let mut seen_rid_p = vec![false; if p { d.rid.len() } else { 0 }];
        let mut fs_full = FirstSeen::new(if a { d.rfull.len() } else { 0 });
        let mut fs_rid = FirstSeen::new(if a { d.rid.len() } else { 0 });

        let mut sh: Vec<SweepHit> = Vec::new();
        let hits = &d.hits;
        for (group_index, group) in hits.chunk_by(|x, y| x.gname == y.gname).enumerate() {
            let base = (group.as_ptr() as usize - hits.as_ptr() as usize) / std::mem::size_of::<Hit>();
            sh.clear();
            sh.extend(group.iter().map(|h| SweepHit {
                start: h.start,
                end: h.end,
                div: h.div,
                truthy: h.div_truthy,
            }));
            let totals = sweep(&sh, |r: Run| {
                let wi = base + r.winner;
                let h = &hits[wi];
                let len = (r.b - r.a + 1) as u64;
                let ni = h.name as usize;
                if p {
                    g.name_nr[ni] += len;
                    if r.double {
                        g.name_dbl[ni] += len;
                    }
                    if !seen_full_p[h.rfull as usize] {
                        seen_full_p[h.rfull as usize] = true;
                        g.name_cnt_nr[ni] += 1;
                    }
                    if !seen_rid_p[h.rid as usize] {
                        seen_rid_p[h.rid as usize] = true;
                        g.name_cnt_tot[ni] += 1;
                    }
                    if !d.is_align {
                        g.median_runs[ni].push(MedianRun {
                            group: group_index,
                            a: r.a,
                            b: r.b,
                            values: [h.div_s, h.del_s, h.ins_s],
                        });
                    }
                    *g.hist[0].entry((h.name, h.div_s)).or_insert(0) += len;
                    *g.hist[1].entry((h.name, h.del_s)).or_insert(0) += len;
                    *g.hist[2].entry((h.name, h.ins_s)).or_insert(0) += len;
                }
                if a {
                    let ty = types[wi];
                    g.age_nr[ty as usize] += len;
                    if r.double {
                        g.age_dbl[ty as usize] += len;
                    }
                    fs_full.push(h.rfull, ty, r.a, r.b);
                    fs_rid.push(h.rid, ty, r.a, r.b);
                }
                if let Some(Some(col)) = cols.get(wi).copied() {
                    let row = g.land_name[ni].get_or_insert_with(|| vec![0; ncols]);
                    let rrow = if rage {
                        Some(g.land_rage.get_or_insert_with(|| vec![0; ncols]))
                    } else {
                        None
                    };
                    if let Some(c) = col {
                        row[c] += len;
                        if let Some(rr) = rrow {
                            rr[c] += len;
                        }
                    }
                }
            });
            g.tot += totals.tot;
            g.double += totals.double;
            if a {
                fs_full.flush(&mut g.age_cnt_nr);
                fs_rid.flush(&mut g.age_cnt_tot);
            }
        }
        g
    }
}

// ---------------------------------------------------------------------------
// Output
// ---------------------------------------------------------------------------

fn push_str(out: &mut Vec<u8>, s: &str) {
    out.extend_from_slice(s.as_bytes());
}

/// get_avg_med_from_list for one value list given as a weighted histogram.
fn avg_med(entries: &mut [(f64, u32, u64)], vals: &Interner, with_median: bool) -> (String, String) {
    // Histogram entries follow hash iteration order, not position order.
    // Distinct strings tied numerically need the tied_median interval override.
    entries.sort_by(|x, y| x.0.partial_cmp(&y.0).unwrap_or(std::cmp::Ordering::Equal));
    let count: u64 = entries.iter().map(|e| e.2).sum();
    let total: f64 = entries.iter().map(|e| e.0 * e.2 as f64).sum();
    let avg = fmt(total / count as f64);
    if !with_median {
        return (avg, "nd".into());
    }
    let at = |k: u64| -> usize {
        let mut acc = 0u64;
        for (i, e) in entries.iter().enumerate() {
            acc += e.2;
            if k < acc {
                return i;
            }
        }
        entries.len() - 1
    };
    let med = if count % 2 == 1 {
        String::from_utf8_lossy(&vals.strs[entries[at(count / 2)].1 as usize]).into_owned()
    } else {
        fmt((entries[at(count / 2)].0 + entries[at(count / 2 - 1)].0) / 2.0)
    };
    (avg, med)
}

/// Resolve the odd median's stable numeric tie using group order and lexical
/// position order. Counting decimal-prefix intersections avoids expanding bases.
fn tied_median(
    entries: &[(f64, u32, u64)],
    runs: &[MedianRun],
    vals: &Interner,
    column: usize,
) -> Option<String> {
    let count: u64 = entries.iter().map(|e| e.2).sum();
    if count.is_multiple_of(2) {
        return None;
    }
    let mut rank = count / 2;
    let value = entries.iter().find_map(|e| {
        if rank < e.2 {
            Some(e.0)
        } else {
            rank -= e.2;
            None
        }
    })?;
    let ties: Vec<_> = entries.iter().filter(|e| e.0 == value).collect();
    if ties.len() < 2 {
        return None;
    }
    rank = count / 2
        - entries
            .iter()
            .filter(|e| e.0 < value)
            .map(|e| e.2)
            .sum::<u64>();
    let mut ranges = Vec::new();
    for group in runs.chunk_by(|x, y| x.group == y.group) {
        ranges.clear();
        for r in group {
            let si = r.values[column];
            if ties.iter().any(|e| e.1 == si) {
                ranges.push((r.a as i128, r.b as i128, si));
            }
        }
        let n: u64 = ranges.iter().map(|r| (r.1 - r.0 + 1) as u64).sum();
        if rank >= n {
            rank -= n;
            continue;
        }
        let si = lexical_value_at(&ranges, rank);
        return Some(String::from_utf8_lossy(&vals.strs[si as usize]).into_owned());
    }
    unreachable!("median rank is covered by winning intervals")
}

fn lexical_value_at(ranges: &[(i128, i128, u32)], mut rank: u64) -> u32 {
    // Negative keys sort before zero and positives; after '-' magnitudes use
    // the same decimal trie as positive keys.
    let negative_count: u64 = ranges
        .iter()
        .map(|r| (r.1.min(-1) - r.0 + 1).max(0) as u64)
        .sum();
    let negative = rank < negative_count;
    let ranges: Vec<_> = if negative {
        ranges
            .iter()
            .filter(|r| r.0 < 0)
            .map(|r| (-r.1.min(-1), -r.0, r.2))
            .collect()
    } else {
        rank -= negative_count;
        ranges
            .iter()
            .filter(|r| r.1 >= 0)
            .map(|r| (r.0.max(0), r.1, r.2))
            .collect()
    };
    if let Some(r) = ranges.iter().find(|r| r.0 == 0) {
        if rank == 0 {
            return r.2;
        }
        rank -= 1;
    }
    let max = ranges.iter().map(|r| r.1).max().unwrap();
    let mut prefix = 0;
    loop {
        for digit in if prefix == 0 { 1..=9 } else { 0..=9 } {
            let candidate = prefix * 10 + digit;
            let mut lo = candidate;
            let mut scale = 1;
            let mut n = 0u64;
            while lo <= max {
                n += ranges
                    .iter()
                    .map(|r| (r.1.min(lo + scale - 1) - r.0.max(lo) + 1).max(0) as u64)
                    .sum::<u64>();
                lo *= 10;
                scale *= 10;
            }
            if rank >= n {
                rank -= n;
                continue;
            }
            prefix = candidate;
            break;
        }
        if let Some(r) = ranges.iter().find(|r| r.0 <= prefix && prefix <= r.1) {
            if rank == 0 {
                return r.2;
            }
            rank -= 1;
        }
    }
}

impl Ctx {
    fn print_outputs(&mut self, f: &str, d: &FileData, g: &Agg) -> R<()> {
        let nr = g.tot as i64 - g.double as i64;
        if self.p {
            self.print_parsed(f, d, g, nr)?;
        }
        if self.a {
            self.print_age(f, g, nr)?;
        }
        if self.land.is_some() {
            self.print_landscape(f, d, g)?;
        }
        Ok(())
    }

    fn print_parsed(&mut self, f: &str, d: &FileData, g: &Agg, nr: i64) -> R<()> {
        let allrep = format!("{f}.parseRM.all-repeats.tab");
        // :761 print_parsed_summary
        let fname = filename(f);
        let mut faname = strip_rm_ext(fname).unwrap_or(fname).to_string();
        if !looks_fasta(&faname) {
            faname.push_str(".fa");
        }
        let fakey = faname.into_bytes();
        let glen = match self.genlen.get(&fakey) {
            Some(v) if truthy(v) => v.clone(),
            _ => {
                self.genlen.insert(fakey.clone(), b"nd".to_vec());
                b"nd".to_vec()
            }
        };
        let gl_nd = glen == b"nd";
        let gl = num(&glen);
        let per = |v: f64, mult: f64| -> R<String> {
            if gl_nd {
                return Ok("nd".into());
            }
            if gl == 0.0 {
                die!("Illegal division by zero at parseRM.pl line 777.\n");
            }
            Ok(fmt(v / gl * mult))
        };
        let mut out = Vec::with_capacity(1 << 16);
        push_str(&mut out, "#This file gives some summary info about the masking, total, by class and by family\n");
        push_str(&mut out, &format!("#(all repeats are detailed in the file {allrep})\n"));
        push_str(&mut out, "#Overlap or double corresponds to DNA fragments that are masked by several elements.\n");
        push_str(&mut out, "#These amounts need to be subtracted in order to get more accurate TE amount.\n");
        push_str(&mut out, "#If a .align file was parsed, these amounts will be much higher than for the associated .out\n");
        push_str(&mut out, "#Note that overlaps may not be parsed correctly if names are not formatted consistently (Rname#Rclass/Rfam)\n");
        push_str(&mut out, "\n#TOTAL:\n");
        push_str(&mut out, "#nt_total_in_genome\tnt_masked-minus-double\t%_masked\tFYI:nt_masked_double\n");
        out.extend_from_slice(&glen);
        push_str(&mut out, &format!("\t{}\t{}\t{}\n", nr, per(nr as f64, 100.0)?, g.double));

        let mut by_class: FxHashMap<&[u8], (u64, u64)> = FxHashMap::default();
        let mut by_fam: FxHashMap<&[u8], (u64, u64)> = FxHashMap::default();
        for ni in 0..d.names.len() {
            if g.name_nr[ni] == 0 {
                continue;
            }
            let e = by_class.entry(&g.class_of[ni]).or_default();
            e.0 += g.name_nr[ni];
            e.1 += g.name_dbl[ni];
            let e = by_fam.entry(&g.fam_of[ni]).or_default();
            e.0 += g.name_nr[ni];
            e.1 += g.name_dbl[ni];
        }
        for (title, head, map) in [
            ("\n#BY CLASS\n", "#class", &by_class),
            ("\n#BY FAMILY\n", "#family", &by_fam),
        ] {
            push_str(&mut out, title);
            push_str(&mut out, &format!("{head}\tnt_masked-minus-double\t%_masked\tFYI:nt_masked_double\n"));
            let mut keys: Vec<_> = map.iter().collect();
            keys.sort();
            for (k, &(nt, db)) in keys {
                out.extend_from_slice(k);
                push_str(&mut out, &format!("\t{}\t{}\t{}\n", nt, per(nt as f64, 100.0)?, db));
            }
        }
        write(&format!("{f}.parseRM.summary.tab"), &out, "sub print_parsed_summary")?;

        // :811 print_parsed_allrep
        let mut per_name: Vec<[Vec<(f64, u32, u64)>; 3]> = (0..d.names.len()).map(|_| Default::default()).collect();
        for (k, hist) in g.hist.iter().enumerate() {
            for (&(ni, si), &w) in hist {
                per_name[ni as usize][k].push((num(&d.vals.strs[si as usize]), si, w));
            }
        }
        let mut out = Vec::with_capacity(1 << 16);
        push_str(&mut out, "#Rname\tRclass\tRfam\tRlen\tFRG_NB_all\tFRG_NB_Reconstructed_repeats");
        push_str(&mut out, "\tLEN_MASKED_NR\tAVG_%DIV\tMED_%DIV\tAVG_%DEL\tMED_%DEL\tAVG_%INS\tMED_%INS\tAVG_LEN_MASKED\t%_GENOME\n\n");
        let mut order: Vec<usize> = (0..d.names.len()).filter(|&ni| g.name_nr[ni] > 0).collect();
        order.sort_by(|&x, &y| d.names.strs[x].cmp(&d.names.strs[y]));
        for ni in order {
            let name = &d.names.strs[ni];
            let rlen = match self.liblen.get(&lc(name)) {
                Some(Some(v)) => String::from_utf8_lossy(v).into_owned(),
                _ => "nd".into(),
            };
            let [div, del, ins] = &mut per_name[ni];
            let med = !d.is_align;
            let (diva, mut divm) = avg_med(div, &d.vals, med);
            let (dela, mut delm) = avg_med(del, &d.vals, med);
            let (insa, mut insm) = avg_med(ins, &d.vals, med);
            if med {
                for (column, entries, median) in [
                    (0, &*div, &mut divm),
                    (1, &*del, &mut delm),
                    (2, &*ins, &mut insm),
                ] {
                    if let Some(ordered) = tied_median(entries, &g.median_runs[ni], &d.vals, column)
                    {
                        *median = ordered;
                    }
                }
            }
            let len = g.name_nr[ni];
            let avglen = fmt(len as f64 / nr as f64);
            let perlen = per(len as f64, 1.0)?;
            out.extend_from_slice(name);
            out.push(b'\t');
            out.extend_from_slice(&g.class_of[ni]);
            out.push(b'\t');
            out.extend_from_slice(&g.fam_of[ni]);
            push_str(
                &mut out,
                &format!(
                    "\t{rlen}\t{}\t{}\t{len}\t{diva}\t{divm}\t{dela}\t{delm}\t{insa}\t{insm}\t{avglen}\t{perlen}\n",
                    g.name_cnt_tot[ni], g.name_cnt_nr[ni]
                ),
            );
        }
        write(&allrep, &out, "sub print_parsed_allrep")
    }

    /// :1017 print_age (+ `$IN.splitage_all.tab` with -dir)
    fn print_age(&self, f: &str, g: &Agg, nr: i64) -> R<()> {
        let mut out = Vec::new();
        push_str(&mut out, "#nr_masked = amount of masked nt to consider\n");
        push_str(&mut out, "#tot_masked = total amount of masked nt, including redundant maskins / overlaps\n");
        push_str(&mut out, "#double_masked = total amount of nt that were masked by at least 2 repeats = redundant masking / overlaps\n\n");
        push_str(&mut out, "#Input_file\tnt_masked-minus-double\ttot_masked\tnt_masked_double\tAgeCat");
        push_str(&mut out, "\tCounts_this_age\tnr_Counts_this_age\tnr_masked_this_age\t%nr_masked_this_age\n\n");
        for (ti, name) in AGE_TYPES.iter().enumerate() {
            let masked = g.age_nr[ti];
            if masked == 0 {
                continue;
            }
            let na = |v: u64| if v > 0 { v.to_string() } else { "na".into() };
            push_str(
                &mut out,
                &format!(
                    "{f}\t{nr}\t{}\t{}\t{name}\t{}\t{}\t{masked}\t{}\n",
                    g.tot,
                    g.double,
                    na(g.age_cnt_tot[ti]),
                    na(g.age_cnt_nr[ti]),
                    fmt(masked as f64 / nr as f64 * 100.0)
                ),
            );
        }
        write(&format!("{f}.splitage.tab"), &out, "sub print_age")?;
        if self.o.dir {
            write(&format!("{}.splitage_all.tab", self.input), &out, "main")?;
        }
        Ok(())
    }

    /// :895 print_landscape
    fn print_landscape(&self, f: &str, d: &FileData, g: &Agg) -> R<()> {
        let land = self.land.as_ref().expect("land");
        let (n, label) = if self.my { ("My", "Million years") } else { ("Div", "% of divergence") };
        let ncols = land.header.len();
        let mut cols_hdr = String::new();
        for (i, i2) in &land.header {
            cols_hdr.push_str(&format!("\t[{i};{i2}["));
        }
        cols_hdr.push_str("\n\n");

        let mut rname: Vec<(Vec<u8>, &Vec<u64>)> = Vec::new();
        let mut rclass: FxHashMap<Vec<u8>, Vec<u64>> = FxHashMap::default();
        let mut rfam: FxHashMap<Vec<u8>, Vec<u64>> = FxHashMap::default();
        for (ni, row) in g.land_name.iter().enumerate() {
            let Some(row) = row else { continue };
            let (c, fam) = (&g.class_of[ni], &g.fam_of[ni]);
            let mut k = d.names.strs[ni].clone();
            k.push(b'\t');
            k.extend_from_slice(c);
            k.push(b'\t');
            k.extend_from_slice(fam);
            rname.push((k, row));
            let add = |acc: &mut Vec<u64>| acc.iter_mut().zip(row).for_each(|(x, y)| *x += y);
            add(rclass.entry(c.clone()).or_insert_with(|| vec![0; ncols]));
            let mut fk = c.clone();
            fk.push(b'\t');
            fk.extend_from_slice(fam);
            add(rfam.entry(fk).or_insert_with(|| vec![0; ncols]));
        }
        let emit = |ty: &str, header: &str, mut rows: Vec<(Vec<u8>, &Vec<u64>)>| -> R<()> {
            if rows.is_empty() {
                return Ok(());
            }
            rows.sort();
            let mut out = Vec::new();
            push_str(&mut out, header);
            push_str(&mut out, &cols_hdr);
            for (k, row) in rows {
                out.extend_from_slice(&k);
                for v in row {
                    push_str(&mut out, &format!("\t{v}"));
                }
                out.push(b'\n');
            }
            write(&format!("{f}.landscape.{n}.{ty}.tab"), &out, "sub print_landscape")
        };
        emit("Rname", &format!("\t\t\t{label} bins:\nRname\tRclass\tRfam"), rname)?;
        emit("Rclass", &format!("\t{label} bins:\nRclass"), rclass.iter().map(|(k, v)| (k.clone(), v)).collect())?;
        emit("Rfam", &format!("\t\t{label} bins:\nRclass\tRfam"), rfam.iter().map(|(k, v)| (k.clone(), v)).collect())?;
        if let Some(rr) = &g.land_rage {
            emit("Rage", &format!("\t{label} bins:\nLineage"), vec![(REF_KEY.as_bytes().to_vec(), rr)])?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lexical_median_selection_matches_position_enumeration() {
        for offset in [-25, 0, 90, 990] {
            let ranges = [(offset, offset + 12, 0), (offset + 13, offset + 25, 1)];
            let mut positions: Vec<_> = ranges
                .iter()
                .flat_map(|&(a, b, si)| (a..=b).map(move |p| (p.to_string(), si)))
                .collect();
            positions.sort();
            for (rank, (_, want)) in positions.iter().enumerate() {
                assert_eq!(
                    lexical_value_at(&ranges, rank as u64),
                    *want,
                    "offset {offset}, rank {rank}"
                );
            }
        }
    }

    #[test]
    fn median_ties_preserve_scaffold_group_order() {
        let mut vals = Interner::default();
        let decimal = vals.id(b"5.0");
        let integer = vals.id(b"5");
        let entries = [(5.0, decimal, 2), (5.0, integer, 1)];
        let runs = [
            MedianRun {
                group: 0,
                a: 100,
                b: 100,
                values: [decimal; 3],
            },
            MedianRun {
                group: 1,
                a: 1,
                b: 1,
                values: [integer; 3],
            },
            MedianRun {
                group: 2,
                a: 10,
                b: 10,
                values: [decimal; 3],
            },
        ];
        assert_eq!(tied_median(&entries, &runs, &vals, 0).as_deref(), Some("5"));
    }

    #[test]
    fn parse_name_cuts_like_perl_split() {
        assert_eq!(parse_name(b"AB-#-C"), b"AB");
        assert_eq!(parse_name(b"X-#"), b"X");
        assert_eq!(parse_name(b"a-"), b"a-");
        assert_eq!(parse_name(b"L1PA2"), b"L1PA2");
        assert_eq!(parse_name(b""), b"");
    }

    #[test]
    fn fasta_and_extension_helpers() {
        assert!(looks_fasta("genome.fa"));
        assert!(looks_fasta("genome.fasta"));
        assert!(!looks_fasta("genome"));
        assert_eq!(strip_rm_ext("a.out.b.out"), Some("a.out.b"));
        assert_eq!(strip_rm_ext("x.align"), Some("x"));
    }
}
