//! Helpers that reproduce Perl scalar semantics used by parseRM.pl:
//! string truthiness, numification, number stringification and `split`.

/// Perl `\s` (ASCII): space, \t, \n, \v, \f, \r.
#[inline]
pub fn is_space(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

/// Truthiness of a Perl string: false iff "" or "0".
#[inline]
pub fn truthy(s: &[u8]) -> bool {
    !(s.is_empty() || s == b"0")
}

/// `$s =~ /\w/` on a byte string (no locale, no utf8).
#[inline]
pub fn has_word(s: &[u8]) -> bool {
    s.iter().any(|&b| b.is_ascii_alphanumeric() || b == b'_')
}

/// Perl numification of a string: the longest leading decimal number, 0 otherwise.
pub fn num(s: &[u8]) -> f64 {
    let n = s.len();
    let mut i = 0;
    while i < n && is_space(s[i]) {
        i += 1;
    }
    let start = i;
    if i < n && (s[i] == b'+' || s[i] == b'-') {
        i += 1;
    }
    let digits_start = i;
    while i < n && s[i].is_ascii_digit() {
        i += 1;
    }
    let mut ndig = i - digits_start;
    if i < n && s[i] == b'.' {
        i += 1;
        let f = i;
        while i < n && s[i].is_ascii_digit() {
            i += 1;
        }
        ndig += i - f;
    }
    if ndig == 0 {
        let rest = s[digits_start..].to_ascii_lowercase();
        let neg = digits_start > start && s[start] == b'-';
        let v = if rest.starts_with(b"inf") {
            f64::INFINITY
        } else if rest.starts_with(b"nan") {
            f64::NAN
        } else {
            return 0.0;
        };
        return if neg { -v } else { v };
    }
    let mut end = i;
    if i < n && (s[i] == b'e' || s[i] == b'E') {
        let mut j = i + 1;
        if j < n && (s[j] == b'+' || s[j] == b'-') {
            j += 1;
        }
        let e = j;
        while j < n && s[j].is_ascii_digit() {
            j += 1;
        }
        if j > e {
            end = j;
        }
    }
    // The slice is ASCII by construction.
    std::str::from_utf8(&s[start..end])
        .ok()
        .and_then(|t| t.parse::<f64>().ok())
        .unwrap_or(0.0)
}

/// Perl stringification of an NV (`%.15g`, zero printed as "0").
pub fn fmt(v: f64) -> String {
    if v == 0.0 {
        return "0".to_string();
    }
    if v.is_nan() {
        return "NaN".to_string();
    }
    if v.is_infinite() {
        return if v > 0.0 { "Inf" } else { "-Inf" }.to_string();
    }
    const P: i32 = 15;
    let e_repr = format!("{:.*e}", (P - 1) as usize, v);
    let (mant, exp) = e_repr.split_once('e').expect("exponent");
    let x: i32 = exp.parse().expect("exponent value");
    if !(-4..P).contains(&x) {
        let mant = strip_zeros(mant);
        let sign = if x < 0 { '-' } else { '+' };
        format!("{}e{}{:02}", mant, sign, x.abs())
    } else {
        let prec = (P - 1 - x) as usize;
        strip_zeros(&format!("{:.*}", prec, v)).to_string()
    }
}

fn strip_zeros(s: &str) -> &str {
    if s.contains('.') {
        s.trim_end_matches('0').trim_end_matches('.')
    } else {
        s
    }
}

/// Perl `split(/\s+/, $s)`: a leading empty field when `$s` starts with
/// whitespace, trailing empty fields removed.
pub fn split_ws(s: &[u8]) -> Vec<&[u8]> {
    let n = s.len();
    let mut out = Vec::with_capacity(16);
    let mut i = 0;
    if n > 0 && is_space(s[0]) {
        out.push(&s[0..0]);
        while i < n && is_space(s[i]) {
            i += 1;
        }
    }
    while i < n {
        let st = i;
        while i < n && !is_space(s[i]) {
            i += 1;
        }
        out.push(&s[st..i]);
        while i < n && is_space(s[i]) {
            i += 1;
        }
    }
    while out.last().is_some_and(|f| f.is_empty()) {
        out.pop();
    }
    out
}

/// Perl `split(/c/, $s)` for a single literal byte: trailing empty fields removed.
pub fn split_byte(s: &[u8], c: u8) -> Vec<&[u8]> {
    let mut out: Vec<&[u8]> = s.split(|&b| b == c).collect();
    while out.last().is_some_and(|f| f.is_empty()) {
        out.pop();
    }
    out
}

/// Iterate over lines like `while (<$fh>) { chomp; ... }`.
pub fn lines(buf: &[u8]) -> impl Iterator<Item = &[u8]> {
    let mut pos = 0;
    std::iter::from_fn(move || {
        if pos >= buf.len() {
            return None;
        }
        let rest = &buf[pos..];
        match memchr::memchr(b'\n', rest) {
            Some(i) => {
                pos += i + 1;
                Some(&rest[..i])
            }
            None => {
                pos = buf.len();
                Some(rest)
            }
        }
    })
}

/// `s/(.*)\/.*$/$1/` or "." (parseRM.pl `get_path`).
pub fn get_path(file: &str) -> String {
    match file.rfind('/') {
        Some(i) => file[..i].to_string(),
        None => ".".to_string(),
    }
}

/// `s/.*\/(.*)$/$1/` (parseRM.pl `filename`).
pub fn filename(file: &str) -> &str {
    match file.rfind('/') {
        Some(i) => &file[i + 1..],
        None => file,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fmt_matches_perl() {
        assert_eq!(fmt(0.1 + 0.2), "0.3");
        assert_eq!(fmt(-0.0), "0");
        assert_eq!(fmt(1e21), "1e+21");
        assert_eq!(fmt(1e15), "1e+15");
        assert_eq!(fmt(1.0 / 3e7), "3.33333333333333e-08");
        assert_eq!(fmt(10.0), "10");
        assert_eq!(fmt(0.9375), "0.9375");
        assert_eq!(fmt(1.0 / 7.0), "0.142857142857143");
        assert_eq!(fmt(123456789012345.0), "123456789012345");
        assert_eq!(fmt(0.0001), "0.0001");
        assert_eq!(fmt(0.00001), "1e-05");
        assert_eq!(fmt(-2.5), "-2.5");
        assert_eq!(fmt(10.7325632730505), "10.7325632730505");
    }

    #[test]
    fn num_prefix() {
        assert_eq!(num(b"12.3"), 12.3);
        assert_eq!(num(b"  7abc"), 7.0);
        assert_eq!(num(b"(100)"), 0.0);
        assert_eq!(num(b"-1.5e3x"), -1500.0);
        assert_eq!(num(b"1e"), 1.0);
        assert_eq!(num(b".5"), 0.5);
        assert_eq!(num(b""), 0.0);
    }

    #[test]
    fn splits() {
        assert_eq!(split_ws(b" a  b "), vec![&b""[..], b"a", b"b"]);
        assert_eq!(split_ws(b"   "), Vec::<&[u8]>::new());
        assert_eq!(split_byte(b"DNA/", b'/'), vec![&b"DNA"[..]]);
        assert_eq!(split_byte(b"/x", b'/'), vec![&b""[..], b"x"]);
        assert_eq!(split_byte(b"", b'#'), Vec::<&[u8]>::new());
    }

    #[test]
    fn truthiness() {
        assert!(!truthy(b"0"));
        assert!(!truthy(b""));
        assert!(truthy(b"0.0"));
        assert!(truthy(b"00"));
    }
}
