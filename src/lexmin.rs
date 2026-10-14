//! Lexicographically smallest decimal representation in an integer range.
//!
//! parseRM.pl walks positions with `sort keys %LOWDIV` (string order); the
//! first position seen for a fragment ID decides which age category gets the
//! fragment count. With intervals instead of positions we need the
//! string-minimal position of each interval.

fn pow10(k: u32) -> i128 {
    10i128.pow(k)
}

fn digits(v: i64) -> u32 {
    let mut d = 1;
    let mut x = v / 10;
    while x > 0 {
        d += 1;
        x /= 10;
    }
    d
}

/// Does some integer in `[a, b]` have a decimal representation starting with `prefix`?
fn has_prefix_in(prefix: i128, plen: u32, a: i64, b: i64, maxlen: u32) -> bool {
    for len in plen..=maxlen {
        let scale = pow10(len - plen);
        let lo = prefix * scale;
        let hi = lo + scale - 1;
        if lo.max(a as i128) <= hi.min(b as i128) {
            return true;
        }
    }
    false
}

/// String-minimal integer in `[a, b]` (requires `a <= b`). Negative ranges
/// fall back to plain enumeration (never produced by RepeatMasker).
pub fn lexmin(a: i64, b: i64) -> String {
    debug_assert!(a <= b);
    if a < 0 {
        return (a..=b.min(a.saturating_add(100_000)))
            .map(|v| v.to_string())
            .min()
            .unwrap_or_default();
    }
    if a == 0 {
        return "0".to_string();
    }
    let maxlen = digits(b);
    let mut prefix: i128 = 0;
    let mut plen: u32 = 0;
    loop {
        if plen > 0 && prefix >= a as i128 && prefix <= b as i128 {
            return prefix.to_string();
        }
        let first = if plen == 0 { 1 } else { 0 };
        let mut found = false;
        for d in first..=9 {
            let q = prefix * 10 + d;
            if has_prefix_in(q, plen + 1, a, b, maxlen) {
                prefix = q;
                plen += 1;
                found = true;
                break;
            }
        }
        assert!(found, "lexmin: empty range");
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn brute_force() {
        for a in 0..400i64 {
            for b in a..(a + 260) {
                let want = (a..=b).map(|v| v.to_string()).min().unwrap();
                assert_eq!(lexmin(a, b), want, "range {a}..{b}");
            }
        }
        assert_eq!(lexmin(95, 1203), "100");
        assert_eq!(lexmin(2_000_001, 2_999_999), "2000001");
        assert_eq!(lexmin(999_999_998, 1_000_000_001), "1000000000");
    }
}
