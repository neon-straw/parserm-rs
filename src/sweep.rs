//! Sweep-line replacement for the per-nucleotide loop of `parseRM_table`.
//!
//! Perl, for every position `n` of every hit (in file order):
//! ```text
//! $DOUBLE{$n}++ if ($LOWDIV{$n}); $TOT{double}++ if ($LOWDIV{$n});
//! if (! $LOWDIV{$n} || $LOWDIV{$n} > $DIV) { $LOWDIV{$n} = $DIV; winner = hit }
//! ```
//! Between two consecutive hit boundaries the set of covering hits is fixed,
//! so the outcome is the same for every position of that elementary segment.
//! When every covering div is truthy the winner is the minimal div (ties: the
//! first hit in file order) and each position gets `depth - 1` double
//! arrivals. When a falsy div (numeric 0, string "0") is involved, the
//! arrival sequence is replayed in file order to reproduce the quirk.

use std::collections::BTreeSet;

#[derive(Clone, Copy, Debug)]
pub struct SweepHit {
    pub start: i64,
    pub end: i64,
    pub div: f64,
    pub truthy: bool,
}

/// Maximal run of positions `[a, b]` with the same winner and the same
/// "position was double masked" flag.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Run {
    pub a: i64,
    pub b: i64,
    pub winner: usize,
    pub double: bool,
}

#[derive(Default, Debug, Clone, Copy, PartialEq, Eq)]
pub struct Totals {
    pub tot: u64,
    pub double: u64,
}

/// Total order on f64 consistent with Perl's numeric `>` for non-NaN values
/// (-0.0 and 0.0 compare equal).
#[inline]
fn order_key(v: f64) -> u64 {
    let v = if v == 0.0 { 0.0 } else { v };
    let b = v.to_bits();
    if b >> 63 == 1 {
        !b
    } else {
        b | (1 << 63)
    }
}

/// Process one scaffold group. `emit` receives runs in increasing position order.
pub fn sweep(hits: &[SweepHit], mut emit: impl FnMut(Run)) -> Totals {
    let mut totals = Totals::default();
    let mut events: Vec<(i64, u32, bool)> = Vec::with_capacity(hits.len() * 2);
    for (i, h) in hits.iter().enumerate() {
        if h.start <= h.end {
            totals.tot += (h.end - h.start + 1) as u64;
            events.push((h.start, i as u32, true));
            events.push((h.end + 1, i as u32, false));
        }
    }
    events.sort_unstable_by_key(|e| e.0);

    let mut by_div: BTreeSet<(u64, u32)> = BTreeSet::new();
    let mut by_idx: BTreeSet<u32> = BTreeSet::new();
    let mut falsy = 0usize;
    let mut pending: Option<Run> = None;

    let mut k = 0;
    while k < events.len() {
        let p = events[k].0;
        while k < events.len() && events[k].0 == p {
            let (_, i, add) = events[k];
            let h = &hits[i as usize];
            if add {
                by_div.insert((order_key(h.div), i));
                by_idx.insert(i);
                falsy += usize::from(!h.truthy);
            } else {
                by_div.remove(&(order_key(h.div), i));
                by_idx.remove(&i);
                falsy -= usize::from(!h.truthy);
            }
            k += 1;
        }
        if k == events.len() || by_idx.is_empty() {
            continue;
        }
        let q = events[k].0 - 1;
        let (winner, arrivals) = if falsy == 0 {
            let &(_, w) = by_div.iter().next().expect("non-empty");
            (w as usize, (by_idx.len() - 1) as u64)
        } else {
            let mut cur: Option<usize> = None;
            let mut c = 0u64;
            for &i in &by_idx {
                let i = i as usize;
                match cur {
                    Some(w) if hits[w].truthy => {
                        c += 1;
                        if hits[w].div > hits[i].div {
                            cur = Some(i);
                        }
                    }
                    _ => cur = Some(i),
                }
            }
            (cur.expect("non-empty"), c)
        };
        totals.double += arrivals * (q - p + 1) as u64;
        let double = arrivals > 0;
        match &mut pending {
            Some(r) if r.winner == winner && r.double == double && r.b + 1 == p => r.b = q,
            _ => {
                if let Some(r) = pending.take() {
                    emit(r);
                }
                pending = Some(Run {
                    a: p,
                    b: q,
                    winner,
                    double,
                });
            }
        }
    }
    if let Some(r) = pending {
        emit(r);
    }
    totals
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeMap;

    /// Literal port of the Perl per-position loop.
    fn naive(hits: &[SweepHit]) -> (BTreeMap<i64, (usize, bool)>, Totals) {
        let mut low: BTreeMap<i64, usize> = BTreeMap::new();
        let mut dbl: BTreeMap<i64, u64> = BTreeMap::new();
        let mut t = Totals::default();
        for (i, h) in hits.iter().enumerate() {
            let mut n = h.start;
            while n <= h.end {
                t.tot += 1;
                let cur = low.get(&n).copied();
                let cur_true = cur.is_some_and(|w| hits[w].truthy);
                if cur_true {
                    *dbl.entry(n).or_default() += 1;
                    t.double += 1;
                }
                if !cur_true || hits[cur.unwrap()].div > h.div {
                    low.insert(n, i);
                }
                n += 1;
            }
        }
        let out = low
            .into_iter()
            .map(|(n, w)| (n, (w, dbl.get(&n).copied().unwrap_or(0) > 0)))
            .collect();
        (out, t)
    }

    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            self.0 >> 33
        }
    }

    #[test]
    fn sweep_equals_naive() {
        let mut rng = Lcg(42);
        for round in 0..3000 {
            let nh = 1 + (rng.next() % 12) as usize;
            let span = 5 + rng.next() % 60;
            let hits: Vec<SweepHit> = (0..nh)
                .map(|_| {
                    let s = (rng.next() % span) as i64 + 1;
                    let len = (rng.next() % 25) as i64 - 2;
                    let dv = (rng.next() % 6) as f64;
                    let falsy_mode = round % 3 == 0;
                    SweepHit {
                        start: s,
                        end: s + len,
                        div: dv,
                        truthy: if falsy_mode { dv != 0.0 } else { true },
                    }
                })
                .collect();
            let (want, want_t) = naive(&hits);
            let mut got = BTreeMap::new();
            let mut prev_b = i64::MIN;
            let t = sweep(&hits, |r| {
                assert!(r.a > prev_b);
                prev_b = r.b;
                for n in r.a..=r.b {
                    got.insert(n, (r.winner, r.double));
                }
            });
            assert_eq!(t, want_t, "totals, round {round}: {hits:?}");
            assert_eq!(got, want, "positions, round {round}: {hits:?}");
        }
    }
}
