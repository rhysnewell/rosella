use rand::{Rng, SeedableRng, rngs::StdRng};

/// Hartigan and Hartigan (1985) on a weighted sample: the minorant runs through the step
/// below each point and the majorant through the step above, which unit weights hide.
pub fn dip(values: &[f64], weights: &[f64]) -> f64 {
    let (x, w) = sorted_unique(values, weights);
    weighted_dip(&x, &w)
}

pub fn exceeds_null(values: &[f64], weights: &[f64], draws: usize, seed: u64) -> bool {
    let observed = dip(values, weights);
    if observed <= 0.0 {
        return false;
    }
    let mut rng = StdRng::seed_from_u64(seed);
    let mut positions = vec![0.0; values.len()];
    for _ in 0..draws {
        for position in positions.iter_mut() {
            *position = rng.random::<f64>();
        }
        if dip(&positions, weights) >= observed {
            return false;
        }
    }
    true
}

fn sorted_unique(values: &[f64], weights: &[f64]) -> (Vec<f64>, Vec<f64>) {
    let mut order = (0..values.len()).collect::<Vec<_>>();
    order.sort_by(|a, b| values[*a].total_cmp(&values[*b]));
    let mut x = Vec::with_capacity(values.len());
    let mut w = Vec::with_capacity(values.len());
    for index in order {
        if weights[index] <= 0.0 {
            continue;
        }
        if x.last().is_some_and(|last: &f64| *last == values[index]) {
            *w.last_mut().expect("x and w grow together") += weights[index];
        } else {
            x.push(values[index]);
            w.push(weights[index]);
        }
    }
    (x, w)
}

fn weighted_dip(x: &[f64], w: &[f64]) -> f64 {
    let n = x.len();
    if n < 2 {
        return 0.0;
    }
    let mut xs = vec![0.0; n + 1];
    xs[1..].copy_from_slice(x);
    let mut lower = vec![0.0; n + 1];
    let mut upper = vec![0.0; n + 1];
    let mut total = 0.0;
    for j in 1..=n {
        lower[j] = total;
        total += w[j - 1];
        upper[j] = total;
    }
    if total <= 0.0 {
        return 0.0;
    }

    let mut mn = vec![0usize; n + 1];
    mn[1] = 1;
    for j in 2..=n {
        mn[j] = j - 1;
        loop {
            let mnj = mn[j];
            let mnmnj = mn[mnj];
            if mnj == 1
                || (xs[j] - xs[mnj]) * (lower[mnj] - lower[mnmnj])
                    < (xs[mnj] - xs[mnmnj]) * (lower[j] - lower[mnj])
            {
                break;
            }
            mn[j] = mnmnj;
        }
    }
    let mut mj = vec![0usize; n + 1];
    mj[n] = n;
    for k in (1..n).rev() {
        mj[k] = k + 1;
        loop {
            let mjk = mj[k];
            let mjmjk = mj[mjk];
            if mjk == n
                || (xs[k] - xs[mjk]) * (upper[mjk] - upper[mjmjk])
                    < (xs[mjk] - xs[mjmjk]) * (upper[k] - upper[mjk])
            {
                break;
            }
            mj[k] = mjmjk;
        }
    }

    let mut low = 1;
    let mut high = n;
    let mut dip = 0.0f64;
    let mut gcm = vec![0usize; n + 2];
    let mut lcm = vec![0usize; n + 2];
    while low < high {
        gcm[1] = high;
        let mut i = 1;
        while gcm[i] > low {
            gcm[i + 1] = mn[gcm[i]];
            i += 1;
        }
        let l_gcm = i;
        let mut ig = l_gcm;
        let mut ix = l_gcm - 1;

        lcm[1] = low;
        let mut i = 1;
        while lcm[i] < high {
            lcm[i + 1] = mj[lcm[i]];
            i += 1;
        }
        let l_lcm = i;
        let mut ih = l_lcm;
        let mut iv = 2;

        let mut d = 0.0f64;
        if l_gcm != 2 || l_lcm != 2 {
            loop {
                let gcmix = gcm[ix];
                let lcmiv = lcm[iv];
                if gcmix > lcmiv {
                    let gcmi1 = gcm[ix + 1];
                    let line = lower[gcmi1]
                        + (xs[lcmiv] - xs[gcmi1]) * (lower[gcmix] - lower[gcmi1])
                            / (xs[gcmix] - xs[gcmi1]);
                    let dx = upper[lcmiv] - line;
                    iv += 1;
                    if dx >= d {
                        d = dx;
                        ig = ix + 1;
                        ih = iv - 1;
                    }
                } else {
                    let lcmiv1 = lcm[iv - 1];
                    let line = upper[lcmiv1]
                        + (xs[gcmix] - xs[lcmiv1]) * (upper[lcmiv] - upper[lcmiv1])
                            / (xs[lcmiv] - xs[lcmiv1]);
                    let dx = line - lower[gcmix];
                    ix = ix.saturating_sub(1);
                    if dx >= d {
                        d = dx;
                        ig = ix + 1;
                        ih = iv;
                    }
                }
                if ix < 1 {
                    ix = 1;
                }
                if iv > l_lcm {
                    iv = l_lcm;
                }
                if gcm[ix] == lcm[iv] {
                    break;
                }
            }
        }
        if d < dip {
            break;
        }

        let mut dip_l = 0.0f64;
        for j in ig..l_gcm {
            let jb = gcm[j + 1];
            let je = gcm[j];
            let slope = (lower[je] - lower[jb]) / (xs[je] - xs[jb]);
            for jj in jb..=je {
                let t = upper[jj] - (lower[jb] + (xs[jj] - xs[jb]) * slope);
                dip_l = dip_l.max(t);
            }
        }
        let mut dip_u = 0.0f64;
        for j in ih..l_lcm {
            let jb = lcm[j];
            let je = lcm[j + 1];
            let slope = (upper[je] - upper[jb]) / (xs[je] - xs[jb]);
            for jj in jb..=je {
                let t = (upper[jb] + (xs[jj] - xs[jb]) * slope) - lower[jj];
                dip_u = dip_u.max(t);
            }
        }
        dip = dip.max(dip_l.max(dip_u));
        if low == gcm[ig] && high == lcm[ih] {
            break;
        }
        low = gcm[ig];
        high = lcm[ih];
    }
    dip / (2.0 * total)
}
