// A byte per value turns most pairs away on a quarter of the reads. The integer dot is exact and
// every rounding is bounded, so the bound never sits under the dot the full rows give.
pub(super) struct Coded {
    width: usize,
    codes: Vec<i8>,
    rows: Vec<Code>,
    rounding: f64,
}

#[derive(Clone, Copy)]
struct Code {
    scale: f64,
    kept: f64,
    lost: f64,
    norm: f64,
}

const STEPS: f64 = 127.0;

impl Coded {
    pub(super) fn of(tnf: &[f32], width: usize) -> Self {
        // Higham's bound on a recursive f32 dot of `width` terms, with room for the f64 sums here.
        let unit = (width as f64) * f64::from(f32::EPSILON) / 2.0;
        let rows = tnf.len().checked_div(width).unwrap_or(0);
        let mut coded = Self {
            width,
            codes: Vec::with_capacity(tnf.len()),
            rows: Vec::with_capacity(rows),
            rounding: unit / (1.0 - unit) + 1e-9,
        };
        for row in tnf.chunks_exact(width.max(1)).take(rows) {
            coded.push(row);
        }
        coded
    }

    fn push(&mut self, row: &[f32]) {
        let peak = row
            .iter()
            .fold(0.0f64, |peak, value| peak.max(f64::from(*value).abs()));
        let scale = peak / STEPS;
        let start = self.codes.len();
        self.codes.extend(row.iter().map(|value| match scale > 0.0 {
            true => (f64::from(*value) / scale).round().clamp(-STEPS, STEPS) as i8,
            false => 0,
        }));
        let (mut kept, mut lost, mut norm) = (0.0, 0.0, 0.0);
        for (value, code) in row.iter().zip(&self.codes[start..]) {
            let (value, coded) = (f64::from(*value), scale * f64::from(*code));
            kept += coded * coded;
            lost += (value - coded) * (value - coded);
            norm += value * value;
        }
        self.rows.push(Code {
            scale,
            kept: kept.sqrt(),
            lost: lost.sqrt(),
            norm: norm.sqrt(),
        });
    }

    pub(super) fn highest_dot(&self, a: usize, b: usize) -> f64 {
        let (x, y) = (self.rows[a], self.rows[b]);
        let exact = self
            .codes_of(a)
            .iter()
            .zip(self.codes_of(b))
            .map(|(p, q)| i32::from(*p) * i32::from(*q))
            .sum::<i32>();
        x.scale * y.scale * f64::from(exact)
            + x.kept * y.lost
            + y.kept * x.lost
            + x.lost * y.lost
            + self.rounding * x.norm * y.norm
    }

    fn codes_of(&self, row: usize) -> &[i8] {
        &self.codes[row * self.width..(row + 1) * self.width]
    }
}
