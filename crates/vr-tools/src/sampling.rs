//! Walker/Vose alias-table source sampling.
//!
//! [`AliasTable`] implements Walker's method with Vose's construction
//! (Wuttke 2013), including the reverted small/large index ordering and the
//! `prob == 1` drain for numerically degenerate leftovers.
//! [`MeshSourceSampler`] applies it to [`MeshTallyData`]: ANALOG samples the
//! (volume-weighted) tally totals, UNIFORM samples space uniformly, and USER
//! consumes an external density array. Only voxel-level sampling is provided;
//! subvoxel/cell-fraction modes are out of scope.

use mcnp_io::meshtal::MeshTallyData;

use crate::Error;

/// Walker/Vose alias table over a discrete PDF.
///
/// Build once with [`AliasTable::new`], then draw with
/// [`AliasTable::sample`] using two uniform random numbers in `[0, 1)`.
#[derive(Debug, Clone, PartialEq)]
pub struct AliasTable {
    n: usize,
    /// Normalized input PDF, stored for reference/inspection.
    pdf: Vec<f64>,
    /// Per-bin cut probability (`<= 1`) as built by Vose's algorithm.
    prob: Vec<f64>,
    /// Per-bin alias index.
    alias: Vec<usize>,
}

impl AliasTable {
    /// Build an alias table from a non-negative PDF.
    ///
    /// The PDF is normalized internally so the standalone API is safe with
    /// unnormalized input (the underlying construction assumes a unit sum).
    pub fn new(pdf: &[f64]) -> Result<Self, Error> {
        if pdf.is_empty() {
            return Err(Error::EmptyPdf);
        }
        let mut sum = 0.0;
        for (i, &p) in pdf.iter().enumerate() {
            if !p.is_finite() {
                return Err(Error::NonFinitePdf { index: i });
            }
            if p < 0.0 {
                return Err(Error::NegativePdf { index: i, value: p });
            }
            sum += p;
        }
        if sum <= 0.0 {
            return Err(Error::ZeroSumPdf);
        }

        let n = pdf.len();
        // Normalized PDF: also what `pdf()` reports, since this is what the
        // table actually samples.
        let normalized_pdf: Vec<f64> = pdf.iter().map(|&x| x / sum).collect();
        let mut p = normalized_pdf.clone();

        // Scale so the mean bin probability is exactly 1.
        for x in p.iter_mut() {
            *x *= n as f64;
        }

        // Separate index lists for small and large probabilities. As in the
        // Wuttke implementation, indices are visited in reverted order.
        let mut small = Vec::with_capacity(n);
        let mut large = Vec::with_capacity(n);
        let mut i = n;
        while i > 0 {
            i -= 1;
            if p[i] < 1.0 {
                small.push(i);
            } else {
                large.push(i);
            }
        }

        let mut prob = vec![0.0; n];
        let mut alias = vec![0usize; n];
        while !small.is_empty() && !large.is_empty() {
            let a = small.pop().expect("small non-empty");
            let g = large.pop().expect("large non-empty");
            prob[a] = p[a];
            alias[a] = g;
            p[g] += p[a] - 1.0;
            if p[g] < 1.0 {
                small.push(g);
            } else {
                large.push(g);
            }
        }
        while let Some(g) = large.pop() {
            prob[g] = 1.0;
        }
        // Can only happen through numeric instability.
        while let Some(a) = small.pop() {
            prob[a] = 1.0;
        }

        Ok(Self {
            n,
            pdf: normalized_pdf,
            prob,
            alias,
        })
    }

    /// Draw one index using two uniforms in `[0, 1)`.
    ///
    /// Mirrors `sample_pdf(rand1, rand2)`: pick column `n * rand1`, return it
    /// when `rand2 < prob[column]`, else the column's alias. Values outside
    /// `[0, 1)` saturate instead of reading out of bounds (the C++ relies on
    /// caller discipline here).
    pub fn sample(&self, r1: f64, r2: f64) -> usize {
        let mut i = (self.n as f64 * r1) as usize;
        if i >= self.n {
            i = self.n - 1;
        }
        if r2 < self.prob[i] {
            i
        } else {
            self.alias[i]
        }
    }

    /// The normalized PDF the table samples from.
    pub fn pdf(&self) -> &[f64] {
        &self.pdf
    }

    /// Number of bins.
    pub fn len(&self) -> usize {
        self.n
    }

    /// Whether the table has no bins.
    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// Exact outcome probabilities implied by the table:
    /// `P(i) = (prob[i] + sum over j with alias[j] == i of (1 - prob[j])) / n`.
    #[cfg(test)]
    fn exact_probabilities(&self) -> Vec<f64> {
        let mut out = vec![0.0; self.n];
        for c in 0..self.n {
            let base = self.prob[c].min(1.0);
            out[c] += base / self.n as f64;
            if base < 1.0 {
                out[self.alias[c]] += (1.0 - base) / self.n as f64;
            }
        }
        out
    }
}

/// Source-sampling bias mode (values 0–2 at voxel level).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Mode {
    /// Sample where particles are born: PDF ∝ total source strength per voxel.
    #[default]
    Analog,
    /// Sample uniformly in space: PDF ∝ voxel volume.
    Uniform,
    /// Sample from an external (unnormalized) density array.
    User,
}

/// One sampled particle birth site (voxel resolution).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SampledVoxel {
    /// Flat volume-element index, x slowest → z fastest.
    pub index: usize,
    /// X cell index.
    pub i: usize,
    /// Y cell index.
    pub j: usize,
    /// Z cell index.
    pub k: usize,
    /// Birth weight: 1.0 for analog sampling, otherwise
    /// `analog_pdf[voxel] / biased_pdf[voxel]`.
    pub weight: f64,
}

/// Mesh source sampler over [`MeshTallyData`] energy-integrated totals.
///
/// The PDF lives over volume elements only (no energy dimension):
/// densities are multiplied by cell volumes before
/// normalization, since structured meshes may have unequal cells.
#[derive(Debug, Clone, PartialEq)]
pub struct MeshSourceSampler {
    dims: [usize; 3],
    num_ves: usize,
    mode: Mode,
    table: AliasTable,
    /// Birth weight per bin: `pdf / bias_pdf` elementwise (all ones in analog).
    biased_weights: Vec<f64>,
}

impl MeshSourceSampler {
    /// Build a sampler over `tally`'s totals in `mode`.
    ///
    /// For [`Mode::User`], `user_pdf` supplies one unnormalized density value
    /// per volume element (length must equal `tally.num_ves()`); it is
    /// scaled by cell volumes and normalized, like a bias tag.
    pub fn new(tally: &MeshTallyData, mode: Mode, user_pdf: Option<&[f64]>) -> Result<Self, Error> {
        let num_ves = tally.num_ves();
        if num_ves == 0 {
            return Err(Error::EmptyTally);
        }

        let volumes = cell_volumes(tally);
        let analog_pdf: Vec<f64> = tally
            .total_result
            .iter()
            .zip(volumes.iter())
            .map(|(q, v)| q.abs() * v)
            .collect();
        if analog_pdf.iter().all(|&x| x == 0.0) {
            return Err(Error::ZeroSumPdf);
        }

        let bias_pdf: Vec<f64> = match mode {
            Mode::Analog => analog_pdf.clone(),
            Mode::Uniform => volumes.clone(),
            Mode::User => {
                let user = user_pdf.ok_or(Error::EmptyPdf)?;
                if user.len() != num_ves {
                    return Err(Error::LengthMismatch {
                        expected: num_ves,
                        got: user.len(),
                    });
                }
                user.iter()
                    .zip(volumes.iter())
                    .map(|(q, v)| q.abs() * v)
                    .collect()
            }
        };

        let normalized_analog = normalized(&analog_pdf).ok_or(Error::ZeroSumPdf)?;
        let normalized_bias = normalized(&bias_pdf).ok_or(Error::ZeroSumPdf)?;
        let biased_weights: Vec<f64> = match mode {
            Mode::Analog => vec![1.0; num_ves],
            _ => normalized_analog
                .iter()
                .zip(normalized_bias.iter())
                .map(|(&a, &b)| if b > 0.0 { a / b } else { 1.0 })
                .collect(),
        };

        Ok(Self {
            dims: tally.dims(),
            num_ves,
            mode,
            table: AliasTable::new(&normalized_bias)?,
            biased_weights,
        })
    }

    /// Sample a birth voxel with two uniforms in `[0, 1)`.
    pub fn sample(&self, r1: f64, r2: f64) -> SampledVoxel {
        let index = self.table.sample(r1, r2);
        let nz = self.dims[2];
        let ny = self.dims[1];
        let k = index % nz;
        let j = (index / nz) % ny;
        let i = index / (nz * ny);
        SampledVoxel {
            index,
            i,
            j,
            k,
            weight: self.biased_weights[index],
        }
    }

    /// The bias mode this sampler was constructed with.
    pub fn mode(&self) -> Mode {
        self.mode
    }

    /// Number of voxels in the sampling domain.
    pub fn num_voxels(&self) -> usize {
        self.num_ves
    }

    /// The underlying alias table (biased/analog PDF included).
    pub fn table(&self) -> &AliasTable {
        &self.table
    }
}

fn cell_volumes(tally: &MeshTallyData) -> Vec<f64> {
    let d = tally.dims();
    let mut vols = Vec::with_capacity(tally.num_ves());
    for i in 0..d[0] {
        let dx = tally.x_bounds[i + 1] - tally.x_bounds[i];
        for j in 0..d[1] {
            let dy = tally.y_bounds[j + 1] - tally.y_bounds[j];
            for k in 0..d[2] {
                let dz = tally.z_bounds[k + 1] - tally.z_bounds[k];
                vols.push(dx * dy * dz);
            }
        }
    }
    vols
}

fn normalized(pdf: &[f64]) -> Option<Vec<f64>> {
    let sum: f64 = pdf.iter().sum();
    if sum.is_nan() || sum <= 0.0 {
        return None;
    }
    Some(pdf.iter().map(|&x| x / sum).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::magic::magic;

    /// Park–Miller minimal standard LCG; deterministic, returns values in
    /// (0, 1) exclusive.
    struct Lcg(u64);

    impl Lcg {
        fn next_f64(&mut self) -> f64 {
            self.0 = (16807 * self.0) % 2147483647;
            self.0 as f64 / 2147483647.0
        }
    }

    fn approx(actual: f64, want: f64, tol: f64) {
        assert!(
            (actual - want).abs() <= tol * want.abs().max(1e-30),
            "expected {want}, got {actual}"
        );
    }

    #[test]
    fn uniform_two_bin_pdf_builds_degenerate_table() {
        // [0.5, 0.5]: after scaling both bins equal 1, no aliasing needed.
        let t = AliasTable::new(&[0.5, 0.5]).unwrap();
        assert_eq!(t.pdf(), &[0.5, 0.5]);
        assert_eq!(t.prob, vec![1.0, 1.0]);
        assert_eq!(t.len(), 2);
        assert!(!t.is_empty());
    }

    #[test]
    fn skewed_pdf_exact_outcome_probabilities_match_input() {
        // The alias method is exact by construction: reconstructing P(i) from
        // (prob, alias) must reproduce the normalized PDF to machine precision.
        for pdf in [
            vec![0.9, 0.06, 0.03, 0.01],
            vec![0.5, 0.5],
            vec![0.01, 0.01, 0.96],
            vec![1.0, 2.0, 3.0, 4.0, 5.0],
        ] {
            let t = AliasTable::new(&pdf).unwrap();
            let exact = t.exact_probabilities();
            for (e, w) in exact.iter().zip(t.pdf().iter()) {
                approx(*e, *w, 1e-12);
            }
        }
    }

    #[test]
    fn skewed_nine_bin_table_exact_probabilities() {
        let pdf = [0.9, 0.02, 0.02, 0.02, 0.01, 0.01, 0.008, 0.007, 0.005];
        let t = AliasTable::new(&pdf).unwrap();
        assert_eq!(t.len(), 9);
        for (c, (&p, &a)) in t.prob.iter().zip(t.alias.iter()).enumerate() {
            assert!((0.0..=1.0).contains(&p), "prob[{c}] = {p} outside [0, 1]");
            assert!(a < 9, "alias[{c}] = {a} out of range");
        }
        let exact = t.exact_probabilities();
        for (&e, &w) in exact.iter().zip(t.pdf().iter()) {
            approx(e, w, 1e-12);
        }
    }

    #[test]
    fn sample_saturates_out_of_range_uniforms() {
        let t = AliasTable::new(&[0.25, 0.75]).unwrap();
        // r1 = 1.0 would read index n in the C++; here it clamps to n - 1.
        // Column 1 has prob 1.0, so any finite r2 keeps index 1.
        assert_eq!(t.sample(1.0, 0.999), 1);
        // Negative r1 clamps to column 0 (prob 0.5); small r2 keeps it,
        // large r2 takes its alias.
        assert_eq!(t.sample(-3.0, 0.25), 0);
        assert_eq!(t.sample(-3.0, 0.75), t.alias[0]);
    }

    #[test]
    fn bad_pdfs_rejected() {
        assert_eq!(AliasTable::new(&[]), Err(Error::EmptyPdf));
        assert_eq!(
            AliasTable::new(&[-0.1, 0.5]),
            Err(Error::NegativePdf {
                index: 0,
                value: -0.1
            })
        );
        assert!(matches!(
            AliasTable::new(&[f64::NAN, 0.5]),
            Err(Error::NonFinitePdf { index: 0 })
        ));
        assert_eq!(AliasTable::new(&[0.0, 0.0]), Err(Error::ZeroSumPdf));
    }

    #[test]
    fn exhaustive_sampling_statistics_match_pdf() {
        // 100k deterministic draws against the task's skewed [0.9, ...] pdf;
        // every realized frequency within ~1% absolute of its probability.
        let pdf = [0.9, 0.06, 0.03, 0.01];
        let t = AliasTable::new(&pdf).unwrap();
        let mut rng = Lcg(123_456_789);
        let n_draws = 100_000;
        let mut counts = [0usize; 4];
        for _ in 0..n_draws {
            let idx = t.sample(rng.next_f64(), rng.next_f64());
            counts[idx] += 1;
        }
        for (i, &c) in counts.iter().enumerate() {
            let freq = c as f64 / n_draws as f64;
            assert!(
                (freq - pdf[i]).abs() < 0.01,
                "bin {i}: freq {freq} vs pdf {}",
                pdf[i]
            );
        }
    }

    #[test]
    fn lcg_coverage_hits_every_column_and_alias_branch() {
        // With 16 bins, both outcomes of every column get exercised across
        // 100k draws; all bins receive samples.
        let pdf = [1.0f64; 16];
        let t = AliasTable::new(&pdf).unwrap();
        let mut rng = Lcg(987_654_321);
        let mut seen = [false; 16];
        for _ in 0..100_000 {
            seen[t.sample(rng.next_f64(), rng.next_f64())] = true;
        }
        assert!(seen.iter().all(|&s| s));
    }

    fn fixture_tally(name: &str, num: u32) -> MeshTallyData {
        let path = format!(
            "{}/../../fixtures/mcnp/meshtal/{name}",
            env!("CARGO_MANIFEST_DIR")
        );
        let m = mcnp_io::meshtal::Meshtal::from_file(path).unwrap();
        m.tallies[&num].clone()
    }

    #[test]
    fn fixture_sampler_analog_pdf_proportional_to_totals_times_volume() {
        let t = fixture_tally("mcnp_meshtal_single_meshtal.txt", 4);
        let s = MeshSourceSampler::new(&t, Mode::Analog, None).unwrap();
        assert_eq!(s.mode(), Mode::Analog);
        assert_eq!(s.num_voxels(), t.num_ves());

        let d = t.dims();
        let mut expect = Vec::with_capacity(t.num_ves());
        for ve in 0..t.num_ves() {
            let k = ve % d[2];
            let j = (ve / d[2]) % d[1];
            let i = ve / (d[2] * d[1]);
            let vol = (t.x_bounds[i + 1] - t.x_bounds[i])
                * (t.y_bounds[j + 1] - t.y_bounds[j])
                * (t.z_bounds[k + 1] - t.z_bounds[k]);
            expect.push(t.total_result[ve] * vol);
        }
        let norm: f64 = expect.iter().sum();
        let got = s.table().pdf();
        for (&g, &w) in got.iter().zip(expect.iter()) {
            approx(g, w / norm, 1e-12);
        }
        // Analog births carry unit weight regardless of draw.
        let mut rng = Lcg(555);
        for _ in 0..100 {
            let sv = s.sample(rng.next_f64(), rng.next_f64());
            approx(sv.weight, 1.0, 1e-15);
        }
    }

    #[test]
    fn magic_bounds_proportional_to_analog_pdf_on_equal_cells() {
        // On an equal-volume mesh both quantities are proportional to the
        // cell flux: MAGIC lower bound = flux / (2 * max_flux) (no nulls at
        // default tolerance here), analog pdf = flux / sum(flux).
        let t = MeshTallyData {
            tally_number: 1,
            particle: mcnp_io::meshtal::ParticleKind::Neutron,
            dose_response: false,
            x_bounds: vec![0.0, 1.0],
            y_bounds: vec![0.0, 1.0],
            z_bounds: vec![0.0, 1.0, 2.0, 3.0],
            e_bounds: vec![0.0, 1.0],
            column_idx: Default::default(),
            result: vec![vec![2.0], vec![1.0], vec![4.0]],
            rel_error: vec![vec![0.1], vec![0.1], vec![0.1]],
            total_result: vec![2.0, 1.0, 4.0],
            total_rel_error: vec![0.1, 0.1, 0.1],
        };
        let ww = magic(&t).unwrap();
        let s = MeshSourceSampler::new(&t, Mode::Analog, None).unwrap();
        let pdf = s.table().pdf();

        let ratio = ww.lower_bounds_ww[0] / pdf[0];
        for (&w, &p) in ww.lower_bounds_ww[1..].iter().zip(pdf[1..].iter()) {
            approx(w / p, ratio, 1e-12);
        }
        // Spot-check both normalizations independently.
        approx(pdf[2], 4.0 / 7.0, 1e-12);
        approx(ww.lower_bounds_ww[2], 0.5, 1e-12);
    }

    #[test]
    fn uniform_mode_weights_cells_by_volume() {
        // Two stacked x-cells of widths 1 and 3 with densities 10 and 1:
        // analog PDF ∝ strength (density × volume) = [10, 3],
        // uniform PDF ∝ volume = [1, 3].
        let t = MeshTallyData {
            tally_number: 1,
            particle: mcnp_io::meshtal::ParticleKind::Neutron,
            dose_response: false,
            x_bounds: vec![0.0, 1.0, 4.0],
            y_bounds: vec![0.0, 1.0],
            z_bounds: vec![0.0, 1.0],
            e_bounds: vec![0.0, 1.0],
            column_idx: Default::default(),
            result: vec![vec![10.0], vec![1.0]],
            rel_error: vec![vec![0.0], vec![0.0]],
            total_result: vec![10.0, 1.0],
            total_rel_error: vec![0.0, 0.0],
        };
        let uni = MeshSourceSampler::new(&t, Mode::Uniform, None).unwrap();
        let ana = MeshSourceSampler::new(&t, Mode::Analog, None).unwrap();

        let mut rng = Lcg(42_424_242);
        let (mut u0, mut a0) = (0u32, 0u32);
        for _ in 0..100_000 {
            if uni.sample(rng.next_f64(), rng.next_f64()).index == 0 {
                u0 += 1;
            }
            if ana.sample(rng.next_f64(), rng.next_f64()).index == 0 {
                a0 += 1;
            }
        }
        // Analog hits cell 0 w.p. 10/13; uniform only 1/4.
        approx(u0 as f64 / 100_000.0, 0.25, 0.02);
        approx(a0 as f64 / 100_000.0, 10.0 / 13.0, 0.02);

        // Biased birth weights = analog pdf / biased pdf:
        // narrow cell (10/13)/(1/4) = 40/13, wide cell (3/13)/(3/4) = 4/13.
        approx(uni.biased_weights[0], 40.0 / 13.0, 1e-12);
        approx(uni.biased_weights[1], 4.0 / 13.0, 1e-12);
        approx(ana.biased_weights[0], 1.0, 1e-12);
    }

    #[test]
    fn user_mode_respects_external_pdf_and_reports_weights() {
        let t = MeshTallyData {
            tally_number: 1,
            particle: mcnp_io::meshtal::ParticleKind::Photon,
            dose_response: false,
            x_bounds: vec![0.0, 1.0, 2.0],
            y_bounds: vec![0.0, 1.0],
            z_bounds: vec![0.0, 1.0],
            e_bounds: vec![0.0, 1.0],
            column_idx: Default::default(),
            result: vec![vec![10.0], vec![0.001]],
            rel_error: vec![vec![0.0], vec![0.0]],
            total_result: vec![10.0, 0.001],
            total_rel_error: vec![0.0, 0.0],
        };
        // Bias everything into the weakly-populated second cell.
        let s = MeshSourceSampler::new(&t, Mode::User, Some(&[0.0, 5.0])).unwrap();
        assert_eq!(s.mode(), Mode::User);
        assert_eq!(s.table().pdf(), &[0.0, 1.0]);

        let mut rng = Lcg(777_777);
        for _ in 0..1_000 {
            let sv = s.sample(rng.next_f64(), rng.next_f64());
            assert_eq!(sv.index, 1);
            assert_eq!((sv.i, sv.j, sv.k), (1, 0, 0));
            // Birth weight = analog pdf / biased pdf ≈ 9.999e-5: births are
            // concentrated where particles rarely stream, so each carries a
            // small weight.
            approx(sv.weight, 0.001 / 10.001, 1e-9);
        }
    }

    #[test]
    fn sampled_indices_decompose_consistently() {
        let t = fixture_tally("mcnp_meshtal_single_meshtal.txt", 4);
        let s = MeshSourceSampler::new(&t, Mode::Analog, None).unwrap();
        let d = t.dims();
        let mut rng = Lcg(20_260_824);
        for _ in 0..10_000 {
            let sv = s.sample(rng.next_f64(), rng.next_f64());
            assert!(sv.i < d[0] && sv.j < d[1] && sv.k < d[2], "{sv:?}");
            assert_eq!((sv.i * d[1] + sv.j) * d[2] + sv.k, sv.index);
            assert_eq!(t.ve_index(sv.i, sv.j, sv.k), sv.index);
        }
    }

    #[test]
    fn user_pdf_length_mismatch_rejected() {
        let t = fixture_tally("mcnp_meshtal_single_meshtal.txt", 4);
        assert_eq!(
            MeshSourceSampler::new(&t, Mode::User, Some(&[1.0; 3])),
            Err(Error::LengthMismatch {
                expected: 45,
                got: 3
            })
        );
    }

    #[test]
    fn all_zero_totals_rejected() {
        let mut t = fixture_tally("mcnp_meshtal_single_meshtal.txt", 4);
        t.total_result.iter_mut().for_each(|v| *v = 0.0);
        assert_eq!(
            MeshSourceSampler::new(&t, Mode::Analog, None),
            Err(Error::ZeroSumPdf)
        );
    }
}
