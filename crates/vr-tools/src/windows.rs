//! Weight-window emission for OpenMC and Serpent over [`MagicOutput`].
//!
//! Both emitters are pure formatting steps: they reorder
//! [`MagicOutput::lower_bounds_ww`] into the target code's flat layout and
//! wrap it with the mesh/energy structure taken from the source
//! [`MeshTallyData`]. No new math lives here.
//!
//! # Pinned format spellings (public docs)
//!
//! - **OpenMC `settings.xml`** — `<mesh>` element per the settings
//!   specification §3.29 (`id` attribute; `dimension`, `lower_left`,
//!   `upper_right` sub-elements) and `<weight_windows>` per §3.66, spelled
//!   exactly as the OpenMC Python writer emits and the C++ XML reader parses
//!   it: `id` attribute plus `mesh`, `particle_type`, `energy_bounds` (eV),
//!   `lower_ww_bounds`, `upper_ww_bounds`, `survival_ratio`, `max_split`, and
//!   `weight_cutoff` sub-elements. Flat bound order follows the C++ bounds
//!   tensor layout `(energy_bin, mesh_bin)` with mesh bins x-fastest — i.e.
//!   per energy group, all volume elements in x-slowest/z-fastest mesh order
//!   transposed to x-fastest order, groups outermost. The emitted fragment is
//!   pasted inside the existing `<settings>` root next to the user's other
//!   elements.
//! - **Serpent `wwin`** — the user guide §2.2.8.3 and the input-syntax manual
//!   (`wwin` card) define `wwin NAME wf FILE FMT` with `FMT = 1` for a
//!   Serpent-generated mesh and `FMT = 2` for the MCNP WWINP text format.
//!   Only FMT = 2 is emitted: its layout is pinned by public documentation
//!   (MCNP user manual) and by the in-workspace reader in
//!   [`nucleide_mcnp_io::wwinp`], which the tests re-parse through. The
//!   Serpent-native FMT = 1 layout is not publicly specified; requesting it
//!   is a loud error (there is deliberately no API that emits it), and the
//!   returned [`SerpentWwin::card`] pins `wf "<file>" 2`.
//!
//! # Loud-error boundary
//!
//! A [`MagicOutput`] has no well-formed spelling when any lower bound is
//! non-finite or negative (both formats treat non-positive windows as inert,
//! but a negative bound has no meaning), when the energy upper bounds are not
//! strictly increasing positive values, when the mesh bounds are degenerate,
//! or when the requested OpenMC tuning parameters fall outside the ranges
//! the OpenMC reader enforces. The Serpent `wwin` card tokens (`name`,
//! `file`) must be single whitespace-and-quote-free words, and the OpenMC
//! scaling steps (MeV→eV energies, `lower × upper_bound_ratio`) must stay
//! finite after formatting. Each case is a named [`crate::Error`]
//! variant; nothing is silently dropped or clamped.

use nucleide_mcnp_io::meshtal::{MeshTallyData, ParticleKind};
use nucleide_mcnp_io::wwinp::Wwinp;

use crate::{Error, MagicOutput, Result};

/// OpenMC emission settings, defaulting to the values OpenMC itself uses.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OpenMcOptions {
    /// ID written to the `<mesh id="...">` attribute (referenced by
    /// `<weight_windows><mesh>`).
    pub mesh_id: u32,
    /// ID written to the `<weight_windows id="...">` attribute.
    pub window_id: u32,
    /// Ratio of upper to lower weight-window bounds. MAGIC produces lower
    /// bounds only; OpenMC requires both, so the uppers are synthesized as
    /// `upper = lower * ratio`. OpenMC's own WWINP importer uses 5.0.
    pub upper_bound_ratio: f64,
    /// Survival weight over lower bound; must be greater than 1 and less
    /// than [`OpenMcOptions::upper_bound_ratio`] (OpenMC reader range).
    pub survival_ratio: f64,
    /// Maximum split count; must be at least 2 (OpenMC reader range).
    pub max_split: u32,
    /// Russian-roulette weight cutoff; must be in `(0, 1]` (OpenMC reader
    /// range).
    pub weight_cutoff: f64,
}

impl Default for OpenMcOptions {
    fn default() -> Self {
        Self {
            mesh_id: 1,
            window_id: 1,
            upper_bound_ratio: 5.0,
            survival_ratio: 3.0,
            max_split: 10,
            weight_cutoff: 1.0e-38,
        }
    }
}

/// Emitted OpenMC weight-window fragment.
#[derive(Debug, Clone, PartialEq)]
pub struct OpenMcWeightWindows {
    /// `<mesh>` + `<weight_windows>` elements to paste into `settings.xml`.
    pub xml: String,
    /// Machine-readable drift notes (synthesized structure, null cells).
    pub notes: Vec<String>,
}

/// Emitted Serpent weight-window file (MCNP WWINP text, read via `wf ... 2`).
#[derive(Debug, Clone, PartialEq)]
pub struct SerpentWwin {
    /// WWINP-format file content (the `.wwd` file body).
    pub text: String,
    /// Input card referencing the file: `wwin <name> wf "<file>" 2`.
    pub card: String,
    /// Machine-readable drift notes (format pinning, synthesized structure).
    pub notes: Vec<String>,
}

/// Emit MAGIC lower bounds as an OpenMC `settings.xml` fragment.
pub fn emit_openmc_weight_windows(
    output: &MagicOutput,
    tally: &MeshTallyData,
    options: &OpenMcOptions,
) -> Result<OpenMcWeightWindows> {
    let dims = validate(output, tally)?;
    validate_openmc_options(options)?;

    let groups = output.groups_per_ve;
    let particle = match tally.particle {
        ParticleKind::Neutron => "neutron",
        ParticleKind::Photon => "photon",
    };
    // Energy bounds in eV with the implicit zero lower edge, matching the
    // convention OpenMC's own WWINP importer applies to upper-bound-only
    // inputs. Meshtal energies are MeV.
    let energy_bounds: Vec<f64> = std::iter::once(0.0)
        .chain(output.e_upper_bounds.iter().map(|e| e * 1.0e6))
        .collect();
    let lower = flat_bounds_xfastest(output, dims);
    let upper: Vec<f64> = lower
        .iter()
        .map(|v| v * options.upper_bound_ratio)
        .collect();
    // Finite validated inputs can still overflow at the scaling steps
    // (MeV→eV, lower × ratio); a non-finite value in the XML would silently
    // break the OpenMC reader, so the emission is loud instead.
    if let Some(&v) = energy_bounds.iter().find(|v| !v.is_finite()) {
        return Err(Error::BadEmissionOption {
            option: "energy_bounds",
            value: v.to_string(),
            detail: "MeV-to-eV scaling overflowed to a non-finite value",
        });
    }
    if let Some(&v) = upper.iter().find(|v| !v.is_finite()) {
        return Err(Error::BadEmissionOption {
            option: "upper_ww_bounds",
            value: v.to_string(),
            detail: "lower * upper_bound_ratio overflowed to a non-finite value",
        });
    }

    let mut xml = String::new();
    xml += &format!("  <mesh id=\"{}\">\n", options.mesh_id);
    xml += &format!(
        "    <dimension>{} {} {}</dimension>\n",
        dims[0], dims[1], dims[2]
    );
    xml += &format!(
        "    <lower_left>{} {} {}</lower_left>\n",
        xml_f64(tally.x_bounds[0]),
        xml_f64(tally.y_bounds[0]),
        xml_f64(tally.z_bounds[0])
    );
    xml += &format!(
        "    <upper_right>{} {} {}</upper_right>\n",
        xml_f64(*tally.x_bounds.last().expect("validated")),
        xml_f64(*tally.y_bounds.last().expect("validated")),
        xml_f64(*tally.z_bounds.last().expect("validated"))
    );
    xml += "  </mesh>\n";
    xml += &format!("  <weight_windows id=\"{}\">\n", options.window_id);
    xml += &format!("    <mesh>{}</mesh>\n", options.mesh_id);
    xml += &format!("    <particle_type>{particle}</particle_type>\n");
    xml += &format!(
        "    <energy_bounds>{}</energy_bounds>\n",
        join_f64(&energy_bounds)
    );
    xml += &format!(
        "    <lower_ww_bounds>{}</lower_ww_bounds>\n",
        join_f64(&lower)
    );
    xml += &format!(
        "    <upper_ww_bounds>{}</upper_ww_bounds>\n",
        join_f64(&upper)
    );
    xml += &format!(
        "    <survival_ratio>{}</survival_ratio>\n",
        xml_f64(options.survival_ratio)
    );
    xml += &format!("    <max_split>{}</max_split>\n", options.max_split);
    xml += &format!(
        "    <weight_cutoff>{}</weight_cutoff>\n",
        xml_f64(options.weight_cutoff)
    );
    xml += "  </weight_windows>\n";

    let mut notes = vec![format!(
        "upper-bounds-synthesized: upper = lower * {} (MAGIC produces lower bounds only; OpenMC requires both)",
        options.upper_bound_ratio
    )];
    if groups > 1 {
        notes.push(format!(
            "energy-bounds-mev-to-ev: {} group upper bounds converted MeV -> eV with implicit 0 eV lower edge",
            groups
        ));
    }
    notes.extend(null_notes(output));

    Ok(OpenMcWeightWindows { xml, notes })
}

/// Validate one token interpolated into the Serpent `wwin` card. The card
/// is assembled as `wwin {name} wf "{file}" 2`, so a token holding
/// whitespace or a quote would corrupt the deck (an injected newline starts
/// a new card). Only a single whitespace-and-quote-free word is accepted.
fn check_wwin_token(option: &'static str, token: &str) -> Result<()> {
    if token.is_empty()
        || token
            .chars()
            .any(|c| c.is_whitespace() || c == '"' || c == '\'')
    {
        return Err(Error::BadEmissionOption {
            option,
            value: token.to_string(),
            detail: "must be a single token without whitespace or quotes",
        });
    }
    Ok(())
}

/// Emit MAGIC lower bounds as a Serpent-readable weight-window file.
///
/// The file is written in the MCNP WWINP text spelling, which Serpent reads
/// through `wwin <name> wf "<file>" 2` (user guide §2.2.8.3). The
/// Serpent-native FMT = 1 layout is not publicly documented and is never
/// emitted. WWINP carries only window lower bounds and the energies/mesh in
/// MeV/cm, so no unit conversion is applied.
pub fn emit_serpent_wwin(
    output: &MagicOutput,
    tally: &MeshTallyData,
    name: &str,
    file: &str,
) -> Result<SerpentWwin> {
    check_wwin_token("name", name)?;
    check_wwin_token("file", file)?;
    let dims = validate(output, tally)?;
    let groups = output.groups_per_ve;

    let mut cm = Vec::with_capacity(3);
    let mut fm = Vec::with_capacity(3);
    let mut bounds = Vec::with_capacity(3);
    for axis in [&tally.x_bounds, &tally.y_bounds, &tally.z_bounds] {
        let (c, f) = coarse_fine(axis);
        cm.push(c);
        fm.push(f);
        bounds.push(axis.clone());
    }
    let nc = [cm[0].len() as u32, cm[1].len() as u32, cm[2].len() as u32];
    let nf = [dims[0] as u32, dims[1] as u32, dims[2] as u32];
    let nft = nf.iter().map(|v| u64::from(*v)).product();

    let nve = dims[0] * dims[1] * dims[2];
    let mut particle_windows: Vec<Vec<f64>> = vec![vec![0.0; nve]; groups];
    for_each_bound_xfastest(output, dims, |g, ve, w| particle_windows[g][ve] = w);
    let energies = output.e_upper_bounds.clone();

    // MCNP convention (see the nucleide-mcnp-io WWINP fixtures): neutron-only
    // files set ni = 1 with a single energy-group slot; photon-only files set
    // ni = 2 and declare an empty neutron slot, with photon data in the
    // second slot.
    let (ni, ne, e, ww) = match tally.particle {
        ParticleKind::Neutron => (
            1,
            vec![groups as u32],
            vec![energies],
            vec![particle_windows],
        ),
        ParticleKind::Photon => (
            2,
            vec![0, groups as u32],
            vec![energies],
            vec![particle_windows],
        ),
    };

    let wwinp = Wwinp {
        ni,
        nr: 10,
        ne,
        nf,
        nft,
        origin: [tally.x_bounds[0], tally.y_bounds[0], tally.z_bounds[0]],
        nc,
        nwg: 1,
        date_time: String::new(),
        cm,
        fm,
        bounds,
        e,
        ww,
    };
    let text = wwinp.to_text().map_err(|e| Error::Wwinp(e.to_string()))?;

    let mut notes = vec![
        "format-pinned: Serpent reads this file via `wwin <name> wf \"<file>\" 2` (MCNP WWINP \
         spelling, user guide 2.2.8.3); the Serpent-native FMT=1 wwd layout is not publicly \
         documented and is not emitted"
            .to_string(),
        "lower-bounds-only: WWINP carries lower bounds; Serpent derives upper bounds from \
         importances at load time (set wwb LB=0.5 UB=2 defaults)"
            .to_string(),
    ];
    notes.extend(null_notes(output));

    Ok(SerpentWwin {
        text,
        card: format!("wwin {name} wf \"{file}\" 2"),
        notes,
    })
}

/// Validate the shared input structure; returns `[nx, ny, nz]` on success.
fn validate(output: &MagicOutput, tally: &MeshTallyData) -> Result<[usize; 3]> {
    let axes = [&tally.x_bounds, &tally.y_bounds, &tally.z_bounds];
    for (axis, b) in axes.iter().enumerate() {
        if b.len() < 2 {
            return Err(Error::BadMeshBounds { axis, index: 0 });
        }
        for (i, &v) in b.iter().enumerate() {
            if !v.is_finite() {
                return Err(Error::BadMeshBounds { axis, index: i });
            }
            if i > 0 && v <= b[i - 1] {
                return Err(Error::BadMeshBounds { axis, index: i });
            }
        }
    }
    let dims = [
        tally.x_bounds.len() - 1,
        tally.y_bounds.len() - 1,
        tally.z_bounds.len() - 1,
    ];
    let groups = output.groups_per_ve;
    if groups == 0 || output.e_upper_bounds.len() != groups {
        return Err(Error::BadEnergyBounds { index: 0 });
    }
    for (g, &e) in output.e_upper_bounds.iter().enumerate() {
        if !e.is_finite() || e <= 0.0 {
            return Err(Error::BadEnergyBounds { index: g });
        }
        if g > 0 && e <= output.e_upper_bounds[g - 1] {
            return Err(Error::BadEnergyBounds { index: g });
        }
    }
    let expected = dims[0] * dims[1] * dims[2] * groups;
    if output.lower_bounds_ww.len() != expected {
        return Err(Error::LengthMismatch {
            expected,
            got: output.lower_bounds_ww.len(),
        });
    }
    for (i, &w) in output.lower_bounds_ww.iter().enumerate() {
        if !w.is_finite() {
            return Err(Error::NonFiniteWindow { index: i });
        }
        if w < 0.0 {
            return Err(Error::NegativeWindow { index: i, value: w });
        }
    }
    Ok(dims)
}

fn validate_openmc_options(options: &OpenMcOptions) -> Result<()> {
    if options.upper_bound_ratio <= 1.0 {
        return Err(Error::BadEmissionOption {
            option: "upper_bound_ratio",
            value: options.upper_bound_ratio.to_string(),
            detail: "must be greater than 1",
        });
    }
    if !(options.survival_ratio > 1.0 && options.survival_ratio < options.upper_bound_ratio) {
        return Err(Error::BadEmissionOption {
            option: "survival_ratio",
            value: options.survival_ratio.to_string(),
            detail: "must be greater than 1 and less than upper_bound_ratio",
        });
    }
    if options.max_split < 2 {
        return Err(Error::BadEmissionOption {
            option: "max_split",
            value: options.max_split.to_string(),
            detail: "must be at least 2",
        });
    }
    if !(options.weight_cutoff > 0.0 && options.weight_cutoff <= 1.0) {
        return Err(Error::BadEmissionOption {
            option: "weight_cutoff",
            value: options.weight_cutoff.to_string(),
            detail: "must be in (0, 1]",
        });
    }
    Ok(())
}

/// Visit every `(group, x-fastest ve index, bound)` triple in the flat
/// x-fastest window layout both target codes use: group outermost, then k,
/// j, i with x fastest. [`MagicOutput::lower_bounds_ww`] stores the
/// transposed (z-fastest) layout with groups innermost, so this walks it in
/// one sequential pass (volume-element rows are contiguous) and scatters
/// each row across the `groups` destination regions.
fn for_each_bound_xfastest(
    output: &MagicOutput,
    dims: [usize; 3],
    mut emit: impl FnMut(usize, usize, f64),
) {
    let (nx, ny, nz) = (dims[0], dims[1], dims[2]);
    let groups = output.groups_per_ve;
    for i in 0..nx {
        for j in 0..ny {
            for k in 0..nz {
                let row = &output.lower_bounds_ww[((i * ny + j) * nz + k) * groups..][..groups];
                let ve = (k * ny + j) * nx + i;
                for (g, &w) in row.iter().enumerate() {
                    emit(g, ve, w);
                }
            }
        }
    }
}

/// Flat x-fastest window vector for all energy groups: group `g` outermost,
/// then k, j, i with x fastest — the ordering both OpenMC and WWINP use.
fn flat_bounds_xfastest(output: &MagicOutput, dims: [usize; 3]) -> Vec<f64> {
    let nve = dims[0] * dims[1] * dims[2];
    let mut flat = vec![0.0; nve * output.groups_per_ve];
    for_each_bound_xfastest(output, dims, |g, ve, w| flat[g * nve + ve] = w);
    flat
}

/// Group equal-width runs of cells into WWINP coarse bins. WWINP block 2
/// stores only coarse boundaries plus a fine-bin count, with fine bins
/// defined by uniform interpolation, so cells may share a coarse bin only
/// when their widths are exactly equal (bitwise) — anything else gets its
/// own coarse bin and reproduces the grid exactly.
fn coarse_fine(bounds: &[f64]) -> (Vec<f64>, Vec<f64>) {
    debug_assert!(bounds.len() >= 2);
    let mut cm = Vec::new();
    let mut fm = Vec::new();
    let mut run_start = 0usize;
    let mut run_w = bounds[1] - bounds[0];
    for k in 1..bounds.len() - 1 {
        let w = bounds[k + 1] - bounds[k];
        if w != run_w {
            cm.push(bounds[k]);
            fm.push((k - run_start) as f64);
            run_start = k;
            run_w = w;
        }
    }
    cm.push(bounds[bounds.len() - 1]);
    fm.push((bounds.len() - 1 - run_start) as f64);
    (cm, fm)
}

fn null_notes(output: &MagicOutput) -> Vec<String> {
    let nulled = output.lower_bounds_ww.iter().filter(|&&w| w == 0.0).count();
    if nulled == 0 {
        Vec::new()
    } else {
        vec![format!(
            "null-cells: {nulled} of {} lower bounds are 0.0 (MAGIC null value); non-positive \
             windows are inert in OpenMC and Serpent treats the cell importance as infinite \
             (set wwb) — review whether null_value=0 is intended",
            output.lower_bounds_ww.len()
        )]
    }
}

/// Format a float for XML text nodes: shortest round-trip decimal, switching
/// to scientific notation outside `[1e-6, 1e21)` so tiny values like the
/// default `1e-38` weight cutoff stay compact. Any spelling `std::stod`
/// (OpenMC) parses is acceptable; fixtures pin this one.
fn xml_f64(v: f64) -> String {
    XmlF64(v).to_string()
}

/// [`std::fmt::Display`] adapter spelling [`xml_f64`] without an
/// intermediate heap `String` per value.
struct XmlF64(f64);

impl std::fmt::Display for XmlF64 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let v = self.0;
        if v == 0.0 {
            return f.write_str("0.0");
        }
        if (1.0e-6..1.0e21).contains(&v.abs()) {
            write!(f, "{v}")
        } else {
            write!(f, "{v:e}")
        }
    }
}

fn join_f64(values: &[f64]) -> String {
    use std::fmt::Write as _;
    let mut out = String::with_capacity(values.len() * 4);
    for (i, &v) in values.iter().enumerate() {
        if i > 0 {
            out.push(' ');
        }
        write!(out, "{}", XmlF64(v)).expect("write to String");
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::magic::MagicOutput;

    /// Synthetic 2×2×1 neutron tally with uniform 1.0/2.0/2.0 cm cells.
    fn sample_tally() -> MeshTallyData {
        MeshTallyData {
            tally_number: 4,
            particle: ParticleKind::Neutron,
            dose_response: false,
            x_bounds: vec![0.0, 1.0, 2.0],
            y_bounds: vec![-1.0, 1.0, 3.0],
            z_bounds: vec![10.0, 12.0],
            e_bounds: vec![0.0, 0.5, 1.0],
            column_idx: Default::default(),
            result: Vec::new(),
            rel_error: Vec::new(),
            total_result: Vec::new(),
            total_rel_error: Vec::new(),
        }
    }

    /// Two-group windows, values chosen dyadic so `{:13.5E}` round-trips
    /// exactly through the WWINP reader and `upper = 5*lower` is exact.
    /// `lower_bounds_ww` is ve-major (z-fastest ve, groups innermost).
    fn sample_output() -> MagicOutput {
        MagicOutput {
            lower_bounds_ww: vec![
                0.5, 0.25, // ve0 = (i0, j0, k0)
                0.25, 0.5, // ve1 = (i0, j1, k0)
                0.125, 0.75, // ve2 = (i1, j0, k0)
                0.0625, 1.0, // ve3 = (i1, j1, k0)
            ],
            groups_per_ve: 2,
            scale_factors: vec![1.0, 1.0],
            e_upper_bounds: vec![0.5, 1.0],
            ww_tag_name: "ww_n".to_string(),
            e_upper_bounds_tag_name: "n_e_upper_bounds".to_string(),
        }
    }

    /// The x-fastest flat layout both target codes use, group-outermost.
    const XFASTEST_FLAT: [f64; 8] = [0.5, 0.125, 0.25, 0.0625, 0.25, 0.75, 0.5, 1.0];

    fn floats_of(tag: &str, xml: &str) -> Vec<f64> {
        let open = format!("<{tag}>");
        let line = xml
            .lines()
            .find(|l| l.trim_start().starts_with(&open))
            .unwrap_or_else(|| panic!("missing <{tag}> in:\n{xml}"));
        line.trim()
            .strip_prefix(&open)
            .and_then(|s| s.strip_suffix(&format!("</{tag}>")))
            .expect("malformed element")
            .split_whitespace()
            .map(|t| t.parse::<f64>().unwrap())
            .collect()
    }

    #[test]
    fn openmc_per_group_matches_golden_xml() {
        let out = emit_openmc_weight_windows(
            &sample_output(),
            &sample_tally(),
            &OpenMcOptions::default(),
        )
        .unwrap();
        let expected = "  <mesh id=\"1\">\n\
            \x20   <dimension>2 2 1</dimension>\n\
            \x20   <lower_left>0.0 -1 10</lower_left>\n\
            \x20   <upper_right>2 3 12</upper_right>\n\
            \x20 </mesh>\n\
            \x20 <weight_windows id=\"1\">\n\
            \x20   <mesh>1</mesh>\n\
            \x20   <particle_type>neutron</particle_type>\n\
            \x20   <energy_bounds>0.0 500000 1000000</energy_bounds>\n\
            \x20   <lower_ww_bounds>0.5 0.125 0.25 0.0625 0.25 0.75 0.5 1</lower_ww_bounds>\n\
            \x20   <upper_ww_bounds>2.5 0.625 1.25 0.3125 1.25 3.75 2.5 5</upper_ww_bounds>\n\
            \x20   <survival_ratio>3</survival_ratio>\n\
            \x20   <max_split>10</max_split>\n\
            \x20   <weight_cutoff>1e-38</weight_cutoff>\n\
            \x20 </weight_windows>\n";
        assert_eq!(out.xml, expected);
    }

    #[test]
    fn openmc_reparse_recovers_windows_and_structure() {
        let out = emit_openmc_weight_windows(
            &sample_output(),
            &sample_tally(),
            &OpenMcOptions::default(),
        )
        .unwrap();
        // Structural assertions on our own emitted text (test-local parsing,
        // no new reader crate).
        assert_eq!(
            floats_of("lower_ww_bounds", &out.xml),
            XFASTEST_FLAT.to_vec()
        );
        let upper: Vec<f64> = XFASTEST_FLAT.iter().map(|v| v * 5.0).collect();
        assert_eq!(floats_of("upper_ww_bounds", &out.xml), upper);
        assert_eq!(
            floats_of("energy_bounds", &out.xml),
            vec![0.0, 500000.0, 1.0e6]
        );
        assert_eq!(floats_of("dimension", &out.xml), vec![2.0, 2.0, 1.0]);
        assert_eq!(floats_of("lower_left", &out.xml), vec![0.0, -1.0, 10.0]);
        assert_eq!(floats_of("upper_right", &out.xml), vec![2.0, 3.0, 12.0]);
        assert!(out.xml.contains("<mesh id=\"1\">"));
        assert!(out.xml.contains("<weight_windows id=\"1\">"));
        assert!(out.xml.contains("<particle_type>neutron</particle_type>"));
        assert!(out
            .notes
            .iter()
            .any(|n| n.starts_with("upper-bounds-synthesized")));
    }

    #[test]
    fn openmc_total_mode_emits_single_energy_bin() {
        let mut output = sample_output();
        output.groups_per_ve = 1;
        output.lower_bounds_ww = vec![0.5, 0.25, 0.125, 0.0625];
        output.e_upper_bounds = vec![1.0];
        let out = emit_openmc_weight_windows(&output, &sample_tally(), &OpenMcOptions::default())
            .unwrap();
        assert_eq!(floats_of("energy_bounds", &out.xml), vec![0.0, 1.0e6]);
        assert_eq!(
            floats_of("lower_ww_bounds", &out.xml),
            vec![0.5, 0.125, 0.25, 0.0625]
        );
    }

    #[test]
    fn openmc_post_scaling_overflow_is_loud() {
        // Finite, validated inputs can still overflow at the scaling steps:
        // huge MeV group bounds x1e6 and huge lower bounds x the ratio.
        let mut output = sample_output();
        output.e_upper_bounds = vec![1.0e300, 1.0e308];
        assert!(matches!(
            emit_openmc_weight_windows(&output, &sample_tally(), &OpenMcOptions::default()),
            Err(Error::BadEmissionOption {
                option: "energy_bounds",
                ..
            })
        ));
        let mut output = sample_output();
        output.lower_bounds_ww = vec![1.0e308; 8];
        assert!(matches!(
            emit_openmc_weight_windows(&output, &sample_tally(), &OpenMcOptions::default()),
            Err(Error::BadEmissionOption {
                option: "upper_ww_bounds",
                ..
            })
        ));
    }

    #[test]
    fn serpent_wwin_rejects_card_token_injection() {
        // `name`/`file` interpolate into `wwin {name} wf "{file}" 2`; a
        // whitespace or quote character would corrupt the deck.
        for (name, file, option) in [
            ("ww1\nsrc 2 1 1 1", "mesh.wwd", "name"),
            ("ww1 src", "mesh.wwd", "name"),
            ("ww1", "mesh.wwd\nsrc 2", "file"),
            ("ww1", "mesh\"wwd", "file"),
        ] {
            let err = emit_serpent_wwin(&sample_output(), &sample_tally(), name, file).unwrap_err();
            assert!(
                matches!(err, Error::BadEmissionOption { option: o, .. } if o == option),
                "{name:?}/{file:?}: {err}"
            );
        }
        // Plain tokens still emit the pinned card.
        let out = emit_serpent_wwin(&sample_output(), &sample_tally(), "ww1", "mesh.wwd").unwrap();
        assert_eq!(out.card, "wwin ww1 wf \"mesh.wwd\" 2");
    }

    #[test]
    fn openmc_photon_spells_photon_type() {
        let mut tally = sample_tally();
        tally.particle = ParticleKind::Photon;
        let out = emit_openmc_weight_windows(&sample_output(), &tally, &OpenMcOptions::default())
            .unwrap();
        assert!(out.xml.contains("<particle_type>photon</particle_type>"));
    }

    #[test]
    fn openmc_option_validation_is_loud() {
        let bad = OpenMcOptions {
            survival_ratio: 1.0,
            ..Default::default()
        };
        assert!(matches!(
            emit_openmc_weight_windows(&sample_output(), &sample_tally(), &bad),
            Err(Error::BadEmissionOption {
                option: "survival_ratio",
                ..
            })
        ));
        let bad = OpenMcOptions {
            survival_ratio: 6.0,
            ..Default::default()
        };
        assert!(matches!(
            emit_openmc_weight_windows(&sample_output(), &sample_tally(), &bad),
            Err(Error::BadEmissionOption {
                option: "survival_ratio",
                ..
            })
        ));
        let bad = OpenMcOptions {
            max_split: 1,
            ..Default::default()
        };
        assert!(matches!(
            emit_openmc_weight_windows(&sample_output(), &sample_tally(), &bad),
            Err(Error::BadEmissionOption {
                option: "max_split",
                ..
            })
        ));
        let bad = OpenMcOptions {
            weight_cutoff: 2.0,
            ..Default::default()
        };
        assert!(matches!(
            emit_openmc_weight_windows(&sample_output(), &sample_tally(), &bad),
            Err(Error::BadEmissionOption {
                option: "weight_cutoff",
                ..
            })
        ));
    }

    #[test]
    fn serpent_matches_golden_wwinp_text() {
        let out = emit_serpent_wwin(&sample_output(), &sample_tally(), "ww1", "mesh.wwd").unwrap();
        assert_eq!(out.card, "wwin ww1 wf \"mesh.wwd\" 2");
        // The emitter must reproduce exactly what assembling the same Wwinp
        // through the canonical mcnp-io writer produces.
        let expected = Wwinp {
            ni: 1,
            nr: 10,
            ne: vec![2],
            nf: [2, 2, 1],
            nft: 4,
            origin: [0.0, -1.0, 10.0],
            nc: [1, 1, 1],
            nwg: 1,
            date_time: String::new(),
            cm: vec![vec![2.0], vec![3.0], vec![12.0]],
            fm: vec![vec![2.0], vec![2.0], vec![1.0]],
            bounds: vec![vec![0.0, 1.0, 2.0], vec![-1.0, 1.0, 3.0], vec![10.0, 12.0]],
            e: vec![vec![0.5, 1.0]],
            ww: vec![vec![
                vec![0.5, 0.125, 0.25, 0.0625],
                vec![0.25, 0.75, 0.5, 1.0],
            ]],
        };
        assert_eq!(out.text, expected.to_text().unwrap());
    }

    #[test]
    fn serpent_reparse_through_wwinp_reader_is_exact() {
        let out = emit_serpent_wwin(&sample_output(), &sample_tally(), "ww1", "mesh.wwd").unwrap();
        // Real reader round trip: bounds, energies, and every window value
        // must come back exactly (fixture values are 13.5E-exact).
        let w = Wwinp::parse(&out.text).unwrap();
        assert_eq!(w.ni, 1);
        assert_eq!(w.nr, 10);
        assert_eq!(w.ne, vec![2]);
        assert_eq!(w.nf, [2, 2, 1]);
        assert_eq!(w.nft, 4);
        assert_eq!(w.origin, [0.0, -1.0, 10.0]);
        assert_eq!(
            w.bounds,
            vec![vec![0.0, 1.0, 2.0], vec![-1.0, 1.0, 3.0], vec![10.0, 12.0]]
        );
        assert_eq!(w.e, vec![vec![0.5, 1.0]]);
        // x-fastest order, group rows: this is the transpose-sensitive check.
        assert_eq!(w.ww[0][0], vec![0.5, 0.125, 0.25, 0.0625]);
        assert_eq!(w.ww[0][1], vec![0.25, 0.75, 0.5, 1.0]);
    }

    #[test]
    fn serpent_nonuniform_mesh_reproduces_bounds() {
        let mut tally = sample_tally();
        // Widths 1, 2, 1 with no exactly-uniform run beyond single cells:
        // every cell becomes its own coarse bin and the grid is reproduced
        // exactly by the block-2 stream.
        tally.x_bounds = vec![0.0, 1.0, 3.0, 4.0];
        let mut output = sample_output();
        output.lower_bounds_ww = vec![
            0.5, 0.25, 0.25, 0.5, 0.125, 0.75, 0.0625, 1.0, 0.5, 0.25, 0.25, 0.5,
        ];
        output.e_upper_bounds = vec![0.5, 1.0];
        let out = emit_serpent_wwin(&output, &tally, "ww1", "m.wwd").unwrap();
        let w = Wwinp::parse(&out.text).unwrap();
        assert_eq!(w.nf, [3, 2, 1]);
        assert_eq!(w.nc, [3, 1, 1]);
        assert_eq!(w.fm[0], vec![1.0, 1.0, 1.0]);
        assert_eq!(w.cm[0], vec![1.0, 3.0, 4.0]);
        assert_eq!(w.bounds[0], vec![0.0, 1.0, 3.0, 4.0]);
        assert_eq!(w.ww[0][0].len(), 6);
    }

    #[test]
    fn serpent_photon_uses_two_slot_header() {
        let mut tally = sample_tally();
        tally.particle = ParticleKind::Photon;
        let out = emit_serpent_wwin(&sample_output(), &tally, "ww1", "m.wwd").unwrap();
        let w = Wwinp::parse(&out.text).unwrap();
        assert_eq!(w.ni, 2);
        assert_eq!(w.ne, vec![0, 2]);
        assert_eq!(w.ww.len(), 1);
        assert_eq!(w.ww[0][0], vec![0.5, 0.125, 0.25, 0.0625]);
        assert!(out.notes.iter().any(|n| n.starts_with("format-pinned")));
    }

    #[test]
    fn serpent_total_mode_single_group() {
        let mut output = sample_output();
        output.groups_per_ve = 1;
        output.lower_bounds_ww = vec![0.5, 0.25, 0.125, 0.0625];
        output.e_upper_bounds = vec![1.0];
        let out = emit_serpent_wwin(&output, &sample_tally(), "ww1", "m.wwd").unwrap();
        let w = Wwinp::parse(&out.text).unwrap();
        assert_eq!(w.ne, vec![1]);
        assert_eq!(w.e, vec![vec![1.0]]);
        assert_eq!(w.ww[0].len(), 1);
    }

    #[test]
    fn magic_end_to_end_emission_recovers_fixture_windows() {
        // The full pipeline on the real meshtal fixture: MAGIC per-group →
        // both emitters → reader re-parse of the Serpent file and flat-order
        // spot checks of the OpenMC fragment.
        let m = nucleide_mcnp_io::meshtal::Meshtal::from_file(format!(
            "{}/../../fixtures/mcnp/meshtal/mcnp_meshtal_single_meshtal.txt",
            env!("CARGO_MANIFEST_DIR")
        ))
        .unwrap();
        let t = &m.tallies[&4];
        let output = crate::magic_with(
            t,
            crate::MagicSelection::PerGroup,
            crate::MagicParams::default(),
        )
        .unwrap();

        let serpent = emit_serpent_wwin(&output, t, "ww1", "mesh.wwd").unwrap();
        let w = Wwinp::parse(&serpent.text).unwrap();
        let dims = t.dims();
        assert_eq!(dims, [3, 5, 3]);
        assert_eq!(w.nf, [3, 5, 3]);
        assert_eq!(w.nft, 45);
        assert_eq!(w.ne, vec![3]);
        assert_eq!(w.e[0], output.e_upper_bounds);
        // Cell 22 (z-fastest ve) holds the max of every group; decompose
        // ve 22 = (i*ny + j)*nz + k and recompose x-fastest.
        let (i, j, k) = (
            22 / (dims[1] * dims[2]),
            (22 / dims[2]) % dims[1],
            22 % dims[2],
        );
        let ve_w = (k * dims[1] + j) * dims[0] + i;
        assert_eq!(w.ww[0].len(), 3);
        for group in &w.ww[0] {
            assert!((group[ve_w] - 0.5).abs() < 1e-12);
        }

        let openmc = emit_openmc_weight_windows(&output, t, &OpenMcOptions::default()).unwrap();
        let lower = floats_of("lower_ww_bounds", &openmc.xml);
        assert_eq!(lower.len(), 45 * 3);
        // Same cell through the OpenMC fragment (group-outermost flat list).
        for g in 0..3 {
            assert!((lower[g * 45 + ve_w] - 0.5).abs() < 1e-12);
        }
    }

    #[test]
    fn negative_window_is_rejected() {
        let mut output = sample_output();
        output.lower_bounds_ww[3] = -0.25;
        assert!(matches!(
            emit_serpent_wwin(&output, &sample_tally(), "ww1", "m.wwd"),
            Err(Error::NegativeWindow {
                index: 3,
                value: -0.25
            })
        ));
        assert!(matches!(
            emit_openmc_weight_windows(&output, &sample_tally(), &OpenMcOptions::default()),
            Err(Error::NegativeWindow {
                index: 3,
                value: -0.25
            })
        ));
    }

    #[test]
    fn non_finite_window_is_rejected() {
        let mut output = sample_output();
        output.lower_bounds_ww[0] = f64::NAN;
        assert!(matches!(
            emit_serpent_wwin(&output, &sample_tally(), "ww1", "m.wwd"),
            Err(Error::NonFiniteWindow { index: 0 })
        ));
    }

    #[test]
    fn bad_energy_bounds_are_rejected() {
        let mut output = sample_output();
        output.e_upper_bounds = vec![1.0, 0.5];
        assert!(matches!(
            emit_serpent_wwin(&output, &sample_tally(), "ww1", "m.wwd"),
            Err(Error::BadEnergyBounds { index: 1 })
        ));
        output.e_upper_bounds = vec![0.0, 1.0];
        assert!(matches!(
            emit_serpent_wwin(&output, &sample_tally(), "ww1", "m.wwd"),
            Err(Error::BadEnergyBounds { index: 0 })
        ));
        output.e_upper_bounds = vec![f64::INFINITY, 1.0];
        assert!(matches!(
            emit_openmc_weight_windows(&output, &sample_tally(), &OpenMcOptions::default()),
            Err(Error::BadEnergyBounds { index: 0 })
        ));
    }

    #[test]
    fn bad_mesh_bounds_are_rejected() {
        let mut tally = sample_tally();
        tally.x_bounds = vec![2.0, 1.0];
        assert!(matches!(
            emit_serpent_wwin(&sample_output(), &tally, "ww1", "m.wwd"),
            Err(Error::BadMeshBounds { axis: 0, index: 1 })
        ));
        tally.x_bounds = vec![0.0, f64::NAN];
        assert!(matches!(
            emit_serpent_wwin(&sample_output(), &tally, "ww1", "m.wwd"),
            Err(Error::BadMeshBounds { axis: 0, index: 1 })
        ));
        tally.x_bounds = vec![0.0];
        assert!(matches!(
            emit_serpent_wwin(&sample_output(), &tally, "ww1", "m.wwd"),
            Err(Error::BadMeshBounds { axis: 0, index: 0 })
        ));
    }

    #[test]
    fn length_mismatch_is_rejected() {
        let mut output = sample_output();
        output.lower_bounds_ww.pop();
        assert!(matches!(
            emit_serpent_wwin(&output, &sample_tally(), "ww1", "m.wwd"),
            Err(Error::LengthMismatch {
                expected: 8,
                got: 7
            })
        ));
    }

    #[test]
    fn null_cells_produce_drift_notes() {
        let mut output = sample_output();
        output.lower_bounds_ww[1] = 0.0;
        output.lower_bounds_ww[5] = 0.0;
        let serpent = emit_serpent_wwin(&output, &sample_tally(), "ww1", "m.wwd").unwrap();
        let openmc =
            emit_openmc_weight_windows(&output, &sample_tally(), &OpenMcOptions::default())
                .unwrap();
        let note = |notes: &[String]| {
            notes
                .iter()
                .find(|n| n.starts_with("null-cells"))
                .unwrap()
                .clone()
        };
        assert!(note(&serpent.notes).contains("2 of 8"));
        assert!(note(&openmc.notes).contains("2 of 8"));
        // And the zeros must actually appear in the outputs.
        assert!(serpent.text.contains("0.00000E+00"));
        assert!(openmc.xml.contains(" 0.0 "));
    }
}
