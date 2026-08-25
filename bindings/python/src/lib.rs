//! Python bindings (`nucleide._internal`).
//!
//! Thin facade only: all logic lives in workspace crates so the Rust API
//! stays usable without Python. Type stubs live in `python/nucleide/_internal.pyi`.

use std::collections::BTreeMap;

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;

use nuclei::NuclideId;

/// Package version, re-exported to Python.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

fn wrap_nucid_err(e: nuclei::Error) -> PyErr {
    PyValueError::new_err(e.to_string())
}

// ---------------------------------------------------------------------------
// Nuclide naming
// ---------------------------------------------------------------------------

/// A nuclide identifier (canonical nucid integer + naming conversions).
#[pyclass(name = "Nuclide")]
struct PyNuclide {
    inner: NuclideId,
}

#[pymethods]
impl PyNuclide {
    /// Create from a name such as "U235" or "Am242_m1".
    #[new]
    fn new(name: &str) -> PyResult<Self> {
        NuclideId::from_name(name)
            .map(|inner| Self { inner })
            .map_err(wrap_nucid_err)
    }

    /// GNDS-style name ("U235", "Am242_m1").
    #[getter]
    fn name(&self) -> String {
        self.inner.to_name()
    }

    /// Raw nucid integer.
    #[getter]
    fn nucid(&self) -> u32 {
        self.inner.nucid()
    }

    /// ZZAAAM form (922350 for U-235).
    #[getter]
    fn zzaaam(&self) -> u32 {
        self.inner.zzaaam()
    }

    /// Atomic number.
    #[getter]
    fn z(&self) -> u32 {
        self.inner.z()
    }

    /// Mass number.
    #[getter]
    fn a(&self) -> u32 {
        self.inner.a()
    }

    /// Metastable state index (0 = ground).
    #[getter]
    fn state(&self) -> u32 {
        self.inner.state()
    }

    /// MCNP ZAID integer.
    #[getter]
    fn zaid(&self) -> u32 {
        nuclei::dialects::to_zaid(self.inner)
    }

    /// zzllaaam form ("U-235").
    #[getter]
    fn zzllaaam(&self) -> String {
        nuclei::dialects::zzllaaam(self.inner)
    }

    /// Serpent-style name ("U-235").
    #[getter]
    fn serpent(&self) -> String {
        nuclei::dialects::serpent(self.inner)
    }

    /// NIST-style name.
    #[getter]
    fn nist(&self) -> String {
        nuclei::dialects::nist(self.inner)
    }

    /// Cinder integer id.
    #[getter]
    fn cinder(&self) -> u32 {
        nuclei::dialects::to_cinder(self.inner)
    }

    /// ALARA name ("u:235").
    #[getter]
    fn alara(&self) -> String {
        nuclei::dialects::alara(self.inner)
    }

    /// SZA integer.
    #[getter]
    fn sza(&self) -> u32 {
        nuclei::dialects::to_sza(self.inner)
    }

    /// FLUKA element-isotope name; raises ValueError if unavailable.
    fn fluka(&self) -> PyResult<&'static str> {
        nuclei::dialects::id_to_fluka(self.inner).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Atomic mass in u (AME2020), or None if unknown.
    #[getter]
    fn mass(&self) -> Option<f64> {
        nuclei::data::atomic_mass(self.inner.nucid())
    }

    /// Natural abundance fraction, or None.
    #[getter]
    fn abundance(&self) -> Option<f64> {
        nuclei::data::natural_abundance(self.inner.nucid())
    }

    fn __repr__(&self) -> String {
        format!("Nuclide({})", self.inner.to_name())
    }
}

/// Parse a MCNP ZAID integer into a Nuclide.
#[pyfunction]
fn from_zaid(zaid: u32) -> PyResult<PyNuclide> {
    nuclei::dialects::from_zaid(zaid)
        .map(|inner| PyNuclide { inner })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

fn lookup(key: &Bound<'_, PyAny>, f: impl Fn(u32) -> Option<f64>) -> PyResult<Option<f64>> {
    if let Ok(nucid) = key.extract::<u32>() {
        return Ok(f(nucid));
    }
    if let Ok(name) = key.extract::<&str>() {
        let id = NuclideId::from_name(name).map_err(wrap_nucid_err)?;
        return Ok(f(id.nucid()));
    }
    Err(PyTypeError::new_err("expected int nucid or str name"))
}

/// Atomic mass in u for a nucid integer or name string.
#[pyfunction]
fn atomic_mass(key: &Bound<'_, PyAny>) -> PyResult<Option<f64>> {
    lookup(key, nuclei::data::atomic_mass)
}

/// Natural abundance fraction for a nucid integer or name string.
#[pyfunction]
fn natural_abundance(key: &Bound<'_, PyAny>) -> PyResult<Option<f64>> {
    lookup(key, nuclei::data::natural_abundance)
}

/// A particle species with cross-code name translations.
#[pyclass(name = "Particle")]
struct PyParticle {
    inner: nuclei::particles::ParticleId,
}

#[pymethods]
impl PyParticle {
    /// Create from any alias ("n", "neutron", "gamma", PDC int, ...).
    #[new]
    fn new(spec: &Bound<'_, PyAny>) -> PyResult<Self> {
        let inner = if let Ok(pdc) = spec.extract::<i32>() {
            nuclei::particles::ParticleId::from_pdc(pdc)
                .ok_or_else(|| PyValueError::new_err(format!("unknown PDC code {pdc}")))?
        } else if let Ok(s) = spec.extract::<&str>() {
            s.parse::<nuclei::particles::ParticleId>()
                .map_err(|e| PyValueError::new_err(e.to_string()))?
        } else {
            return Err(PyTypeError::new_err("expected str alias or int PDC"));
        };
        Ok(Self { inner })
    }

    #[getter]
    fn name(&self) -> &'static str {
        self.inner.name()
    }

    #[getter]
    fn describe(&self) -> &'static str {
        self.inner.describe()
    }

    fn mcnp(&self) -> Option<&'static str> {
        self.inner.mcnp()
    }
    fn mcnp6(&self) -> Option<&'static str> {
        self.inner.mcnp6()
    }
    fn fluka(&self) -> Option<&'static str> {
        self.inner.fluka()
    }
    fn geant4(&self) -> Option<&'static str> {
        self.inner.geant4()
    }

    fn __repr__(&self) -> String {
        format!("Particle('{}')", self.inner.name())
    }
}

/// Resolve a reaction name/MT/id string to its numeric id.
#[pyfunction]
fn rxname_id(name: &str) -> PyResult<u32> {
    nuclei::rxname::name_to_id(name).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Canonical short name for a reaction id.
#[pyfunction]
fn rxname_name(id: u32) -> Option<&'static str> {
    nuclei::rxname::id_to_name(id)
}

/// ENDF MT number for a reaction id (0 if none registered).
#[pyfunction]
fn rxname_mt(id: u32) -> i32 {
    nuclei::rxname::id_to_mt(id)
}

// ---------------------------------------------------------------------------
// MCNP file I/O
// ---------------------------------------------------------------------------

fn io_err(e: mcnp_io::xsdir::Error) -> PyErr {
    PyValueError::new_err(e.to_string())
}
fn m_err<T>(r: Result<T, impl std::fmt::Display>) -> PyResult<T> {
    r.map_err(|e| PyValueError::new_err(e.to_string()))
}

/// One xsdir directory entry.
#[pyclass(name = "XsdirTable")]
struct PyXsdirTable {
    inner: mcnp_io::xsdir::XsdirTable,
}

#[pymethods]
impl PyXsdirTable {
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    fn awr(&self) -> f64 {
        self.inner.awr
    }
    #[getter]
    fn filename(&self) -> &str {
        &self.inner.filename
    }
    #[getter]
    fn filetype(&self) -> i64 {
        self.inner.filetype
    }
    #[getter]
    fn address(&self) -> i64 {
        self.inner.address
    }
    #[getter]
    fn tablelength(&self) -> i64 {
        self.inner.tablelength
    }
    #[getter]
    fn temperature(&self) -> Option<f64> {
        self.inner.temperature
    }
    #[getter]
    fn ptable(&self) -> bool {
        self.inner.ptable
    }
    /// ZAID text before the '.'.
    fn zaid(&self) -> &str {
        self.inner.zaid()
    }
    /// Serpent directory-entry line.
    fn to_serpent(&self, directory: &str) -> PyResult<String> {
        m_err(self.inner.to_serpent(directory))
    }
    fn __repr__(&self) -> String {
        format!("<XsdirTable: {}>", self.inner.name)
    }
}

/// Parsed xsdir index file.
#[pyclass(name = "Xsdir")]
struct PyXsdir {
    inner: mcnp_io::xsdir::Xsdir,
}

#[pymethods]
impl PyXsdir {
    #[getter]
    fn datapath(&self) -> Option<&str> {
        self.inner.datapath.as_deref()
    }
    /// Atomic weight ratios keyed by zaid integer.
    #[getter]
    fn awr(&self) -> BTreeMap<u32, f64> {
        self.inner.awr.clone()
    }
    /// Directory entries in file order.
    #[getter]
    fn tables(&self) -> Vec<PyXsdirTable> {
        self.inner
            .tables
            .iter()
            .map(|t| PyXsdirTable { inner: t.clone() })
            .collect()
    }
    /// Tables whose name contains `name`.
    fn find_table(&self, name: &str) -> Vec<PyXsdirTable> {
        self.inner
            .find_table(name)
            .into_iter()
            .map(|t| PyXsdirTable { inner: t.clone() })
            .collect()
    }
    /// Distinct nuclides referenced by the entries.
    fn nucs(&self) -> Vec<u32> {
        self.inner.nucs().iter().map(|n| n.nucid()).collect()
    }
}

/// Parse an MCNP xsdir file.
#[pyfunction]
fn read_xsdir(path: &str) -> PyResult<PyXsdir> {
    mcnp_io::xsdir::Xsdir::from_file(path)
        .map(|inner| PyXsdir { inner })
        .map_err(io_err)
}

/// One fmesh4 tally from a meshtal file.
#[pyclass(name = "MeshTally")]
struct PyMeshTally {
    inner: mcnp_io::meshtal::MeshTallyData,
}

#[pymethods]
impl PyMeshTally {
    #[getter]
    fn tally_number(&self) -> u32 {
        self.inner.tally_number
    }
    /// 'n', 'p', ...
    #[getter]
    fn particle(&self) -> char {
        self.inner.particle.letter()
    }
    #[getter]
    fn dose_response(&self) -> bool {
        self.inner.dose_response
    }
    #[getter]
    fn x_bounds(&self) -> Vec<f64> {
        self.inner.x_bounds.clone()
    }
    #[getter]
    fn y_bounds(&self) -> Vec<f64> {
        self.inner.y_bounds.clone()
    }
    #[getter]
    fn z_bounds(&self) -> Vec<f64> {
        self.inner.z_bounds.clone()
    }
    #[getter]
    fn e_bounds(&self) -> Vec<f64> {
        self.inner.e_bounds.clone()
    }
    /// [nx, ny, nz] cell counts.
    fn dims(&self) -> [usize; 3] {
        self.inner.dims()
    }
    fn num_ves(&self) -> usize {
        self.inner.num_ves()
    }
    fn num_e_groups(&self) -> usize {
        self.inner.num_e_groups()
    }
    /// All-group results for cell (i,j,k): [result_per_group, error_per_group].
    fn cell(&self, i: usize, j: usize, k: usize) -> (Vec<f64>, Vec<f64>) {
        let (r, e) = self.inner.cell(i, j, k);
        (r.to_vec(), e.to_vec())
    }
    /// Energy-integrated totals for cell (i,j,k).
    fn cell_total(&self, i: usize, j: usize, k: usize) -> (f64, f64) {
        self.inner.cell_total(i, j, k)
    }
    /// Full results array `[ve][group]`.
    #[getter]
    fn result(&self) -> Vec<Vec<f64>> {
        self.inner.result.clone()
    }
    /// Full relative-error array `[ve][group]`.
    #[getter]
    fn rel_error(&self) -> Vec<Vec<f64>> {
        self.inner.rel_error.clone()
    }
    /// Per-cell energy-integrated totals.
    #[getter]
    fn total_result(&self) -> Vec<f64> {
        self.inner.total_result.clone()
    }
}

/// Parsed meshtal file.
#[pyclass(name = "Meshtal")]
struct PyMeshtal {
    inner: mcnp_io::meshtal::Meshtal,
}

#[pymethods]
impl PyMeshtal {
    #[getter]
    fn version(&self) -> &str {
        &self.inner.version
    }
    #[getter]
    fn ld(&self) -> &str {
        &self.inner.ld
    }
    #[getter]
    fn title(&self) -> &str {
        &self.inner.title
    }
    #[getter]
    fn histories(&self) -> u64 {
        self.inner.histories
    }
    /// Tallies keyed by fmesh4 number.
    #[getter]
    fn tallies(&self) -> BTreeMap<u32, PyMeshTally> {
        self.inner
            .tallies
            .iter()
            .map(|(k, v)| (*k, PyMeshTally { inner: v.clone() }))
            .collect()
    }
}

/// Parse an MCNP meshtal file.
#[pyfunction]
fn read_meshtal(path: &str) -> PyResult<PyMeshtal> {
    m_err(mcnp_io::meshtal::Meshtal::from_file(path).map(|inner| PyMeshtal { inner }))
}

/// Parsed WWINP weight-window file.
#[pyclass(name = "Wwinp")]
struct PyWwinp {
    inner: mcnp_io::wwinp::Wwinp,
}

#[pymethods]
impl PyWwinp {
    #[getter]
    fn ni(&self) -> u32 {
        self.inner.ni
    }
    #[getter]
    fn nr(&self) -> u32 {
        self.inner.nr
    }
    #[getter]
    fn ne(&self) -> Vec<u32> {
        self.inner.ne.clone()
    }
    #[getter]
    fn nf(&self) -> [u32; 3] {
        self.inner.nf
    }
    #[getter]
    fn origin(&self) -> [f64; 3] {
        self.inner.origin
    }
    #[getter]
    fn nc(&self) -> [u32; 3] {
        self.inner.nc
    }
    /// Coarse boundaries per dimension.
    #[getter]
    fn cm(&self) -> Vec<Vec<f64>> {
        self.inner.cm.clone()
    }
    /// Expanded spatial bounds per dimension.
    #[getter]
    fn bounds(&self) -> Vec<Vec<f64>> {
        self.inner.bounds.clone()
    }
    /// Energy upper bounds per particle present.
    #[getter]
    fn e(&self) -> Vec<Vec<f64>> {
        self.inner.e.clone()
    }
    /// Lower bounds for one group: ww_row(particle, group) -> list[nve].
    fn ww_row(&self, particle: usize, group: usize) -> Vec<f64> {
        self.inner.ww[particle][group].clone()
    }
    /// Lower-bound vector for one volume element across groups.
    fn ww_column(&self, particle: usize, ve: usize) -> Vec<f64> {
        self.inner.ww_column(particle, ve)
    }
}

/// Parse an MCNP WWINP weight-window file.
#[pyfunction]
fn read_wwinp(path: &str) -> PyResult<PyWwinp> {
    m_err(mcnp_io::wwinp::Wwinp::from_file(path).map(|inner| PyWwinp { inner }))
}

/// Parsed MCTAL kcode data.
#[pyclass(name = "Mctal")]
struct PyMctal {
    inner: mcnp_io::mctal::Mctal,
}

#[pymethods]
impl PyMctal {
    #[getter]
    fn code_name(&self) -> &str {
        &self.inner.code_name
    }
    #[getter]
    fn comment(&self) -> &str {
        &self.inner.comment
    }
    #[getter]
    fn n_histories(&self) -> u64 {
        self.inner.n_histories
    }
    #[getter]
    fn n_cycles(&self) -> usize {
        self.inner.n_cycles
    }
    #[getter]
    fn n_inactive(&self) -> usize {
        self.inner.n_inactive
    }
    #[getter]
    fn vars_per_cycle(&self) -> usize {
        self.inner.vars_per_cycle
    }
    #[getter]
    fn k_col(&self) -> Vec<f64> {
        self.inner.k_col.clone()
    }
    #[getter]
    fn k_abs(&self) -> Vec<f64> {
        self.inner.k_abs.clone()
    }
    #[getter]
    fn k_path(&self) -> Vec<f64> {
        self.inner.k_path.clone()
    }
    #[getter]
    fn prompt_life_col(&self) -> Vec<f64> {
        self.inner.prompt_life_col.clone()
    }
    #[getter]
    fn prompt_life_path(&self) -> Vec<f64> {
        self.inner.prompt_life_path.clone()
    }
    /// Running averages (empty unless vars_per_cycle == 19); each entry is a
    /// dict of the averaged pairs plus cycle_histories/fom.
    #[getter]
    fn averages(&self) -> Vec<BTreeMap<String, f64>> {
        self.inner
            .averages
            .iter()
            .map(|a| {
                let mut m = BTreeMap::new();
                m.insert("avg_k_col".into(), a.avg_k_col.0);
                m.insert("avg_k_col_stdev".into(), a.avg_k_col.1);
                m.insert("avg_k_abs".into(), a.avg_k_abs.0);
                m.insert("avg_k_abs_stdev".into(), a.avg_k_abs.1);
                m.insert("avg_k_path".into(), a.avg_k_path.0);
                m.insert("avg_k_path_stdev".into(), a.avg_k_path.1);
                m.insert("avg_k_combined".into(), a.avg_k_combined.0);
                m.insert("avg_k_combined_stdev".into(), a.avg_k_combined.1);
                m.insert("avg_k_combined_active".into(), a.avg_k_combined_active.0);
                m.insert(
                    "avg_k_combined_active_stdev".into(),
                    a.avg_k_combined_active.1,
                );
                m.insert("prompt_life_combined".into(), a.prompt_life_combined.0);
                m.insert(
                    "prompt_life_combined_stdev".into(),
                    a.prompt_life_combined.1,
                );
                m.insert("cycle_histories".into(), a.cycle_histories);
                m.insert("fom".into(), a.fom);
                m
            })
            .collect()
    }
}

/// Parse an MCNP MCTAL file (kcode subset, upstream parity).
#[pyfunction]
fn read_mctal(path: &str) -> PyResult<PyMctal> {
    m_err(mcnp_io::mctal::Mctal::from_file(path).map(|inner| PyMctal { inner }))
}

/// Parsed SSW surface-source file.
#[pyclass(name = "SurfSrc")]
struct PySurfSrc {
    inner: mcnp_io::surfsrc::SurfSrc,
}

#[pymethods]
impl PySurfSrc {
    #[getter]
    fn kod(&self) -> String {
        self.inner.header.kod.trim_end().to_string()
    }
    #[getter]
    fn ver(&self) -> String {
        self.inner.header.ver.trim_end().to_string()
    }
    #[getter]
    fn np1(&self) -> i64 {
        self.inner.header.np1
    }
    #[getter]
    fn nrss(&self) -> i64 {
        self.inner.header.nrss
    }
    #[getter]
    fn ncrd(&self) -> i32 {
        self.inner.header.ncrd
    }
    #[getter]
    fn njsw(&self) -> i32 {
        self.inner.header.njsw
    }
    #[getter]
    fn niss(&self) -> i64 {
        self.inner.header.niss
    }
    /// Formatted header block.
    fn print_header(&self) -> String {
        self.inner.header.print_header()
    }
    /// Track records as dicts of named fields.
    fn tracks(&self) -> PyResult<Vec<BTreeMap<String, f64>>> {
        let tracks = self
            .inner
            .read_tracklist()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(tracks
            .iter()
            .map(|t| {
                let mut d = BTreeMap::new();
                d.insert("nps".into(), t.nps);
                d.insert("bitarray".into(), t.bitarray);
                d.insert("wgt".into(), t.wgt);
                d.insert("erg".into(), t.erg);
                d.insert("tme".into(), t.tme);
                d.insert("x".into(), t.x);
                d.insert("y".into(), t.y);
                d.insert("z".into(), t.z);
                d.insert("u".into(), t.u);
                d.insert("v".into(), t.v);
                d.insert("cs".into(), t.cs);
                d.insert("w".into(), t.w);
                d
            })
            .collect())
    }
}

/// Read an MCNP SSW surface-source file (header eagerly; tracks on demand).
#[pyfunction]
fn read_ssw(path: &str) -> PyResult<PySurfSrc> {
    mcnp_io::surfsrc::SurfSrc::open(path)
        .map(|inner| PySurfSrc { inner })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Detected PTRAC layout: 0 = i4 little-endian, 1 = i8 little-endian.
#[pyclass(name = "PtracFile")]
struct PyPtracFile {
    inner: mcnp_io::ptrac::PtracFile,
}

#[pymethods]
impl PyPtracFile {
    #[getter]
    fn problem_title(&self) -> &str {
        &self.inner.problem_title
    }
    /// 0 for i4, 1 for i8.
    #[getter]
    fn width_code(&self) -> u8 {
        match self.inner.format {
            mcnp_io::ptrac::Format::I4LittleEndian => 0,
            mcnp_io::ptrac::Format::I8LittleEndian => 1,
        }
    }
    /// Variable counts per event type as {nps,src,bnk,sur,col,ter}.
    #[getter]
    fn variable_nums(&self) -> BTreeMap<String, usize> {
        let v = &self.inner.variable_nums;
        let mut m = BTreeMap::new();
        m.insert("nps".into(), v.nps);
        m.insert("src".into(), v.src);
        m.insert("bnk".into(), v.bnk);
        m.insert("sur".into(), v.sur);
        m.insert("col".into(), v.col);
        m.insert("ter".into(), v.ter);
        m
    }
    /// All events as dicts: {'event_type': int, '<var>': float, ...}.
    fn events(&self) -> PyResult<Vec<BTreeMap<String, f64>>> {
        let events = self
            .inner
            .events()
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(events
            .iter()
            .map(|ev| {
                let mut d = BTreeMap::new();
                d.insert("event_type".to_string(), ev.event_type as f64);
                for (n, v) in ev.iter() {
                    d.insert(n.to_string(), v);
                }
                d
            })
            .collect())
    }
}

/// Read an MCNP PTRAC event file.
#[pyfunction]
fn read_ptrac(path: &str) -> PyResult<PyPtracFile> {
    mcnp_io::ptrac::PtracFile::open(path)
        .map(|inner| PyPtracFile { inner })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

// ---------------------------------------------------------------------------
// Depletion / CRAM
// ---------------------------------------------------------------------------

/// A parsed depletion chain (XML format).
#[pyclass(name = "Chain")]
struct PyChain {
    inner: std::sync::Arc<depletion::Chain>,
}

#[pymethods]
impl PyChain {
    /// Nuclide names in chain order.
    #[getter]
    fn nuclides(&self) -> Vec<String> {
        self.inner.nuclides.iter().map(|n| n.name.clone()).collect()
    }

    fn index_of(&self, name: &str) -> Option<usize> {
        self.inner.index_of(name)
    }
}

/// Parse a depletion-chain XML file.
#[pyfunction]
fn read_chain(path: &str) -> PyResult<PyChain> {
    depletion::Chain::from_file(path)
        .map(|inner| PyChain {
            inner: std::sync::Arc::new(inner),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// One-group reaction rates keyed by "NuclideName:reaction".
type RateMap = BTreeMap<String, f64>;

fn split_rates(rates: &RateMap, chain: &depletion::Chain) -> PyResult<depletion::ReactionRates> {
    let mut out = depletion::ReactionRates::new();
    for (key, v) in rates {
        let (nuc, rx) = key.split_once(':').ok_or_else(|| {
            PyValueError::new_err(format!("rate key `{key}` must be `Name:reaction`"))
        })?;
        let idx = chain
            .index_of(nuc)
            .ok_or_else(|| PyValueError::new_err(format!("rate for unknown nuclide `{nuc}`")))?;
        out.insert((idx, rx.to_string()), *v);
    }
    Ok(out)
}

/// Solve one depletion step with IPF CRAM.
///
/// `n0` maps nuclide names to initial atom counts; `rates` maps
/// `"Name:(n,gamma)"`-style keys to one-group rates [1/s]; `dt` is the step
/// length in seconds; `order` is 16 or 48.
#[pyfunction]
#[pyo3(signature = (chain, n0, dt, rates=None, order=48))]
fn deplete(
    chain: &PyChain,
    n0: BTreeMap<String, f64>,
    dt: f64,
    rates: Option<RateMap>,
    order: u8,
) -> PyResult<BTreeMap<String, f64>> {
    let order = match order {
        16 => depletion::Order::Order16,
        48 => depletion::Order::Order48,
        other => {
            return Err(PyValueError::new_err(format!(
                "unsupported CRAM order {other}"
            )))
        }
    };
    let rates = split_rates(rates.as_ref().unwrap_or(&BTreeMap::new()), &chain.inner)?;
    let sys = depletion::DepletionSystem::build((*chain.inner).clone(), &rates)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    depletion::deplete(&sys, order, &n0, dt)
        .map(|r| r.atoms)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

// ---------------------------------------------------------------------------
// Serpent / FLUKA / variance reduction + writers
// ---------------------------------------------------------------------------

/// Parse a Serpent .m output file ("res", "dep", or "det") into a nested
/// Python dict.
#[pyfunction]
fn read_serpent(path: &str, kind: &str) -> PyResult<Py<PyAny>> {
    let text = std::fs::read_to_string(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
    let table = match kind {
        "res" => serpent_io::parse_res(&text),
        "dep" => serpent_io::parse_dep(&text),
        "det" => serpent_io::parse_det(&text),
        other => {
            return Err(PyValueError::new_err(format!(
                "kind must be res|dep|det, got `{other}`"
            )))
        }
    }
    .map_err(|e| PyValueError::new_err(e.to_string()))?;
    fn entry_to_py(py: Python<'_>, e: &serpent_io::Entry) -> Py<PyAny> {
        use serpent_io::Entry as E;
        match e {
            E::Scalar(serpent_io::Value::Num(n)) => {
                n.into_pyobject(py).unwrap().unbind().into_any()
            }
            E::Scalar(serpent_io::Value::Str(s)) => {
                s.into_pyobject(py).unwrap().unbind().into_any()
            }
            E::Vector(vs) => vs
                .iter()
                .map(|v| match v {
                    serpent_io::Value::Num(n) => n.into_pyobject(py).unwrap().unbind().into_any(),
                    serpent_io::Value::Str(s) => s.into_pyobject(py).unwrap().unbind().into_any(),
                })
                .collect::<Vec<_>>()
                .into_pyobject(py)
                .unwrap()
                .unbind()
                .into_any(),
            E::Matrix(m) => {
                let rows: Vec<Py<PyAny>> = m
                    .to_rows_f64()
                    .iter()
                    .map(|row| row.into_pyobject(py).unwrap().unbind().into_any())
                    .collect();
                rows.into_pyobject(py).unwrap().unbind().into_any()
            }
        }
    }
    Ok(Python::attach(|py| {
        let dict = pyo3::types::PyDict::new(py);
        for (k, e) in table.iter() {
            dict.set_item(k, entry_to_py(py, e)).ok();
        }
        dict.into_any().unbind()
    }))
}

/// One FLUKA USRBIN detector.
#[pyclass(name = "UsrbinTally")]
struct PyUsrbinTally {
    inner: fluka_io::usrbin::UsrbinTally,
}

#[pymethods]
impl PyUsrbinTally {
    #[getter]
    fn name(&self) -> &str {
        &self.inner.name
    }
    #[getter]
    fn particle(&self) -> &str {
        &self.inner.particle
    }
    #[getter]
    fn nx(&self) -> usize {
        self.inner.x_info.bins
    }
    #[getter]
    fn ny(&self) -> usize {
        self.inner.y_info.bins
    }
    #[getter]
    fn nz(&self) -> usize {
        self.inner.z_info.bins
    }
    #[getter]
    fn x_bounds(&self) -> Vec<f64> {
        self.inner.x_bounds.clone()
    }
    #[getter]
    fn y_bounds(&self) -> Vec<f64> {
        self.inner.y_bounds.clone()
    }
    #[getter]
    fn z_bounds(&self) -> Vec<f64> {
        self.inner.z_bounds.clone()
    }
    /// Scored values, x slowest -> z fastest.
    #[getter]
    fn data(&self) -> Vec<f64> {
        self.inner.part_data.clone()
    }
    /// Statistical errors, same layout as `data`.
    #[getter]
    fn error(&self) -> Vec<f64> {
        self.inner.error_data.clone()
    }
    fn dims(&self) -> [usize; 3] {
        [self.nx(), self.ny(), self.nz()]
    }
}

/// Parse all USRBIN tallies from a FLUKA .lis file.
#[pyfunction]
fn read_usrbin(path: &str) -> PyResult<Vec<PyUsrbinTally>> {
    let tallies = fluka_io::usrbin::read_usrbin_file(path)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(tallies
        .into_iter()
        .map(|inner| PyUsrbinTally { inner })
        .collect())
}

/// MAGIC weight-window output.
#[pyclass(name = "MagicOutput")]
struct PyMagicOutput {
    inner: vr_tools::magic::MagicOutput,
}

#[pymethods]
impl PyMagicOutput {
    /// Flat lower bounds ([ve] in total mode, [ve*g+g] per-group).
    #[getter]
    fn lower_bounds_ww(&self) -> Vec<f64> {
        self.inner.lower_bounds_ww.clone()
    }
    #[getter]
    fn groups_per_ve(&self) -> usize {
        self.inner.groups_per_ve
    }
    #[getter]
    fn scale_factors(&self) -> Vec<f64> {
        self.inner.scale_factors.clone()
    }
    #[getter]
    fn e_upper_bounds(&self) -> Vec<f64> {
        self.inner.e_upper_bounds.clone()
    }
    #[getter]
    fn ww_tag_name(&self) -> &str {
        &self.inner.ww_tag_name
    }
}

/// Generate MAGIC weight-window lower bounds from a meshtal tally.
#[pyfunction]
#[pyo3(signature = (tally, per_group=false, tolerance=0.5))]
fn magic(tally: &PyMeshTally, per_group: bool, tolerance: f64) -> PyResult<PyMagicOutput> {
    let selection = if per_group {
        vr_tools::magic::MagicSelection::PerGroup
    } else {
        vr_tools::magic::MagicSelection::Total
    };
    let params = vr_tools::magic::MagicParams {
        tolerance,
        ..Default::default()
    };
    vr_tools::magic::magic_with(&tally.inner, selection, params)
        .map(|inner| PyMagicOutput { inner })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Walker alias table for discrete sampling.
#[pyclass(name = "AliasTable")]
struct PyAliasTable {
    inner: vr_tools::sampling::AliasTable,
}

#[pymethods]
impl PyAliasTable {
    /// Build from a probability density (normalized internally).
    #[new]
    fn new(pdf: Vec<f64>) -> PyResult<Self> {
        vr_tools::sampling::AliasTable::new(&pdf)
            .map(|inner| PyAliasTable { inner })
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }
    /// Sample an index from two uniform random numbers.
    fn sample(&self, r1: f64, r2: f64) -> usize {
        self.inner.sample(r1, r2)
    }
    #[getter]
    fn pdf(&self) -> Vec<f64> {
        self.inner.pdf().to_vec()
    }
    fn __len__(&self) -> usize {
        self.inner.len()
    }
}

/// Mesh source sampler over a meshtal tally (ANALOG/UNIFORM/USER modes).
#[pyclass(name = "MeshSourceSampler")]
struct PyMeshSourceSampler {
    inner: vr_tools::sampling::MeshSourceSampler,
}

#[pymethods]
impl PyMeshSourceSampler {
    /// mode: "analog" | "uniform" | "user" (user requires user_pdf).
    #[new]
    #[pyo3(signature = (tally, mode, user_pdf=None))]
    fn new(tally: &PyMeshTally, mode: &str, user_pdf: Option<Vec<f64>>) -> PyResult<Self> {
        let user = if matches!(mode, "user") {
            Some(user_pdf.ok_or_else(|| PyValueError::new_err("user mode needs user_pdf"))?)
        } else {
            None
        };
        let m = match mode {
            "analog" => vr_tools::sampling::Mode::Analog,
            "uniform" => vr_tools::sampling::Mode::Uniform,
            "user" => vr_tools::sampling::Mode::User,
            other => {
                return Err(PyValueError::new_err(format!(
                    "mode must be analog|uniform|user, got `{other}`"
                )))
            }
        };
        vr_tools::sampling::MeshSourceSampler::new(&tally.inner, m, user.as_deref())
            .map(|inner| PyMeshSourceSampler { inner })
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }
    /// Sample a voxel; returns dict(index, i, j, k, weight).
    fn sample(&self, r1: f64, r2: f64) -> BTreeMap<String, f64> {
        let s = self.inner.sample(r1, r2);
        let mut d = BTreeMap::new();
        d.insert("index".into(), s.index as f64);
        d.insert("i".into(), s.i as f64);
        d.insert("j".into(), s.j as f64);
        d.insert("k".into(), s.k as f64);
        d.insert("weight".into(), s.weight);
        d
    }
}

/// Write a SurfSrc file back to disk. `tracks` defaults to re-reading the
/// original file's tracks.
#[pyfunction]
#[pyo3(signature = (ssw, path, tracks=None))]
fn write_ssw(
    ssw: &PySurfSrc,
    path: &str,
    tracks: Option<Vec<BTreeMap<String, f64>>>,
) -> PyResult<()> {
    let header = ssw.inner.header.clone();
    let track_data: Vec<mcnp_io::surfsrc::TrackData> = match tracks {
        Some(dict_tracks) => dict_tracks
            .iter()
            .map(|d| {
                let g = |k: &str| d.get(k).copied().unwrap_or(0.0);
                let mut record = vec![0.0f64; mcnp_io::surfsrc::TrackData::RECORD_WIDTH];
                record[0] = g("nps");
                record[1] = g("bitarray");
                record[2] = g("wgt");
                record[3] = g("erg");
                record[4] = g("tme");
                record[5] = g("x");
                record[6] = g("y");
                record[7] = g("z");
                record[8] = g("u");
                record[9] = g("v");
                record[10] = g("cs");
                mcnp_io::surfsrc::TrackData::from_record(record)
            })
            .collect(),
        None => ssw
            .inner
            .read_tracklist()
            .map_err(|e| PyValueError::new_err(e.to_string()))?,
    };
    let mut f = std::fs::File::create(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
    mcnp_io::surfsrc::write_to(&mut f, &header, &track_data)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Generate MCNP input-deck text from a structured mesh.
#[pyfunction]
fn mesh_to_geom(
    x_bounds: Vec<f64>,
    y_bounds: Vec<f64>,
    z_bounds: Vec<f64>,
    cell_materials: Vec<Option<(String, f64)>>,
    title_card: &str,
) -> String {
    let opts = mcnp_io::deck::DeckOptions {
        title_card: title_card.to_string(),
        frac_type: mcnp_io::deck::FracType::Mass,
    };
    mcnp_io::deck::mesh_to_geom(&x_bounds, &y_bounds, &z_bounds, &cell_materials, &opts)
}

// ---------------------------------------------------------------------------
// Data accessors, input parsing, enrichment, materials
// ---------------------------------------------------------------------------

/// Half-life [s] for a nucid integer or name string.
#[pyfunction]
fn half_life(key: &Bound<'_, PyAny>) -> PyResult<Option<f64>> {
    lookup(key, nuclei::data::half_life)
}

/// Decay constant lambda = ln2 / t_half [1/s].
#[pyfunction]
fn decay_constant(key: &Bound<'_, PyAny>) -> PyResult<Option<f64>> {
    lookup(key, nuclei::data::decay_constant)
}

/// Neutron-capture Q value computed from AME2020 masses [MeV].
#[pyfunction]
fn q_value_capture(key: &Bound<'_, PyAny>) -> PyResult<Option<f64>> {
    lookup(key, nuclei::data::q_value_neutron_capture)
}

/// Alpha-decay Q value from AME2020 masses [MeV].
#[pyfunction]
fn q_value_alpha(key: &Bound<'_, PyAny>) -> PyResult<Option<f64>> {
    lookup(key, nuclei::data::q_value_alpha)
}

/// Parse MCNP material cards from an input deck.
/// Returns a list of dicts: {number, fractions: {NuclideName: frac},
/// fraction_type: "atom"|"mass", density, comments}.
#[pyfunction]
fn read_inp(path: &str) -> PyResult<Vec<BTreeMap<String, Py<PyAny>>>> {
    let mats = mcnp_io::inp::materials_from_file(path)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Python::attach(|py| {
        Ok(mats
            .into_iter()
            .map(|m| {
                let mut d = BTreeMap::new();
                d.insert(
                    "number".to_string(),
                    m.number.into_pyobject(py).unwrap().unbind().into_any(),
                );
                let fr: BTreeMap<String, f64> = m
                    .fractions
                    .iter()
                    .map(|(id, f)| (id.to_name(), *f))
                    .collect();
                d.insert(
                    "fractions".to_string(),
                    fr.into_pyobject(py).unwrap().unbind().into_any(),
                );
                d.insert(
                    "fraction_type".to_string(),
                    match m.fraction_type {
                        mcnp_io::inp::FracKind::Atom => "atom",
                        mcnp_io::inp::FracKind::Mass => "mass",
                    }
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
                );
                d.insert(
                    "density".to_string(),
                    m.density.into_pyobject(py).unwrap().unbind().into_any(),
                );
                d.insert(
                    "comments".to_string(),
                    m.comments
                        .join(" ")
                        .into_pyobject(py)
                        .unwrap()
                        .unbind()
                        .into_any(),
                );
                d
            })
            .collect())
    })
}

fn comp_to_material(comp: BTreeMap<String, f64>) -> PyResult<material::Material> {
    let mut mat = material::Material::new();
    for (name, grams) in &comp {
        let id = nuclei::NuclideId::from_name(name)
            .map_err(|e| PyValueError::new_err(format!("`{name}`: {e}")))?;
        mat.add_nuclide(id, *grams);
    }
    Ok(mat)
}

/// Expand a chemical formula into a natural-isotope composition dict
/// ({nuclide_name: atom_fraction}) using AME2020 masses + abundances.
#[pyfunction]
fn from_formula(formula: &str) -> PyResult<BTreeMap<String, f64>> {
    use material::AbundanceProvider;
    let parsed =
        material::parse_formula(formula).map_err(|e| PyValueError::new_err(e.to_string()))?;
    // Build a temporary element-count material then expand via abundances:
    let mut nat = Vec::new();
    for (z, count) in &parsed {
        if let Some(isotopes) = material::NaturalAbundances.natural_isotopes(*z) {
            for (id, frac) in isotopes {
                nat.push((id, frac * count));
            }
        }
    }
    let total: f64 = nat.iter().map(|(_, c)| c).sum();
    if total <= 0.0 {
        return Err(PyValueError::new_err("empty formula expansion"));
    }
    let mut out: BTreeMap<String, f64> = BTreeMap::new();
    for (id, atoms) in nat {
        *out.entry(id.to_name()).or_insert(0.0) += atoms / total;
    }
    Ok(out)
}

/// Activity [Bq] per nuclide plus whole-material specific activity.
/// Returns {name: Bq} entries and "specific" = Bq/g of the composition.
#[pyfunction]
fn activity(comp: BTreeMap<String, f64>) -> PyResult<BTreeMap<String, f64>> {
    let mat = comp_to_material(comp)?;
    let analytics = material::Analytics {
        masses: &material::Ame2020,
        decays: &material::ChainDecays,
    };
    let per_nuc = mat
        .activity(&analytics)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let specific = mat
        .specific_activity(&analytics)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let mut out: BTreeMap<String, f64> = per_nuc
        .into_iter()
        .map(|(id, v)| (id.to_name(), v))
        .collect();
    out.insert("specific".to_string(), specific);
    Ok(out)
}

/// Enrichment cascade with numeric multicomponent solving.
#[pyclass(name = "Cascade")]
struct PyCascade {
    inner: std::sync::Mutex<enrichment::Cascade>,
}

#[pymethods]
impl PyCascade {
    /// Natural-uranium default cascade (alpha=1.05, Mstar=236, j=U235, k=U238).
    #[staticmethod]
    fn default_uranium() -> Self {
        Self {
            inner: std::sync::Mutex::new(enrichment::default_uranium_cascade()),
        }
    }

    /// Solve via the numeric fixed-point + secant scheme in place.
    #[pyo3(signature = (tolerance=None, max_iterations=None))]
    fn solve(&self, tolerance: Option<f64>, max_iterations: Option<u32>) -> PyResult<()> {
        let tol = tolerance.unwrap_or(enrichment::DEFAULT_TOLERANCE);
        let iters = max_iterations.unwrap_or(enrichment::DEFAULT_MAX_ITER);
        let mut c = self
            .inner
            .lock()
            .map_err(|_| PyValueError::new_err("cascade lock poisoned"))?;
        *c = enrichment::solve_numeric(&c, tol, iters)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(())
    }

    #[getter]
    fn swu_per_feed(&self) -> f64 {
        self.inner.lock().unwrap().swu_per_feed
    }
    #[getter]
    fn swu_per_prod(&self) -> f64 {
        self.inner.lock().unwrap().swu_per_prod
    }
    /// Product assay of the key component j.
    #[getter]
    fn x_prod_j(&self) -> f64 {
        self.inner.lock().unwrap().x_prod_j
    }
    #[getter]
    fn x_tail_j(&self) -> f64 {
        self.inner.lock().unwrap().x_tail_j
    }
    /// Separative work per product [kg SWU/kg] from the key assays.
    fn separative_work_per_product(&self) -> f64 {
        let c = self.inner.lock().unwrap();
        enrichment::swu_per_prod(c.x_feed_j, c.x_prod_j, c.x_tail_j)
    }

    fn __repr__(&self) -> String {
        let c = self.inner.lock().unwrap();
        format!(
            "Cascade(alpha={}, Mstar={}, x_prod_j={:.5})",
            c.alpha, c.Mstar, c.x_prod_j
        )
    }
}

/// PNNL/DOE Materials Compendium library (411 named materials).
#[pyclass(name = "MaterialsCompendium")]
struct PyMaterialsCompendium {
    inner: material::MaterialsLibrary,
}

#[pymethods]
impl PyMaterialsCompendium {
    /// Load from the official MaterialsCompendium.json.
    #[staticmethod]
    fn load(path: &str) -> PyResult<Self> {
        material::MaterialsLibrary::from_file(path)
            .map(|inner| PyMaterialsCompendium { inner })
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// All display names in file order.
    fn names(&self) -> Vec<String> {
        self.inner.names().into_iter().map(String::from).collect()
    }

    /// Case-insensitive lookup by name; returns
    /// {name, mat_num, density, fractions: {ZAID: weight_fraction}} or None.
    /// With as_material=True fractions are keyed by nuclide name instead.
    #[pyo3(signature = (name, as_material=false))]
    #[allow(clippy::type_complexity)]
    fn get(&self, name: &str, as_material: bool) -> PyResult<Option<BTreeMap<String, Py<PyAny>>>> {
        let entry = match self.inner.get(name) {
            Some(e) => e,
            None => return Ok(None),
        };
        // Material conversion needs no GIL; do it before attaching.
        let named_fractions = if as_material {
            Some(
                entry
                    .to_material()
                    .map_err(|e| PyValueError::new_err(e.to_string()))?,
            )
        } else {
            None
        };

        Ok(Python::attach(|py| {
            let mut d: BTreeMap<String, Py<PyAny>> = BTreeMap::new();
            d.insert(
                "name".into(),
                entry
                    .name
                    .as_str()
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            d.insert(
                "mat_num".into(),
                entry.mat_num.into_pyobject(py).unwrap().unbind().into_any(),
            );
            d.insert(
                "density".into(),
                entry.density.into_pyobject(py).unwrap().unbind().into_any(),
            );
            match &named_fractions {
                Some(mat) => {
                    let fr: BTreeMap<String, f64> =
                        mat.comp.iter().map(|(id, g)| (id.to_name(), *g)).collect();
                    d.insert(
                        "fractions".into(),
                        fr.into_pyobject(py).unwrap().unbind().into_any(),
                    );
                }
                None => {
                    let fr = entry.weight_fractions();
                    d.insert(
                        "fractions".into(),
                        fr.into_pyobject(py).unwrap().unbind().into_any(),
                    );
                }
            }
            Some(d)
        }))
    }
}

/// Python module entry point.
#[pymodule]
fn _internal(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(version, m)?)?;
    m.add_function(wrap_pyfunction!(from_zaid, m)?)?;
    m.add_function(wrap_pyfunction!(atomic_mass, m)?)?;
    m.add_function(wrap_pyfunction!(natural_abundance, m)?)?;
    m.add_function(wrap_pyfunction!(rxname_id, m)?)?;
    m.add_function(wrap_pyfunction!(rxname_name, m)?)?;
    m.add_function(wrap_pyfunction!(rxname_mt, m)?)?;
    m.add_function(wrap_pyfunction!(read_xsdir, m)?)?;
    m.add_function(wrap_pyfunction!(read_meshtal, m)?)?;
    m.add_function(wrap_pyfunction!(read_wwinp, m)?)?;
    m.add_function(wrap_pyfunction!(read_mctal, m)?)?;
    m.add_function(wrap_pyfunction!(read_ssw, m)?)?;
    m.add_function(wrap_pyfunction!(read_ptrac, m)?)?;
    m.add_function(wrap_pyfunction!(read_chain, m)?)?;
    m.add_function(wrap_pyfunction!(deplete, m)?)?;
    m.add_function(wrap_pyfunction!(read_serpent, m)?)?;
    m.add_function(wrap_pyfunction!(read_usrbin, m)?)?;
    m.add_function(wrap_pyfunction!(magic, m)?)?;
    m.add_function(wrap_pyfunction!(write_ssw, m)?)?;
    m.add_function(wrap_pyfunction!(mesh_to_geom, m)?)?;
    m.add_function(wrap_pyfunction!(half_life, m)?)?;
    m.add_function(wrap_pyfunction!(decay_constant, m)?)?;
    m.add_function(wrap_pyfunction!(q_value_capture, m)?)?;
    m.add_function(wrap_pyfunction!(q_value_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(read_inp, m)?)?;
    m.add_function(wrap_pyfunction!(from_formula, m)?)?;
    m.add_function(wrap_pyfunction!(activity, m)?)?;
    m.add_class::<PyNuclide>()?;
    m.add_class::<PyParticle>()?;
    m.add_class::<PyXsdir>()?;
    m.add_class::<PyXsdirTable>()?;
    m.add_class::<PyMeshtal>()?;
    m.add_class::<PyMeshTally>()?;
    m.add_class::<PyWwinp>()?;
    m.add_class::<PyMctal>()?;
    m.add_class::<PySurfSrc>()?;
    m.add_class::<PyPtracFile>()?;
    m.add_class::<PyChain>()?;
    m.add_class::<PyUsrbinTally>()?;
    m.add_class::<PyMagicOutput>()?;
    m.add_class::<PyAliasTable>()?;
    m.add_class::<PyMeshSourceSampler>()?;
    m.add_class::<PyCascade>()?;
    m.add_class::<PyMaterialsCompendium>()?;
    Ok(())
}
