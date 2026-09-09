//! Python bindings (`nucleide._internal`).
//!
//! Thin facade only: all logic lives in workspace crates so the Rust API
//! stays usable without Python. Type stubs live in `python/nucleide/_internal.pyi`.

use std::collections::BTreeMap;
use std::str::FromStr;

use pyo3::exceptions::{PyTypeError, PyValueError};
use pyo3::prelude::*;

use nucleide_nuclei::NuclideId;

/// Package version, re-exported to Python.
#[pyfunction]
fn version() -> &'static str {
    env!("CARGO_PKG_VERSION")
}

fn wrap_nucid_err(e: nucleide_nuclei::Error) -> PyErr {
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
        nucleide_nuclei::dialects::to_zaid(self.inner)
    }

    /// zzllaaam form ("U-235").
    #[getter]
    fn zzllaaam(&self) -> String {
        nucleide_nuclei::dialects::zzllaaam(self.inner)
    }

    /// Serpent-style name ("U-235").
    #[getter]
    fn serpent(&self) -> String {
        nucleide_nuclei::dialects::serpent(self.inner)
    }

    /// NIST-style name.
    #[getter]
    fn nist(&self) -> String {
        nucleide_nuclei::dialects::nist(self.inner)
    }

    /// Cinder integer id.
    #[getter]
    fn cinder(&self) -> u32 {
        nucleide_nuclei::dialects::to_cinder(self.inner)
    }

    /// ALARA name ("u:235").
    #[getter]
    fn alara(&self) -> String {
        nucleide_nuclei::dialects::alara(self.inner)
    }

    /// SZA integer.
    #[getter]
    fn sza(&self) -> u32 {
        nucleide_nuclei::dialects::to_sza(self.inner)
    }

    /// FLUKA element-isotope name; raises ValueError if unavailable.
    fn fluka(&self) -> PyResult<&'static str> {
        nucleide_nuclei::dialects::id_to_fluka(self.inner)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Atomic mass in u (AME2020), or None if unknown.
    #[getter]
    fn mass(&self) -> Option<f64> {
        nucleide_nuclei::data::atomic_mass(self.inner.nucid())
    }

    /// Natural abundance fraction, or None.
    #[getter]
    fn abundance(&self) -> Option<f64> {
        nucleide_nuclei::data::natural_abundance(self.inner.nucid())
    }

    fn __repr__(&self) -> String {
        format!("Nuclide({})", self.inner.to_name())
    }
}

/// Parse a MCNP ZAID integer into a Nuclide.
#[pyfunction]
fn from_zaid(zaid: u32) -> PyResult<PyNuclide> {
    nucleide_nuclei::dialects::from_zaid(zaid)
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
    lookup(key, nucleide_nuclei::data::atomic_mass)
}

/// Natural abundance fraction for a nucid integer or name string.
#[pyfunction]
fn natural_abundance(key: &Bound<'_, PyAny>) -> PyResult<Option<f64>> {
    lookup(key, nucleide_nuclei::data::natural_abundance)
}

/// A particle species with cross-code name translations.
#[pyclass(name = "Particle")]
struct PyParticle {
    inner: nucleide_nuclei::particles::ParticleId,
}

#[pymethods]
impl PyParticle {
    /// Create from any alias ("n", "neutron", "gamma", PDC int, ...).
    #[new]
    fn new(spec: &Bound<'_, PyAny>) -> PyResult<Self> {
        let inner = if let Ok(pdc) = spec.extract::<i32>() {
            nucleide_nuclei::particles::ParticleId::from_pdc(pdc)
                .ok_or_else(|| PyValueError::new_err(format!("unknown PDC code {pdc}")))?
        } else if let Ok(s) = spec.extract::<&str>() {
            s.parse::<nucleide_nuclei::particles::ParticleId>()
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
    nucleide_nuclei::rxname::name_to_id(name).map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Canonical short name for a reaction id.
#[pyfunction]
fn rxname_name(id: u32) -> Option<&'static str> {
    nucleide_nuclei::rxname::id_to_name(id)
}

/// ENDF MT number for a reaction id (0 if none registered).
#[pyfunction]
fn rxname_mt(id: u32) -> i32 {
    nucleide_nuclei::rxname::id_to_mt(id)
}

// ---------------------------------------------------------------------------
// MCNP file I/O
// ---------------------------------------------------------------------------

fn io_err(e: nucleide_mcnp_io::xsdir::Error) -> PyErr {
    PyValueError::new_err(e.to_string())
}
fn m_err<T>(r: Result<T, impl std::fmt::Display>) -> PyResult<T> {
    r.map_err(|e| PyValueError::new_err(e.to_string()))
}

/// One xsdir directory entry.
#[pyclass(name = "XsdirTable")]
struct PyXsdirTable {
    inner: nucleide_mcnp_io::xsdir::XsdirTable,
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
    inner: nucleide_mcnp_io::xsdir::Xsdir,
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
    nucleide_mcnp_io::xsdir::Xsdir::from_file(path)
        .map(|inner| PyXsdir { inner })
        .map_err(io_err)
}

/// One fmesh4 tally from a meshtal file.
#[pyclass(name = "MeshTally")]
struct PyMeshTally {
    inner: nucleide_mcnp_io::meshtal::MeshTallyData,
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
    /// Per-cell energy-integrated total relative errors.
    #[getter]
    fn total_rel_error(&self) -> Vec<f64> {
        self.inner.total_rel_error.clone()
    }
    /// Full results + relative errors as nested lists (plain copy).
    ///
    /// Zero-copy `result_array()` via NumPy stays deferred (see the Stream C
    /// module note); use this until the ndarray/NumPy bridge lands.
    fn to_list(&self) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
        (self.inner.result.clone(), self.inner.rel_error.clone())
    }
    /// Per-cell energy-integrated totals + errors as flat lists (plain copy).
    fn totals_list(&self) -> (Vec<f64>, Vec<f64>) {
        (
            self.inner.total_result.clone(),
            self.inner.total_rel_error.clone(),
        )
    }
}

/// Parsed meshtal file.
#[pyclass(name = "Meshtal")]
struct PyMeshtal {
    inner: nucleide_mcnp_io::meshtal::Meshtal,
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
    m_err(nucleide_mcnp_io::meshtal::Meshtal::from_file(path).map(|inner| PyMeshtal { inner }))
}

/// Parsed WWINP weight-window file.
#[pyclass(name = "Wwinp")]
struct PyWwinp {
    inner: nucleide_mcnp_io::wwinp::Wwinp,
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
    m_err(nucleide_mcnp_io::wwinp::Wwinp::from_file(path).map(|inner| PyWwinp { inner }))
}

/// Parsed MCTAL kcode data.
#[pyclass(name = "Mctal")]
struct PyMctal {
    inner: nucleide_mcnp_io::mctal::Mctal,
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
    m_err(nucleide_mcnp_io::mctal::Mctal::from_file(path).map(|inner| PyMctal { inner }))
}

/// Parsed SSW surface-source file.
#[pyclass(name = "SurfSrc")]
struct PySurfSrc {
    inner: nucleide_mcnp_io::surfsrc::SurfSrc,
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
    nucleide_mcnp_io::surfsrc::SurfSrc::open(path)
        .map(|inner| PySurfSrc { inner })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Detected PTRAC layout: 0 = i4 little-endian, 1 = i8 little-endian.
#[pyclass(name = "PtracFile")]
struct PyPtracFile {
    inner: nucleide_mcnp_io::ptrac::PtracFile,
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
            nucleide_mcnp_io::ptrac::Format::I4LittleEndian => 0,
            nucleide_mcnp_io::ptrac::Format::I8LittleEndian => 1,
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
    nucleide_mcnp_io::ptrac::PtracFile::open(path)
        .map(|inner| PyPtracFile { inner })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

// ---------------------------------------------------------------------------
// Depletion / CRAM
// ---------------------------------------------------------------------------

/// A parsed depletion chain (XML format).
#[pyclass(name = "Chain")]
struct PyChain {
    inner: std::sync::Arc<nucleide_depletion::Chain>,
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
    nucleide_depletion::Chain::from_file(path)
        .map(|inner| PyChain {
            inner: std::sync::Arc::new(inner),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// One-group reaction rates keyed by "NuclideName:reaction".
type RateMap = BTreeMap<String, f64>;

/// Pre-built depletion system for repeated CRAM solves.
#[pyclass(name = "DepletionSystem")]
struct PyDepletionSystem {
    inner: std::sync::Arc<nucleide_depletion::DepletionSystem>,
}

#[pymethods]
impl PyDepletionSystem {
    /// Solve one depletion step with the pre-built system.
    #[pyo3(signature = (n0, dt, order=48))]
    fn solve(
        &self,
        n0: BTreeMap<String, f64>,
        dt: f64,
        order: u8,
    ) -> PyResult<BTreeMap<String, f64>> {
        let order = parse_order(order)?;
        nucleide_depletion::deplete(&self.inner, order, &n0, dt)
            .map(|r| r.atoms)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Solve one depletion step using pre-built index vectors.
    ///
    /// `n0` and the returned vector are in chain index order; this avoids the
    /// name-to-index mapping overhead of `solve()` for tight timing loops.
    #[pyo3(signature = (n0, dt, order=48))]
    fn solve_vec(&self, n0: Vec<f64>, dt: f64, order: u8) -> PyResult<Vec<f64>> {
        let order = parse_order(order)?;
        nucleide_depletion::cram(&self.inner, order, &n0, dt)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }
}

/// Build a reusable depletion system from a chain and reaction rates.
#[pyfunction]
fn build_depletion_system(chain: &PyChain, rates: RateMap) -> PyResult<PyDepletionSystem> {
    let rs = split_rates(&rates, &chain.inner)?;
    nucleide_depletion::DepletionSystem::build((*chain.inner).clone(), &rs)
        .map(|sys| PyDepletionSystem {
            inner: std::sync::Arc::new(sys),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

fn parse_order(order: u8) -> PyResult<nucleide_depletion::Order> {
    match order {
        16 => Ok(nucleide_depletion::Order::Order16),
        48 => Ok(nucleide_depletion::Order::Order48),
        other => Err(PyValueError::new_err(format!(
            "unsupported CRAM order {other}"
        ))),
    }
}

fn split_rates(
    rates: &RateMap,
    chain: &nucleide_depletion::Chain,
) -> PyResult<nucleide_depletion::ReactionRates> {
    let mut out = nucleide_depletion::ReactionRates::new();
    for (key, v) in rates {
        let (nuc, rx) = key.split_once(':').ok_or_else(|| {
            PyValueError::new_err(format!("rate key `{key}` must be `Name:reaction`"))
        })?;
        let idx = chain
            .index_of(nuc)
            .ok_or_else(|| PyValueError::new_err(format!("rate for unknown nuclide `{nuc}`")))?;
        out.entry(idx).or_default().insert(rx.to_string(), *v);
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
        16 => nucleide_depletion::Order::Order16,
        48 => nucleide_depletion::Order::Order48,
        other => {
            return Err(PyValueError::new_err(format!(
                "unsupported CRAM order {other}"
            )))
        }
    };
    let rates = split_rates(rates.as_ref().unwrap_or(&BTreeMap::new()), &chain.inner)?;
    let sys = nucleide_depletion::DepletionSystem::build((*chain.inner).clone(), &rates)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    nucleide_depletion::deplete(&sys, order, &n0, dt)
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
        "res" => nucleide_serpent_io::parse_res(&text),
        "dep" => nucleide_serpent_io::parse_dep(&text),
        "det" => nucleide_serpent_io::parse_det(&text),
        other => {
            return Err(PyValueError::new_err(format!(
                "kind must be res|dep|det, got `{other}`"
            )))
        }
    }
    .map_err(|e| PyValueError::new_err(e.to_string()))?;
    fn entry_to_py(py: Python<'_>, e: &nucleide_serpent_io::Entry) -> Py<PyAny> {
        use nucleide_serpent_io::Entry as E;
        match e {
            E::Scalar(nucleide_serpent_io::Value::Num(n)) => {
                n.into_pyobject(py).unwrap().unbind().into_any()
            }
            E::Scalar(nucleide_serpent_io::Value::Str(s)) => {
                s.into_pyobject(py).unwrap().unbind().into_any()
            }
            E::Vector(vs) => vs
                .iter()
                .map(|v| match v {
                    nucleide_serpent_io::Value::Num(n) => {
                        n.into_pyobject(py).unwrap().unbind().into_any()
                    }
                    nucleide_serpent_io::Value::Str(s) => {
                        s.into_pyobject(py).unwrap().unbind().into_any()
                    }
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
    inner: nucleide_fluka_io::usrbin::UsrbinTally,
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
    let tallies = nucleide_fluka_io::usrbin::read_usrbin_file(path)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(tallies
        .into_iter()
        .map(|inner| PyUsrbinTally { inner })
        .collect())
}

/// MAGIC weight-window output.
#[pyclass(name = "MagicOutput")]
struct PyMagicOutput {
    inner: nucleide_vr_tools::magic::MagicOutput,
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
        nucleide_vr_tools::magic::MagicSelection::PerGroup
    } else {
        nucleide_vr_tools::magic::MagicSelection::Total
    };
    let params = nucleide_vr_tools::magic::MagicParams {
        tolerance,
        ..Default::default()
    };
    nucleide_vr_tools::magic::magic_with(&tally.inner, selection, params)
        .map(|inner| PyMagicOutput { inner })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Walker alias table for discrete sampling.
#[pyclass(name = "AliasTable")]
struct PyAliasTable {
    inner: nucleide_vr_tools::sampling::AliasTable,
}

#[pymethods]
impl PyAliasTable {
    /// Build from a probability density (normalized internally).
    #[new]
    fn new(pdf: Vec<f64>) -> PyResult<Self> {
        nucleide_vr_tools::sampling::AliasTable::new(&pdf)
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
    inner: nucleide_vr_tools::sampling::MeshSourceSampler,
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
            "analog" => nucleide_vr_tools::sampling::Mode::Analog,
            "uniform" => nucleide_vr_tools::sampling::Mode::Uniform,
            "user" => nucleide_vr_tools::sampling::Mode::User,
            other => {
                return Err(PyValueError::new_err(format!(
                    "mode must be analog|uniform|user, got `{other}`"
                )))
            }
        };
        nucleide_vr_tools::sampling::MeshSourceSampler::new(&tally.inner, m, user.as_deref())
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
    let track_data: Vec<nucleide_mcnp_io::surfsrc::TrackData> = match tracks {
        Some(dict_tracks) => dict_tracks
            .iter()
            .map(|d| {
                let g = |k: &str| d.get(k).copied().unwrap_or(0.0);
                let mut record = vec![0.0f64; nucleide_mcnp_io::surfsrc::TrackData::RECORD_WIDTH];
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
                nucleide_mcnp_io::surfsrc::TrackData::from_record(record)
            })
            .collect(),
        None => ssw
            .inner
            .read_tracklist()
            .map_err(|e| PyValueError::new_err(e.to_string()))?,
    };
    let mut f = std::fs::File::create(path).map_err(|e| PyValueError::new_err(e.to_string()))?;
    nucleide_mcnp_io::surfsrc::write_to(&mut f, &header, &track_data)
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
    let opts = nucleide_mcnp_io::deck::DeckOptions {
        title_card: title_card.to_string(),
        frac_type: nucleide_mcnp_io::deck::FracType::Mass,
    };
    nucleide_mcnp_io::deck::mesh_to_geom(&x_bounds, &y_bounds, &z_bounds, &cell_materials, &opts)
}

// ---------------------------------------------------------------------------
// ALARA I/O (thin glue over `alara-io`; solver stays out of scope)
// ---------------------------------------------------------------------------

/// Parse an ALARA input deck into plain Python containers.
///
/// Returns a dict with `block_kinds` (list[str] in file order), `geometry`
/// (str | None), `mixtures` (list of {name, entries}), `fluxes` (list of
/// {name, file, scale, skip, format}), `cooling_times_s` (list[float]),
/// `schedules`, `pulse_histories`, `outputs`, and `truncation`.
#[pyfunction]
fn alara_parse_deck(py: Python<'_>, text: &str) -> PyResult<Py<PyAny>> {
    let owned = text.to_owned();
    let deck = py
        .detach(move || nucleide_alara_io::AlaraDeck::parse(&owned))
        .map_err(ala_err)?;
    Ok(deck_to_py(py, &deck))
}

fn ala_err(e: nucleide_alara_io::Error) -> PyErr {
    PyValueError::new_err(e.to_string())
}

fn deck_to_py(py: Python<'_>, deck: &nucleide_alara_io::AlaraDeck) -> Py<PyAny> {
    use pyo3::types::PyDict;
    let out = PyDict::new(py);
    let block_kinds: Vec<&str> = deck.block_kinds();
    out.set_item("block_kinds", block_kinds).ok();
    out.set_item("geometry", deck.geometry.as_ref().map(|g| g.kind.clone()))
        .ok();
    let mixtures: Vec<Py<PyAny>> = deck.mixtures.iter().map(|m| mixture_to_py(py, m)).collect();
    out.set_item("mixtures", mixtures).ok();
    let fluxes: Vec<Py<PyAny>> = deck.fluxes.iter().map(|f| fluxdef_to_py(py, f)).collect();
    out.set_item("fluxes", fluxes).ok();
    out.set_item(
        "cooling_times_s",
        deck.cooling
            .as_ref()
            .map(|c| c.times_s.clone())
            .unwrap_or_default(),
    )
    .ok();
    let schedules: Vec<Py<PyAny>> = deck
        .schedules
        .iter()
        .map(|s| {
            let d = PyDict::new(py);
            let items: Vec<Vec<String>> = s.items.iter().map(|it| it.tokens.clone()).collect();
            d.set_item("name", &s.name).ok();
            d.set_item("items", items).ok();
            d.into_any().unbind()
        })
        .collect();
    out.set_item("schedules", schedules).ok();
    let histories: Vec<Py<PyAny>> = deck
        .pulse_histories
        .iter()
        .map(|h| {
            let d = PyDict::new(py);
            let levels: Vec<Py<PyAny>> = h
                .levels
                .iter()
                .map(|l| {
                    let e = PyDict::new(py);
                    e.set_item("pulses", l.pulses).ok();
                    e.set_item("delay_s", l.delay_s).ok();
                    e.into_any().unbind()
                })
                .collect();
            d.set_item("name", &h.name).ok();
            d.set_item("levels", levels).ok();
            d.into_any().unbind()
        })
        .collect();
    out.set_item("pulse_histories", histories).ok();
    let outputs: Vec<Py<PyAny>> = deck
        .outputs
        .iter()
        .map(|o| {
            let d = PyDict::new(py);
            d.set_item("resolution", &o.resolution).ok();
            d.set_item("entries", o.entries.clone()).ok();
            d.into_any().unbind()
        })
        .collect();
    out.set_item("outputs", outputs).ok();
    out.set_item("truncation", deck.truncation.as_ref().map(|t| t.tolerance))
        .ok();
    out.into_any().unbind()
}

fn mixture_to_py(py: Python<'_>, mix: &nucleide_alara_io::deck::Mixture) -> Py<PyAny> {
    use pyo3::types::PyDict;
    let entries: Vec<Py<PyAny>> = mix
        .entries
        .iter()
        .map(|e| mixture_entry_to_py(py, e))
        .collect();
    let d = PyDict::new(py);
    d.set_item("name", &mix.name).ok();
    d.set_item("entries", entries).ok();
    d.into_any().unbind()
}

fn mixture_entry_to_py(py: Python<'_>, entry: &nucleide_alara_io::deck::MixtureEntry) -> Py<PyAny> {
    use nucleide_alara_io::deck::MixtureEntry as E;
    use pyo3::types::PyDict;
    let d = PyDict::new(py);
    match entry {
        E::Material {
            name,
            rel_density,
            vol_fraction,
        } => {
            d.set_item("kind", "material").ok();
            d.set_item("name", name).ok();
            d.set_item("rel_density", *rel_density).ok();
            d.set_item("vol_fraction", *vol_fraction).ok();
        }
        E::Element {
            symbol,
            rel_density,
            vol_fraction,
        } => {
            d.set_item("kind", "element").ok();
            d.set_item("symbol", symbol).ok();
            d.set_item("rel_density", *rel_density).ok();
            d.set_item("vol_fraction", *vol_fraction).ok();
        }
        E::Like {
            mixture,
            rel_density,
        } => {
            d.set_item("kind", "like").ok();
            d.set_item("mixture", mixture).ok();
            d.set_item("rel_density", *rel_density).ok();
        }
        E::Target { target_kind, name } => {
            d.set_item("kind", "target").ok();
            d.set_item("target_kind", target_kind).ok();
            d.set_item("name", name).ok();
        }
    }
    d.into_any().unbind()
}

fn fluxdef_to_py(py: Python<'_>, flux: &nucleide_alara_io::deck::FluxDef) -> Py<PyAny> {
    use pyo3::types::PyDict;
    let d = PyDict::new(py);
    d.set_item("name", &flux.name).ok();
    d.set_item("file", &flux.file).ok();
    d.set_item("scale", flux.scale).ok();
    d.set_item("skip", flux.skip).ok();
    d.set_item("format", &flux.format).ok();
    d.into_any().unbind()
}

/// Parse an ALARA default-format group-flux file into plain containers.
///
/// Returns a dict with `name`, `groups_per_interval`, `num_intervals`,
/// `totals` (per-interval sums), `total` (grand sum), and `intervals`.
#[pyfunction]
fn alara_parse_flux(py: Python<'_>, text: &str, name: &str) -> PyResult<Py<PyAny>> {
    let owned_text = text.to_owned();
    let owned_name = name.to_owned();
    let spectra = py
        .detach(move || nucleide_alara_io::FluxSpectra::parse(&owned_name, &owned_text))
        .map_err(ala_err)?;
    use pyo3::types::PyDict;
    let d = PyDict::new(py);
    d.set_item("name", spectra.name.clone()).ok();
    d.set_item("groups_per_interval", spectra.groups_per_interval)
        .ok();
    d.set_item("num_intervals", spectra.num_intervals()).ok();
    let totals: Vec<f64> = spectra.intervals.iter().map(|iv| iv.iter().sum()).collect();
    d.set_item("totals", totals).ok();
    d.set_item("total", spectra.total()).ok();
    d.set_item("intervals", spectra.intervals.clone()).ok();
    Ok(d.into_any().unbind())
}

/// Parse an ALARA activation-output listing into a list of row dicts.
///
/// Each row carries the 11 `ResponseRow` fields as plain floats/strings/ints:
/// `time_s`, `time_label`, `nuclide`, `half_life_s`, `run_lbl`, `block`,
/// `block_name`, `block_num`, `variable`, `var_unit`, `value`.
#[pyfunction]
fn alara_parse_output(
    py: Python<'_>,
    text: &str,
    run_lbl: &str,
) -> PyResult<Vec<BTreeMap<String, Py<PyAny>>>> {
    let owned_text = text.to_owned();
    let owned_lbl = run_lbl.to_owned();
    let rows = py
        .detach(move || {
            nucleide_alara_io::output::ResponseFrame::parse(&owned_text, &owned_lbl).map(|f| f.rows)
        })
        .map_err(ala_err)?;
    Ok(rows
        .iter()
        .map(|r| {
            let mut d = BTreeMap::new();
            d.insert(
                "time_s".to_string(),
                r.time_s.into_pyobject(py).unwrap().unbind().into_any(),
            );
            d.insert(
                "time_label".to_string(),
                r.time_label
                    .clone()
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            d.insert(
                "nuclide".to_string(),
                r.nuclide
                    .clone()
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            d.insert(
                "half_life_s".to_string(),
                r.half_life_s.into_pyobject(py).unwrap().unbind().into_any(),
            );
            d.insert(
                "run_lbl".to_string(),
                r.run_lbl
                    .clone()
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            d.insert(
                "block".to_string(),
                r.block
                    .as_str()
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            d.insert(
                "block_name".to_string(),
                r.block_name
                    .clone()
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            d.insert(
                "block_num".to_string(),
                r.block_num.into_pyobject(py).unwrap().unbind().into_any(),
            );
            d.insert(
                "variable".to_string(),
                r.variable
                    .as_str()
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            d.insert(
                "var_unit".to_string(),
                r.var_unit
                    .clone()
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            d.insert(
                "value".to_string(),
                r.value.into_pyobject(py).unwrap().unbind().into_any(),
            );
            d
        })
        .collect())
}

/// Expand a deck's schedule hierarchy into flat irradiation/cooling steps.
///
/// Choice: takes deck text (plus optional top schedule name) instead of JSON
/// schedule/history blobs, so callers reuse the already-parsed deck blocks
/// without a parallel JSON schema. Returns a list of
/// {duration_s, flux, is_cooling} dicts.
#[pyfunction]
#[pyo3(signature = (deck_text, top=None))]
fn alara_expand_schedule(
    py: Python<'_>,
    deck_text: &str,
    top: Option<&str>,
) -> PyResult<Vec<BTreeMap<String, Py<PyAny>>>> {
    let owned_text = deck_text.to_owned();
    let owned_top = top.map(str::to_owned);
    let steps = py
        .detach(move || expand_deck_schedules(&owned_text, owned_top.as_deref()))
        .map_err(PyValueError::new_err)?;
    Ok(steps
        .into_iter()
        .map(|s| {
            let mut d = BTreeMap::new();
            let cooling = s.is_cooling();
            d.insert(
                "duration_s".to_string(),
                s.duration_s.into_pyobject(py).unwrap().unbind().into_any(),
            );
            d.insert(
                "flux".to_string(),
                s.flux
                    .clone()
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            d.insert(
                "is_cooling".to_string(),
                pyo3::types::PyBool::new(py, cooling)
                    .to_owned()
                    .into_any()
                    .unbind(),
            );
            d
        })
        .collect())
}

fn expand_deck_schedules(
    deck_text: &str,
    top: Option<&str>,
) -> Result<Vec<nucleide_alara_io::FlatStep>, String> {
    let deck = nucleide_alara_io::AlaraDeck::parse(deck_text).map_err(|e| e.to_string())?;
    let mut scheds = Vec::with_capacity(deck.schedules.len());
    for raw in &deck.schedules {
        let mut items = Vec::with_capacity(raw.items.len());
        for entry in &raw.items {
            items.push(
                parse_deck_sched_item(&entry.tokens)
                    .map_err(|m| format!("schedule `{}` line {}: {m}", raw.name, entry.line))?,
            );
        }
        scheds.push(nucleide_alara_io::schedule::ScheduleDef {
            name: raw.name.clone(),
            items,
        });
    }
    let histories: Vec<nucleide_alara_io::schedule::PulseHistory> = deck
        .pulse_histories
        .iter()
        .map(|h| nucleide_alara_io::schedule::PulseHistory {
            name: h.name.clone(),
            levels: h
                .levels
                .iter()
                .map(|l| nucleide_alara_io::schedule::PulseLevel {
                    count: l.pulses,
                    delay_s: l.delay_s,
                })
                .collect(),
        })
        .collect();
    match top {
        Some(name) => {
            nucleide_alara_io::expand_from(name, &scheds, &histories).map_err(|e| e.to_string())
        }
        None => nucleide_alara_io::expand(&scheds, &histories).map_err(|e| e.to_string()),
    }
}

fn parse_deck_sched_item(tokens: &[String]) -> Result<nucleide_alara_io::SchedItem, String> {
    match tokens {
        [op_text, op_unit, flux, history, delay_text, delay_unit] => {
            let op: f64 = op_text
                .parse()
                .map_err(|_| format!("expected operating time, found `{op_text}`"))?;
            let delay: f64 = delay_text
                .parse()
                .map_err(|_| format!("expected delay, found `{delay_text}`"))?;
            let op_time_s =
                nucleide_alara_io::parse_time_to_seconds(op, op_unit).map_err(|e| e.to_string())?;
            let delay_s = nucleide_alara_io::parse_time_to_seconds(delay, delay_unit)
                .map_err(|e| e.to_string())?;
            Ok(nucleide_alara_io::SchedItem::Pulse {
                op_time_s,
                flux: flux.clone(),
                history: history.clone(),
                delay_s,
            })
        }
        [name, history, delay_text, delay_unit] => {
            let delay: f64 = delay_text
                .parse()
                .map_err(|_| format!("expected delay, found `{delay_text}`"))?;
            let delay_s = nucleide_alara_io::parse_time_to_seconds(delay, delay_unit)
                .map_err(|e| e.to_string())?;
            Ok(nucleide_alara_io::SchedItem::SubSchedule {
                name: name.clone(),
                history: history.clone(),
                delay_s,
            })
        }
        _ => Err(format!(
            "expected 4- or 6-token schedule item, found {}",
            tokens.join(" ")
        )),
    }
}

// ---------------------------------------------------------------------------
// Data accessors, input parsing, enrichment, materials
// ---------------------------------------------------------------------------

/// Half-life [s] for a nucid integer or name string.
#[pyfunction]
fn half_life(key: &Bound<'_, PyAny>) -> PyResult<Option<f64>> {
    lookup(key, nucleide_nuclei::data::half_life)
}

/// Decay constant lambda = ln2 / t_half [1/s].
#[pyfunction]
fn decay_constant(key: &Bound<'_, PyAny>) -> PyResult<Option<f64>> {
    lookup(key, nucleide_nuclei::data::decay_constant)
}

/// Neutron-capture Q value computed from AME2020 masses [MeV].
#[pyfunction]
fn q_value_capture(key: &Bound<'_, PyAny>) -> PyResult<Option<f64>> {
    lookup(key, nucleide_nuclei::data::q_value_neutron_capture)
}

/// Alpha-decay Q value from AME2020 masses [MeV].
#[pyfunction]
fn q_value_alpha(key: &Bound<'_, PyAny>) -> PyResult<Option<f64>> {
    lookup(key, nucleide_nuclei::data::q_value_alpha)
}

/// Parse MCNP material cards from an input deck.
/// Returns a list of dicts: {number, fractions: {NuclideName: frac},
/// fraction_type: "atom"|"mass", density, comments}.
#[pyfunction]
fn read_inp(path: &str) -> PyResult<Vec<BTreeMap<String, Py<PyAny>>>> {
    let mats = nucleide_mcnp_io::inp::materials_from_file(path)
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
                        nucleide_mcnp_io::inp::FracKind::Atom => "atom",
                        nucleide_mcnp_io::inp::FracKind::Mass => "mass",
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

fn comp_to_material(comp: BTreeMap<String, f64>) -> PyResult<nucleide_material::Material> {
    let mut mat = nucleide_material::Material::new();
    for (name, grams) in &comp {
        let id = nucleide_nuclei::NuclideId::from_name(name)
            .map_err(|e| PyValueError::new_err(format!("`{name}`: {e}")))?;
        mat.add_nuclide(id, *grams);
    }
    Ok(mat)
}

/// Expand a chemical formula into a natural-isotope composition dict
/// ({nuclide_name: atom_fraction}) using AME2020 masses + abundances.
#[pyfunction]
fn from_formula(formula: &str) -> PyResult<BTreeMap<String, f64>> {
    use nucleide_material::AbundanceProvider;
    let parsed = nucleide_material::parse_formula(formula)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    // Build a temporary element-count material then expand via abundances:
    let mut nat = Vec::new();
    for (z, count) in &parsed {
        if let Some(isotopes) = nucleide_material::NaturalAbundances.natural_isotopes(*z) {
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
    let analytics = nucleide_material::Analytics {
        masses: &nucleide_material::Ame2020,
        decays: &nucleide_material::ChainDecays,
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

/// Serialize a composition dictionary to a `<material>` XML fragment.
#[pyfunction]
fn to_xml(comp: BTreeMap<String, f64>, name: &str, density: f64, units: &str) -> PyResult<String> {
    let mat = comp_to_material(comp)?;
    mat.to_xml(name, density, units)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Enrichment cascade with numeric multicomponent solving.
#[pyclass(name = "Cascade")]
struct PyCascade {
    inner: std::sync::Mutex<nucleide_enrichment::Cascade>,
}

#[pymethods]
impl PyCascade {
    /// Natural-uranium default cascade (alpha=1.05, Mstar=236, j=U235, k=U238).
    #[staticmethod]
    fn default_uranium() -> Self {
        Self {
            inner: std::sync::Mutex::new(nucleide_enrichment::default_uranium_cascade()),
        }
    }

    /// Build a cascade from full parameters. `mat_feed` is a dict of
    /// nuclide-name strings to mass fractions.
    #[new]
    #[allow(non_snake_case)]
    #[allow(clippy::too_many_arguments)]
    fn new(
        alpha: f64,
        Mstar: f64,
        j: u32,
        k: u32,
        N: f64,
        M: f64,
        x_feed_j: f64,
        x_prod_j: f64,
        x_tail_j: f64,
        mat_feed: BTreeMap<String, f64>,
    ) -> PyResult<Self> {
        let mut feed = BTreeMap::new();
        for (name, frac) in mat_feed {
            let id = NuclideId::from_name(&name).map_err(wrap_nucid_err)?;
            feed.insert(id, frac);
        }
        let casc = nucleide_enrichment::Cascade {
            alpha,
            Mstar,
            j: NuclideId::from_nucid(j),
            k: NuclideId::from_nucid(k),
            N,
            M,
            x_feed_j,
            x_prod_j,
            x_tail_j,
            mat_feed: nucleide_enrichment::Stream::with_total_mass(feed, 1.0),
            mat_prod: nucleide_enrichment::Stream::new(),
            mat_tail: nucleide_enrichment::Stream::new(),
            l_t_per_feed: 0.0,
            swu_per_feed: 0.0,
            swu_per_prod: 0.0,
        };
        Ok(Self {
            inner: std::sync::Mutex::new(casc),
        })
    }

    /// Solve via the numeric fixed-point + secant scheme in place.
    #[pyo3(signature = (tolerance=None, max_iterations=None))]
    fn solve(&self, tolerance: Option<f64>, max_iterations: Option<u32>) -> PyResult<()> {
        let tol = tolerance.unwrap_or(nucleide_enrichment::DEFAULT_TOLERANCE);
        let iters = max_iterations.unwrap_or(nucleide_enrichment::DEFAULT_MAX_ITER);
        let mut c = self
            .inner
            .lock()
            .map_err(|_| PyValueError::new_err("cascade lock poisoned"))?;
        *c = nucleide_enrichment::solve_numeric(&c, tol, iters)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(())
    }

    /// Solve and optimize `M*` for a multicomponent feed in place.
    #[pyo3(signature = (tolerance=None, max_iterations=None))]
    fn solve_multicomponent(
        &self,
        tolerance: Option<f64>,
        max_iterations: Option<u32>,
    ) -> PyResult<()> {
        let tol = tolerance.unwrap_or(nucleide_enrichment::DEFAULT_TOLERANCE);
        let iters = max_iterations.unwrap_or(nucleide_enrichment::DEFAULT_MAX_ITER);
        let mut c = self
            .inner
            .lock()
            .map_err(|_| PyValueError::new_err("cascade lock poisoned"))?;
        *c = nucleide_enrichment::multicomponent(&c, tol, iters)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(())
    }

    #[getter]
    fn alpha(&self) -> f64 {
        self.inner.lock().unwrap().alpha
    }
    #[getter]
    #[allow(non_snake_case)]
    fn Mstar(&self) -> f64 {
        self.inner.lock().unwrap().Mstar
    }
    #[getter]
    #[allow(non_snake_case)]
    fn N(&self) -> f64 {
        self.inner.lock().unwrap().N
    }
    #[getter]
    #[allow(non_snake_case)]
    fn M(&self) -> f64 {
        self.inner.lock().unwrap().M
    }
    #[getter]
    fn x_feed_j(&self) -> f64 {
        self.inner.lock().unwrap().x_feed_j
    }
    #[getter]
    fn x_prod_j(&self) -> f64 {
        self.inner.lock().unwrap().x_prod_j
    }
    #[getter]
    fn x_tail_j(&self) -> f64 {
        self.inner.lock().unwrap().x_tail_j
    }
    #[getter]
    fn l_t_per_feed(&self) -> f64 {
        self.inner.lock().unwrap().l_t_per_feed
    }
    #[getter]
    fn swu_per_feed(&self) -> f64 {
        self.inner.lock().unwrap().swu_per_feed
    }
    #[getter]
    fn swu_per_prod(&self) -> f64 {
        self.inner.lock().unwrap().swu_per_prod
    }
    /// Feed composition as {nuclide_name: mass_fraction}.
    #[getter]
    fn mat_feed(&self) -> BTreeMap<String, f64> {
        self.inner
            .lock()
            .unwrap()
            .mat_feed
            .comp
            .iter()
            .map(|(id, frac)| (id.to_name(), *frac))
            .collect()
    }
    /// Product composition as {nuclide_name: mass_fraction}.
    #[getter]
    fn mat_prod(&self) -> BTreeMap<String, f64> {
        self.inner
            .lock()
            .unwrap()
            .mat_prod
            .comp
            .iter()
            .map(|(id, frac)| (id.to_name(), *frac))
            .collect()
    }
    /// Tails composition as {nuclide_name: mass_fraction}.
    #[getter]
    fn mat_tail(&self) -> BTreeMap<String, f64> {
        self.inner
            .lock()
            .unwrap()
            .mat_tail
            .comp
            .iter()
            .map(|(id, frac)| (id.to_name(), *frac))
            .collect()
    }
    /// Separative work per product [kg SWU/kg] from the key assays.
    fn separative_work_per_product(&self) -> f64 {
        let c = self.inner.lock().unwrap();
        nucleide_enrichment::swu_per_prod(c.x_feed_j, c.x_prod_j, c.x_tail_j)
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
    inner: nucleide_material::MaterialsLibrary,
}

#[pymethods]
impl PyMaterialsCompendium {
    /// Load from the official MaterialsCompendium.json.
    #[staticmethod]
    fn load(path: &str) -> PyResult<Self> {
        nucleide_material::MaterialsLibrary::from_file(path)
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

// ---------------------------------------------------------------------------
// CCCC I/O (thin glue over `cccc-io`; no solver)
// ---------------------------------------------------------------------------

/// Parse ISOTXS text into plain Python containers.
///
/// Returns a dict with `nuclides` (list of {label, zaid, groups, total_xs}
/// in file order).
#[pyfunction]
fn isotxs_parse(py: Python<'_>, text: &str) -> PyResult<Py<PyAny>> {
    let owned = text.to_owned();
    let lib = py
        .detach(move || nucleide_cccc_io::IsotxsLib::parse(&owned))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(isotxs_to_py(py, &lib))
}

fn isotxs_to_py(py: Python<'_>, lib: &nucleide_cccc_io::IsotxsLib) -> Py<PyAny> {
    use pyo3::types::PyDict;
    let out = PyDict::new(py);
    let nuclides: Vec<Py<PyAny>> = lib
        .nuclides
        .iter()
        .map(|n| {
            let d = PyDict::new(py);
            d.set_item("label", &n.label).ok();
            d.set_item("zaid", &n.zaid).ok();
            d.set_item("groups", n.groups).ok();
            d.set_item("total_xs", n.total_xs.clone()).ok();
            d.into_any().unbind()
        })
        .collect();
    out.set_item("nuclides", nuclides).ok();
    out.into_any().unbind()
}

/// Parse an RTFLUX/ATFLUX/RZFLUX flux file into plain containers.
///
/// `kind` selects the expected header keyword (`rtflux`|`atflux`|`rzflux`,
/// case-insensitive). Returns a dict with `kind`, `groups`, `per_point`,
/// `npoints`, `values`, and `total`.
#[pyfunction]
#[pyo3(signature = (text, kind="rtflux"))]
fn rtflux_parse(py: Python<'_>, text: &str, kind: &str) -> PyResult<Py<PyAny>> {
    let flux_kind = match kind.to_ascii_lowercase().as_str() {
        "rtflux" => nucleide_cccc_io::rtflux::FluxKind::Rtflux,
        "atflux" => nucleide_cccc_io::rtflux::FluxKind::Atflux,
        "rzflux" => nucleide_cccc_io::rtflux::FluxKind::Rzflux,
        other => {
            return Err(PyValueError::new_err(format!(
                "kind must be rtflux|atflux|rzflux, got `{other}`"
            )))
        }
    };
    let owned = text.to_owned();
    let flux = py
        .detach(move || nucleide_cccc_io::FluxFile::parse(flux_kind, &owned))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    use pyo3::types::PyDict;
    let d = PyDict::new(py);
    d.set_item("kind", flux.kind.keyword()).ok();
    d.set_item("groups", flux.groups).ok();
    d.set_item("per_point", flux.per_point).ok();
    d.set_item("npoints", flux.npoints()).ok();
    d.set_item("values", flux.values.clone()).ok();
    d.set_item("total", flux.total()).ok();
    Ok(d.into_any().unbind())
}

fn partisn_deck_from_dict(
    deck: &Bound<'_, pyo3::types::PyDict>,
) -> PyResult<nucleide_cccc_io::PartisnDeck> {
    let title: String = match deck.get_item("title")? {
        Some(v) => v
            .extract()
            .map_err(|_| PyValueError::new_err("partisn deck `title` must be str"))?,
        None => return Err(PyValueError::new_err("partisn deck missing `title`")),
    };
    let dim: u8 = match deck.get_item("dim")? {
        Some(v) => v
            .extract()
            .map_err(|_| PyValueError::new_err("partisn deck `dim` must be 1, 2, or 3"))?,
        None => return Err(PyValueError::new_err("partisn deck missing `dim`")),
    };
    let zones_value = match deck.get_item("zones")? {
        Some(v) => v,
        None => return Err(PyValueError::new_err("partisn deck missing `zones`")),
    };
    let zone_dicts: Vec<Bound<'_, pyo3::types::PyDict>> = zones_value
        .extract()
        .map_err(|_| PyValueError::new_err("partisn deck `zones` must be a list of dicts"))?;
    let mut zones = Vec::with_capacity(zone_dicts.len());
    for z in &zone_dicts {
        let id: u32 = match z.get_item("id")? {
            Some(v) => v
                .extract()
                .map_err(|_| PyValueError::new_err("partisn zone `id` must be int"))?,
            None => return Err(PyValueError::new_err("partisn zone missing `id`")),
        };
        let material: String = match z.get_item("material")? {
            Some(v) => v
                .extract()
                .map_err(|_| PyValueError::new_err("partisn zone `material` must be str"))?,
            None => return Err(PyValueError::new_err("partisn zone missing `material`")),
        };
        let isotxs_labels: Vec<String> = match z.get_item("isotxs_labels")? {
            Some(v) => v.extract().map_err(|_| {
                PyValueError::new_err("partisn zone `isotxs_labels` must be a list of str")
            })?,
            None => {
                return Err(PyValueError::new_err(
                    "partisn zone missing `isotxs_labels`",
                ))
            }
        };
        let density: f64 = match z.get_item("density")? {
            Some(v) => v
                .extract()
                .map_err(|_| PyValueError::new_err("partisn zone `density` must be float"))?,
            None => return Err(PyValueError::new_err("partisn zone missing `density`")),
        };
        zones.push(nucleide_cccc_io::partisn::PartisnZone {
            id,
            material,
            isotxs_labels,
            density,
        });
    }
    let source: Option<String> = match deck.get_item("source")? {
        Some(v) if v.is_none() => None,
        Some(v) => Some(
            v.extract()
                .map_err(|_| PyValueError::new_err("partisn deck `source` must be str or None"))?,
        ),
        None => None,
    };
    Ok(nucleide_cccc_io::PartisnDeck {
        title,
        dim,
        zones,
        source,
    })
}

/// Render a PARTISN deck dict to PARTISN input text.
///
/// Deck shape: {title: str, dim: 1|2|3, zones: [{id, material,
/// isotxs_labels, density}], source: str | None}.
#[pyfunction]
fn partisn_render(py: Python<'_>, deck: &Bound<'_, pyo3::types::PyDict>) -> PyResult<String> {
    let rust_deck = partisn_deck_from_dict(deck)?;
    Ok(py.detach(move || rust_deck.render()))
}

/// Validate a PARTISN deck dict against ISOTXS text.
///
/// Raises `ValueError` when `dim` is not 1/2/3 or a zone names an ISOTXS
/// label absent from the library.
#[pyfunction]
fn partisn_validate(
    py: Python<'_>,
    deck: &Bound<'_, pyo3::types::PyDict>,
    isotxs_text: &str,
) -> PyResult<()> {
    let rust_deck = partisn_deck_from_dict(deck)?;
    let owned = isotxs_text.to_owned();
    let lib = py
        .detach(move || nucleide_cccc_io::IsotxsLib::parse(&owned))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    rust_deck
        .validate(&lib)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

// ---------------------------------------------------------------------------
// FISPACT-II output (thin glue over `fispact-io`; reuses ResponseFrame)
// ---------------------------------------------------------------------------

fn fispact_row_to_map(
    py: Python<'_>,
    r: &nucleide_alara_io::output::ResponseRow,
) -> BTreeMap<String, Py<PyAny>> {
    let mut d = BTreeMap::new();
    d.insert(
        "time_s".to_string(),
        r.time_s.into_pyobject(py).unwrap().unbind().into_any(),
    );
    d.insert(
        "time_label".to_string(),
        r.time_label
            .clone()
            .into_pyobject(py)
            .unwrap()
            .unbind()
            .into_any(),
    );
    d.insert(
        "nuclide".to_string(),
        r.nuclide
            .clone()
            .into_pyobject(py)
            .unwrap()
            .unbind()
            .into_any(),
    );
    d.insert(
        "half_life_s".to_string(),
        r.half_life_s.into_pyobject(py).unwrap().unbind().into_any(),
    );
    d.insert(
        "run_lbl".to_string(),
        r.run_lbl
            .clone()
            .into_pyobject(py)
            .unwrap()
            .unbind()
            .into_any(),
    );
    d.insert(
        "block".to_string(),
        r.block
            .as_str()
            .into_pyobject(py)
            .unwrap()
            .unbind()
            .into_any(),
    );
    d.insert(
        "block_name".to_string(),
        r.block_name
            .clone()
            .into_pyobject(py)
            .unwrap()
            .unbind()
            .into_any(),
    );
    d.insert(
        "block_num".to_string(),
        r.block_num.into_pyobject(py).unwrap().unbind().into_any(),
    );
    d.insert(
        "variable".to_string(),
        r.variable
            .as_str()
            .into_pyobject(py)
            .unwrap()
            .unbind()
            .into_any(),
    );
    d.insert(
        "var_unit".to_string(),
        r.var_unit
            .clone()
            .into_pyobject(py)
            .unwrap()
            .unbind()
            .into_any(),
    );
    d.insert(
        "value".to_string(),
        r.value.into_pyobject(py).unwrap().unbind().into_any(),
    );
    d
}

/// Parse a FISPACT-II inventory listing into a list of row dicts.
///
/// Each row carries the 11 `ResponseRow` fields as plain floats/strings/ints:
/// `time_s`, `time_label`, `nuclide`, `half_life_s`, `run_lbl`, `block`,
/// `block_name`, `block_num`, `variable`, `var_unit`, `value`.
#[pyfunction]
fn fispact_parse_output(
    py: Python<'_>,
    text: &str,
    run_lbl: &str,
) -> PyResult<Vec<BTreeMap<String, Py<PyAny>>>> {
    let owned_text = text.to_owned();
    let owned_lbl = run_lbl.to_owned();
    let rows = py
        .detach(move || {
            nucleide_fispact_io::parse_to_frame(&owned_text, &owned_lbl).map(|f| f.rows)
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(rows.iter().map(|r| fispact_row_to_map(py, r)).collect())
}

// ---------------------------------------------------------------------------
// ORIGEN TAPE readers (thin glue over `origen-io`; scoped TAPE5/6/9)
// ---------------------------------------------------------------------------

/// Parse ORIGEN TAPE5 input-echo text into plain containers.
///
/// Returns a dict with `titles` (list[str]), `irradiation_steps`
/// (list of {flux, days}), and `materials` (list of {name, entries:
/// [{nuclide, grams}]}).
#[pyfunction]
fn origen_parse_tape5(py: Python<'_>, text: &str) -> PyResult<Py<PyAny>> {
    let owned = text.to_owned();
    let tape = py
        .detach(move || nucleide_origen_io::Tape5::parse(&owned))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    use pyo3::types::PyDict;
    let out = PyDict::new(py);
    out.set_item("titles", tape.titles.clone()).ok();
    let steps: Vec<Py<PyAny>> = tape
        .irradiation_steps
        .iter()
        .map(|s| {
            let d = PyDict::new(py);
            d.set_item("flux", s.flux).ok();
            d.set_item("days", s.days).ok();
            d.into_any().unbind()
        })
        .collect();
    out.set_item("irradiation_steps", steps).ok();
    let materials: Vec<Py<PyAny>> = tape
        .materials
        .iter()
        .map(|m| {
            let d = PyDict::new(py);
            d.set_item("name", &m.name).ok();
            let entries: Vec<Py<PyAny>> = m
                .grams
                .iter()
                .map(|(nuclide, grams)| {
                    let e = PyDict::new(py);
                    e.set_item("nuclide", nuclide).ok();
                    e.set_item("grams", *grams).ok();
                    e.into_any().unbind()
                })
                .collect();
            d.set_item("entries", entries).ok();
            d.into_any().unbind()
        })
        .collect();
    out.set_item("materials", materials).ok();
    Ok(out.into_any().unbind())
}

/// Parse ORIGEN TAPE6 output-inventory text into plain containers.
///
/// Returns a dict with `records` (list of {nuclide, grams, activity_bq} in
/// file order) and `total_activity` (sum over records).
#[pyfunction]
fn origen_parse_tape6(py: Python<'_>, text: &str) -> PyResult<Py<PyAny>> {
    let owned = text.to_owned();
    let tape = py
        .detach(move || nucleide_origen_io::Tape6::parse(&owned))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    use pyo3::types::PyDict;
    let out = PyDict::new(py);
    let records: Vec<Py<PyAny>> = tape
        .records
        .iter()
        .map(|r| {
            let d = PyDict::new(py);
            d.set_item("nuclide", &r.nuclide).ok();
            d.set_item("grams", r.grams).ok();
            d.set_item("activity_bq", r.activity_bq).ok();
            d.into_any().unbind()
        })
        .collect();
    out.set_item("records", records).ok();
    out.set_item("total_activity", tape.total_activity()).ok();
    Ok(out.into_any().unbind())
}

/// Parse ORIGEN TAPE9 decay-constant text into a list of row dicts.
///
/// Each entry is {nuclide, decay_const} in file order.
#[pyfunction]
fn origen_parse_tape9(py: Python<'_>, text: &str) -> PyResult<Vec<BTreeMap<String, Py<PyAny>>>> {
    let owned = text.to_owned();
    let entries = py
        .detach(move || nucleide_origen_io::Tape9Entry::parse(&owned))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(entries
        .iter()
        .map(|e| {
            let mut d = BTreeMap::new();
            d.insert(
                "nuclide".to_string(),
                e.nuclide
                    .clone()
                    .into_pyobject(py)
                    .unwrap()
                    .unbind()
                    .into_any(),
            );
            d.insert(
                "decay_const".to_string(),
                e.decay_const.into_pyobject(py).unwrap().unbind().into_any(),
            );
            d
        })
        .collect())
}

// ---------------------------------------------------------------------------
// R2S workflow builder (thin glue over `r2s`; no transport/activation solve)
// ---------------------------------------------------------------------------

fn r2s_workflow_to_py(py: Python<'_>, workflow: &nucleide_r2s::R2sWorkflow) -> Py<PyAny> {
    use pyo3::types::PyDict;
    let out = PyDict::new(py);
    let steps: Vec<Py<PyAny>> = workflow
        .steps
        .iter()
        .map(|s| {
            let d = PyDict::new(py);
            d.set_item("zone", &s.zone).ok();
            d.set_item("flux", &s.flux).ok();
            d.into_any().unbind()
        })
        .collect();
    out.set_item("steps", steps).ok();
    out.set_item("cooling_s", workflow.cooling_s.clone()).ok();
    out.set_item("top_schedule", &workflow.top_schedule).ok();
    out.into_any().unbind()
}

fn r2s_workflow_from_dict(
    workflow: &Bound<'_, pyo3::types::PyDict>,
) -> PyResult<nucleide_r2s::R2sWorkflow> {
    let steps_value = match workflow.get_item("steps")? {
        Some(v) => v,
        None => return Err(PyValueError::new_err("r2s workflow missing `steps`")),
    };
    let step_dicts: Vec<Bound<'_, pyo3::types::PyDict>> = steps_value
        .extract()
        .map_err(|_| PyValueError::new_err("r2s workflow `steps` must be a list of dicts"))?;
    let mut steps = Vec::with_capacity(step_dicts.len());
    for s in &step_dicts {
        let zone: String = match s.get_item("zone")? {
            Some(v) => v
                .extract()
                .map_err(|_| PyValueError::new_err("r2s step `zone` must be str"))?,
            None => return Err(PyValueError::new_err("r2s step missing `zone`")),
        };
        let flux: String = match s.get_item("flux")? {
            Some(v) => v
                .extract()
                .map_err(|_| PyValueError::new_err("r2s step `flux` must be str"))?,
            None => return Err(PyValueError::new_err("r2s step missing `flux`")),
        };
        steps.push(nucleide_r2s::R2sStep { zone, flux });
    }
    let cooling_s: Vec<f64> = match workflow.get_item("cooling_s")? {
        Some(v) => v.extract().map_err(|_| {
            PyValueError::new_err("r2s workflow `cooling_s` must be a list of float")
        })?,
        None => return Err(PyValueError::new_err("r2s workflow missing `cooling_s`")),
    };
    let top_schedule: String = match workflow.get_item("top_schedule")? {
        Some(v) => v
            .extract()
            .map_err(|_| PyValueError::new_err("r2s workflow `top_schedule` must be str"))?,
        None => return Err(PyValueError::new_err("r2s workflow missing `top_schedule`")),
    };
    Ok(nucleide_r2s::R2sWorkflow {
        steps,
        cooling_s,
        top_schedule,
    })
}

/// Derive an R2S workflow summary from an ALARA deck.
///
/// Returns a dict with `steps` (list of {zone, flux}), `cooling_s`
/// (list[float]), and `top_schedule` (str).
#[pyfunction]
fn r2s_from_deck(py: Python<'_>, deck_text: &str) -> PyResult<Py<PyAny>> {
    let owned = deck_text.to_owned();
    let workflow = py
        .detach(move || {
            let deck = nucleide_alara_io::AlaraDeck::parse(&owned)
                .map_err(|e| nucleide_r2s::Error::Invalid(e.to_string()))?;
            nucleide_r2s::R2sWorkflow::from_deck(&deck)
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(r2s_workflow_to_py(py, &workflow))
}

/// Validate an R2S workflow dict against an ALARA deck.
///
/// Raises `ValueError` when a step zone/flux is unknown or cooling histories
/// are missing.
#[pyfunction]
fn r2s_validate(
    py: Python<'_>,
    workflow: &Bound<'_, pyo3::types::PyDict>,
    deck_text: &str,
) -> PyResult<()> {
    let rust_workflow = r2s_workflow_from_dict(workflow)?;
    let owned = deck_text.to_owned();
    let deck = py
        .detach(move || nucleide_alara_io::AlaraDeck::parse(&owned))
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    rust_workflow
        .validate_against(&deck)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Expand an ALARA deck's irradiation hierarchy into flat steps via R2S.
///
/// Returns a list of {duration_s, flux, is_cooling} dicts. When `top` is
/// given it overrides the workflow's discovered top schedule.
#[pyfunction]
#[pyo3(signature = (deck_text, top=None))]
fn r2s_expand(
    py: Python<'_>,
    deck_text: &str,
    top: Option<&str>,
) -> PyResult<Vec<BTreeMap<String, Py<PyAny>>>> {
    let owned_text = deck_text.to_owned();
    let owned_top = top.map(str::to_owned);
    let steps = py
        .detach(move || {
            let deck = nucleide_alara_io::AlaraDeck::parse(&owned_text)
                .map_err(|e| nucleide_r2s::Error::Invalid(e.to_string()))?;
            let mut workflow = nucleide_r2s::R2sWorkflow::from_deck(&deck)?;
            if let Some(top) = owned_top {
                workflow.top_schedule = top;
            }
            workflow.expand(&deck, &[])
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(steps
        .into_iter()
        .map(|s| {
            let mut d = BTreeMap::new();
            let cooling = s.is_cooling();
            d.insert(
                "duration_s".to_string(),
                s.duration_s.into_pyobject(py).unwrap().unbind().into_any(),
            );
            d.insert(
                "flux".to_string(),
                s.flux.into_pyobject(py).unwrap().unbind().into_any(),
            );
            d.insert(
                "is_cooling".to_string(),
                pyo3::types::PyBool::new(py, cooling)
                    .to_owned()
                    .into_any()
                    .unbind(),
            );
            d
        })
        .collect())
}

/// Assemble a uniform-split photon source summary for `zone`.
///
/// Parses an ALARA activation-output listing, sums shutdown
/// `SpecificActivity` over the zone's nuclide rows (skipping `total`
/// aggregates), and splits the total uniformly over `groups` energy groups.
/// Returns a dict with `zone`, `groups` (list[float]), and `total`.
///
/// Approximation: the uniform split preserves only the total shutdown
/// strength; real decay photons follow the nuclide- and energy-dependent
/// lines in ALARA `.photonSrc` spectra.
#[pyfunction]
fn r2s_assemble(
    py: Python<'_>,
    output_text: &str,
    run_lbl: &str,
    zone: &str,
    groups: usize,
) -> PyResult<Py<PyAny>> {
    let owned_text = output_text.to_owned();
    let owned_lbl = run_lbl.to_owned();
    let owned_zone = zone.to_owned();
    let source = py
        .detach(move || {
            let frame = nucleide_alara_io::output::ResponseFrame::parse(&owned_text, &owned_lbl)
                .map_err(|e| nucleide_r2s::Error::Invalid(e.to_string()))?;
            Ok::<_, nucleide_r2s::Error>(nucleide_r2s::photon::assemble(
                &frame,
                &owned_zone,
                groups,
            ))
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    use pyo3::types::PyDict;
    let out = PyDict::new(py);
    out.set_item("zone", source.zone.clone()).ok();
    out.set_item("groups", source.groups.clone()).ok();
    out.set_item("total", source.total()).ok();
    Ok(out.into_any().unbind())
}

// ---------------------------------------------------------------------------
// 0.3.0 series driver, data accessors, list helpers
// ---------------------------------------------------------------------------
//
// Thin facade only (core tables and integrators live in `nucleide-nuclei` /
// `nucleide-material` / `nucleide-depletion`; nothing duplicated here):
//
// - `deplete_series` wraps the core `integrate` series (`predictor`/`cecm`/
//   `cf4`), omitting the core `t = 0` row so there is one output per step.
// - `simple_xs` / `scattering_length` / `decay_energy` / `decay_heat` are
//   thin wrappers over the vendored TSV tables + material analytics.
// - `MeshTally::to_list` / `totals_list` are plain-copy helpers. The zero-copy
//   NumPy bridge (`result_array()`, roadmap "ndarray/NumPy zero-copy") stays
//   DEFERRED: adding the `numpy` crate was judged too risky for this change
//   (native build + abi3 version matching), so no `numpy` dependency is
//   introduced here.

/// Supported `deplete_series` integrators (core `Integrator` variants).
fn parse_integrator(name: &str) -> PyResult<nucleide_depletion::Integrator> {
    use nucleide_depletion::Integrator as I;
    if name.eq_ignore_ascii_case("predictor") {
        return Ok(I::Predictor);
    }
    if name.eq_ignore_ascii_case("cecm") {
        return Ok(I::Cecm);
    }
    if name.eq_ignore_ascii_case("cf4") {
        return Ok(I::Cf4);
    }
    Err(PyValueError::new_err(format!(
        "unsupported integrator `{name}` (supported: predictor, cecm, cf4)"
    )))
}

/// Solve a multi-step depletion series with the chosen core integrator.
///
/// Thin wrapper over `nucleide_depletion::integrate`: one [`Step`] per `dt`
/// (per-step `rates`/`rates_list`, `None` meaning decay-only), `n0` keyed by
/// nuclide name. Returns a dict with `times` (cumulative seconds, one entry
/// per step — the core `t = 0` initial row is omitted so `atoms[k]`
/// matches a single `deplete` call over `dts[k]`), `atoms`, `activity`
/// ([Bq]), and `decay_heat` ([W] per nuclide via the shared chain →
/// ENDF/B-VII.1 → 0.0 energy resolution).
#[pyfunction]
#[pyo3(signature = (chain, n0, dts, rates=None, rates_list=None, integrator="predictor", order=48))]
#[allow(clippy::too_many_arguments)]
fn deplete_series(
    chain: &PyChain,
    n0: BTreeMap<String, f64>,
    dts: Vec<f64>,
    rates: Option<RateMap>,
    rates_list: Option<Vec<Option<RateMap>>>,
    integrator: &str,
    order: u8,
) -> PyResult<Py<PyAny>> {
    use nucleide_depletion::{DepletionSystem, ReactionRates, Step};
    let integrator = parse_integrator(integrator)?;
    let order = parse_order(order)?;
    if let Some(list) = &rates_list {
        if list.len() != dts.len() {
            return Err(PyValueError::new_err(format!(
                "rates_list has {} entries but dts has {}",
                list.len(),
                dts.len()
            )));
        }
    }
    if dts.is_empty() {
        return Err(PyValueError::new_err("dts must not be empty"));
    }
    // Atom vector in chain order; unknown names fail loudly like `deplete`.
    let mut n0_vec = vec![0.0; chain.inner.len()];
    for (name, value) in &n0 {
        let idx = chain.inner.index_of(name).ok_or_else(|| {
            PyValueError::new_err(format!("unknown nuclide `{name}` for this chain"))
        })?;
        n0_vec[idx] = *value;
    }
    let empty = BTreeMap::new();
    let mut steps = Vec::with_capacity(dts.len());
    for (i, dt) in dts.iter().enumerate() {
        let step_rates = rates_list
            .as_ref()
            .and_then(|list| list[i].as_ref())
            .or(rates.as_ref())
            .unwrap_or(&empty);
        let rs = split_rates(step_rates, &chain.inner)?;
        steps.push(Step::new(*dt, rs));
    }
    // Template system: `integrate` rebuilds the matrix per step from the
    // chain + step rates; the template's own rates are unused.
    let template = DepletionSystem::build((*chain.inner).clone(), &ReactionRates::new())
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    // NOTE: plain (GIL held) call by design, matching the other CRAM
    // bindings; batch sizes here are small.
    let series = nucleide_depletion::integrate(&template, &n0_vec, &steps, integrator, order)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let names: Vec<&str> = template
        .chain
        .nuclides
        .iter()
        .map(|nuc| nuc.name.as_str())
        .collect();
    let keyed = |rows: &[Vec<f64>]| -> Vec<BTreeMap<String, f64>> {
        rows.iter()
            .map(|row| {
                names
                    .iter()
                    .zip(row)
                    .map(|(name, v)| ((*name).to_string(), *v))
                    .collect()
            })
            .collect()
    };
    // Skip the t = 0 initial row: one output entry per requested step.
    let atoms = keyed(&series.atoms[1..]);
    let activity = keyed(&series.activity[1..]);
    let decay_heat = keyed(&series.decay_heat[1..]);
    let times = series.times[1..].to_vec();
    Ok(Python::attach(|py| {
        use pyo3::types::PyDict;
        let out = PyDict::new(py);
        out.set_item("times", &times).ok();
        out.set_item("atoms", &atoms).ok();
        out.set_item("activity", &activity).ok();
        out.set_item("decay_heat", &decay_heat).ok();
        out.into_any().unbind()
    }))
}

/// Thermal/fast cross sections [barn] for a nuclide name.
///
/// Screening-level values from the `nucleide-nuclei` table (thermal 2200 m/s
/// total + 14-MeV total); `None` for nuclides outside the table.
#[pyfunction]
fn simple_xs(name: &str) -> PyResult<Option<(f64, f64)>> {
    NuclideId::from_name(name).map_err(wrap_nucid_err)?;
    Ok(nucleide_nuclei::data::simple_xs_by_name(name))
}

/// Coherent scattering length [fm] for a nuclide name.
///
/// First element of the `nucleide-nuclei` (coherent, incoherent) pair;
/// `None` for nuclides outside the table.
#[pyfunction]
fn scattering_length(name: &str) -> PyResult<Option<f64>> {
    NuclideId::from_name(name).map_err(wrap_nucid_err)?;
    Ok(nucleide_nuclei::data::scattering_length_by_name(name).map(|(b_coh, _)| b_coh))
}

/// Mean decay energy per disintegration [MeV] for a nuclide name.
///
/// Screening-level placeholder values (NOT ENSDF); `None` when unknown.
#[pyfunction]
fn decay_energy(name: &str) -> PyResult<Option<f64>> {
    NuclideId::from_name(name).map_err(wrap_nucid_err)?;
    Ok(nucleide_nuclei::data::decay_energy_mev_by_name(name))
}

/// Decay heat [W] of a composition dict ({nuclide name: grams}).
///
/// Screening-level estimate via `Material::total_decay_heat` (Ame2020 masses,
/// ENDF/B-VIII.0 decay constants, placeholder decay energies). Errors when a
/// nuclide lacks mass, decay, or energy data.
#[pyfunction]
fn decay_heat(comp: BTreeMap<String, f64>) -> PyResult<f64> {
    let mat = comp_to_material(comp)?;
    let analytics = nucleide_material::Analytics {
        masses: &nucleide_material::Ame2020,
        decays: &nucleide_material::ChainDecays,
    };
    mat.total_decay_heat(&analytics, &nucleide_material::DecayEnergies)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

fn parse_dose_pathway(s: &str) -> PyResult<nucleide_material::DosePathway> {
    nucleide_material::DosePathway::parse(s)
        .ok_or_else(|| PyValueError::new_err(format!("unknown dose pathway `{s}`")))
}

fn parse_dose_source(s: &str) -> PyResult<nucleide_material::DoseSource> {
    nucleide_material::DoseSource::parse(s)
        .ok_or_else(|| PyValueError::new_err(format!("unknown dose source `{s}`")))
}

/// Raw dose factor for a nuclide name, pathway, and source.
///
/// Pathway is one of `air`/`soil`/`ingest`/`inhale` (`ext_air`/`ext_soil`
/// aliases accepted); source is one of `EPA`/`DOE`/`GENII` (default `EPA`,
/// matching PyNE source id 0). Returns `None` when the nuclide has no row;
/// GENII/DOE air resolve to `-1.0` (PyNE missing-air sentinel).
#[pyfunction]
#[pyo3(signature = (name, pathway, source="EPA"))]
fn dose_factor(name: &str, pathway: &str, source: &str) -> PyResult<Option<f64>> {
    NuclideId::from_name(name).map_err(wrap_nucid_err)?;
    let p = parse_dose_pathway(pathway)?;
    let s = parse_dose_source(source)?;
    Ok(nucleide_nuclei::data::dose_factor_by_name(name, p, s))
}

/// Total dose per gram of a composition dict ({nuclide name: grams}).
///
/// Thin wrapper over `Material::total_dose_per_g` (Ame2020 masses,
/// ENDF/B-VIII.0 decay constants, HNF-5636/PyNE dose factors). Pathway is one
/// of `air`/`soil`/`ingest`/`inhale`; source is `EPA`/`DOE`/`GENII` (default
/// `EPA`). Units follow the table: air `mrem/h per g per m^3`, soil
/// `mrem/h per g per m^2`, ingest/inhale `mrem per g`. Screening-level only —
/// not for safety decisions. Errors when a nuclide lacks mass, decay, or
/// dose data (including `-1` GENII/DOE air sentinels).
#[pyfunction]
#[pyo3(signature = (comp, pathway, source="EPA"))]
fn dose_per_g(comp: BTreeMap<String, f64>, pathway: &str, source: &str) -> PyResult<f64> {
    let mat = comp_to_material(comp)?;
    let analytics = nucleide_material::Analytics {
        masses: &nucleide_material::Ame2020,
        decays: &nucleide_material::ChainDecays,
    };
    let p = parse_dose_pathway(pathway)?;
    let s = parse_dose_source(source)?;
    mat.total_dose_per_g(&analytics, &nucleide_material::DoseFactors, p, s)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

// ---------------------------------------------------------------------------
// 0.3.0 Tier 1: deck round-trip, decay inventories, ARMI dialects, checks
// ---------------------------------------------------------------------------

/// A parsed MCNP input deck with format-preserving write-back.
#[pyclass(name = "DeckProblem")]
struct PyDeckProblem {
    inner: std::sync::Mutex<nucleide_mcnp_io::problem::DeckProblem>,
}

fn deck_cell_dict(cell: &nucleide_mcnp_io::cell::CellCard) -> BTreeMap<String, String> {
    let mut d = BTreeMap::new();
    d.insert("num".to_string(), cell.num.to_string());
    d.insert("mat".to_string(), cell.mat.to_string());
    d.insert(
        "dens".to_string(),
        cell.dens.map(|v| v.to_string()).unwrap_or_default(),
    );
    d.insert("geom".to_string(), cell.geom.render());
    d.insert("params".to_string(), cell.params.join(" "));
    d
}

#[pymethods]
impl PyDeckProblem {
    /// Parse a deck from text.
    #[staticmethod]
    fn loads(text: &str) -> PyResult<Self> {
        nucleide_mcnp_io::problem::parse_deck(text)
            .map(|inner| Self {
                inner: std::sync::Mutex::new(inner),
            })
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Message (first) line.
    #[getter]
    fn message(&self) -> String {
        self.inner.lock().unwrap().message.clone()
    }

    /// Title card (second line).
    #[getter]
    fn title(&self) -> String {
        self.inner.lock().unwrap().title.clone()
    }

    /// Cell cards as `{num, mat, dens, geom, params}` dicts (`dens` is `""`
    /// for void cells).
    #[getter]
    fn cells(&self) -> Vec<BTreeMap<String, String>> {
        self.inner
            .lock()
            .unwrap()
            .cells
            .iter()
            .map(deck_cell_dict)
            .collect()
    }

    /// Surface cards as `{num, reflecting, transform, kind, coeffs}` dicts
    /// (`transform` is `""` when absent).
    #[getter]
    fn surfs(&self) -> Vec<BTreeMap<String, String>> {
        self.inner
            .lock()
            .unwrap()
            .surfs
            .iter()
            .map(|s| {
                let mut d = BTreeMap::new();
                d.insert("num".to_string(), s.num.to_string());
                d.insert("reflecting".to_string(), s.reflecting.to_string());
                d.insert(
                    "transform".to_string(),
                    s.transform.map(|v| v.to_string()).unwrap_or_default(),
                );
                d.insert("kind".to_string(), s.kind.keyword().to_string());
                d.insert(
                    "coeffs".to_string(),
                    s.coeffs
                        .iter()
                        .map(|v| v.to_string())
                        .collect::<Vec<_>>()
                        .join(" "),
                );
                d
            })
            .collect()
    }

    /// Material numbers in file order.
    #[getter]
    fn material_numbers(&self) -> Vec<u32> {
        self.inner
            .lock()
            .unwrap()
            .materials
            .iter()
            .map(|m| m.number)
            .collect()
    }

    /// Data-card names in file order (`MODE`, `M1`, `KCODE`, ...).
    #[getter]
    fn data_names(&self) -> Vec<String> {
        self.inner
            .lock()
            .unwrap()
            .data
            .iter()
            .map(|d| d.name.clone())
            .collect()
    }

    /// Serialize back to MCNP input text (byte-identical when unedited).
    fn dumps(&self) -> String {
        nucleide_mcnp_io::problem::write_deck(&self.inner.lock().unwrap())
    }

    /// Set a cell's density (re-renders that card canonically).
    fn set_cell_density(&self, cell: u32, dens: f64) -> PyResult<()> {
        self.inner
            .lock()
            .unwrap()
            .set_cell_density(cell, dens)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Set a cell's material number (re-renders that card canonically).
    fn set_cell_material(&self, cell: u32, mat: u32) -> PyResult<()> {
        self.inner
            .lock()
            .unwrap()
            .set_cell_material(cell, mat)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }
}

/// Parse an MCNP input deck file into a [`PyDeckProblem`].
#[pyfunction]
fn read_deck(path: &str) -> PyResult<PyDeckProblem> {
    nucleide_mcnp_io::problem::parse_deck_file(path)
        .map(|inner| PyDeckProblem {
            inner: std::sync::Mutex::new(inner),
        })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Parse MCNP input deck text into a [`PyDeckProblem`].
#[pyfunction]
fn parse_deck(text: &str) -> PyResult<PyDeckProblem> {
    PyDeckProblem::loads(text)
}

/// A unit-aware decay inventory over a depletion chain.
#[pyclass(name = "Inventory")]
struct PyInventory {
    chain: std::sync::Arc<nucleide_depletion::Chain>,
    atoms: BTreeMap<String, f64>,
}

fn inventory_sys(
    chain: &nucleide_depletion::Chain,
    rates: &RateMap,
) -> PyResult<nucleide_depletion::DepletionSystem> {
    let rs = split_rates(rates, chain)?;
    nucleide_depletion::DepletionSystem::build(chain.clone(), &rs)
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

fn parse_quantity_unit(unit: &str) -> PyResult<nucleide_depletion::QuantityUnit> {
    nucleide_depletion::QuantityUnit::from_str(unit)
        .map_err(|e| PyValueError::new_err(format!("{e:?}")))
}

#[pymethods]
impl PyInventory {
    /// Build from quantities in `units` (atom counts, `Bq`/`Ci` activity,
    /// `g`/`kg` mass, `mol`, ... — see `QuantityUnit`).
    #[new]
    #[pyo3(signature = (chain, comp, units="atoms"))]
    fn new(chain: &PyChain, comp: BTreeMap<String, f64>, units: &str) -> PyResult<Self> {
        let unit = parse_quantity_unit(units)?;
        let sys = inventory_sys(&chain.inner, &BTreeMap::new())?;
        let inv = nucleide_depletion::DecayInventory::from_units(&comp, unit, &sys)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        Ok(Self {
            chain: chain.inner.clone(),
            atoms: inv.atoms,
        })
    }

    /// Atom counts by nuclide name.
    fn numbers(&self) -> BTreeMap<String, f64> {
        self.atoms.clone()
    }

    /// Decay over `dt` in `time_unit` (`s`, `m`, `h`, `d`, `y`); optional
    /// one-group `rates` (`"Name:reaction"` keys) and CRAM `order`.
    #[pyo3(signature = (dt, time_unit="s", rates=None, order=48))]
    fn decay(&self, dt: f64, time_unit: &str, rates: Option<RateMap>, order: u8) -> PyResult<Self> {
        let order = parse_order(order)?;
        let unit = nucleide_depletion::inventory::time_unit_from_str(time_unit)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let seconds = dt * unit.as_seconds();
        let empty = BTreeMap::new();
        let step_rates = rates.as_ref().unwrap_or(&empty);
        let template = inventory_sys(&self.chain, step_rates)?;
        // Route through the core series: predictor over one step equals the
        // single CRAM solve, and rates/order stay honored.
        let steps = vec![nucleide_depletion::Step::new(
            seconds,
            split_rates(step_rates, &self.chain)?,
        )];
        let series = nucleide_depletion::integrate(
            &template,
            &chain_vec(&self.chain, &self.atoms)?,
            &steps,
            nucleide_depletion::Integrator::Predictor,
            order,
        )
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
        let names: Vec<String> = self.chain.nuclides.iter().map(|n| n.name.clone()).collect();
        let atoms = names
            .iter()
            .zip(series.atoms.last().cloned().unwrap_or_default())
            .map(|(n, v)| (n.clone(), v))
            .collect();
        Ok(Self {
            chain: self.chain.clone(),
            atoms,
        })
    }

    /// Activity per nuclide in `units`.
    fn activities(&self, units: &str) -> PyResult<BTreeMap<String, f64>> {
        let unit = parse_quantity_unit(units)?;
        let sys = inventory_sys(&self.chain, &BTreeMap::new())?;
        let inv = nucleide_depletion::DecayInventory {
            atoms: self.atoms.clone(),
        };
        inv.activities(&sys, unit)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Mass per nuclide in `units`.
    fn masses(&self, units: &str) -> PyResult<BTreeMap<String, f64>> {
        let unit = parse_quantity_unit(units)?;
        let inv = nucleide_depletion::DecayInventory {
            atoms: self.atoms.clone(),
        };
        inv.masses(unit)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Moles per nuclide in `units`.
    fn moles(&self, units: &str) -> PyResult<BTreeMap<String, f64>> {
        let unit = parse_quantity_unit(units)?;
        let inv = nucleide_depletion::DecayInventory {
            atoms: self.atoms.clone(),
        };
        inv.moles(unit)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Activity fractions by nuclide name.
    fn activity_fractions(&self) -> PyResult<BTreeMap<String, f64>> {
        let sys = inventory_sys(&self.chain, &BTreeMap::new())?;
        let inv = nucleide_depletion::DecayInventory {
            atoms: self.atoms.clone(),
        };
        inv.activity_fractions(&sys)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Mass fractions by nuclide name.
    fn mass_fractions(&self) -> PyResult<BTreeMap<String, f64>> {
        let inv = nucleide_depletion::DecayInventory {
            atoms: self.atoms.clone(),
        };
        inv.mass_fractions()
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Mole fractions by nuclide name.
    fn mole_fractions(&self) -> BTreeMap<String, f64> {
        nucleide_depletion::DecayInventory {
            atoms: self.atoms.clone(),
        }
        .mole_fractions()
    }

    /// Human-readable half-lives (`"3.2 d"`, `"stable"`, `"unknown"`).
    fn half_lives_readable(&self) -> BTreeMap<String, String> {
        nucleide_depletion::DecayInventory {
            atoms: self.atoms.clone(),
        }
        .half_lives_readable()
    }

    /// Add two inventories (atom counts sum).
    fn add(&self, other: &Self) -> Self {
        let a = nucleide_depletion::DecayInventory {
            atoms: self.atoms.clone(),
        };
        let b = nucleide_depletion::DecayInventory {
            atoms: other.atoms.clone(),
        };
        Self {
            chain: self.chain.clone(),
            atoms: a.add(&b).atoms,
        }
    }

    /// Subtract (clamped at zero).
    fn sub(&self, other: &Self) -> Self {
        let a = nucleide_depletion::DecayInventory {
            atoms: self.atoms.clone(),
        };
        let b = nucleide_depletion::DecayInventory {
            atoms: other.atoms.clone(),
        };
        Self {
            chain: self.chain.clone(),
            atoms: a.sub(&b).atoms,
        }
    }

    /// Scale by a scalar.
    fn mul(&self, scalar: f64) -> Self {
        let a = nucleide_depletion::DecayInventory {
            atoms: self.atoms.clone(),
        };
        Self {
            chain: self.chain.clone(),
            atoms: a.mul(scalar).atoms,
        }
    }

    /// Divide by a scalar.
    fn div(&self, scalar: f64) -> Self {
        let a = nucleide_depletion::DecayInventory {
            atoms: self.atoms.clone(),
        };
        Self {
            chain: self.chain.clone(),
            atoms: a.div(scalar).atoms,
        }
    }

    /// Serialize as `nuclide,atoms` CSV rows.
    fn to_csv(&self) -> String {
        nucleide_depletion::DecayInventory {
            atoms: self.atoms.clone(),
        }
        .to_csv()
    }

    /// Parse `to_csv` output back into an inventory over `chain`.
    #[staticmethod]
    fn from_csv(chain: &PyChain, text: &str) -> PyResult<Self> {
        // Validate names against the chain (core from_csv is chain-free).
        let inv = nucleide_depletion::DecayInventory::from_csv(text)
            .map_err(|e| PyValueError::new_err(e.to_string()))?;
        for name in inv.atoms.keys() {
            if chain.inner.index_of(name).is_none() {
                return Err(PyValueError::new_err(format!(
                    "unknown nuclide `{name}` for this chain"
                )));
            }
        }
        Ok(Self {
            chain: chain.inner.clone(),
            atoms: inv.atoms,
        })
    }
}

/// Atom vector in chain order for an inventory map (unknown names error).
fn chain_vec(
    chain: &nucleide_depletion::Chain,
    atoms: &BTreeMap<String, f64>,
) -> PyResult<Vec<f64>> {
    let mut vec = vec![0.0; chain.len()];
    for (name, value) in atoms {
        let idx = chain.index_of(name).ok_or_else(|| {
            PyValueError::new_err(format!("unknown nuclide `{name}` for this chain"))
        })?;
        vec[idx] = *value;
    }
    Ok(vec)
}

/// Time-integrated decays per nuclide over one step (chain order → names).
#[pyfunction]
#[pyo3(signature = (chain, n0, dt, rates=None))]
fn cumulative_decays(
    chain: &PyChain,
    n0: BTreeMap<String, f64>,
    dt: f64,
    rates: Option<RateMap>,
) -> PyResult<BTreeMap<String, f64>> {
    let empty = BTreeMap::new();
    let sys = inventory_sys(&chain.inner, rates.as_ref().unwrap_or(&empty))?;
    let vec = chain_vec(&chain.inner, &n0)?;
    let out = nucleide_depletion::cumulative_decays(&sys, &vec, dt)
        .map_err(|e| PyValueError::new_err(e.to_string()))?;
    Ok(chain
        .inner
        .nuclides
        .iter()
        .zip(out)
        .map(|(nuc, v)| (nuc.name.clone(), v))
        .collect())
}

/// `(child, branching_ratio, decay_mode)` triples for a chain nuclide.
#[pyfunction]
fn progeny(chain: &PyChain, name: &str) -> Vec<(String, f64, String)> {
    nucleide_depletion::progeny(&chain.inner, name)
}

/// Branching fraction from parent to child, if the decay exists.
#[pyfunction]
fn branching_fraction(chain: &PyChain, parent: &str, child: &str) -> Option<f64> {
    nucleide_depletion::branching_fraction(&chain.inner, parent, child)
}

/// Decay-mode label from parent to child, if the decay exists.
#[pyfunction]
fn decay_mode(chain: &PyChain, parent: &str, child: &str) -> Option<String> {
    nucleide_depletion::decay_mode(&chain.inner, parent, child)
}

/// `(parent, child, branching_ratio, decay_mode)` edges of a chain.
#[pyfunction]
fn chain_edges(chain: &PyChain) -> Vec<(String, String, f64, String)> {
    nucleide_depletion::chain_edges(&chain.inner)
}

/// Parse an ARMI nuclide label (`nU235`, `92235`, ...) into a [`PyNuclide`].
#[pyfunction]
fn armi_to_nucid(name: &str) -> PyResult<PyNuclide> {
    nucleide_nuclei::armi::armi_name_to_nucid(name)
        .map(|inner| PyNuclide { inner })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Render a nuclide in ARMI database-label form.
#[pyfunction]
fn nucid_to_armi(nuclide: &PyNuclide) -> String {
    nucleide_nuclei::armi::nucid_to_armi_label(nuclide.inner)
}

/// Parse an MCC3-style nuclide label into a [`PyNuclide`].
#[pyfunction]
fn mcc3_to_nucid(name: &str) -> PyResult<PyNuclide> {
    nucleide_nuclei::armi::mcc3_to_nucid(name)
        .map(|inner| PyNuclide { inner })
        .map_err(|e| PyValueError::new_err(e.to_string()))
}

/// Truncated-label collisions in a composition at DIF3D/MC2 widths.
///
/// `comp` maps nuclide names to grams; `widths` defaults to `[6, 8]`.
/// Returns `[{truncated, width, members}]`.
#[pyfunction]
#[pyo3(signature = (comp, widths=None))]
fn check_labels(
    comp: BTreeMap<String, f64>,
    widths: Option<Vec<usize>>,
) -> PyResult<Vec<BTreeMap<String, Py<PyAny>>>> {
    let mat = comp_to_material(comp)?;
    let widths = widths.unwrap_or_else(|| nucleide_material::DEFAULT_WIDTHS.to_vec());
    let collisions = nucleide_material::check_labels(&mat, &widths);
    Python::attach(|py| {
        Ok(collisions
            .into_iter()
            .map(|c| {
                let mut d = BTreeMap::new();
                d.insert(
                    "truncated".to_string(),
                    c.truncated.into_pyobject(py).unwrap().unbind().into_any(),
                );
                d.insert(
                    "width".to_string(),
                    c.width.into_pyobject(py).unwrap().unbind().into_any(),
                );
                let members: Vec<String> = c.members.iter().map(|id| id.to_name()).collect();
                d.insert(
                    "members".to_string(),
                    members.into_pyobject(py).unwrap().unbind().into_any(),
                );
                d
            })
            .collect())
    })
}

/// Conservation audit of a composition: `[{kind, detail}]` (empty = clean).
#[pyfunction]
fn audit_material(comp: BTreeMap<String, f64>) -> PyResult<Vec<BTreeMap<String, String>>> {
    let mat = comp_to_material(comp)?;
    Ok(nucleide_material::audit(&mat, &nucleide_material::Ame2020)
        .into_iter()
        .map(|issue| {
            let mut d = BTreeMap::new();
            d.insert("kind".to_string(), format!("{:?}", issue.kind));
            d.insert("detail".to_string(), issue.detail);
            d
        })
        .collect())
}

/// Emit one composition through all five code dialects (MCNP, Serpent, FLUKA,
/// ALARA, PARTISN). Returns `{code: card_text}`.
///
/// `comp` maps nuclide names to grams; `density` is mass density [g/cm³] for
/// dialects that need one (falls back to none — Serpent/FLUKA/PARTISN error
/// without it).
#[pyfunction]
#[pyo3(signature = (comp, name, density=None, mcnp_number=1, xs_suffix="80c", serpent_lib="03c", fluka_fid=1, partisn_zone=1))]
#[allow(clippy::too_many_arguments)]
fn emit_cards(
    comp: BTreeMap<String, f64>,
    name: &str,
    density: Option<f64>,
    mcnp_number: u32,
    xs_suffix: &str,
    serpent_lib: &str,
    fluka_fid: u32,
    partisn_zone: u32,
) -> PyResult<BTreeMap<String, String>> {
    let (emitted, _) = emit_drift_inner(
        comp,
        name,
        density,
        mcnp_number,
        xs_suffix,
        serpent_lib,
        fluka_fid,
        partisn_zone,
    )?;
    Ok(emitted
        .into_iter()
        .map(|e| (e.code.to_string(), e.text))
        .collect())
}

/// Mass-drift report for one composition across all five code dialects.
/// Returns `[{code, mass_in, mass_out, rel_drift, dropped: [{nuclide, mass,
/// reason}], reparsed}]`.
#[pyfunction]
#[pyo3(signature = (comp, name, density=None, mcnp_number=1, xs_suffix="80c", serpent_lib="03c", fluka_fid=1, partisn_zone=1))]
#[allow(clippy::too_many_arguments)]
fn emit_drift_table(
    comp: BTreeMap<String, f64>,
    name: &str,
    density: Option<f64>,
    mcnp_number: u32,
    xs_suffix: &str,
    serpent_lib: &str,
    fluka_fid: u32,
    partisn_zone: u32,
) -> PyResult<Vec<BTreeMap<String, Py<PyAny>>>> {
    let (_, table) = emit_drift_inner(
        comp,
        name,
        density,
        mcnp_number,
        xs_suffix,
        serpent_lib,
        fluka_fid,
        partisn_zone,
    )?;
    Python::attach(|py| {
        Ok(table
            .rows
            .into_iter()
            .map(|r| {
                let mut d = BTreeMap::new();
                d.insert(
                    "code".to_string(),
                    r.code
                        .to_string()
                        .into_pyobject(py)
                        .unwrap()
                        .unbind()
                        .into_any(),
                );
                d.insert(
                    "mass_in".to_string(),
                    r.mass_in.into_pyobject(py).unwrap().unbind().into_any(),
                );
                d.insert(
                    "mass_out".to_string(),
                    r.mass_out.into_pyobject(py).unwrap().unbind().into_any(),
                );
                d.insert(
                    "rel_drift".to_string(),
                    r.rel_drift.into_pyobject(py).unwrap().unbind().into_any(),
                );
                let dropped: Vec<BTreeMap<String, Py<PyAny>>> = r
                    .dropped
                    .into_iter()
                    .map(|x| {
                        let mut dd = BTreeMap::new();
                        dd.insert(
                            "nuclide".to_string(),
                            x.id.to_name()
                                .into_pyobject(py)
                                .unwrap()
                                .unbind()
                                .into_any(),
                        );
                        dd.insert(
                            "mass".to_string(),
                            x.mass.into_pyobject(py).unwrap().unbind().into_any(),
                        );
                        dd.insert(
                            "reason".to_string(),
                            x.reason.into_pyobject(py).unwrap().unbind().into_any(),
                        );
                        dd
                    })
                    .collect();
                d.insert(
                    "dropped".to_string(),
                    dropped.into_pyobject(py).unwrap().unbind().into_any(),
                );
                d.insert(
                    "reparsed".to_string(),
                    pyo3::types::PyBool::new(py, r.reparsed)
                        .to_owned()
                        .into_any()
                        .unbind(),
                );
                d
            })
            .collect())
    })
}

#[allow(clippy::too_many_arguments)]
fn emit_drift_inner(
    comp: BTreeMap<String, f64>,
    name: &str,
    density: Option<f64>,
    mcnp_number: u32,
    xs_suffix: &str,
    serpent_lib: &str,
    fluka_fid: u32,
    partisn_zone: u32,
) -> PyResult<(Vec<nucleide_emit::Emitted>, nucleide_emit::DriftTable)> {
    let mut mat = comp_to_material(comp)?;
    mat.set_density(density);
    let mut opts = nucleide_emit::EmitOptions::new(name);
    opts.mcnp_number = mcnp_number;
    opts.xs_suffix = xs_suffix.to_string();
    opts.serpent_lib = serpent_lib.to_string();
    opts.fluka_fid = fluka_fid;
    opts.partisn_zone = partisn_zone;
    nucleide_emit::emit_drift(&mat, &opts).map_err(|e| PyValueError::new_err(e.to_string()))
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
    m.add_function(wrap_pyfunction!(build_depletion_system, m)?)?;
    m.add_function(wrap_pyfunction!(deplete, m)?)?;
    m.add_function(wrap_pyfunction!(deplete_series, m)?)?;
    m.add_function(wrap_pyfunction!(simple_xs, m)?)?;
    m.add_function(wrap_pyfunction!(scattering_length, m)?)?;
    m.add_function(wrap_pyfunction!(decay_energy, m)?)?;
    m.add_function(wrap_pyfunction!(decay_heat, m)?)?;
    m.add_function(wrap_pyfunction!(dose_factor, m)?)?;
    m.add_function(wrap_pyfunction!(dose_per_g, m)?)?;
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
    m.add_function(wrap_pyfunction!(to_xml, m)?)?;
    m.add_function(wrap_pyfunction!(alara_parse_deck, m)?)?;
    m.add_function(wrap_pyfunction!(alara_parse_flux, m)?)?;
    m.add_function(wrap_pyfunction!(alara_parse_output, m)?)?;
    m.add_function(wrap_pyfunction!(alara_expand_schedule, m)?)?;
    m.add_function(wrap_pyfunction!(isotxs_parse, m)?)?;
    m.add_function(wrap_pyfunction!(rtflux_parse, m)?)?;
    m.add_function(wrap_pyfunction!(partisn_render, m)?)?;
    m.add_function(wrap_pyfunction!(partisn_validate, m)?)?;
    m.add_function(wrap_pyfunction!(fispact_parse_output, m)?)?;
    m.add_function(wrap_pyfunction!(origen_parse_tape5, m)?)?;
    m.add_function(wrap_pyfunction!(origen_parse_tape6, m)?)?;
    m.add_function(wrap_pyfunction!(origen_parse_tape9, m)?)?;
    m.add_function(wrap_pyfunction!(r2s_from_deck, m)?)?;
    m.add_function(wrap_pyfunction!(r2s_validate, m)?)?;
    m.add_function(wrap_pyfunction!(r2s_expand, m)?)?;
    m.add_function(wrap_pyfunction!(r2s_assemble, m)?)?;
    m.add_function(wrap_pyfunction!(parse_deck, m)?)?;
    m.add_function(wrap_pyfunction!(read_deck, m)?)?;
    m.add_function(wrap_pyfunction!(cumulative_decays, m)?)?;
    m.add_function(wrap_pyfunction!(progeny, m)?)?;
    m.add_function(wrap_pyfunction!(branching_fraction, m)?)?;
    m.add_function(wrap_pyfunction!(decay_mode, m)?)?;
    m.add_function(wrap_pyfunction!(chain_edges, m)?)?;
    m.add_function(wrap_pyfunction!(armi_to_nucid, m)?)?;
    m.add_function(wrap_pyfunction!(nucid_to_armi, m)?)?;
    m.add_function(wrap_pyfunction!(mcc3_to_nucid, m)?)?;
    m.add_function(wrap_pyfunction!(check_labels, m)?)?;
    m.add_function(wrap_pyfunction!(audit_material, m)?)?;
    m.add_function(wrap_pyfunction!(emit_cards, m)?)?;
    m.add_function(wrap_pyfunction!(emit_drift_table, m)?)?;
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
    m.add_class::<PyDepletionSystem>()?;
    m.add_class::<PyUsrbinTally>()?;
    m.add_class::<PyMagicOutput>()?;
    m.add_class::<PyAliasTable>()?;
    m.add_class::<PyMeshSourceSampler>()?;
    m.add_class::<PyCascade>()?;
    m.add_class::<PyMaterialsCompendium>()?;
    m.add_class::<PyDeckProblem>()?;
    m.add_class::<PyInventory>()?;
    Ok(())
}
