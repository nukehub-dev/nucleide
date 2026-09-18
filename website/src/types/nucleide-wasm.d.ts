// Hand-written type declarations for the wasm-pack generated Nucleide module.
// The module is loaded dynamically at runtime from /wasm/nucleide_wasm.js.

export interface NuclideInfo {
  name: string;
  nucid: number;
  z: number;
  a: number;
  state: number;
  zzaaam: number;
  zaid: number;
  zzllaaam: string;
  serpent: string;
  nist: string;
  cinder: number;
  alara: string;
  sza: number;
  mass?: number;
  abundance?: number;
}

export interface WasmNuclide {
  toObject(): NuclideInfo;
  name: string;
  nucid: number;
  z: number;
  a: number;
  state: number;
  zzaaam: number;
  zaid: number;
  zzllaaam: string;
  serpent: string;
  nist: string;
  cinder: number;
  alara: string;
  sza: number;
  mass?: number;
  abundance?: number;
  fluka(): string;
}

export interface WasmMaterial {
  mass: number;
  density?: number;
  weightFractions(): Record<string, number>;
  atomFractions(): Record<string, number>;
  toXml(name: string, density: number): string;
}

export interface CascadeResult {
  alpha: number;
  Mstar: number;
  feedAssay: number;
  productAssay: number;
  tailsAssay: number;
  stagesEnriching: number;
  stagesStripping: number;
  swuPerFeed: number;
  swuPerProduct: number;
  productPerFeed: number;
  tailsPerFeed: number;
  feed: Record<string, number>;
  product: Record<string, number>;
  tails: Record<string, number>;
}

export interface StagePointJson {
  stage: number;
  assayJ: number;
}

export interface WasmCascade {
  solve(): void;
  solveMulticomponent(): void;
  toObject(): CascadeResult;
  stageProfile(): StagePointJson[];
  alpha: number;
  feedAssay: number;
  productAssay: number;
  tailsAssay: number;
  stagesEnriching: number;
  stagesStripping: number;
  swuPerFeed: number;
  swuPerProduct: number;
}

export interface WasmChain {
  nuclides(): string[];
}

export interface DeckCellJson {
  num: number;
  mat: number;
  dens?: number;
  geom: string;
  params: string[];
}

export interface DeckSurfJson {
  num: number;
  reflecting: boolean;
  transform?: number;
  periodic?: number;
  kind: string;
  coeffs: number[];
}

export interface ModeJson {
  particles: string[];
}

export interface TransformJson {
  number: number;
  displacement: [number, number, number];
  rotation: number[];
  inDegrees: boolean;
  mainToAux: boolean;
  hidden: boolean;
}

export interface UniverseJson {
  number: number;
  cells: number[];
  notTruncated: number[];
}

export interface LatticeJson {
  cell: number;
  lattice: number;
}

export interface FillJson {
  cell: number;
  kind: "single" | "matrix";
  universe?: number;
  minIndex?: [number, number, number];
  maxIndex?: [number, number, number];
  universes?: (number | null)[];
  transform?: number;
  hiddenTransform?: number[];
  inDegrees: boolean;
}

export interface ImportanceJson {
  cell: number;
  particle: string;
  value: number;
}

export interface VolumeJson {
  cell: number;
  volume: number;
}

export interface TallyJson {
  number: number;
  type: number;
  particles: string[];
  entries: string[];
  fm?: string[];
  eBins?: string[];
}

export interface WasmDeckProblem {
  dumps(): string;
  message: string;
  title: string;
  cells(): DeckCellJson[];
  surfs(): DeckSurfJson[];
  materialNumbers(): number[];
  dataNames(): string[];
  mode(): ModeJson;
  transforms(): TransformJson[];
  universes(): UniverseJson[];
  lattices(): LatticeJson[];
  fills(): FillJson[];
  importances(): ImportanceJson[];
  volumes(): VolumeJson[];
  tallies(): TallyJson[];
  cellInventory(): Record<number, number>;
  validate(): void;
  validationNotes(): string[];
  setCellDensity(cell: number, dens: number): void;
  setCellMaterial(cell: number, mat: number): void;
  setMode(particles: string[]): void;
  setCellUniverse(cell: number, universe: number, notTruncated?: boolean): void;
  setCellLattice(cell: number, lattice?: number): void;
  setCellFill(cell: number, universe: number): void;
}

export interface WasmInventory {
  numbers(): Record<string, number>;
  decay(
    dt: number,
    timeUnit?: string,
    rates?: Record<string, number>,
    order?: number,
    method?: string,
  ): WasmInventory;
  activities(units: string): Record<string, number>;
  masses(units: string): Record<string, number>;
  moles(units: string): Record<string, number>;
  activityFractions(): Record<string, number>;
  massFractions(): Record<string, number>;
  moleFractions(): Record<string, number>;
  halfLivesReadable(): Record<string, string>;
  add(other: WasmInventory): WasmInventory;
  sub(other: WasmInventory): WasmInventory;
  mul(s: number): WasmInventory;
  div(s: number): WasmInventory;
  toCsv(): string;
}

export interface DepleteSeriesResult {
  times: number[];
  atoms: Record<string, number>[];
  activity: Record<string, number>[];
  decay_heat: Record<string, number>[];
}

export interface McnpMaterialJson {
  number: number;
  fractions: Record<string, number>;
  fractionType: "atom" | "mass";
  density?: number;
  comments: string[];
}

export interface DecayBranchJson {
  progeny: string;
  branching_fraction: number;
  mode: string;
}

export interface XsdirTableJson {
  name: string;
  zaid: string;
  serpentType?: string;
  awr: number;
  filename: string;
  fileType: number;
  temperature?: number;
  metastable?: boolean;
}

export interface XsdirSummary {
  datapath?: string;
  awrCount: number;
  tableCount: number;
  tables: XsdirTableJson[];
}

export interface MeshTallySummary {
  tallyNumber: number;
  particle: string;
  doseResponse: boolean;
  dims: number[];
  numVes: number;
  numEGroups: number;
  xBounds: number[];
  yBounds: number[];
  zBounds: number[];
  eBounds: number[];
  result: number[][];
  relError: number[][];
  totalResult: number[];
  totalRelError: number[];
}

export interface MeshtalSummary {
  version: string;
  title: string;
  histories: number;
  tallyCount: number;
  tallies: Record<string, MeshTallySummary>;
}

export interface WwinpSummary {
  ni: number;
  nr: number;
  ne: number[];
  nf: number[];
  origin: number[];
  nc: number[];
  bounds: number[][];
  e: number[][];
  ww: number[][][];
}

export interface MctalTallySummary {
  number: number;
  particleType: number;
  detectorType: number | null;
  bins: number[];
  pairs: number;
  total: number;
}

export interface MctalSummary {
  codeName: string;
  codeVersion: string;
  nHistories: number;
  tallyNums: number[];
  npert: string | null;
  tallies: MctalTallySummary[];
  nCycles: number;
}

export interface MagicSummary {
  lowerBoundsWw: number[];
  groupsPerVe: number;
  scaleFactors: number[];
  eUpperBounds: number[];
  wwTagName: string;
  eUpperBoundsTagName: string;
}

export interface SampledVoxelSummary {
  index: number;
  i: number;
  j: number;
  k: number;
  weight: number;
}

export interface AlaraMixtureJson {
  name: string;
  entries: string[];
}

export interface AlaraFluxJson {
  name: string;
  file: string;
  scale: number;
  skip: number;
  format: string;
}

export interface AlaraScheduleJson {
  name: string;
  items: string[][];
}

export interface AlaraDeckSummary {
  block_kinds: string[];
  mixtures: AlaraMixtureJson[];
  fluxes: AlaraFluxJson[];
  cooling_times_s: number[];
  schedules: AlaraScheduleJson[];
}

export interface ResponseRowJson {
  time_s: number;
  time_label: string;
  nuclide: string;
  half_life_s: number;
  run_lbl: string;
  block: string;
  block_name: string;
  block_num: number;
  variable: string;
  var_unit: string;
  value: number;
}

export interface AlaraOutputSummary {
  rows: ResponseRowJson[];
  variables: string[];
  blocks: string[];
}

export interface FispactOutputSummary {
  rows: ResponseRowJson[];
  variables: string[];
}

export interface R2sStepJson {
  zone: string;
  flux: string;
}

export interface R2sSummary {
  steps: R2sStepJson[];
  cooling_s: number[];
  top_schedule: string;
  total_s: number;
}

export interface IsotxsNuclideJson {
  label: string;
  zaid: string;
  groups: number;
  total_xs: number[];
}

export interface IsotxsSummary {
  nuclides: IsotxsNuclideJson[];
  groups: number;
}

export interface PartisnZone {
  id: number;
  material: string;
  isotxs_labels: string[];
  density: number;
}

export interface PartisnDeck {
  title: string;
  dim: number;
  zones: PartisnZone[];
  source?: string | null;
}

export interface SerpentVariableJson {
  name: string;
  kind: "scalar" | "vector" | "matrix";
  shape: string;
  value?: number | string;
}

export interface SerpentResSummary {
  variable_count: number;
  version?: string;
  title?: string;
  keff?: number[];
  /** `[mean, err]` per burnup block — the full IMP_KEFF matrix. */
  keff_history?: [number, number][];
  variables: SerpentVariableJson[];
}

export interface SerpentDepSummary {
  variable_count: number;
  nuclides: string[];
  zai: number[];
  variables: SerpentVariableJson[];
}

export interface SerpentDetSpectrum {
  name: string;
  /** Energy-bin midpoint column of the matching DET<name>E grid (empty when
   *  the grid does not pair 1:1 with the value rows). */
  energy_mid: number[];
  /** Per-bin tally values. */
  values: number[];
  /** Per-bin relative errors. */
  errors: number[];
}

export interface SerpentDetSummary {
  variable_count: number;
  detectors: string[];
  spectra: SerpentDetSpectrum[];
  variables: SerpentVariableJson[];
}

export interface UsrbinTallyJson {
  name: string;
  particle: string;
  coord_sys: string;
  dims: number[];
  x_bounds: number[];
  y_bounds: number[];
  z_bounds: number[];
  part_data: number[];
  error_data: number[];
}

export interface UsrbinSummary {
  tally_count: number;
  tallies: UsrbinTallyJson[];
}

export interface OrigenTape5StepJson {
  flux: number;
  days: number;
}

export interface OrigenTape5MaterialJson {
  name: string;
  entries: { nuclide: string; grams: number }[];
}

export interface OrigenTape5Summary {
  titles: string[];
  steps: OrigenTape5StepJson[];
  materials: OrigenTape5MaterialJson[];
}

export interface OrigenTape6RecordJson {
  nuclide: string;
  grams: number;
  activity_bq: number;
}

export interface OrigenTape6Summary {
  total_activity_bq: number;
  records: OrigenTape6RecordJson[];
}

export interface OrigenTape9Summary {
  entries: { nuclide: string; decay_const: number }[];
}

export interface DroppedJson {
  nuclide: string;
  mass: number;
  reason: string;
}

export interface DriftRowJson {
  code: string;
  massIn: number;
  massOut: number;
  relDrift: number;
  dropped: DroppedJson[];
  reparsed: boolean;
}

export interface EmitOpts {
  mcnpNumber?: number;
  xsSuffix?: string;
  serpentLib?: string;
  flukaFid?: number;
  partisnZone?: number;
}

export interface SnapshotZoneJson {
  id: string;
  volumeCm3: number;
  zbottomCm?: number;
  ztopCm?: number;
  material?: string;
  xsType?: string;
  temperatureC?: number;
  composition: Record<string, number>;
  flux?: string;
}

export interface SnapshotFluxJson {
  name: string;
  file: string;
  scale: number;
}

export interface SnapshotInputJson {
  zones: SnapshotZoneJson[];
  fluxDefs: SnapshotFluxJson[];
  coolingS: number[];
  scheduleText?: string;
  output?: string;
}

export interface SnapshotBundleJson {
  workflow: R2sSummary;
  deck: string;
  decks: string[];
}

export interface VoxelTagsInputJson {
  totals: number[];
  zoneOfVoxel: number[];
  split?: boolean;
}

export interface VoxelTagsSummary {
  n_zones: number;
  n_voxels: number;
  zone_of_voxel: number[];
  source_strength: number[];
  decay_time_s: number[];
  zone_totals: number[];
  total: number;
}

export interface VoxelPhotonGroupJson {
  nuclide: string;
  time_s: number;
  strengths: number[];
}

export interface VoxelPhotonInputJson {
  photonText: string;
  nuclides: string[];
  timeS: number;
}

export interface VoxelPhotonSummary {
  groups: VoxelPhotonGroupJson[];
  sums: number[];
  total: number;
}

export interface CompendiumEntryInfo {
  name: string;
  acronym: string[];
  mat_num: number;
  density: number;
  atom_density: number;
  source: string;
  comment: string[];
  weight_fractions: Record<string, number>;
}

export interface KineticsTransientResult {
  times: number[];
  n: number[];
  promptJump: number | null;
  betaTotal: number;
}

export interface TritiumBreakthroughResult {
  times: number[];
  fluxOverJss: number[];
  tLag: number;
  jss: number;
}

export interface SpectroscopySmoothResult {
  smoothed: number[];
  gross: number;
  background: number;
  net: number;
}

export interface UqSampleResult {
  samples: number[][];
  sampleMean: number[];
  sampleCov: number[][];
  method: string;
  minEigen: number | null;
  maxEigen: number | null;
}

export interface McplParticleJson {
  ekin: number;
  position: [number, number, number];
  direction: [number, number, number];
  time: number;
  weight: number;
  pdgcode: number;
  userflags: number;
}

export interface McplSummary {
  version: number;
  nparticles: number;
  srcname: string;
  comments: string[];
  hasUserflags: boolean;
  hasPolarisation: boolean;
  doublePrec: boolean;
  universalPdgcode?: number;
  universalWeight?: number;
  blobs: { key: string; len: number }[];
  particles: McplParticleJson[];
  truncated: boolean;
}

export interface SpeSummary {
  spec_name: string;
  channels: number;
  startChan: number;
  liveTime: number;
  realTime: number;
  detId: string;
  energyFit: number[];
  counts: number[];
  ebins: number[];
  truncated: boolean;
}

export interface RtfluxSummary {
  kind: string;
  groups: number;
  npoints: number;
  values: number[];
  truncated: boolean;
}

export interface CusumResult {
  alarmed: boolean;
  statistic: number;
  mean: number;
  std: number;
  count: number;
}

export interface FusionSpectrumMoments {
  reaction: string;
  nominalMeV: number;
  meanMeV: number;
  sigmaMeV: number;
  mono: boolean;
}

export interface FusionParticleJson {
  positionCm: [number, number, number];
  direction: [number, number, number];
  energyMeV: number;
  weight: number;
}

export interface FusionSampleResult {
  kind: string;
  count: number;
  particles: FusionParticleJson[];
}

export interface FusionRingSpec {
  kind: "ring";
  radiusCm: number;
  heightCm: number;
  reaction: string;
  tiKev: number;
  n: number;
  seed: number;
}

export interface FusionPointSpec {
  kind: "point";
  xCm: number;
  yCm: number;
  zCm: number;
  reaction: string;
  tiKev: number;
  n: number;
  seed: number;
}

export interface FusionParametricSpec {
  kind: "parametric";
  majorRadiusCm: number;
  minorRadiusCm: number;
  elongation: number;
  triangularity: number;
  shafranovFactorCm: number;
  mode: string;
  fuel: string;
  centreDensityM3: number;
  densityPeaking: number;
  pedestalDensityM3: number;
  separatrixDensityM3: number;
  centreTempKev: number;
  tempPeaking: number;
  tempBeta: number;
  pedestalTempKev: number;
  separatrixTempKev: number;
  pedestalRadiusCm: number;
  fuelDeuterium?: number;
  fuelTritium?: number;
  tailFraction?: number;
  tailTempKev?: number;
  speciesDeuteriumKev?: number;
  speciesTritiumKev?: number;
  n: number;
  seed: number;
}

export interface LatticePointSpec {
  positionCm: [number, number, number];
  rate: number;
  tiKev: number;
}

export interface LatticeSourceSpec {
  points: LatticePointSpec[];
  reaction: string;
  fieldPeriods?: number;
  baseAngle?: number;
  n: number;
  seed: number;
}

export interface LatticeSampleResult {
  kind: string;
  totalStrength: number;
  count: number;
  particles: FusionParticleJson[];
}

export interface EcrhScalars {
  gyrofrequencyGhz: number;
  resonantFieldT: number;
  relativisticFieldT?: number;
  o1CutoffDensityM3: number;
  x1CutoffDensityM3?: number;
}

export interface EcrhAccessibility {
  resonantFieldT: number;
  relativisticFieldT?: number;
  o1CutoffDensityM3: number;
  resonanceM: number[];
  o1CutoffM: number[];
  x1CutoffM: number[];
}

export interface TbrScalars {
  rawTbr: number;
  effectiveTbr: number;
  margin: number;
  meets?: boolean;
  burnGPerDay: number;
  surplusGPerDay: number;
}

export interface CoilFlux {
  fastFlux: number;
  fastFluence: number;
}

export interface CoilLifetime {
  seconds: number;
  limiting?: number;
}

export interface SubletActivityEntry {
  nuclide: string;
  activityBq: number;
  irt: number;
  alphaFrac?: number;
}

export interface SubletActivity {
  totalBq: number;
  alphaBq: number;
  betaBq: number;
  gammaBq: number;
  exTritiumBq: number;
}

export interface SubletHeatEntry {
  nuclide: string;
  activityBq: number;
  eAlphaEv: number;
  eBetaEv: number;
  eGammaEv: number;
}

export interface SubletHeat {
  alphaKw: number;
  betaKw: number;
  gammaKw: number;
  totalKw: number;
  exTritiumKw: number;
}

export interface SubletHazardEntry {
  nuclide: string;
  activityBq: number;
  coeffSvPerBq: number;
}

export interface SubletHazard {
  totalSv: number;
  exTritiumSv: number;
}

export interface SubletTransportEntry {
  nuclide: string;
  activityBq: number;
  a2Tbq: number;
}

export interface SubletTransport {
  ratio: number;
  totalBq: number;
  effectiveA2Tbq: number;
}

export interface SubletIaeaEntry {
  nuclide: string;
  activityBq: number;
  limitBqPerKg: number;
}

export interface SubletIaea {
  index: number;
  clearanceClass: string;
  maxFraction: number;
  maxNuclide: string | null;
}

export interface SubletDoseGroup {
  intensity: number;
  muAir: number;
  mu: number;
}

export interface SubletDose {
  doseSvPerH: number;
}

export interface SubletPointDose {
  doseSvPerH: number;
  distanceUsedM: number;
  clamped: boolean;
}

export interface RatioUq {
  metric: string;
  nominal: number;
  mean: number;
  std: number;
  expected: number;
  analyticStd: number;
  q16: number;
  q50: number;
  q84: number;
  expectedQ16: number;
  expectedQ50: number;
  expectedQ84: number;
  quantilesPassed: boolean;
  k: number;
  n: number;
  seed: number;
  passed: boolean;
}

export type IndataValue = { Int: number } | { Float: number } | { Bool: boolean } | { Str: string };

export interface IndataDoc {
  scalars: Record<string, IndataValue>;
  indexed: Record<string, { index: number[]; value: IndataValue }[]>;
}

export interface WallLoadResult {
  ntheta: number;
  nzeta: number;
  loads: number[][];
  total: number;
}

export type FusionSourceSpec = FusionRingSpec | FusionPointSpec | FusionParametricSpec;

export interface FusionCardsSpec {
  kind: "ring" | "point" | "parametric";
  radiusCm?: number;
  heightCm?: number;
  xCm?: number;
  yCm?: number;
  zCm?: number;
  reaction?: string;
  tiKev?: number;
  majorRadiusCm?: number;
  minorRadiusCm?: number;
  elongation?: number;
  triangularity?: number;
  shafranovFactorCm?: number;
  mode?: string;
  fuel?: string;
  centreDensityM3?: number;
  densityPeaking?: number;
  pedestalDensityM3?: number;
  separatrixDensityM3?: number;
  centreTempKev?: number;
  tempPeaking?: number;
  tempBeta?: number;
  pedestalTempKev?: number;
  separatrixTempKev?: number;
  pedestalRadiusCm?: number;
  nBins?: number;
  mcnpVersion?: number;
}

export interface FusionCardsResult {
  mcnp: string;
  serpent: string;
}

export interface SandiiSolutionJson {
  spectrum: number[];
  rates: number[];
  rateFactors: number[];
  iterations: number;
  tolerance: number;
  maxRelChange: number;
}

export interface ClearanceTableEntry {
  nuclide: string;
  limitBqG: number;
}

export interface ClearanceTableJson {
  count: number;
  entries: ClearanceTableEntry[];
}

export interface SumOfFractionsJson {
  sum: number;
  class: "satisfied" | "exceeded";
  maxFraction: number;
  maxNuclide: string | null;
}

export interface WasmMaterialsCompendium {
  len: number;
  is_empty: boolean;
  names(): string[];
  get(name: string): CompendiumEntryInfo;
}

export interface WasmApi {
  default: () => Promise<void>;
  WasmNuclide: {
    new (name: string): WasmNuclide;
    fromZzaaam(v: number): WasmNuclide;
    fromNucid(s: string): WasmNuclide;
  };
  WasmMaterial: {
    new (formula: string): WasmMaterial;
    fromAtomFrac(atoms: Record<string, number>): WasmMaterial;
    mixByMass(parts: { formula: string; fraction: number }[]): WasmMaterial;
  };
  WasmCascade: {
    defaultUranium(): WasmCascade;
    new (config: unknown): WasmCascade;
  };
  WasmChain: {
    fromXml(xml: string): WasmChain;
  };
  WasmDeckProblem: {
    fromText(text: string): WasmDeckProblem;
  };
  WasmInventory: {
    new (chain: WasmChain, comp: Record<string, number>, units?: string): WasmInventory;
    fromCsv(chain: WasmChain, text: string): WasmInventory;
  };
  cumulativeDecays(
    chain: WasmChain,
    n0: Record<string, number>,
    dt: number,
    rates?: Record<string, number>,
  ): Record<string, number>;
  progeny(chain: WasmChain, name: string): [string, number, string][];
  branchingFraction(chain: WasmChain, parent: string, child: string): number | undefined;
  decayMode(chain: WasmChain, parent: string, child: string): string | undefined;
  chainEdges(chain: WasmChain): [string, string, number, string][];
  WasmMaterialsCompendium: {
    fromJson(text: string): WasmMaterialsCompendium;
  };
  atomicMass(key: string): number | undefined;
  naturalAbundance(key: string): number | undefined;
  halfLife(key: string): number | undefined;
  decayConstant(key: string): number | undefined;
  qValueCapture(key: string): number | undefined;
  qValueAlpha(key: string): number | undefined;
  normalize_nuclide(name: string): string;
  decay_branches(key: string): DecayBranchJson[];
  decay_branch_fraction(parent: string, progeny: string): number | undefined;
  deplete(
    chain: WasmChain,
    n0: Record<string, number>,
    dt: number,
    rates: Record<string, number>,
    order: number,
    method?: string,
  ): Record<string, number>;
  depleteSeries(
    chain: WasmChain,
    n0: Record<string, number>,
    dts: number[],
    rates: Record<string, number>,
    integrator: string,
    order: number,
    method?: string,
  ): DepleteSeriesResult;
  parseMcnpMaterials(text: string): McnpMaterialJson[];
  parseXsdir(text: string): XsdirSummary;
  parseMeshtal(text: string): MeshtalSummary;
  parseWwinp(text: string): WwinpSummary;
  parseMctal(text: string): MctalSummary;
  magicBounds(
    meshtalText: string,
    tallyNumber: number,
    selection: "total" | "perGroup",
    tolerance: number,
    nullValue: number,
  ): MagicSummary;
  aliasTableSample(pdf: number[], r1: number, r2: number): number;
  meshSourceSample(
    meshtalText: string,
    tallyNumber: number,
    mode: "analog" | "uniform",
    r1: number,
    r2: number,
  ): SampledVoxelSummary;
  parseAlaraDeck(text: string): AlaraDeckSummary;
  parseAlaraOutput(text: string, runLbl: string): AlaraOutputSummary;
  parseFispactOutput(text: string, runLbl: string): FispactOutputSummary;
  parseOrigenTape5(text: string): OrigenTape5Summary;
  parseOrigenTape6(text: string): OrigenTape6Summary;
  parseOrigenTape9(text: string): OrigenTape9Summary;
  r2sFromDeck(text: string): R2sSummary;
  r2sFromSnapshot(snapshot: SnapshotInputJson): SnapshotBundleJson;
  voxelTagsFromTotals(input: VoxelTagsInputJson): VoxelTagsSummary;
  voxelPhotonSums(input: VoxelPhotonInputJson): VoxelPhotonSummary;
  parseIsotxs(text: string): IsotxsSummary;
  partisnRender(deck: PartisnDeck): string;
  partisnValidate(deck: PartisnDeck, isotxsText: string): void;
  parseSerpentRes(text: string): SerpentResSummary;
  parseSerpentDep(text: string): SerpentDepSummary;
  parseSerpentDet(text: string): SerpentDetSummary;
  parseUsrbin(text: string): UsrbinSummary;
  emitCards(
    comp: Record<string, number>,
    name: string,
    density?: number,
    opts?: EmitOpts,
  ): Record<string, string>;
  emitDriftTable(
    comp: Record<string, number>,
    name: string,
    density?: number,
    opts?: EmitOpts,
  ): DriftRowJson[];
  emitArmiCards(
    comp: Record<string, number>,
    name: string,
    density?: number,
    opts?: EmitOpts,
  ): Record<string, string>;
  emitArmiDriftTable(
    comp: Record<string, number>,
    name: string,
    density?: number,
    opts?: EmitOpts,
  ): DriftRowJson[];
  doseFactor(name: string, pathway: string, source?: string): number | undefined;
  dosePerGram(comp: Record<string, number>, pathway: string, source?: string): number;
  kineticsTransient(
    betas: number[],
    lambdas: number[],
    lambdaGen: number,
    tStep: number,
    rhoInit: number,
    rhoFinal: number,
    times: number[],
    n0: number,
  ): KineticsTransientResult;
  tritiumBreakthrough(
    length: number,
    diffusivity: number,
    c0: number,
    times: number[],
  ): TritiumBreakthroughResult;
  spectroscopySmooth(
    counts: number[],
    method: string,
    c1: number,
    c2: number,
  ): SpectroscopySmoothResult;
  parseLinesTsv(text: string): [number, number][];
  energyBins(channels: number[], fit: number[]): number[];
  detectorEfficiency(energyMev: number, coeff: number[], effFit: number): number;
  parseDollarSpe(text: string): SpeSummary;
  parsePlainSpe(text: string): SpeSummary;
  inhourRho(betas: number[], lambdas: number[], lambdaGen: number, omega: number): number;
  stablePeriod(betas: number[], lambdas: number[], lambdaGen: number, rho: number): number;
  promptJump(nBefore: number, rhoBefore: number, rhoAfter: number, betaTotal: number): number;
  parseRtflux(text: string, kind: string): RtfluxSummary;
  materialSeparate(
    comp: Record<string, number>,
    effs: Record<string, number>,
  ): { product: Record<string, number>; tails: Record<string, number> };
  materialBlend(parts: { comp: Record<string, number>; ratio: number }[]): Record<string, number>;
  cusumDetect(series: number[], k?: number, h?: number, startup?: number): CusumResult;
  readMcpl(bytes: Uint8Array): McplSummary;
  writeMcpl(header: unknown, particles: unknown): Uint8Array;
  ssw2mcpl(sswBytes: Uint8Array, surfs: number[], kinds: string[], options?: unknown): Uint8Array;
  mcpl2ssw(mcplBytes: Uint8Array, referenceSswBytes: Uint8Array, surface?: number): Uint8Array;
  uqSample(mean: number[], cov: number[][], n: number, seed: number): UqSampleResult;
  sampleLhs(mean: number[], cov: number[][], n: number, seed: number): UqSampleResult;
  damageNrtDpa(flux: number[], response: number[], bounds: number[], seconds: number): number;
  damageArcDpa(flux: number[], response: number[], bounds: number[], seconds: number): number;
  damageGasAppm(flux: number[], response: number[], bounds: number[], seconds: number): number;
  damageHeDpaRatio(
    flux: number[],
    heResponse: number[],
    damageResponse: number[],
    bounds: number[],
    seconds: number,
  ): number;
  damageEnergy(tEv: number, target: string): number;
  nrtDisplacements(tEv: number, edEv: number, target: string): number;
  arcEfficiency(tDamEv: number, edEv: number, bArc: number, cArc: number): number;
  fusionSpectrumMoments(reaction: string, tiKev: number): FusionSpectrumMoments;
  fusionReactivity(reaction: string, tiKev: number): number;
  sampleFusionSource(spec: FusionSourceSpec): FusionSampleResult;
  emitFusionSourceCards(spec: FusionCardsSpec): FusionCardsResult;
  ecrhScalars(frequencyGhz: number, harmonic: number, bT: number, teKev?: number): EcrhScalars;
  ecrhAccessibility(
    sM: number[],
    bT: number[],
    neM3: number[],
    frequencyGhz: number,
    harmonic: number,
    teKev?: number,
  ): EcrhAccessibility;
  sampleLatticeSource(spec: LatticeSourceSpec): LatticeSampleResult;
  tbrScalars(
    tritonsBred: number,
    sourceNeutrons: number,
    portFractions: number[],
    fusionPowerMw: number,
    requiredTbr?: number,
  ): TbrScalars;
  coilFastFlux(flux: number[], bounds: number[], thresholdMev: number, seconds: number): CoilFlux;
  coilLifetime(limits: number[], rates: number[]): CoilLifetime;
  subletActivity(entries: SubletActivityEntry[]): SubletActivity;
  subletDecayHeat(entries: SubletHeatEntry[]): SubletHeat;
  subletIngestionHazard(entries: SubletHazardEntry[]): SubletHazard;
  subletInhalationHazard(entries: SubletHazardEntry[]): SubletHazard;
  subletTransportRatio(entries: SubletTransportEntry[]): SubletTransport;
  subletIaeaClearance(totalMassKg: number, entries: SubletIaeaEntry[]): SubletIaea;
  subletDoseSlab(activityBqPerKg: number, groups: SubletDoseGroup[]): SubletDose;
  subletDosePoint(
    activityBqPerKg: number,
    sourceMassKg: number,
    distanceM: number,
    groups: SubletDoseGroup[],
  ): SubletPointDose;
  subletDoseMixtureMu(fractions: number[], elementMus: number[][]): number[];
  damageHeDpaRatioUq(
    flux: number[],
    heResponse: number[],
    damageResponse: number[],
    bounds: number[],
    seconds: number,
    mean: number[],
    cov: number[][],
    n: number,
    seed: number,
    k: number,
  ): RatioUq;
  parseIndata(text: string): IndataDoc;
  wallLoad(
    sEdges: number[],
    ntheta: number,
    nzeta: number,
    nfp: number,
    birth: number[],
    jacobian: number[],
  ): WallLoadResult;
  unfoldForwardFold(response: number[][], spectrum: number[]): number[];
  sandiiSolve(
    response: number[][],
    rates: number[],
    guess: number[],
    tolerance?: number,
    maxIterations?: number,
  ): SandiiSolutionJson;
  euClearanceTable(): ClearanceTableJson;
  clearanceIndex(inventory: Record<string, number>): number;
  clearanceSumOfFractions(inventory: Record<string, number>): SumOfFractionsJson;
}
