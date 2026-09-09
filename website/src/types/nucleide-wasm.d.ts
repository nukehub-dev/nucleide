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
  r2sFromDeck(text: string): R2sSummary;
  parseIsotxs(text: string): IsotxsSummary;
}
