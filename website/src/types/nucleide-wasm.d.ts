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
  ): Record<string, number>;
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
}
