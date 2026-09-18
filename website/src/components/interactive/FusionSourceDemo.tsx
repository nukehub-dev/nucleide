import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type {
  FusionCardsResult,
  FusionParametricSpec,
  FusionRingSpec,
  FusionSampleResult,
  FusionSpectrumMoments,
  LatticeSampleResult,
  LatticeSourceSpec,
} from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";

// Tokamak midplane ring at R = 300 cm, T_i = 20 keV D-T: the crate's example
// configuration (synthetic source geometry, no evaluated data).
const DEFAULT_RADIUS = "300";
const DEFAULT_HEIGHT = "0";
const DEFAULT_TI = "20";
const DEFAULT_N = "2000";
const DEFAULT_SEED = "42";
const MAX_SAMPLES = 5000;

type Reaction = "dt" | "dd";

type SourceKind = "ring" | "lattice" | "parametric";

// Synthetic five-node base-sector cloud: an arc at R = 300 cm spanning
// ±4° with unit rates (no evaluated data, no machine geometry).
const LATTICE_PRESET: { positionCm: [number, number, number]; rate: number }[] = [
  { positionCm: [299.27, -20.92, 0], rate: 1 },
  { positionCm: [299.82, -10.47, 0], rate: 1 },
  { positionCm: [300, 0, 0], rate: 1 },
  { positionCm: [299.82, 10.47, 0], rate: 1 },
  { positionCm: [299.27, 20.92, 0], rate: 1 },
];
const DEFAULT_PERIODS = "5";
const DEFAULT_BASE_ANGLE = "0";

// Synthetic ITER-ish H-mode parametric preset (hand round numbers, no
// machine data): fixed Miller geometry and profiles, editable core below.
const PARAM_PRESET = {
  majorRadiusCm: 620,
  minorRadiusCm: 200,
  elongation: 1.85,
  triangularity: 0.35,
  shafranovFactorCm: 15,
  mode: "H",
  fuel: "dt",
  centreDensityM3: 1.2e20,
  densityPeaking: 1.1,
  pedestalDensityM3: 4.0e19,
  separatrixDensityM3: 3.0e19,
  tempPeaking: 2.5,
  tempBeta: 2.0,
  pedestalTempKev: 4.0,
  separatrixTempKev: 0.1,
  pedestalRadiusCm: 150,
};
const DEFAULT_PARAM_TI = "20";
const DEFAULT_FUEL_D = "0.5";
const DEFAULT_FUEL_T = "0.5";
const DEFAULT_TAIL_FRAC = "0.05";
const DEFAULT_TAIL_TEMP = "60";

const REACTIONS: { value: Reaction; label: string }[] = [
  { value: "dt", label: "D-T (14.1 MeV)" },
  { value: "dd", label: "D-D (2.45 MeV)" },
];

export function FusionSourceDemo() {
  const { wasm, ready, error } = useWasm();
  const [sourceKind, setSourceKind] = useState<SourceKind>("ring");
  const [reaction, setReaction] = useState<Reaction>("dt");
  const [radiusText, setRadiusText] = useState(DEFAULT_RADIUS);
  const [heightText, setHeightText] = useState(DEFAULT_HEIGHT);
  const [tiText, setTiText] = useState(DEFAULT_TI);
  const [nText, setNText] = useState(DEFAULT_N);
  const [seedText, setSeedText] = useState(DEFAULT_SEED);
  const [sample, setSample] = useState<FusionSampleResult | null>(null);
  const [lattice, setLattice] = useState<LatticeSampleResult | null>(null);
  const [periodsText, setPeriodsText] = useState(DEFAULT_PERIODS);
  const [baseAngleText, setBaseAngleText] = useState(DEFAULT_BASE_ANGLE);
  const [paramTiText, setParamTiText] = useState(DEFAULT_PARAM_TI);
  const [fuelDText, setFuelDText] = useState(DEFAULT_FUEL_D);
  const [fuelTText, setFuelTText] = useState(DEFAULT_FUEL_T);
  const [tailFracText, setTailFracText] = useState(DEFAULT_TAIL_FRAC);
  const [tailTempText, setTailTempText] = useState(DEFAULT_TAIL_TEMP);
  const [moments, setMoments] = useState<FusionSpectrumMoments | null>(null);
  const [cards, setCards] = useState<FusionCardsResult | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function parseEntry(text: string, label: string): number {
    const v = parseFloat(text);
    if (!Number.isFinite(v)) throw new Error(`bad ${label} \`${text}\``);
    return v;
  }

  function parseCountSeed(): { n: number; seed: number } {
    const n = parseEntry(nText, "n");
    if (!Number.isInteger(n) || n < 1 || n > MAX_SAMPLES)
      throw new Error(`bad n \`${nText}\` (expected an integer in [1, ${MAX_SAMPLES}])`);
    const seed = parseEntry(seedText, "seed");
    if (!Number.isInteger(seed) || seed < 0 || seed >= 2 ** 53)
      throw new Error(`bad seed \`${seedText}\` (expected a non-negative integer below 2^53)`);
    return { n, seed };
  }

  function buildSpec(): FusionRingSpec {
    const { n, seed } = parseCountSeed();
    const radiusCm = parseEntry(radiusText, "radius");
    if (radiusCm <= 0) throw new Error(`bad radius \`${radiusText}\` (expected > 0)`);
    const heightCm = parseEntry(heightText, "height");
    const tiKev = parseEntry(tiText, "ion temperature");
    if (tiKev < 0) throw new Error(`bad ion temperature \`${tiText}\` (expected >= 0)`);
    return {
      kind: "ring",
      radiusCm,
      heightCm,
      reaction,
      tiKev,
      n,
      seed,
    };
  }

  function buildLatticeSpec(): LatticeSourceSpec {
    const { n, seed } = parseCountSeed();
    const tiKev = parseEntry(tiText, "ion temperature");
    if (tiKev < 0) throw new Error(`bad ion temperature \`${tiText}\` (expected >= 0)`);
    const fieldPeriods = parseEntry(periodsText, "field periods");
    if (!Number.isInteger(fieldPeriods) || fieldPeriods < 1)
      throw new Error(`bad field periods \`${periodsText}\` (expected an integer >= 1)`);
    const baseAngle = parseEntry(baseAngleText, "base angle");
    return {
      points: LATTICE_PRESET.map((p) => ({ ...p, tiKev })),
      reaction,
      fieldPeriods,
      baseAngle,
      n,
      seed,
    };
  }

  function buildParametricSpec(): FusionParametricSpec {
    const { n, seed } = parseCountSeed();
    const centreTempKev = parseEntry(paramTiText, "centre ion temperature");
    if (centreTempKev < 0)
      throw new Error(`bad centre ion temperature \`${paramTiText}\` (expected >= 0)`);
    const fuelDeuterium = parseEntry(fuelDText, "fuel D");
    const fuelTritium = parseEntry(fuelTText, "fuel T");
    const fracBlank = tailFracText.trim() === "";
    const tempBlank = tailTempText.trim() === "";
    let tailFraction: number | undefined;
    let tailTempKev: number | undefined;
    if (!fracBlank || !tempBlank) {
      tailFraction = parseEntry(tailFracText, "tail fraction");
      tailTempKev = parseEntry(tailTempText, "tail temperature");
    }
    return {
      kind: "parametric",
      ...PARAM_PRESET,
      centreTempKev,
      fuelDeuterium,
      fuelTritium,
      tailFraction,
      tailTempKev,
      n,
      seed,
    };
  }

  function resetOutputs() {
    setSample(null);
    setLattice(null);
    setMoments(null);
    setCards(null);
  }

  function run() {
    if (!wasm) return;
    try {
      if (sourceKind === "lattice") {
        const spec = buildLatticeSpec();
        setLattice(wasm.sampleLatticeSource(spec));
        setSample(null);
        setMoments(
          wasm.fusionSpectrumMoments(spec.reaction, parseEntry(tiText, "ion temperature")),
        );
        setCards(null);
        setLocalError(null);
        return;
      }
      if (sourceKind === "parametric") {
        const spec = buildParametricSpec();
        setSample(wasm.sampleFusionSource(spec));
        setLattice(null);
        setMoments(null);
        setCards(null);
        setLocalError(null);
        return;
      }
      const spec = buildSpec();
      const out = wasm.sampleFusionSource(spec);
      setSample(out);
      setLattice(null);
      setMoments(wasm.fusionSpectrumMoments(spec.reaction, spec.tiKev));
      setCards(wasm.emitFusionSourceCards({ ...spec, nBins: 21 }));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      resetOutputs();
    }
  }

  const displayError = error ?? localError;

  return (
    <div className="rounded-xl border border-border/50 bg-background p-4 space-y-4">
      {!ready && <p className="text-sm text-muted-foreground">Loading the live demo…</p>}
      {displayError && (
        <div className="flex items-start justify-between gap-2 rounded-lg border border-red-200 bg-red-50 p-3 text-sm text-red-700 dark:border-red-900 dark:bg-red-950 dark:text-red-300">
          <span>Demo error: {displayError}</span>
          <button
            onClick={clearError}
            className="font-bold leading-none"
            aria-label="Dismiss error"
          >
            ×
          </button>
        </div>
      )}

      {ready && (
        <>
          <div className="flex flex-wrap gap-2" role="group" aria-label="Source geometry">
            {(
              [
                { value: "ring", label: "Tokamak ring" },
                { value: "lattice", label: "3D lattice" },
                { value: "parametric", label: "Parametric + tail" },
              ] as const
            ).map((k) => (
              <Button
                key={k.value}
                onClick={() => {
                  setSourceKind(k.value);
                  resetOutputs();
                  clearError();
                }}
                variant={sourceKind === k.value ? "default" : "outline"}
                size="sm"
              >
                {k.label}
              </Button>
            ))}
          </div>

          <div className="flex flex-wrap gap-2" role="group" aria-label="Fusion reaction">
            {REACTIONS.map((r) => (
              <Button
                key={r.value}
                onClick={() => {
                  setReaction(r.value);
                  resetOutputs();
                  clearError();
                }}
                variant={reaction === r.value ? "default" : "outline"}
                size="sm"
              >
                {r.label}
              </Button>
            ))}
          </div>

          {sourceKind === "parametric" ? (
            <div className="space-y-3">
              <p className="text-sm text-muted-foreground">
                Synthetic H-mode parametric plasma (fixed Miller geometry and profiles, no machine
                data): a D/T fuel mixture with an optional deuterium hot-tail fraction. Leave the
                tail fraction blank for a pure Maxwellian mix.
              </p>
              <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-5">
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Centre Tᵢ [keV]</Label>
                  <Input
                    value={paramTiText}
                    onChange={(e) => {
                      setParamTiText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Fuel D fraction</Label>
                  <Input
                    value={fuelDText}
                    onChange={(e) => {
                      setFuelDText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Fuel T fraction</Label>
                  <Input
                    value={fuelTText}
                    onChange={(e) => {
                      setFuelTText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Tail fraction (blank = none)</Label>
                  <Input
                    value={tailFracText}
                    onChange={(e) => {
                      setTailFracText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Tail T [keV]</Label>
                  <Input
                    value={tailTempText}
                    onChange={(e) => {
                      setTailTempText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Particles n (max {MAX_SAMPLES})</Label>
                  <Input
                    value={nText}
                    onChange={(e) => {
                      setNText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Seed</Label>
                  <Input
                    value={seedText}
                    onChange={(e) => {
                      setSeedText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
              </div>
            </div>
          ) : sourceKind === "ring" ? (
            <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-5">
              <div className="space-y-1">
                <Label className="flex min-h-10 items-end">Ring radius [cm]</Label>
                <Input
                  value={radiusText}
                  onChange={(e) => {
                    setRadiusText(e.target.value);
                    clearError();
                  }}
                  className="font-mono text-xs"
                />
              </div>
              <div className="space-y-1">
                <Label className="flex min-h-10 items-end">Ring height [cm]</Label>
                <Input
                  value={heightText}
                  onChange={(e) => {
                    setHeightText(e.target.value);
                    clearError();
                  }}
                  className="font-mono text-xs"
                />
              </div>
              <div className="space-y-1">
                <Label className="flex min-h-10 items-end">Ion temperature [keV]</Label>
                <Input
                  value={tiText}
                  onChange={(e) => {
                    setTiText(e.target.value);
                    clearError();
                  }}
                  className="font-mono text-xs"
                />
              </div>
              <div className="space-y-1">
                <Label className="flex min-h-10 items-end">Particles n (max {MAX_SAMPLES})</Label>
                <Input
                  value={nText}
                  onChange={(e) => {
                    setNText(e.target.value);
                    clearError();
                  }}
                  className="font-mono text-xs"
                />
              </div>
              <div className="space-y-1">
                <Label className="flex min-h-10 items-end">Seed</Label>
                <Input
                  value={seedText}
                  onChange={(e) => {
                    setSeedText(e.target.value);
                    clearError();
                  }}
                  className="font-mono text-xs"
                />
              </div>
            </div>
          ) : (
            <div className="space-y-3">
              <p className="text-sm text-muted-foreground">
                Five-node base-sector cloud on an arc at R = 300 cm (±4°, unit rates, shared ion
                temperature below), replicated by field-period symmetry — the stellarator-style
                lattice spelling. Synthetic preset, no machine geometry.
              </p>
              <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-5">
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Ion temperature [keV]</Label>
                  <Input
                    value={tiText}
                    onChange={(e) => {
                      setTiText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Field periods</Label>
                  <Input
                    value={periodsText}
                    onChange={(e) => {
                      setPeriodsText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Base angle [rad]</Label>
                  <Input
                    value={baseAngleText}
                    onChange={(e) => {
                      setBaseAngleText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Particles n (max {MAX_SAMPLES})</Label>
                  <Input
                    value={nText}
                    onChange={(e) => {
                      setNText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Seed</Label>
                  <Input
                    value={seedText}
                    onChange={(e) => {
                      setSeedText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
              </div>
            </div>
          )}

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Sample source</Button>
          </div>

          {moments && (
            <p className="text-sm">
              Spectrum: <span className="font-mono">{moments.reaction}</span> line at{" "}
              <span className="font-mono">{moments.nominalMeV.toFixed(4)} MeV</span>, Ballabio mean{" "}
              <span className="font-mono">{moments.meanMeV.toFixed(4)} MeV</span>, σ{" "}
              <span className="font-mono">
                {moments.mono ? "0 (monoenergetic)" : `${moments.sigmaMeV.toExponential(4)} MeV`}
              </span>
            </p>
          )}

          {sample && (
            <div className="space-y-2">
              <p className="text-sm font-medium">
                Sampled birth energies ({sample.count} particles)
              </p>
              <Plotly
                aspect="video"
                data={[
                  {
                    type: "histogram",
                    name: "energy",
                    x: sample.particles.map((p) => p.energyMeV),
                    nbinsx: 60,
                  },
                ]}
                layout={{
                  xaxis: { title: { text: "Birth energy (MeV)" }, type: "linear" },
                  yaxis: { title: { text: "Count" }, type: "linear" },
                  margin: { t: 16, r: 24, b: 48, l: 64 },
                  legend: { orientation: "h", y: -0.25 },
                }}
              />
              <p className="text-xs text-muted-foreground">
                {sourceKind === "parametric"
                  ? "Birth positions follow the parametric density profile over the plasma volume; directions are isotropic. Identical inputs reproduce the identical stream (seeded)."
                  : "Birth positions sit on the ring (fixed radius and height, uniform azimuth); directions are isotropic. Identical inputs reproduce the identical stream (seeded)."}
              </p>
            </div>
          )}

          {lattice && (
            <div className="space-y-2">
              <p className="text-sm font-medium">
                Lattice sample ({lattice.count} particles, total strength{" "}
                <span className="font-mono">{lattice.totalStrength.toPrecision(6)}</span>)
              </p>
              <Plotly
                aspect="video"
                data={[
                  {
                    type: "histogram",
                    name: "energy",
                    x: lattice.particles.map((p) => p.energyMeV),
                    nbinsx: 60,
                  },
                ]}
                layout={{
                  xaxis: { title: { text: "Birth energy (MeV)" }, type: "linear" },
                  yaxis: { title: { text: "Count" }, type: "linear" },
                  margin: { t: 16, r: 24, b: 48, l: 64 },
                  legend: { orientation: "h", y: -0.25 },
                }}
              />
              <p className="text-xs text-muted-foreground">
                Birth positions follow the symmetry-expanded cloud (base sector replicated by
                field-period rotation); directions are isotropic. Total strength is the rate sum
                scaled by the field-period count.
              </p>
            </div>
          )}

          {cards && (
            <div className="space-y-2">
              <p className="text-sm font-medium">Emitted source cards</p>
              <p className="text-sm">MCNP SDEF:</p>
              <pre className="overflow-x-auto rounded-lg border border-border/50 bg-muted/30 p-3 font-mono text-xs whitespace-pre-wrap">
                {cards.mcnp}
              </pre>
              <p className="text-sm">Serpent src:</p>
              <pre className="overflow-x-auto rounded-lg border border-border/50 bg-muted/30 p-3 font-mono text-xs whitespace-pre-wrap">
                {cards.serpent}
              </pre>
            </div>
          )}
        </>
      )}
    </div>
  );
}
