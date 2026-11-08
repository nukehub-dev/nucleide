import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { KdeFitResult, MagicSummary, SampledVoxelSummary } from "../../types/nucleide-wasm";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Select } from "@nukehub/docs-kit/components/ui/Select";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";

const DEFAULT_MESHTAL = `mcnp version 5.mpi ld=00000000
Demo
Number of histories used for normalizing tallies = 1000

Mesh Tally Number 4
This is a neutron mesh tally.

Tally bin boundaries:
X direction: 0.0 1.0
Y direction: 0.0 1.0
Z direction: 0.0 1.0
Energy bin boundaries: 0.0 1.0

Energy X Y Z Result Rel Error
1.00000E+00 5.00000E-01 5.00000E-01 5.00000E-01 1.00000E+00 1.00000E-01`;

const SELECTION_OPTIONS = [
  { value: "total", label: "Total" },
  { value: "perGroup", label: "Per-group" },
];

const SAMPLE_MODE_OPTIONS = [
  { value: "analog", label: "Analog" },
  { value: "uniform", label: "Uniform" },
];

const KDE_RULE_OPTIONS = [
  { value: "silverman", label: "Silverman" },
  { value: "fixed", label: "Fixed" },
  { value: "adaptive", label: "Adaptive" },
];

// Synthetic 1-D hand vector {0, 0.2, 0.4, 10} (one row per line): a tight
// cluster plus a sparse tail, so adaptive kNN bandwidths visibly differ per
// sample while Silverman gives one global width.
const DEFAULT_KDE_SAMPLES = "0\n0.2\n0.4\n10";

type KdeRule = "silverman" | "fixed" | "adaptive";

interface KdeResult {
  fit: KdeFitResult;
  draw: number[];
  xs: number[];
  ys: number[];
  firstValues: number[];
  mean: number;
  std: number;
  pdfAtMean: number;
}

export function VrDemo() {
  const { wasm, ready, error } = useWasm();
  const [meshtalText, setMeshtalText] = useState(DEFAULT_MESHTAL);
  const [tallyNumber, setTallyNumber] = useState(4);
  const [selection, setSelection] = useState<"total" | "perGroup">("total");
  const [tolerance, setTolerance] = useState(0.5);
  const [nullValue, setNullValue] = useState(0.0);
  const [magic, setMagic] = useState<MagicSummary | null>(null);

  const [pdfInput, setPdfInput] = useState("1 3 2 1");
  const [aliasResult, setAliasResult] = useState<number | null>(null);

  const [sample, setSample] = useState<SampledVoxelSummary | null>(null);
  const [sampleMode, setSampleMode] = useState<"analog" | "uniform">("analog");

  const [kdeSamplesText, setKdeSamplesText] = useState(DEFAULT_KDE_SAMPLES);
  const [kdeRule, setKdeRule] = useState<KdeRule>("silverman");
  const [kdeKText, setKdeKText] = useState("1");
  const [kdeScaleText, setKdeScaleText] = useState("1");
  const [kdeWidthsText, setKdeWidthsText] = useState("1");
  const [kdeResult, setKdeResult] = useState<KdeResult | null>(null);

  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function runMagic() {
    if (!wasm) return;
    try {
      if (!Number.isInteger(tallyNumber)) throw new Error(`bad tally number \`${tallyNumber}\``);
      if (!Number.isFinite(tolerance)) throw new Error(`bad tolerance \`${tolerance}\``);
      if (!Number.isFinite(nullValue)) throw new Error(`bad null value \`${nullValue}\``);
      const out = wasm.magicBounds(meshtalText, tallyNumber, selection, tolerance, nullValue);
      setMagic(out);
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setMagic(null);
    }
  }

  function runAlias() {
    if (!wasm) return;
    try {
      const tokens = pdfInput.split(/\s+/).filter((t) => t.length > 0);
      if (tokens.length === 0) throw new Error(`bad PDF \`${pdfInput}\``);
      const pdf = tokens.map((t) => {
        const v = parseFloat(t);
        if (!Number.isFinite(v)) throw new Error(`bad PDF value \`${t}\``);
        return v;
      });
      const r1 = Math.random();
      const r2 = Math.random();
      setAliasResult(wasm.aliasTableSample(pdf, r1, r2));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setAliasResult(null);
    }
  }

  function runSample() {
    if (!wasm) return;
    try {
      if (!Number.isInteger(tallyNumber)) throw new Error(`bad tally number \`${tallyNumber}\``);
      const r1 = Math.random();
      const r2 = Math.random();
      setSample(wasm.meshSourceSample(meshtalText, tallyNumber, sampleMode, r1, r2));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setSample(null);
    }
  }

  function parseKdeSamples(text: string): number[][] {
    const lines = text
      .split("\n")
      .map((line) => line.trim())
      .filter((line) => line.length > 0);
    if (lines.length === 0) throw new Error(`bad KDE samples \`${text}\``);
    const rows = lines.map((line) =>
      line.split(/\s+/).map((token) => {
        const value = parseFloat(token);
        if (!Number.isFinite(value)) throw new Error(`bad KDE sample value \`${token}\``);
        return value;
      }),
    );
    const dims = rows[0].length;
    for (const row of rows) {
      if (row.length !== dims) {
        throw new Error(
          `bad KDE samples \`${text}\` (expected every row to have ${dims} value(s))`,
        );
      }
    }
    return rows;
  }

  function runKde() {
    if (!wasm) return;
    try {
      const samples = parseKdeSamples(kdeSamplesText);
      const base: {
        samples: number[][];
        rule: KdeRule;
        k?: number;
        scale?: number;
        widths?: number[];
      } = { samples, rule: kdeRule };
      if (kdeRule === "adaptive") {
        const k = parseFloat(kdeKText);
        if (!Number.isInteger(k) || k < 1)
          throw new Error(`bad k \`${kdeKText}\` (expected an integer >= 1)`);
        const scale = parseFloat(kdeScaleText);
        if (!Number.isFinite(scale) || scale <= 0)
          throw new Error(`bad scale \`${kdeScaleText}\` (expected > 0)`);
        base.k = k;
        base.scale = scale;
      } else if (kdeRule === "fixed") {
        const tokens = kdeWidthsText.split(/\s+/).filter((t) => t.length > 0);
        if (tokens.length === 0) throw new Error(`bad widths \`${kdeWidthsText}\``);
        base.widths = tokens.map((token) => {
          const value = parseFloat(token);
          if (!Number.isFinite(value) || value <= 0)
            throw new Error(`bad width \`${token}\` (expected > 0)`);
          return value;
        });
      }
      const fit = wasm.kdeFit(base);
      const dims = fit.dims;
      const draw = wasm.kdeDraw({ ...base, u: 0.5, normals: Array(dims).fill(0) });
      const meanVec = Array.from({ length: dims }, (_, j) =>
        samples.reduce((sum, row) => sum + row[j], 0),
      ).map((sum) => sum / samples.length);
      const pdfAtMean = wasm.kdePdf({ ...base, point: meanVec });
      // Density curve along the first dimension (other dimensions held at
      // their sample means) plus the fitted-sample histogram.
      const firstValues = samples.map((row) => row[0]);
      const mean = meanVec[0];
      const variance =
        samples.length > 1
          ? firstValues.reduce((sum, v) => sum + (v - mean) ** 2, 0) / (samples.length - 1)
          : 0;
      const std = Math.sqrt(variance);
      const lo = Math.min(...firstValues);
      const hi = Math.max(...firstValues);
      const pad = hi - lo || 1;
      const loX = lo - 0.2 * pad;
      const hiX = hi + 0.2 * pad;
      const steps = 120;
      const xs = Array.from({ length: steps + 1 }, (_, i) => loX + ((hiX - loX) * i) / steps);
      const ys = xs.map((x) =>
        wasm.kdePdf({ ...base, point: meanVec.map((v, j) => (j === 0 ? x : v)) }),
      );
      setKdeResult({ fit, draw, xs, ys, firstValues, mean, std, pdfAtMean });
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setKdeResult(null);
    }
  }

  const displayError = error ?? localError;

  return (
    <div className="rounded-xl border border-border/50 bg-background p-4 space-y-6">
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
          <div className="space-y-2">
            <Label>Mesh tally text</Label>
            <Textarea
              value={meshtalText}
              onChange={(e) => {
                setMeshtalText(e.target.value);
                clearError();
              }}
              className="font-mono text-xs"
            />
          </div>

          <div className="grid gap-3 sm:grid-cols-4">
            <div className="space-y-1">
              <Label>Tally</Label>
              <Input
                type="number"
                value={tallyNumber}
                onChange={(e) => {
                  setTallyNumber(parseInt(e.target.value));
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label>Selection</Label>
              <Select
                value={selection}
                onChange={(v) => {
                  setSelection(v as "total" | "perGroup");
                  clearError();
                }}
                options={SELECTION_OPTIONS}
                className="min-w-[140px]"
              />
            </div>
            <div className="space-y-1">
              <Label>Tolerance</Label>
              <Input
                type="number"
                step="0.1"
                value={tolerance}
                onChange={(e) => {
                  setTolerance(parseFloat(e.target.value));
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label>Null value</Label>
              <Input
                type="number"
                step="any"
                value={nullValue}
                onChange={(e) => {
                  setNullValue(parseFloat(e.target.value));
                  clearError();
                }}
              />
            </div>
          </div>

          <Button onClick={runMagic}>Generate MAGIC bounds</Button>

          {magic && (
            <div className="space-y-3 text-sm">
              <p>Tag: {magic.wwTagName}</p>
              <p>Groups per voxel: {magic.groupsPerVe}</p>
              <p>Lower bounds: [{magic.lowerBoundsWw.map((v) => v.toExponential(3)).join(", ")}]</p>
              <p>Scale factors: [{magic.scaleFactors.map((v) => v.toExponential(3)).join(", ")}]</p>
              <MagicBoundsChart magic={magic} />
            </div>
          )}

          <div className="space-y-2 border-t border-border/50 pt-4">
            <Label>Alias-table sampling</Label>
            <Input
              type="text"
              value={pdfInput}
              onChange={(e) => {
                setPdfInput(e.target.value);
                clearError();
              }}
              className="font-mono"
            />
            <Button onClick={runAlias}>Sample index</Button>
            {aliasResult !== null && <p className="font-mono">Sampled index: {aliasResult}</p>}
          </div>

          <div className="space-y-2 border-t border-border/50 pt-4">
            <Label>Mesh source sampling</Label>
            <div className="flex flex-wrap items-end gap-2">
              <Select
                value={sampleMode}
                onChange={(v) => {
                  setSampleMode(v as "analog" | "uniform");
                  clearError();
                }}
                options={SAMPLE_MODE_OPTIONS}
                className="min-w-[140px]"
              />
              <Button onClick={runSample}>Sample voxel</Button>
            </div>
            {sample && (
              <p className="font-mono">
                index={sample.index} i={sample.i} j={sample.j} k={sample.k} weight=
                {sample.weight.toFixed(4)}
              </p>
            )}
          </div>

          <div className="space-y-2 border-t border-border/50 pt-4">
            <Label>KDE source resampling</Label>
            <p className="text-sm text-muted-foreground">
              Fits a Gaussian KDE over caller-supplied particle rows (synthetic
              {` {0, 0.2, 0.4, 10}`} hand vector by default) and deterministically resamples with
              caller-supplied uniforms and normals — no RNG inside. Silverman gives one global width
              per dimension, fixed takes explicit widths, and adaptive gives one width per sample
              from k-nearest-neighbour distances.
            </p>
            <div className="space-y-1">
              <Label>KDE samples (one row per line)</Label>
              <Textarea
                value={kdeSamplesText}
                onChange={(e) => {
                  setKdeSamplesText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="grid gap-3 sm:grid-cols-4">
              <div className="space-y-1">
                <Label>Bandwidth rule</Label>
                <Select
                  value={kdeRule}
                  onChange={(v) => {
                    setKdeRule(v as KdeRule);
                    clearError();
                  }}
                  options={KDE_RULE_OPTIONS}
                  className="min-w-[140px]"
                />
              </div>
              {kdeRule === "adaptive" && (
                <>
                  <div className="space-y-1">
                    <Label>k (adaptive)</Label>
                    <Input
                      type="number"
                      value={kdeKText}
                      onChange={(e) => {
                        setKdeKText(e.target.value);
                        clearError();
                      }}
                    />
                  </div>
                  <div className="space-y-1">
                    <Label>Scale (adaptive)</Label>
                    <Input
                      type="number"
                      step="any"
                      value={kdeScaleText}
                      onChange={(e) => {
                        setKdeScaleText(e.target.value);
                        clearError();
                      }}
                    />
                  </div>
                </>
              )}
              {kdeRule === "fixed" && (
                <div className="space-y-1">
                  <Label>Widths (fixed)</Label>
                  <Input
                    type="text"
                    value={kdeWidthsText}
                    onChange={(e) => {
                      setKdeWidthsText(e.target.value);
                      clearError();
                    }}
                    className="font-mono"
                  />
                </div>
              )}
            </div>
            <Button onClick={runKde}>Fit KDE</Button>
            {kdeResult && (
              <div className="space-y-3 text-sm">
                <p>
                  KDE bandwidths: [
                  {kdeResult.fit.bandwidths.map((v) => v.toExponential(3)).join(", ")}]
                </p>
                <p className="font-mono">
                  Sample bandwidths:{" "}
                  {kdeResult.fit.sampleBandwidths
                    .map((row) => `[${row.map((v) => v.toExponential(3)).join(", ")}]`)
                    .join(" ")}
                </p>
                <p>
                  KDE moments: mean={kdeResult.mean.toFixed(4)} std=
                  {kdeResult.std.toExponential(4)} pdf(mean)=
                  {kdeResult.pdfAtMean.toExponential(4)}
                </p>
                <p className="font-mono">
                  Draw (u=0.5, normals=0): [{kdeResult.draw.map((v) => v.toFixed(4)).join(", ")}]
                </p>
                <KdeChart result={kdeResult} />
              </div>
            )}
          </div>
        </>
      )}
    </div>
  );
}

function KdeChart({ result }: { result: KdeResult }) {
  return (
    <Plotly
      aspect="video"
      data={[
        {
          type: "histogram",
          name: "samples",
          x: result.firstValues,
          nbinsx: 10,
        },
        {
          type: "scatter",
          mode: "lines",
          name: "KDE pdf",
          x: result.xs,
          y: result.ys,
        },
      ]}
      layout={{
        title: { text: "KDE density (samples + pdf)" },
        xaxis: { title: { text: "Value (dimension 1)" }, type: "linear" },
        yaxis: { title: { text: "Count / density" }, type: "linear" },
        margin: { t: 16, r: 24, b: 48, l: 64 },
        legend: { orientation: "h", y: -0.25 },
      }}
    />
  );
}

function MagicBoundsChart({ magic }: { magic: MagicSummary }) {
  const n = magic.eUpperBounds.length;
  const groups = magic.eUpperBounds.map((_, i) => i);
  // A log axis over a single point (the default demo tally has one energy
  // group) collapses its ticks into an unreadable stack, and a one-point
  // "line" renders nothing. Pad explicitly and mark lone points instead.
  const positives = magic.eUpperBounds.filter((v) => Number.isFinite(v) && v > 0);
  const logAxis = positives.length ? { type: "log" as const, dtick: 1 } : {};
  return (
    <Plotly
      aspect="video"
      data={[
        {
          type: "scatter",
          mode: n === 1 ? "markers" : "lines",
          name: "Energy upper bounds",
          x: groups,
          y: magic.eUpperBounds,
          line: { shape: "hv" },
          marker: { size: 8 },
        },
      ]}
      layout={{
        title: { text: `Energy upper bounds (${n} group${n === 1 ? "" : "s"})` },
        xaxis: {
          title: { text: "Energy group" },
          range: [-0.5, Math.max(n - 0.5, 0.5)],
        },
        yaxis: { title: { text: "Upper bound (MeV)" }, ...logAxis },
        margin: { t: 16, r: 16, b: 48, l: 64 },
      }}
    />
  );
}
