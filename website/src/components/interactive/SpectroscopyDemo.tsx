import { useId, useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { SpeSummary, SpectroscopySmoothResult, WasmApi } from "../../types/nucleide-wasm";
import { DataTable } from "@nukehub/docs-kit/components/mdx/DataTable";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Select } from "@nukehub/docs-kit/components/ui/Select";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";

// Defaults trace to synthetic sources only (no laboratory data):
// - oracle: fixtures/spectroscopy/smooth_oracle.json counts
//   [2, 5, 1, 6, 3, 8, 4] with peak window c1=2, c2=5 (the worked example
//   in docs/theory/spectroscopy.mdx: gross 10, bg 76/6, net -16/6).
// - peak: hand-picked small numbers (flat ~12 background with one
//   Gaussian-style bump peaking at channel 6).
const ORACLE_COUNTS = "2, 5, 1, 6, 3, 8, 4";
const PEAK_COUNTS = "12, 11, 13, 12, 30, 80, 150, 80, 30, 12, 11, 13, 12";
const EXAMPLE_TSV = "# synthetic pair\n0.662 2.0\n1.170 1.0\n";
const EXAMPLE_DOLLAR_SPE = `$SPEC_ID:
SYNTH
$SPEC_REM:
DET# 7
DETDESC# SYNDET
$DATE_MEA:
02/03/2026 09:15:00
$MEAS_TIM:
90 100
$DATA:
0 3
1
2
4
8
$MCA_CAL:
3
1.5 2.0 0.5 keV
$SHAPE_CAL:
3
0.7 0.0005 0.0000002
`;
const EXAMPLE_PLAIN_SPE = `Spectrum name:  SYNTH-PLAIN
Detector ID:  7
Detector description: SYNDET
Real Time:  100.5
Live Time:  90.25
Acquisition start date:  03-Feb-2026
Acquisition start time:  09:15:00
Starting channel number:  0
Number of channels:  4
Energy Fit:  1.5  2.0  0.5
FWHM Fit:  0.7  0.0005  0.0000002
SPECTRUM

     0:    1.0
     1:    2.0
     2:    4.0
     3:    8.0
`;
const CALIBRATION_CHANNELS = [0, 1, 2];

const METHOD_OPTIONS = [
  { value: "rect3", label: "Rectangular m=3 (E1)" },
  { value: "rect5", label: "Rectangular m=5 (E1)" },
  { value: "rect7", label: "Rectangular m=7 (E1)" },
  { value: "five-point", label: "Five-point (E2)" },
];

type Method = "rect3" | "rect5" | "rect7" | "five-point";

const ACTIVITY_FIELDS = [
  { name: "netCounts", label: "Net peak counts", unit: "counts" },
  { name: "countingUncertainty", label: "Counting uncertainty", unit: "counts" },
  { name: "liveTime", label: "Live time", unit: "s" },
  { name: "peakEfficiency", label: "Peak efficiency", unit: "fraction" },
  { name: "emissionProbability", label: "Emission probability", unit: "fraction" },
] as const;

type ActivityField = (typeof ACTIVITY_FIELDS)[number]["name"];
type EstimateActivity = WasmApi["spectroscopyGammaLineActivity"];

function GammaLineActivityForm({ estimate }: { estimate: EstimateActivity }) {
  const id = useId();
  const [values, setValues] = useState<Record<ActivityField, string>>({
    netCounts: "120",
    countingUncertainty: "12",
    liveTime: "60",
    peakEfficiency: "0.1",
    emissionProbability: "0.5",
  });
  const [result, setResult] = useState<ReturnType<EstimateActivity> | null>(null);
  const [error, setError] = useState<string | null>(null);

  function run() {
    setResult(null);
    setError(null);
    try {
      setResult(
        estimate(
          Number(values.netCounts),
          Number(values.countingUncertainty),
          Number(values.liveTime),
          Number(values.peakEfficiency),
          Number(values.emissionProbability),
        ),
      );
    } catch (cause) {
      setError(cause instanceof Error ? cause.message : String(cause));
    }
  }

  return (
    <form
      aria-labelledby={`${id}-heading`}
      className="space-y-4 border-t border-border/50 pt-4"
      onSubmit={(event) => {
        event.preventDefault();
        run();
      }}
    >
      <div className="space-y-1">
        <h3 id={`${id}-heading`} className="text-sm font-medium">
          Gamma-line activity
        </h3>
        <p className="text-sm text-muted-foreground">
          Use a background-subtracted peak area and its standard uncertainty. Enter the full-energy
          peak efficiency and emission probability as fractions greater than 0 and at most 1.
        </p>
      </div>
      <div className="grid gap-3 sm:grid-cols-2">
        {ACTIVITY_FIELDS.map(({ name, label, unit }) => (
          <div key={name} className="grid content-end gap-1">
            <Label htmlFor={`${id}-${name}`}>{label}</Label>
            <div className="flex items-center gap-2">
              <Input
                id={`${id}-${name}`}
                type="number"
                step="any"
                required
                value={values[name]}
                onChange={(event) => {
                  setValues({ ...values, [name]: event.target.value });
                  setResult(null);
                  setError(null);
                }}
                aria-describedby={`${id}-${name}-unit`}
                className="min-w-0"
              />
              <span
                id={`${id}-${name}-unit`}
                className="w-14 shrink-0 text-xs text-muted-foreground"
              >
                {unit}
              </span>
            </div>
          </div>
        ))}
      </div>
      <Button type="submit">Estimate activity</Button>
      {error && (
        <p
          role="alert"
          className="rounded-lg border border-red-200 bg-red-50 p-3 text-sm text-red-700 dark:border-red-900 dark:bg-red-950 dark:text-red-300"
        >
          {error}
        </p>
      )}
      {result && (
        <output className="flex flex-wrap items-baseline gap-x-3 gap-y-1 rounded-lg bg-muted p-3">
          <span className="text-sm text-muted-foreground">Estimated activity</span>
          <span className="break-all font-mono text-lg">
            {Number(result.activityBq.toPrecision(8))} ±{" "}
            {Number(result.countingUncertaintyBq.toPrecision(8))} Bq
          </span>
        </output>
      )}
    </form>
  );
}

function SpectroscopyExamples({ wasm }: { wasm: WasmApi }) {
  const [lines, setLines] = useState<ReturnType<WasmApi["parseLinesTsv"]> | null>(null);
  const [calibration, setCalibration] = useState<{ energies: number[]; efficiency: number } | null>(
    null,
  );
  const [spe, setSpe] = useState<{
    format: "Dollar SPE" | "Plain SPE";
    summary: SpeSummary;
  } | null>(null);
  const [error, setError] = useState<string | null>(null);

  function run(action: () => void) {
    setError(null);
    try {
      action();
    } catch (cause) {
      setError(cause instanceof Error ? cause.message : String(cause));
    }
  }

  return (
    <section className="space-y-4 border-t border-border/50 pt-4">
      <div className="space-y-1">
        <h3 className="text-sm font-medium">Decay-line TSV, calibration, and SPE readers</h3>
        <p className="text-sm text-muted-foreground">
          Run the parsers or calibration on built-in sample data.
        </p>
      </div>
      <div className="flex flex-wrap gap-2">
        <Button
          variant="outline"
          onClick={() => run(() => setLines(wasm.parseLinesTsv(EXAMPLE_TSV)))}
        >
          Parse TSV
        </Button>
        <Button
          variant="outline"
          onClick={() =>
            run(() =>
              setCalibration({
                energies: wasm.energyBins(CALIBRATION_CHANNELS, [1.5, 2.0, 0.5]),
                efficiency: wasm.detectorEfficiency(1, [-2.81861504261204, -0.727352820018942], 1),
              }),
            )
          }
        >
          Calibrate
        </Button>
        <Button
          variant="outline"
          onClick={() =>
            run(() =>
              setSpe({ format: "Dollar SPE", summary: wasm.parseDollarSpe(EXAMPLE_DOLLAR_SPE) }),
            )
          }
        >
          Parse dollar SPE
        </Button>
        <Button
          variant="outline"
          onClick={() =>
            run(() =>
              setSpe({ format: "Plain SPE", summary: wasm.parsePlainSpe(EXAMPLE_PLAIN_SPE) }),
            )
          }
        >
          Parse plain SPE
        </Button>
      </div>
      {error && (
        <p role="alert" className="text-sm text-red-700 dark:text-red-300">
          {error}
        </p>
      )}
      <div aria-live="polite" className="grid gap-4 empty:sr-only lg:grid-cols-2">
        {lines && (
          <div className="min-w-0 space-y-2">
            <h4 className="text-sm font-medium">Decay lines</h4>
            <DataTable
              columns={[
                { key: "energy", header: "Energy (MeV)", align: "right" },
                { key: "intensity", header: "Relative intensity", align: "right" },
              ]}
              data={lines.map(([energy, intensity]) => ({ energy: energy.toFixed(3), intensity }))}
              className="tabular-nums"
            />
          </div>
        )}
        {calibration && (
          <div className="min-w-0 space-y-2">
            <h4 className="text-sm font-medium">Energy calibration</h4>
            <DataTable
              columns={[
                { key: "channel", header: "Channel" },
                { key: "energy", header: "Calibrated energy (keV)", align: "right" },
              ]}
              data={calibration.energies.map((energy, i) => ({
                channel: CALIBRATION_CHANNELS[i],
                energy: energy.toFixed(1),
              }))}
              className="tabular-nums"
              caption={
                <>
                  Detector efficiency at 1 MeV{" "}
                  <strong className="text-foreground">
                    {(calibration.efficiency * 100).toFixed(3)}%
                  </strong>
                </>
              }
            />
          </div>
        )}
        {spe && (
          <div className="space-y-2 rounded-xl border border-border/50 bg-card p-4 lg:col-span-2">
            <div className="flex flex-wrap items-baseline justify-between gap-x-3 gap-y-1">
              <h4 className="text-sm font-medium">{spe.format} summary</h4>
              <span className="break-all font-mono text-sm">{spe.summary.spec_name}</span>
            </div>
            <dl className="grid grid-cols-2 gap-3 text-sm sm:grid-cols-4">
              {[
                ["Detector ID", spe.summary.detId],
                ["Channels", spe.summary.channels],
                ["Live time", `${spe.summary.liveTime} s`],
                ["Real time", `${spe.summary.realTime} s`],
              ].map(([label, value]) => (
                <div key={label}>
                  <dt className="text-xs text-muted-foreground">{label}</dt>
                  <dd className="font-mono tabular-nums">{value}</dd>
                </div>
              ))}
            </dl>
          </div>
        )}
      </div>
    </section>
  );
}

export function SpectroscopyDemo() {
  const { wasm, ready, error } = useWasm();
  const [countsText, setCountsText] = useState(ORACLE_COUNTS);
  const [method, setMethod] = useState<Method>("rect5");
  const [c1Text, setC1Text] = useState("2");
  const [c2Text, setC2Text] = useState("5");
  const [counts, setCounts] = useState<number[] | null>(null);
  const [result, setResult] = useState<SpectroscopySmoothResult | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function loadOracleCounts() {
    setCountsText(ORACLE_COUNTS);
    setC1Text("2");
    setC2Text("5");
    setResult(null);
    setCounts(null);
    clearError();
  }

  function loadSyntheticPeak() {
    setCountsText(PEAK_COUNTS);
    setC1Text("4");
    setC2Text("8");
    setResult(null);
    setCounts(null);
    clearError();
  }

  function parseCounts(): number[] {
    const values = countsText
      .split(/[\s,]+/)
      .filter((s) => s.length > 0)
      .map((s) => parseFloat(s));
    if (values.length === 0) throw new Error("counts are empty");
    for (const v of values) {
      if (!Number.isFinite(v) || v < 0)
        throw new Error(`bad count \`${v}\` (expected finite counts >= 0)`);
    }
    return values;
  }

  function parseChannel(text: string, label: string): number {
    const v = parseInt(text, 10);
    if (!Number.isInteger(v)) throw new Error(`bad ${label} \`${text}\` (expected integer)`);
    return v;
  }

  function run() {
    if (!wasm) return;
    try {
      const parsed = parseCounts();
      const c1 = parseChannel(c1Text, "c1");
      const c2 = parseChannel(c2Text, "c2");
      const out = wasm.spectroscopySmooth(parsed, method, c1, c2);
      setCounts(parsed);
      setResult(out);
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setResult(null);
      setCounts(null);
    }
  }

  const displayError = error ?? localError;
  const channels = counts ? counts.map((_, i) => i) : [];

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
          <div className="flex flex-wrap gap-2">
            <Button variant="outline" size="sm" onClick={loadOracleCounts}>
              Oracle counts
            </Button>
            <Button variant="outline" size="sm" onClick={loadSyntheticPeak}>
              Synthetic peak
            </Button>
          </div>

          <div className="space-y-2">
            <Label>Counts per channel (comma/space separated)</Label>
            <Textarea
              value={countsText}
              onChange={(e) => {
                setCountsText(e.target.value);
                clearError();
              }}
              className="font-mono text-sm"
            />
          </div>

          <div className="grid gap-3 sm:grid-cols-3">
            <div className="space-y-1">
              <Label>Smoother</Label>
              <Select
                value={method}
                onChange={(v) => {
                  setMethod(v as Method);
                  clearError();
                }}
                options={METHOD_OPTIONS}
                className="min-w-[180px]"
              />
            </div>
            <div className="space-y-1">
              <Label>Peak start c1 (inclusive)</Label>
              <Input
                value={c1Text}
                onChange={(e) => {
                  setC1Text(e.target.value);
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label>Peak stop c2 (exclusive)</Label>
              <Input
                value={c2Text}
                onChange={(e) => {
                  setC2Text(e.target.value);
                  clearError();
                }}
              />
            </div>
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Smooth spectrum</Button>
          </div>

          {result && (
            <div className="space-y-2">
              <table className="w-full text-sm">
                <thead>
                  <tr className="border-b border-border">
                    <th className="py-1 text-left">Quantity</th>
                    <th className="py-1 text-right">Counts</th>
                  </tr>
                </thead>
                <tbody>
                  <tr className="border-b border-border/50">
                    <td className="py-1 font-mono">Gross counts</td>
                    <td className="py-1 text-right font-mono">{result.gross.toFixed(4)}</td>
                  </tr>
                  <tr className="border-b border-border/50">
                    <td className="py-1 font-mono">Background</td>
                    <td className="py-1 text-right font-mono">{result.background.toFixed(4)}</td>
                  </tr>
                  <tr className="border-b border-border/50">
                    <td className="py-1 font-mono">Net counts</td>
                    <td className="py-1 text-right font-mono">{result.net.toFixed(4)}</td>
                  </tr>
                </tbody>
              </table>
              <Plotly
                aspect="video"
                data={[
                  {
                    type: "scatter",
                    mode: "lines+markers",
                    name: "Raw counts",
                    x: channels,
                    y: counts ?? [],
                  },
                  {
                    type: "scatter",
                    mode: "lines",
                    name: "Smoothed",
                    x: channels,
                    y: result.smoothed,
                  },
                ]}
                layout={{
                  xaxis: { title: { text: "Channel" }, type: "linear" },
                  yaxis: { title: { text: "Counts" }, type: "linear" },
                  margin: { t: 16, r: 24, b: 48, l: 64 },
                  legend: { orientation: "h", y: -0.25 },
                }}
              />
            </div>
          )}
          {wasm && (
            <>
              <GammaLineActivityForm estimate={wasm.spectroscopyGammaLineActivity} />
              <SpectroscopyExamples wasm={wasm} />
            </>
          )}
        </>
      )}
    </div>
  );
}
