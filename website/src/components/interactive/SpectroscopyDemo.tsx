import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { SpeSummary, SpectroscopySmoothResult } from "../../types/nucleide-wasm";
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

const METHOD_OPTIONS = [
  { value: "rect3", label: "Rectangular m=3 (E1)" },
  { value: "rect5", label: "Rectangular m=5 (E1)" },
  { value: "rect7", label: "Rectangular m=7 (E1)" },
  { value: "five-point", label: "Five-point (E2)" },
];

type Method = "rect3" | "rect5" | "rect7" | "five-point";

export function SpectroscopyDemo() {
  const { wasm, ready, error } = useWasm();
  const [countsText, setCountsText] = useState(ORACLE_COUNTS);
  const [method, setMethod] = useState<Method>("rect5");
  const [c1Text, setC1Text] = useState("2");
  const [c2Text, setC2Text] = useState("5");
  const [counts, setCounts] = useState<number[] | null>(null);
  const [result, setResult] = useState<SpectroscopySmoothResult | null>(null);
  const [tsvOut, setTsvOut] = useState<string | null>(null);
  const [calibOut, setCalibOut] = useState<string | null>(null);
  const [spe, setSpe] = useState<SpeSummary | null>(null);
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
      {!ready && <p className="text-sm text-muted-foreground">Loading Nucleide WASM…</p>}
      {displayError && (
        <div className="flex items-start justify-between gap-2 rounded-lg border border-red-200 bg-red-50 p-3 text-sm text-red-700 dark:border-red-900 dark:bg-red-950 dark:text-red-300">
          <span>WASM error: {displayError}</span>
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

          <div className="space-y-2 border-t border-border/50 pt-4">
            <p className="text-sm font-medium">Decay-line TSV, calibration, and SPE readers</p>
            <div className="flex flex-wrap gap-2">
              <Button
                onClick={() => {
                  if (!wasm) return;
                  try {
                    const rows = wasm.parseLinesTsv("# synthetic pair\n0.662 2.0\n1.170 1.0\n");
                    setTsvOut(
                      `TSV rows: ${rows.length}, first energy ${(rows[0][0] as number).toFixed(3)} MeV`,
                    );
                    clearError();
                  } catch (e) {
                    setLocalError(e instanceof Error ? e.message : String(e));
                  }
                }}
              >
                Parse TSV
              </Button>
              <Button
                onClick={() => {
                  if (!wasm) return;
                  try {
                    const ebins = wasm.energyBins([0, 1, 2], [1.5, 2.0, 0.5]) as number[];
                    const eff = wasm.detectorEfficiency(
                      1.0,
                      [-2.81861504261204, -0.727352820018942],
                      1,
                    );
                    setCalibOut(
                      `ebins: ${ebins.map((v) => v.toFixed(1)).join(", ")}, eff(1 MeV) = ${eff.toExponential(3)}`,
                    );
                    clearError();
                  } catch (e) {
                    setLocalError(e instanceof Error ? e.message : String(e));
                  }
                }}
              >
                Calibrate
              </Button>
              <Button
                onClick={() => {
                  if (!wasm) return;
                  try {
                    setSpe(
                      wasm.parseDollarSpe(
                        "$SPEC_ID:\nSYNTH\n$SPEC_REM:\nDET# 7\nDETDESC# SYNDET\n$DATE_MEA:\n02/03/2026 09:15:00\n$MEAS_TIM:\n90 100\n$DATA:\n0 3\n1\n2\n4\n8\n$MCA_CAL:\n3\n1.5 2.0 0.5 keV\n$SHAPE_CAL:\n3\n0.7 0.0005 0.0000002\n",
                      ),
                    );
                    clearError();
                  } catch (e) {
                    setLocalError(e instanceof Error ? e.message : String(e));
                  }
                }}
              >
                Parse dollar SPE
              </Button>
              <Button
                onClick={() => {
                  if (!wasm) return;
                  try {
                    setSpe(
                      wasm.parsePlainSpe(
                        "Spectrum name:  SYNTH-PLAIN\nDetector ID:  7\nDetector description: SYNDET\nReal Time:  100.5\nLive Time:  90.25\nAcquisition start date:  03-Feb-2026\nAcquisition start time:  09:15:00\nStarting channel number:  0\nNumber of channels:  4\nEnergy Fit:  1.5  2.0  0.5\nFWHM Fit:  0.7  0.0005  0.0000002\nSPECTRUM\n\n     0:    1.0\n     1:    2.0\n     2:    4.0\n     3:    8.0\n",
                      ),
                    );
                    clearError();
                  } catch (e) {
                    setLocalError(e instanceof Error ? e.message : String(e));
                  }
                }}
              >
                Parse plain SPE
              </Button>
            </div>
            {tsvOut && <p className="text-sm">{tsvOut}</p>}
            {calibOut && <p className="text-sm">{calibOut}</p>}
            {spe && (
              <p className="text-sm">
                SPE {spe.spec_name}: channels {spe.channels}, live {spe.liveTime} s, real{" "}
                {spe.realTime} s
              </p>
            )}
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
        </>
      )}
    </div>
  );
}
