import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { TritiumBreakthroughResult } from "../../types/nucleide-wasm";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";

// Preset traces to synthetic sources only (no library data):
// fixtures/tritium/g2_timelag.json (D = 1e-9 m²/s, L = 1e-3 m,
// c0 = 1.0 mol/m³, t_lag = 166.67 s).
const DEFAULTS = {
  length: "0.001",
  diffusivity: "1e-9",
  c0: "1",
  tEnd: "1000",
  steps: "200",
};

const MAX_STEPS = 500;

export function TritiumBreakthrough() {
  const { wasm, ready, error } = useWasm();
  const [lengthText, setLengthText] = useState(DEFAULTS.length);
  const [diffusivityText, setDiffusivityText] = useState(DEFAULTS.diffusivity);
  const [c0Text, setC0Text] = useState(DEFAULTS.c0);
  const [tEndText, setTEndText] = useState(DEFAULTS.tEnd);
  const [stepsText, setStepsText] = useState(DEFAULTS.steps);
  const [result, setResult] = useState<TritiumBreakthroughResult | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function parseScalar(text: string, label: string): number {
    const v = parseFloat(text);
    if (!Number.isFinite(v)) throw new Error(`bad ${label} \`${text}\``);
    return v;
  }

  function run() {
    if (!wasm) return;
    try {
      const length = parseScalar(lengthText, "length");
      const diffusivity = parseScalar(diffusivityText, "diffusivity");
      const c0 = parseScalar(c0Text, "c0");
      const tEnd = parseScalar(tEndText, "t_end");
      const steps = Math.min(Math.max(Math.floor(parseScalar(stepsText, "steps")), 2), MAX_STEPS);
      if (length <= 0) throw new Error(`bad length \`${lengthText}\` (expected > 0)`);
      if (diffusivity <= 0)
        throw new Error(`bad diffusivity \`${diffusivityText}\` (expected > 0)`);
      if (c0 < 0) throw new Error(`bad c0 \`${c0Text}\` (expected >= 0)`);
      if (tEnd <= 0) throw new Error(`bad t_end \`${tEndText}\` (expected > 0)`);
      const times = Array.from({ length: steps }, (_, k) => (tEnd * (k + 1)) / steps);
      const out = wasm.tritiumBreakthrough(length, diffusivity, c0, times);
      setResult(out);
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setResult(null);
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
          <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-4">
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Slab length L [m]</Label>
              <Input
                value={lengthText}
                onChange={(e) => {
                  setLengthText(e.target.value);
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Diffusivity D [m²/s]</Label>
              <Input
                value={diffusivityText}
                onChange={(e) => {
                  setDiffusivityText(e.target.value);
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Entry concentration c0 [mol/m³]</Label>
              <Input
                value={c0Text}
                onChange={(e) => {
                  setC0Text(e.target.value);
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">End time [s]</Label>
              <Input
                value={tEndText}
                onChange={(e) => {
                  setTEndText(e.target.value);
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Output steps (max {MAX_STEPS})</Label>
              <Input
                value={stepsText}
                onChange={(e) => {
                  setStepsText(e.target.value);
                  clearError();
                }}
              />
            </div>
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Run breakthrough</Button>
          </div>

          {result && (
            <div className="space-y-2">
              <p className="text-sm">
                Time lag t_lag = L²/6D ={" "}
                <span className="font-mono">{result.tLag.toExponential(4)} s</span> (steady flux
                J_ss = <span className="font-mono">{result.jss.toExponential(4)} mol/m²/s</span>)
              </p>
              <p className="text-sm">
                Final J/J_ss({result.times[result.times.length - 1].toExponential(2)} s) ={" "}
                <span className="font-mono">
                  {result.fluxOverJss[result.fluxOverJss.length - 1].toExponential(4)}
                </span>
              </p>
              <Plotly
                aspect="video"
                data={[
                  {
                    type: "scatter",
                    mode: "lines",
                    name: "J/J_ss",
                    x: result.times,
                    y: result.fluxOverJss,
                  },
                ]}
                layout={{
                  xaxis: { title: { text: "Time (s)" }, type: "log" },
                  yaxis: { title: { text: "Normalized outlet flux J/J_ss" }, type: "linear" },
                  margin: { t: 16, r: 24, b: 48, l: 64 },
                  legend: { orientation: "h", y: -0.25 },
                  shapes: [
                    {
                      type: "line",
                      x0: result.tLag,
                      x1: result.tLag,
                      y0: 0,
                      y1: 1,
                      line: { dash: "dash", width: 1 },
                    },
                  ],
                }}
              />
            </div>
          )}
        </>
      )}
    </div>
  );
}
