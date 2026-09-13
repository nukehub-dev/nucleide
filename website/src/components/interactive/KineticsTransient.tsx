import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { KineticsTransientResult } from "../../types/nucleide-wasm";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";

// Presets trace to synthetic sources only (no library data):
// - one-group: fixtures/kinetics/step_oracle.json (beta=0.0065,
//   lambda=0.08 1/s, Lambda=1e-4 s, step 0 -> 0.002 at t=1 s, n0=1).
// - six-group: fixtures/kinetics/inhour_check.json (hand-picked round
//   numbers, betas sum to 0.0065, Lambda=1e-5 s) with the same step shape.
const ONE_GROUP = {
  betas: "0.0065",
  lambdas: "0.08",
  lambdaGen: "0.0001",
  tStep: "1",
  rhoInit: "0",
  rhoFinal: "0.002",
  tEnd: "10",
  steps: "200",
  n0: "1",
};

const SIX_GROUP = {
  betas: "0.00021, 0.00141, 0.00127, 0.00255, 0.00074, 0.00032",
  lambdas: "0.01, 0.03, 0.1, 0.3, 1.0, 3.0",
  lambdaGen: "0.00001",
  tStep: "1",
  rhoInit: "0",
  rhoFinal: "0.002",
  tEnd: "100",
  steps: "200",
  n0: "1",
};

type PresetName = "one-group" | "six-group";

const PRESETS: Record<PresetName, typeof ONE_GROUP> = {
  "one-group": ONE_GROUP,
  "six-group": SIX_GROUP,
};

const MAX_STEPS = 500;

export function KineticsTransient() {
  const { wasm, ready, error } = useWasm();
  const [preset, setPreset] = useState<PresetName>("one-group");
  const [betasText, setBetasText] = useState(ONE_GROUP.betas);
  const [lambdasText, setLambdasText] = useState(ONE_GROUP.lambdas);
  const [lambdaGenText, setLambdaGenText] = useState(ONE_GROUP.lambdaGen);
  const [tStepText, setTStepText] = useState(ONE_GROUP.tStep);
  const [rhoInitText, setRhoInitText] = useState(ONE_GROUP.rhoInit);
  const [rhoFinalText, setRhoFinalText] = useState(ONE_GROUP.rhoFinal);
  const [tEndText, setTEndText] = useState(ONE_GROUP.tEnd);
  const [stepsText, setStepsText] = useState(ONE_GROUP.steps);
  const [n0Text, setN0Text] = useState(ONE_GROUP.n0);
  const [result, setResult] = useState<KineticsTransientResult | null>(null);
  const [inhourOut, setInhourOut] = useState<string | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function applyPreset(name: PresetName) {
    const p = PRESETS[name];
    setPreset(name);
    setBetasText(p.betas);
    setLambdasText(p.lambdas);
    setLambdaGenText(p.lambdaGen);
    setTStepText(p.tStep);
    setRhoInitText(p.rhoInit);
    setRhoFinalText(p.rhoFinal);
    setTEndText(p.tEnd);
    setStepsText(p.steps);
    setN0Text(p.n0);
    setResult(null);
    clearError();
  }

  function parseArray(text: string, label: string): number[] {
    const values = text
      .split(/[\s,]+/)
      .filter((s) => s.length > 0)
      .map((s) => parseFloat(s));
    if (values.length === 0) throw new Error(`${label} is empty`);
    for (const v of values) {
      if (!Number.isFinite(v)) throw new Error(`bad ${label} entry \`${text}\``);
    }
    return values;
  }

  function parseScalar(text: string, label: string): number {
    const v = parseFloat(text);
    if (!Number.isFinite(v)) throw new Error(`bad ${label} \`${text}\``);
    return v;
  }

  function run() {
    if (!wasm) return;
    try {
      const betas = parseArray(betasText, "betas");
      const lambdas = parseArray(lambdasText, "lambdas");
      const lambdaGen = parseScalar(lambdaGenText, "Lambda");
      const tStep = parseScalar(tStepText, "t_step");
      const rhoInit = parseScalar(rhoInitText, "rho_init");
      const rhoFinal = parseScalar(rhoFinalText, "rho_final");
      const tEnd = parseScalar(tEndText, "t_end");
      const n0 = parseScalar(n0Text, "n0");
      const steps = Math.min(Math.max(Math.floor(parseScalar(stepsText, "steps")), 2), MAX_STEPS);
      if (tEnd <= 0) throw new Error(`bad t_end \`${tEndText}\` (expected > 0)`);
      if (tStep < 0 || tStep >= tEnd)
        throw new Error(`bad t_step \`${tStepText}\` (expected 0 <= t_step < t_end)`);
      if (n0 < 0) throw new Error(`bad n0 \`${n0Text}\` (expected >= 0)`);
      const times = Array.from({ length: steps }, (_, k) => (tEnd * (k + 1)) / steps);
      const out = wasm.kineticsTransient(
        betas,
        lambdas,
        lambdaGen,
        tStep,
        rhoInit,
        rhoFinal,
        times,
        n0,
      );
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
            <Button
              variant={preset === "one-group" ? "default" : "outline"}
              size="sm"
              onClick={() => applyPreset("one-group")}
            >
              One-group preset
            </Button>
            <Button
              variant={preset === "six-group" ? "default" : "outline"}
              size="sm"
              onClick={() => applyPreset("six-group")}
            >
              Six-group preset
            </Button>
          </div>

          <div className="grid gap-3 sm:grid-cols-2">
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Betas (delayed fractions)</Label>
              <Input
                value={betasText}
                onChange={(e) => {
                  setBetasText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Lambdas (decay constants [1/s])</Label>
              <Input
                value={lambdasText}
                onChange={(e) => {
                  setLambdasText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
          </div>

          <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-4">
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Lambda (gen. time [s])</Label>
              <Input
                value={lambdaGenText}
                onChange={(e) => {
                  setLambdaGenText(e.target.value);
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Step time [s]</Label>
              <Input
                value={tStepText}
                onChange={(e) => {
                  setTStepText(e.target.value);
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Rho init [Δk]</Label>
              <Input
                value={rhoInitText}
                onChange={(e) => {
                  setRhoInitText(e.target.value);
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Rho final [Δk]</Label>
              <Input
                value={rhoFinalText}
                onChange={(e) => {
                  setRhoFinalText(e.target.value);
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
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Initial n0</Label>
              <Input
                value={n0Text}
                onChange={(e) => {
                  setN0Text(e.target.value);
                  clearError();
                }}
              />
            </div>
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Run transient</Button>
          </div>

          <div className="space-y-2 border-t border-border/50 pt-4">
            <p className="text-sm font-medium">Inhour and prompt jump (E3–E4 scalars)</p>
            <div className="flex flex-wrap gap-2">
              <Button
                onClick={() => {
                  if (!wasm) return;
                  try {
                    const betas = betasText
                      .split(/[\s,]+/)
                      .filter(Boolean)
                      .map(Number);
                    const lambdas = lambdasText
                      .split(/[\s,]+/)
                      .filter(Boolean)
                      .map(Number);
                    const rho = wasm.inhourRho(betas, lambdas, parseFloat(lambdaGenText), 0.05);
                    const period = wasm.stablePeriod(
                      betas,
                      lambdas,
                      parseFloat(lambdaGenText),
                      0.002,
                    );
                    const jump = wasm.promptJump(
                      1.0,
                      0.0,
                      0.002,
                      betas.reduce((a, b) => a + b, 0),
                    );
                    setInhourOut(
                      `rho(0.05) = ${rho.toExponential(4)}, stable period = ${period.toExponential(4)} s, prompt jump = ${jump.toFixed(4)}`,
                    );
                    clearError();
                  } catch (e) {
                    setLocalError(e instanceof Error ? e.message : String(e));
                    setInhourOut(null);
                  }
                }}
              >
                Inhour analysis
              </Button>
            </div>
            {inhourOut && <p className="text-sm">{inhourOut}</p>}
          </div>

          {result && (
            <div className="space-y-2">
              <p className="text-sm">
                Prompt jump n⁺ ={" "}
                <span className="font-mono">
                  {result.promptJump === null
                    ? "— (ρ⁺ ≥ β, no finite jump)"
                    : result.promptJump.toExponential(4)}
                </span>{" "}
                (β total = <span className="font-mono">{result.betaTotal.toExponential(4)}</span>)
              </p>
              <p className="text-sm">
                Final n({result.times[result.times.length - 1].toExponential(2)} s) ={" "}
                <span className="font-mono">{result.n[result.n.length - 1].toExponential(4)}</span>
              </p>
              <Plotly
                aspect="video"
                data={[
                  {
                    type: "scatter",
                    mode: "lines",
                    name: "n(t)",
                    x: result.times,
                    y: result.n,
                  },
                ]}
                layout={{
                  xaxis: { title: { text: "Time (s)" }, type: "linear" },
                  yaxis: { title: { text: "Neutron population n(t)" }, type: "linear" },
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
