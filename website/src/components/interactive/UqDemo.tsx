import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { UqSampleResult } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";

// Defaults trace to a synthetic source only (no evaluated data):
// fixtures/uq/cov_2x2.json (mean [1, 2], variances 0.25/0.16,
// covariance 0.10, seed 20260913).
const DEFAULT_MEAN = "1, 2";
const DEFAULT_COV = "0.25, 0.1, 0.1, 0.16";
const DEFAULT_N = "2000";
const DEFAULT_SEED = "20260913";

const MAX_SAMPLES = 5000;

function fmt(v: number): string {
  return Number.isFinite(v) ? v.toPrecision(5) : String(v);
}

export function UqDemo() {
  const { wasm, ready, error } = useWasm();
  const [meanText, setMeanText] = useState(DEFAULT_MEAN);
  const [covText, setCovText] = useState(DEFAULT_COV);
  const [nText, setNText] = useState(DEFAULT_N);
  const [seedText, setSeedText] = useState(DEFAULT_SEED);
  const [mean, setMean] = useState<number[] | null>(null);
  const [cov, setCov] = useState<number[][] | null>(null);
  const [result, setResult] = useState<UqSampleResult | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function parseArray(text: string, label: string): number[] {
    const values = text
      .split(/[\s,;]+/)
      .filter((s) => s.length > 0)
      .map((s) => parseFloat(s));
    if (values.length === 0) throw new Error(`${label} is empty`);
    for (const v of values) {
      if (!Number.isFinite(v)) throw new Error(`bad ${label} entry \`${text}\``);
    }
    return values;
  }

  function run() {
    if (!wasm) return;
    try {
      const meanValues = parseArray(meanText, "mean");
      const dim = meanValues.length;
      const covValues = parseArray(covText, "cov");
      if (covValues.length !== dim * dim)
        throw new Error(
          `bad cov \`${covText}\` (expected ${dim * dim} row-major entries for dim ${dim})`,
        );
      const covMatrix = Array.from({ length: dim }, (_, i) =>
        covValues.slice(i * dim, (i + 1) * dim),
      );
      const n = parseInt(nText, 10);
      if (!Number.isInteger(n) || n < 2 || n > MAX_SAMPLES)
        throw new Error(`bad n \`${nText}\` (expected an integer in [2, ${MAX_SAMPLES}])`);
      const seed = parseInt(seedText, 10);
      if (!Number.isInteger(seed) || seed < 0)
        throw new Error(`bad seed \`${seedText}\` (expected a non-negative integer)`);
      const out = wasm.uqSample(meanValues, covMatrix, n, seed);
      setMean(meanValues);
      setCov(covMatrix);
      setResult(out);
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setMean(null);
      setCov(null);
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
          <div className="grid gap-3 sm:grid-cols-2">
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Mean (comma-separated)</Label>
              <Input
                value={meanText}
                onChange={(e) => {
                  setMeanText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Covariance (row-major)</Label>
              <Input
                value={covText}
                onChange={(e) => {
                  setCovText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Draws n (max {MAX_SAMPLES})</Label>
              <Input
                value={nText}
                onChange={(e) => {
                  setNText(e.target.value);
                  clearError();
                }}
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
              />
            </div>
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Run sampling</Button>
          </div>

          {result && mean && cov && (
            <div className="space-y-3">
              <p className="text-sm">
                Sampling method: <span className="font-mono">{result.method}</span> (
                {result.samples.length} draws, seed <span className="font-mono">{seedText}</span>)
              </p>
              <div>
                <p className="text-sm font-medium">Sample mean vs input mean</p>
                <table className="mt-1 text-sm">
                  <thead>
                    <tr className="text-left text-muted-foreground">
                      <th className="pr-4 font-medium">Dim</th>
                      <th className="pr-4 font-medium">Input</th>
                      <th className="font-medium">Sample</th>
                    </tr>
                  </thead>
                  <tbody className="font-mono text-xs">
                    {mean.map((m, i) => (
                      <tr key={i}>
                        <td className="pr-4">{i}</td>
                        <td className="pr-4">{fmt(m)}</td>
                        <td>{fmt(result.sampleMean[i])}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              <div>
                <p className="text-sm font-medium">Sample covariance vs input covariance</p>
                <table className="mt-1 text-sm">
                  <thead>
                    <tr className="text-left text-muted-foreground">
                      <th className="pr-4 font-medium">Entry</th>
                      <th className="pr-4 font-medium">Input</th>
                      <th className="font-medium">Sample</th>
                    </tr>
                  </thead>
                  <tbody className="font-mono text-xs">
                    {cov.flatMap((row, i) =>
                      row.map((c, j) => (
                        <tr key={`${i}-${j}`}>
                          <td className="pr-4">
                            [{i},{j}]
                          </td>
                          <td className="pr-4">{fmt(c)}</td>
                          <td>{fmt(result.sampleCov[i][j])}</td>
                        </tr>
                      )),
                    )}
                  </tbody>
                </table>
              </div>
            </div>
          )}
        </>
      )}
    </div>
  );
}
