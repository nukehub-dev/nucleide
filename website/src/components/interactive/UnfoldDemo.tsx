import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { SandiiSolutionJson } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";

// Determined synthetic unfolding problem (no evaluated data): 4 detectors x
// 6 energy groups. Each detector row is a Gaussian bump centered on a
// different group (plus a small floor so every group carries some weight);
// the bump matrix is the crate-test construction. Default measured rates are
// the forward fold of the synthetic truth spectrum [3, 5, 2, 0.5, 0.2, 0.05].
const RESPONSE: number[][] = [
  [1.001, 0.64218039, 0.17001332, 0.019315639, 0.0018159878, 0.0010149453],
  [0.29196046, 0.82175481, 0.95281678, 0.45478877, 0.089943576, 0.008166975],
  [0.008166975, 0.089943576, 0.45478877, 0.95281678, 0.82175481, 0.29196046],
  [0.0010149453, 0.0018159878, 0.019315639, 0.17001332, 0.64218039, 1.001],
];
const GROUP_MIDPOINTS_MEV = [3.831e-6, 5.623e-5, 8.254e-4, 1.212e-2, 0.1778, 2.61];
const GUESS = [1, 1, 1, 1, 1, 1];
const DEFAULT_RATES = ["6.5640", "7.1361", "2.0392", "0.3142"];
const DEFAULT_TOLERANCE = "1e-3";

function groupLabel(mid: number): string {
  return mid.toPrecision(2);
}

export function UnfoldDemo() {
  const { wasm, ready, error } = useWasm();
  const [ratesText, setRatesText] = useState(DEFAULT_RATES);
  const [toleranceText, setToleranceText] = useState(DEFAULT_TOLERANCE);
  const [solution, setSolution] = useState<SandiiSolutionJson | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function run() {
    if (!wasm) return;
    try {
      const rates = ratesText.map((t, i) => {
        const v = parseFloat(t);
        if (!Number.isFinite(v)) throw new Error(`bad measured rate ${i + 1} \`${t}\``);
        if (v < 0) throw new Error(`measured rate ${i + 1} must be >= 0`);
        return v;
      });
      const tolerance = parseFloat(toleranceText);
      if (!Number.isFinite(tolerance) || tolerance <= 0)
        throw new Error(`bad tolerance \`${toleranceText}\` (expected > 0)`);
      const out = wasm.sandiiSolve(RESPONSE, rates, GUESS, tolerance, 500);
      setSolution(out);
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setSolution(null);
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
          <p className="text-sm text-muted-foreground">
            Response matrix (rows = detectors, columns = energy groups; values are response
            strengths on a synthetic Gaussian-bump system, not evaluated data):
          </p>
          <table className="text-sm">
            <thead>
              <tr className="text-left text-muted-foreground">
                <th className="pr-4 font-medium">Detector</th>
                {GROUP_MIDPOINTS_MEV.map((m) => (
                  <th key={m} className="pr-3 font-mono text-xs font-medium">
                    {groupLabel(m)}
                  </th>
                ))}
              </tr>
            </thead>
            <tbody className="font-mono text-xs">
              {RESPONSE.map((row, i) => (
                <tr key={i}>
                  <td className="pr-4">{i + 1}</td>
                  {row.map((v, j) => (
                    <td key={j} className="pr-3">
                      {v.toPrecision(3)}
                    </td>
                  ))}
                </tr>
              ))}
            </tbody>
          </table>
          <p className="text-xs text-muted-foreground">
            Column headers are group midpoints in MeV (log-spaced 10⁻⁶–10 MeV). The guess spectrum
            is one flat value per group; SAND-II adjusts it multiplicatively until the fold
            reproduces the measured rates.
          </p>

          <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-5">
            {DEFAULT_RATES.map((_, i) => (
              <div key={i} className="space-y-1">
                <Label className="flex min-h-10 items-end">Measured rate, detector {i + 1}</Label>
                <Input
                  value={ratesText[i]}
                  onChange={(e) => {
                    setRatesText(ratesText.map((t, j) => (j === i ? e.target.value : t)));
                    clearError();
                  }}
                  className="font-mono text-xs"
                />
              </div>
            ))}
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Tolerance</Label>
              <Input
                value={toleranceText}
                onChange={(e) => {
                  setToleranceText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Unfold spectrum</Button>
          </div>

          {solution && (
            <div className="space-y-3">
              <p className="text-sm">
                Converged in <span className="font-mono">{solution.iterations}</span> adjustments
                (final max per-group relative change{" "}
                <span className="font-mono">{solution.maxRelChange.toExponential(3)}</span>,
                tolerance <span className="font-mono">{solution.tolerance}</span>).
              </p>
              <Plotly
                aspect="video"
                data={[
                  {
                    type: "bar",
                    name: "guess",
                    x: GROUP_MIDPOINTS_MEV.map(groupLabel),
                    y: GUESS,
                  },
                  {
                    type: "bar",
                    name: "unfolded",
                    x: GROUP_MIDPOINTS_MEV.map(groupLabel),
                    y: solution.spectrum,
                  },
                ]}
                layout={{
                  xaxis: { title: { text: "Group midpoint (MeV)" }, type: "category" },
                  yaxis: { title: { text: "Spectrum per group" }, type: "linear" },
                  barmode: "group",
                  margin: { t: 16, r: 24, b: 48, l: 64 },
                  legend: { orientation: "h", y: -0.25 },
                }}
              />
              <div>
                <p className="text-sm font-medium">Measured vs folded rates</p>
                <table className="mt-1 text-sm">
                  <thead>
                    <tr className="text-left text-muted-foreground">
                      <th className="pr-4 font-medium">Detector</th>
                      <th className="pr-4 font-medium">Measured</th>
                      <th className="pr-4 font-medium">Folded</th>
                      <th className="font-medium">N/c ratio</th>
                    </tr>
                  </thead>
                  <tbody className="font-mono text-xs">
                    {solution.rates.map((c, i) => (
                      <tr key={i}>
                        <td className="pr-4">{i + 1}</td>
                        <td className="pr-4">{parseFloat(ratesText[i]).toPrecision(5)}</td>
                        <td className="pr-4">{c.toPrecision(5)}</td>
                        <td>{solution.rateFactors[i].toPrecision(5)}</td>
                      </tr>
                    ))}
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
