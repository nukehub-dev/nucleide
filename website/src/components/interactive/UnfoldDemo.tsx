import { useState, type SubmitEvent } from "react";
import { useWasm } from "../../lib/wasm";
import type {
  CorrelatedSolutionJson,
  SandiiSolutionJson,
  WasmApi,
} from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";

// Synthetic unfolding problem: 4 detectors x
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

const CORRELATED_RESPONSE = [
  [1, 0],
  [0, 1],
];
const CORRELATED_GUESS = [1, 1];
const CORRELATED_DAMPING = 1;

function groupLabel(mid: number): string {
  return mid.toPrecision(2);
}

function inputNumber(text: string, label: string): number {
  const value = Number(text);
  if (!text.trim() || !Number.isFinite(value)) throw new Error(`Enter a number for ${label}.`);
  return value;
}

function CorrelatedUnfoldPanel({ wasm }: { wasm: WasmApi }) {
  const [ratesText, setRatesText] = useState(["3", "2"]);
  const [correlationText, setCorrelationText] = useState("0.5");
  const [fit, setFit] = useState<{
    independent: CorrelatedSolutionJson;
    correlated: CorrelatedSolutionJson;
  } | null>(null);
  const [error, setError] = useState<string | null>(null);

  function clearFit() {
    setFit(null);
    setError(null);
  }

  function run(event: SubmitEvent<HTMLFormElement>) {
    event.preventDefault();
    try {
      const rates = ratesText.map((text, i) =>
        inputNumber(text, `detector ${i + 1}'s measured rate`),
      );
      const correlation = inputNumber(correlationText, "the correlation coefficient");
      const correlated = wasm.unfoldStayslCorrelated(
        CORRELATED_RESPONSE,
        rates,
        [
          [1, correlation],
          [correlation, 1],
        ],
        CORRELATED_GUESS,
        CORRELATED_DAMPING,
      );
      const independent = wasm.unfoldStayslCorrelated(
        CORRELATED_RESPONSE,
        rates,
        [
          [1, 0],
          [0, 1],
        ],
        CORRELATED_GUESS,
        CORRELATED_DAMPING,
      );
      setFit({ independent, correlated });
      setError(null);
    } catch (err) {
      setError(err instanceof Error ? err.message : String(err));
      setFit(null);
    }
  }

  return (
    <section
      aria-labelledby="correlated-unfold-title"
      className="space-y-3 border-t border-border/50 pt-4"
    >
      <h3 id="correlated-unfold-title" className="text-base font-semibold">
        Fit correlated measurements
      </h3>
      <p className="text-sm text-muted-foreground">
        Each of the two synthetic detectors responds only to its matching group. The guess is 1 per
        group, both rates have unit variance, and damping is 1.
      </p>
      <form onSubmit={run} className="space-y-3">
        <div className="grid gap-3 sm:grid-cols-3">
          {ratesText.map((text, i) => (
            <div key={i} className="flex flex-col justify-end gap-1">
              <Label htmlFor={`correlated-rate-${i}`}>Detector {i + 1} measured rate</Label>
              <Input
                id={`correlated-rate-${i}`}
                inputMode="decimal"
                value={text}
                onChange={(event) => {
                  setRatesText(ratesText.map((rate, j) => (j === i ? event.target.value : rate)));
                  clearFit();
                }}
                className="font-mono"
              />
            </div>
          ))}
          <div className="flex flex-col justify-end gap-1">
            <Label htmlFor="unfold-correlation">Correlation coefficient</Label>
            <Input
              id="unfold-correlation"
              inputMode="decimal"
              value={correlationText}
              onChange={(event) => {
                setCorrelationText(event.target.value);
                clearFit();
              }}
              aria-describedby="unfold-correlation-range"
              className="font-mono"
            />
          </div>
        </div>
        <p id="unfold-correlation-range" className="text-xs text-muted-foreground">
          Rates use arbitrary units. A coefficient of 0 gives independent errors, and the allowed
          range is strictly between −1 and 1.
        </p>
        <Button type="submit">Fit spectrum</Button>
      </form>

      {error && (
        <p
          role="alert"
          className="break-words rounded-lg border border-red-200 bg-red-50 p-3 text-sm text-red-700 dark:border-red-900 dark:bg-red-950 dark:text-red-300"
        >
          Fit error: {error}
        </p>
      )}
      <div aria-live="polite">
        {fit && (
          <div className="space-y-3">
            <table className="w-full text-left text-sm sm:w-auto">
              <caption className="pb-1 text-left font-medium">Adjusted spectrum</caption>
              <thead className="text-muted-foreground">
                <tr>
                  <th scope="col" className="pr-3 font-medium">
                    Group
                  </th>
                  <th scope="col" className="pr-3 font-medium">
                    Independent errors
                  </th>
                  <th scope="col" className="font-medium">
                    Correlated errors
                  </th>
                </tr>
              </thead>
              <tbody>
                {fit.correlated.spectrum.map((value, i) => (
                  <tr key={i} className="font-mono tabular-nums">
                    <th scope="row" className="pr-3 font-normal">
                      {i + 1}
                    </th>
                    <td className="pr-3">{fit.independent.spectrum[i].toFixed(4)}</td>
                    <td>{value.toFixed(4)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
            <dl className="grid gap-2 rounded-lg bg-muted/30 p-3 text-sm sm:grid-cols-3">
              <div>
                <dt className="text-muted-foreground">Measurement χ²</dt>
                <dd className="font-mono">{fit.correlated.measurementChi2.toPrecision(4)}</dd>
              </div>
              <div>
                <dt className="text-muted-foreground">Damping penalty</dt>
                <dd className="font-mono">{fit.correlated.dampingPenalty.toPrecision(4)}</dd>
              </div>
              <div>
                <dt className="text-muted-foreground">Optimality error</dt>
                <dd className="font-mono">{fit.correlated.optimality.toExponential(2)}</dd>
              </div>
            </dl>
          </div>
        )}
      </div>
    </section>
  );
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
      {wasm && <CorrelatedUnfoldPanel wasm={wasm} />}
    </div>
  );
}
