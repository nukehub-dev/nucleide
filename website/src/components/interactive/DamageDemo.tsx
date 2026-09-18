import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { CoilFlux, CoilLifetime } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";

// Three-group synthetic grid (MeV) with hand-picked round flux/response
// values — the same shape as the in-crate fold tests (no evaluated data):
// flux is the per-group integrated flux [n/cm^2/s], responses are group
// cross sections [barns] (dpa damage and He gas production separately).
const GROUP_BOUNDS = [0.0, 0.1, 1.0, 20.0];
const DEFAULT_FLUX = ["1e13", "1e12", "1e11"];
const DEFAULT_DPA_RESPONSE = ["50", "100", "200"];
const DEFAULT_HE_RESPONSE = ["0.2", "0.4", "0.6"];
const DEFAULT_SECONDS = "3";
const DEFAULT_TARGET = "Fe56";
const DEFAULT_ED = "40";
const DEFAULT_B_ARC = "-0.55";
const DEFAULT_C_ARC = "0.3";
const DEFAULT_THRESHOLD = "0.1";
const DEFAULT_LIMITS = "1e20, 5e20";
const DEFAULT_RATES = "1e12, 1e12";

const GROUP_LABELS = GROUP_BOUNDS.slice(0, 3).map((lo, i) => `${lo}-${GROUP_BOUNDS[i + 1]} MeV`);

interface DamageMetrics {
  nrtDpa: number;
  arcDpa: number;
  heAppm: number;
  heDpaRatio: number | null;
}

interface ArcCurve {
  tEv: number[];
  xi: number[];
}

function fmt(v: number): string {
  return Number.isFinite(v) ? v.toPrecision(5) : String(v);
}

export function DamageDemo() {
  const { wasm, ready, error } = useWasm();
  const [fluxText, setFluxText] = useState(DEFAULT_FLUX);
  const [dpaText, setDpaText] = useState(DEFAULT_DPA_RESPONSE);
  const [heText, setHeText] = useState(DEFAULT_HE_RESPONSE);
  const [secondsText, setSecondsText] = useState(DEFAULT_SECONDS);
  const [targetText, setTargetText] = useState(DEFAULT_TARGET);
  const [edText, setEdText] = useState(DEFAULT_ED);
  const [bArcText, setBArcText] = useState(DEFAULT_B_ARC);
  const [cArcText, setCArcText] = useState(DEFAULT_C_ARC);
  const [metrics, setMetrics] = useState<DamageMetrics | null>(null);
  const [curve, setCurve] = useState<ArcCurve | null>(null);
  const [thresholdText, setThresholdText] = useState(DEFAULT_THRESHOLD);
  const [limitsText, setLimitsText] = useState(DEFAULT_LIMITS);
  const [ratesText, setRatesText] = useState(DEFAULT_RATES);
  const [coilFlux, setCoilFlux] = useState<CoilFlux | null>(null);
  const [coilLife, setCoilLife] = useState<CoilLifetime | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function parseEntry(text: string, label: string): number {
    const v = parseFloat(text);
    if (!Number.isFinite(v)) throw new Error(`bad ${label} \`${text}\``);
    return v;
  }

  function run() {
    if (!wasm) return;
    try {
      const flux = fluxText.map((t, i) => parseEntry(t, `flux group ${i + 1}`));
      const dpaResponse = dpaText.map((t, i) => parseEntry(t, `dpa response group ${i + 1}`));
      const heResponse = heText.map((t, i) => parseEntry(t, `He response group ${i + 1}`));
      for (const [values, label] of [
        [flux, "flux"],
        [dpaResponse, "dpa response"],
        [heResponse, "He response"],
      ] as const) {
        if (values.some((v) => v < 0)) throw new Error(`${label} entries must be >= 0`);
      }
      const seconds = parseEntry(secondsText, "seconds");
      if (seconds <= 0) throw new Error(`bad seconds \`${secondsText}\` (expected > 0)`);
      const target = targetText.trim();
      const edEv = parseEntry(edText, "E_d");
      if (edEv <= 0) throw new Error(`bad E_d \`${edText}\` (expected > 0)`);
      const bArc = parseEntry(bArcText, "b_arc");
      const cArc = parseEntry(cArcText, "c_arc");

      const nrtDpa = wasm.damageNrtDpa(flux, dpaResponse, GROUP_BOUNDS, seconds);
      const arcDpa = wasm.damageArcDpa(flux, dpaResponse, GROUP_BOUNDS, seconds);
      const heAppm = wasm.damageGasAppm(flux, heResponse, GROUP_BOUNDS, seconds);
      let heDpaRatio: number | null;
      try {
        heDpaRatio = wasm.damageHeDpaRatio(flux, heResponse, dpaResponse, GROUP_BOUNDS, seconds);
      } catch {
        heDpaRatio = null; // exactly zero dpa: reported as "—", not an error
      }

      // arc efficiency curve xi(T_dam(T)) over a PKA-energy sweep (log axis).
      const tEv: number[] = [];
      const xi: number[] = [];
      for (let k = 0; k < 60; k++) {
        const t = 10 ** (1 + (k * 6) / 59); // 10 eV .. 10^7 eV
        const tDam = wasm.damageEnergy(t, target);
        tEv.push(t);
        xi.push(wasm.arcEfficiency(tDam, edEv, bArc, cArc));
      }

      setMetrics({ nrtDpa, arcDpa, heAppm, heDpaRatio });
      setCurve({ tEv, xi });

      // Coil fast-flux / weakest-link lifetime over the same group grid.
      const thresholdMeV = parseEntry(thresholdText, "fast threshold");
      const limits = limitsText.split(",").map((t, i) => parseEntry(t.trim(), `limit ${i + 1}`));
      const rates = ratesText.split(",").map((t, i) => parseEntry(t.trim(), `rate ${i + 1}`));
      setCoilFlux(wasm.coilFastFlux(flux, GROUP_BOUNDS, thresholdMeV, seconds));
      setCoilLife(wasm.coilLifetime(limits, rates));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setMetrics(null);
      setCurve(null);
      setCoilFlux(null);
      setCoilLife(null);
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
            Three energy groups on a fixed grid ({GROUP_LABELS.join(", ")}). Flux is the per-group
            integrated flux [n/cm²/s]; responses are group cross sections [barns] condensed
            consistently with the flux spectrum. The fold is piecewise-constant per group: metric =
            seconds · Σ_g flux[g]·response[g] (with the barns→cm² and per-atom→appm scale factors).
          </p>

          <div className="grid gap-3 sm:grid-cols-3">
            {GROUP_LABELS.map((label, i) => (
              <div key={label} className="space-y-1">
                <Label className="flex min-h-10 items-end">Flux, {label} [n/cm²/s]</Label>
                <Input
                  value={fluxText[i]}
                  onChange={(e) => {
                    setFluxText(fluxText.map((t, j) => (j === i ? e.target.value : t)));
                    clearError();
                  }}
                  className="font-mono text-xs"
                />
              </div>
            ))}
          </div>
          <div className="grid gap-3 sm:grid-cols-3">
            {GROUP_LABELS.map((label, i) => (
              <div key={label} className="space-y-1">
                <Label className="flex min-h-10 items-end">dpa response, {label} [barns]</Label>
                <Input
                  value={dpaText[i]}
                  onChange={(e) => {
                    setDpaText(dpaText.map((t, j) => (j === i ? e.target.value : t)));
                    clearError();
                  }}
                  className="font-mono text-xs"
                />
              </div>
            ))}
          </div>
          <div className="grid gap-3 sm:grid-cols-3">
            {GROUP_LABELS.map((label, i) => (
              <div key={label} className="space-y-1">
                <Label className="flex min-h-10 items-end">He response, {label} [barns]</Label>
                <Input
                  value={heText[i]}
                  onChange={(e) => {
                    setHeText(heText.map((t, j) => (j === i ? e.target.value : t)));
                    clearError();
                  }}
                  className="font-mono text-xs"
                />
              </div>
            ))}
          </div>

          <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-5">
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Exposure [s]</Label>
              <Input
                value={secondsText}
                onChange={(e) => {
                  setSecondsText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Target nuclide</Label>
              <Input
                value={targetText}
                onChange={(e) => {
                  setTargetText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">E_d [eV]</Label>
              <Input
                value={edText}
                onChange={(e) => {
                  setEdText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">b_arc</Label>
              <Input
                value={bArcText}
                onChange={(e) => {
                  setBArcText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">c_arc</Label>
              <Input
                value={cArcText}
                onChange={(e) => {
                  setCArcText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Compute damage metrics</Button>
          </div>

          <p className="text-sm text-muted-foreground">
            Coil fast-flux and weakest-link lifetime use the same group grid: the fast flux sums
            groups above the threshold, and each limit/rate pair is one ageing channel in shared
            units (comma-separated lists must match in length).
          </p>

          <div className="grid gap-3 sm:grid-cols-3">
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Fast threshold [MeV]</Label>
              <Input
                value={thresholdText}
                onChange={(e) => {
                  setThresholdText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Channel limits</Label>
              <Input
                value={limitsText}
                onChange={(e) => {
                  setLimitsText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Channel rates [/s]</Label>
              <Input
                value={ratesText}
                onChange={(e) => {
                  setRatesText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
          </div>

          {metrics && (
            <div className="space-y-2">
              <p className="text-sm font-medium">Spectral-fold results</p>
              <table className="text-sm">
                <tbody className="font-mono text-xs">
                  <tr>
                    <td className="pr-4">NRT-dpa</td>
                    <td>{fmt(metrics.nrtDpa)}</td>
                  </tr>
                  <tr>
                    <td className="pr-4">arc-dpa</td>
                    <td>{fmt(metrics.arcDpa)}</td>
                  </tr>
                  <tr>
                    <td className="pr-4">He production</td>
                    <td>{fmt(metrics.heAppm)} appm</td>
                  </tr>
                  <tr>
                    <td className="pr-4">He/dpa ratio</td>
                    <td>
                      {metrics.heDpaRatio === null
                        ? "— (zero dpa)"
                        : `${fmt(metrics.heDpaRatio)} appm/dpa`}
                    </td>
                  </tr>
                </tbody>
              </table>
              <p className="text-xs text-muted-foreground">
                NRT-dpa and arc-dpa use the same fold here; in practice the response cross sections
                embed the respective damage function (NRT 1975 vs Nordlund et al. 2018). arc-dpa
                &le; NRT-dpa once the efficiency ξ drops below 1.
              </p>
            </div>
          )}

          {curve && (
            <div className="space-y-2">
              <p className="text-sm font-medium">arc-dpa efficiency ξ vs PKA energy</p>
              <Plotly
                aspect="video"
                data={[
                  {
                    type: "scatter",
                    mode: "lines",
                    name: "ξ(T_dam(T))",
                    x: curve.tEv,
                    y: curve.xi,
                  },
                ]}
                layout={{
                  xaxis: { title: { text: "PKA energy T (eV)" }, type: "log" },
                  yaxis: { title: { text: "arc efficiency ξ" }, type: "linear" },
                  margin: { t: 16, r: 24, b: 48, l: 64 },
                  legend: { orientation: "h", y: -0.25 },
                }}
              />
            </div>
          )}

          {coilFlux && coilLife && (
            <div className="space-y-2">
              <p className="text-sm font-medium">Coil fast-flux and lifetime</p>
              <table className="text-sm">
                <tbody className="font-mono text-xs">
                  <tr>
                    <td className="pr-4">Fast flux</td>
                    <td>{fmt(coilFlux.fastFlux)} n/cm²/s</td>
                  </tr>
                  <tr>
                    <td className="pr-4">Fast fluence</td>
                    <td>{fmt(coilFlux.fastFluence)} n/cm²</td>
                  </tr>
                  <tr>
                    <td className="pr-4">Weakest-link life</td>
                    <td>
                      {Number.isFinite(coilLife.seconds)
                        ? `${fmt(coilLife.seconds)} s (channel ${coilLife.limiting})`
                        : "infinite (no ageing channel)"}
                    </td>
                  </tr>
                </tbody>
              </table>
            </div>
          )}
        </>
      )}
    </div>
  );
}
