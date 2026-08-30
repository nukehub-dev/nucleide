import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Select } from "@nukehub/docs-kit/components/ui/Select";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";

const BASE = import.meta.env.BASE_URL.endsWith("/")
  ? import.meta.env.BASE_URL
  : `${import.meta.env.BASE_URL}/`;
const CHAIN_SAMPLE_URL = `${BASE}data/chain_simple.xml`;

const DEFAULT_CHAIN = `<?xml version="1.0"?>
<depletion_chain>
  <nuclide name="I135" half_life="2.36520E+04">
    <decay type="beta" target="Xe135" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="Xe135" half_life="3.29040E+04">
    <decay type="beta" target="Cs135" branching_ratio="1.0"/>
  </nuclide>
  <nuclide name="Cs135" />
</depletion_chain>`;

const ORDER_OPTIONS = [
  { value: "16", label: "16" },
  { value: "48", label: "48" },
];

const MAX_BURNUP_STEPS = 200;
const DEFAULT_BURNUP_STEPS = 50;

export function DepletionStep() {
  const { wasm, ready, error } = useWasm();
  const [xml, setXml] = useState(DEFAULT_CHAIN);
  const [n0Input, setN0Input] = useState("I135 1e15\nXe135 0\nCs135 0");
  const [dt, setDt] = useState(86400.0);
  const [order, setOrder] = useState<16 | 48>(48);
  const [result, setResult] = useState<Record<string, number> | null>(null);
  const [burnupSteps, setBurnupSteps] = useState(DEFAULT_BURNUP_STEPS);
  const [burnup, setBurnup] = useState<BurnupCurve | null>(null);
  const [burnupBusy, setBurnupBusy] = useState(false);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  async function loadSampleChain() {
    try {
      const res = await fetch(CHAIN_SAMPLE_URL);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const sample = await res.text();
      setXml(sample);
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
    }
  }

  function parseN0(): Record<string, number> {
    const n0: Record<string, number> = {};
    for (const line of n0Input.split("\n")) {
      const [name, count] = line.trim().split(/\s+/);
      if (name && count) n0[name] = parseFloat(count);
    }
    return n0;
  }

  function run() {
    if (!wasm) return;
    try {
      const chain = wasm.WasmChain.fromXml(xml);
      const n0 = parseN0();
      const out = wasm.deplete(chain, n0, dt, {}, order);
      setResult(out);
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setResult(null);
    }
  }

  function runBurnup() {
    if (!wasm) return;
    const steps = Math.min(Math.max(1, burnupSteps), MAX_BURNUP_STEPS);
    setBurnupSteps(steps);
    setBurnupBusy(true);
    try {
      const chain = wasm.WasmChain.fromXml(xml);
      let n = parseN0();
      const times: number[] = [0];
      const series: Record<string, number[]> = {};
      for (const name of Object.keys(n)) {
        series[name] = [n[name]];
      }

      for (let i = 1; i <= steps; i++) {
        n = wasm.deplete(chain, n, dt, {}, order);
        times.push(i * dt);
        for (const [name, value] of Object.entries(n)) {
          if (!series[name]) series[name] = [];
          series[name].push(value);
        }
      }

      // Align series length (new nuclides may appear mid-curve).
      for (const values of Object.values(series)) {
        while (values.length < times.length) {
          values.push(0);
        }
      }

      setBurnup({ times, series });
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setBurnup(null);
    } finally {
      setBurnupBusy(false);
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
          <div className="space-y-2">
            <div className="flex flex-wrap items-center justify-between gap-2">
              <Label>Depletion chain XML</Label>
              <Button variant="outline" size="sm" onClick={loadSampleChain}>
                Load sample chain
              </Button>
            </div>
            <Textarea
              value={xml}
              onChange={(e) => {
                setXml(e.target.value);
                clearError();
              }}
              className="font-mono text-xs"
            />
          </div>

          <div className="space-y-2">
            <Label>Initial atom counts</Label>
            <p className="text-xs text-muted-foreground">One "Nuclide count" per line.</p>
            <Textarea
              value={n0Input}
              onChange={(e) => {
                setN0Input(e.target.value);
                clearError();
              }}
              className="font-mono text-sm"
            />
          </div>

          <div className="grid gap-3 sm:grid-cols-2">
            <div className="space-y-1">
              <Label>Time step (s)</Label>
              <Input
                type="number"
                step="any"
                value={dt}
                onChange={(e) => {
                  setDt(parseFloat(e.target.value));
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label>CRAM order</Label>
              <Select
                value={String(order)}
                onChange={(v) => {
                  setOrder(parseInt(v) as 16 | 48);
                  clearError();
                }}
                options={ORDER_OPTIONS}
                className="min-w-[100px]"
              />
            </div>
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Deplete</Button>
            <Button variant="outline" onClick={runBurnup} disabled={burnupBusy}>
              {burnupBusy ? "Running burnup curve…" : "Burnup curve"}
            </Button>
          </div>

          <div className="grid gap-3 sm:grid-cols-2">
            <div className="space-y-1">
              <Label>Burnup steps (max {MAX_BURNUP_STEPS})</Label>
              <Input
                type="number"
                min={1}
                max={MAX_BURNUP_STEPS}
                value={burnupSteps}
                onChange={(e) => {
                  setBurnupSteps(parseInt(e.target.value));
                  clearError();
                }}
              />
            </div>
          </div>

          {result && Object.keys(result).length === 0 && (
            <p className="text-sm text-muted-foreground">No nuclides in result.</p>
          )}

          {result && Object.keys(result).length > 0 && (
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b border-border">
                  <th className="py-1 text-left">Nuclide</th>
                  <th className="py-1 text-right">Atom count</th>
                </tr>
              </thead>
              <tbody>
                {Object.entries(result).map(([nuc, count]) => (
                  <tr key={nuc} className="border-b border-border/50">
                    <td className="py-1 font-mono">{nuc}</td>
                    <td className="py-1 text-right font-mono">{count.toExponential(4)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          )}

          {burnup && <BurnupPlot burnup={burnup} />}
        </>
      )}
    </div>
  );
}

interface BurnupCurve {
  times: number[];
  series: Record<string, number[]>;
}

function BurnupPlot({ burnup }: { burnup: BurnupCurve }) {
  const traces = Object.entries(burnup.series).map(([name, values]) => ({
    type: "scatter" as const,
    mode: "lines" as const,
    name,
    x: burnup.times,
    y: values,
  }));

  return (
    <Plotly
      aspect="video"
      data={traces}
      layout={{
        xaxis: { title: { text: "Time (s)" }, type: "linear" },
        yaxis: { title: { text: "Atom count" }, type: "log" },
        margin: { t: 16, r: 16, b: 48, l: 64 },
        legend: { orientation: "h", y: -0.25 },
      }}
    />
  );
}
