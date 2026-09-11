import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type {
  SerpentDepSummary,
  SerpentDetSpectrum,
  SerpentDetSummary,
  SerpentResSummary,
  SerpentVariableJson,
} from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";
import { DataTable } from "@nukehub/docs-kit/components/mdx/DataTable";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";

const BASE = import.meta.env.BASE_URL.endsWith("/")
  ? import.meta.env.BASE_URL
  : `${import.meta.env.BASE_URL}/`;

type ParserMode = "res" | "dep" | "det";

const SAMPLE_URLS: Record<ParserMode, string> = {
  res: `${BASE}data/serpent_res_sample.m`,
  dep: `${BASE}data/serpent_dep_sample.m`,
  det: `${BASE}data/serpent_det_sample.m`,
};

const DEFAULT_RES = `% Minimal _res.m-style block:
VERSION  (idx, [1: 13]) = 'Serpent 1.1.7' ;
TITLE    (idx, [1:  8]) = 'Untitled' ;
POP      (idx, 1)        = 5000 ;
CYCLES   (idx, 1)        = 130 ;
IMP_KEFF (idx, [1:  2])  = [ 1.24207 0.00053 ] ;`;

const DEFAULT_DEP = `% Minimal _dep.m-style inventory:
ZAI = [
922350
922380
];
NAMES = [
'U235'
'U238'
];
MAT_fuel_ADENS = [ 1.0000E-03 2.0000E-02 % U235
9.9000E-01 9.8000E-01 ]; % U238
MAT_fuel_BURNUP = [ 0.00000E+00 1.00874E-01 ];`;

const DEFAULT_DET = `% Minimal _det.m-style detector:
DETphi_EBINS = 2;
DETphi_VALS = 2;
DETphi = [
1 1 1 1 1 1 1 1 1 1 2.00638E-03 0.01824 2833
2 2 1 1 1 1 1 1 1 1 7.26203E-03 0.01195 6989
];
DETphiE = [ 0.0 1.0 2.0 % low bound
1.0 2.0 3.0 ]; % high bound`;

const DEFAULTS: Record<ParserMode, string> = {
  res: DEFAULT_RES,
  dep: DEFAULT_DEP,
  det: DEFAULT_DET,
};

async function fetchSample(url: string): Promise<string> {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  return res.text();
}

// WASM `keff` holds up to 2 entries from the first IMP_KEFF matrix row, but a
// well-formed row can have a single column — guard every access.
function formatKeff(keff: number[] | undefined): string {
  if (!keff || keff.length === 0) return " not present";
  const mean = ` ${keff[0].toFixed(5)}`;
  return keff.length >= 2 ? `${mean} ± ${keff[1].toExponential(2)}` : mean;
}

function KeffConvergenceChart({ history }: { history: [number, number][] }) {
  return (
    <Plotly
      aspect="video"
      data={[
        {
          type: "scatter",
          mode: "lines+markers",
          name: "IMP_KEFF",
          x: history.map((_, i) => i + 1),
          y: history.map(([mean]) => mean),
          error_y: {
            type: "data",
            array: history.map(([, err]) => err),
            visible: true,
          },
        },
      ]}
      layout={{
        xaxis: { title: { text: "Cycle (burnup block)" }, type: "linear" },
        yaxis: { title: { text: "k-eff (IMP_KEFF mean ± σ)" }, type: "linear" },
        showlegend: false,
        margin: { t: 16, r: 24, b: 48, l: 64 },
      }}
    />
  );
}

function DetSpectrumChart({ spectra }: { spectra: SerpentDetSpectrum[] }) {
  const [selected, setSelected] = useState<string | null>(null);
  const active = spectra.find((s) => s.name === selected) ?? spectra[0];
  if (!active || active.values.length === 0) return null;
  const hasEnergy = active.energy_mid.length === active.values.length;
  const x = hasEnergy ? active.energy_mid : active.values.map((_, i) => i + 1);
  return (
    <div className="space-y-2">
      {spectra.length > 1 && (
        <div className="flex flex-wrap gap-2">
          {spectra.map((s) => (
            <Button
              key={s.name}
              onClick={() => setSelected(s.name)}
              variant={s.name === active.name ? "default" : "outline"}
              size="sm"
            >
              {s.name}
            </Button>
          ))}
        </div>
      )}
      <Plotly
        aspect="video"
        data={[
          {
            type: "scatter",
            mode: "lines",
            name: active.name,
            x,
            y: active.values,
            error_y: {
              type: "data",
              array: active.values.map((v, i) => v * active.errors[i]),
              visible: true,
            },
          },
        ]}
        layout={{
          xaxis: {
            title: { text: hasEnergy ? "Energy (MeV)" : "Bin" },
            type: hasEnergy ? "log" : "linear",
          },
          yaxis: { title: { text: "Tally value" }, type: "log" },
          showlegend: false,
          margin: { t: 16, r: 24, b: 48, l: 64 },
        }}
      />
    </div>
  );
}

export function SerpentParser() {
  const { wasm, ready, error } = useWasm();
  const [mode, setMode] = useState<ParserMode>("res");
  const [text, setText] = useState(DEFAULT_RES);
  const [res, setRes] = useState<SerpentResSummary | null>(null);
  const [dep, setDep] = useState<SerpentDepSummary | null>(null);
  const [det, setDet] = useState<SerpentDetSummary | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);
  const [loadingSample, setLoadingSample] = useState(false);

  function clearError() {
    setLocalError(null);
  }

  function clearResults() {
    setRes(null);
    setDep(null);
    setDet(null);
  }

  function selectMode(next: ParserMode) {
    setMode(next);
    setText(DEFAULTS[next]);
    clearResults();
    clearError();
  }

  async function loadSample(which: ParserMode) {
    setLoadingSample(true);
    try {
      const sample = await fetchSample(SAMPLE_URLS[which]);
      setMode(which);
      setText(sample);
      clearResults();
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
    } finally {
      setLoadingSample(false);
    }
  }

  function run() {
    if (!wasm) return;
    clearResults();
    try {
      switch (mode) {
        case "res":
          setRes(wasm.parseSerpentRes(text));
          break;
        case "dep":
          setDep(wasm.parseSerpentDep(text));
          break;
        case "det":
          setDet(wasm.parseSerpentDet(text));
          break;
      }
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
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
            {(["res", "dep", "det"] as ParserMode[]).map((m) => (
              <Button
                key={m}
                onClick={() => selectMode(m)}
                variant={mode === m ? "default" : "outline"}
                size="sm"
              >
                {m}
              </Button>
            ))}
          </div>

          <Textarea
            value={text}
            onChange={(e) => {
              setText(e.target.value);
              clearError();
            }}
            className="font-mono text-xs"
          />

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Parse</Button>
            <Button variant="outline" onClick={() => loadSample(mode)} disabled={loadingSample}>
              {loadingSample ? "Loading…" : `Load sample ${mode}`}
            </Button>
          </div>

          {res && (
            <div className="space-y-3">
              <div className="grid gap-2 text-sm sm:grid-cols-3">
                <p>Version: {res.version ?? "unknown"}</p>
                <p>Variables: {res.variable_count}</p>
                <p>IMP_KEFF (block 1):{formatKeff(res.keff)}</p>
              </div>
              {res.title && <p className="text-sm">Title: {res.title}</p>}
              {res.keff_history && <KeffConvergenceChart history={res.keff_history} />}
              <VariablesTable variables={res.variables} />
            </div>
          )}

          {dep && (
            <div className="space-y-3">
              <div className="grid gap-2 text-sm sm:grid-cols-3">
                <p>Nuclides: {dep.nuclides.length}</p>
                <p>ZAI entries: {dep.zai.length}</p>
                <p>Variables: {dep.variable_count}</p>
              </div>
              {dep.nuclides.length > 0 && (
                <DataTable
                  data={dep.nuclides.map((name, i) => ({
                    nuclide: <span className="font-mono">{name}</span>,
                    zai: dep.zai[i] !== undefined ? String(dep.zai[i]) : "",
                  }))}
                  columns={[
                    { key: "nuclide", header: "Nuclide" },
                    { key: "zai", header: "ZAI", align: "right" },
                  ]}
                  pagination
                  pageSize={10}
                />
              )}
              <VariablesTable variables={dep.variables} />
            </div>
          )}

          {det && (
            <div className="space-y-3">
              <div className="grid gap-2 text-sm sm:grid-cols-2">
                <p>Detectors: {det.detectors.join(", ") || "none"}</p>
                <p>Variables: {det.variable_count}</p>
              </div>
              {det.spectra.length > 0 && <DetSpectrumChart spectra={det.spectra} />}
              <VariablesTable variables={det.variables} />
            </div>
          )}
        </>
      )}
    </div>
  );
}

function VariablesTable({ variables }: { variables: SerpentVariableJson[] }) {
  return (
    <div className="space-y-2">
      <p className="text-sm font-medium">Variables ({variables.length})</p>
      <DataTable
        data={variables.map((v) => ({
          name: <span className="font-mono">{v.name}</span>,
          kind: v.kind,
          shape: v.shape,
          value: v.value === undefined ? "" : String(v.value),
        }))}
        columns={[
          { key: "name", header: "Name" },
          { key: "kind", header: "Kind" },
          { key: "shape", header: "Shape", align: "right" },
          { key: "value", header: "Scalar value" },
        ]}
        pagination
        pageSize={10}
      />
    </div>
  );
}
