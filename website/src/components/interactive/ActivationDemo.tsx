import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type {
  AlaraDeckSummary,
  AlaraOutputSummary,
  FispactOutputSummary,
  OrigenTape5Summary,
  OrigenTape6Summary,
  OrigenTape9Summary,
  R2sSummary,
  SnapshotBundleJson,
  SnapshotInputJson,
} from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";
import { DataTable } from "@nukehub/docs-kit/components/mdx/DataTable";

type ActivationMode =
  | "alara-deck"
  | "alara-output"
  | "fispact"
  | "origen-tape5"
  | "origen-tape6"
  | "origen-tape9"
  | "r2s"
  | "r2s-snapshot";

const DEFAULT_DECK = `geometry rectangular
mat_loading
inner_zone inner_mix
end
mixture inner_mix
material WATER 1.0 1.0
end
flux flux_1 data/fluxin1 1.0 1 default
schedule 1_year
1 y flux_1 steady_state 0 s
end
pulsehistory steady_state
1 0 s
end
cooling
1 d
end`;

const DEFAULT_ALARA_OUTPUT = `***Number Density [atoms/cm3]
Interval #1 (Zone: inner_zone)
isotope t_1/2(s) shutdown 1 d
h-1 stable 6.6354E+21 6.0000E+21
total none 1.0000E+22 9.0000E+21`;

const DEFAULT_FISPACT = `SYNTHETIC INVENTORY DEMO RUN
TIME = 0.0 SECS
NUCLIDE ATOMS ACTIVITY HEAT
h-3 1.0000E+20 1.2000E+09 3.5000E-03
total 1.3000E+20 5.5000E+09 7.6600E-02`;

const DEFAULT_ORIGEN_TAPE5 = `# SYNTHETIC ORIGEN TAPE5 sample - simplified input echo for parser tests (not real ORIGEN output).
SYNTHETIC PWR PIN - TAPE5 SAMPLE
CASE 1 - BASE DEPLETION
FLUX= 3.0e13 DAYS= 100.0
FLUX= 0.0 DAYS= 30.0
MAT fuel
U235 10.5
U238 1000.0
Pu239 0.25
MAT clad
Zr90 50.0`;

const DEFAULT_ORIGEN_TAPE6 = `# SYNTHETIC ORIGEN TAPE6 sample - simplified inventory for parser tests (not real ORIGEN output).
U235 10.2 8.16e5
U238 995.0 1.23e4
Pu239 0.245 5.63e8
Cs137 0.012 3.84e10`;

const DEFAULT_ORIGEN_TAPE9 = `# SYNTHETIC ORIGEN TAPE9 sample - simplified decay constants for parser tests (not real ORIGEN data).
U235 3.1209e-17
U238 4.9161e-18
Pu239 9.1105e-13

# short-lived fission products and activation nuclides below
Cs137 7.3217e-10
Co60 4.1674e-09`;

const DEFAULT_SNAPSHOT = `{
  "zones": [
    {"id": "zone1", "volumeCm3": 1000, "composition": {"U235": 0.001, "U238": 0.02}},
    {"id": "zone2", "volumeCm3": 500, "composition": {"H1": 0.06, "O16": 0.03}}
  ],
  "fluxDefs": [{"name": "flux_1", "file": "flux1", "scale": 1.0}],
  "coolingS": [86400]
}`;

const DEFAULTS: Record<ActivationMode, string> = {
  "alara-deck": DEFAULT_DECK,
  "alara-output": DEFAULT_ALARA_OUTPUT,
  fispact: DEFAULT_FISPACT,
  "origen-tape5": DEFAULT_ORIGEN_TAPE5,
  "origen-tape6": DEFAULT_ORIGEN_TAPE6,
  "origen-tape9": DEFAULT_ORIGEN_TAPE9,
  r2s: DEFAULT_DECK,
  "r2s-snapshot": DEFAULT_SNAPSHOT,
};

const MODES: { value: ActivationMode; label: string }[] = [
  { value: "alara-deck", label: "ALARA deck" },
  { value: "alara-output", label: "ALARA output" },
  { value: "fispact", label: "FISPACT output" },
  { value: "origen-tape5", label: "ORIGEN TAPE5" },
  { value: "origen-tape6", label: "ORIGEN TAPE6" },
  { value: "origen-tape9", label: "ORIGEN TAPE9" },
  { value: "r2s", label: "R2S workflow" },
  { value: "r2s-snapshot", label: "R2S snapshot" },
];

export function ActivationDemo() {
  const { wasm, ready, error } = useWasm();
  const [mode, setMode] = useState<ActivationMode>("alara-deck");
  const [text, setText] = useState(DEFAULTS["alara-deck"]);
  const [runLbl, setRunLbl] = useState("demo");
  const [deck, setDeck] = useState<AlaraDeckSummary | null>(null);
  const [alaraOutput, setAlaraOutput] = useState<AlaraOutputSummary | null>(null);
  const [fispact, setFispact] = useState<FispactOutputSummary | null>(null);
  const [tape5, setTape5] = useState<OrigenTape5Summary | null>(null);
  const [tape6, setTape6] = useState<OrigenTape6Summary | null>(null);
  const [tape9, setTape9] = useState<OrigenTape9Summary | null>(null);
  const [r2s, setR2s] = useState<R2sSummary | null>(null);
  const [snapshot, setSnapshot] = useState<SnapshotBundleJson | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function clearResults() {
    setDeck(null);
    setAlaraOutput(null);
    setFispact(null);
    setTape5(null);
    setTape6(null);
    setTape9(null);
    setR2s(null);
    setSnapshot(null);
  }

  function selectMode(next: ActivationMode) {
    setMode(next);
    setText(DEFAULTS[next]);
    clearResults();
    clearError();
  }

  function run() {
    if (!wasm) return;
    clearResults();
    try {
      switch (mode) {
        case "alara-deck":
          setDeck(wasm.parseAlaraDeck(text));
          break;
        case "alara-output":
          setAlaraOutput(wasm.parseAlaraOutput(text, runLbl));
          break;
        case "fispact":
          setFispact(wasm.parseFispactOutput(text, runLbl));
          break;
        case "origen-tape5":
          setTape5(wasm.parseOrigenTape5(text));
          break;
        case "origen-tape6":
          setTape6(wasm.parseOrigenTape6(text));
          break;
        case "origen-tape9":
          setTape9(wasm.parseOrigenTape9(text));
          break;
        case "r2s":
          setR2s(wasm.r2sFromDeck(text));
          break;
        case "r2s-snapshot": {
          let parsed: unknown;
          try {
            parsed = JSON.parse(text);
          } catch {
            throw new Error("bad snapshot JSON (expected zones/fluxDefs/coolingS)");
          }
          setSnapshot(wasm.r2sFromSnapshot(parsed as SnapshotInputJson));
          break;
        }
      }
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
    }
  }

  const displayError = error ?? localError;
  const needsRunLbl = mode === "alara-output" || mode === "fispact";

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
            {MODES.map((m) => (
              <Button
                key={m.value}
                onClick={() => selectMode(m.value)}
                variant={mode === m.value ? "default" : "outline"}
                size="sm"
              >
                {m.label}
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

          {needsRunLbl && (
            <div className="w-48 space-y-1">
              <Label>Run label</Label>
              <Input
                value={runLbl}
                onChange={(e) => {
                  setRunLbl(e.target.value);
                  clearError();
                }}
              />
            </div>
          )}

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Parse</Button>
          </div>

          {deck && (
            <div className="space-y-3">
              <p className="text-sm">Blocks: {deck.block_kinds.join(", ")}</p>
              <div className="space-y-2">
                <p className="text-sm font-medium">Mixtures ({deck.mixtures.length})</p>
                <DataTable
                  data={deck.mixtures.map((m) => ({
                    name: <span className="font-mono">{m.name}</span>,
                    entries: m.entries.join("; "),
                  }))}
                  columns={[
                    { key: "name", header: "Mixture" },
                    { key: "entries", header: "Entries" },
                  ]}
                />
              </div>
              <div className="space-y-2">
                <p className="text-sm font-medium">Fluxes ({deck.fluxes.length})</p>
                <DataTable
                  data={deck.fluxes.map((f) => ({
                    name: <span className="font-mono">{f.name}</span>,
                    file: <span className="font-mono">{f.file}</span>,
                    scale: f.scale.toExponential(2),
                    format: f.format,
                  }))}
                  columns={[
                    { key: "name", header: "Flux" },
                    { key: "file", header: "File" },
                    { key: "scale", header: "Scale", align: "right" },
                    { key: "format", header: "Format" },
                  ]}
                />
              </div>
              <p className="text-sm">
                Cooling times (s): {deck.cooling_times_s.join(", ") || "none"}
              </p>
              <div className="space-y-2">
                <p className="text-sm font-medium">Schedules ({deck.schedules.length})</p>
                <DataTable
                  data={deck.schedules.map((s) => ({
                    name: <span className="font-mono">{s.name}</span>,
                    items: s.items.map((tokens) => tokens.join(" ")).join(" | "),
                  }))}
                  columns={[
                    { key: "name", header: "Schedule" },
                    { key: "items", header: "Items" },
                  ]}
                />
              </div>
            </div>
          )}

          {alaraOutput && (
            <OutputTable
              rows={alaraOutput.rows}
              variables={alaraOutput.variables}
              blocks={alaraOutput.blocks}
            />
          )}

          {fispact && (
            <OutputTable rows={fispact.rows} variables={fispact.variables} blocks={null} />
          )}

          {tape5 && (
            <div className="space-y-3">
              <div className="grid gap-2 text-sm sm:grid-cols-3">
                <p>Titles: {tape5.titles.length}</p>
                <p>Steps: {tape5.steps.length}</p>
                <p>Materials: {tape5.materials.length}</p>
              </div>
              {tape5.titles.length > 0 && (
                <p className="text-sm">Title cards: {tape5.titles.join(" | ")}</p>
              )}
              <DataTable
                data={tape5.steps.map((s) => ({
                  flux: s.flux.toExponential(3),
                  days: s.days.toFixed(1),
                }))}
                columns={[
                  { key: "flux", header: "Flux (n/cm²/s)", align: "right" },
                  { key: "days", header: "Days", align: "right" },
                ]}
              />
              <div className="space-y-2">
                <p className="text-sm font-medium">Materials ({tape5.materials.length})</p>
                <DataTable
                  data={tape5.materials.flatMap((m) =>
                    m.entries.map((e) => ({
                      material: <span className="font-mono">{m.name}</span>,
                      nuclide: <span className="font-mono">{e.nuclide}</span>,
                      grams: e.grams.toExponential(4),
                    })),
                  )}
                  columns={[
                    { key: "material", header: "Material" },
                    { key: "nuclide", header: "Nuclide" },
                    { key: "grams", header: "Grams", align: "right" },
                  ]}
                  pagination
                  pageSize={10}
                />
              </div>
            </div>
          )}

          {tape6 && (
            <div className="space-y-3">
              <p className="text-sm">
                Records: {tape6.records.length}, total activity:{" "}
                {tape6.total_activity_bq.toExponential(4)} Bq
              </p>
              <DataTable
                data={tape6.records.map((r) => ({
                  nuclide: <span className="font-mono">{r.nuclide}</span>,
                  grams: r.grams.toExponential(4),
                  activity: r.activity_bq.toExponential(4),
                }))}
                columns={[
                  { key: "nuclide", header: "Nuclide" },
                  { key: "grams", header: "Grams", align: "right" },
                  { key: "activity", header: "Activity (Bq)", align: "right" },
                ]}
              />
            </div>
          )}

          {tape9 && (
            <div className="space-y-3">
              <p className="text-sm">Entries: {tape9.entries.length}</p>
              <DataTable
                data={tape9.entries.map((e) => ({
                  nuclide: <span className="font-mono">{e.nuclide}</span>,
                  lambda: e.decay_const.toExponential(4),
                }))}
                columns={[
                  { key: "nuclide", header: "Nuclide" },
                  { key: "lambda", header: "λ (s⁻¹)", align: "right" },
                ]}
              />
            </div>
          )}

          {r2s && (
            <div className="space-y-3">
              <div className="grid gap-2 text-sm sm:grid-cols-3">
                <p>Top schedule: {r2s.top_schedule}</p>
                <p>Total time (s): {r2s.total_s.toExponential(4)}</p>
                <p>Cooling (s): {r2s.cooling_s.join(", ") || "none"}</p>
              </div>
              <DataTable
                data={r2s.steps.map((s) => ({
                  zone: <span className="font-mono">{s.zone}</span>,
                  flux: <span className="font-mono">{s.flux}</span>,
                }))}
                columns={[
                  { key: "zone", header: "Zone" },
                  { key: "flux", header: "Flux" },
                ]}
              />
            </div>
          )}

          {snapshot && (
            <div className="space-y-3">
              <div className="grid gap-2 text-sm sm:grid-cols-3">
                <p>Top schedule: {snapshot.workflow.top_schedule}</p>
                <p>Total time (s): {snapshot.workflow.total_s.toExponential(4)}</p>
                <p>Decks: {snapshot.decks.length}</p>
              </div>
              <DataTable
                data={snapshot.workflow.steps.map((s) => ({
                  zone: <span className="font-mono">{s.zone}</span>,
                  flux: <span className="font-mono">{s.flux}</span>,
                }))}
                columns={[
                  { key: "zone", header: "Zone" },
                  { key: "flux", header: "Flux" },
                ]}
              />
              <div className="space-y-1">
                <Label>Template deck</Label>
                <pre className="overflow-x-auto rounded-lg border border-border/50 bg-muted/30 p-3 font-mono text-xs whitespace-pre-wrap">
                  {snapshot.deck}
                </pre>
              </div>
            </div>
          )}
        </>
      )}
    </div>
  );
}

interface OutputRow {
  time_s: number;
  time_label: string;
  nuclide: string;
  variable: string;
  var_unit: string;
  value: number;
  block: string;
  block_name: string;
}

function OutputTable({
  rows,
  variables,
  blocks,
}: {
  rows: OutputRow[];
  variables: string[];
  blocks: string[] | null;
}) {
  return (
    <div className="space-y-3">
      <p className="text-sm">
        Rows: {rows.length}, variables: {variables.join(", ")}
        {blocks !== null && `, blocks: ${blocks.join(", ")}`}
      </p>
      <DataTable
        data={rows.map((r) => ({
          nuclide: <span className="font-mono">{r.nuclide}</span>,
          variable: r.variable,
          time: r.time_label,
          block: r.block_name,
          value: r.value.toExponential(4),
          unit: r.var_unit,
        }))}
        columns={[
          { key: "nuclide", header: "Nuclide" },
          { key: "variable", header: "Variable" },
          { key: "time", header: "Cooling" },
          { key: "block", header: "Block" },
          { key: "value", header: "Value", align: "right" },
          { key: "unit", header: "Unit" },
        ]}
        pagination
        pageSize={10}
      />
    </div>
  );
}
