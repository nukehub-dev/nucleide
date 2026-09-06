import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type {
  AlaraDeckSummary,
  AlaraOutputSummary,
  FispactOutputSummary,
  R2sSummary,
} from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";
import { DataTable } from "@nukehub/docs-kit/components/mdx/DataTable";

type ActivationMode = "alara-deck" | "alara-output" | "fispact" | "r2s";

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

const DEFAULTS: Record<ActivationMode, string> = {
  "alara-deck": DEFAULT_DECK,
  "alara-output": DEFAULT_ALARA_OUTPUT,
  fispact: DEFAULT_FISPACT,
  r2s: DEFAULT_DECK,
};

const MODES: { value: ActivationMode; label: string }[] = [
  { value: "alara-deck", label: "ALARA deck" },
  { value: "alara-output", label: "ALARA output" },
  { value: "fispact", label: "FISPACT output" },
  { value: "r2s", label: "R2S workflow" },
];

export function ActivationDemo() {
  const { wasm, ready, error } = useWasm();
  const [mode, setMode] = useState<ActivationMode>("alara-deck");
  const [text, setText] = useState(DEFAULTS["alara-deck"]);
  const [runLbl, setRunLbl] = useState("demo");
  const [deck, setDeck] = useState<AlaraDeckSummary | null>(null);
  const [alaraOutput, setAlaraOutput] = useState<AlaraOutputSummary | null>(null);
  const [fispact, setFispact] = useState<FispactOutputSummary | null>(null);
  const [r2s, setR2s] = useState<R2sSummary | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function clearResults() {
    setDeck(null);
    setAlaraOutput(null);
    setFispact(null);
    setR2s(null);
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
        case "r2s":
          setR2s(wasm.r2sFromDeck(text));
          break;
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
