import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { IsotxsSummary, RtfluxSummary } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";
import { DataTable } from "@nukehub/docs-kit/components/mdx/DataTable";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";
import { Select } from "@nukehub/docs-kit/components/ui/Select";

const DEFAULT_ISOTXS = `ISOTXS 2
NUCLIDE U235 92235 2
1.1 2.2
NUCLIDE PU239 94239 2
4.4 5.5`;

const DEFAULT_RTFLUX = `RTFLUX 2 3
1.0 2.0 3.0
4.0 5.0 6.0`;

export function DeterministicDemo() {
  const { wasm, ready, error } = useWasm();
  const [text, setText] = useState(DEFAULT_ISOTXS);
  const [summary, setSummary] = useState<IsotxsSummary | null>(null);
  const [rtfluxText, setRtfluxText] = useState(DEFAULT_RTFLUX);
  const [rtfluxKind, setRtfluxKind] = useState("rtflux");
  const [rtflux, setRtflux] = useState<RtfluxSummary | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function run() {
    if (!wasm) return;
    setSummary(null);
    try {
      setSummary(wasm.parseIsotxs(text));
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
          </div>

          {summary && (
            <div className="space-y-3">
              <p className="text-sm">
                Nuclides ({summary.nuclides.length}), energy groups: {summary.groups}
              </p>
              <DataTable
                data={summary.nuclides.map((n) => ({
                  label: <span className="font-mono">{n.label}</span>,
                  zaid: <span className="font-mono">{n.zaid}</span>,
                  groups: n.groups,
                  totalXs: n.total_xs.map((v) => v.toFixed(2)).join(", "),
                }))}
                columns={[
                  { key: "label", header: "Nuclide" },
                  { key: "zaid", header: "ZAID" },
                  { key: "groups", header: "Groups", align: "right" },
                  { key: "totalXs", header: "Total xs" },
                ]}
                pagination
                pageSize={10}
              />
              <Plotly
                aspect="video"
                data={summary.nuclides.map((n) => ({
                  type: "bar",
                  name: n.label,
                  x: n.total_xs.map((_, g) => `g${g + 1}`),
                  y: n.total_xs,
                }))}
                layout={{
                  barmode: "group",
                  xaxis: { title: { text: "Energy group" } },
                  yaxis: { title: { text: "Total xs" } },
                  margin: { t: 16, r: 16, b: 48, l: 64 },
                  legend: { orientation: "h", y: -0.25 },
                }}
              />
            </div>
          )}

          <div className="space-y-2 border-t border-border/50 pt-4">
            <p className="text-sm font-medium">RTFLUX fluxes (PARTISN stays Python-only)</p>
            <Textarea
              value={rtfluxText}
              onChange={(e) => {
                setRtfluxText(e.target.value);
                clearError();
              }}
              className="font-mono text-xs"
            />
            <div className="flex flex-wrap items-end gap-2">
              <Select
                value={rtfluxKind}
                onChange={(v) => {
                  setRtfluxKind(v);
                  clearError();
                }}
                options={[
                  { value: "rtflux", label: "RTFLUX" },
                  { value: "atflux", label: "ATFLUX" },
                  { value: "rzflux", label: "RZFLUX" },
                ]}
              />
              <Button
                onClick={() => {
                  if (!wasm) return;
                  try {
                    setRtflux(wasm.parseRtflux(rtfluxText, rtfluxKind));
                    clearError();
                  } catch (e) {
                    setLocalError(e instanceof Error ? e.message : String(e));
                    setRtflux(null);
                  }
                }}
              >
                Parse RTFLUX
              </Button>
            </div>
            {rtflux && (
              <p className="text-sm">
                Flux kind: <span className="font-mono">{rtflux.kind}</span>, groups: {rtflux.groups}
                , points: {rtflux.npoints}
              </p>
            )}
          </div>
        </>
      )}
    </div>
  );
}
