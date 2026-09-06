import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { IsotxsSummary } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";
import { DataTable } from "@nukehub/docs-kit/components/mdx/DataTable";

const DEFAULT_ISOTXS = `ISOTXS 2
NUCLIDE U235 92235 2
1.1 2.2
NUCLIDE PU239 94239 2
4.4 5.5`;

export function DeterministicDemo() {
  const { wasm, ready, error } = useWasm();
  const [text, setText] = useState(DEFAULT_ISOTXS);
  const [summary, setSummary] = useState<IsotxsSummary | null>(null);
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
            </div>
          )}
        </>
      )}
    </div>
  );
}
