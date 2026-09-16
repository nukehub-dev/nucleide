import { useEffect, useRef, useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { SumOfFractionsJson } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";

// Default inventory sits exactly on the screening boundary: 0.06/0.1 (Co-60)
// + 0.03/0.1 (Cs-137) + 10/100 (Ni-63) = 0.6 + 0.3 + 0.1 = 1.0 — every edit
// tips the class one way or the other.
interface Row {
  id: number;
  nuclide: string;
  activity: string;
}

const DEFAULT_ROWS: Row[] = [
  { id: 1, nuclide: "Co60", activity: "0.06" },
  { id: 2, nuclide: "Cs137", activity: "0.03" },
  { id: 3, nuclide: "Ni63", activity: "10" },
];

export function ClearanceDemo() {
  const { wasm, ready, error } = useWasm();
  const [rows, setRows] = useState<Row[]>(DEFAULT_ROWS);
  const [limits, setLimits] = useState<Record<string, number>>({});
  const [outcome, setOutcome] = useState<SumOfFractionsJson | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);
  const nextId = useRef(4);

  useEffect(() => {
    if (!wasm || !ready) return;
    const table = wasm.euClearanceTable();
    const map: Record<string, number> = {};
    for (const entry of table.entries) {
      map[entry.nuclide] = entry.limitBqG;
    }
    setLimits(map);
  }, [wasm, ready]);

  function clearError() {
    setLocalError(null);
  }

  function normalize(name: string): string | null {
    if (!wasm) return null;
    try {
      return wasm.normalize_nuclide(name);
    } catch {
      return null;
    }
  }

  function updateRow(id: number, patch: Partial<Row>) {
    setRows(rows.map((r) => (r.id === id ? { ...r, ...patch } : r)));
    clearError();
  }

  function run() {
    if (!wasm) return;
    try {
      const inventory: Record<string, number> = {};
      for (const row of rows) {
        const name = row.nuclide.trim();
        if (name.length === 0) continue;
        const activity = parseFloat(row.activity);
        if (!Number.isFinite(activity) || activity < 0)
          throw new Error(`bad activity for ${name} \`${row.activity}\` (expected >= 0)`);
        const canonical = normalize(name);
        if (canonical === null) throw new Error(`unknown nuclide \`${name}\``);
        inventory[canonical] = activity;
      }
      if (Object.keys(inventory).length === 0)
        throw new Error("inventory is empty (add at least one nuclide row)");
      setOutcome(wasm.clearanceSumOfFractions(inventory));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setOutcome(null);
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
            Activities are massic activities (Bq/g) screened against the default EU 2013/59/Euratom
            Annex VII Table A (activity-concentration limits for solid material, in Bq/g). The
            screening criterion is the sum of fractions Σ A_i / CL_i ≤ 1 (boundary included).
          </p>

          <div className="space-y-2">
            {rows.map((row) => {
              const canonical = normalize(row.nuclide);
              const limit = canonical !== null ? limits[canonical] : undefined;
              const activity = parseFloat(row.activity);
              const fraction =
                limit !== undefined && Number.isFinite(activity) ? activity / limit : null;
              return (
                <div key={row.id} className="grid gap-3 sm:grid-cols-[1fr_1fr_1fr_auto]">
                  <div className="space-y-1">
                    <Label className="flex min-h-10 items-end">Nuclide</Label>
                    <Input
                      value={row.nuclide}
                      onChange={(e) => updateRow(row.id, { nuclide: e.target.value })}
                      className="font-mono text-xs"
                    />
                  </div>
                  <div className="space-y-1">
                    <Label className="flex min-h-10 items-end">Activity [Bq/g]</Label>
                    <Input
                      value={row.activity}
                      onChange={(e) => updateRow(row.id, { activity: e.target.value })}
                      className="font-mono text-xs"
                    />
                  </div>
                  <div className="space-y-1">
                    <Label className="flex min-h-10 items-end">Limit / fraction</Label>
                    <p className="flex min-h-10 items-center font-mono text-xs">
                      {limit !== undefined
                        ? `${limit} Bq/g / ${fraction === null ? "—" : fraction.toPrecision(3)}`
                        : canonical === null && row.nuclide.trim().length > 0
                          ? "not in table"
                          : "—"}
                    </p>
                  </div>
                  <div className="flex items-end">
                    <Button
                      variant="outline"
                      size="sm"
                      onClick={() => {
                        setRows(rows.filter((r) => r.id !== row.id));
                        clearError();
                      }}
                    >
                      Remove
                    </Button>
                  </div>
                </div>
              );
            })}
          </div>

          <div className="flex flex-wrap gap-2">
            <Button
              variant="outline"
              size="sm"
              onClick={() => {
                setRows([...rows, { id: nextId.current, nuclide: "", activity: "" }]);
                nextId.current += 1;
                clearError();
              }}
            >
              Add nuclide
            </Button>
            <Button onClick={run}>Screen inventory</Button>
          </div>

          {outcome && (
            <div className="space-y-2">
              <p className="text-sm font-medium">Sum-of-fractions screening</p>
              <table className="text-sm">
                <tbody className="font-mono text-xs">
                  <tr>
                    <td className="pr-4">Sum of fractions</td>
                    <td>{outcome.sum.toPrecision(5)}</td>
                  </tr>
                  <tr>
                    <td className="pr-4">Screening class</td>
                    <td>{outcome.class}</td>
                  </tr>
                  <tr>
                    <td className="pr-4">Dominant contributor</td>
                    <td>
                      {outcome.maxNuclide === null
                        ? "—"
                        : `${outcome.maxNuclide} (fraction ${outcome.maxFraction.toPrecision(3)})`}
                    </td>
                  </tr>
                </tbody>
              </table>
              <p className="text-xs text-muted-foreground">
                This is screening arithmetic, not a compliance decision: real clearance requires the
                governing regulatory table, material bookkeeping, and the national transposition of
                the directive.
              </p>
            </div>
          )}
        </>
      )}
    </div>
  );
}
