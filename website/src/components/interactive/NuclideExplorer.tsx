import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { NuclideInfo } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";

export function NuclideExplorer() {
  const { wasm, ready, error } = useWasm();
  const [name, setName] = useState("U235");
  const [info, setInfo] = useState<NuclideInfo | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function lookup() {
    if (!wasm) return;
    try {
      const nuc = new wasm.WasmNuclide(name);
      setInfo(nuc.toObject());
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setInfo(null);
    }
  }

  const displayError = error ?? localError;

  return (
    <div className="rounded-xl border border-border/50 bg-background p-4 space-y-3">
      {!ready && <p className="text-sm text-muted-foreground">Loading Nucleide WASM…</p>}
      {ready && (
        <>
          <div className="flex flex-wrap items-end gap-2">
            <div className="space-y-1">
              <Label htmlFor="nuclide-name">Nuclide</Label>
              <Input
                id="nuclide-name"
                type="text"
                value={name}
                onChange={(e) => {
                  setName(e.target.value);
                  setLocalError(null);
                }}
                className="w-32"
              />
            </div>
            <Button onClick={lookup}>Look up</Button>
          </div>

          {displayError && (
            <div className="flex items-start justify-between gap-2 rounded-lg border border-red-200 bg-red-50 p-3 text-sm text-red-700 dark:border-red-900 dark:bg-red-950 dark:text-red-300">
              <span>WASM error: {displayError}</span>
              <button
                onClick={() => setLocalError(null)}
                className="font-bold leading-none"
                aria-label="Dismiss error"
              >
                ×
              </button>
            </div>
          )}

          {info && (
            <table className="w-full text-sm">
              <tbody>
                {Object.entries(info).map(([key, value]) => (
                  <tr key={key} className="border-b border-border/50">
                    <td className="py-1 font-medium">{key}</td>
                    <td className="py-1 text-right font-mono">
                      {typeof value === "number" ? value.toString() : (value ?? "—")}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          )}
        </>
      )}
    </div>
  );
}
