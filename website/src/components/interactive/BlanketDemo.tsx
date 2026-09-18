import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { TbrScalars } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";

// Stellaris-adjacent synthetic preset (round caller tallies, no plant data).
const DEFAULT_BRED = "1.12";
const DEFAULT_SOURCE = "1.0";
const DEFAULT_PORTS = "0.03, 0.02";
const DEFAULT_POWER = "500";
const DEFAULT_REQUIRED = "1.05";

export function BlanketDemo() {
  const { wasm, ready, error } = useWasm();
  const [bredText, setBredText] = useState(DEFAULT_BRED);
  const [sourceText, setSourceText] = useState(DEFAULT_SOURCE);
  const [portsText, setPortsText] = useState(DEFAULT_PORTS);
  const [powerText, setPowerText] = useState(DEFAULT_POWER);
  const [requiredText, setRequiredText] = useState(DEFAULT_REQUIRED);
  const [scalars, setScalars] = useState<TbrScalars | null>(null);
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
      const tritonsBred = parseEntry(bredText, "tritons bred");
      const sourceNeutrons = parseEntry(sourceText, "source neutrons");
      const portFractions =
        portsText.trim() === ""
          ? []
          : portsText.split(",").map((t, i) => {
              const v = parseEntry(t.trim(), `port ${i + 1}`);
              if (v < 0 || v >= 1)
                throw new Error(`bad port ${i + 1} \`${t.trim()}\` (expected in [0, 1))`);
              return v;
            });
      const fusionPowerMw = parseEntry(powerText, "fusion power");
      const requiredTbr =
        requiredText.trim() === "" ? undefined : parseEntry(requiredText, "required TBR");
      setScalars(
        wasm.tbrScalars(tritonsBred, sourceNeutrons, portFractions, fusionPowerMw, requiredTbr),
      );
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setScalars(null);
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
            TBR bookkeeping over caller tallies: raw TBR, the multiplicative port-penalty haircut,
            the breeding margin, and the daily tritium burn and surplus. Ports are a comma-separated
            list of fractional coverage losses. All numbers are synthetic.
          </p>

          <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-5">
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Tritons bred</Label>
              <Input
                value={bredText}
                onChange={(e) => {
                  setBredText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Source neutrons</Label>
              <Input
                value={sourceText}
                onChange={(e) => {
                  setSourceText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Port fractions</Label>
              <Input
                value={portsText}
                onChange={(e) => {
                  setPortsText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Fusion power [MW]</Label>
              <Input
                value={powerText}
                onChange={(e) => {
                  setPowerText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Required TBR</Label>
              <Input
                value={requiredText}
                onChange={(e) => {
                  setRequiredText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Compute TBR bookkeeping</Button>
          </div>

          {scalars && (
            <div className="space-y-2">
              <p className="text-sm font-medium">Bookkeeping results</p>
              <table className="text-sm">
                <tbody className="font-mono text-xs">
                  <tr>
                    <td className="pr-4">Raw TBR</td>
                    <td>{scalars.rawTbr.toPrecision(6)}</td>
                  </tr>
                  <tr>
                    <td className="pr-4">Effective TBR</td>
                    <td>{scalars.effectiveTbr.toPrecision(6)}</td>
                  </tr>
                  <tr>
                    <td className="pr-4">Breeding margin</td>
                    <td>{scalars.margin.toPrecision(6)}</td>
                  </tr>
                  <tr>
                    <td className="pr-4">Meets requirement</td>
                    <td>
                      {scalars.meets === undefined
                        ? "— (no threshold)"
                        : scalars.meets
                          ? "yes"
                          : "no"}
                    </td>
                  </tr>
                  <tr>
                    <td className="pr-4">Tritium burn</td>
                    <td>{scalars.burnGPerDay.toPrecision(6)} g/day</td>
                  </tr>
                  <tr>
                    <td className="pr-4">Net surplus</td>
                    <td>{scalars.surplusGPerDay.toPrecision(6)} g/day</td>
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
