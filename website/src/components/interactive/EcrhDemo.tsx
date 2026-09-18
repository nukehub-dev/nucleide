import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { EcrhScalars } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";

// ITER-class 170 GHz fundamental at 6 T: the crate's worked example
// (synthetic beamline inputs, no machine data).
const DEFAULT_FREQ = "170";
const DEFAULT_HARMONIC = "1";
const DEFAULT_B = "6";
const DEFAULT_TE = "10";

export function EcrhDemo() {
  const { wasm, ready, error } = useWasm();
  const [freqText, setFreqText] = useState(DEFAULT_FREQ);
  const [harmonicText, setHarmonicText] = useState(DEFAULT_HARMONIC);
  const [bText, setBText] = useState(DEFAULT_B);
  const [teText, setTeText] = useState(DEFAULT_TE);
  const [scalars, setScalars] = useState<EcrhScalars | null>(null);
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
      const frequencyGhz = parseEntry(freqText, "frequency");
      const harmonic = parseEntry(harmonicText, "harmonic");
      const bT = parseEntry(bText, "field");
      const teKev = teText.trim() === "" ? undefined : parseEntry(teText, "electron temperature");
      setScalars(wasm.ecrhScalars(frequencyGhz, harmonic, bT, teKev));
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
            Closed-form ECRH accessibility scalars: the cold cyclotron resonance, the relativistic
            shift at finite electron temperature, and the O1/X1 cut-off densities. Leave the
            temperature blank for the cold result.
          </p>

          <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-4">
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Frequency [GHz]</Label>
              <Input
                value={freqText}
                onChange={(e) => {
                  setFreqText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Harmonic (1–10)</Label>
              <Input
                value={harmonicText}
                onChange={(e) => {
                  setHarmonicText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Field [T]</Label>
              <Input
                value={bText}
                onChange={(e) => {
                  setBText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Electron temp [keV]</Label>
              <Input
                value={teText}
                onChange={(e) => {
                  setTeText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Compute ECRH access</Button>
          </div>

          {scalars && (
            <div className="space-y-2">
              <p className="text-sm font-medium">Accessibility scalars</p>
              <table className="text-sm">
                <tbody className="font-mono text-xs">
                  <tr>
                    <td className="pr-4">Gyrofrequency</td>
                    <td>{scalars.gyrofrequencyGhz.toPrecision(6)} GHz</td>
                  </tr>
                  <tr>
                    <td className="pr-4">Cold resonant field</td>
                    <td>{scalars.resonantFieldT.toPrecision(6)} T</td>
                  </tr>
                  <tr>
                    <td className="pr-4">Relativistic field</td>
                    <td>
                      {scalars.relativisticFieldT === undefined
                        ? "— (no temperature)"
                        : `${scalars.relativisticFieldT.toPrecision(6)} T`}
                    </td>
                  </tr>
                  <tr>
                    <td className="pr-4">O1 cut-off density</td>
                    <td>{scalars.o1CutoffDensityM3.toExponential(4)} m⁻³</td>
                  </tr>
                  <tr>
                    <td className="pr-4">X1 cut-off density</td>
                    <td>
                      {scalars.x1CutoffDensityM3 === undefined
                        ? "— (evanescent)"
                        : `${scalars.x1CutoffDensityM3.toExponential(4)} m⁻³`}
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
