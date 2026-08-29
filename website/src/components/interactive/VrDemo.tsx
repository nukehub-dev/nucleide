import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { MagicSummary, SampledVoxelSummary } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Select } from "@nukehub/docs-kit/components/ui/Select";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";

const DEFAULT_MESHTAL = `mcnp version 5.mpi ld=00000000
Demo
Number of histories used for normalizing tallies = 1000

Mesh Tally Number 4
This is a neutron mesh tally.

Tally bin boundaries:
X direction: 0.0 1.0
Y direction: 0.0 1.0
Z direction: 0.0 1.0
Energy bin boundaries: 0.0 1.0

Energy X Y Z Result Rel Error
1.00000E+00 5.00000E-01 5.00000E-01 5.00000E-01 1.00000E+00 1.00000E-01`;

const SELECTION_OPTIONS = [
  { value: "total", label: "Total" },
  { value: "perGroup", label: "Per-group" },
];

const SAMPLE_MODE_OPTIONS = [
  { value: "analog", label: "Analog" },
  { value: "uniform", label: "Uniform" },
];

export function VrDemo() {
  const { wasm, ready, error } = useWasm();
  const [meshtalText, setMeshtalText] = useState(DEFAULT_MESHTAL);
  const [tallyNumber, setTallyNumber] = useState(4);
  const [selection, setSelection] = useState<"total" | "perGroup">("total");
  const [tolerance, setTolerance] = useState(0.5);
  const [nullValue, setNullValue] = useState(0.0);
  const [magic, setMagic] = useState<MagicSummary | null>(null);

  const [pdfInput, setPdfInput] = useState("1 3 2 1");
  const [aliasResult, setAliasResult] = useState<number | null>(null);

  const [sample, setSample] = useState<SampledVoxelSummary | null>(null);
  const [sampleMode, setSampleMode] = useState<"analog" | "uniform">("analog");

  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function runMagic() {
    if (!wasm) return;
    try {
      const out = wasm.magicBounds(meshtalText, tallyNumber, selection, tolerance, nullValue);
      setMagic(out);
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setMagic(null);
    }
  }

  function runAlias() {
    if (!wasm) return;
    try {
      const pdf = pdfInput.split(/\s+/).map(parseFloat).filter(Number.isFinite);
      const r1 = Math.random();
      const r2 = Math.random();
      setAliasResult(wasm.aliasTableSample(pdf, r1, r2));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setAliasResult(null);
    }
  }

  function runSample() {
    if (!wasm) return;
    try {
      const r1 = Math.random();
      const r2 = Math.random();
      setSample(wasm.meshSourceSample(meshtalText, tallyNumber, sampleMode, r1, r2));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setSample(null);
    }
  }

  const displayError = error ?? localError;

  return (
    <div className="rounded-xl border border-border/50 bg-background p-4 space-y-6">
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
            <Label>Mesh tally text</Label>
            <Textarea
              value={meshtalText}
              onChange={(e) => {
                setMeshtalText(e.target.value);
                clearError();
              }}
              className="font-mono text-xs"
            />
          </div>

          <div className="grid gap-3 sm:grid-cols-4">
            <div className="space-y-1">
              <Label>Tally</Label>
              <Input
                type="number"
                value={tallyNumber}
                onChange={(e) => {
                  setTallyNumber(parseInt(e.target.value));
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label>Selection</Label>
              <Select
                value={selection}
                onChange={(v) => {
                  setSelection(v as "total" | "perGroup");
                  clearError();
                }}
                options={SELECTION_OPTIONS}
                className="min-w-[140px]"
              />
            </div>
            <div className="space-y-1">
              <Label>Tolerance</Label>
              <Input
                type="number"
                step="0.1"
                value={tolerance}
                onChange={(e) => {
                  setTolerance(parseFloat(e.target.value));
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label>Null value</Label>
              <Input
                type="number"
                step="any"
                value={nullValue}
                onChange={(e) => {
                  setNullValue(parseFloat(e.target.value));
                  clearError();
                }}
              />
            </div>
          </div>

          <Button onClick={runMagic}>Generate MAGIC bounds</Button>

          {magic && (
            <div className="space-y-2 text-sm">
              <p>Tag: {magic.wwTagName}</p>
              <p>Groups per voxel: {magic.groupsPerVe}</p>
              <p>Lower bounds: [{magic.lowerBoundsWw.map((v) => v.toExponential(3)).join(", ")}]</p>
              <p>Scale factors: [{magic.scaleFactors.map((v) => v.toExponential(3)).join(", ")}]</p>
            </div>
          )}

          <div className="space-y-2 border-t border-border/50 pt-4">
            <Label>Alias-table sampling</Label>
            <Input
              type="text"
              value={pdfInput}
              onChange={(e) => {
                setPdfInput(e.target.value);
                clearError();
              }}
              className="font-mono"
            />
            <Button onClick={runAlias}>Sample index</Button>
            {aliasResult !== null && <p className="font-mono">Sampled index: {aliasResult}</p>}
          </div>

          <div className="space-y-2 border-t border-border/50 pt-4">
            <Label>Mesh source sampling</Label>
            <div className="flex flex-wrap items-end gap-2">
              <Select
                value={sampleMode}
                onChange={(v) => {
                  setSampleMode(v as "analog" | "uniform");
                  clearError();
                }}
                options={SAMPLE_MODE_OPTIONS}
                className="min-w-[140px]"
              />
              <Button onClick={runSample}>Sample voxel</Button>
            </div>
            {sample && (
              <p className="font-mono">
                index={sample.index} i={sample.i} j={sample.j} k={sample.k} weight=
                {sample.weight.toFixed(4)}
              </p>
            )}
          </div>
        </>
      )}
    </div>
  );
}
