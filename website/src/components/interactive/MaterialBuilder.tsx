import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";

export function MaterialBuilder() {
  const { wasm, ready, error } = useWasm();
  const [formula, setFormula] = useState("H2O");
  const [mixInput, setMixInput] = useState("UO2 0.9\nH2O 0.1");
  const [xmlName, setXmlName] = useState("water");
  const [xmlDensity, setXmlDensity] = useState(1.0);
  const [atomFracs, setAtomFracs] = useState<Record<string, number> | null>(null);
  const [weightFracs, setWeightFracs] = useState<Record<string, number> | null>(null);
  const [mixedFracs, setMixedFracs] = useState<Record<string, number> | null>(null);
  const [xml, setXml] = useState<string | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function runFractions() {
    if (!wasm) return;
    try {
      const mat = new wasm.WasmMaterial(formula);
      setAtomFracs(mat.atomFractions());
      setWeightFracs(mat.weightFractions());
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setAtomFracs(null);
      setWeightFracs(null);
    }
  }

  function runMix() {
    if (!wasm) return;
    try {
      const parts = mixInput
        .split("\n")
        .map((line) => line.trim())
        .filter(Boolean)
        .map((line) => {
          const [f, r] = line.split(/\s+/);
          return { formula: f, fraction: parseFloat(r) };
        });
      const mixed = wasm.WasmMaterial.mixByMass(parts);
      setMixedFracs(mixed.atomFractions());
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setMixedFracs(null);
    }
  }

  function runXml() {
    if (!wasm) return;
    try {
      const mat = new wasm.WasmMaterial(formula);
      setXml(mat.toXml(xmlName, xmlDensity));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setXml(null);
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
          <div className="space-y-2">
            <div className="flex flex-wrap items-end gap-2">
              <div className="space-y-1">
                <Label htmlFor="mat-formula">Formula</Label>
                <Input
                  id="mat-formula"
                  type="text"
                  value={formula}
                  onChange={(e) => {
                    setFormula(e.target.value);
                    clearError();
                  }}
                />
              </div>
              <Button onClick={runFractions}>Compute fractions</Button>
            </div>
            <div className="grid gap-4 md:grid-cols-2">
              <FractionTable title="Atom fractions" data={atomFracs} />
              <FractionTable title="Weight fractions" data={weightFracs} />
            </div>
          </div>

          <div className="space-y-2">
            <Label>Mix by mass</Label>
            <p className="text-xs text-muted-foreground">
              One formula and relative mass fraction per line.
            </p>
            <Textarea
              value={mixInput}
              onChange={(e) => {
                setMixInput(e.target.value);
                clearError();
              }}
              className="font-mono text-xs"
            />
            <Button onClick={runMix}>Mix</Button>
            <FractionTable title="Mixed atom fractions" data={mixedFracs} />
          </div>

          <div className="space-y-2">
            <div className="flex flex-wrap items-end gap-2">
              <Label>XML export</Label>
              <Input
                type="text"
                value={xmlName}
                onChange={(e) => {
                  setXmlName(e.target.value);
                  clearError();
                }}
                placeholder="name"
              />
              <Input
                type="number"
                step="0.01"
                value={xmlDensity}
                onChange={(e) => {
                  setXmlDensity(parseFloat(e.target.value));
                  clearError();
                }}
                className="w-28"
              />
              <Button onClick={runXml}>To XML</Button>
            </div>
            {xml && (
              <pre className="max-h-48 overflow-auto rounded-lg border border-border bg-muted p-2 text-xs">
                {xml}
              </pre>
            )}
          </div>
        </>
      )}
    </div>
  );
}

function FractionTable({
  title,
  data,
}: {
  title: string;
  data: Record<string, number> | null;
}) {
  if (!data) return null;
  return (
    <div>
      <p className="text-sm font-medium">{title}</p>
      <table className="w-full text-sm">
        <thead>
          <tr className="border-b border-border">
            <th className="py-1 text-left">Nuclide</th>
            <th className="py-1 text-right">Fraction</th>
          </tr>
        </thead>
        <tbody>
          {Object.entries(data).map(([nuc, frac]) => (
            <tr key={nuc} className="border-b border-border/50">
              <td className="py-1 font-mono">{nuc}</td>
              <td className="py-1 text-right font-mono">{frac.toExponential(4)}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
}
