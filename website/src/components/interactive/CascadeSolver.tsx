import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { CascadeResult } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";

export function CascadeSolver() {
  const { wasm, ready, error } = useWasm();
  const [config, setConfig] = useState({
    alpha: 1.05,
    Mstar: 236.5,
    N: 30,
    M: 10,
    feedAssay: 0.0072,
    productAssay: 0.05,
    tailsAssay: 0.0025,
  });
  const [result, setResult] = useState<CascadeResult | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function run(solveMode: "solve" | "multicomponent") {
    if (!wasm) return;
    try {
      const c = new wasm.WasmCascade({
        ...config,
        enrichingKey: "U235",
        strippingKey: "U238",
        feed: {
          U234: 0.000055,
          U235: config.feedAssay,
          U238: 1.0 - config.feedAssay - 0.000055,
        },
      });
      if (solveMode === "solve") {
        c.solve();
      } else {
        c.solveMulticomponent();
      }
      setResult(c.toObject());
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setResult(null);
    }
  }

  function update<K extends keyof typeof config>(key: K, value: typeof config[K]) {
    setConfig((prev) => ({ ...prev, [key]: value }));
    setLocalError(null);
  }

  const displayError = error ?? localError;

  return (
    <div className="rounded-xl border border-border/50 bg-background p-4 space-y-4">
      {!ready && <p className="text-sm text-muted-foreground">Loading Nucleide WASM…</p>}
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

      {ready && (
        <>
          <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-3">
            <NumberField label="Stage factor α" value={config.alpha} onChange={(v) => update("alpha", v)} />
            <NumberField label="M*" value={config.Mstar} onChange={(v) => update("Mstar", v)} />
            <NumberField label="Feed assay" value={config.feedAssay} onChange={(v) => update("feedAssay", v)} />
            <NumberField label="Product assay" value={config.productAssay} onChange={(v) => update("productAssay", v)} />
            <NumberField label="Tails assay" value={config.tailsAssay} onChange={(v) => update("tailsAssay", v)} />
            <NumberField label="N (guess)" value={config.N} onChange={(v) => update("N", v)} />
            <NumberField label="M (guess)" value={config.M} onChange={(v) => update("M", v)} />
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={() => run("solve")}>Solve cascade</Button>
            <Button variant="outline" onClick={() => run("multicomponent")}>
              Optimize M*
            </Button>
          </div>

          {result && (
            <div className="space-y-3">
              <table className="w-full text-sm">
                <tbody>
                  {[
                    ["Feed assay", result.feedAssay],
                    ["Product assay", result.productAssay],
                    ["Tails assay", result.tailsAssay],
                    ["Enriching stages", result.stagesEnriching],
                    ["Stripping stages", result.stagesStripping],
                    ["SWU / feed", result.swuPerFeed],
                    ["SWU / product", result.swuPerProduct],
                    ["Product / feed", result.productPerFeed],
                    ["Tails / feed", result.tailsPerFeed],
                  ].map(([label, value]) => (
                    <tr key={label} className="border-b border-border/50">
                      <td className="py-1">{label}</td>
                      <td className="py-1 text-right font-mono">
                        {typeof value === "number" ? value.toFixed(6) : value}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>

              <div className="grid gap-4 md:grid-cols-3">
                <StreamTable title="Feed" data={result.feed} />
                <StreamTable title="Product" data={result.product} />
                <StreamTable title="Tails" data={result.tails} />
              </div>
            </div>
          )}
        </>
      )}
    </div>
  );
}

function NumberField({
  label,
  value,
  onChange,
}: {
  label: string;
  value: number;
  onChange: (v: number) => void;
}) {
  return (
    <div className="space-y-1">
      <Label>{label}</Label>
      <Input
        type="number"
        step="any"
        value={value}
        onChange={(e) => onChange(parseFloat(e.target.value))}
      />
    </div>
  );
}

function StreamTable({ title, data }: { title: string; data: Record<string, number> }) {
  return (
    <div>
      <p className="text-sm font-medium">{title}</p>
      <table className="w-full text-sm">
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
