import { useEffect, useRef, useState } from "react";

// Load the wasm-pack generated module from the site root so it works both in
// dev and production. Astro's BASE_URL may or may not end with a slash depending
// on the command, so normalize it before appending the wasm path.
const BASE = import.meta.env.BASE_URL.endsWith("/")
  ? import.meta.env.BASE_URL
  : `${import.meta.env.BASE_URL}/`;
const WASM_URL = `${BASE}wasm/nucleide_wasm.js`;

interface WasmMaterial {
  atomFractions: () => Record<string, number>;
}

interface WasmCascade {
  solve: () => void;
  feedAssay: number;
  productAssay: number;
  tailsAssay: number;
  stagesEnriching: number;
  stagesStripping: number;
  swuPerFeed: number;
  swuPerProduct: number;
}

interface WasmModule {
  default: () => Promise<void>;
  WasmMaterial: new (formula: string) => WasmMaterial;
  WasmCascade: {
    defaultUranium: () => WasmCascade;
  };
}

interface WasmPlaygroundProps {
  kind: "material" | "cascade";
}

export function WasmPlayground({ kind }: WasmPlaygroundProps) {
  const wasmRef = useRef<WasmModule | null>(null);
  const [ready, setReady] = useState(false);
  const [error, setError] = useState<string | null>(null);

  // material state
  const [formula, setFormula] = useState("H2O");
  const [materialResult, setMaterialResult] = useState<Record<string, number> | null>(null);

  // cascade state
  const [cascadeResult, setCascadeResult] = useState<{
    feedAssay: number;
    productAssay: number;
    tailsAssay: number;
    stagesEnriching: number;
    stagesStripping: number;
    swuPerFeed: number;
    swuPerProduct: number;
  } | null>(null);

  useEffect(() => {
    let cancelled = false;

    async function load() {
      try {
        const mod = (await import(/* @vite-ignore */ WASM_URL)) as WasmModule;
        if (cancelled) return;
        await mod.default();
        if (cancelled) return;
        wasmRef.current = mod;
        setReady(true);
      } catch (e) {
        if (cancelled) return;
        setError(e instanceof Error ? e.message : String(e));
      }
    }

    load();
    return () => {
      cancelled = true;
    };
  }, []);

  function runMaterial() {
    const wasm = wasmRef.current;
    if (!wasm) return;
    try {
      const mat = new wasm.WasmMaterial(formula);
      const fractions = mat.atomFractions();
      setMaterialResult(fractions as Record<string, number>);
      setError(null);
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    }
  }

  function runCascade() {
    const wasm = wasmRef.current;
    if (!wasm) return;
    try {
      const cascade = wasm.WasmCascade.defaultUranium();
      cascade.solve();
      setCascadeResult({
        feedAssay: cascade.feedAssay,
        productAssay: cascade.productAssay,
        tailsAssay: cascade.tailsAssay,
        stagesEnriching: cascade.stagesEnriching,
        stagesStripping: cascade.stagesStripping,
        swuPerFeed: cascade.swuPerFeed,
        swuPerProduct: cascade.swuPerProduct,
      });
      setError(null);
    } catch (e) {
      setError(e instanceof Error ? e.message : String(e));
    }
  }

  if (error) {
    return (
      <div className="rounded-xl border border-red-200 bg-red-50 p-4 text-sm text-red-700 dark:border-red-900 dark:bg-red-950 dark:text-red-300">
        WASM error: {error}
      </div>
    );
  }

  return (
    <div className="rounded-xl border border-border/50 bg-background p-4">
      {!ready && <p className="text-sm text-muted-foreground">Loading Nucleide WASM…</p>}

      {ready && kind === "material" && (
        <div className="space-y-3">
          <div className="flex flex-wrap items-center gap-2">
            <label htmlFor="wasm-formula" className="text-sm font-medium">
              Formula
            </label>
            <input
              id="wasm-formula"
              type="text"
              value={formula}
              onChange={(e) => setFormula(e.target.value)}
              className="rounded-md border border-border bg-background px-2 py-1 text-sm"
            />
            <button
              onClick={runMaterial}
              className="rounded-md bg-primary px-3 py-1 text-sm font-medium text-primary-foreground"
            >
              Compute atom fractions
            </button>
          </div>

          {materialResult && (
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b border-border">
                  <th className="py-1 text-left">Nuclide</th>
                  <th className="py-1 text-right">Atom fraction</th>
                </tr>
              </thead>
              <tbody>
                {Object.entries(materialResult).map(([nuc, frac]) => (
                  <tr key={nuc} className="border-b border-border/50">
                    <td className="py-1 font-mono">{nuc}</td>
                    <td className="py-1 text-right font-mono">{frac.toExponential(4)}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          )}
        </div>
      )}

      {ready && kind === "cascade" && (
        <div className="space-y-3">
          <button
            onClick={runCascade}
            className="rounded-md bg-primary px-3 py-1 text-sm font-medium text-primary-foreground"
          >
            Solve default uranium cascade
          </button>

          {cascadeResult && (
            <table className="w-full text-sm">
              <tbody>
                {[
                  ["Feed assay", cascadeResult.feedAssay],
                  ["Product assay", cascadeResult.productAssay],
                  ["Tails assay", cascadeResult.tailsAssay],
                  ["Enriching stages", cascadeResult.stagesEnriching],
                  ["Stripping stages", cascadeResult.stagesStripping],
                  ["SWU / feed", cascadeResult.swuPerFeed],
                  ["SWU / product", cascadeResult.swuPerProduct],
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
          )}
        </div>
      )}
    </div>
  );
}
