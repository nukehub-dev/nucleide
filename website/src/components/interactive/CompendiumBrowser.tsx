import { useEffect, useRef, useState } from "react";
import { useWasm } from "../../lib/wasm";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";
import { Combobox } from "@nukehub/docs-kit/components/ui/Combobox";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import type { CompendiumEntryInfo, WasmMaterialsCompendium } from "../../types/nucleide-wasm";

const BASE = import.meta.env.BASE_URL.endsWith("/")
  ? import.meta.env.BASE_URL
  : `${import.meta.env.BASE_URL}/`;
const DATA_URL = `${BASE}data/MaterialsCompendium.json`;
const DEFAULT_MATERIAL = "Water, Liquid";
const TOP_N = 12;

export function CompendiumBrowser() {
  const { wasm, ready, error } = useWasm();
  const libRef = useRef<WasmMaterialsCompendium | null>(null);
  const [names, setNames] = useState<string[]>([]);
  const [count, setCount] = useState<number | null>(null);
  const [selected, setSelected] = useState<string | null>(null);
  const [entry, setEntry] = useState<CompendiumEntryInfo | null>(null);
  const [loadingData, setLoadingData] = useState(false);
  const [localError, setLocalError] = useState<string | null>(null);

  // Fetch the compendium JSON and parse it with the WASM build once ready.
  useEffect(() => {
    if (!ready || !wasm) return;
    let cancelled = false;
    setLoadingData(true);

    fetch(DATA_URL)
      .then((r) => {
        if (!r.ok) throw new Error(`failed to fetch compendium (HTTP ${r.status})`);
        return r.text();
      })
      .then((text) => {
        if (cancelled) return;
        const lib = wasm.WasmMaterialsCompendium.fromJson(text);
        const all = lib.names();
        libRef.current = lib;
        setNames(all);
        setCount(lib.len);
        setSelected(all.includes(DEFAULT_MATERIAL) ? DEFAULT_MATERIAL : (all[0] ?? null));
        setLoadingData(false);
      })
      .catch((e) => {
        if (cancelled) return;
        setLocalError(e instanceof Error ? e.message : String(e));
        setLoadingData(false);
      });

    return () => {
      cancelled = true;
    };
  }, [ready, wasm]);

  useEffect(() => {
    if (!selected || !libRef.current) return;
    try {
      setEntry(libRef.current.get(selected));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setEntry(null);
    }
  }, [selected]);

  const displayError = error ?? localError;

  const fractions = entry ? Object.entries(entry.weight_fractions) : [];
  fractions.sort((a, b) => b[1] - a[1]);
  const top = fractions.slice(0, TOP_N);
  const rest = fractions.slice(TOP_N);
  const restTotal = rest.reduce((sum, [, v]) => sum + v, 0);
  const chartLabels = top.map(([nuc]) => nuc);
  const chartValues = top.map(([, v]) => v);
  if (rest.length > 0) {
    chartLabels.push(`other (${rest.length})`);
    chartValues.push(restTotal);
  }

  return (
    <div className="rounded-xl border border-border/50 bg-background p-4 space-y-4">
      {!ready && <p className="text-sm text-muted-foreground">Loading Nucleide WASM…</p>}
      {ready && loadingData && (
        <p className="text-sm text-muted-foreground">Loading Materials Compendium…</p>
      )}

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

      {ready && count !== null && (
        <>
          <div className="space-y-1">
            <Label>Material ({count} materials loaded)</Label>
            <Combobox
              value={selected ?? undefined}
              onChange={setSelected}
              options={names.map((n) => ({ value: n, label: n }))}
              placeholder="Select a material"
              searchPlaceholder="Search materials…"
            />
          </div>

          {entry && (
            <>
              <p className="text-sm text-muted-foreground">
                MatNum {entry.mat_num} · density {entry.density} g/cm³
                {entry.acronym.length > 0 && <> · {entry.acronym.join(", ")}</>}
                {entry.source && <> · source: {entry.source}</>}
              </p>

              {chartValues.length > 0 && (
                <Plotly
                  aspect="video"
                  data={[
                    {
                      type: "bar",
                      orientation: "h",
                      x: chartValues.slice().reverse(),
                      y: chartLabels.slice().reverse(),
                    },
                  ]}
                  layout={{
                    xaxis: { title: { text: "Weight fraction" } },
                    margin: { t: 16, r: 16, b: 48, l: 96 },
                  }}
                />
              )}

              <div className="max-h-64 overflow-auto">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="border-b border-border">
                      <th className="py-1 text-left">Nuclide</th>
                      <th className="py-1 text-right">Weight fraction</th>
                    </tr>
                  </thead>
                  <tbody>
                    {fractions.map(([nuc, frac]) => (
                      <tr key={nuc} className="border-b border-border/50">
                        <td className="py-1 font-mono">{nuc}</td>
                        <td className="py-1 text-right font-mono">{frac.toExponential(4)}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </>
          )}
        </>
      )}
    </div>
  );
}
