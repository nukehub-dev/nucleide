import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { IndataDoc, WallLoadResult } from "../../types/nucleide-wasm";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Label } from "@nukehub/docs-kit/components/ui/Label";

// Minimal STELLOPT-style namelist (synthetic values, grammar only).
const DEFAULT_INDATA = `&INDATA
  NFP = 5
  NCURR = 0
  LBOUND = .TRUE.
  AM(0) = 1.0D0
/`;

function formatValue(value: IndataDoc["scalars"][string]): string {
  if ("Int" in value) return `${value.Int}`;
  if ("Float" in value) return `${value.Float}`;
  if ("Bool" in value) return value.Bool ? ".TRUE." : ".FALSE.";
  return `'${value.Str}'`;
}

export function EquilibDemo() {
  const { wasm, ready, error } = useWasm();
  const [indataText, setIndataText] = useState(DEFAULT_INDATA);
  const [doc, setDoc] = useState<IndataDoc | null>(null);
  const [wall, setWall] = useState<WallLoadResult | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function run() {
    if (!wasm) return;
    try {
      setDoc(wasm.parseIndata(indataText));
      // Uniform unit-density field on a 2×3×2 grid: the wall-load
      // conservation shape in miniature (synthetic, no solver output).
      const sEdges = [0.0, 0.5, 1.0];
      const birth = new Array(2 * 3 * 2).fill(1.0);
      const jacobian = new Array(2 * 3 * 2).fill(1.0);
      setWall(wasm.wallLoad(sEdges, 3, 2, 5, birth, jacobian));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setDoc(null);
      setWall(null);
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
            Equilibrium data layer: parse a STELLOPT-style{" "}
            <span className="font-mono">&INDATA</span> block (grammar only — nothing is solved) and
            map a uniform birth-rate field onto the wall. Edit the namelist and re-run; malformed
            input is a loud error.
          </p>

          <div className="space-y-1">
            <Label className="flex min-h-10 items-end">INDATA text</Label>
            <textarea
              value={indataText}
              onChange={(e) => {
                setIndataText(e.target.value);
                clearError();
              }}
              rows={6}
              spellCheck={false}
              className="w-full rounded-lg border border-border/50 bg-muted/30 p-3 font-mono text-xs"
            />
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Parse equilibrium inputs</Button>
          </div>

          {doc && (
            <div className="space-y-2">
              <p className="text-sm font-medium">
                Parsed namelist ({Object.keys(doc.scalars).length} scalars,{" "}
                {Object.keys(doc.indexed).length} indexed)
              </p>
              <table className="text-sm">
                <tbody className="font-mono text-xs">
                  {Object.entries(doc.scalars).map(([name, value]) => (
                    <tr key={name}>
                      <td className="pr-4">{name}</td>
                      <td>{formatValue(value)}</td>
                    </tr>
                  ))}
                  {Object.entries(doc.indexed).map(([name, rows]) => (
                    <tr key={name}>
                      <td className="pr-4">
                        {name}({rows.map((r) => r.index.join(",")).join("; ")})
                      </td>
                      <td>{rows.map((r) => formatValue(r.value)).join("; ")}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}

          {wall && (
            <div className="space-y-2">
              <p className="text-sm font-medium">
                Wall-load map (uniform field, one-field-period total{" "}
                <span className="font-mono">{wall.total.toPrecision(6)}</span>)
              </p>
              <table className="text-sm">
                <tbody className="font-mono text-xs">
                  {wall.loads.map((row, j) => (
                    <tr key={j}>
                      <td className="pr-4">θ{j}</td>
                      {row.map((v, k) => (
                        <td key={k} className="pr-4">
                          {v.toPrecision(4)}
                        </td>
                      ))}
                    </tr>
                  ))}
                </tbody>
              </table>
              <p className="text-xs text-muted-foreground">
                Every wall node carries the same load here because the birth field and Jacobian are
                uniform; the total conserves the discrete births.
              </p>
              <Plotly
                aspect="video"
                data={[
                  {
                    type: "bar",
                    name: "wall load",
                    x: wall.loads.flatMap((row, j) => row.map((_, k) => `θ${j}ζ${k}`)),
                    y: wall.loads.flat(),
                  },
                ]}
                layout={{
                  xaxis: { title: { text: "Wall node" }, type: "category" },
                  yaxis: { title: { text: "Load" }, type: "linear" },
                  margin: { t: 16, r: 24, b: 48, l: 64 },
                  legend: { orientation: "h", y: -0.25 },
                }}
              />
            </div>
          )}
        </>
      )}
    </div>
  );
}
