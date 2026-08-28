import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type {
  McnpMaterialJson,
  MeshtalSummary,
  WwinpSummary,
  XsdirSummary,
} from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";

type ParserMode = "materials" | "xsdir" | "meshtal" | "wwinp";

const DEFAULTS: Record<ParserMode, string> = {
  materials: `c Test deck
m1 92235 -0.04 92238 -0.96 $ fuel
m2 1001 2 8016 1 $ water`,
  xsdir: `DATAPATH=/tmp
atomic weight ratios
1001 0.999167
directory
1001.00c 0.999167 h1 0 1 0 0`,
  meshtal: `mcnp version 5.mpi ld=00000000
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
1.00000E+00 5.00000E-01 5.00000E-01 5.00000E-01 1.00000E+00 1.00000E-01`,
  wwinp: `1 1 1 10
7
1 1 1 0 0 0
1 1 1 1
0 1 1 1
0 1 1 1
0 1 1 1
0.1 0.2 0.5 1 2 5 10
1
2
3
4
5
6
7`,
};

export function McnpParser() {
  const { wasm, ready, error } = useWasm();
  const [mode, setMode] = useState<ParserMode>("materials");
  const [text, setText] = useState(DEFAULTS[mode]);
  const [materials, setMaterials] = useState<McnpMaterialJson[] | null>(null);
  const [xsdir, setXsdir] = useState<XsdirSummary | null>(null);
  const [meshtal, setMeshtal] = useState<MeshtalSummary | null>(null);
  const [wwinp, setWwinp] = useState<WwinpSummary | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function selectMode(next: ParserMode) {
    setMode(next);
    setText(DEFAULTS[next]);
    setMaterials(null);
    setXsdir(null);
    setMeshtal(null);
    setWwinp(null);
    clearError();
  }

  function updateText(next: string) {
    setText(next);
    clearError();
  }

  function run() {
    if (!wasm) return;
    setMaterials(null);
    setXsdir(null);
    setMeshtal(null);
    setWwinp(null);
    try {
      switch (mode) {
        case "materials":
          setMaterials(wasm.parseMcnpMaterials(text));
          break;
        case "xsdir":
          setXsdir(wasm.parseXsdir(text));
          break;
        case "meshtal":
          setMeshtal(wasm.parseMeshtal(text));
          break;
        case "wwinp":
          setWwinp(wasm.parseWwinp(text));
          break;
      }
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
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
          <div className="flex flex-wrap gap-2">
            {(["materials", "xsdir", "meshtal", "wwinp"] as ParserMode[]).map((m) => (
              <Button
                key={m}
                onClick={() => selectMode(m)}
                variant={mode === m ? "default" : "outline"}
                size="sm"
              >
                {m}
              </Button>
            ))}
          </div>

          <Textarea
            value={text}
            onChange={(e) => updateText(e.target.value)}
            rows={8}
            autoResize={false}
            className="font-mono text-xs"
          />

          <Button onClick={run}>Parse</Button>

          {materials && (
            <div className="space-y-3">
              {materials.map((m) => (
                <div key={m.number} className="rounded-md border border-border/50 p-3">
                  <p className="text-sm font-medium">Material {m.number}</p>
                  {m.density && <p className="text-xs text-muted-foreground">Density: {m.density} g/cm³</p>}
                  <p className="text-xs text-muted-foreground">Type: {m.fractionType}</p>
                  <table className="w-full text-sm">
                    <tbody>
                      {Object.entries(m.fractions).map(([nuc, frac]) => (
                        <tr key={nuc} className="border-b border-border/50">
                          <td className="py-1 font-mono">{nuc}</td>
                          <td className="py-1 text-right font-mono">{frac.toExponential(4)}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              ))}
            </div>
          )}

          {xsdir && (
            <div className="space-y-2 text-sm">
              <p>Datapath: {xsdir.datapath ?? "none"}</p>
              <p>AWR entries: {xsdir.awrCount}</p>
              <p>Table count: {xsdir.tableCount}</p>
              <div className="max-h-64 overflow-auto rounded-md border border-border">
                <table className="w-full text-sm">
                  <thead className="bg-muted">
                    <tr>
                      <th className="px-2 py-1 text-left">Name</th>
                      <th className="px-2 py-1 text-left">ZAID</th>
                      <th className="px-2 py-1 text-right">AWR</th>
                      <th className="px-2 py-1 text-left">File</th>
                    </tr>
                  </thead>
                  <tbody>
                    {xsdir.tables.slice(0, 50).map((t) => (
                      <tr key={t.name} className="border-b border-border/50">
                        <td className="px-2 py-1 font-mono">{t.name}</td>
                        <td className="px-2 py-1 font-mono">{t.zaid}</td>
                        <td className="px-2 py-1 text-right">{t.awr.toFixed(4)}</td>
                        <td className="px-2 py-1">{t.filename}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          )}

          {meshtal && (
            <div className="space-y-2 text-sm">
              <p>Version: {meshtal.version}</p>
              <p>Title: {meshtal.title}</p>
              <p>Histories: {meshtal.histories.toExponential(4)}</p>
              <p>Tallies: {meshtal.tallyCount}</p>
              {Object.entries(meshtal.tallies).map(([num, t]) => (
                <div key={num} className="rounded-md border border-border/50 p-3">
                  <p className="font-medium">Tally {num}</p>
                  <p className="text-xs text-muted-foreground">
                    Particle: {t.particle}, dims: {t.dims.join("×")}, energy groups: {t.numEGroups}
                  </p>
                </div>
              ))}
            </div>
          )}

          {wwinp && (
            <div className="space-y-2 text-sm">
              <p>ni: {wwinp.ni}, nr: {wwinp.nr}</p>
              <p>nf: [{wwinp.nf.join(", ")}]</p>
              <p>nc: [{wwinp.nc.join(", ")}]</p>
              <p>origin: [{wwinp.origin.join(", ")}]</p>
              <p>Energy groups: {wwinp.ne.join(", ")}</p>
            </div>
          )}
        </>
      )}
    </div>
  );
}
