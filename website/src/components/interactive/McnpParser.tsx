import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type {
  McnpMaterialJson,
  MeshTallySummary,
  MeshtalSummary,
  WwinpSummary,
  XsdirSummary,
  XsdirTableJson,
} from "../../types/nucleide-wasm";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Select } from "@nukehub/docs-kit/components/ui/Select";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";
import { DataTable } from "@nukehub/docs-kit/components/mdx/DataTable";

const BASE = import.meta.env.BASE_URL.endsWith("/")
  ? import.meta.env.BASE_URL
  : `${import.meta.env.BASE_URL}/`;
const MESHTAL_SAMPLE_URL = `${BASE}data/meshtal_sample.txt`;
const XSDIR_SAMPLE_URL = `${BASE}data/xsdir_sample.txt`;

const BASE_XSDIR = `DATAPATH=/tmp
atomic weight ratios
1001 0.999167
directory
1001.00c 0.999167 h1 0 1 0 0`;

type ParserMode = "materials" | "xsdir" | "meshtal" | "wwinp";

const DEFAULTS: Record<ParserMode, string> = {
  materials: `c Test deck
m1 92235 -0.04 92238 -0.96 $ fuel
m2 1001 2 8016 1 $ water`,
  xsdir: BASE_XSDIR,
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

async function fetchSample(url: string): Promise<string> {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  return res.text();
}

export function McnpParser() {
  const { wasm, ready, error } = useWasm();
  const [mode, setMode] = useState<ParserMode>("materials");
  const [text, setText] = useState(DEFAULTS[mode]);
  const [materials, setMaterials] = useState<McnpMaterialJson[] | null>(null);
  const [xsdir, setXsdir] = useState<XsdirSummary | null>(null);
  const [meshtal, setMeshtal] = useState<MeshtalSummary | null>(null);
  const [wwinp, setWwinp] = useState<WwinpSummary | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);
  const [loadingSample, setLoadingSample] = useState(false);

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

  async function loadMeshtalSample() {
    setLoadingSample(true);
    try {
      const sample = await fetchSample(MESHTAL_SAMPLE_URL);
      setText(sample);
      setMode("meshtal");
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
    } finally {
      setLoadingSample(false);
    }
  }

  async function loadXsdirSample() {
    setLoadingSample(true);
    try {
      const sample = await fetchSample(XSDIR_SAMPLE_URL);
      setText(sample);
      setMode("xsdir");
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
    } finally {
      setLoadingSample(false);
    }
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
            className="font-mono text-xs"
          />

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Parse</Button>
            {mode === "meshtal" && (
              <Button variant="outline" onClick={loadMeshtalSample} disabled={loadingSample}>
                {loadingSample ? "Loading…" : "Load sample meshtal"}
              </Button>
            )}
            {mode === "xsdir" && (
              <Button variant="outline" onClick={loadXsdirSample} disabled={loadingSample}>
                {loadingSample ? "Loading…" : "Load sample xsdir"}
              </Button>
            )}
          </div>

          {materials && (
            <div className="space-y-4">
              {materials.map((m) => (
                <div key={m.number} className="space-y-2">
                  <div className="flex flex-wrap items-baseline gap-2 text-sm">
                    <span className="font-medium">Material {m.number}</span>
                    {m.density && (
                      <span className="text-xs text-muted-foreground">
                        Density: {m.density} g/cm³
                      </span>
                    )}
                    <span className="text-xs text-muted-foreground">Type: {m.fractionType}</span>
                  </div>
                  <DataTable
                    data={Object.entries(m.fractions).map(([nuclide, fraction]) => ({
                      nuclide,
                      fraction: fraction.toExponential(4),
                    }))}
                    columns={[
                      { key: "nuclide", header: "Nuclide" },
                      { key: "fraction", header: "Fraction", align: "right" },
                    ]}
                  />
                </div>
              ))}
            </div>
          )}

          {xsdir && (
            <div className="space-y-3">
              <div className="grid gap-2 text-sm sm:grid-cols-3">
                <p>Datapath: {xsdir.datapath ?? "none"}</p>
                <p>AWR entries: {xsdir.awrCount}</p>
                <p>Table count: {xsdir.tableCount}</p>
              </div>
              <DataTable
                data={xsdir.tables.map((t) => ({
                  name: <span className="font-mono">{t.name}</span>,
                  zaid: <span className="font-mono">{t.zaid}</span>,
                  awr: t.awr.toFixed(4),
                  filename: <span className="font-mono">{t.filename}</span>,
                }))}
                columns={[
                  { key: "name", header: "Name" },
                  { key: "zaid", header: "ZAID" },
                  { key: "awr", header: "AWR", align: "right" },
                  { key: "filename", header: "File" },
                ]}
                pagination
                pageSize={10}
              />
              <XsdirHistogram tables={xsdir.tables} />
            </div>
          )}

          {meshtal && (
            <div className="space-y-3">
              <div className="grid gap-2 text-sm sm:grid-cols-3">
                <p>Version: {meshtal.version}</p>
                <p>Title: {meshtal.title}</p>
                <p>Histories: {meshtal.histories.toExponential(4)}</p>
              </div>
              <DataTable
                data={Object.entries(meshtal.tallies).map(([num, t]) => ({
                  number: num,
                  particle: t.particle,
                  dims: t.dims.join("×"),
                  energyGroups: t.numEGroups,
                }))}
                columns={[
                  { key: "number", header: "Tally" },
                  { key: "particle", header: "Particle" },
                  { key: "dims", header: "Dimensions" },
                  { key: "energyGroups", header: "Energy groups", align: "right" },
                ]}
              />
              <MeshtalHeatmap tallies={meshtal.tallies} />
            </div>
          )}

          {wwinp && (
            <div className="space-y-2 text-sm">
              <p>
                ni: {wwinp.ni}, nr: {wwinp.nr}
              </p>
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

function MeshtalHeatmap({ tallies }: { tallies: Record<string, MeshTallySummary> }) {
  const tallyKeys = Object.keys(tallies);
  const [selectedKey, setSelectedKey] = useState<string>(tallyKeys[0] ?? "");
  const [field, setField] = useState<"result" | "relError">("result");

  if (tallyKeys.length === 0) return null;

  const tally = tallies[selectedKey];
  if (!tally) return null;

  const [nx, ny, nz] = tally.dims;
  if (nx * ny * nz === 0) return null;

  // Mid-slice through z.
  const k = Math.floor(nz / 2);
  const values: number[][] = [];
  for (let j = 0; j < ny; j++) {
    const row: number[] = [];
    for (let i = 0; i < nx; i++) {
      const ve = i + nx * (j + ny * k);
      row.push(field === "result" ? tally.totalResult[ve] : tally.totalRelError[ve]);
    }
    values.push(row);
  }

  const xMids = midpoints(tally.xBounds);
  const yMids = midpoints(tally.yBounds);
  const { z, tickvals, ticktext } = logTransform(values);

  return (
    <div className="space-y-2">
      <div className="flex flex-wrap items-end gap-4">
        <div className="w-28 space-y-1">
          <Label>Tally</Label>
          <Select
            value={selectedKey}
            onChange={setSelectedKey}
            options={tallyKeys.map((k) => ({ value: k, label: k }))}
          />
        </div>
        <div className="w-36 space-y-1">
          <Label>Field</Label>
          <Select
            value={field}
            onChange={(v) => setField(v as "result" | "relError")}
            options={[
              { value: "result", label: "Result" },
              { value: "relError", label: "Rel. error" },
            ]}
          />
        </div>
      </div>
      <Plotly
        aspect="video"
        data={[
          {
            type: "heatmap",
            x: xMids,
            y: yMids,
            z,
            zsmooth: false,
            colorscale: "Viridis",
            colorbar: {
              title: { text: field === "result" ? "Result (log₁₀)" : "Rel. error" },
              tickvals,
              ticktext,
            },
          },
        ]}
        layout={{
          xaxis: { title: { text: "x" } },
          yaxis: { title: { text: "y" } },
          margin: { t: 40, r: 16, b: 48, l: 48 },
          title: { text: `Mesh tally ${selectedKey}, z-slice ${k}` },
        }}
      />
    </div>
  );
}

function XsdirHistogram({ tables }: { tables: XsdirTableJson[] }) {
  if (tables.length === 0) return null;
  const byElement: Record<string, number> = {};
  for (const t of tables) {
    const z = parseInt(t.zaid.slice(0, -3), 10) || 0;
    const key = z > 0 ? elementSymbol(z) : "other";
    byElement[key] = (byElement[key] ?? 0) + 1;
  }
  const labels = Object.keys(byElement).sort((a, b) => byElement[b] - byElement[a]);
  const counts = labels.map((k) => byElement[k]);

  return (
    <Plotly
      aspect="video"
      data={[
        {
          type: "bar",
          x: labels,
          y: counts,
        },
      ]}
      layout={{
        xaxis: { title: { text: "Element" } },
        yaxis: { title: { text: "Table count" } },
        margin: { t: 16, r: 16, b: 48, l: 48 },
      }}
    />
  );
}

function elementSymbol(z: number): string {
  const symbols = [
    "H",
    "He",
    "Li",
    "Be",
    "B",
    "C",
    "N",
    "O",
    "F",
    "Ne",
    "Na",
    "Mg",
    "Al",
    "Si",
    "P",
    "S",
    "Cl",
    "Ar",
    "K",
    "Ca",
    "Sc",
    "Ti",
    "V",
    "Cr",
    "Mn",
    "Fe",
    "Co",
    "Ni",
    "Cu",
    "Zn",
    "Ga",
    "Ge",
    "As",
    "Se",
    "Br",
    "Kr",
    "Rb",
    "Sr",
    "Y",
    "Zr",
    "Nb",
    "Mo",
    "Tc",
    "Ru",
    "Rh",
    "Pd",
    "Ag",
    "Cd",
    "In",
    "Sn",
    "Sb",
    "Te",
    "I",
    "Xe",
    "Cs",
    "Ba",
    "La",
    "Ce",
    "Pr",
    "Nd",
    "Pm",
    "Sm",
    "Eu",
    "Gd",
    "Tb",
    "Dy",
    "Ho",
    "Er",
    "Tm",
    "Yb",
    "Lu",
    "Hf",
    "Ta",
    "W",
    "Re",
    "Os",
    "Ir",
    "Pt",
    "Au",
    "Hg",
    "Tl",
    "Pb",
    "Bi",
    "Po",
    "At",
    "Rn",
    "Fr",
    "Ra",
    "Ac",
    "Th",
    "Pa",
    "U",
    "Np",
    "Pu",
    "Am",
    "Cm",
    "Bk",
    "Cf",
    "Es",
    "Fm",
    "Md",
    "No",
    "Lr",
  ];
  return symbols[z - 1] ?? `Z${z}`;
}

function midpoints(bounds: number[]): number[] {
  const out: number[] = [];
  for (let i = 0; i < bounds.length - 1; i++) {
    out.push((bounds[i] + bounds[i + 1]) / 2);
  }
  return out;
}

function logTransform(z: number[][]): { z: number[][]; tickvals: number[]; ticktext: string[] } {
  const flat = z.flat();
  const positive = flat.filter((v) => v > 0);
  if (positive.length === 0) {
    return { z, tickvals: [], ticktext: [] };
  }
  const min = Math.min(...positive);
  const max = Math.max(...flat);
  const minExp = Math.floor(Math.log10(min));
  const maxExp = Math.ceil(Math.log10(max));
  const tickvals: number[] = [];
  const ticktext: string[] = [];
  for (let e = minExp; e <= maxExp; e++) {
    tickvals.push(e);
    ticktext.push(`1e${e}`);
  }
  const logMin = minExp - 1;
  return {
    z: z.map((row) => row.map((v) => (v > 0 ? Math.log10(v) : logMin))),
    tickvals,
    ticktext,
  };
}
