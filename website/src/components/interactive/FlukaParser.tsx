import { useMemo, useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { UsrbinSummary, UsrbinTallyJson } from "../../types/nucleide-wasm";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Select } from "@nukehub/docs-kit/components/ui/Select";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";
import { DataTable } from "@nukehub/docs-kit/components/mdx/DataTable";
import { logTransform, midpoints } from "./McnpParser";

const BASE = import.meta.env.BASE_URL.endsWith("/")
  ? import.meta.env.BASE_URL
  : `${import.meta.env.BASE_URL}/`;
const USRBIN_SAMPLE_URL = `${BASE}data/usrbin_sample.lis`;

// One-bin Cartesian tally with the same record layout as real `.lis` output.
const DEFAULT_USRBIN = `1
   Cartesian binning n.   1  "mini    " , generalized particle n.   8
      X coordinate: from -1.0000E+00 to  1.0000E+00 cm,     1 bins ( 2.0000E+00 cm wide)
      Y coordinate: from -1.0000E+00 to  1.0000E+00 cm,     1 bins ( 2.0000E+00 cm wide)
      Z coordinate: from -1.0000E+00 to  1.0000E+00 cm,     1 bins ( 2.0000E+00 cm wide)
      Data follow in a matrix A(ix,iy,iz), format (1(5x,1p,10(1x,e11.4)))

       1.0000E+00

      Percentage errors follow in a matrix A(ix,iy,iz), format (1(5x,1p,10(1x,e11.4)))

       0.0000E+00`;

async function fetchSample(url: string): Promise<string> {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  return res.text();
}

export function FlukaParser() {
  const { wasm, ready, error } = useWasm();
  const [text, setText] = useState(DEFAULT_USRBIN);
  const [summary, setSummary] = useState<UsrbinSummary | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);
  const [loadingSample, setLoadingSample] = useState(false);

  function clearError() {
    setLocalError(null);
  }

  async function loadSample() {
    setLoadingSample(true);
    try {
      setText(await fetchSample(USRBIN_SAMPLE_URL));
      setSummary(null);
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
    } finally {
      setLoadingSample(false);
    }
  }

  function run() {
    if (!wasm) return;
    try {
      setSummary(wasm.parseUsrbin(text));
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
          <Textarea
            value={text}
            onChange={(e) => {
              setText(e.target.value);
              clearError();
            }}
            className="font-mono text-xs"
          />

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Parse</Button>
            <Button variant="outline" onClick={loadSample} disabled={loadingSample}>
              {loadingSample ? "Loading…" : "Load sample USRBIN"}
            </Button>
          </div>

          {summary && (
            <div className="space-y-3">
              <p className="text-sm">Tallies: {summary.tally_count}</p>
              <DataTable
                data={summary.tallies.map((t) => ({
                  name: <span className="font-mono">{t.name}</span>,
                  particle: t.particle,
                  coordSys: t.coord_sys,
                  dims: t.dims.join("×"),
                  maxRelError: percent(t),
                }))}
                columns={[
                  { key: "name", header: "Tally" },
                  { key: "particle", header: "Particle" },
                  { key: "coordSys", header: "Coord. system" },
                  { key: "dims", header: "Bins" },
                  { key: "maxRelError", header: "Max rel. error", align: "right" },
                ]}
              />
              <UsrbinHeatmap tallies={summary.tallies} />
            </div>
          )}
        </>
      )}
    </div>
  );
}

function percent(t: UsrbinTallyJson): string {
  // Explicit loop: spreading a large mesh into Math.max overflows the
  // browser's argument limit (~65k–125k args) on realistic USRBIN meshes.
  if (t.error_data.length === 0) return "";
  let max = -Infinity;
  for (const v of t.error_data) {
    if (v > max) max = v;
  }
  return `${max.toFixed(2)}%`;
}

function UsrbinHeatmap({ tallies }: { tallies: UsrbinTallyJson[] }) {
  const [selectedName, setSelectedName] = useState<string>(tallies[0]?.name ?? "");
  const [field, setField] = useState<"part_data" | "error_data">("part_data");

  const tally = tallies.find((t) => t.name === selectedName) ?? tallies[0];

  const slice = useMemo(() => {
    if (!tally) return null;
    const [nx, ny, nz] = tally.dims;
    if (nx * ny * nz === 0) return null;
    // FLUKA orders volume elements x slowest → z fastest:
    // ve = (i * ny + j) * nz + k. Take a mid-slice through z.
    const k = Math.floor(nz / 2);
    const data = field === "part_data" ? tally.part_data : tally.error_data;
    const values: number[][] = [];
    for (let i = 0; i < nx; i++) {
      const row: number[] = [];
      for (let j = 0; j < ny; j++) {
        row.push(data[(i * ny + j) * nz + k]);
      }
      values.push(row);
    }
    const xMids = midpoints(tally.x_bounds);
    const yMids = midpoints(tally.y_bounds);
    const { z, tickvals, ticktext } = logTransform(values);
    return { k, xMids, yMids, z, tickvals, ticktext };
  }, [tally, field]);

  if (!tally || !slice) return null;
  const { k, xMids, yMids, z, tickvals, ticktext } = slice;

  return (
    <div className="space-y-2">
      <div className="flex flex-wrap items-end gap-4">
        <div className="w-28 space-y-1">
          <Label>Tally</Label>
          <Select
            value={tally.name}
            onChange={setSelectedName}
            options={tallies.map((t) => ({ value: t.name, label: t.name }))}
          />
        </div>
        <div className="w-44 space-y-1">
          <Label>Field</Label>
          <Select
            value={field}
            onChange={(v) => setField(v as "part_data" | "error_data")}
            options={[
              { value: "part_data", label: "Track-length data" },
              { value: "error_data", label: "Percentage error" },
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
              title: {
                text: field === "part_data" ? "Tally (log₁₀)" : "Rel. error",
              },
              tickvals,
              ticktext,
            },
          },
        ]}
        layout={{
          xaxis: { title: { text: "x" } },
          yaxis: { title: { text: "y" } },
          margin: { t: 40, r: 16, b: 48, l: 48 },
          title: { text: `USRBIN ${tally.name}, z-slice ${k}` },
        }}
      />
    </div>
  );
}
