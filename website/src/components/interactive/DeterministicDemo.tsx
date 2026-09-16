import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { IsotxsSummary, PartisnDeck, RtfluxSummary } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";
import { DataTable } from "@nukehub/docs-kit/components/mdx/DataTable";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";
import { Select } from "@nukehub/docs-kit/components/ui/Select";

const DEFAULT_ISOTXS = `ISOTXS 2
NUCLIDE U235 92235 2
1.1 2.2
NUCLIDE PU239 94239 2
4.4 5.5`;

const DEFAULT_RTFLUX = `RTFLUX 2 3
1.0 2.0 3.0
4.0 5.0 6.0`;

interface PartisnZoneRow {
  id: string;
  material: string;
  density: string;
  labels: string;
}

const DEFAULT_ZONES: PartisnZoneRow[] = [
  { id: "1", material: "fuel", density: "10", labels: "U235" },
  { id: "2", material: "blanket", density: "5", labels: "PU239" },
];

export function DeterministicDemo() {
  const { wasm, ready, error } = useWasm();
  const [text, setText] = useState(DEFAULT_ISOTXS);
  const [summary, setSummary] = useState<IsotxsSummary | null>(null);
  const [rtfluxText, setRtfluxText] = useState(DEFAULT_RTFLUX);
  const [rtfluxKind, setRtfluxKind] = useState("rtflux");
  const [rtflux, setRtflux] = useState<RtfluxSummary | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  const [partisnTitle, setPartisnTitle] = useState("synthetic slab");
  const [partisnDim, setPartisnDim] = useState("1");
  const [partisnZones, setPartisnZones] = useState<PartisnZoneRow[]>(DEFAULT_ZONES);
  const [partisnSource, setPartisnSource] = useState("isotropic");
  const [partisnIsotxsText, setPartisnIsotxsText] = useState(DEFAULT_ISOTXS);
  const [partisnLabels, setPartisnLabels] = useState<string | null>(null);
  const [partisnPreview, setPartisnPreview] = useState<string | null>(null);
  const [partisnValid, setPartisnValid] = useState<string | null>(null);
  const [partisnInvalid, setPartisnInvalid] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function run() {
    if (!wasm) return;
    setSummary(null);
    try {
      setSummary(wasm.parseIsotxs(text));
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
    }
  }

  function buildPartisnDeck(): PartisnDeck {
    const title = partisnTitle.trim();
    if (title === "") throw new Error("PARTISN deck title is empty");
    const dim = parseInt(partisnDim, 10);
    if (![1, 2, 3].includes(dim))
      throw new Error(`bad PARTISN dim \`${partisnDim}\` (expected 1, 2, or 3)`);
    if (partisnZones.length === 0) throw new Error("PARTISN deck has no zones");
    const zones = partisnZones.map((row, i) => {
      const id = parseInt(row.id, 10);
      if (!Number.isInteger(id) || id < 0)
        throw new Error(`bad PARTISN zone ${i + 1} id \`${row.id}\``);
      if (row.material.trim() === "") throw new Error(`PARTISN zone ${i + 1} material is empty`);
      const density = parseFloat(row.density);
      if (!Number.isFinite(density))
        throw new Error(`bad PARTISN zone ${i + 1} density \`${row.density}\``);
      const isotxs_labels = row.labels
        .split(/[\s,]+/)
        .map((s) => s.trim())
        .filter((s) => s !== "");
      if (isotxs_labels.length === 0) throw new Error(`PARTISN zone ${i + 1} has no ISOTXS labels`);
      return { id, material: row.material.trim(), isotxs_labels, density };
    });
    const source = partisnSource.trim();
    return {
      title,
      dim,
      zones,
      source: source === "" ? null : source,
    };
  }

  function renderPartisn() {
    if (!wasm) return;
    setPartisnPreview(null);
    try {
      setPartisnPreview(wasm.partisnRender(buildPartisnDeck()));
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
    }
  }

  function validatePartisn() {
    if (!wasm) return;
    setPartisnValid(null);
    setPartisnInvalid(null);
    try {
      const deck = buildPartisnDeck();
      wasm.partisnValidate(deck, partisnIsotxsText);
      setPartisnValid(`PARTISN deck valid: ${deck.zones.length} zone(s), DIM ${deck.dim}`);
    } catch (e) {
      setPartisnInvalid(e instanceof Error ? e.message : String(e));
    }
  }

  function parsePartisnIsotxs() {
    if (!wasm) return;
    try {
      const lib = wasm.parseIsotxs(partisnIsotxsText);
      setPartisnLabels(lib.nuclides.map((n) => n.label).join(", "));
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setPartisnLabels(null);
    }
  }

  function setZone(index: number, patch: Partial<PartisnZoneRow>) {
    setPartisnZones((rows) => rows.map((row, i) => (i === index ? { ...row, ...patch } : row)));
  }

  const displayError = error ?? localError;

  // Per-point profile over the capped `values` (file order: `groups` values
  // per spatial point): one grouped-bar trace per energy group, following the
  // ISOTXS grouped-bar precedent above.
  const rtfluxProfile =
    rtflux && rtflux.groups > 0 && rtflux.values.length >= rtflux.groups
      ? (() => {
          const npoints = Math.floor(rtflux.values.length / rtflux.groups);
          const points = Array.from({ length: npoints }, (_, p) => `p${p + 1}`);
          const traces = Array.from({ length: rtflux.groups }, (_, g) => ({
            type: "bar" as const,
            name: `g${g + 1}`,
            x: points,
            y: Array.from({ length: npoints }, (_, p) => rtflux.values[p * rtflux.groups + g]),
          }));
          return { npoints, traces };
        })()
      : null;

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
          </div>

          {summary && (
            <div className="space-y-3">
              <p className="text-sm">
                Nuclides ({summary.nuclides.length}), energy groups: {summary.groups}
              </p>
              <DataTable
                data={summary.nuclides.map((n) => ({
                  label: <span className="font-mono">{n.label}</span>,
                  zaid: <span className="font-mono">{n.zaid}</span>,
                  groups: n.groups,
                  totalXs: n.total_xs.map((v) => v.toFixed(2)).join(", "),
                }))}
                columns={[
                  { key: "label", header: "Nuclide" },
                  { key: "zaid", header: "ZAID" },
                  { key: "groups", header: "Groups", align: "right" },
                  { key: "totalXs", header: "Total xs" },
                ]}
                pagination
                pageSize={10}
              />
              <Plotly
                aspect="video"
                data={summary.nuclides.map((n) => ({
                  type: "bar",
                  name: n.label,
                  x: n.total_xs.map((_, g) => `g${g + 1}`),
                  y: n.total_xs,
                }))}
                layout={{
                  barmode: "group",
                  xaxis: { title: { text: "Energy group" } },
                  yaxis: { title: { text: "Total xs" } },
                  margin: { t: 16, r: 16, b: 48, l: 64 },
                  legend: { orientation: "h", y: -0.25 },
                }}
              />
            </div>
          )}

          <div className="space-y-2 border-t border-border/50 pt-4">
            <p className="text-sm font-medium">RTFLUX fluxes (PARTISN preview below)</p>
            <Textarea
              value={rtfluxText}
              onChange={(e) => {
                setRtfluxText(e.target.value);
                clearError();
              }}
              className="font-mono text-xs"
            />
            <div className="flex flex-wrap items-end gap-2">
              <Select
                value={rtfluxKind}
                onChange={(v) => {
                  setRtfluxKind(v);
                  clearError();
                }}
                options={[
                  { value: "rtflux", label: "RTFLUX" },
                  { value: "atflux", label: "ATFLUX" },
                  { value: "rzflux", label: "RZFLUX" },
                ]}
              />
              <Button
                onClick={() => {
                  if (!wasm) return;
                  try {
                    setRtflux(wasm.parseRtflux(rtfluxText, rtfluxKind));
                    clearError();
                  } catch (e) {
                    setLocalError(e instanceof Error ? e.message : String(e));
                    setRtflux(null);
                  }
                }}
              >
                Parse RTFLUX
              </Button>
            </div>
            {rtflux && (
              <div className="space-y-3">
                <p className="text-sm">
                  Flux kind: <span className="font-mono">{rtflux.kind}</span>, groups:{" "}
                  {rtflux.groups}, points: {rtflux.npoints}
                </p>
                {rtfluxProfile && (
                  <>
                    <p className="text-sm">
                      Flux profile across {rtfluxProfile.npoints} points × {rtflux.groups} groups
                    </p>
                    <Plotly
                      aspect="video"
                      data={rtfluxProfile.traces}
                      layout={{
                        barmode: "group",
                        xaxis: { title: { text: "Spatial point" } },
                        yaxis: { title: { text: "Flux" } },
                        margin: { t: 16, r: 16, b: 48, l: 64 },
                        legend: { orientation: "h", y: -0.25 },
                      }}
                    />
                    {rtflux.truncated && (
                      <p className="text-xs text-muted-foreground">
                        Values capped at the {rtflux.values.length}-value demo cap; the chart shows
                        the leading points only.
                      </p>
                    )}
                  </>
                )}
              </div>
            )}
          </div>

          <div className="space-y-3 border-t border-border/50 pt-4">
            <p className="text-sm font-medium">PARTISN writer (preview below)</p>
            <p className="text-xs text-muted-foreground">
              Structured deck dict with exact keys: <span className="font-mono">title</span>,{" "}
              <span className="font-mono">dim</span>, <span className="font-mono">zones</span> (each{" "}
              <span className="font-mono">id</span>, <span className="font-mono">material</span>,{" "}
              <span className="font-mono">isotxs_labels</span>,{" "}
              <span className="font-mono">density</span>), optional{" "}
              <span className="font-mono">source</span> — mirroring the Python{" "}
              <span className="font-mono">partisn_render</span>/
              <span className="font-mono">partisn_validate</span> shape. Labels validate against the
              pasted ISOTXS library below.
            </p>

            <div className="grid gap-3 sm:grid-cols-3">
              <div className="space-y-1">
                <Label htmlFor="partisn-title">Deck title</Label>
                <Input
                  id="partisn-title"
                  value={partisnTitle}
                  onChange={(e) => {
                    setPartisnTitle(e.target.value);
                    clearError();
                  }}
                />
              </div>
              <div className="space-y-1">
                <Label>Dimension</Label>
                <Select
                  value={partisnDim}
                  onChange={(v) => {
                    setPartisnDim(v);
                    clearError();
                  }}
                  options={[
                    { value: "1", label: "1" },
                    { value: "2", label: "2" },
                    { value: "3", label: "3" },
                  ]}
                />
              </div>
              <div className="space-y-1">
                <Label htmlFor="partisn-source">Source (optional)</Label>
                <Input
                  id="partisn-source"
                  value={partisnSource}
                  onChange={(e) => {
                    setPartisnSource(e.target.value);
                    clearError();
                  }}
                  placeholder="Empty omits SOURCE"
                />
              </div>
            </div>

            <div className="space-y-3">
              {partisnZones.map((zone, i) => (
                <div key={i} className="grid gap-3 sm:grid-cols-4">
                  <div className="space-y-1">
                    <Label htmlFor={`partisn-zone-${i}-id`}>Zone {i + 1} id</Label>
                    <Input
                      id={`partisn-zone-${i}-id`}
                      value={zone.id}
                      onChange={(e) => {
                        setZone(i, { id: e.target.value });
                        clearError();
                      }}
                    />
                  </div>
                  <div className="space-y-1">
                    <Label htmlFor={`partisn-zone-${i}-material`}>Zone {i + 1} material</Label>
                    <Input
                      id={`partisn-zone-${i}-material`}
                      value={zone.material}
                      onChange={(e) => {
                        setZone(i, { material: e.target.value });
                        clearError();
                      }}
                    />
                  </div>
                  <div className="space-y-1">
                    <Label htmlFor={`partisn-zone-${i}-density`}>Zone {i + 1} density</Label>
                    <Input
                      id={`partisn-zone-${i}-density`}
                      value={zone.density}
                      onChange={(e) => {
                        setZone(i, { density: e.target.value });
                        clearError();
                      }}
                    />
                  </div>
                  <div className="space-y-1">
                    <Label htmlFor={`partisn-zone-${i}-labels`}>Zone {i + 1} labels</Label>
                    <Input
                      id={`partisn-zone-${i}-labels`}
                      value={zone.labels}
                      onChange={(e) => {
                        setZone(i, { labels: e.target.value });
                        clearError();
                      }}
                      placeholder="space/comma separated"
                    />
                  </div>
                </div>
              ))}
              <div className="flex flex-wrap gap-2">
                <Button
                  variant="outline"
                  size="sm"
                  onClick={() => {
                    setPartisnZones((rows) => [
                      ...rows,
                      {
                        id: String(rows.length + 1),
                        material: "shield",
                        density: "1",
                        labels: "U235",
                      },
                    ]);
                    clearError();
                  }}
                >
                  Add zone
                </Button>
                {partisnZones.length > 1 && (
                  <Button
                    variant="outline"
                    size="sm"
                    onClick={() => {
                      setPartisnZones((rows) => rows.slice(0, -1));
                      clearError();
                    }}
                  >
                    Remove zone
                  </Button>
                )}
              </div>
            </div>

            <div className="space-y-2">
              <Label htmlFor="partisn-isotxs">ISOTXS library (pasted)</Label>
              <Textarea
                id="partisn-isotxs"
                value={partisnIsotxsText}
                onChange={(e) => {
                  setPartisnIsotxsText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
              <div className="flex flex-wrap gap-2">
                <Button variant="outline" size="sm" onClick={parsePartisnIsotxs}>
                  Parse PARTISN ISOTXS
                </Button>
              </div>
              {partisnLabels && (
                <p className="text-sm">
                  ISOTXS labels: <span className="font-mono">{partisnLabels}</span>
                </p>
              )}
            </div>

            <div className="flex flex-wrap gap-2">
              <Button onClick={renderPartisn}>Render PARTISN</Button>
              <Button onClick={validatePartisn}>Validate PARTISN</Button>
            </div>

            {partisnPreview && (
              <div className="space-y-1">
                <p className="text-sm font-medium">PARTISN preview</p>
                <pre className="overflow-x-auto rounded-lg border border-border/50 bg-muted/30 p-3 font-mono text-xs whitespace-pre-wrap">
                  {partisnPreview}
                </pre>
              </div>
            )}
            {partisnValid && (
              <p className="text-sm text-green-700 dark:text-green-300">{partisnValid}</p>
            )}
            {partisnInvalid && (
              <p className="text-sm text-red-700 dark:text-red-300">
                PARTISN validation error: {partisnInvalid}
              </p>
            )}
          </div>
        </>
      )}
    </div>
  );
}
