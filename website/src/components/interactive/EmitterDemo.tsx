import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { DriftRowJson, EmitOpts } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Select } from "@nukehub/docs-kit/components/ui/Select";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";
import { DataTable } from "@nukehub/docs-kit/components/mdx/DataTable";

type KeyDialect = "gnds" | "armi";

const GNDS_DEFAULT = `{"U235": 5, "U238": 95}`;
const ARMI_DEFAULT = `{"nU235": 5, "nU238": 95}`;

const DEFAULTS: Record<KeyDialect, string> = {
  gnds: GNDS_DEFAULT,
  armi: ARMI_DEFAULT,
};

const DIALECTS: { value: KeyDialect; label: string }[] = [
  { value: "gnds", label: "GNDS keys" },
  { value: "armi", label: "ARMI keys" },
];

const PRESET_OPTIONS = [
  { value: "umetal-gnds", label: "U metal (GNDS)" },
  { value: "uox-gnds", label: "U oxide (GNDS)" },
  { value: "umetal-armi", label: "U metal (ARMI)" },
];

const PRESETS: Record<string, { dialect: KeyDialect; comp: string }> = {
  "umetal-gnds": { dialect: "gnds", comp: GNDS_DEFAULT },
  "uox-gnds": { dialect: "gnds", comp: `{"U235": 1.0, "U238": 20.0, "O16": 3.0}` },
  "umetal-armi": { dialect: "armi", comp: ARMI_DEFAULT },
};

export function EmitterDemo() {
  const { wasm, ready, error } = useWasm();
  const [dialect, setDialect] = useState<KeyDialect>("gnds");
  const [compText, setCompText] = useState(GNDS_DEFAULT);
  const [preset, setPreset] = useState("umetal-gnds");
  const [name, setName] = useState("umetal");
  const [density, setDensity] = useState("19.1");
  const [mcnpNumber, setMcnpNumber] = useState("1");
  const [xsSuffix, setXsSuffix] = useState("80c");
  const [serpentLib, setSerpentLib] = useState("03c");
  const [flukaFid, setFlukaFid] = useState("1");
  const [partisnZone, setPartisnZone] = useState("1");
  const [cards, setCards] = useState<Record<string, string> | null>(null);
  const [drift, setDrift] = useState<DriftRowJson[] | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function selectDialect(next: KeyDialect) {
    setDialect(next);
    setCompText(DEFAULTS[next]);
    setCards(null);
    setDrift(null);
    clearError();
  }

  function applyPreset(value: string) {
    const found = PRESETS[value];
    if (!found) return;
    setPreset(value);
    setDialect(found.dialect);
    setCompText(found.comp);
    setCards(null);
    setDrift(null);
    clearError();
  }

  function parseComp(): Record<string, number> {
    let parsed: unknown;
    try {
      parsed = JSON.parse(compText);
    } catch {
      throw new Error(`bad composition JSON (expected e.g. ${DEFAULTS[dialect]})`);
    }
    if (typeof parsed !== "object" || parsed === null || Array.isArray(parsed)) {
      throw new Error("composition must be a JSON object mapping nuclide keys to grams");
    }
    const comp: Record<string, number> = {};
    for (const [key, value] of Object.entries(parsed as Record<string, unknown>)) {
      if (typeof value !== "number" || !Number.isFinite(value) || value < 0) {
        throw new Error(
          `bad mass \`${String(value)}\` for \`${key}\` (expected finite grams >= 0)`,
        );
      }
      comp[key] = value;
    }
    if (Object.keys(comp).length === 0) throw new Error("composition is empty");
    return comp;
  }

  function parseName(): string {
    if (name.trim() === "") throw new Error("material name is empty");
    return name.trim();
  }

  function parseDensity(): number | undefined {
    if (density.trim() === "") return undefined;
    const v = parseFloat(density);
    if (!Number.isFinite(v) || v <= 0)
      throw new Error(`bad density \`${density}\` (expected finite g/cm³ > 0)`);
    return v;
  }

  function parseUint(raw: string, label: string): number | undefined {
    if (raw.trim() === "") return undefined;
    const v = parseInt(raw, 10);
    if (!Number.isInteger(v) || v < 0) throw new Error(`bad ${label} \`${raw}\``);
    return v;
  }

  function parseOpts(): EmitOpts {
    const opts: EmitOpts = {};
    const mcnp = parseUint(mcnpNumber, "MCNP number");
    if (mcnp !== undefined) opts.mcnpNumber = mcnp;
    if (xsSuffix.trim() !== "") opts.xsSuffix = xsSuffix.trim();
    if (serpentLib.trim() !== "") opts.serpentLib = serpentLib.trim();
    const fluka = parseUint(flukaFid, "FLUKA fid");
    if (fluka !== undefined) opts.flukaFid = fluka;
    const partisn = parseUint(partisnZone, "PARTISN zone");
    if (partisn !== undefined) opts.partisnZone = partisn;
    return opts;
  }

  function run() {
    if (!wasm) return;
    try {
      const comp = parseComp();
      const matName = parseName();
      const dens = parseDensity();
      const opts = parseOpts();
      if (dialect === "armi") {
        setCards(wasm.emitArmiCards(comp, matName, dens, opts));
        setDrift(wasm.emitArmiDriftTable(comp, matName, dens, opts));
      } else {
        setCards(wasm.emitCards(comp, matName, dens, opts));
        setDrift(wasm.emitDriftTable(comp, matName, dens, opts));
      }
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setCards(null);
      setDrift(null);
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
            {DIALECTS.map((d) => (
              <Button
                key={d.value}
                onClick={() => selectDialect(d.value)}
                variant={dialect === d.value ? "default" : "outline"}
                size="sm"
              >
                {d.label}
              </Button>
            ))}
          </div>

          <div className="max-w-xs space-y-1">
            <Label>Preset composition</Label>
            <Select value={preset} onChange={applyPreset} options={PRESET_OPTIONS} />
          </div>

          <div className="space-y-2">
            <Label>Composition (nuclide keys → grams, JSON)</Label>
            <Textarea
              value={compText}
              onChange={(e) => {
                setCompText(e.target.value);
                clearError();
              }}
              placeholder={DEFAULTS[dialect]}
              className="font-mono text-xs"
            />
          </div>

          <div className="grid gap-3 sm:grid-cols-3">
            <div className="space-y-1">
              <Label>Material name</Label>
              <Input
                value={name}
                onChange={(e) => {
                  setName(e.target.value);
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label>Density (g/cm³)</Label>
              <Input
                value={density}
                onChange={(e) => {
                  setDensity(e.target.value);
                  clearError();
                }}
                placeholder="Empty omits"
              />
            </div>
            <div className="space-y-1">
              <Label>MCNP number</Label>
              <Input
                value={mcnpNumber}
                onChange={(e) => {
                  setMcnpNumber(e.target.value);
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label>MCNP xs suffix</Label>
              <Input
                value={xsSuffix}
                onChange={(e) => {
                  setXsSuffix(e.target.value);
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label>Serpent lib</Label>
              <Input
                value={serpentLib}
                onChange={(e) => {
                  setSerpentLib(e.target.value);
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label>FLUKA fid</Label>
              <Input
                value={flukaFid}
                onChange={(e) => {
                  setFlukaFid(e.target.value);
                  clearError();
                }}
              />
            </div>
            <div className="space-y-1">
              <Label>PARTISN zone</Label>
              <Input
                value={partisnZone}
                onChange={(e) => {
                  setPartisnZone(e.target.value);
                  clearError();
                }}
              />
            </div>
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Emit</Button>
          </div>

          {cards && (
            <div className="space-y-3">
              {Object.entries(cards).map(([code, text]) => (
                <div key={code} className="space-y-1">
                  <p className="text-sm font-medium">{code}</p>
                  <pre className="overflow-x-auto rounded-lg border border-border/50 bg-muted/30 p-3 font-mono text-xs whitespace-pre-wrap">
                    {text}
                  </pre>
                </div>
              ))}
            </div>
          )}

          {drift && (
            <div className="space-y-2">
              <p className="text-sm font-medium">Mass-drift report</p>
              <DataTable
                data={drift.map((row) => ({
                  code: <span className="font-mono">{row.code}</span>,
                  massIn: row.massIn.toExponential(4),
                  massOut: row.massOut.toExponential(4),
                  relDrift: row.relDrift.toExponential(4),
                  dropped:
                    row.dropped.length > 0
                      ? row.dropped
                          .map((d) => `${d.nuclide} ${d.mass.toExponential(2)} (${d.reason})`)
                          .join("; ")
                      : "—",
                  reparsed: row.reparsed ? "yes" : "no",
                }))}
                columns={[
                  { key: "code", header: "Code" },
                  { key: "massIn", header: "Mass in (g)", align: "right" },
                  { key: "massOut", header: "Mass out (g)", align: "right" },
                  { key: "relDrift", header: "Rel. drift", align: "right" },
                  { key: "dropped", header: "Dropped" },
                  { key: "reparsed", header: "Reparsed" },
                ]}
              />
            </div>
          )}
        </>
      )}
    </div>
  );
}
