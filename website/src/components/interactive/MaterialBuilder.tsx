import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Select } from "@nukehub/docs-kit/components/ui/Select";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";

const PATHWAY_OPTIONS = [
  { value: "air", label: "Air" },
  { value: "soil", label: "Soil" },
  { value: "ingest", label: "Ingest" },
  { value: "inhale", label: "Inhale" },
];

const SOURCE_OPTIONS = [
  { value: "EPA", label: "EPA" },
  { value: "DOE", label: "DOE" },
  { value: "GENII", label: "GENII" },
];

const DOSE_UNITS: Record<string, string> = {
  air: "mrem/h per g per m³",
  soil: "mrem/h per g per m²",
  ingest: "mrem per g",
  inhale: "mrem per g",
};

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
  const [pathway, setPathway] = useState("air");
  const [source, setSource] = useState("EPA");
  const [dose, setDose] = useState<number | null>(null);
  const [sepOut, setSepOut] = useState<string | null>(null);
  const [blendOut, setBlendOut] = useState<string | null>(null);
  const [cusumOut, setCusumOut] = useState<string | null>(null);
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

  function runDose() {
    if (!wasm) return;
    try {
      const mat = new wasm.WasmMaterial(formula);
      setDose(wasm.dosePerGram(mat.weightFractions(), pathway, source));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setDose(null);
    }
  }

  const displayError = error ?? localError;

  return (
    <div className="rounded-xl border border-border/50 bg-background p-4 space-y-4">
      {!ready && <p className="text-sm text-muted-foreground">Loading Nucleide WASM…</p>}

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
            {atomFracs && weightFracs && (
              <CompositionChart atomFracs={atomFracs} weightFracs={weightFracs} />
            )}
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
              // clip-path clips the painted scrollbar to the rounded shape too
              // (border-radius alone doesn't clip classic scrollbars). The 1px
              // expansion keeps clip AA from shaving the border's outer edge.
              <pre className="max-h-48 overflow-auto rounded-lg border border-border bg-muted p-2 text-xs [clip-path:inset(-1px_round_calc(0.5rem+1px))]">
                {xml}
              </pre>
            )}
          </div>

          <div className="space-y-2 border-t border-border/50 pt-4">
            <p className="text-sm font-medium">Separator, blender, and CUSUM</p>
            <div className="flex flex-wrap gap-2">
              <Button
                onClick={() => {
                  if (!wasm) return;
                  try {
                    const out = wasm.materialSeparate(
                      { U235: 7.2, U238: 992.8 },
                      { U235: 0.9, U238: 0.05 },
                    );
                    const p = out.product["U235"] ?? 0;
                    setSepOut(`product U235 = ${p.toFixed(2)} g (tails hold the rest)`);
                    clearError();
                  } catch (e) {
                    setLocalError(e instanceof Error ? e.message : String(e));
                  }
                }}
              >
                Separate
              </Button>
              <Button
                onClick={() => {
                  if (!wasm) return;
                  try {
                    const out = wasm.materialBlend([
                      { comp: { U235: 7.2, U238: 992.8 }, ratio: 1 },
                      { comp: { H1: 111.9, O16: 888.1 }, ratio: 2 },
                    ]);
                    setBlendOut(`blended nuclides: ${Object.keys(out).length}`);
                    clearError();
                  } catch (e) {
                    setLocalError(e instanceof Error ? e.message : String(e));
                  }
                }}
              >
                Blend
              </Button>
              <Button
                onClick={() => {
                  if (!wasm) return;
                  try {
                    const out = wasm.cusumDetect([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2]);
                    setCusumOut(
                      `CUSUM alarmed: ${out.alarmed ? "yes" : "no"}, statistic = ${out.statistic.toFixed(4)}`,
                    );
                    clearError();
                  } catch (e) {
                    setLocalError(e instanceof Error ? e.message : String(e));
                  }
                }}
              >
                Detect shift
              </Button>
            </div>
            {sepOut && <p className="text-sm">{sepOut}</p>}
            {blendOut && <p className="text-sm">{blendOut}</p>}
            {cusumOut && <p className="text-sm">{cusumOut}</p>}
          </div>

          <div className="space-y-2">
            <p className="text-sm font-medium">Dose per gram (screening only)</p>
            <p className="text-xs text-muted-foreground">
              Screening-level dose per gram of the current formula from the vendored HNF-5636 tables
              — not for safety decisions. Fully stable compositions such as H2O return 0 (no
              activity to convert). Compositions with nuclides outside the tables report an inline
              error instead of a value.
            </p>
            <div className="flex flex-wrap items-end gap-2">
              <div className="space-y-1 min-w-28">
                <Label>Pathway</Label>
                <Select
                  value={pathway}
                  onChange={(v) => {
                    setPathway(v);
                    // Air dose factors exist for EPA only (DOE/GENII rows are
                    // -1 sentinels); keep the source valid when switching.
                    if (v === "air") setSource("EPA");
                    clearError();
                  }}
                  options={PATHWAY_OPTIONS}
                />
              </div>
              <div className="space-y-1 min-w-32">
                <Label>Source</Label>
                <Select
                  value={source}
                  onChange={(v) => {
                    setSource(v);
                    clearError();
                  }}
                  options={SOURCE_OPTIONS.map((o) => ({
                    ...o,
                    disabled: pathway === "air" && o.value !== "EPA",
                  }))}
                />
              </div>
              <Button onClick={runDose}>Compute dose</Button>
            </div>
            {pathway === "air" && (
              <p className="text-xs text-muted-foreground">
                Air factors are EPA-only in the HNF tables, so the source is locked to EPA.
              </p>
            )}
            {dose !== null && (
              <p className="text-sm">
                Dose per gram: <span className="font-mono">{dose.toExponential(4)}</span>{" "}
                <span className="text-muted-foreground">{DOSE_UNITS[pathway]}</span>
              </p>
            )}
          </div>
        </>
      )}

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
    </div>
  );
}

function FractionTable({ title, data }: { title: string; data: Record<string, number> | null }) {
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

function CompositionChart({
  atomFracs,
  weightFracs,
}: {
  atomFracs: Record<string, number>;
  weightFracs: Record<string, number>;
}) {
  const keys = Array.from(new Set([...Object.keys(atomFracs), ...Object.keys(weightFracs)]));
  keys.sort((a, b) => {
    const maxA = Math.max(atomFracs[a] ?? 0, weightFracs[a] ?? 0);
    const maxB = Math.max(atomFracs[b] ?? 0, weightFracs[b] ?? 0);
    return maxB - maxA;
  });

  const atomValues = keys.map((k) => atomFracs[k] ?? 0);
  const weightValues = keys.map((k) => weightFracs[k] ?? 0);
  // Leave headroom on the value axis so the longest bar stays clear of the
  // modebar overlay in the top-right corner.
  const maxFraction = Math.max(0, ...atomValues, ...weightValues);
  const xMax = Math.min(1, maxFraction * 1.25);

  return (
    <Plotly
      aspect="video"
      data={[
        {
          type: "bar",
          orientation: "h",
          name: "Atom fraction",
          x: atomValues.slice().reverse(),
          y: keys.slice().reverse(),
        },
        {
          type: "bar",
          orientation: "h",
          name: "Weight fraction",
          x: weightValues.slice().reverse(),
          y: keys.slice().reverse(),
        },
      ]}
      layout={{
        barmode: "group",
        xaxis: { title: { text: "Fraction" }, range: [0, xMax] },
        yaxis: { title: { text: "Nuclide" } },
        margin: { t: 16, r: 16, b: 48, l: 64 },
        legend: { orientation: "h", y: -0.2 },
      }}
    />
  );
}
