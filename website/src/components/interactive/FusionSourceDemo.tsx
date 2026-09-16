import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type {
  FusionCardsResult,
  FusionRingSpec,
  FusionSampleResult,
  FusionSpectrumMoments,
} from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Plotly } from "@nukehub/docs-kit/components/mdx/PlotlyClient";

// Tokamak midplane ring at R = 300 cm, T_i = 20 keV D-T: the crate's example
// configuration (synthetic source geometry, no evaluated data).
const DEFAULT_RADIUS = "300";
const DEFAULT_HEIGHT = "0";
const DEFAULT_TI = "20";
const DEFAULT_N = "2000";
const DEFAULT_SEED = "42";
const MAX_SAMPLES = 5000;

type Reaction = "dt" | "dd";

const REACTIONS: { value: Reaction; label: string }[] = [
  { value: "dt", label: "D-T (14.1 MeV)" },
  { value: "dd", label: "D-D (2.45 MeV)" },
];

export function FusionSourceDemo() {
  const { wasm, ready, error } = useWasm();
  const [reaction, setReaction] = useState<Reaction>("dt");
  const [radiusText, setRadiusText] = useState(DEFAULT_RADIUS);
  const [heightText, setHeightText] = useState(DEFAULT_HEIGHT);
  const [tiText, setTiText] = useState(DEFAULT_TI);
  const [nText, setNText] = useState(DEFAULT_N);
  const [seedText, setSeedText] = useState(DEFAULT_SEED);
  const [sample, setSample] = useState<FusionSampleResult | null>(null);
  const [moments, setMoments] = useState<FusionSpectrumMoments | null>(null);
  const [cards, setCards] = useState<FusionCardsResult | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);

  function clearError() {
    setLocalError(null);
  }

  function parseEntry(text: string, label: string): number {
    const v = parseFloat(text);
    if (!Number.isFinite(v)) throw new Error(`bad ${label} \`${text}\``);
    return v;
  }

  function buildSpec(): FusionRingSpec {
    const n = parseEntry(nText, "n");
    if (!Number.isInteger(n) || n < 1 || n > MAX_SAMPLES)
      throw new Error(`bad n \`${nText}\` (expected an integer in [1, ${MAX_SAMPLES}])`);
    const seed = parseEntry(seedText, "seed");
    if (!Number.isInteger(seed) || seed < 0 || seed >= 2 ** 53)
      throw new Error(`bad seed \`${seedText}\` (expected a non-negative integer below 2^53)`);
    const radiusCm = parseEntry(radiusText, "radius");
    if (radiusCm <= 0) throw new Error(`bad radius \`${radiusText}\` (expected > 0)`);
    const heightCm = parseEntry(heightText, "height");
    const tiKev = parseEntry(tiText, "ion temperature");
    if (tiKev < 0) throw new Error(`bad ion temperature \`${tiText}\` (expected >= 0)`);
    return {
      kind: "ring",
      radiusCm,
      heightCm,
      reaction,
      tiKev,
      n,
      seed,
    };
  }

  function run() {
    if (!wasm) return;
    try {
      const spec = buildSpec();
      const out = wasm.sampleFusionSource(spec);
      setSample(out);
      setMoments(wasm.fusionSpectrumMoments(spec.reaction, spec.tiKev));
      setCards(wasm.emitFusionSourceCards({ ...spec, nBins: 21 }));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setSample(null);
      setMoments(null);
      setCards(null);
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
          <div className="flex flex-wrap gap-2" role="group" aria-label="Fusion reaction">
            {REACTIONS.map((r) => (
              <Button
                key={r.value}
                onClick={() => {
                  setReaction(r.value);
                  setSample(null);
                  setMoments(null);
                  setCards(null);
                  clearError();
                }}
                variant={reaction === r.value ? "default" : "outline"}
                size="sm"
              >
                {r.label}
              </Button>
            ))}
          </div>

          <div className="grid gap-3 sm:grid-cols-2 lg:grid-cols-5">
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Ring radius [cm]</Label>
              <Input
                value={radiusText}
                onChange={(e) => {
                  setRadiusText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Ring height [cm]</Label>
              <Input
                value={heightText}
                onChange={(e) => {
                  setHeightText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Ion temperature [keV]</Label>
              <Input
                value={tiText}
                onChange={(e) => {
                  setTiText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Particles n (max {MAX_SAMPLES})</Label>
              <Input
                value={nText}
                onChange={(e) => {
                  setNText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Seed</Label>
              <Input
                value={seedText}
                onChange={(e) => {
                  setSeedText(e.target.value);
                  clearError();
                }}
                className="font-mono text-xs"
              />
            </div>
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Sample source</Button>
          </div>

          {moments && (
            <p className="text-sm">
              Spectrum: <span className="font-mono">{moments.reaction}</span> line at{" "}
              <span className="font-mono">{moments.nominalMeV.toFixed(4)} MeV</span>, Ballabio mean{" "}
              <span className="font-mono">{moments.meanMeV.toFixed(4)} MeV</span>, σ{" "}
              <span className="font-mono">
                {moments.mono ? "0 (monoenergetic)" : `${moments.sigmaMeV.toExponential(4)} MeV`}
              </span>
            </p>
          )}

          {sample && (
            <div className="space-y-2">
              <p className="text-sm font-medium">
                Sampled birth energies ({sample.count} particles)
              </p>
              <Plotly
                aspect="video"
                data={[
                  {
                    type: "histogram",
                    name: "energy",
                    x: sample.particles.map((p) => p.energyMeV),
                    nbinsx: 60,
                  },
                ]}
                layout={{
                  xaxis: { title: { text: "Birth energy (MeV)" }, type: "linear" },
                  yaxis: { title: { text: "Count" }, type: "linear" },
                  margin: { t: 16, r: 24, b: 48, l: 64 },
                  legend: { orientation: "h", y: -0.25 },
                }}
              />
              <p className="text-xs text-muted-foreground">
                Birth positions sit on the ring (fixed radius and height, uniform azimuth);
                directions are isotropic. Identical inputs reproduce the identical stream (seeded).
              </p>
            </div>
          )}

          {cards && (
            <div className="space-y-2">
              <p className="text-sm font-medium">Emitted source cards</p>
              <p className="text-sm">MCNP SDEF:</p>
              <pre className="overflow-x-auto rounded-lg border border-border/50 bg-muted/30 p-3 font-mono text-xs whitespace-pre-wrap">
                {cards.mcnp}
              </pre>
              <p className="text-sm">Serpent src:</p>
              <pre className="overflow-x-auto rounded-lg border border-border/50 bg-muted/30 p-3 font-mono text-xs whitespace-pre-wrap">
                {cards.serpent}
              </pre>
            </div>
          )}
        </>
      )}
    </div>
  );
}
