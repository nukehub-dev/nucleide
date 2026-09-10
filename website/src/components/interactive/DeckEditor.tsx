import { useState } from "react";
import { useWasm } from "../../lib/wasm";
import type {
  DeckCellJson,
  FillJson,
  ImportanceJson,
  LatticeJson,
  TallyJson,
  TransformJson,
  UniverseJson,
  VolumeJson,
  WasmDeckProblem,
  ModeJson,
} from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Checkbox } from "@nukehub/docs-kit/components/ui/Checkbox";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { Textarea } from "@nukehub/docs-kit/components/ui/Textarea";
import { DataTable } from "@nukehub/docs-kit/components/mdx/DataTable";

const BASE = import.meta.env.BASE_URL.endsWith("/")
  ? import.meta.env.BASE_URL
  : `${import.meta.env.BASE_URL}/`;
const DECK_L3_SAMPLE_URL = `${BASE}data/deck_l3.txt`;

// Byte-identical to fixtures/mcnp/inp/deck_minimal.txt (synthetic sample deck).
const DECK_MINIMAL = `Minimal pin-cell message
Minimal pin-cell deck
c -- cell cards --
1 1 -10.0 -1 imp:n=1 $ fuel pin
2 2 -1.0 1 -2 imp:n=1 $ coolant
3 0 2 imp:n=0 $ outside world

c -- surface cards --
1 so 10.0
2 so 11.0
3 px 5.0

mode n
m1 92235 0.04 92238 0.96 $ enriched uranium
m2 1001 2.0 8016 1.0 $ water
`;

interface DeckSnapshot {
  message: string;
  title: string;
  cells: DeckCellJson[];
  materialNumbers: number[];
  dataNames: string[];
  mode: ModeJson;
  transforms: TransformJson[];
  universes: UniverseJson[];
  lattices: LatticeJson[];
  fills: FillJson[];
  importances: ImportanceJson[];
  volumes: VolumeJson[];
  tallies: TallyJson[];
  inventory: Record<number, number>;
  notes: string[];
  verdict: string;
  verdictOk: boolean;
  dumps: string;
}

function snapshot(deck: WasmDeckProblem): DeckSnapshot {
  let verdict = "✓ valid";
  let verdictOk = true;
  try {
    deck.validate();
  } catch (e) {
    verdict = e instanceof Error ? e.message : String(e);
    verdictOk = false;
  }
  return {
    message: deck.message,
    title: deck.title,
    cells: deck.cells(),
    materialNumbers: deck.materialNumbers(),
    dataNames: deck.dataNames(),
    mode: deck.mode(),
    transforms: deck.transforms(),
    universes: deck.universes(),
    lattices: deck.lattices(),
    fills: deck.fills(),
    importances: deck.importances(),
    volumes: deck.volumes(),
    tallies: deck.tallies(),
    inventory: deck.cellInventory(),
    notes: deck.validationNotes(),
    verdict,
    verdictOk,
    dumps: deck.dumps(),
  };
}

export function DeckEditor() {
  const { wasm, ready, error } = useWasm();
  const [text, setText] = useState(DECK_MINIMAL);
  const [deck, setDeck] = useState<WasmDeckProblem | null>(null);
  const [snap, setSnap] = useState<DeckSnapshot | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);
  const [loadingSample, setLoadingSample] = useState(false);

  const [cellNum, setCellNum] = useState("1");
  const [density, setDensity] = useState("-7.0");
  const [material, setMaterial] = useState("1");
  const [universe, setUniverse] = useState("0");
  const [notTruncated, setNotTruncated] = useState(false);
  const [lattice, setLattice] = useState("");
  const [fill, setFill] = useState("0");
  const [modeInput, setModeInput] = useState("n p");

  function clearError() {
    setLocalError(null);
  }

  async function loadL3Sample() {
    setLoadingSample(true);
    try {
      const res = await fetch(DECK_L3_SAMPLE_URL);
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      setText(await res.text());
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
      const parsed = wasm.WasmDeckProblem.fromText(text);
      setDeck(parsed);
      setSnap(snapshot(parsed));
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setDeck(null);
      setSnap(null);
    }
  }

  function applyEdit(edit: (d: WasmDeckProblem) => void) {
    if (!wasm || !deck) return;
    try {
      edit(deck);
      setSnap(snapshot(deck));
      clearError();
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
    }
  }

  function revalidate() {
    if (!deck) return;
    setSnap(snapshot(deck));
  }

  function parseCell(): number {
    const n = parseInt(cellNum, 10);
    if (!Number.isInteger(n) || n < 1) throw new Error(`bad cell number \`${cellNum}\``);
    return n;
  }

  function parseDensity(): number {
    const v = parseFloat(density);
    if (!Number.isFinite(v)) throw new Error(`bad density \`${density}\``);
    return v;
  }

  function parseMaterial(): number {
    const n = parseInt(material, 10);
    if (!Number.isInteger(n)) throw new Error(`bad material number \`${material}\``);
    return n;
  }

  function parseUniverse(): number {
    const n = parseInt(universe, 10);
    if (!Number.isInteger(n)) throw new Error(`bad universe \`${universe}\``);
    return n;
  }

  function parseLattice(): number | undefined {
    if (lattice.trim() === "") return undefined;
    const n = parseInt(lattice, 10);
    if (!Number.isInteger(n)) throw new Error(`bad lattice \`${lattice}\``);
    return n;
  }

  function parseFill(): number {
    const n = parseInt(fill, 10);
    if (!Number.isInteger(n)) throw new Error(`bad fill universe \`${fill}\``);
    return n;
  }

  function parseMode(): string[] {
    const parts = modeInput.split(/\s+/).filter((p) => p.length > 0);
    if (parts.length === 0) throw new Error(`bad mode particles \`${modeInput}\``);
    return parts;
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
          <div className="space-y-2">
            <Label>MCNP deck text</Label>
            <Textarea
              value={text}
              onChange={(e) => {
                setText(e.target.value);
                clearError();
              }}
              className="font-mono text-xs"
            />
          </div>

          <div className="flex flex-wrap gap-2">
            <Button onClick={run}>Parse</Button>
            <Button variant="outline" onClick={loadL3Sample} disabled={loadingSample}>
              {loadingSample ? "Loading…" : "Load L3 sample"}
            </Button>
            {deck && (
              <Button variant="outline" onClick={revalidate}>
                Validate
              </Button>
            )}
          </div>

          {snap && (
            <>
              <div className="grid gap-2 text-sm sm:grid-cols-2">
                <p>
                  <span className="font-medium">Title:</span> {snap.title}
                </p>
                <p>
                  <span className="font-medium">Validation:</span>{" "}
                  <span className={snap.verdictOk ? "text-green-700" : "text-red-700"}>
                    {snap.verdictOk ? "✓ valid" : snap.verdict}
                  </span>
                </p>
              </div>

              <DataTable
                data={snap.cells.map((c) => ({
                  cell: `Cell ${c.num}`,
                  mat: c.mat,
                  dens: c.dens ?? "—",
                  geom: <span className="font-mono">{c.geom}</span>,
                  params: <span className="font-mono">{c.params.join(" ")}</span>,
                }))}
                columns={[
                  { key: "cell", header: "Cell" },
                  { key: "mat", header: "Mat", align: "right" },
                  { key: "dens", header: "Dens", align: "right" },
                  { key: "geom", header: "Geometry" },
                  { key: "params", header: "Params" },
                ]}
              />

              <div className="grid gap-2 text-sm sm:grid-cols-2">
                <p>
                  <span className="font-medium">Materials:</span>{" "}
                  {snap.materialNumbers.length > 0 ? snap.materialNumbers.join(", ") : "none"}
                </p>
                <p>
                  <span className="font-medium">Data cards:</span>{" "}
                  {snap.dataNames.filter((n) => n !== "").join(", ") || "none"}
                </p>
                <p>
                  <span className="font-medium">Mode:</span> {snap.mode.particles.join(" ")}
                </p>
                <p>
                  <span className="font-medium">L3 counts:</span> {snap.transforms.length}{" "}
                  transforms, {snap.universes.length} universes, {snap.lattices.length} lattices,{" "}
                  {snap.fills.length} fills, {snap.importances.length} importances,{" "}
                  {snap.volumes.length} volumes, {snap.tallies.length} tallies
                </p>
              </div>

              <DataTable
                data={Object.entries(snap.inventory).map(([cell, mat]) => ({
                  cell: `Cell ${cell}`,
                  material: mat,
                }))}
                columns={[
                  { key: "cell", header: "Cell" },
                  { key: "material", header: "Material", align: "right" },
                ]}
              />

              {snap.notes.length > 0 && (
                <div className="space-y-1 text-sm">
                  <p className="font-medium">Validation notes (not errors):</p>
                  <ul className="list-disc pl-5 text-muted-foreground">
                    {snap.notes.map((n) => (
                      <li key={n}>{n}</li>
                    ))}
                  </ul>
                </div>
              )}

              <div className="space-y-3 rounded-lg border border-border/50 p-3">
                <div className="grid items-end gap-2 sm:grid-cols-4">
                  <div className="space-y-1">
                    <Label>Cell #</Label>
                    <Input
                      value={cellNum}
                      onChange={(e) => {
                        setCellNum(e.target.value);
                        clearError();
                      }}
                    />
                  </div>
                  <div className="space-y-1">
                    <Label>Density</Label>
                    <Input
                      value={density}
                      onChange={(e) => {
                        setDensity(e.target.value);
                        clearError();
                      }}
                    />
                  </div>
                  <div className="space-y-1">
                    <Label>Material</Label>
                    <Input
                      value={material}
                      onChange={(e) => {
                        setMaterial(e.target.value);
                        clearError();
                      }}
                    />
                  </div>
                </div>
                <div className="flex flex-wrap gap-2">
                  <Button
                    size="sm"
                    variant="outline"
                    onClick={() => applyEdit((d) => d.setCellDensity(parseCell(), parseDensity()))}
                  >
                    Set density
                  </Button>
                  <Button
                    size="sm"
                    variant="outline"
                    onClick={() =>
                      applyEdit((d) => d.setCellMaterial(parseCell(), parseMaterial()))
                    }
                  >
                    Set material
                  </Button>
                </div>

                <div className="grid items-end gap-2 sm:grid-cols-4">
                  <div className="space-y-1">
                    <Label>Universe</Label>
                    <Input
                      value={universe}
                      onChange={(e) => {
                        setUniverse(e.target.value);
                        clearError();
                      }}
                    />
                  </div>
                  <div className="space-y-1">
                    <Label>Lattice (empty clears)</Label>
                    <Input
                      value={lattice}
                      onChange={(e) => {
                        setLattice(e.target.value);
                        clearError();
                      }}
                    />
                  </div>
                  <div className="space-y-1">
                    <Label>Fill universe</Label>
                    <Input
                      value={fill}
                      onChange={(e) => {
                        setFill(e.target.value);
                        clearError();
                      }}
                    />
                  </div>
                  <Checkbox
                    id="deck-not-truncated"
                    checked={notTruncated}
                    onCheckedChange={setNotTruncated}
                  >
                    Not truncated (U=-n)
                  </Checkbox>
                </div>
                <div className="flex flex-wrap gap-2">
                  <Button
                    size="sm"
                    variant="outline"
                    onClick={() =>
                      applyEdit((d) =>
                        d.setCellUniverse(parseCell(), parseUniverse(), notTruncated),
                      )
                    }
                  >
                    Set universe
                  </Button>
                  <Button
                    size="sm"
                    variant="outline"
                    onClick={() => applyEdit((d) => d.setCellLattice(parseCell(), parseLattice()))}
                  >
                    Set lattice
                  </Button>
                  <Button
                    size="sm"
                    variant="outline"
                    onClick={() => applyEdit((d) => d.setCellFill(parseCell(), parseFill()))}
                  >
                    Set fill
                  </Button>
                </div>

                <div className="grid items-end gap-2 sm:grid-cols-2">
                  <div className="space-y-1">
                    <Label>Mode particles (space-separated)</Label>
                    <Input
                      value={modeInput}
                      onChange={(e) => {
                        setModeInput(e.target.value);
                        clearError();
                      }}
                    />
                  </div>
                  <div>
                    <Button
                      size="sm"
                      variant="outline"
                      onClick={() => applyEdit((d) => d.setMode(parseMode()))}
                    >
                      Set mode
                    </Button>
                  </div>
                </div>
              </div>

              <div className="space-y-1">
                <Label>Deck text (dumps)</Label>
                <pre className="overflow-x-auto rounded-lg border border-border/50 bg-muted/30 p-3 font-mono text-xs whitespace-pre-wrap">
                  {snap.dumps}
                </pre>
              </div>
            </>
          )}
        </>
      )}
    </div>
  );
}
