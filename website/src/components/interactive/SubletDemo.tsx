import { useRef, useState } from "react";
import { useWasm } from "../../lib/wasm";
import type {
  SubletDose,
  SubletHazard,
  SubletIaea,
  SubletPointDose,
  SubletTransport,
} from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";

// Radiological-totals family over caller inventories (pure arithmetic, never
// vendored data): S4/S5 committed hazards, S6 transport ratio, S7 IAEA
// clearance index, and the S3 slab/point gamma dose. All coefficients, limits,
// levels, and attenuation tables are caller inputs typed in below.
type Tab = "hazards" | "transport" | "iaea" | "dose";

const TABS: { value: Tab; label: string }[] = [
  { value: "hazards", label: "S4/S5 hazards" },
  { value: "transport", label: "S6 transport" },
  { value: "iaea", label: "S7 IAEA clearance" },
  { value: "dose", label: "S3 gamma dose" },
];

interface InventoryRow {
  id: number;
  nuclide: string;
  activity: string;
  coeff: string;
}

const DEFAULT_HAZARD_ROWS: InventoryRow[] = [
  { id: 1, nuclide: "Co60", activity: "10", coeff: "3" },
  { id: 2, nuclide: "H3", activity: "5", coeff: "2" },
];

const DEFAULT_TRANSPORT_ROWS: InventoryRow[] = [
  { id: 1, nuclide: "Co60", activity: "1e12", coeff: "1" },
  { id: 2, nuclide: "H3", activity: "2e12", coeff: "1" },
];

const DEFAULT_IAEA_ROWS: InventoryRow[] = [
  { id: 1, nuclide: "Co60", activity: "10", coeff: "10" },
  { id: 2, nuclide: "H3", activity: "5", coeff: "5" },
];

interface DoseGroupRow {
  id: number;
  intensity: string;
  muAir: string;
  mu: string;
}

const DEFAULT_DOSE_GROUPS: DoseGroupRow[] = [{ id: 1, intensity: "0.5", muAir: "1", mu: "2" }];

function fmt(v: number): string {
  return Number.isFinite(v) ? v.toPrecision(5) : String(v);
}

export function SubletDemo() {
  const { wasm, ready, error } = useWasm();
  const [tab, setTab] = useState<Tab>("hazards");
  const [hazardRows, setHazardRows] = useState<InventoryRow[]>(DEFAULT_HAZARD_ROWS);
  const [transportRows, setTransportRows] = useState<InventoryRow[]>(DEFAULT_TRANSPORT_ROWS);
  const [iaeaRows, setIaeaRows] = useState<InventoryRow[]>(DEFAULT_IAEA_ROWS);
  const [massText, setMassText] = useState("2");
  const [ingestion, setIngestion] = useState<SubletHazard | null>(null);
  const [inhalation, setInhalation] = useState<SubletHazard | null>(null);
  const [transport, setTransport] = useState<SubletTransport | null>(null);
  const [iaea, setIaea] = useState<SubletIaea | null>(null);
  const [activityText, setActivityText] = useState("2");
  const [sourceMassText, setSourceMassText] = useState("1");
  const [distanceText, setDistanceText] = useState("1");
  const [doseGroups, setDoseGroups] = useState<DoseGroupRow[]>(DEFAULT_DOSE_GROUPS);
  const [slab, setSlab] = useState<SubletDose | null>(null);
  const [point, setPoint] = useState<SubletPointDose | null>(null);
  const [fractionsText, setFractionsText] = useState("0.25, 0.75");
  const [elementsText, setElementsText] = useState("4, 2\n8, 6");
  const [mixture, setMixture] = useState<number[] | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);
  const nextId = useRef(10);

  function clearError() {
    setLocalError(null);
  }

  function parseRow(row: InventoryRow, coeffLabel: string) {
    const nuclide = row.nuclide.trim();
    if (nuclide.length === 0) throw new Error("nuclide name is empty");
    const activity = parseFloat(row.activity);
    if (!Number.isFinite(activity) || activity < 0)
      throw new Error(`bad activity for ${nuclide} \`${row.activity}\` (expected >= 0)`);
    const coeff = parseFloat(row.coeff);
    if (!Number.isFinite(coeff) || coeff < 0)
      throw new Error(`bad ${coeffLabel} for ${nuclide} \`${row.coeff}\` (expected >= 0)`);
    return { nuclide, activityBq: activity, coeff };
  }

  function updateRows(
    rows: InventoryRow[],
    setRows: (r: InventoryRow[]) => void,
    id: number,
    patch: Partial<InventoryRow>,
  ) {
    setRows(rows.map((r) => (r.id === id ? { ...r, ...patch } : r)));
    clearError();
  }

  function addRow(rows: InventoryRow[], setRows: (r: InventoryRow[]) => void) {
    setRows([...rows, { id: nextId.current, nuclide: "", activity: "", coeff: "" }]);
    nextId.current += 1;
    clearError();
  }

  function removeRow(rows: InventoryRow[], setRows: (r: InventoryRow[]) => void, id: number) {
    setRows(rows.filter((r) => r.id !== id));
    clearError();
  }

  function rowInputs(
    rows: InventoryRow[],
    setRows: (r: InventoryRow[]) => void,
    coeffLabel: string,
    coeffUnit: string,
  ) {
    return (
      <div className="space-y-2">
        {rows.map((row) => (
          <div key={row.id} className="grid gap-3 sm:grid-cols-[1fr_1fr_1fr_auto]">
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Nuclide</Label>
              <Input
                value={row.nuclide}
                onChange={(e) => updateRows(rows, setRows, row.id, { nuclide: e.target.value })}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">Activity [Bq]</Label>
              <Input
                value={row.activity}
                onChange={(e) => updateRows(rows, setRows, row.id, { activity: e.target.value })}
                className="font-mono text-xs"
              />
            </div>
            <div className="space-y-1">
              <Label className="flex min-h-10 items-end">
                {coeffLabel} [{coeffUnit}]
              </Label>
              <Input
                value={row.coeff}
                onChange={(e) => updateRows(rows, setRows, row.id, { coeff: e.target.value })}
                className="font-mono text-xs"
              />
            </div>
            <div className="flex items-end">
              <Button variant="outline" size="sm" onClick={() => removeRow(rows, setRows, row.id)}>
                Remove
              </Button>
            </div>
          </div>
        ))}
        <div className="flex flex-wrap gap-2">
          <Button variant="outline" size="sm" onClick={() => addRow(rows, setRows)}>
            Add nuclide
          </Button>
        </div>
      </div>
    );
  }

  function runHazards() {
    if (!wasm) return;
    try {
      const ing = hazardRows.map((r) => {
        const p = parseRow(r, "ingestion coefficient");
        return { nuclide: p.nuclide, activityBq: p.activityBq, coeffSvPerBq: p.coeff };
      });
      const inh = hazardRows.map((r) => {
        const p = parseRow(r, "inhalation coefficient");
        return { nuclide: p.nuclide, activityBq: p.activityBq, coeffSvPerBq: p.coeff };
      });
      // One shared row set carries both coefficient columns in the demo: the
      // typed value feeds ingestion, and inhalation reuses it. Callers with
      // distinct e_ing/e_inh columns pass them separately in Python/Rust.
      setIngestion(wasm.subletIngestionHazard(ing));
      setInhalation(wasm.subletInhalationHazard(inh));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setIngestion(null);
      setInhalation(null);
    }
  }

  function runTransport() {
    if (!wasm) return;
    try {
      const entries = transportRows.map((r) => {
        const p = parseRow(r, "A2 limit");
        if (!(p.coeff > 0)) throw new Error(`bad A2 limit for ${p.nuclide} (expected > 0 TBq)`);
        return { nuclide: p.nuclide, activityBq: p.activityBq, a2Tbq: p.coeff };
      });
      setTransport(wasm.subletTransportRatio(entries));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setTransport(null);
    }
  }

  function runIaea() {
    if (!wasm) return;
    try {
      const mass = parseFloat(massText);
      if (!Number.isFinite(mass) || mass <= 0)
        throw new Error(`bad total mass \`${massText}\` (expected > 0 kg)`);
      const entries = iaeaRows.map((r) => {
        const p = parseRow(r, "IAEA level");
        if (!(p.coeff > 0)) throw new Error(`bad IAEA level for ${p.nuclide} (expected > 0 Bq/kg)`);
        return { nuclide: p.nuclide, activityBq: p.activityBq, limitBqPerKg: p.coeff };
      });
      setIaea(wasm.subletIaeaClearance(mass, entries));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setIaea(null);
    }
  }

  function runDose() {
    if (!wasm) return;
    try {
      const activity = parseFloat(activityText);
      if (!Number.isFinite(activity) || activity < 0)
        throw new Error(`bad specific activity \`${activityText}\` (expected >= 0 Bq/kg)`);
      const mass = parseFloat(sourceMassText);
      if (!Number.isFinite(mass) || mass < 0)
        throw new Error(`bad source mass \`${sourceMassText}\` (expected >= 0 kg)`);
      const distance = parseFloat(distanceText);
      if (!Number.isFinite(distance))
        throw new Error(`bad distance \`${distanceText}\` (expected finite, clamps at 0.3 m)`);
      const groups = doseGroups.map((g, i) => {
        const intensity = parseFloat(g.intensity);
        const muAir = parseFloat(g.muAir);
        const mu = parseFloat(g.mu);
        if (!Number.isFinite(intensity) || intensity < 0)
          throw new Error(`bad intensity in group ${i + 1} (expected >= 0)`);
        if (!Number.isFinite(muAir) || muAir < 0)
          throw new Error(`bad muAir in group ${i + 1} (expected >= 0)`);
        if (!Number.isFinite(mu) || mu < 0)
          throw new Error(`bad mu in group ${i + 1} (expected >= 0)`);
        return { intensity, muAir, mu };
      });
      setSlab(wasm.subletDoseSlab(activity, groups));
      setPoint(wasm.subletDosePoint(activity, mass, distance, groups));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setSlab(null);
      setPoint(null);
    }
  }

  function runMixture() {
    if (!wasm) return;
    try {
      const fractions = fractionsText.split(",").map((t, i) => {
        const v = parseFloat(t.trim());
        if (!Number.isFinite(v) || v < 0) throw new Error(`bad fraction ${i + 1} (expected >= 0)`);
        return v;
      });
      const rows = elementsText
        .split("\n")
        .map((l) => l.trim())
        .filter((l) => l.length > 0)
        .map((line, j) =>
          line.split(",").map((t) => {
            const v = parseFloat(t.trim());
            if (!Number.isFinite(v) || v < 0)
              throw new Error(`bad coefficient in element row ${j + 1} (expected >= 0)`);
            return v;
          }),
        );
      setMixture(wasm.subletDoseMixtureMu(fractions, rows));
      setLocalError(null);
    } catch (e) {
      setLocalError(e instanceof Error ? e.message : String(e));
      setMixture(null);
    }
  }

  // The He/dpa ratio-UQ surface lives in the damage demo; this component
  // covers the S3–S7 inventory-totals family only.

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
          <div className="flex flex-wrap gap-2" role="group" aria-label="Radiological family">
            {TABS.map((t) => (
              <Button
                key={t.value}
                onClick={() => {
                  setTab(t.value);
                  clearError();
                }}
                variant={tab === t.value ? "default" : "outline"}
                size="sm"
              >
                {t.label}
              </Button>
            ))}
          </div>

          {tab === "hazards" && (
            <div className="space-y-3">
              <p className="text-sm text-muted-foreground">
                Committed ingestion/inhalation dose Σ A_i·e_i [Sv] over caller 50-year committed
                coefficients (Sv/Bq, never vendored), with the excluding-tritium companion. The demo
                reuses one coefficient column for both kernels; the Python/Rust API takes e_ing and
                e_inh separately.
              </p>
              {rowInputs(hazardRows, setHazardRows, "Coefficient", "Sv/Bq")}
              <div className="flex flex-wrap gap-2">
                <Button onClick={runHazards}>Compute hazards</Button>
              </div>
              {ingestion && inhalation && (
                <table className="text-sm">
                  <tbody className="font-mono text-xs">
                    <tr>
                      <td className="pr-4">Ingestion total</td>
                      <td>
                        {fmt(ingestion.totalSv)} Sv (ex-tritium {fmt(ingestion.exTritiumSv)})
                      </td>
                    </tr>
                    <tr>
                      <td className="pr-4">Inhalation total</td>
                      <td>
                        {fmt(inhalation.totalSv)} Sv (ex-tritium {fmt(inhalation.exTritiumSv)})
                      </td>
                    </tr>
                  </tbody>
                </table>
              )}
            </div>
          )}

          {tab === "transport" && (
            <div className="space-y-3">
              <p className="text-sm text-muted-foreground">
                Transport Bq/A₂ ratio Σ A_i/(A2,i·C₂) (dimensionless, C₂ = 10¹² TBq→Bq) with the
                effective A₂ it defines. Limits are caller inputs, never vendored.
              </p>
              {rowInputs(transportRows, setTransportRows, "A2 limit", "TBq")}
              <div className="flex flex-wrap gap-2">
                <Button onClick={runTransport}>Compute transport ratio</Button>
              </div>
              {transport && (
                <table className="text-sm">
                  <tbody className="font-mono text-xs">
                    <tr>
                      <td className="pr-4">Total Bq/A₂</td>
                      <td>{fmt(transport.ratio)}</td>
                    </tr>
                    <tr>
                      <td className="pr-4">Total activity</td>
                      <td>{fmt(transport.totalBq)} Bq</td>
                    </tr>
                    <tr>
                      <td className="pr-4">Effective A₂</td>
                      <td>{fmt(transport.effectiveA2Tbq)} TBq</td>
                    </tr>
                  </tbody>
                </table>
              )}
            </div>
          )}

          {tab === "iaea" && (
            <div className="space-y-3">
              <p className="text-sm text-muted-foreground">
                IAEA clearance index Σ A_i/(M·L_i) (dimensionless) with the ≤ 1 screening class
                (boundary included) and dominant contributor. Levels are caller inputs, never
                vendored. Screening arithmetic only, never a compliance decision.
              </p>
              <div className="grid gap-3 sm:grid-cols-2">
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Total mass [kg]</Label>
                  <Input
                    value={massText}
                    onChange={(e) => {
                      setMassText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
              </div>
              {rowInputs(iaeaRows, setIaeaRows, "IAEA level", "Bq/kg")}
              <div className="flex flex-wrap gap-2">
                <Button onClick={runIaea}>Compute IAEA clearance</Button>
              </div>
              {iaea && (
                <table className="text-sm">
                  <tbody className="font-mono text-xs">
                    <tr>
                      <td className="pr-4">Clearance index</td>
                      <td>{fmt(iaea.index)}</td>
                    </tr>
                    <tr>
                      <td className="pr-4">Screening class</td>
                      <td>{iaea.clearanceClass}</td>
                    </tr>
                    <tr>
                      <td className="pr-4">Dominant contributor</td>
                      <td>
                        {iaea.maxNuclide === null
                          ? "—"
                          : `${iaea.maxNuclide} (fraction ${fmt(iaea.maxFraction)})`}
                      </td>
                    </tr>
                  </tbody>
                </table>
              )}
            </div>
          )}

          {tab === "dose" && (
            <div className="space-y-3">
              <p className="text-sm text-muted-foreground">
                Slab dose C·B/2·Σ μa/μm·Sγ and point dose C·Σ μa/(4πr²)·e^(−μr)·m_s·Sγ [Sv/h] over
                caller gamma groups (yields plus air/mixture attenuation, never vendored). Point
                distances below 0.3 m clamp loudly. The mixture fold builds μm from elemental values
                below.
              </p>
              <div className="grid gap-3 sm:grid-cols-3">
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Specific activity [Bq/kg]</Label>
                  <Input
                    value={activityText}
                    onChange={(e) => {
                      setActivityText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Source mass [kg]</Label>
                  <Input
                    value={sourceMassText}
                    onChange={(e) => {
                      setSourceMassText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Distance [m]</Label>
                  <Input
                    value={distanceText}
                    onChange={(e) => {
                      setDistanceText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
              </div>
              <div className="space-y-2">
                {doseGroups.map((g) => (
                  <div key={g.id} className="grid gap-3 sm:grid-cols-[1fr_1fr_1fr_auto]">
                    <div className="space-y-1">
                      <Label className="flex min-h-10 items-end">Yield Iᵢ</Label>
                      <Input
                        value={g.intensity}
                        onChange={(e) =>
                          setDoseGroups(
                            doseGroups.map((x) =>
                              x.id === g.id ? { ...x, intensity: e.target.value } : x,
                            ),
                          )
                        }
                        className="font-mono text-xs"
                      />
                    </div>
                    <div className="space-y-1">
                      <Label className="flex min-h-10 items-end">μa</Label>
                      <Input
                        value={g.muAir}
                        onChange={(e) =>
                          setDoseGroups(
                            doseGroups.map((x) =>
                              x.id === g.id ? { ...x, muAir: e.target.value } : x,
                            ),
                          )
                        }
                        className="font-mono text-xs"
                      />
                    </div>
                    <div className="space-y-1">
                      <Label className="flex min-h-10 items-end">μm</Label>
                      <Input
                        value={g.mu}
                        onChange={(e) =>
                          setDoseGroups(
                            doseGroups.map((x) =>
                              x.id === g.id ? { ...x, mu: e.target.value } : x,
                            ),
                          )
                        }
                        className="font-mono text-xs"
                      />
                    </div>
                    <div className="flex items-end">
                      <Button
                        variant="outline"
                        size="sm"
                        onClick={() => {
                          setDoseGroups(doseGroups.filter((x) => x.id !== g.id));
                          clearError();
                        }}
                      >
                        Remove
                      </Button>
                    </div>
                  </div>
                ))}
                <div className="flex flex-wrap gap-2">
                  <Button
                    variant="outline"
                    size="sm"
                    onClick={() => {
                      setDoseGroups([
                        ...doseGroups,
                        { id: nextId.current, intensity: "", muAir: "", mu: "" },
                      ]);
                      nextId.current += 1;
                      clearError();
                    }}
                  >
                    Add group
                  </Button>
                  <Button onClick={runDose}>Compute gamma dose</Button>
                </div>
              </div>
              {slab && point && (
                <table className="text-sm">
                  <tbody className="font-mono text-xs">
                    <tr>
                      <td className="pr-4">Slab dose</td>
                      <td>{fmt(slab.doseSvPerH)} Sv/h</td>
                    </tr>
                    <tr>
                      <td className="pr-4">Point dose</td>
                      <td>
                        {fmt(point.doseSvPerH)} Sv/h at {fmt(point.distanceUsedM)} m
                        {point.clamped ? " (clamped)" : ""}
                      </td>
                    </tr>
                  </tbody>
                </table>
              )}
              <p className="text-sm text-muted-foreground">
                Mixture fold μ_m = Σ f_j·μ_j: fractions (comma-separated, summing to 1) weight every
                group of each element row (one row per line, comma-separated groups).
              </p>
              <div className="grid gap-3 sm:grid-cols-2">
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Mass fractions</Label>
                  <Input
                    value={fractionsText}
                    onChange={(e) => {
                      setFractionsText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
                <div className="space-y-1">
                  <Label className="flex min-h-10 items-end">Element rows</Label>
                  <Input
                    value={elementsText}
                    onChange={(e) => {
                      setElementsText(e.target.value);
                      clearError();
                    }}
                    className="font-mono text-xs"
                  />
                </div>
              </div>
              <div className="flex flex-wrap gap-2">
                <Button onClick={runMixture}>Fold mixture</Button>
              </div>
              {mixture && (
                <p className="font-mono text-xs">
                  Mixture μm: {mixture.map((v) => fmt(v)).join(", ")}
                </p>
              )}
            </div>
          )}
        </>
      )}
    </div>
  );
}
