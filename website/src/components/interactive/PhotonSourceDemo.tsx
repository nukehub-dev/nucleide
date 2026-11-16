import { useEffect, useState } from "react";
import { useWasm } from "../../lib/wasm";
import type {
  PhotonSourceEvaluateResult,
  PhotonSourceInputJson,
  WasmApi,
} from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { DataTable } from "@nukehub/docs-kit/components/mdx/DataTable";

// Synthetic snapshots: three declared cooling times, two voxels, three
// groups. The browser never computes decay physics; it selects one declared
// snapshot and lets the shared Rust source contract validate and normalize it.
const EDGES_EV = [0, 1e5, 5e5, 2e6, 1e7];
const GROUP_LABELS = ["0–0.1", "0.1–0.5", "0.5–2", "2–10"];
const VOLUMES_CM3 = [1, 1, 1, 1];
const TIMES = [
  {
    coolingS: 0,
    strengths: [
      [420, 180, 40, 4],
      [160, 360, 90, 8],
      [80, 140, 380, 50],
      [20, 60, 220, 300],
    ],
  },
  {
    coolingS: 3600,
    strengths: [
      [300, 150, 35, 3.5],
      [120, 270, 80, 7],
      [60, 110, 300, 45],
      [15, 45, 180, 270],
    ],
  },
  {
    coolingS: 7200,
    strengths: [
      [210, 110, 28, 3],
      [90, 200, 65, 6],
      [45, 85, 230, 38],
      [12, 35, 140, 220],
    ],
  },
];

function initialText(): string[][] {
  return TIMES[0].strengths.map((row) => row.map(String));
}

interface CoolingSelection {
  requested: number[];
  indices: number[];
}

export function PhotonSourceDemo() {
  const { wasm, ready, error } = useWasm();
  const [timeIndex, setTimeIndex] = useState(0);
  const [strengthsText, setStrengthsText] = useState<string[][]>(initialText);
  const [result, setResult] = useState<PhotonSourceEvaluateResult | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);
  const [requestedText, setRequestedText] = useState("3600");
  const [selection, setSelection] = useState<CoolingSelection | null>(null);
  const [selectionError, setSelectionError] = useState<string | null>(null);

  function buildInput(text: string[][]): PhotonSourceInputJson {
    const strengths = text.map((row, voxel) =>
      row.map((cell, group) => {
        if (cell.trim() === "") {
          throw new Error(`voxel ${voxel + 1} group ${group + 1} is empty`);
        }
        const value = Number(cell);
        if (!Number.isFinite(value)) {
          throw new Error(`voxel ${voxel + 1} group ${group + 1} is not a finite number`);
        }
        return value;
      }),
    );
    return {
      meshShape: [2, 2, 1],
      lowerLeftCm: [0, 0, 0],
      upperRightCm: [2, 2, 1],
      voxelVolumesCm3: VOLUMES_CM3,
      edgesEv: EDGES_EV,
      coolingS: TIMES[timeIndex].coolingS,
      strengthsPhotonsPerS: strengths,
    };
  }

  function evaluate(api: WasmApi, text: string[][]) {
    try {
      const next = api.photonSourceEvaluate(buildInput(text));
      setResult(next);
      setLocalError(null);
    } catch (e) {
      setResult(null);
      setLocalError(e instanceof Error ? e.message : String(e));
    }
  }

  useEffect(() => {
    if (wasm) evaluate(wasm, initialText());
    // The first evaluation deliberately runs once per loaded module; later
    // edits are explicit so a malformed cell remains visible.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [wasm]);

  function selectTime(index: number) {
    setTimeIndex(index);
    setStrengthsText(TIMES[index].strengths.map((row) => row.map(String)));
    setSelection(null);
    setSelectionError(null);
  }

  function updateCell(voxel: number, group: number, value: string) {
    setStrengthsText((current) =>
      current.map((row, i) =>
        i === voxel ? row.map((cell, g) => (g === group ? value : cell)) : row,
      ),
    );
  }

  function checkSelection() {
    if (!wasm) return;
    try {
      const requested = requestedText
        .split(",")
        .map((part) => part.trim())
        .filter((part) => part.length > 0)
        .map((part) => Number(part));
      if (requested.length === 0) {
        throw new Error("enter at least one cooling offset");
      }
      if (requested.some((time) => !Number.isFinite(time) || time < 0)) {
        throw new Error("cooling offsets must be finite nonnegative seconds");
      }
      const indices = wasm.photonSourceSelectCooling(
        TIMES.map((entry) => entry.coolingS),
        requested,
      );
      setSelection({ requested, indices });
      setSelectionError(null);
    } catch (e) {
      setSelection(null);
      setSelectionError(e instanceof Error ? e.message : String(e));
    }
  }

  const total = result?.total_photons_per_s ?? null;
  const peakVoxelTotal = Math.max(0, ...(result?.voxel_totals_photons_per_s ?? []), 1e-30);

  return (
    <div className="not-prose my-6 space-y-6 rounded-xl border p-4 sm:p-6">
      <div className="flex flex-wrap items-center gap-2">
        <span className="text-sm font-medium">Cooling time</span>
        {TIMES.map((entry, index) => (
          <Button
            key={entry.coolingS}
            variant={index === timeIndex ? "default" : "outline"}
            size="sm"
            onClick={() => selectTime(index)}
          >
            {entry.coolingS} s
          </Button>
        ))}
        <span className="text-xs text-muted-foreground">
          declared snapshots, no decay calculation in the browser
        </span>
      </div>

      <div
        className="grid gap-6"
        style={{ gridTemplateColumns: "repeat(auto-fit, minmax(min(100%, 380px), 1fr))" }}
      >
        <div className="min-w-0 space-y-3">
          <div>
            <p className="text-sm font-semibold">Voxel spectra</p>
            <p className="text-xs text-muted-foreground">
              Integrated photons/s. Regular 2×2×1 mesh (four 1 cm³ cells); group edges{" "}
              {EDGES_EV.join(", ")} eV.
            </p>
          </div>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="text-left text-xs text-muted-foreground">
                  <th className="py-1 pr-3">Voxel</th>
                  {GROUP_LABELS.map((label) => (
                    <th key={label} className="whitespace-nowrap py-1 pr-1 text-center text-[10px]">
                      <span className="block">{label}</span>
                      <span className="block font-normal text-muted-foreground">MeV</span>
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {strengthsText.map((row, voxel) => (
                  <tr key={voxel}>
                    <td className="whitespace-nowrap py-1 pr-2 font-mono text-xs">V{voxel}</td>
                    {row.map((cell, group) => (
                      <td key={group} className="py-1 pr-2">
                        <Input
                          aria-label={`voxel ${voxel + 1} group ${group + 1} photons per second`}
                          className="h-8 w-14 px-1 text-center font-mono text-xs"
                          inputMode="decimal"
                          value={cell}
                          onChange={(event) => updateCell(voxel, group, event.target.value)}
                        />
                      </td>
                    ))}
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <div className="flex items-center gap-3">
            <Button onClick={() => wasm && evaluate(wasm, strengthsText)} disabled={!wasm}>
              Evaluate source
            </Button>
            {total !== null && (
              <span className="text-sm">
                Total <span className="font-mono">{total.toFixed(3)}</span> photons/s
              </span>
            )}
          </div>
          {error && <p className="text-sm text-destructive">{error}</p>}
          {localError && <p className="text-sm text-destructive">{localError}</p>}
          {!ready && !error && (
            <p className="text-sm text-muted-foreground">Loading the live demo…</p>
          )}
        </div>

        <div className="min-w-0 space-y-4">
          <div>
            <p className="text-sm font-semibold">Checks</p>
            {result ? (
              <p className="text-xs text-muted-foreground">
                {result.valid ? "valid source" : "invalid source"} · {result.issue_count} issue
                {result.issue_count === 1 ? "" : "s"}
                {result.truncated ? " (diagnostics capped)" : ""}
              </p>
            ) : (
              <p className="text-xs text-muted-foreground">Evaluate to run the shared validator.</p>
            )}
          </div>
          {result && result.issues.length > 0 && (
            <ul className="space-y-2">
              {result.issues.map((issue, index) => (
                <li key={`${issue.code}-${issue.location}-${index}`} className="text-xs">
                  <span
                    className={
                      issue.severity === "error"
                        ? "font-semibold text-destructive"
                        : issue.severity === "warning"
                          ? "font-semibold text-amber-600"
                          : "font-semibold text-muted-foreground"
                    }
                  >
                    {issue.severity}
                  </span>{" "}
                  <code>{issue.code}</code> at <code>{issue.location}</code>: {issue.message}{" "}
                  <span className="text-muted-foreground">{issue.remedy}</span>
                </li>
              ))}
            </ul>
          )}
          {result?.sampling_error && (
            <p className="text-xs text-destructive">{result.sampling_error}</p>
          )}
          {result && result.voxel_totals_photons_per_s.length > 0 && (
            <div>
              <p className="text-sm font-semibold">Per-voxel totals</p>
              <div className="mt-2 space-y-1.5">
                {result.voxel_totals_photons_per_s.map((total, voxel) => (
                  <div key={voxel} className="flex items-center gap-2 text-xs">
                    <span className="w-5 font-mono">V{voxel}</span>
                    <div className="h-2 flex-1 overflow-hidden rounded bg-muted">
                      <div
                        className="h-full rounded bg-primary"
                        style={{
                          width: `${Math.min(100, (total / peakVoxelTotal) * 100)}%`,
                        }}
                      />
                    </div>
                    <span className="w-24 text-right font-mono">{total.toFixed(1)} photons/s</span>
                  </div>
                ))}
              </div>
            </div>
          )}
          {result?.joint_probabilities && (
            <div>
              <p className="text-sm font-semibold">Joint voxel/group probabilities</p>
              <p className="text-xs text-muted-foreground">
                Sampling uses these pairs directly; independent spatial and energy marginals would
                invent cross-pairs that carry no source.
              </p>
              <DataTable
                className="mt-2"
                data={result.joint_probabilities.map((row, voxel) => {
                  const record: Record<string, string> = { voxel: `V${voxel}` };
                  row.forEach((probability, group) => {
                    record[`g${group}`] = `${(probability * 100).toFixed(2)}%`;
                  });
                  return record;
                })}
                columns={[
                  { key: "voxel", header: "Voxel", align: "left" as const },
                  ...GROUP_LABELS.map((label, group) => ({
                    key: `g${group}`,
                    header: (
                      <span>
                        <span className="block">{label}</span>
                        <span className="block text-[10px] font-normal text-muted-foreground">
                          MeV
                        </span>
                      </span>
                    ),
                    align: "right" as const,
                  })),
                ]}
              />
            </div>
          )}
        </div>
      </div>

      <div className="space-y-2 border-t pt-4">
        <p className="text-sm font-semibold">Exact cooling selection</p>
        <div className="flex flex-wrap items-center gap-1.5 text-xs text-muted-foreground">
          <span>Available offsets:</span>
          {TIMES.map((entry) => (
            <button
              key={entry.coolingS}
              type="button"
              onClick={() => {
                setRequestedText(String(entry.coolingS));
                setSelection(null);
                setSelectionError(null);
              }}
              className="rounded-md border px-2 py-0.5 font-mono text-[11px] transition-colors hover:bg-muted focus-visible:outline-none focus-visible:ring-2 focus-visible:ring-ring/50"
            >
              {entry.coolingS} s
            </button>
          ))}
        </div>
        <p className="text-xs text-muted-foreground">
          Requested offsets must match exactly; duplicates and near matches fail.
        </p>
        <div className="flex flex-wrap items-end gap-2">
          <div className="min-w-0 flex-1 space-y-1 sm:flex-none">
            <Label htmlFor="cooling-request" className="text-xs">
              Requested seconds (comma separated)
            </Label>
            <Input
              id="cooling-request"
              className="h-8 w-full font-mono text-xs sm:w-56"
              value={requestedText}
              onChange={(event) => {
                setRequestedText(event.target.value);
                setSelection(null);
                setSelectionError(null);
              }}
            />
          </div>
          <Button variant="outline" size="sm" onClick={checkSelection} disabled={!wasm}>
            Check selection
          </Button>
        </div>
        {selectionError && (
          <div
            role="alert"
            className="rounded-lg border border-destructive/40 bg-destructive/5 p-3 text-xs"
          >
            <p className="font-medium text-destructive">No exact match</p>
            <p className="mt-1 text-destructive">{selectionError}</p>
            <p className="mt-1 text-muted-foreground">
              Pick one of the available offsets above; near matches are not substituted.
            </p>
          </div>
        )}
        {selection && (
          <div
            role="status"
            aria-live="polite"
            className="rounded-lg border border-emerald-500/30 bg-emerald-500/5 p-3 text-xs"
          >
            <p className="font-medium text-emerald-700 dark:text-emerald-400">Exact match</p>
            <ul className="mt-2 space-y-1">
              {selection.requested.map((time, position) => {
                const index = selection.indices[position];
                const matchedTime = TIMES[index]?.coolingS;
                return (
                  <li key={`${time}-${position}`} className="flex flex-wrap items-center gap-1.5">
                    <span className="font-mono">{time} s</span>
                    <span className="text-muted-foreground">→</span>
                    <span className="rounded border bg-background px-1.5 py-0.5 font-mono text-[11px]">
                      snapshot {index}
                    </span>
                    <span className="text-muted-foreground">({matchedTime ?? "?"} s)</span>
                  </li>
                );
              })}
            </ul>
          </div>
        )}
      </div>
    </div>
  );
}
