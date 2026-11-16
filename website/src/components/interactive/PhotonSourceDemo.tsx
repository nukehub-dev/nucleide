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

// Synthetic snapshots: three declared cooling times, two voxels, three
// groups. The browser never computes decay physics; it selects one declared
// snapshot and lets the shared Rust source contract validate and normalize it.
const EDGES_EV = [0, 1, 3, 10];
const VOLUMES_CM3 = [1, 1];
const TIMES = [
  {
    coolingS: 0,
    strengths: [
      [100, 0, 20],
      [0, 300, 0],
    ],
  },
  {
    coolingS: 3600,
    strengths: [
      [71, 0, 14],
      [0, 212, 0],
    ],
  },
  {
    coolingS: 7200,
    strengths: [
      [50, 0, 10],
      [0, 150, 0],
    ],
  },
];

function initialText(): string[][] {
  return TIMES[0].strengths.map((row) => row.map(String));
}

export function PhotonSourceDemo() {
  const { wasm, ready, error } = useWasm();
  const [timeIndex, setTimeIndex] = useState(0);
  const [strengthsText, setStrengthsText] = useState<string[][]>(initialText);
  const [result, setResult] = useState<PhotonSourceEvaluateResult | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);
  const [requestedText, setRequestedText] = useState("3600");
  const [selection, setSelection] = useState<number[] | null>(null);

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
      meshShape: [2, 1, 1],
      lowerLeftCm: [0, 0, 0],
      upperRightCm: [2, 1, 1],
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
      setSelection(
        wasm.photonSourceSelectCooling(
          TIMES.map((entry) => entry.coolingS),
          requested,
        ),
      );
      setLocalError(null);
    } catch (e) {
      setSelection(null);
      setLocalError(e instanceof Error ? e.message : String(e));
    }
  }

  const total = result?.total_photons_per_s ?? null;

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

      <div className="grid gap-6 lg:grid-cols-2">
        <div className="min-w-0 space-y-3">
          <div>
            <p className="text-sm font-semibold">Voxel spectra</p>
            <p className="text-xs text-muted-foreground">
              Integrated photons/s. Regular mesh (2, 1, 1); group edges {EDGES_EV.join(", ")} eV;
              voxel volumes {VOLUMES_CM3.join(", ")} cm³.
            </p>
          </div>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="text-left text-xs text-muted-foreground">
                  <th className="py-1 pr-3">Voxel</th>
                  {EDGES_EV.slice(0, -1).map((low, group) => (
                    <th key={low} className="py-1 pr-3">
                      {low}–{EDGES_EV[group + 1]} eV
                    </th>
                  ))}
                </tr>
              </thead>
              <tbody>
                {strengthsText.map((row, voxel) => (
                  <tr key={voxel}>
                    <td className="py-1 pr-3 font-mono text-xs">
                      {voxel} · {VOLUMES_CM3[voxel]} cm³
                    </td>
                    {row.map((cell, group) => (
                      <td key={group} className="py-1 pr-3">
                        <Input
                          aria-label={`voxel ${voxel + 1} group ${group + 1} photons per second`}
                          className="h-8 w-24 font-mono text-xs"
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
          {result?.joint_probabilities && (
            <div>
              <p className="text-sm font-semibold">Joint voxel/group probabilities</p>
              <p className="text-xs text-muted-foreground">
                Sampling uses these pairs directly; independent spatial and energy marginals would
                invent cross-pairs that carry no source.
              </p>
              <div className="mt-2 overflow-x-auto">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="text-left text-xs text-muted-foreground">
                      <th className="py-1 pr-3">Voxel</th>
                      {EDGES_EV.slice(0, -1).map((low, group) => (
                        <th key={low} className="py-1 pr-3">
                          group {group}
                        </th>
                      ))}
                    </tr>
                  </thead>
                  <tbody>
                    {result.joint_probabilities.map((row, voxel) => (
                      <tr key={voxel}>
                        <td className="py-1 pr-3 font-mono text-xs">{voxel}</td>
                        {row.map((probability, group) => (
                          <td key={group} className="py-1 pr-3 font-mono text-xs">
                            {(probability * 100).toFixed(2)}%
                          </td>
                        ))}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
          )}
        </div>
      </div>

      <div className="space-y-2 border-t pt-4">
        <p className="text-sm font-semibold">Exact cooling selection</p>
        <p className="text-xs text-muted-foreground">
          Available offsets: {TIMES.map((entry) => entry.coolingS).join(", ")} s. Requested offsets
          must match exactly; duplicates and near matches fail.
        </p>
        <div className="flex flex-wrap items-end gap-2">
          <div className="space-y-1">
            <Label htmlFor="cooling-request" className="text-xs">
              Requested seconds (comma separated)
            </Label>
            <Input
              id="cooling-request"
              className="h-8 w-56 font-mono text-xs"
              value={requestedText}
              onChange={(event) => setRequestedText(event.target.value)}
            />
          </div>
          <Button variant="outline" size="sm" onClick={checkSelection} disabled={!wasm}>
            Check selection
          </Button>
          {selection && (
            <span className="font-mono text-xs">indices: [{selection.join(", ")}]</span>
          )}
        </div>
      </div>
    </div>
  );
}
