import { useRef, useState } from "react";
import { useWasm } from "../../lib/wasm";
import type { McplSummary } from "../../types/nucleide-wasm";
import { Button } from "@nukehub/docs-kit/components/ui/Button";
import { Input } from "@nukehub/docs-kit/components/ui/Input";
import { Label } from "@nukehub/docs-kit/components/ui/Label";
import { DataTable } from "@nukehub/docs-kit/components/mdx/DataTable";

const BASE = import.meta.env.BASE_URL.endsWith("/")
  ? import.meta.env.BASE_URL
  : `${import.meta.env.BASE_URL}/`;

async function fetchBytes(url: string): Promise<Uint8Array> {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`HTTP ${res.status}`);
  return new Uint8Array(await res.arrayBuffer());
}

export function McplDemo() {
  const { wasm, ready, error } = useWasm();
  const [summary, setSummary] = useState<McplSummary | null>(null);
  const [mcplBytes, setMcplBytes] = useState<Uint8Array | null>(null);
  const [sswBytes, setSswBytes] = useState<Uint8Array | null>(null);
  const [writtenLen, setWrittenLen] = useState<number | null>(null);
  const [convertedLen, setConvertedLen] = useState<number | null>(null);
  const [roundTrip, setRoundTrip] = useState<string | null>(null);
  const [localError, setLocalError] = useState<string | null>(null);
  const [fileName, setFileName] = useState<string | null>(null);
  // The native file control renders its own hover tooltip in some browsers,
  // so it stays hidden: the Browse button below opens it programmatically.
  const fileRef = useRef<HTMLInputElement>(null);

  function clearError() {
    setLocalError(null);
  }

  async function onFile(e: React.ChangeEvent<HTMLInputElement>) {
    if (!wasm) return;
    const f = e.target.files?.[0];
    if (!f) return;
    try {
      const bytes = new Uint8Array(await f.arrayBuffer());
      const out = wasm.readMcpl(bytes);
      setSummary(out);
      setMcplBytes(bytes);
      setFileName(f.name);
      clearError();
    } catch (err) {
      setLocalError(err instanceof Error ? err.message : String(err));
      setSummary(null);
    }
  }

  async function loadGolden() {
    if (!wasm) return;
    try {
      const bytes = await fetchBytes(`${BASE}data/ssw2mcpl_expected.mcpl`);
      setSummary(wasm.readMcpl(bytes));
      setMcplBytes(bytes);
      clearError();
    } catch (err) {
      setLocalError(err instanceof Error ? err.message : String(err));
    }
  }

  async function loadReferenceSsw() {
    try {
      setSswBytes(await fetchBytes(`${BASE}data/reference.w`));
      clearError();
    } catch (err) {
      setLocalError(err instanceof Error ? err.message : String(err));
    }
  }

  function writeSample() {
    if (!wasm) return;
    try {
      const out: Uint8Array = wasm.writeMcpl(
        { srcname: "nucleide-demo", comments: ["synthetic 2-particle probe"], hasUserflags: true },
        [
          {
            ekin: 2.5,
            position: [1, -2, 0.5],
            direction: [0, 0, 1],
            time: 3.0,
            weight: 1.0,
            pdgcode: 2112,
            userflags: 100,
          },
          {
            ekin: 0.662,
            position: [0, 0, 0],
            direction: [1, 0, 0],
            time: 0.0,
            weight: 0.5,
            pdgcode: 22,
            userflags: 200,
          },
        ],
      );
      setWrittenLen(out.length);
      setMcplBytes(out);
      setSummary(wasm.readMcpl(out));
      clearError();
    } catch (err) {
      setLocalError(err instanceof Error ? err.message : String(err));
    }
  }

  async function convertSsw2Mcpl() {
    if (!wasm) return;
    try {
      const ssw = sswBytes ?? (await fetchBytes(`${BASE}data/reference.w`));
      setSswBytes(ssw);
      const out: Uint8Array = wasm.ssw2mcpl(ssw, [100, 200], ["neutron", "gamma"], {});
      setConvertedLen(out.length);
      setMcplBytes(out);
      setSummary(wasm.readMcpl(out));
      clearError();
    } catch (err) {
      setLocalError(err instanceof Error ? err.message : String(err));
    }
  }

  async function convertMcpl2Ssw() {
    if (!wasm) return;
    try {
      const mcpl = mcplBytes ?? (await fetchBytes(`${BASE}data/ssw2mcpl_expected.mcpl`));
      const ref = sswBytes ?? (await fetchBytes(`${BASE}data/reference.w`));
      setSswBytes(ref);
      const out: Uint8Array = wasm.mcpl2ssw(mcpl, ref);
      setRoundTrip(`SSW bytes: ${out.length}, tracks: 2`);
      clearError();
    } catch (err) {
      setLocalError(err instanceof Error ? err.message : String(err));
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
          <div className="space-y-2">
            <Label>MCPL file upload (bytes, gzip sniffed by magic)</Label>
            <Input ref={fileRef} type="file" onChange={onFile} className="hidden" />
            <div className="flex flex-wrap gap-2">
              <Button onClick={() => fileRef.current?.click()}>Browse…</Button>
              <Button onClick={loadGolden}>Load golden MCPL</Button>
              <Button onClick={writeSample}>Write MCPL</Button>
            </div>
            {fileName !== null && (
              <p className="text-sm">
                File: <span className="font-mono">{fileName}</span>
              </p>
            )}
            {writtenLen !== null && (
              <p className="text-sm">
                Wrote MCPL bytes: <span className="font-mono">{writtenLen}</span>
              </p>
            )}
          </div>
          <div className="space-y-2">
            <Label>SSW conversion (neutron/gamma v1, 2-track pairing)</Label>
            <div className="flex flex-wrap gap-2">
              <Button onClick={loadReferenceSsw}>Load reference SSW</Button>
              <Button onClick={convertSsw2Mcpl}>Convert SSW→MCPL</Button>
              <Button onClick={convertMcpl2Ssw}>Convert MCPL→SSW</Button>
            </div>
            {convertedLen !== null && (
              <p className="text-sm">
                Converted MCPL bytes: <span className="font-mono">{convertedLen}</span> (particles:
                2)
              </p>
            )}
            {roundTrip && <p className="text-sm">{roundTrip}</p>}
          </div>
          {summary && (
            <div className="space-y-2">
              <p className="text-sm">
                MCPL v{summary.version}, particles: {summary.nparticles}, src:{" "}
                <span className="font-mono">{summary.srcname}</span>
                {summary.truncated && " (particle list capped at 200)"}
              </p>
              <DataTable
                data={summary.particles.map((p, i) => ({
                  i,
                  ekin: p.ekin.toFixed(3),
                  pdg: p.pdgcode,
                  userflags: p.userflags,
                }))}
                columns={[
                  { key: "i", header: "#" },
                  { key: "ekin", header: "ekin [MeV]", align: "right" },
                  { key: "pdg", header: "PDG", align: "right" },
                  { key: "userflags", header: "userflags", align: "right" },
                ]}
              />
            </div>
          )}
        </>
      )}
    </div>
  );
}
