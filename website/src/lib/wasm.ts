import { useEffect, useRef, useState } from "react";
import type { WasmApi } from "../types/nucleide-wasm";

const BASE = import.meta.env.BASE_URL.endsWith("/")
  ? import.meta.env.BASE_URL
  : `${import.meta.env.BASE_URL}/`;
const WASM_URL = `${BASE}wasm/nucleide_wasm.js`;
const WASM_BG_URL = `${BASE}wasm/nucleide_wasm_bg.wasm`;

export interface UseWasmResult {
  wasm: WasmApi | null;
  ready: boolean;
  error: string | null;
}

/// Load the wasm-pack ES module without routing it through Vite's module
/// pipeline: `website/public/wasm/` is copied as-is (never transformed), so a
/// direct `import(WASM_URL)` fails under `astro dev` once the file exists.
/// Instead the JS is fetched as text, executed from a blob URL, and
/// initialized with an explicit `.wasm` URL (plain static fetch, which dev
/// and preview both serve from `public/`).
export async function loadWasmModule<T>(jsUrl: string, bgUrl: string): Promise<T> {
  const res = await fetch(jsUrl);
  if (!res.ok) {
    throw new Error(
      `wasm module not found at ${jsUrl} (HTTP ${res.status}); run \`npm run build:wasm\` from website/ first`,
    );
  }
  const src = await res.text();
  const blobUrl = URL.createObjectURL(new Blob([src], { type: "text/javascript" }));
  try {
    const mod = (await import(/* @vite-ignore */ blobUrl)) as T & {
      default: (path?: { module_or_path: string }) => Promise<void>;
    };
    await mod.default({ module_or_path: bgUrl });
    return mod;
  } finally {
    URL.revokeObjectURL(blobUrl);
  }
}

export function useWasm(): UseWasmResult {
  const wasmRef = useRef<WasmApi | null>(null);
  const [ready, setReady] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;

    async function load() {
      try {
        const mod = await loadWasmModule<WasmApi>(WASM_URL, WASM_BG_URL);
        if (cancelled) return;
        wasmRef.current = mod;
        setReady(true);
      } catch (e) {
        if (cancelled) return;
        setError(e instanceof Error ? e.message : String(e));
      }
    }

    load();
    return () => {
      cancelled = true;
    };
  }, []);

  return { wasm: wasmRef.current, ready, error };
}
