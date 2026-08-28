import { useEffect, useRef, useState } from "react";
import type { WasmApi } from "../types/nucleide-wasm";

const BASE = import.meta.env.BASE_URL.endsWith("/")
  ? import.meta.env.BASE_URL
  : `${import.meta.env.BASE_URL}/`;
const WASM_URL = `${BASE}wasm/nucleide_wasm.js`;

export interface UseWasmResult {
  wasm: WasmApi | null;
  ready: boolean;
  error: string | null;
}

export function useWasm(): UseWasmResult {
  const wasmRef = useRef<WasmApi | null>(null);
  const [ready, setReady] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    let cancelled = false;

    async function load() {
      try {
        const mod = (await import(/* @vite-ignore */ WASM_URL)) as WasmApi;
        if (cancelled) return;
        await mod.default();
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
