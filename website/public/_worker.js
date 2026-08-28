/* eslint-disable no-undef -- Cloudflare Workers runtime provides Request/Response/Headers/URL globals */
/**
 * Cloudflare Pages advanced-mode Worker (Pages "Functions" advanced mode).
 *
 * Placing this file at dist/_worker.js (delivered here from public/_worker.js
 * by the Astro build) switches Pages into advanced mode, where this Worker
 * runs for *every* request before static assets are served.
 *
 * It performs HTTP content negotiation: if the client prefers
 * `text/markdown` (per RFC 9110 §12.5.1 — q-values, specificity, q=0
 * rejections), it serves the pre-built .md sibling that the
 * markdown-negotiation Astro integration emits at build time — same
 * canonical URL, markdown content-type. Otherwise it falls through to the
 * static asset server (`env.ASSETS.fetch`).
 *
 * Output is annotated with `Vary: Accept` so edge caches key correctly.
 */

const PRODUCES = ["text/html", "text/markdown"];

/**
 * @param {string} header
 * @returns {Array<{type: string; q: number; specificity: number}>}
 */
function parseAccept(header) {
  if (!header) return [];
  return header
    .split(",")
    .map((raw) => raw.trim())
    .filter(Boolean)
    .map((raw) => {
      const parts = raw.split(";").map((s) => s.trim());
      const type = parts[0].toLowerCase();
      let q = 1;
      for (const param of parts.slice(1)) {
        const eq = param.indexOf("=");
        const name = (eq === -1 ? param : param.slice(0, eq)).trim().toLowerCase();
        const value = eq === -1 ? "" : param.slice(eq + 1).trim();
        if (name === "q") {
          const parsed = Number(value);
          if (!Number.isNaN(parsed)) q = Math.max(0, Math.min(1, parsed));
        }
      }
      const specificity = type === "*/*" ? 0 : type.endsWith("/*") ? 1 : 2;
      return { type, q, specificity };
    });
}

/**
 * @param {{type: string; q: number; specificity: number}} entry
 * @param {string} candidate
 */
function matches(entry, candidate) {
  if (entry.type === "*/*") return true;
  if (entry.type.endsWith("/*")) return candidate.startsWith(entry.type.slice(0, -1));
  return entry.type === candidate;
}

/**
 * @param {string | null} header
 * @returns {string | null}
 */
function preferredType(header) {
  if (!header) return PRODUCES[0];
  const entries = parseAccept(header);
  if (entries.length === 0) return PRODUCES[0];

  let best = null;
  let bestQ = -1;
  let bestPosition = Infinity;

  for (let c = 0; c < PRODUCES.length; c++) {
    const candidate = PRODUCES[c];
    // Most specific matching range wins per type (RFC 9110 §12.5.1).
    let matched = null;
    let matchedPosition = Infinity;
    for (let i = 0; i < entries.length; i++) {
      const e = entries[i];
      if (!matches(e, candidate)) continue;
      if (
        matched === null ||
        e.specificity > matched.specificity ||
        (e.specificity === matched.specificity && i < matchedPosition)
      ) {
        matched = e;
        matchedPosition = i;
      }
    }
    if (matched === null) continue;
    if (matched.q <= 0) continue; // explicit rejection ("q=0")

    // Across candidates: highest q wins, tie-break on client order.
    if (matched.q > bestQ || (matched.q === bestQ && matchedPosition < bestPosition)) {
      bestQ = matched.q;
      bestPosition = matchedPosition;
      best = candidate;
    }
  }

  return best;
}

/**
 * Build candidate .md asset paths for a given request pathname.
 * Handles Astro's `directory` build format (e.g. /about/ -> dist/about/index.html).
 * @param {string} pathname
 * @returns {string[]}
 */
function mdCandidates(pathname) {
  const decoded = decodeURIComponent(pathname);
  if (decoded === "/" || decoded === "") return ["/index.md"];
  if (decoded.endsWith("/")) {
    const base = decoded.slice(0, -1);
    return [`${decoded}index.md`, `${base}.md`];
  }
  return [`${decoded}.md`, `${decoded}/index.md`];
}

/**
 * @param {Headers} headers
 */
function appendVaryAccept(headers) {
  const existing = headers.get("Vary");
  if (!existing) {
    headers.set("Vary", "Accept");
    return;
  }
  const tokens = existing.split(",").map((s) => s.trim().toLowerCase());
  if (!tokens.includes("accept")) headers.set("Vary", `${existing}, Accept`);
}

export default {
  /** @param {Request} request @param {{ ASSETS: { fetch: (req: Request) => Promise<Response> } }} env */
  async fetch(request, env) {
    const chosen = preferredType(request.headers.get("accept"));

    if (chosen === "text/markdown") {
      const url = new URL(request.url);
      for (const candidate of mdCandidates(url.pathname)) {
        try {
          const mdReq = new Request(new URL(candidate, url).toString());
          const mdRes = await env.ASSETS.fetch(mdReq);
          if (mdRes.ok) {
            const body = await mdRes.text();
            const headers = new Headers({
              "Content-Type": "text/markdown; charset=utf-8",
              "Vary": "Accept",
              "Cache-Control": "public, max-age=3600, s-maxage=86400",
            });
            return new Response(body, { status: 200, headers });
          }
        } catch {
          // asset not found — try next candidate
        }
      }
      // No .md available; fall through to the HTML page.
    }

    const res = await env.ASSETS.fetch(request);
    // Annotate the HTML response so caches differentiate by Accept.
    const outHeaders = new Headers(res.headers);
    appendVaryAccept(outHeaders);
    return new Response(res.body, { status: res.status, statusText: res.statusText, headers: outHeaders });
  },
};