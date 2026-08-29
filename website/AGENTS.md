# Website NAD

## Purpose

This document owns the Nucleide documentation website under `website/`. It
builds an Astro + React + Tailwind site from `docs/` and exposes interactive
tutorials through a `wasm-bindgen` build in `bindings/wasm`.

## Ownership

This doc owns the website build, preview, sync, and end-to-end test workflow.
It covers `website/` and its generated artifacts (`website/public/wasm`,
`website/public/data`, `website/dist`, `website/src/content/docs`).
Parent-level Rust/Python/verification rules remain in the root `AGENTS.md`.

## Local Contracts

- **Docs source of truth**: `docs/`. The site syncs content into
  `website/src/content/docs` via `npm run sync-docs`.
- **WASM source of truth**: `bindings/wasm/`. Rebuild with `npm run build:wasm`
  after any Rust change that affects the WASM bindings.
- **Runtime data**: `fixtures/data/` is the source of truth for data files the
  interactive tutorials load at runtime (currently the Materials Compendium).
  `npm run sync-data` stages them, minified, into the git-ignored
  `website/public/data/`; it runs automatically via `predev`/`prebuild`.
- **Shared UI**: `@nukehub/docs-kit` provides layout, navigation, and
  markdown-negotiation integration. The dynamic favicon is customized through
  the kit's `SiteConfig.faviconPaths` field in `src/data/site.ts`.
- **Base path**: The site is configured with `base: "/nucleide"`. Preview and
  tests must account for this prefix.

## Work Guidance

After changing docs, components, or the WASM bindings:

1. `cd website`
2. `npm run build:wasm` — rebuild the WASM package if `bindings/wasm` changed.
3. `npm run check` — run Astro type checks.
4. `npm run build` — build the static site.
5. `npm run test:e2e:ci` — run the Playwright smoke tests.

## Verification

Run these in order from `website/`:

```bash
npm run build:wasm
npm run check
npm run build
npm run test:e2e:ci
```

Notes:

- `npm run test:e2e` runs the full Playwright suite locally with the HTML
  reporter.
- `npm run test:e2e:ci` runs only the Chromium project and is the command
  used in CI.
- Interactive tutorial components use `client:visible`; the E2E tests scroll
  each demo into view so Astro hydrates it and the WASM module loads.

## Child NAD Index

No children yet.
