// Copies runtime data files into public/data/ for the interactive tutorials.
// The Materials Compendium stays in fixtures/ (single source of truth); this
// script minifies and stages it so the docs site can serve it statically.
import { mkdir, copyFile, readFile, writeFile } from "node:fs/promises";
import { fileURLToPath } from "node:url";
import path from "node:path";

const here = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(here, "..", "..");
const src = path.join(repoRoot, "fixtures", "data");
const dst = path.join(here, "..", "public", "data");

await mkdir(dst, { recursive: true });

const raw = await readFile(path.join(src, "MaterialsCompendium.json"), "utf8");
const minified = JSON.stringify(JSON.parse(raw));
await writeFile(path.join(dst, "MaterialsCompendium.json"), minified);
await copyFile(
  path.join(src, "MaterialsCompendium.LICENSE"),
  path.join(dst, "MaterialsCompendium.LICENSE"),
);

console.log(
  `sync-data: staged MaterialsCompendium.json (${(minified.length / 1e6).toFixed(1)} MB) + LICENSE`,
);
