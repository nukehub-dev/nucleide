// Copies runtime data files into public/data/ for the interactive tutorials.
// Fixtures stay under fixtures/ (single source of truth); this script minifies
// and stages them so the docs site can serve them statically.
import { mkdir, copyFile, readFile, readdir, writeFile } from "node:fs/promises";
import { fileURLToPath } from "node:url";
import path from "node:path";

const here = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(here, "..", "..");
const fixturesRoot = path.join(repoRoot, "fixtures");
const dst = path.join(here, "..", "public", "data");

await mkdir(dst, { recursive: true });

// Materials Compendium: minified JSON + LICENSE.
const raw = await readFile(path.join(fixturesRoot, "data", "MaterialsCompendium.json"), "utf8");
const minified = JSON.stringify(JSON.parse(raw));
await writeFile(path.join(dst, "MaterialsCompendium.json"), minified);
await copyFile(
  path.join(fixturesRoot, "data", "MaterialsCompendium.LICENSE"),
  path.join(dst, "MaterialsCompendium.LICENSE"),
);

// MCNP sample files consumed by the interactive tutorials.
await copyFile(
  path.join(fixturesRoot, "mcnp", "meshtal", "mcnp_meshtal_single_meshtal.txt"),
  path.join(dst, "meshtal_sample.txt"),
);
await copyFile(
  path.join(fixturesRoot, "mcnp", "xsdir", "dummy_xsdir"),
  path.join(dst, "xsdir_sample.txt"),
);

// Depletion chain sample.
await copyFile(
  path.join(fixturesRoot, "depletion", "chain_simple.xml"),
  path.join(dst, "chain_simple.xml"),
);

// MCNP full-deck samples consumed by the deck-editor tutorial.
await copyFile(
  path.join(fixturesRoot, "mcnp", "inp", "deck_minimal.txt"),
  path.join(dst, "deck_minimal.txt"),
);
await copyFile(
  path.join(fixturesRoot, "mcnp", "inp", "deck_l3.txt"),
  path.join(dst, "deck_l3.txt"),
);

// Serpent output samples consumed by the serpent-io tutorial.
await copyFile(
  path.join(fixturesRoot, "serpent", "serp2_res.m"),
  path.join(dst, "serpent_res_sample.m"),
);
await copyFile(
  path.join(fixturesRoot, "serpent", "sample2_dep.m"),
  path.join(dst, "serpent_dep_sample.m"),
);
await copyFile(
  path.join(fixturesRoot, "serpent", "sample_det.m"),
  path.join(dst, "serpent_det_sample.m"),
);

// FLUKA USRBIN sample consumed by the fluka-io tutorial.
await copyFile(
  path.join(fixturesRoot, "fluka", "fluka_usrbin_single.lis"),
  path.join(dst, "usrbin_sample.lis"),
);

// Theory-page SVG figures (docs/theory/figures/ is the source of truth).
const figuresSrc = path.join(repoRoot, "docs", "theory", "figures");
const figuresDst = path.join(here, "..", "public", "theory", "figures");
await mkdir(figuresDst, { recursive: true });
const figures = (await readdir(figuresSrc)).filter((f) => f.endsWith(".svg"));
for (const f of figures) {
  await copyFile(path.join(figuresSrc, f), path.join(figuresDst, f));
}

console.log(
  `sync-data: staged MaterialsCompendium.json (${(minified.length / 1e6).toFixed(1)} MB) + LICENSE, ` +
    `meshtal_sample.txt, xsdir_sample.txt, chain_simple.xml, deck_minimal.txt, deck_l3.txt, ` +
    `serpent_res_sample.m, serpent_dep_sample.m, serpent_det_sample.m, usrbin_sample.lis, ` +
    `${figures.length} theory figures`,
);
