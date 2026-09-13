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
await copyFile(
  path.join(fixturesRoot, "mcnp", "mctal", "synthetic_tally_bodies.mctal"),
  path.join(dst, "mctal_sample.mctal"),
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
// sample_res.m (multi-burnup-block) so the k-eff convergence chart draws a
// real line; serp2_res.m has a single block and would chart one point.
await copyFile(
  path.join(fixturesRoot, "serpent", "sample_res.m"),
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

// Activation output samples consumed by the activation tutorial.
await copyFile(
  path.join(fixturesRoot, "alara", "output", "sample2.out"),
  path.join(dst, "alara_output_sample.out"),
);
await copyFile(
  path.join(fixturesRoot, "fispact", "inventory.fis"),
  path.join(dst, "fispact_inventory_sample.fis"),
);

// FLUKA USRBIN sample consumed by the fluka-io tutorial.
await copyFile(
  path.join(fixturesRoot, "fluka", "fluka_usrbin_single.lis"),
  path.join(dst, "usrbin_sample.lis"),
);

// MCPL/SSW samples consumed by the mcpl-io tutorial (binary staging:
// golden MCPL output + synthetic SSW reference; gzip sniffed by magic).
await copyFile(
  path.join(fixturesRoot, "mcpl", "ssw_conversion", "ssw2mcpl_expected.mcpl"),
  path.join(dst, "ssw2mcpl_expected.mcpl"),
);
await copyFile(
  path.join(fixturesRoot, "mcpl", "ssw_conversion", "reference.w"),
  path.join(dst, "reference.w"),
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
    `meshtal_sample.txt, xsdir_sample.txt, mctal_sample.mctal, chain_simple.xml, deck_minimal.txt, deck_l3.txt, ` +
    `serpent_res_sample.m, serpent_dep_sample.m, serpent_det_sample.m, ` +
    `alara_output_sample.out, fispact_inventory_sample.fis, usrbin_sample.lis, ` +
    `ssw2mcpl_expected.mcpl, reference.w, ` +
    `${figures.length} theory figures`,
);
