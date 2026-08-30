import { test, expect } from "@playwright/test";

async function waitForWasmReady(page) {
  // Pages can host several WASM demos, each rendering its own loader and
  // hydrating on visibility (client:visible). Scroll the page end-to-end so
  // every demo hydrates, then wait until no loader remains.
  const loaders = page.locator("text=Loading Nucleide WASM…");
  await loaders.first().waitFor({ state: "attached" });
  try {
    await loaders.first().scrollIntoViewIfNeeded();
  } catch {
    // The loader may have already hydrated and been removed; that's fine.
  }
  await page.evaluate(() => window.scrollTo(0, document.body.scrollHeight));
  await expect(loaders).toHaveCount(0, { timeout: 15_000 });
  await page.evaluate(() => window.scrollTo(0, 0));
}

async function assertNoKatexErrors(page) {
  await expect(page.locator(".katex-error")).toHaveCount(0);
}

async function assertNoWasmError(page) {
  await expect(page.getByText("WASM error:", { exact: false })).not.toBeVisible();
}

const INTERACTIVE_PAGES = [
  {
    path: "tutorials/interactive/nuclides",
    button: "Look up",
    output: "text=cinder",
  },
  {
    path: "tutorials/interactive/materials",
    button: "Compute fractions",
    output: "text=Atom fractions",
    // A known data cell: guards against map-shaped WASM values rendering as
    // empty tables (Object.entries on a JS Map yields no rows).
    cell: "H1",
    chart: { selector: ".js-plotly-plot" },
  },
  {
    path: "tutorials/interactive/enrichment",
    button: "Solve cascade",
    output: "text=Enriching stages",
    chart: { button: "Solve cascade", selector: ".js-plotly-plot" },
  },
  {
    path: "tutorials/interactive/depletion",
    button: "Deplete",
    output: "text=Atom count",
    chart: { button: "Burnup curve", selector: ".js-plotly-plot" },
  },
  {
    path: "tutorials/interactive/mcnp-io",
    button: "Parse",
    output: "text=Material 1",
    chart: { actions: ["meshtal", "Load sample meshtal", "Parse"], selector: ".js-plotly-plot" },
  },
  {
    path: "tutorials/interactive/variance-reduction",
    button: "Generate MAGIC bounds",
    output: "text=Groups per voxel:",
    chart: { selector: ".js-plotly-plot" },
  },
];

for (const { path, button, output, cell, chart } of INTERACTIVE_PAGES) {
  test(`interactive demo: ${path}`, async ({ page }) => {
    await page.goto(path);
    await waitForWasmReady(page);
    await assertNoKatexErrors(page);

    await page.getByRole("button", { name: button }).click();

    await assertNoWasmError(page);
    await expect(page.locator(output).first()).toBeVisible();
    if (cell) {
      await expect(page.locator("tbody td", { hasText: cell }).first()).toBeVisible();
    }
    if (path === "tutorials/interactive/materials") {
      // The compendium browser hydrates on visibility (client:visible), then
      // fetches public/data/MaterialsCompendium.json and parses it with WASM.
      await page
        .getByRole("heading", { name: "Browse the Materials Compendium" })
        .scrollIntoViewIfNeeded();
      await expect(page.getByText("materials loaded")).toBeVisible();
    }

    if (chart) {
      if (chart.actions) {
        for (const action of chart.actions) {
          await page.getByRole("button", { name: action }).click();
        }
      } else if (chart.button) {
        await page.getByRole("button", { name: chart.button }).click();
      }
      await assertNoWasmError(page);
      await expect(page.locator(chart.selector).first()).toBeVisible({ timeout: 10_000 });
    }
  });
}

test("theory page has no KaTeX errors", async ({ page }) => {
  await page.goto("theory/variance-reduction");
  await assertNoKatexErrors(page);
});
