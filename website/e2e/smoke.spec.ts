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

interface ExtraStep {
  button: string;
  output?: string;
  fill?: { label: string; text: string };
  chart?: string;
}

interface InteractivePage {
  path: string;
  button: string;
  output: string;
  cell?: string;
  paste?: string;
  chart?: { button?: string; actions?: string[]; selector: string };
  extraSteps?: ExtraStep[];
}

const INTERACTIVE_PAGES: InteractivePage[] = [
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
    extraSteps: [
      { button: "Mix", output: "text=Mixed atom fractions" },
      { button: "To XML", output: "text=density" },
      {
        button: "Compute dose",
        fill: { label: "Formula", text: "U" },
        output: "text=Dose per gram",
      },
    ],
  },
  {
    path: "tutorials/interactive/enrichment",
    button: "Solve cascade",
    output: "text=Enriching stages",
    chart: { button: "Solve cascade", selector: ".js-plotly-plot" },
    extraSteps: [{ button: "Optimize M*", output: "text=Enriching stages" }],
  },
  {
    path: "tutorials/interactive/depletion",
    button: "Deplete",
    output: "text=Atom count",
    chart: { button: "Burnup curve", selector: ".js-plotly-plot" },
    extraSteps: [{ button: "Load sample chain" }, { button: "Deplete", output: "text=Atom count" }],
  },
  {
    path: "tutorials/interactive/mcnp-io",
    button: "Parse",
    output: "text=Material 1",
    chart: {
      actions: ["meshtal", "Load sample meshtal", "Parse", "xsdir", "Load sample xsdir", "Parse"],
      selector: ".js-plotly-plot",
    },
    extraSteps: [{ button: "wwinp" }, { button: "Parse", output: "text=ni:" }],
  },
  {
    path: "tutorials/interactive/serpent-io",
    button: "Parse",
    output: "text=Variables:",
    chart: { button: "Parse", selector: ".js-plotly-plot" },
    extraSteps: [
      { button: "dep" },
      { button: "Load sample dep" },
      { button: "Parse", output: "text=Nuclides:" },
      { button: "det" },
      { button: "Load sample det" },
      { button: "Parse", output: "text=Detectors:", chart: ".js-plotly-plot" },
    ],
  },
  {
    path: "tutorials/interactive/fluka-io",
    button: "Parse",
    output: "text=Tallies:",
    chart: { actions: ["Load sample USRBIN", "Parse"], selector: ".js-plotly-plot" },
  },
  {
    path: "tutorials/interactive/variance-reduction",
    button: "Generate MAGIC bounds",
    output: "text=Groups per voxel:",
    chart: { button: "Generate MAGIC bounds", selector: ".js-plotly-plot" },
    extraSteps: [
      { button: "Sample index", output: "text=Sampled index:" },
      { button: "Sample voxel", output: "text=index=" },
    ],
  },
  {
    path: "tutorials/interactive/activation",
    button: "Parse",
    output: "text=inner_mix",
    paste: `geometry rectangular
mat_loading
inner_zone inner_mix
end
mixture inner_mix
material WATER 1.0 1.0
end
flux flux_1 data/fluxin1 1.0 1 default
schedule 1_year
1 y flux_1 steady_state 0 s
end
pulsehistory steady_state
1 0 s
end
cooling
1 d
end`,
    extraSteps: [
      { button: "ALARA output" },
      { button: "Load sample" },
      { button: "Parse", output: "text=Rows:", chart: ".js-plotly-plot" },
      { button: "FISPACT output" },
      { button: "Load sample" },
      { button: "Parse", output: "text=Rows:", chart: ".js-plotly-plot" },
      { button: "ORIGEN TAPE5" },
      { button: "Parse", output: "text=Materials:" },
      { button: "ORIGEN TAPE6" },
      { button: "Parse", output: "text=total activity:" },
      { button: "ORIGEN TAPE9" },
      { button: "Parse", output: "text=Entries:" },
      { button: "R2S workflow" },
      { button: "Parse", output: "text=Top schedule:" },
      { button: "R2S snapshot" },
      { button: "Parse", output: "text=Top schedule:" },
      { button: "Parse", output: "text=Decks:" },
    ],
  },
  {
    path: "tutorials/interactive/deterministic",
    button: "Parse",
    output: "text=U235",
    paste: `ISOTXS 2
NUCLIDE U235 92235 2
1.1 2.2
NUCLIDE PU239 94239 2
4.4 5.5`,
  },
  {
    path: "tutorials/interactive/deck-editor",
    button: "Parse",
    output: "text=Cell 1",
    extraSteps: [
      { button: "Validate", output: "text=✓ valid" },
      { button: "Set density", output: "text=Cell 1" },
      { button: "Load L3 sample" },
      { button: "Parse", output: "text=Cell 1" },
    ],
  },
  {
    path: "tutorials/interactive/emitter",
    button: "Emit",
    output: "text=MCNP",
    extraSteps: [{ button: "ARMI keys" }, { button: "Emit", output: "text=Mass-drift report" }],
  },
  {
    path: "tutorials/interactive/kinetics",
    button: "Run transient",
    output: "text=Prompt jump",
    chart: { button: "Run transient", selector: ".js-plotly-plot" },
    extraSteps: [
      { button: "Six-group preset" },
      { button: "Run transient", output: "text=Final n" },
    ],
  },
  {
    path: "tutorials/interactive/spectroscopy",
    button: "Smooth spectrum",
    output: "text=Net counts",
    cell: "Gross counts",
    chart: { button: "Smooth spectrum", selector: ".js-plotly-plot" },
    extraSteps: [
      { button: "Synthetic peak" },
      { button: "Smooth spectrum", output: "text=Net counts" },
    ],
  },
];

for (const { path, button, output, cell, chart, paste, extraSteps } of INTERACTIVE_PAGES) {
  test(`interactive demo: ${path}`, async ({ page }) => {
    await page.goto(path);
    await waitForWasmReady(page);
    await assertNoKatexErrors(page);

    if (paste) {
      await page.locator("textarea").first().fill(paste);
    }
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

    if (extraSteps) {
      for (const step of extraSteps) {
        if (step.fill) {
          const input = page.getByLabel(step.fill.label);
          await input.scrollIntoViewIfNeeded();
          await input.fill(step.fill.text);
        }
        const stepButton = page.getByRole("button", { name: step.button, exact: true });
        await stepButton.scrollIntoViewIfNeeded();
        await stepButton.click();
        // Sample loaders disable their button while fetching; wait for the
        // fetch to land before the next step (instant for other buttons).
        await expect(stepButton).toBeEnabled({ timeout: 10_000 });
        await assertNoWasmError(page);
        if (step.output) {
          await expect(page.locator(step.output).first()).toBeVisible();
        }
        if (step.chart) {
          await expect(page.locator(step.chart).first()).toBeVisible({ timeout: 10_000 });
        }
      }
    }
  });
}

test("theory page has no KaTeX errors", async ({ page }) => {
  await page.goto("theory/variance-reduction");
  await assertNoKatexErrors(page);
});
