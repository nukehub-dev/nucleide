import { test, expect } from "@playwright/test";

async function waitForWasmReady(page) {
  const loader = page.locator("text=Loading Nucleide WASM…");
  await loader.waitFor({ state: "attached" });
  try {
    await loader.scrollIntoViewIfNeeded();
  } catch {
    // The loader may have already hydrated and been removed; that's fine.
  }
  await loader.waitFor({ state: "detached" });
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
  },
  {
    path: "tutorials/interactive/enrichment",
    button: "Solve cascade",
    output: "text=Enriching stages",
  },
  {
    path: "tutorials/interactive/depletion",
    button: "Deplete",
    output: "text=Atom count",
  },
  {
    path: "tutorials/interactive/mcnp-io",
    button: "Parse",
    output: "text=Material 1",
  },
  {
    path: "tutorials/interactive/variance-reduction",
    button: "Generate MAGIC bounds",
    output: "text=Groups per voxel:",
  },
];

for (const { path, button, output, cell } of INTERACTIVE_PAGES) {
  test(`interactive demo: ${path}`, async ({ page }) => {
    await page.goto(path);
    await waitForWasmReady(page);
    await assertNoKatexErrors(page);

    await page.getByRole("button", { name: button }).click();

    await assertNoWasmError(page);
    await expect(page.locator(output).first()).toBeVisible();
    if (cell) {
      await expect(
        page.locator("tbody td", { hasText: cell }).first(),
      ).toBeVisible();
    }
  });
}

test("theory page has no KaTeX errors", async ({ page }) => {
  await page.goto("theory/variance-reduction");
  await assertNoKatexErrors(page);
});
