import { defineConfig, devices } from "@playwright/test";

// Distinctive port: 4321 is Astro's default, so an unrelated project's dev
// server can squat on it and be silently reused (`reuseExistingServer`),
// failing every test with 404s. 4387 is unique to this repo.
const e2eBase = "http://localhost:4387";

export default defineConfig({
  testDir: "./e2e",
  fullyParallel: true,
  forbidOnly: !!process.env.CI,
  retries: process.env.CI ? 2 : 0,
  workers: process.env.CI ? 1 : undefined,
  reporter: process.env.CI ? "list" : "html",
  use: {
    baseURL: e2eBase,
    trace: "on-first-retry",
  },
  projects: [
    {
      name: "chromium",
      use: { ...devices["Desktop Chrome"] },
    },
  ],
  webServer: {
    command: "npm run preview -- --port 4387",
    // Use Astro's preview child mode to suppress automatic detachment in
    // agent sessions. Playwright owns the foreground process lifecycle.
    env: { ASTRO_PREVIEW_BACKGROUND: "1" },
    url: e2eBase,
    timeout: 120_000,
    reuseExistingServer: !process.env.CI,
  },
});
