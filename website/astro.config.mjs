import { defineConfig } from "astro/config";
import react from "@astrojs/react";
import mdx from "@astrojs/mdx";
import sitemap from "@astrojs/sitemap";
import tailwindcss from "@tailwindcss/vite";
import remarkMath from "remark-math";
import rehypeKatex from "rehype-katex";
import { SITE } from "./src/data/site.ts";
import markdownNegotiation from "@nukehub/docs-kit/integrations/markdown-negotiation";

// https://astro.build/config
export default defineConfig({
  site: SITE.site,
  base: SITE.base,
  output: "static",
  integrations: [react(), mdx(), sitemap(), markdownNegotiation()],
  markdown: {
    remarkPlugins: [remarkMath],
    rehypePlugins: [rehypeKatex],
    shikiConfig: {
      themes: {
        light: "github-light",
        dark: "github-dark",
      },
      wrap: true,
    },
  },
  vite: {
    plugins: [tailwindcss()],
    build: {
      sourcemap: true,
    },
  },
});