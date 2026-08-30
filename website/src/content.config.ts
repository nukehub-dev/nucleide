import { defineCollection } from "astro:content";
import { glob } from "astro/loaders";
import { z } from "zod";

const referenceSchema = z.object({
  id: z.string(),
  title: z.string(),
  url: z.string().url(),
  source: z.string().optional(),
  date: z.string().optional(),
  authors: z.array(z.string()).optional(),
  type: z.enum(["article", "book", "inproceedings", "techreport", "misc"]).optional(),
  publisher: z.string().optional(),
  doi: z.string().optional(),
  arxiv: z.string().optional(),
  journal: z.string().optional(),
  volume: z.string().optional(),
  issue: z.string().optional(),
  pages: z.string().optional(),
});

const docs = defineCollection({
  loader: glob({ pattern: "**/*.{md,mdx}", base: "./src/content/docs" }),
  schema: z.object({
    title: z.string(),
    description: z.string().optional(),
    // Injected by nukehub-sync-docs: repo-relative source path for EditLink.
    editPath: z.string().optional(),
    sidebar: z
      .object({
        label: z.string().optional(),
        order: z.number().optional(),
      })
      .optional(),
    draft: z.boolean().optional(),
    references: z.array(referenceSchema).default([]),
  }),
});

export const collections = { docs };
