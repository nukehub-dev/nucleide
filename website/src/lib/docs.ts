import type { CollectionEntry } from "astro:content";
import {
  buildSidebar as buildSidebarBase,
  getCategoryLabel,
  type SidebarNode,
  type SidebarSection,
} from "@nukehub/docs-kit/lib/docs";

/**
 * Custom category order for the Nucleide docs sidebar.
 *
 * The docs-kit default puts unknown categories at the end (order 99). We want
 * Theory to sit right after Tutorials and before Reference.
 */
const CATEGORY_ORDER: Record<string, number> = {
  home: 0,
  tutorials: 1,
  theory: 2,
  reference: 3,
  development: 4,
  architecture: 5,
  plan: 6,
  changelog: 7,
};

export function buildSidebar(docs: CollectionEntry<"docs">[]): SidebarSection[] {
  const sections = buildSidebarBase(docs);
  return sections.sort((a, b) => {
    const orderA = CATEGORY_ORDER[a.category] ?? 99;
    const orderB = CATEGORY_ORDER[b.category] ?? 99;
    return orderA - orderB;
  });
}

function flattenNodes(nodes: SidebarNode[]): { title: string; slug: string }[] {
  const result: { title: string; slug: string }[] = [];
  for (const node of nodes) {
    if (node.slug) {
      result.push({ title: node.title, slug: node.slug });
    }
    if (node.children) {
      result.push(...flattenNodes(node.children));
    }
  }
  return result;
}

export interface CommandPalettePage {
  id: string;
  title: string;
  url: string;
  description?: string;
  category: string;
}

export function getFlatPages(
  docs: CollectionEntry<"docs">[],
  nav: { title: string; url: string; newpage?: boolean }[],
): CommandPalettePage[] {
  const pages: CommandPalettePage[] = [];
  const sections = buildSidebar(docs);

  for (const section of sections) {
    for (const item of flattenNodes(section.items)) {
      pages.push({
        id: item.slug,
        title: item.title,
        url: item.slug === "index" ? "" : `${item.slug}/`,
        description: docs.find((d) => d.id.replace(/\.mdx?$/i, "") === item.slug)?.data.description,
        category: getCategoryLabel(section.category),
      });
    }
  }

  for (const item of nav) {
    if (item.newpage) continue;
    pages.push({
      id: `nav-${item.title}`,
      title: item.title,
      url: item.url,
      category: "Navigation",
    });
  }

  return pages;
}

export interface PaginationLink {
  slug: string;
  title: string;
}

export function getPrevNext(
  docs: CollectionEntry<"docs">[],
  currentSlug: string,
): { prev?: PaginationLink; next?: PaginationLink } {
  const sections = buildSidebar(docs);
  const flat: PaginationLink[] = [];

  for (const section of sections) {
    for (const item of flattenNodes(section.items)) {
      flat.push({ slug: item.slug, title: item.title });
    }
  }

  const idx = flat.findIndex((item) => item.slug === currentSlug);
  if (idx === -1) return {};

  return {
    prev: idx > 0 ? flat[idx - 1] : undefined,
    next: idx < flat.length - 1 ? flat[idx + 1] : undefined,
  };
}
