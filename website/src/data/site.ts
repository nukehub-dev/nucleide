import { Logo } from "../components/Logo";
import { type SiteConfig } from "@nukehub/docs-kit";

// Same artwork as Logo.tsx / favicon.svg, scaled from the 0-100 viewBox into
// the 1024x1024 canvas used by the kit's dynamic favicon.
const FAVICON_PATHS =
  '<g fill="currentColor" transform="scale(10.24)">' +
  '<rect x="2.5" y="2.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="22.5" y="2.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="42.5" y="2.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="62.5" y="2.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="2.5" y="22.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="22.5" y="22.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="42.5" y="22.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="82.5" y="22.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="2.5" y="42.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="22.5" y="42.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="62.5" y="42.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="82.5" y="42.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="2.5" y="62.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="42.5" y="62.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="62.5" y="62.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="82.5" y="62.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="22.5" y="82.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="42.5" y="82.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="62.5" y="82.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="82.5" y="82.5" width="15" height="15" rx="3" opacity="0.18" />' +
  '<rect x="2.5" y="82.5" width="15" height="15" rx="3" opacity="0.7" />' +
  '<rect x="22.5" y="62.5" width="15" height="15" rx="3" opacity="0.85" />' +
  '<rect x="62.5" y="22.5" width="15" height="15" rx="3" opacity="0.85" />' +
  '<rect x="82.5" y="2.5" width="15" height="15" rx="3" opacity="0.7" />' +
  '<path fill-rule="evenodd" d="M45.5 42.5 h9 a3 3 0 0 1 3 3 v9 a3 3 0 0 1 -3 3 h-9 a3 3 0 0 1 -3 -3 v-9 a3 3 0 0 1 3 -3 z M45.8 50 a4.2 4.2 0 1 1 8.4 0 a4.2 4.2 0 1 1 -8.4 0 z" />' +
  "</g>";

export const SITE: SiteConfig = {
  name: "Nucleide",
  logoText: "Nucleide",
  description: "A modern Rust toolkit for nuclear-engineering workflow glue.",
  site: "https://nukehub-dev.github.io",
  base: "/nucleide",
  github: "https://github.com/nukehub-dev/nucleide",
  editBranch: "main",
  editPath: "docs/",
  logo: Logo,
  faviconPaths: FAVICON_PATHS,
};
