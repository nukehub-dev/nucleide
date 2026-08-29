import { FlaskConical, Globe, Rss, MessageSquare } from "lucide-react";
import { type NavItem } from "@nukehub/docs-kit";

export const navItems: NavItem[] = [
  {
    title: "NukeHub",
    icon: Globe,
    url: "https://nukehub.org",
    newpage: true,
  },
  {
    title: "NukeBlog",
    icon: Rss,
    url: "https://blog.nukehub.org",
    newpage: true,
  },
  {
    title: "NukeTalk",
    icon: MessageSquare,
    url: "https://talk.nukehub.org",
    newpage: true,
  },
  {
    title: "PyPI",
    icon: FlaskConical,
    url: "https://pypi.org/project/nucleide/",
    newpage: true,
  },
];
