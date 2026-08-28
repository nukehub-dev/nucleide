import { FlaskConical } from "lucide-react";
import { GitHubIcon, type FooterColumn, type FooterLink } from "@nukehub/docs-kit";

export const footerColumns: FooterColumn[] = [
  {
    title: "Project",
    links: [
      {
        title: "Source code",
        url: "https://github.com/nukehub-dev/nucleide",
        icon: GitHubIcon,
        newpage: true,
      },
      {
        title: "PyPI",
        url: "https://pypi.org/project/nucleide/",
        icon: FlaskConical,
        newpage: true,
      },
    ],
  },
  {
    title: "Community",
    links: [
      {
        title: "NukeHub",
        url: "https://nukehub.org",
        newpage: true,
      },
      {
        title: "NukeBlog",
        url: "https://blog.nukehub.org",
        newpage: true,
      },
      {
        title: "NukeTalk",
        url: "https://talk.nukehub.org",
        newpage: true,
      },
    ],
  },
];

export const footerLegal: FooterLink[] = [
  {
    title: "License",
    url: "https://github.com/nukehub-dev/nucleide/blob/main/LICENSE",
    newpage: true,
  },
];
