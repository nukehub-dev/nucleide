import { NuclideExplorer } from "./NuclideExplorer";
import { MaterialBuilder } from "./MaterialBuilder";
import { CompendiumBrowser } from "./CompendiumBrowser";
import { CascadeSolver } from "./CascadeSolver";
import { DepletionStep } from "./DepletionStep";
import { McnpParser } from "./McnpParser";
import { VrDemo } from "./VrDemo";

interface InteractiveDemoProps {
  kind:
    | "nuclides"
    | "materials"
    | "materials-compendium"
    | "enrichment"
    | "depletion"
    | "mcnp-io"
    | "variance-reduction";
}

export function InteractiveDemo({ kind }: InteractiveDemoProps) {
  switch (kind) {
    case "nuclides":
      return <NuclideExplorer />;
    case "materials":
      return <MaterialBuilder />;
    case "materials-compendium":
      return <CompendiumBrowser />;
    case "enrichment":
      return <CascadeSolver />;
    case "depletion":
      return <DepletionStep />;
    case "mcnp-io":
      return <McnpParser />;
    case "variance-reduction":
      return <VrDemo />;
    default:
      return <div className="text-sm text-muted-foreground">Unknown demo kind: {kind}</div>;
  }
}
