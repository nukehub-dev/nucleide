import { NuclideExplorer } from "./NuclideExplorer";
import { MaterialBuilder } from "./MaterialBuilder";
import { CompendiumBrowser } from "./CompendiumBrowser";
import { CascadeSolver } from "./CascadeSolver";
import { DepletionStep } from "./DepletionStep";
import { McnpParser } from "./McnpParser";
import { DeckEditor } from "./DeckEditor";
import { VrDemo } from "./VrDemo";
import { ActivationDemo } from "./ActivationDemo";
import { DeterministicDemo } from "./DeterministicDemo";
import { EmitterDemo } from "./EmitterDemo";
import { KineticsTransient } from "./KineticsTransient";
import { SpectroscopyDemo } from "./SpectroscopyDemo";

interface InteractiveDemoProps {
  kind:
    | "nuclides"
    | "materials"
    | "materials-compendium"
    | "enrichment"
    | "depletion"
    | "mcnp-io"
    | "deck-editor"
    | "variance-reduction"
    | "activation"
    | "deterministic"
    | "emitter"
    | "kinetics"
    | "spectroscopy";
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
    case "deck-editor":
      return <DeckEditor />;
    case "variance-reduction":
      return <VrDemo />;
    case "activation":
      return <ActivationDemo />;
    case "deterministic":
      return <DeterministicDemo />;
    case "emitter":
      return <EmitterDemo />;
    case "kinetics":
      return <KineticsTransient />;
    case "spectroscopy":
      return <SpectroscopyDemo />;
    default:
      return <div className="text-sm text-muted-foreground">Unknown demo kind: {kind}</div>;
  }
}
