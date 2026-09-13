import { NuclideExplorer } from "./NuclideExplorer";
import { MaterialBuilder } from "./MaterialBuilder";
import { CompendiumBrowser } from "./CompendiumBrowser";
import { CascadeSolver } from "./CascadeSolver";
import { DepletionStep } from "./DepletionStep";
import { McnpParser } from "./McnpParser";
import { SerpentParser } from "./SerpentParser";
import { FlukaParser } from "./FlukaParser";
import { DeckEditor } from "./DeckEditor";
import { VrDemo } from "./VrDemo";
import { ActivationDemo } from "./ActivationDemo";
import { DeterministicDemo } from "./DeterministicDemo";
import { EmitterDemo } from "./EmitterDemo";
import { KineticsTransient } from "./KineticsTransient";
import { McplDemo } from "./McplDemo";
import { SpectroscopyDemo } from "./SpectroscopyDemo";
import { UqDemo } from "./UqDemo";

interface InteractiveDemoProps {
  kind:
    | "nuclides"
    | "materials"
    | "materials-compendium"
    | "enrichment"
    | "depletion"
    | "mcnp-io"
    | "serpent-io"
    | "fluka-io"
    | "deck-editor"
    | "variance-reduction"
    | "activation"
    | "deterministic"
    | "emitter"
    | "kinetics"
    | "mcpl-io"
    | "spectroscopy"
    | "uq";
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
    case "serpent-io":
      return <SerpentParser />;
    case "fluka-io":
      return <FlukaParser />;
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
    case "mcpl-io":
      return <McplDemo />;
    case "spectroscopy":
      return <SpectroscopyDemo />;
    case "uq":
      return <UqDemo />;
    default:
      return <div className="text-sm text-muted-foreground">Unknown demo kind: {kind}</div>;
  }
}
