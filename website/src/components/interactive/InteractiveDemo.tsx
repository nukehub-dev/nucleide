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
import { TritiumBreakthrough } from "./TritiumBreakthrough";
import { McplDemo } from "./McplDemo";
import { SpectroscopyDemo } from "./SpectroscopyDemo";
import { UqDemo } from "./UqDemo";
import { DamageDemo } from "./DamageDemo";
import { FusionSourceDemo } from "./FusionSourceDemo";
import { EcrhDemo } from "./EcrhDemo";
import { BlanketDemo } from "./BlanketDemo";
import { EquilibDemo } from "./EquilibDemo";
import { UnfoldDemo } from "./UnfoldDemo";
import { ClearanceDemo } from "./ClearanceDemo";
import { SubletDemo } from "./SubletDemo";

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
    | "tritium"
    | "mcpl-io"
    | "spectroscopy"
    | "uq"
    | "damage"
    | "fusion-sources"
    | "ecrh"
    | "blanket"
    | "equilibrium"
    | "unfold"
    | "clearance"
    | "sublet";
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
    case "tritium":
      return <TritiumBreakthrough />;
    case "mcpl-io":
      return <McplDemo />;
    case "spectroscopy":
      return <SpectroscopyDemo />;
    case "uq":
      return <UqDemo />;
    case "damage":
      return <DamageDemo />;
    case "fusion-sources":
      return <FusionSourceDemo />;
    case "ecrh":
      return <EcrhDemo />;
    case "blanket":
      return <BlanketDemo />;
    case "equilibrium":
      return <EquilibDemo />;
    case "unfold":
      return <UnfoldDemo />;
    case "clearance":
      return <ClearanceDemo />;
    case "sublet":
      return <SubletDemo />;
    default:
      return <div className="text-sm text-muted-foreground">Unknown demo kind: {kind}</div>;
  }
}
