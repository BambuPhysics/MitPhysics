#include "MitPhysics/SelMods/interface/EventCategoryMod.h"
#include "MitAna/DataCont/interface/Types.h"

ClassImp(mithep::EventCategoryMod)

void
mithep::EventCategoryMod::Process()
{
  auto* selectionBits = GetObject<mithep::NFArrBool>(fInputName);

  if (!selectionBits) {
    SendError(kAbortAnalysis, "Process", "No selection bits named " + fInputName + " found.");
    return;
  }

  if (selectionBits->GetEntries() <= fCategory || !selectionBits->At(fCategory))
    SkipEvent();
}
