#include "MitPhysics/SelMods/interface/BadEventsFilterMod.h"

#include "MitAna/DataTree/interface/EvtSelData.h"

ClassImp(mithep::BadEventsFilterMod)

void
mithep::BadEventsFilterMod::Process()
{
  auto* evtSelData = GetObject<mithep::EvtSelData>(fEvtSelDataName);
  if (!evtSelData) {
    SendError(kWarning, "Process", "EvtSelData " + fEvtSelDataName + " not found.");
    return;
  }

  if ((fBitMask & evtSelData->metFiltersWord()) != fBitMask)
    SkipEvent();

  return;
}

void
mithep::BadEventsFilterMod::SetMask(EvtSelFilter f, Bool_t b)
{
  if (b)
    fBitMask |= (1 << f);
  else
    fBitMask &= ~(1 << f);
}
