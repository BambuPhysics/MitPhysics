#include "MitPhysics/SelMods/interface/BadEventsFilterMod.h"

#include "MitAna/DataTree/interface/EvtSelData.h"

ClassImp(mithep::BadEventsFilterMod)

void
mithep::BadEventsFilterMod::SlaveBegin()
{
  if (GetFillHist()) {
    AddTH1(hCounter, "hMETFilterCounter", "Number of events flagged bad", nEvtSelFilters, 0., double(nEvtSelFilters));
  }
}

void
mithep::BadEventsFilterMod::Process()
{
  auto* evtSelData = GetObject<mithep::EvtSelData>(fEvtSelDataName);
  if (!evtSelData) {
    SendError(kWarning, "Process", "EvtSelData " + fEvtSelDataName + " not found.");
    return;
  }

  if (GetFillHist()) {
    Int_t w = evtSelData->metFiltersWord();
    for (unsigned iF = 0; iF != nEvtSelFilters; ++iF) {
      if ((w & 1) == 0)
        hCounter->Fill(iF + 0.5);
      w >>= 1;
    }
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
