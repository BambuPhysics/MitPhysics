//--------------------------------------------------------------------------------------------------
// BadEventsFilterMod
// 
// Filter events according to the EvtSelData (i.e. MET filters).
// Implementation of MET filters is not ideal at the moment and might change in the future.
//
// Authors: Y.Iiyama
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_SELMODS_BADEVENTSFILTERMOD_H
#define MITPHYSICS_SELMODS_BADEVENTSFILTERMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "TH1D.h"

namespace mithep {

  class BadEventsFilterMod : public BaseMod {
  public:
    BadEventsFilterMod(char const* name = "BadEventsFilterMod", char const* title = "MET filters") : BaseMod(name, title) {}
    ~BadEventsFilterMod() {}

    // Must be identical to MitProd/TreeFiller/interface/FillerEvtSelData
    // -> One reason why the implementation is poor. At least this enum should be defined in MitAna/DataTree.
    enum EvtSelFilter {
      kHBHENoiseFilter,
      kECALDeadCellFilter,
      kTrackingFailureFilter,
      kEEBadScFilter,
      kECALaserCorrFilter,
      kManyStripClusFilter,
      kTooManyStripClusFilter,
      kLogErrorTooManyClustersFilter,
      kCSCTightHaloFilter,
      kCSCLooseHaloFilter,
      nEvtSelFilters
    };

    void SetInputName(char const* n) { fEvtSelDataName = n; }

    void SetHBHENoiseFilter(Bool_t b = kTRUE) { SetMask(kHBHENoiseFilter, b); }
    void SetECALDeadCellFilter(Bool_t b = kTRUE) { SetMask(kECALDeadCellFilter, b); }
    void SetTrackingFailureFilter(Bool_t b = kTRUE) { SetMask(kTrackingFailureFilter, b); }
    void SetEEBadScFilter(Bool_t b = kTRUE) { SetMask(kEEBadScFilter, b); }
    void SetECALaserCorrFilter(Bool_t b = kTRUE) { SetMask(kECALaserCorrFilter, b); }
    void SetManyStripClusFilter(Bool_t b = kTRUE) { SetMask(kManyStripClusFilter, b); }
    void SetTooManyStripClusFilter(Bool_t b = kTRUE) { SetMask(kTooManyStripClusFilter, b); }
    void SetLogErrorTooManyClustersFilter(Bool_t b = kTRUE) { SetMask(kLogErrorTooManyClustersFilter, b); }
    void SetCSCTightHaloFilter(Bool_t b = kTRUE) { SetMask(kCSCTightHaloFilter, b); }
    void SetCSCLooseHaloFilter(Bool_t b = kTRUE) { SetMask(kCSCLooseHaloFilter, b); }

  protected:
    void SlaveBegin() override;
    void Process() override;

    void SetMask(EvtSelFilter, Bool_t);

    TString fEvtSelDataName{"EvtSelData"};
    Int_t fBitMask{0};

    TH1D* hCounter{0};

    ClassDef(BadEventsFilterMod, 0)
  };

}

#endif
