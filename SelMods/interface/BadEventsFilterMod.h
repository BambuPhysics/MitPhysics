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

    void SetInputName(char const* n) { fEvtSelDataName = n; }
    void SetLabelTreeName(char const* n) { fLabelTreeName = n; }
    void SetLabelBranchName(char const* n) { fLabelBranchName = n; }

    void SetFilter(char const* name, Bool_t enable = kTRUE);

  protected:
    void SlaveBegin() override;
    void Process() override;

    TString fEvtSelDataName{"EvtSelData"};
    TString fLabelTreeName{"EvtSelBits"};
    TString fLabelBranchName{"FilterLabels"};
    std::vector<std::string> fEnabledFilters{};
    Int_t fBitMask{0};

    TH1D* hCounter{0};

    ClassDef(BadEventsFilterMod, 0)
  };

}

#endif
