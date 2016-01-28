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
#include "MitAna/DataCont/interface/Types.h"
#include "TH1D.h"
#include "TObjArray.h"

namespace mithep {

  class BadEventsFilterMod : public BaseMod {
  public:
    BadEventsFilterMod(char const* name = "BadEventsFilterMod", char const* title = "MET filters");
    ~BadEventsFilterMod() {}

    void SetInputName(char const* n) { fEvtSelDataName = n; }
    void SetLabelTreeName(char const* n) { fLabelTreeName = n; }
    void SetLabelBranchName(char const* n) { fLabelBranchName = n; }
    void SetTaggingMode(Bool_t b = kTRUE) { fTaggingMode = b; }
    void SetInvertDecision(Bool_t b = kTRUE) { fNormalDecision = !b; }
    void SetOutputName(char const* n) { fTagResults.SetName(n); fFilterNames.SetName(TString(n) + "Names"); }

    void SetFilter(char const* name, Bool_t enable = kTRUE);
    void AddEventList(char const* name, char const* fileName);

    char const* GetOutputName() const { return fTagResults.GetName(); }

    struct EventID {
      UInt_t run;
      UInt_t lumi;
      UInt_t event;

      EventID(UInt_t r, UInt_t l, UInt_t e) : run(r), lumi(l), event(e) {}
      bool operator<(EventID const&) const;
    };

    typedef std::set<EventID> EventList;

  protected:
    void SlaveBegin() override;
    void Process() override;
    void BeginRun() override;
    Bool_t Notify() override;

    TString fEvtSelDataName{"EvtSelData"};
    TString fLabelTreeName{"EvtSelBits"};
    TString fLabelBranchName{"FilterLabels"};
    std::vector<std::string> fEnabledFilters{};
    std::map<std::string, EventList> fEventLists{};
    Int_t fBitMask{0};

    TObjArray fFilterNames;
    NFArrBool fTagResults;

    Bool_t fTaggingMode{kFALSE};
    Bool_t fNormalDecision{kTRUE};
    
    UInt_t fReload{0};

    TH1D* hCounter{0};

    ClassDef(BadEventsFilterMod, 0)
  };

}

inline
bool
mithep::BadEventsFilterMod::EventID::operator<(EventID const& rhs) const
{
  if (run < rhs.run)
    return true;

  if (run > rhs.run)
    return false;

  if (lumi < rhs.lumi)
    return true;

  if (lumi > rhs.lumi)
    return false;

  if (event < rhs.event)
    return true;
  else
    return false;
}

#endif
