//--------------------------------------------------------------------------------------------------
// BadPFTrackFilterMod
// 
// Filter events if there is a bad PFMuon or PFChargedHadron. Ported from
// CMSSW cms-met/cmssw:CMSSW_8_0_X-METFilterUpdate to cope with 8_0_X PF ill-reconstruction.
//
// Authors: Y.Iiyama
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_SELMODS_BADPFTRACKFILTERMOD_H
#define MITPHYSICS_SELMODS_BADPFTRACKFILTERMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataCont/interface/Types.h"

namespace mithep {

  class BadPFTrackFilterMod : public BaseMod {
  public:
    BadPFTrackFilterMod(char const* name = "BadPFTackFilterMod", char const* title = "BadPFTrackFilter") : BaseMod(name, title) {}
    ~BadPFTrackFilterMod() {}

    void SetMuonsName(char const* n) { fMuonsName = n; }
    void SetPFCandidatesName(char const* n) { fPFCandidatesName = n; }
    void SetMinPt(Double_t m) { fMinPt = m; }
    void SetMinTrackRelErr(Double_t m) { fMinTrackRelErr = m; }
    void SetMinRelDPt(Double_t m) { fMinRelDPt = m; }
    void SetTaggingMode(Bool_t t) { fTaggingMode = t; }

    void SetOutputName(char const* n) { fOutput.SetName(n); }
    char const* GetOutputName(Int_t = -1) const;

  protected:
    void SlaveBegin() override;
    void SlaveTerminate() override;
    void Process() override;

    TString fMuonsName{Names::gkMuonBrn};
    TString fPFCandidatesName{Names::gkPFCandidatesBrn};
    Double_t fMinPt{100.};
    Double_t fMinTrackRelErr{0.5};
    Double_t fMinRelDPt{-0.5};

    Bool_t fTaggingMode{kFALSE};
    NFArrBool fOutput{};

    ClassDef(BadPFTrackFilterMod, 0)
  };

}

#endif
