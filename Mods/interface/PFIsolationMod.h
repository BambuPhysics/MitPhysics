//--------------------------------------------------------------------------------------------------
// PFIsolationMod
//
// This module runs the IsolationTools and exports the isolation sums as arrays. Just a convenient
// tool to avoid running the loop-heavy isolation too many times.
// The module, IsolationTools::PFIsoDeposit, and related functions were written to replace the
// Electron isolation variables (recomputing isolation sums on the fly) but the CMSSW values could
// not be reproduced exactly, in particular for CH iso. Therefore the isolation variables were brought
// back into the Electron object, and this module is not used at the moment.
//
// Authors: Y.Iiyama
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PFISOLATIONMOD_H
#define MITPHYSICS_MODS_PFISOLATIONMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataCont/interface/Types.h"

#include "TH1D.h"

namespace mithep {

  class PFIsolationMod : public BaseMod {
  public:
    enum IsolationType {
      kChargedHadron,
      kNeutralHadron,
      kPhoton,
      nIsolationTypes
    };

    PFIsolationMod(char const* name = "PFIsolationMod", char const* title = "PF isolation exporter");
    ~PFIsolationMod() {}

    char const* GetChargedHadronIsoName() const { return fOutputName + "CH"; }
    char const* GetNeutralHadronIsoName() const { return fOutputName + "NH"; }
    char const* GetPhotonIsoName() const { return fOutputName + "Ph"; }

    void SetInputName(char const* n) { fCandsName = n; }
    void SetPFCandidatesName(char const* n) { fPFCandidatesName = n; }
    void SetPFNoPileupCandidatesName(char const* n) { fPFNoPileupCandidatesName = n; }
    void SetOutputName(char const* n) { fOutputName = n; }

    void SetMinDR(IsolationType type, Double_t dr) { fMinDR[type] = dr; }
    void SetMaxDR(IsolationType type, Double_t dr) { fMaxDR[type] = dr; }
    void SetMaxDRrho(IsolationType type, Double_t drho) { fMaxDRho[type] = drho; }

  protected:
    void SlaveBegin() override;
    void SlaveTerminate() override;
    void Process() override;

    TString fCandsName;
    TString fPFCandidatesName{"PFCandidates"};
    TString fPFNoPileupCandidatesName{"PFNoPUCandidates"};
    TString fOutputName;

    NFArrDouble fIsoSum[nIsolationTypes];

    Double_t fMinDR[nIsolationTypes]{};
    Double_t fMaxDR[nIsolationTypes]{};
    Double_t fMaxDRho[nIsolationTypes]{};
    Double_t fVetoDRBarrel[nIsolationTypes]{};
    Double_t fVetoDREndcap[nIsolationTypes]{};

    TH1D* fHIso[nIsolationTypes]{};

    ClassDef(PFIsolationMod, 0)
  };

}

#endif
