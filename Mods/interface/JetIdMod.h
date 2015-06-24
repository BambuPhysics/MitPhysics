//--------------------------------------------------------------------------------------------------
// JetIdMod
//
// This module applies jet identification criteria and exports a pointer to a collection
// of "good jet" according to the specified identification scheme.
//
// Authors: S.Xie, Y.Iiyama, S.Narayanan
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_JETIDMOD_H
#define MITPHYSICS_MODS_JETIDMOD_H
#include "MitPhysics/Mods/interface/IdMod.h"

#include "MitAna/DataTree/interface/Jet.h"
#include "MitPhysics/Utils/interface/JetIDMVA.h"

namespace mithep {
  class JetIdMVA; 

  class JetIdMod : public IdMod<mithep::Jet> {
  public:
    JetIdMod(const char* name="JetIdMod",
             const char* title="Jet identification module");

    Bool_t   GetUseJetCorrection() const         { return fUseJetCorrection; }
    Double_t GetJetEEMFractionMinCut() const     { return fJetEEMFractionMinCut; }
    Double_t GetMinChargedHadronFraction() const { return fMinChargedHadronFraction; }
    Double_t GetMaxChargedHadronFraction() const { return fMaxChargedHadronFraction; }
    Double_t GetMinNeutralHadronFraction() const { return fMinNeutralHadronFraction; }
    Double_t GetMaxNeutralHadronFraction() const { return fMaxNeutralHadronFraction; }
    Double_t GetMinChargedEMFraction() const     { return fMinChargedEMFraction; }
    Double_t GetMaxChargedEMFraction() const     { return fMaxChargedEMFraction; }
    Double_t GetMinNeutralEMFraction() const     { return fMinNeutralEMFraction; }
    Double_t GetMaxNeutralEMFraction() const     { return fMaxNeutralEMFraction; }
    Bool_t   GetApplyBetaCut() const             { return fApplyBetaCut; }
    Bool_t   GetApplyPFLooseId() const           { return fApplyPFLooseId; }
    Bool_t   GetApplyMVACut() const              { return fApplyMVACut; }
    Bool_t   GetApplyMVACHS() const              { return fApplyMVACHS; }
    
    JetIDMVA* GetJetIDMVA() const { return fJetIDMVA; }

    void SetJetUseCorrection(Bool_t b)         { fUseJetCorrection = b; }
    void SetJetEEMFractionMinCut(Double_t cut) { fJetEEMFractionMinCut = cut; }
    void SetMinChargedHadronFraction(Double_t m) { fMinChargedHadronFraction = m; }
    void SetMaxChargedHadronFraction(Double_t m) { fMaxChargedHadronFraction = m; }
    void SetMinNeutralHadronFraction(Double_t m) { fMinNeutralHadronFraction = m; }
    void SetMaxNeutralHadronFraction(Double_t m) { fMaxNeutralHadronFraction = m; }
    void SetMinChargedEMFraction(Double_t m)     { fMinChargedEMFraction = m; }
    void SetMaxChargedEMFraction(Double_t m)     { fMaxChargedEMFraction = m; }
    void SetMinNeutralEMFraction(Double_t m)     { fMinNeutralEMFraction = m; }
    void SetMaxNeutralEMFraction(Double_t m)     { fMaxNeutralEMFraction = m; }
    void SetApplyBetaCut(Bool_t b)             { fApplyBetaCut = b; }
    void SetApplyPFLooseId(Bool_t b)           { fApplyPFLooseId = b; }
    void SetApplyMVACut(Bool_t b)              { fApplyMVACut = b; }
    void SetApplyMVACHS(Bool_t b)              { fApplyMVACHS = b; }

    void SetJetIDMVA(JetIDMVA* mva) { fJetIDMVA = mva; }

    void SetChargedHadronFractionRange(Double_t min, Double_t max)
    { fMinChargedHadronFraction = min; fMaxChargedHadronFraction = max; }
    void SetNeutralHadronFraction(Double_t min, Double_t max)
    { fMinNeutralHadronFraction = min; fMaxNeutralHadronFraction = max; }
    void SetChargedEMFraction(Double_t min, Double_t max)
    { fMinChargedEMFraction = min; fMaxChargedEMFraction = max; }
    void SetNeutralEMFraction(Double_t min, Double_t max)
    { fMinNeutralEMFraction = min; fMaxNeutralEMFraction = max; }

  protected:
    enum CutFlow {
      cAll,
      cEta,
      cPt,
      cEEMFraction,
      cChargedHFrac,
      cNeutralHFrac,
      cChargedEMFrac,
      cNeutralEMFrac,
      cPFLooseId,
      cBeta,
      cMVA,
      nCuts
    };

    Bool_t IsGood(mithep::Jet const&) override;
    void IdBegin() override;

    Bool_t   fUseJetCorrection = kTRUE;        //=true then use corrected energy
    Double_t fJetEEMFractionMinCut = 0.01;  //jet Eem fraction min cut for calo jets
    Double_t fMinChargedHadronFraction = 0.;
    Double_t fMaxChargedHadronFraction = 1.;
    Double_t fMinNeutralHadronFraction = 0.;
    Double_t fMaxNeutralHadronFraction = 1.;
    Double_t fMinChargedEMFraction = 0.;
    Double_t fMaxChargedEMFraction = 1.;
    Double_t fMinNeutralEMFraction = 0.;
    Double_t fMaxNeutralEMFraction = 1.;
    Bool_t   fApplyBetaCut = kFALSE;          //=true then apply beta cut
    Bool_t   fApplyPFLooseId = kFALSE;        //=true then apply PF loose ID
    Bool_t   fApplyMVACut = kFALSE;           //=true then apply MVA cut
    Bool_t   fApplyMVACHS = kFALSE;           //=true then apply MVA for CHS

    JetIDMVA* fJetIDMVA = 0;

    ClassDef(JetIdMod, 0) // Jet identification module
  };
}

#endif
