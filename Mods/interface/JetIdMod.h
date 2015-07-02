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

    Bool_t      GetInputIsCorrected() const         { return fCorrectedInput; }
    Bool_t      GetUseL1Correction() const          { return fCorrections.TestBit(mithep::Jet::L1); }
    Bool_t      GetUseL2Correction() const          { return fCorrections.TestBit(mithep::Jet::L2); }
    Bool_t      GetUseL3Correction() const          { return fCorrections.TestBit(mithep::Jet::L3); }
    Double_t    GetJetEEMFractionMinCut() const     { return fJetEEMFractionMinCut; }
    Double_t    GetMinChargedHadronFraction() const { return fMinChargedHadronFraction; }
    Double_t    GetMaxChargedHadronFraction() const { return fMaxChargedHadronFraction; }
    Double_t    GetMinNeutralHadronFraction() const { return fMinNeutralHadronFraction; }
    Double_t    GetMaxNeutralHadronFraction() const { return fMaxNeutralHadronFraction; }
    Double_t    GetMinChargedEMFraction() const     { return fMinChargedEMFraction; }
    Double_t    GetMaxChargedEMFraction() const     { return fMaxChargedEMFraction; }
    Double_t    GetMinNeutralEMFraction() const     { return fMinNeutralEMFraction; }
    Double_t    GetMaxNeutralEMFraction() const     { return fMaxNeutralEMFraction; }
    Bool_t      GetApplyBetaCut() const             { return fApplyBetaCut; }
    Bool_t      GetApplyPFLooseId() const           { return fApplyPFLooseId; }
    UInt_t      GetMVATrainingSet() const           { return fMVATrainingSet; }
    UInt_t      GetMVACutWP() const                 { return fMVACutWP; }
    char const* GetMVAWeightsFile() const           { return fMVAWeightsFile; }
    char const* GetMVACutsFile() const              { return fMVACutsFile; }
    
    JetIDMVA* GetJetIDMVA() const { return fJetIDMVA; }

    void SetInputIsCorrected(Bool_t b)           { fCorrectedInput = b; }
    void SetUseL1Correction(Bool_t b)            { fCorrections.SetBit(mithep::Jet::L1, b); }
    void SetUseL2Correction(Bool_t b)            { fCorrections.SetBit(mithep::Jet::L2, b); }
    void SetUseL3Correction(Bool_t b)            { fCorrections.SetBit(mithep::Jet::L3, b); }
    void SetJetEEMFractionMinCut(Double_t cut)   { fJetEEMFractionMinCut = cut; }
    void SetMinChargedHadronFraction(Double_t m) { fMinChargedHadronFraction = m; }
    void SetMaxChargedHadronFraction(Double_t m) { fMaxChargedHadronFraction = m; }
    void SetMinNeutralHadronFraction(Double_t m) { fMinNeutralHadronFraction = m; }
    void SetMaxNeutralHadronFraction(Double_t m) { fMaxNeutralHadronFraction = m; }
    void SetMinChargedEMFraction(Double_t m)     { fMinChargedEMFraction = m; }
    void SetMaxChargedEMFraction(Double_t m)     { fMaxChargedEMFraction = m; }
    void SetMinNeutralEMFraction(Double_t m)     { fMinNeutralEMFraction = m; }
    void SetMaxNeutralEMFraction(Double_t m)     { fMaxNeutralEMFraction = m; }
    void SetApplyBetaCut(Bool_t b)               { fApplyBetaCut = b; }
    void SetApplyPFLooseId(Bool_t b)             { fApplyPFLooseId = b; }
    void SetMVATrainingSet(UInt_t s)             { fMVATrainingSet = s; }
    void SetMVACutWP(UInt_t w)                   { fMVACutWP = w; }
    void SetMVAWeightsFile(char const* n)        { fMVAWeightsFile = n; }
    void SetMVACutsFile(char const* n)           { fMVACutsFile = n; }

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

    Bool_t   fCorrectedInput = kTRUE;
    BitMask8 fCorrections = BitMask8(7); // default = L1 + L2 + L3
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
    UInt_t   fMVATrainingSet = JetIDMVA::nMVATypes; //JetIDMVA::MVAType
    UInt_t   fMVACutWP = JetIDMVA::kLoose; //JetIDMVA::CutType
    TString  fMVAWeightsFile = "";
    TString  fMVACutsFile = "";

    JetIDMVA* fJetIDMVA = 0;

    ClassDef(JetIdMod, 0) // Jet identification module
  };
}

#endif
