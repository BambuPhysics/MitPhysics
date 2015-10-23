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
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/JetIDMVA.h"

namespace mithep {
  class JetIdMVA; 

  class JetIdMod : public IdMod<mithep::Jet> {
  public:
    JetIdMod(const char* name="JetIdMod",
             const char* title="Jet identification module");

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
    Bool_t      GetPFId() const                     { return fPFId; }
    UInt_t      GetMVATrainingSet() const           { return fMVATrainingSet; }
    UInt_t      GetMVACutWP() const                 { return fMVACutWP; }
    char const* GetMVAWeightsFile(UInt_t idx = 0) const { return fMVAWeightsFile.at(idx); }
    char const* GetMVACutsFile() const              { return fMVACutsFile; }
    Bool_t      GetUseClassicBetaForMVA() const     { return fUseClassicBetaForMVA; }
    
    JetIDMVA* GetJetIDMVA() const { return fJetIDMVA; }

    void SetJetEEMFractionMinCut(Double_t cut)   { fJetEEMFractionMinCut = cut; }
    void SetMinChargedHadronFraction(Double_t m) { fMinChargedHadronFraction = m; }
    void SetMaxChargedHadronFraction(Double_t m) { fMaxChargedHadronFraction = m; }
    void SetMinNeutralHadronFraction(Double_t m) { fMinNeutralHadronFraction = m; }
    void SetMaxNeutralHadronFraction(Double_t m) { fMaxNeutralHadronFraction = m; }
    void SetMinChargedEMFraction(Double_t m)     { fMinChargedEMFraction = m; }
    void SetMaxChargedEMFraction(Double_t m)     { fMaxChargedEMFraction = m; }
    void SetMinNeutralEMFraction(Double_t m)     { fMinNeutralEMFraction = m; }
    void SetMaxNeutralEMFraction(Double_t m)     { fMaxNeutralEMFraction = m; }
    void SetMaxMuonFraction(Double_t m)          { fMaxMuonFraction = m; }
    void SetMinMuonFraction(Double_t m)          { fMinMuonFraction = m; }
    void SetApplyBetaCut(Bool_t b)               { fApplyBetaCut = b; }
    void SetPFId(UInt_t w)                       { fPFId = w; }
    void SetMVATrainingSet(UInt_t s)             { fMVATrainingSet = s; }
    void SetMVACutWP(UInt_t w)                   { fMVACutWP = w; }
    void SetMVAWeightsFile(char const*, UInt_t idx = 0);
    void SetMVACutsFile(char const* n)           { fMVACutsFile = n; }
    void SetUseClassicBetaForMVA(Bool_t b)       { fUseClassicBetaForMVA = b; }
    void SetMinNPFCandidates(UInt_t n)           { fMinNPFCandidates = n; }
    void SetMinNChargedPFCandidates(UInt_t n)    { fMinNChargedPFCandidates = n; }

    void SetJetIDMVA(JetIDMVA* mva) { fJetIDMVA = mva; }

    void SetChargedHadronFractionRange(Double_t min, Double_t max)
    { fMinChargedHadronFraction = min; fMaxChargedHadronFraction = max; }
    void SetNeutralHadronFraction(Double_t min, Double_t max)
    { fMinNeutralHadronFraction = min; fMaxNeutralHadronFraction = max; }
    void SetChargedEMFraction(Double_t min, Double_t max)
    { fMinChargedEMFraction = min; fMaxChargedEMFraction = max; }
    void SetNeutralEMFraction(Double_t min, Double_t max)
    { fMinNeutralEMFraction = min; fMaxNeutralEMFraction = max; }
    void SetMuonFraction(Double_t min, Double_t max)
    { fMinMuonFraction = min; fMaxMuonFraction = max; }

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
      cMuonFrac,
      cNPFCandidates,
      cNChargedPFCandidates,
      cPFId,
      cBeta,
      cMVA,
      nCuts
    };

    Bool_t IsGood(mithep::Jet const&) override;
    void IdBegin() override;
    void IdTerminate() override;

    Double_t fJetEEMFractionMinCut = 0.01;  //jet Eem fraction min cut for calo jets
    Double_t fMinChargedHadronFraction = 0.;
    Double_t fMaxChargedHadronFraction = 1.;
    Double_t fMinNeutralHadronFraction = 0.;
    Double_t fMaxNeutralHadronFraction = 1.;
    Double_t fMinChargedEMFraction = 0.;
    Double_t fMaxChargedEMFraction = 1.;
    Double_t fMinNeutralEMFraction = 0.;
    Double_t fMaxNeutralEMFraction = 1.;
    Double_t fMinMuonFraction = 0.;
    Double_t fMaxMuonFraction = 1.;
    UInt_t   fMinNPFCandidates = 0;
    UInt_t   fMinNChargedPFCandidates = 0;
    Bool_t   fApplyBetaCut = kFALSE;          //=true then apply beta cut
    UInt_t   fPFId = JetTools::nPFIdWorkingPoints; // PF Id working point
    UInt_t   fMVATrainingSet = JetIDMVA::nMVATypes; //JetIDMVA::MVAType
    UInt_t   fMVACutWP = JetIDMVA::kLoose; //JetIDMVA::CutType
    std::vector<TString>  fMVAWeightsFile{""};
    TString  fMVACutsFile = "";
    Bool_t   fUseClassicBetaForMVA = kFALSE; //set to true to replicate CMSSW PU jet ID on MiniAOD

    Bool_t    fOwnJetIDMVA = kFALSE;
    JetIDMVA* fJetIDMVA = 0;

    ClassDef(JetIdMod, 1) // Jet identification module
  };
}

#endif
