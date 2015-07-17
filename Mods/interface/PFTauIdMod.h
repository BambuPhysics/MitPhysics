//--------------------------------------------------------------------------------------------------
// PFTauIdMod
//
// This module applies tau identification criteria and exports a pointer to a collection
// of "good Taus" according to the specified identification scheme.
//
// Authors: Y.Iiyama
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PFTauIDMod_H
#define MITPHYSICS_MODS_PFTauIDMod_H

#include "MitPhysics/Mods/interface/IdMod.h" 
#include "MitAna/DataTree/interface/PFTau.h"

namespace mithep {

  class PFTauIdMod : public IdMod<PFTau> {
  public:
    PFTauIdMod(char const* name = "PFTauIdMod", 
               char const* title = "Tau identification module");

    Long64_t GetDiscriminatorMask() const;
    Bool_t   GetIsHPSSel() const                     { return fIsHPSSel; }
    // cuts for non-HPS taus
    Double_t GetPtLeadChargedHadronPFCandMin() const { return fPtLeadChargedHadronPFCandMin; }
    UInt_t   GetIsoChargedHadronPtSumMax() const     { return fIsoChargedHadronPtSumMax; }
    UInt_t   GetIsoGammaEtSumMax() const             { return fIsoGammaEtSumMax; }
    Double_t GetSignalMassMin() const                { return fSignalMassMin; }
    Double_t GetSignalMassMax() const                { return fSignalMassMax; }

    void AddDiscriminator(UInt_t d, Bool_t invert = kFALSE)
    { AddCutDiscriminator(d, 0.5, !invert); }
    void AddCutDiscriminator(UInt_t d, Double_t c, Bool_t passAbove = kTRUE)
    { fDiscriminators[d] = CutConfig(c, passAbove); }
    void RemoveDiscriminator(UInt_t d)               { fDiscriminators.erase(d); }
    void SetIsHPSSel(Bool_t b)                       { fIsHPSSel = b; }
    void SetPtLeadChargedHadronPFCandMin(Double_t x) { fPtLeadChargedHadronPFCandMin = x; }
    void SetIsoChargedHadronPtSumMax(Double_t x)     { fIsoChargedHadronPtSumMax = x; }
    void SetIsoGammaEtSumMax(Double_t x)             { fIsoGammaEtSumMax = x; }
    void SetSignalMassMin(Double_t x)                { fSignalMassMin = x; }
    void SetSignalMassMax(Double_t x)                { fSignalMassMax = x; }

  protected:
    enum CutFlow {
      cAll,
      cPt,
      cEta,
      cDisc,
      nHPSCuts = cDisc,
      cCharge = cEta + 1,
      cNPF,
      cLeadPF,
      cLeadPt,
      cChargedPtSum,
      cGammaEtSum,
      cNTrk,
      cSystemPt,
      cSystemMass,
      nNonHPSCuts
    };

    Bool_t IsGood(PFTau const&) override;
    void   IdBegin() override;

    typedef std::pair<Double_t, Bool_t> CutConfig;

    std::map<UInt_t, CutConfig> fDiscriminators{};
    Bool_t   fIsHPSSel = kTRUE;               //apply HPS Tau selection
    Double_t fPtLeadChargedHadronPFCandMin = 5.; //min LeadChargedHadronPFCand pt cut
    Double_t fIsoChargedHadronPtSumMax = 2.;       //maximum of Pt iso tracks
    Double_t fIsoGammaEtSumMax = 3.;               //maximum of Pt iso neutrals
    Double_t fSignalMassMin = 0.13;                  //minimum of tau mass
    Double_t fSignalMassMax = 2.;                  //maximum of tau mass
    Bool_t   fIsLooseId = kTRUE;              //apply Loose Tau selection
    
    ClassDef(PFTauIdMod, 0) // Tau identification module
  };
}

#endif
