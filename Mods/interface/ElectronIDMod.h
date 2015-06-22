//--------------------------------------------------------------------------------------------------
// ElectronIDMod
//
// This module applies electron identification criteria and exports a pointer to a collection
// of "good electrons" according to the specified identification scheme.
//
// Authors: S.Xie, C.Loizides, Y.Iiyama
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_ELECTRONIDMOD_H
#define MITPHYSICS_MODS_ELECTRONIDMOD_H

#include "MitPhysics/Mods/interface/IDMod.h"

#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"

namespace mithep {
  class ElectronIDMod : public IDMod {
  public:
    ElectronIDMod(const char* name="ElectronIDMod",
                  const char* title="Electron identification module");

    const char* GetConversionBranchName() const     { return fAuxInputNames[kConversions]; }
    const char* GetNonIsolatedMuonsName() const     { return fAuxInputNames[kNonIsolatedMuons]; }
    const char* GetNonIsolatedElectronsName() const { return fAuxInputNames[kNonIsolatedElectrons]; }
    const char* GetVerticesName() const             { return fAuxInputNames[kVertices]; }
    const char* GetTriggerObjectsName() const       { return fAuxInputNames[kTrigObjects]; }
    const char* GetPileupEnergyDensityName() const  { return fAuxInputNames[kPileupEnergyDensity]; }

    Bool_t   GetApplyConversionFilterType1() const  { return fApplyConvFilterType1; }
    Bool_t   GetApplyConversionFilterType2() const  { return fApplyConvFilterType2; }
    Bool_t   GetApplyNExpectedHitsInnerCut() const  { return fApplyNExpectedHitsInnerCut; }
    Bool_t   GetInvertNExpectedHitsInnerCut() const { return fInvertNExpectedHitsInnerCut; }
    Bool_t   GetApplySpikeRemoval() const           { return fApplySpikeRemoval; }
    Bool_t   GetApplyD0Cut() const                  { return fApplyD0Cut; }
    Bool_t   GetApplyDZCut() const                  { return fApplyDZCut; }
    Bool_t   GetChargeFilter() const                { return fChargeFilter; }
    Bool_t   GetApplyTriggerMatching() const        { return fApplyTriggerMatching; }
    Bool_t   GetApplyEcalSeeded() const             { return fApplyEcalSeeded; }
    Bool_t   GetApplyEcalFiducial() const           { return fApplyEcalFiducial; }
    UInt_t   GetRhoAlgo() const                     { return fRhoAlgo; }
    Int_t    GetWhichVertex() const                 { return fWhichVertex; }
    Double_t GetEtMin() const                       { return fElectronEtMin; }
    Double_t GetIDLikelihoodCut() const             { return fIDLikelihoodCut; }

    ElectronLikelihood* GetLH() const { return fLH; }

    void SetConversionBranchName(const char* n)    { fAuxInputNames[kConversions] = n; }
    void SetNonIsolatedMuonsName(const char* n)    { fAuxInputNames[kNonIsolatedMuons] = n; }
    void SetNonIsolatedElectronsName(const char* n){ fAuxInputNames[kNonIsolatedElectrons] = n; }
    void SetVertexName(const char* n)              { fAuxInputNames[kVertices] = n; }
    void SetTriggerObjectsName(const char* n)      { fAuxInputNames[kTrigObjects] = n; }
    void SetPileupEnergyDensityName(const char* n) { fAuxInputNames[kPileupEnergyDensity] = n; }

    void SetApplyConversionFilterType1(Bool_t b)   { fApplyConvFilterType1 = b; }
    void SetApplyConversionFilterType2(Bool_t b)   { fApplyConvFilterType2 = b; }
    void SetApplyNExpectedHitsInnerCut(Bool_t b)   { fApplyNExpectedHitsInnerCut = b; }
    void SetInvertNExpectedHitsInnerCut(Bool_t b)  { fInvertNExpectedHitsInnerCut = b; }
    void SetApplySpikeRemoval(Bool_t b)            { fApplySpikeRemoval = b; }
    void SetApplyD0Cut(Bool_t b)                   { fApplyD0Cut = b; }
    void SetApplyDZCut(Bool_t b)                   { fApplyDZCut = b; }
    void SetChargeFilter(Bool_t b)                 { fChargeFilter = b; }
    void SetApplyTriggerMatching(Bool_t b)         { fApplyTriggerMatching = b; }
    void SetApplyEcalSeeded(Bool_t b)              { fApplyEcalSeeded = b; }
    void SetApplyEcalFiducial(Bool_t b)            { fApplyEcalFiducial = b; }
    void SetRhoAlgo(UInt_t algo)                   { fRhoAlgo = algo; }
    void SetWhichVertex(Int_t d)                   { fWhichVertex = d; }
    void SetEtMin(Double_t et)                     { fElectronEtMin = et; }
    void SetIDLikelihoodCut(Double_t cut)          { fIDLikelihoodCut = cut; }

    void SetLH(ElectronLikelihood* l)              { fLH = l; }

  protected:
    enum AuxInput {
      kTrigObjects,
      kConversions,
      kVertices,
      kBeamSpot,
      kPFCandidates,
      kPileupEnergyDensity,
      kNonIsolatedMuons,
      kNonIsolatedElectrons,
      nAuxInputs
    };

    enum CutFlow {
      cAll,
      cIsEcalDriven,
      cPt,
      cEta,
      cEt,
      cFiducial,
      cSpikeRemoval,
      cChargeFilter,
      cTriggerMatching,
      cConvFilterType1,
      cConvFilterType2,
      cNExpectedHits,
      cD0,
      cDZ,
      cId,
      cIsolation,
      nCuts
    };

    Bool_t   PassLikelihoodID(Electron const&);
    Bool_t   PassIDCut(Electron const&, TObject const**);
    Bool_t   PassIsolationCut(Electron const&, TObject const**);
    template<class T> void GetAuxInput(AuxInput, TObject const**);

    void Process() override;
    void IdBegin() override;
    
    TString  fAuxInputNames[nAuxInputs];

    Bool_t   fApplyConvFilterType1 = kTRUE;   //whether remove conversions using fit method
    Bool_t   fApplyConvFilterType2 = kFALSE;   //whether remove conversions using DCotTheta method
    Bool_t   fApplyNExpectedHitsInnerCut = kTRUE;
    Bool_t   fInvertNExpectedHitsInnerCut = kFALSE; //whether to invert NExpectedHitsInner cut
    Bool_t   fApplySpikeRemoval = kTRUE;      //whether spike removal
    Bool_t   fApplyD0Cut = kTRUE;             //whether apply d0 cut
    Bool_t   fApplyDZCut = kTRUE;             //whether apply dz cut
    Bool_t   fChargeFilter = kTRUE;           //whether apply GSF and CFT equal requirement
    Bool_t   fApplyTriggerMatching = kFALSE;   //match to hlt electron (default=0)
    Bool_t   fApplyEcalSeeded = kFALSE;        //require ecal seeded flag
    Bool_t   fApplyEcalFiducial = kFALSE;      //apply ecal fiducial cuts on supercluster eta
    UInt_t   fRhoAlgo = mithep::PileupEnergyDensity::kFixedGridFastjetAll;
    Int_t    fWhichVertex = -1;            //vertex to use (-2: beamspot, -1: closest in Z)
    Double_t fElectronEtMin = 0.;          //min pt cut
    Double_t fIDLikelihoodCut = -999.;        //cut value for ID likelihood

    // pointers not owned
    ElectronLikelihood*           fLH = 0;

    ClassDef(ElectronIDMod, 0) // Electron identification module
  };
}

#endif
