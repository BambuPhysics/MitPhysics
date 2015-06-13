//--------------------------------------------------------------------------------------------------
// ElectronIDMod
//
// This module applies electron identification criteria and exports a pointer to a collection
// of "good electrons" according to the specified identification scheme.
//
// See http://indico.cern.ch/contributionDisplay.py?contribId=1&confId=42251
//
// Authors: S.Xie, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_ELECTRONIDMOD_H
#define MITPHYSICS_MODS_ELECTRONIDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/DecayParticleFwd.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"

// for Rho definitons
#include "MitPhysics/Utils/interface/RhoUtilities.h"

namespace mithep {
  class ElectronIDMod : public BaseMod {
  public:
    ElectronIDMod(const char* name="ElectronIDMod",
                  const char* title="Electron identification module");

    const char* GetGoodElectronsName() const       { return fGoodElectronsName; }
    const char* GetOutputName() const              { return GetGoodElectronsName(); }
    const char* GetGoodName() const                { return GetGoodElectronsName(); }
    const char* GetInputName() const               { return fElectronBranchName; }
    const char* GetConversionBranchName() const    { return fConversionBranchName; }
    const char* GetOldMuonsName() const            { return fNonIsolatedMuonsName; }
    const char* GetOldElectronsName() const        { return fNonIsolatedElectronsName ; }
    const char* GetVertexName() const              { return fVertexName; }
    const char* GetPVName() const                  { return GetVertexName(); }
    const char* GetTriggerObjectsName() const      { return fTrigObjectsName; }
    const char* GetPFNoPileUpName() const          { return fPFNoPileUpName; }
    const char* GetPileupEnergyDensityName() const { return fPileupEnergyDensityName; }

    Bool_t   GetApplyConversionFilterType1() const  { return fApplyConvFilterType1; }
    Bool_t   GetApplyConversionFilterType2() const  { return fApplyConvFilterType2; }
    Bool_t   GetApplyNExpectedHitsInnerCut() const  { return fApplyNExpectedHitsInnerCut; }
    Bool_t   GetInvertNExpectedHitsInnerCut() const { return fInvertNExpectedHitsInnerCut; }
    Bool_t   GetCombinedIdCut() const               { return fCombinedIdCut; }
    Bool_t   GetApplySpikeRemoval() const           { return fApplySpikeRemoval; }
    Bool_t   GetApplyD0Cut() const                  { return fApplyD0Cut; }
    Bool_t   GetApplyDZCut() const                  { return fApplyDZCut; }
    Bool_t   GetChargeFilter() const                { return fChargeFilter; }
    Bool_t   GetApplyTriggerMatching() const        { return fApplyTriggerMatching; }
    Bool_t   GetApplyEcalSeeded() const             { return fApplyEcalSeeded; }
    Bool_t   GetApplyEcalFiducial() const           { return fApplyEcalFiducial; }
    Bool_t   GetElectronsFromBranch() const         { return fElectronsFromBranch; }
    UInt_t   GetIDType() const                      { return fElIdType; }
    UInt_t   GetIsoType() const                     { return fElIsoType; }
    UInt_t   GetRhoAlgo() const                     { return fRhoAlgo; }
    Int_t    GetWhichVertex() const                 { return fWhichVertex; }
    Double_t GetPtMin() const                       { return fElectronPtMin; }
    Double_t GetEtMin() const                       { return fElectronEtMin; }
    Double_t GetEtaMax() const                      { return fElectronEtaMax; }
    Double_t GetIDLikelihoodCut() const             { return fIDLikelihoodCut; }
    Double_t GetIntRadius() const                   { return fIntRadius; }

    ElectronLikelihood* GetLH() const { return fLH; }
    ElectronIDMVA* GetMVA() const { return fElectronIDMVA; }

    UInt_t GetRhoType() const { return -1; } /* DEPRECATED */

    void SetGoodElectronsName(const char* n)       { fGoodElectronsName = n; }
    void SetOutputName(const char* n)              { SetGoodElectronsName(n); }
    void SetGoodName(const char* n)                { SetGoodElectronsName(n); }
    void SetInputName(const char* n)               { fElectronBranchName = n; }
    void SetConversionBranchName(const char* n)    { fConversionBranchName = n; }
    void SetOldMuonsName(const char* n)            { fNonIsolatedMuonsName = n; }
    void SetOldElectronsName(const char* n)        { fNonIsolatedElectronsName = n; }
    void SetVertexName(const char* n)              { fVertexName = n; }
    void SetPVName(const char* n)                  { SetVertexName(n); }
    void SetTriggerObjectsName(const char* n)      { fTrigObjectsName = n; }
    void SetPFNoPileUpName(const char* n)          { fPFNoPileUpName = n; }
    void SetPileupEnergyDensityName(const char* n) { fPileupEnergyDensityName = n; }

    void SetApplyConversionFilterType1(Bool_t b)   { fApplyConvFilterType1 = b; }
    void SetApplyConversionFilterType2(Bool_t b)   { fApplyConvFilterType2 = b; }
    void SetApplyNExpectedHitsInnerCut(Bool_t b)   { fApplyNExpectedHitsInnerCut = b; }
    void SetInvertNExpectedHitsInnerCut(Bool_t b)  { fInvertNExpectedHitsInnerCut = b; }
    void SetCombinedIdCut(Bool_t b)                { fCombinedIdCut = b; }
    void SetApplySpikeRemoval(Bool_t b)            { fApplySpikeRemoval = b; }
    void SetApplyD0Cut(Bool_t b);
    void SetApplyDZCut(Bool_t b);
    void SetChargeFilter(Bool_t b)                 { fChargeFilter = b; }
    void SetApplyTriggerMatching(Bool_t b)         { fApplyTriggerMatching = b; }
    void SetApplyEcalSeeded(Bool_t b)              { fApplyEcalSeeded = b; }
    void SetApplyEcalFiducial(Bool_t b)            { fApplyEcalFiducial = b; }
    void SetElectronsFromBranch(Bool_t b)          { fElectronsFromBranch = b; }
    void SetIDType(UInt_t t)                       { fElIdType = static_cast<ElectronTools::EElIdType>(t); }
    void SetIsoType(UInt_t t)                      { fElIsoType = static_cast<ElectronTools::EElIsoType>(t); }
    void SetRhoAlgo(UInt_t algo)                   { fRhoAlgo = algo; }
    void SetWhichVertex(Int_t d)                   { fWhichVertex = d; }
    void SetPtMin(Double_t pt)                     { fElectronPtMin = pt; }
    void SetEtMin(Double_t et)                     { fElectronEtMin = et; }
    void SetEtaMax(Double_t eta)                   { fElectronEtaMax = eta; }
    void SetIDLikelihoodCut(Double_t cut)          { fIDLikelihoodCut = cut; }
    void SetIntRadius(Double_t dr)                 { fIntRadius = dr; }

    void SetLH(ElectronLikelihood* l)              { fLH = l; }
    void SetMVA(ElectronIDMVA* mva)                { fElectronIDMVA = mva; }

    /* DEPRECATED METHODS */
    void SetRhoType(RhoUtilities::RhoType);
    void SetNExpectedHitsInnerCut(Int_t);

  protected:
    void SetCutValues();
    Bool_t   PassLikelihoodID(const Electron *ele) const;
    Bool_t   PassMVAID(const Electron *el) const;
    Double_t EvaluateMVAID(const Electron *el) const;
    Bool_t   PassIDCut(const Electron *el) const;
    Bool_t   PassIsolationCut(const Electron *el) const;

    void Process() override;
    void SlaveBegin() override;
    void Terminate() override;

    TString fGoodElectronsName;      //name of exported "good electrons" col
    TString fElectronBranchName;     //name of electron collection (input)
    TString fConversionBranchName;   //name of electron collection (input)
    TString fVertexName;             //name of vertex collection
    TString fBeamSpotName;           //name of beamspot collection
    TString fTrackName;              //name of track collection
    TString fPFCandidatesName;       //name of pfcandidates collection
    TString fPileupEnergyDensityName;
    TString fPFNoPileUpName;         //name of pfcandidates collection
    TString fTrigObjectsName;        //name of trigger object collection
    TString fNonIsolatedMuonsName;    //name of imported "old muon" collection
    TString fNonIsolatedElectronsName;//name of imported "old electron" collection

    const ElectronCol*            fElectrons;              //!electron collection
    const DecayParticleCol*       fConversions;            //!conversion collection
    const VertexCol*              fVertices;               //!vertices branches
    const BeamSpotCol*            fBeamSpot;               //!beamspot branch
    const TrackCol*               fTracks;                 //!Track branch
    const PFCandidateCol*         fPFCandidates;           //!pfcandidate branch
    const PileupEnergyDensityCol* fPileupEnergyDensity;
    const PFCandidateCol*         fPFNoPileUpCands;        //!pfcandidate branch with PFNoPU
    const MuonCol*                fNonIsolatedMuons;       //!pointer to old muon collection
    const ElectronCol*            fNonIsolatedElectrons;   //!pointer to old electron collection

    ElectronTools::EElIdType  fElIdType;              //!identification scheme
    ElectronTools::EElIsoType fElIsoType;              //!isolation scheme
    Bool_t   fApplyConvFilterType1;   //whether remove conversions using fit method
    Bool_t   fApplyConvFilterType2;   //whether remove conversions using DCotTheta method
    Bool_t   fApplyNExpectedHitsInnerCut;
    Bool_t   fInvertNExpectedHitsInnerCut; //whether to invert NExpectedHitsInner cut
    Bool_t   fCombinedIdCut;          //whether to use full combined id
    Bool_t   fApplySpikeRemoval;      //whether spike removal
    Bool_t   fApplyD0Cut;             //whether apply d0 cut
    Bool_t   fApplyDZCut;             //whether apply dz cut
    Bool_t   fChargeFilter;           //whether apply GSF and CFT equal requirement
    Bool_t   fApplyTriggerMatching;   //match to hlt electron (default=0)
    Bool_t   fApplyEcalSeeded;        //require ecal seeded flag
    Bool_t   fApplyEcalFiducial;      //apply ecal fiducial cuts on supercluster eta
    Bool_t   fElectronsFromBranch;    //where to get input electrons
    UInt_t   fRhoAlgo;
    Int_t    fWhichVertex;            //vertex to use (-2: beamspot, -1: closest in Z)
    Double_t fElectronPtMin;          //min pt cut
    Double_t fElectronEtMin;          //min pt cut
    Double_t fElectronEtaMax;         //max eta cut
    Double_t fIDLikelihoodCut;        //cut value for ID likelihood
    Double_t fIntRadius;              //!min IntRadius cut in pf isolation

    // pointers not owned
    ElectronLikelihood*           fLH;                     //LH
    ElectronIDMVA*                fElectronIDMVA;

    ClassDef(ElectronIDMod, 1) // Electron identification module
  };
}
#endif
