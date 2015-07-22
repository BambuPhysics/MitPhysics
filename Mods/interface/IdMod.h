//--------------------------------------------------------------------------------------------------
// IdMod
//
// Base class for object Id modules
//
// Authors: Y.Iiyama
//--------------------------------------------------------------------------------------------------


#ifndef MITPHYSICS_MODS_IDMOD_H
#define MITPHYSICS_MODS_IDMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataCont/interface/NamedFastArrayBasic.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/TriggerObjectFwd.h"
#include "MitAna/DataTree/interface/DecayParticleFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/BeamSpotFwd.h"
#include "MitAna/DataTree/interface/PFCandidateFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityFwd.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitPhysics/Init/interface/ModNames.h"

#include "TString.h"
#include "TH1D.h"

#include <limits>

namespace mithep {

  template<class O>
  class IdMod : public BaseMod {
  public:
    IdMod(char const* name, char const* title);
    ~IdMod();

    char const* GetInputName() const { return fInputName; }
    char const* GetOutputName() const { return fOutputName; }

    char const* GetTriggerObjectsName() const       { return fAuxInputNames[kTrigObjects]; }
    char const* GetConversionsName() const          { return fAuxInputNames[kConversions]; }
    char const* GetElectronsName() const            { return fAuxInputNames[kElectrons]; }
    char const* GetVerticesName() const             { return fAuxInputNames[kVertices]; }
    char const* GetBeamSpotName() const             { return fAuxInputNames[kBeamSpot]; }
    char const* GetPFCandidatesName() const         { return fAuxInputNames[kPFCandidates]; }
    char const* GetPFNoPileupCandidatesName() const { return fAuxInputNames[kPFNoPileupCandidates]; }
    char const* GetPileupEnergyDensityName() const  { return fAuxInputNames[kPileupEnergyDensity]; }
    char const* GetNonIsolatedMuonsName() const     { return fAuxInputNames[kNonIsolatedMuons]; }
    char const* GetNonIsolatedElectronsName() const { return fAuxInputNames[kNonIsolatedElectrons]; }

    Bool_t GetIsFilterMode() const { return fIsFilterMode; }
    UInt_t GetIdType() const { return fIdType; }
    UInt_t GetIsoType() const { return fIsoType; }
    Double_t GetPtMin() const { return fPtMin; }
    Double_t GetEtaMax() const { return fEtaMax; }
    UInt_t GetMinOutput() const { return fMinOutput; }

    void SetInputName(char const* n) { fInputName = n; }
    void SetOutputName(char const* n) { fOutputName = n; }

    void SetTriggerObjectsName(const char* n)       { fAuxInputNames[kTrigObjects] = n; }
    void SetConversionsName(const char* n)          { fAuxInputNames[kConversions] = n; }
    void SetElectronsName(const char* n)            { fAuxInputNames[kElectrons] = n; }
    void SetVertexName(const char* n)               { fAuxInputNames[kVertices] = n; }
    void SetBeamSpotName(const char* n)             { fAuxInputNames[kBeamSpot] = n; }
    void SetPFCandidatesName(const char* n)         { fAuxInputNames[kPFCandidates] = n; } 
    void SetPFNoPileupCandidatesName(const char* n) { fAuxInputNames[kPFNoPileupCandidates] = n; } 
    void SetPFPileupCandidatesName(const char* n)   { fAuxInputNames[kPFPileupCandidates] = n; } 
    void SetPileupEnergyDensityName(const char* n)  { fAuxInputNames[kPileupEnergyDensity] = n; }
    void SetNonIsolatedMuonsName(const char* n)     { fAuxInputNames[kNonIsolatedMuons] = n; }
    void SetNonIsolatedElectronsName(const char* n) { fAuxInputNames[kNonIsolatedElectrons] = n; }

    void SetIsFilterMode(Bool_t b) { fIsFilterMode = b; }
    void SetIdType(UInt_t t) { fIdType = t; }
    void SetIsoType(UInt_t t) { fIsoType = t; }
    void SetPtMin(Double_t m) { fPtMin = m; }
    void SetEtaMax(Double_t m) { fEtaMax = m; }
    void SetMinOutput(UInt_t m) { fMinOutput = m; }

  protected:
    void SlaveBegin() override;
    void SlaveTerminate() override;
    void Process() override;
    virtual Bool_t IsGood(O const&) = 0;
    virtual void IdBegin() {}
    virtual void IdTerminate() {}

    enum AuxInput {
      kTrigObjects,
      kConversions,
      kElectrons,
      kVertices,
      kBeamSpot,
      kPFCandidates,
      kPFNoPileupCandidates,
      kPFPileupCandidates,
      kPileupEnergyDensity,
      kNonIsolatedMuons,
      kNonIsolatedElectrons,
      nAuxInputs
    };

    template<class T> T const* GetAuxInput(AuxInput);

    mithep::TriggerObjectCol const* GetTrigObjects()
    { return GetAuxInput<mithep::TriggerObjectCol>(kTrigObjects); }
    mithep::DecayParticleCol const* GetConversions()
    { return GetAuxInput<mithep::DecayParticleCol>(kConversions); }
    mithep::ElectronCol const* GetElectrons()
    { return GetAuxInput<mithep::ElectronCol>(kElectrons); }
    mithep::VertexCol const* GetVertices()
    { return GetAuxInput<mithep::VertexCol>(kVertices); }
    mithep::BeamSpotCol const* GetBeamSpot()
    { return GetAuxInput<mithep::BeamSpotCol>(kBeamSpot); }
    mithep::PFCandidateCol const* GetPFCandidates()
    { return GetAuxInput<mithep::PFCandidateCol>(kPFCandidates); }
    mithep::PFCandidateCol const* GetPFNoPileupCandidates()
    { return GetAuxInput<mithep::PFCandidateCol>(kPFNoPileupCandidates); }
    mithep::PFCandidateCol const* GetPFPileupCandidates()
    { return GetAuxInput<mithep::PFCandidateCol>(kPFPileupCandidates); }
    mithep::PileupEnergyDensityCol const* GetPileupEnergyDensity()
    { return GetAuxInput<mithep::PileupEnergyDensityCol>(kPileupEnergyDensity); }
    mithep::MuonCol const* GetNonIsolatedMuons()
    { return GetAuxInput<mithep::MuonCol>(kNonIsolatedMuons); }
    mithep::ElectronCol const* GetNonIsolatedElectrons()
    { return GetAuxInput<mithep::ElectronCol>(kNonIsolatedElectrons); }
    
    // fGoodObjects: Collection of Id'ed objects (filter mode)
    ObjArray<O> fGoodObjects;
    // fFlags: Good / bad flags published in flag mode
    NamedFastArrayBasic<Bool_t> fFlags;

    TH1D* fCutFlow = 0;

    TString        fInputName = ""; // input collection of objects to be Id'ed
    TString        fOutputName = ""; // input collection of objects to be Id'ed
    TString        fAuxInputNames[nAuxInputs] = {};
    TObject const* fAuxInputs[nAuxInputs] = {};

    Bool_t   fIsFilterMode = kTRUE;
    UInt_t   fMinOutput = 0;
    UInt_t   fIdType = 0xffffffff;
    UInt_t   fIsoType = 0xffffffff;
    Double_t fPtMin = 0.;
    Double_t fEtaMax = std::numeric_limits<double>::max();

    ClassDef(IdMod, 0)
  };

  template<class O>
  mithep::IdMod<O>::IdMod(char const* name, char const* title) :
    BaseMod(name, title),
    fGoodObjects(32, TString(name) + "Output"),
    fFlags(32, TString(name) + "Flags")
  {
    fAuxInputNames[kConversions] = Names::gkMvfConversionBrn;
    fAuxInputNames[kElectrons] = Names::gkElectronBrn;
    fAuxInputNames[kVertices] = ModNames::gkGoodVertexesName;
    fAuxInputNames[kBeamSpot] = Names::gkBeamSpotBrn;
    fAuxInputNames[kPFCandidates] = Names::gkPFCandidatesBrn;
    fAuxInputNames[kPileupEnergyDensity] = Names::gkPileupEnergyDensityBrn;
  }
  
  template<class O>
  mithep::IdMod<O>::~IdMod()
  {
  }

  template<class O>
  void
  mithep::IdMod<O>::SlaveBegin()
  {
    if (fIsFilterMode) {
      fGoodObjects.SetName(fOutputName);
      if (!PublishObj(&fGoodObjects))
        SendError(kAbortAnalysis, "SlaveBegin", "Cannot publish output");
    }
    else {
      fFlags.SetName(fOutputName);
      if (!PublishObj(&fFlags))
        SendError(kAbortAnalysis, "SlaveBegin", "Cannot publish output");
    }
  
    AddTH1(fCutFlow, TString(GetName()) + "CutFlow", "Identification cut flow", 1, 0., 1.);
  
    IdBegin();
  }
  
  template<class O>
  void
  mithep::IdMod<O>::SlaveTerminate()
  {
    RetractObj(GetOutputName());
  
    IdTerminate();
  }
  
  template<class O>
  void
  mithep::IdMod<O>::Process()
  {
    for (unsigned iA = 0; iA != nAuxInputs; ++iA)
      fAuxInputs[iA] = 0;

    // Process entries of the tree.
    auto* input = GetObject<mithep::Collection<O> >(fInputName, true);
    if (!input)
      SendError(kAbortAnalysis, "Process", fInputName + " not found");

    if (fIsFilterMode) {
      fGoodObjects.Reset();
    }
    else {
      fFlags.Resize(input->GetEntries());
      for (UInt_t iO = 0; iO != input->GetEntries(); ++iO)
        fFlags.At(iO) = false;
    }

    UInt_t nGoodObjects = 0;
    for (UInt_t iO = 0; iO != input->GetEntries(); ++iO) {
      O const& obj(*input->At(iO));
      if (IsGood(obj)) {
        if (fIsFilterMode)
          fGoodObjects.Add(&obj);
        else
          fFlags.At(iO) = true;

        obj.Mark();
        ++nGoodObjects;
      }
    }
  
    if (nGoodObjects < fMinOutput) {
      SkipEvent();
      return;
    }

    if (fIsFilterMode) {
      // sort according to pt
      fGoodObjects.Sort();
    }
  }

  template<class O>
  template<class T>
  T const*
  mithep::IdMod<O>::GetAuxInput(AuxInput inputCol)
  {
    if (fAuxInputs[inputCol])
      return static_cast<T const*>(fAuxInputs[inputCol]);

    if (inputCol == kTrigObjects) {
      fAuxInputs[kTrigObjects] = GetHLTObjects(fAuxInputNames[kTrigObjects]);
    }
    else {
      if (fAuxInputNames[inputCol].Length() == 0)
        SendError(kAbortAnalysis, "GetAuxInput", "Auxiliary input name is empty");
  
      fAuxInputs[inputCol] = GetObject<T>(fAuxInputNames[inputCol], true);
    }

    if (!fAuxInputs[inputCol])
      SendError(kAbortAnalysis, "GetAuxInput", "Could not retrieve auxiliary input " + fAuxInputNames[inputCol]);

    return static_cast<T const*>(fAuxInputs[inputCol]);
  }

}

#endif
