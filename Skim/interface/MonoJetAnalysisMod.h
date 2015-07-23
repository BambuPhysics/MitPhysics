//--------------------------------------------------------------------------------------------------
// MonoJetAnalysisMod
//
// An analysis module for selecting signal and calibration regions for the MonoJet analysis
// and produces some basic distributions.
//
// Authors: LDM, TJW, YI
//--------------------------------------------------------------------------------------------------
#ifndef MITPHYSICS_SKIM_MONOJETANALYSISMOD_H
#define MITPHYSICS_SKIM_MONOJETANALYSISMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataCont/interface/Types.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"

#include "TString.h"
#include "TH1D.h"

namespace mithep {

  class MonoJetAnalysisMod : public BaseMod {
  public:
    MonoJetAnalysisMod(const char *name  = "MonoJetAnalysisMod", 
                       const char *title = "MoneJet Slection Module");
    ~MonoJetAnalysisMod() {}

    char const* GetCategoryFlagsName() const { return fCategoryFlags.GetName(); }
    
    // set input names
    void SetMetName(const char* n)           { fMetName= n; }
    void SetJetsName(const char* n)          { fJetsName = n; } 
    void SetVetoElectronsName(const char* n) { fVetoElectronsName = n; }
    void SetElectronMaskName(const char* n)  { fElectronMaskName = n; }
    void SetVetoMuonsName(const char* n)     { fVetoMuonsName = n; }
    void SetMuonMaskName(const char* n)      { fMuonMaskName = n; }
    void SetVetoTausName(const char* n)      { fVetoTausName = n; }
    void SetVetoPhotonsName(const char* n)   { fVetoPhotonsName = n; }
    void SetPhotonMaskName(const char* n)    { fPhotonMaskName = n; }

    // set output name
    void SetCategoryFlagsName(const char* n) { fCategoryFlags.SetName(n); }

    // Setting cut values
    void SetCategoryActive(UInt_t c, Bool_t a = kTRUE) { fCategoryActive[c] = a; }
    void AddTriggerName(UInt_t c, const char* n)       { fTriggerNames[c].push_back(n); }
    void SetMaxNumJets(UInt_t c, Int_t n)              { fMaxNumJets[c] = n; }
    void SetMinLeadJetPt(UInt_t c, Double_t x)         { fMinLeadJetPt[c] = x; }
    void SetMaxJetEta(UInt_t c, Double_t x)            { fMaxJetEta[c] = x; }
    void SetMinMetPt(UInt_t c, Double_t x)             { fMinMetPt[c] = x; }
    void SetMinChargedHadronFrac(UInt_t c, Double_t x) { fMinChargedHadronFrac[c] = x; }
    void SetMaxNeutralHadronFrac(UInt_t c, Double_t x) { fMaxNeutralHadronFrac[c] = x; }
    void SetMaxNeutralEmFrac(UInt_t c, Double_t x)     { fMaxNeutralEmFrac[c] = x; }

    enum MonoJetCategory {
      kSignal,
      kDielectron,
      kDimuon,
      kSingleElectron,
      kSingleMuon,
      kPhoton,
      nMonoJetCategories
    };

    static TString fgMonoJetCategories[nMonoJetCategories];

  protected:
    // Standard module methods
    void SlaveBegin() override;
    void BeginRun() override;
    void Process() override;
    void SlaveTerminate() override;

    // names of the input collections
    TString fMetName{"PFType1CorrectedMet"};
    TString fJetsName{ModNames::gkCleanJetsName};
    TString fVetoElectronsName{"VetoElectrons"};
    TString fElectronMaskName{"GoodElectrons"};
    TString fVetoMuonsName{"VetoMuons"};
    TString fMuonMaskName{"GoodMuons"};
    TString fVetoTausName{"GoodTaus"};
    TString fVetoPhotonsName{"VetoPhotons"};
    TString fPhotonMaskName{"GoodPhotons"};
    TString fTriggerBitsName{Names::gkHltBitBrn};

    std::vector<TString> fTriggerNames[nMonoJetCategories]{};
    std::vector<UInt_t> fTriggerIds[nMonoJetCategories]{};

    Bool_t   fCategoryActive[nMonoJetCategories] = {};
    UInt_t   fMaxNumJets[nMonoJetCategories] = {};
    Double_t fMinLeadJetPt[nMonoJetCategories] = {};
    Double_t fMaxJetEta[nMonoJetCategories] = {};
    Double_t fMinMetPt[nMonoJetCategories] = {};
    Double_t fMinChargedHadronFrac[nMonoJetCategories] = {};
    Double_t fMaxNeutralHadronFrac[nMonoJetCategories] = {};
    Double_t fMaxNeutralEmFrac[nMonoJetCategories] = {};

    // Category (signal/calibration regions) bitmask
    mithep::NFArrBool fCategoryFlags{nMonoJetCategories, "MonoJetCategories"};
    
    // Counters
    Int_t fNEventsSelected = 0;

    // Histograms
    TH1D* fCutflow[nMonoJetCategories] = {};   // tally of events passing various cuts
    TH1D* fRealMetPt[nMonoJetCategories] = {}; // met spectrum
    TH1D* fMetPt[nMonoJetCategories] = {};     // met spectrum
    TH1D* fJetPt[nMonoJetCategories] = {};     // jet Pt spectrum
    TH1D* fJetEta[nMonoJetCategories] = {};    // jet eta

    enum Cut {
      cAll,
      cTrigger,
      cNElectrons,
      cNMuons,
      cNTaus,
      cNPhotons,
      cLeadJet,
      cNJets,
      cMet,
      nCuts
    };
    
    ClassDef(MonoJetAnalysisMod,1) // MonJet Selection Module
  };

}
#endif
