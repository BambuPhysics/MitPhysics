//--------------------------------------------------------------------------------------------------
// MonoXSkimMod
//
// An analysis module for selecting signal and calibration regions for the MonoJet analysis
// and produces some basic distributions.
//
// Authors: LDM, TJW, YI
//--------------------------------------------------------------------------------------------------
#ifndef MITPHYSICS_SKIM_MONOXSKIMMOD_H
#define MITPHYSICS_SKIM_MONOXSKIMMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataCont/interface/Types.h"
#include "MitPhysics/Init/interface/ModNames.h"

#include "TString.h"
#include "TH1D.h"

namespace mithep {

  class MonoXSkimMod : public BaseMod {
  public:
    enum EventCategory {
      kMet,
      kDielectron,
      kDimuon,
      kSingleElectron,
      kSingleMuon,
      kPhoton,
      nEventCategories
    };

    enum AnalysisType {
      kMonoJet,
      kMonoPhoton,
      nAnalysisTypes
    };

    MonoXSkimMod(const char *name  = "MonoXSkimMod", 
                 const char *title = "MoneJet Slection Module");
    ~MonoXSkimMod() {}

    char const* GetCategoryFlagsName() const { return fCategoryFlags.GetName(); }
    
    // set input names
    void SetMetName(const char* n)           { fMetName= n; }
    void SetJetsName(const char* n)          { fJetsName = n; } 
    void SetVetoElectronsName(const char* n) { fVetoElectronsName = n; }
    void SetGoodElectronsName(const char* n) { fGoodElectronsName = n; }
    void SetVetoMuonsName(const char* n)     { fVetoMuonsName = n; }
    void SetGoodMuonsName(const char* n)     { fGoodMuonsName = n; }
    void SetVetoTausName(const char* n)      { fVetoTausName = n; }
    void SetVetoPhotonsName(const char* n)   { fVetoPhotonsName = n; }
    void SetGoodPhotonsName(const char* n)   { fGoodPhotonsName = n; }

    // set output name
    void SetCategoryFlagsName(const char* n) { fCategoryFlags.SetName(n); }

    void SetAnalysisType(UInt_t t) { fAnalysisType = t; }

    // Setting cut values
    void SetMinMonoXPt(Double_t m)  { fMinMonoXPt = m; }
    void SetMinMetPt(Double_t m)    { fMinMetPt = m; }
    void SetVetoElectrons(Bool_t v) { fVetoElectrons = v; }
    void SetVetoMuons(Bool_t v)     { fVetoMuons = v; }
    void SetVetoTaus(Bool_t v)      { fVetoTaus = v; }
    void SetVetoPhotons(Bool_t v)   { fVetoPhotons = v; }

    void SetCategoryActive(UInt_t c, Bool_t a = kTRUE) { fCategoryActive[c] = a; }
    void AddTriggerName(UInt_t c, const char* n)       { fTriggerNames[c].push_back(n); }
    void SetMinNumJets(UInt_t c, Int_t n)              { fMinNumJets[c] = n; }
    void SetMaxNumJets(UInt_t c, Int_t n)              { fMaxNumJets[c] = n; }

    void SetIgnoreTrigger(Bool_t i) { fIgnoreTrigger = i; }

    static TString fgEventCategories[nEventCategories];

  protected:
    // Standard module methods
    void SlaveBegin() override;
    void BeginRun() override;
    void Process() override;
    void SlaveTerminate() override;

    // names of the input collections
    TString fMetName{};
    TString fJetsName{};
    TString fVetoElectronsName{};
    TString fGoodElectronsName{};
    TString fVetoMuonsName{};
    TString fGoodMuonsName{};
    TString fVetoTausName{};
    TString fVetoPhotonsName{};
    TString fGoodPhotonsName{};

    std::vector<TString> fTriggerNames[nEventCategories] = {};
    std::vector<UInt_t> fTriggerIds[nEventCategories] = {};

    UInt_t fAnalysisType{kMonoJet};

    Double_t fMinMonoXPt{0.};
    Double_t fMinMetPt{0.};
    Bool_t   fVetoElectrons{kFALSE};
    Bool_t   fVetoMuons{kFALSE};
    Bool_t   fVetoTaus{kFALSE};
    Bool_t   fVetoPhotons{kFALSE};

    Bool_t fCategoryActive[nEventCategories] = {};
    UInt_t fMinNumJets[nEventCategories] = {};
    UInt_t fMaxNumJets[nEventCategories] = {};

    Bool_t fIgnoreTrigger{kFALSE};

    // Category (signal/calibration regions) bitmask
    mithep::NFArrBool fCategoryFlags;
    
    // Counters
    UInt_t fNEventsSelected{0};

    // Histograms
    TH1D* fCutflow[nEventCategories] = {};   // tally of events passing various cuts

    enum Cut {
      cMonoX,
      cTrigger,
      cNElectrons,
      cNMuons,
      cNTaus,
      cNPhotons,
      cNJets,
      cMet,
      nCuts
    };
    
    ClassDef(MonoXSkimMod, 1) // MonJet Selection Module
  };

}
#endif
