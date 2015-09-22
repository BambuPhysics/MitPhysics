//--------------------------------------------------------------------------------------------------
// JetCleaningMod
//
// This Module performs cleaning of jets, ie it removes jets which point 
// in the same direction as a clean isolated electrons.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_JETCLEANINGMOD_H
#define MITPHYSICS_MODS_JETCLEANINGMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/JetCol.h"

namespace mithep {
  class JetCleaningMod : public BaseMod {
  public:
    JetCleaningMod(const char *name="JetCleaningMod", 
                   const char *title="Jet cleaning module");
    ~JetCleaningMod();

    const char *GetCleanElectronsName() const { return fCleanElectronsName; }
    const char *GetCleanMuonsName() const { return fCleanMuonsName; }
    const char *GetCleanJetsName() const { return fCleanJets->GetName(); }
    const char *GetCleanName() const { return GetCleanJetsName(); }
    const char *GetCleanPhotonsName() const { return fCleanPhotonsName; }
    const char *GetCleanTausName() const { return fCleanTausName; }
    const char *GetGoodJetsName() const { return fGoodJetsName; }
    Double_t GetMinDeltaRToElectron() const { return fMinDeltaR[0]; }
    Double_t GetMinDeltaRToMuon() const { return fMinDeltaR[1]; }
    Double_t GetMinDeltaRToTau() const { return fMinDeltaR[2]; }
    Double_t GetMinDeltaRToPhoton() const { return fMinDeltaR[3]; }
    const char *GetOutputName() const { return GetCleanJetsName(); }
    void SetCleanElectronsName(const char *name) { fCleanElectronsName = name; }
    void SetCleanJetsName(const char *name) { fCleanJets->SetName(name); }
    void SetCleanMuonsName(const char *name) { fCleanMuonsName = name; }
    void SetCleanName(const char *name) { SetCleanJetsName(name); }
    void SetCleanPhotonsName(const char *name) { fCleanPhotonsName = name; }
    void SetCleanTausName(const char *name) { fCleanTausName = name; }
    void SetGoodJetsName(const char *name) { fGoodJetsName = name; } 
    void SetMinDeltaRToElectron(Double_t dr) { fMinDeltaR[0] = dr; }
    void SetMinDeltaRToMuon(Double_t dr) { fMinDeltaR[1] = dr; }
    void SetMinDeltaRToTau(Double_t dr) { fMinDeltaR[2] = dr; }
    void SetMinDeltaRToPhoton(Double_t dr) { fMinDeltaR[3] = dr; }
    void SetOutputName(const char *name) { SetCleanJetsName(name); }

  protected:
    void Process() override;
    void SlaveBegin() override;
    void SlaveTerminate() override;
 
    TString fCleanElectronsName; //name of clean electrons (input)
    TString fCleanMuonsName; //name of clean muons (input)
    TString fCleanPhotonsName; //name of clean photons (input)
    TString fCleanTausName; //name of clean taus (input)
    TString fGoodJetsName; //name of good jets (input)
    JetOArr* fCleanJets; //clean jets (output)
    Double_t fMinDeltaR[4]; //e, mu, tau, pho
 
    ClassDef(JetCleaningMod, 1) // Jet cleaning module
  };
}
#endif
