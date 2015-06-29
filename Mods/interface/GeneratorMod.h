//--------------------------------------------------------------------------------------------------
// GeneratorMod
//
// This module collects interesting generator information and publishes collections
// for subsequent modules.
//
// Authors: G.Ceballos
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_GENERATORMOD_H
#define MITPHYSICS_MODS_GENERATORMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MCParticleFwd.h"

class TH1D;
class TH2D;

namespace mithep 
{
  class GeneratorMod : public BaseMod
  {
    public:
      GeneratorMod(const char *name="GeneratorMod", 
                   const char *title="Generator information module");
      ~GeneratorMod();

      Bool_t               GetIsData()           const { return fIsData;           }
      Bool_t               GetApplyISRFilter()   const { return fApplyISRFilter;   }
      Bool_t               GetCopyArrays()	 const { return fCopyArrays;       }	
      const char          *GetMCAllLeptonsName() const { return fMCAllLeptonsName; }	
      const char          *GetMCBosonsName()     const { return fMCBosonsName;     }
      const char          *GetMCISRPhotonsName() const { return fMCISRPhotonsName; }	
      const char          *GetMCLeptonsName()    const { return fMCLeptonsName;    }	
      const char          *GetMCMETName()        const { return fMCMETName;        }	
      const char          *GetMCNeutrinosName()  const { return fMCNeutrinosName;  }
      const char          *GetMCPartName()	 const { return fMCPartName;       } 	
      const char          *GetMCPhotonsName()    const { return fMCPhotonsName;    }
      const char          *GetMCQuarksName()     const { return fMCQuarksName;     }	
      const char          *GetMCRadPhotonsName() const { return fMCRadPhotonsName; }	
      const char          *GetMCTausName()	 const { return fMCTausName;       }	
      const char          *GetMCqqHsName()	 const { return fMCqqHsName;       }	
      Bool_t               GetPrintDebug()	 const { return fPrintDebug;       }	
      void                 SetIsData(Bool_t b)	               { fIsData           = b; }	 
      void                 SetApplyISRFilter(Bool_t b)	       { fApplyISRFilter   = b; }	 
      void                 SetApplyVVFilter(Bool_t b)	       { fApplyVVFilter    = b; }	 
      void                 SetApplyVGFilter(Bool_t b)	       { fApplyVGFilter    = b; }	 
      void                 SetAllowWWEvents(Bool_t b)	       { fAllowWWEvents    = b; }	 
      void                 SetAllowWZEvents(Bool_t b)	       { fAllowWZEvents    = b; }	 
      void                 SetAllowZZEvents(Bool_t b)	       { fAllowZZEvents    = b; }	       
      void                 SetFilterBTEvents(Bool_t b)	       { fFilterBTEvents   = b; }	       
      void                 SetCopyArrays(Bool_t b)	       { fCopyArrays       = b; }	 
      void                 SetEtaLeptonMax(Double_t x)         { fEtaLeptonMax     = x; }     
      void                 SetEtaPhotonMax(Double_t x)         { fEtaPhotonMax     = x; }	  
      void                 SetEtaRadPhotonMax(Double_t x)      { fEtaRadPhotonMax  = x; }	  
      void                 SetMCAllLeptonsName(const char * s) { fMCAllLeptonsName = s; }	
      void                 SetMCBosonsName(const char *s)      { fMCBosonsName     = s; }	
      void                 SetMCISRPhotonsName(const char *s)  { fMCISRPhotonsName = s; }	
      void                 SetMCLeptonsName(const char * s)    { fMCLeptonsName    = s; }	
      void                 SetMCMETName(const char * s)        { fMCMETName        = s; }	
      void                 SetMCNeutrinosName(const char *s)   { fMCNeutrinosName  = s; }	
      void                 SetMCPartName(const char *s)	       { fMCPartName       = s; }	
      void                 SetMCPhotonsName(const char *s)     { fMCPhotonsName    = s; }	
      void                 SetMCQuarksName(const char *s)      { fMCQuarksName     = s; }	
      void                 SetMCRadPhotonsName(const char *s)  { fMCRadPhotonsName = s; }	
      void                 SetMCTausName(const char *s)	       { fMCTausName       = s; }	
      void                 SetMCqqHsName(const char *s)	       { fMCqqHsName       = s; }	
      void                 SetMassMaxCut(Double_t x)	       { fMassMaxCut	   = x; }	 
      void                 SetMassMinCut(Double_t x)	       { fMassMinCut       = x; }	 
      void                 SetPdgIdCut(UInt_t d)	       { fPdgIdCut         = d; }	 
      void                 SetPrintDebug(bool b)               { fPrintDebug       = b; }   
      void                 SetPtLeptonMin(Double_t x)          { fPtLeptonMin      = x; }     
      void                 SetPtPhotonMin(Double_t x)          { fPtPhotonMin      = x; }     
      void                 SetPtRadPhotonMin(Double_t x)       { fPtRadPhotonMin   = x; }     

    protected:
      void                 Process();
      void                 SlaveBegin();

      Bool_t               fIsData;             //=true then it does nothing (def=0)
      Bool_t               fPrintDebug;         //=true then print debug info (def=0)
      Bool_t               fCopyArrays;         //=true then copy array content for skimming  (def=0)
      TString              fMCPartName;         //name of MCParticle branch
      TString              fMCMETName;          //name of met coll
      TString              fMCLeptonsName;      //name of lepton coll (from W/Z/H)
      TString              fMCAllLeptonsName;   //name of lepton coll (all)
      TString              fMCTausName;         //name of tau coll (hadronic decays)
      TString              fMCNeutrinosName;    //name of neutrinos coll
      TString              fMCQuarksName;       //name of quarks coll
      TString              fMCqqHsName;         //name of qqH coll
      TString              fMCBosonsName;       //name of bosons coll
      TString              fMCPhotonsName;      //name of photons coll
      TString              fMCRadPhotonsName;   //name of rad photons coll
      TString              fMCISRPhotonsName;   //name of ISR photons coll
      Double_t             fPtLeptonMin;        //pt min for leptons
      Double_t             fEtaLeptonMax;       //eta max for leptons
      Double_t             fPtPhotonMin;        //pt min for photons
      Double_t             fEtaPhotonMax;       //eta max for photons
      Double_t             fPtRadPhotonMin;     //pt min for rad photons
      Double_t             fEtaRadPhotonMax;    //eta max for rad photons
      UInt_t               fPdgIdCut;           //pdg id for particle used to select on mass (def=0)
      Double_t             fMassMinCut;	        //mass min for given PdgId particle 
      Double_t             fMassMaxCut;	        //mass max for given PdgId particle
      Bool_t               fApplyISRFilter;     //=true then apply ISR filter (def=0)
      Bool_t               fApplyVVFilter;      //=true then apply VV filters (def=0)
      Bool_t               fApplyVGFilter;      //=true then apply VG filters (def=0)
      Bool_t               fAllowWWEvents;            //=true then allow WW events  (def=0)
      Bool_t               fAllowWZEvents;            //=true then allow WZ events  (def=0)
      Bool_t               fAllowZZEvents;            //=true then allow ZZ events  (def=0)
      Bool_t               fFilterBTEvents;           //=true then filter B and T events  (def=0)
      TH1D                *hDGenPtMin;          //!histo for Pt min for leptons from W/Z
      TH1D                *hDGenMet[10];        //!histos for gen MET
      TH1D                *hDGenLeptons[40];    //!histos for W/Z/H leptons
      TH1D		  *hDGenWMuons ;        //!histo for # of W muons only
      TH1D		  *hDGenWElectrons;     //!histo for # of W electrons only
      TH1D                *hDGenAllLeptons[20]; //!histos for all leptons
      TH1D                *hDGenTaus[20];       //!histos for taus
      TH1D                *hDGenNeutrinos[20];  //!histos for neutrinos
      TH1D                *hDGenQuarks[20];     //!histos for quarks
      TH1D                *hDGenWBF[20];        //!histos for WBF
      TH1D                *hDGenBosons[15];     //!histos for bosons
      TH1D                *hDGenPhotons[20];    //!histos for photons
      TH1D                *hDGenRadPhotons[20]; //!histos for rad photons
      TH1D                *hDGenISRPhotons[20]; //!histos for ISR photons
      TH1D                *hDVMass[20];         //!histos for auxiliar MG work
      TH1D                *hDVVMass[50];        //!histos for auxiliar VV work
      TH1D                *hDTopMass[2];        //!histos for auxiliar Top work
      MCParticleArr       *fGenLeptons;         //!copied owning array for skimming
      MCParticleArr       *fGenAllLeptons;      //!copied owning array for skimming
      MCParticleArr       *fGenTaus;            //!copied owning array for skimming
      MCParticleArr       *fGenNeutrinos;       //!copied owning array for skimming
      MCParticleArr       *fGenQuarks;          //!copied owning array for skimming
      MCParticleArr       *fGenqqHs;            //!copied owning array for skimming
      MCParticleArr       *fGenBosons;          //!copied owning array for skimming
      MCParticleArr       *fGenPhotons;         //!copied owning array for skimming
      MCParticleArr       *fGenRadPhotons;      //!copied owning array for skimming
      MCParticleArr       *fGenISRPhotons;      //!copied owning array for skimming

    ClassDef(GeneratorMod, 1) // Module to gather generator information
  };
}
#endif
