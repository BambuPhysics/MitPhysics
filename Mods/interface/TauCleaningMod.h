//--------------------------------------------------------------------------------------------------
// TauCleaningMod
//
// This Module performs cleaning of taus, ie it removes taus which point 
// in the same direction as a clean isolated muons or electrons.
//
// Authors: G.Ceballos
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_TAUCLEANINGMOD_H
#define MITPHYSICS_MODS_TAUCLEANINGMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CaloTauCol.h"

namespace mithep 
{
  class TauCleaningMod : public BaseMod
  {
    public:
      TauCleaningMod(const char *name="TauCleaningMod", 
                     const char *title="Tau cleaning module");
      ~TauCleaningMod();

      const char      *GetCleanElectronsName()   const { return fCleanElectronsName;   }
      const char      *GetCleanMuonsName()       const { return fCleanMuonsName;       }
      const char      *GetCleanName()            const { return GetCleanCaloTausName();}
      const char      *GetCleanCaloTausName()    const { return fCleanCaloTaus->GetName(); }
      const char      *GetGoodTausName()         const { return fGoodTausName;         }
      Double_t         GetMinDeltaRToElectron()  const { return fMinDeltaRToElectron;  }
      Double_t         GetMinDeltaRToMuon()      const { return fMinDeltaRToMuon;      }
      const char      *GetOutputName()           const { return GetCleanCaloTausName();}
      void             SetCleanElectronsName(const char *name)  { fCleanElectronsName = name; }
      void             SetCleanMuonsName(const char *name)      { fCleanMuonsName     = name; }
      void             SetCleanName(const char *name)           { SetCleanCaloTausName(name); }
      void             SetCleanCaloTausName(const char *name)   { fCleanCaloTaus->SetName(name); }
      void             SetGoodTausName(const char *name)        { fGoodTausName        = name;} 
      void             SetMinDeltaRToElectron(Double_t dr)      { fMinDeltaRToElectron = dr;  }
      void             SetMinDeltaRToMuon(Double_t dr)          { fMinDeltaRToMuon     = dr;  }
      void             SetOutputName(const char *name)          { SetCleanCaloTausName(name); }

    protected:
      void             Process();
      void             SlaveBegin();
      void             SlaveTerminate();

      TString          fCleanElectronsName;  //name of clean electrons (input)
      TString          fCleanMuonsName;      //name of clean muons (input)
      TString          fGoodTausName;        //name of good taus (input)
      Double_t         fMinDeltaRToElectron; //delta R threshold for separating electrons/taus
      Double_t         fMinDeltaRToMuon;     //delta R threshold for separating muons/taus
      CaloTauOArr*     fCleanCaloTaus;
   
    ClassDef(TauCleaningMod, 1) // Tau cleaning module
  };
}
#endif
