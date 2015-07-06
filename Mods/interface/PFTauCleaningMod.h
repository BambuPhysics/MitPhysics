//--------------------------------------------------------------------------------------------------
// PFTauCleaningMod
//
// This Module performs cleaning of taus, ie it removes taus which point 
// in the same direction as a clean isolated muons or electrons.
//
// Authors: G.Ceballos
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PFTauCleaningMod_H
#define MITPHYSICS_MODS_PFTauCleaningMod_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PFTauCol.h"

namespace mithep 
{
  class PFTauCleaningMod : public BaseMod
  {
    public:
      PFTauCleaningMod(const char *name="PFTauCleaningMod", 
                        const char *title="Tau cleaning module");
      ~PFTauCleaningMod();

      const char      *GetCleanElectronsName()   const { return fCleanElectronsName;    }
      const char      *GetCleanMuonsName()       const { return fCleanMuonsName;        }
      const char      *GetCleanName()            const { return GetCleanPFTausName();   }
      const char      *GetCleanPFTausName()      const { return fCleanPFTaus->GetName();}
      const char      *GetGoodPFTausName()       const { return fGoodPFTausName;        }
      Double_t         GetMinDeltaRToElectron()  const { return fMinDeltaRToElectron;   }
      Double_t         GetMinDeltaRToMuon()      const { return fMinDeltaRToMuon;       }
      const char      *GetOutputName()           const { return GetCleanPFTausName();}
      void             SetCleanElectronsName(const char *name)  { fCleanElectronsName = name; }
      void             SetCleanMuonsName(const char *name)      { fCleanMuonsName     = name; }
      void             SetCleanName(const char *name)           { SetCleanPFTausName(name);   }
      void             SetCleanPFTausName(const char *name)     { fCleanPFTaus->SetName(name);}
      void             SetGoodPFTausName(const char *name)      { fGoodPFTausName      = name;} 
      void             SetMinDeltaRToElectron(Double_t dr)      { fMinDeltaRToElectron = dr;  }
      void             SetMinDeltaRToMuon(Double_t dr)          { fMinDeltaRToMuon     = dr;  }
      void             SetOutputName(const char *name)          { SetCleanPFTausName(name);   }

    protected:
      void Process();
      void SlaveBegin();
      void SlaveEnd();

      TString          fCleanElectronsName;  //name of clean electrons (input)
      TString          fCleanMuonsName;      //name of clean muons (input)
      TString          fGoodPFTausName;      //name of good taus (input)
      PFTauOArr*       fCleanPFTaus;     //clean taus (output)
      Double_t         fMinDeltaRToElectron; //delta R threshold for separating electrons/taus
      Double_t         fMinDeltaRToMuon;     //delta R threshold for separating muons/taus
   
    ClassDef(PFTauCleaningMod, 1) // Tau cleaning module
  };
}
#endif
