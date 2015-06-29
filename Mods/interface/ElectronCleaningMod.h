//--------------------------------------------------------------------------------------------------
// $Id: ElectronCleaningMod.h,v 1.6 2009/03/23 14:23:06 loizides Exp $
//
// ElectronCleaningMod
//
// This Module performs cleaning of electrons, ie. it removes duplicate objects and good muons 
// from the good electrons.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_ELECTRONCLEANINGMOD_H
#define MITPHYSICS_MODS_ELECTRONCLEANINGMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/ElectronCol.h"

namespace mithep 
{
  class ElectronCleaningMod : public BaseMod
  {
    public:
      ElectronCleaningMod(const char *name="ElectronCleaningMod", 
                          const char *title="Electron cleaning module");

      const char        *GetCleanElectronsName() const { return fCleanElectrons->GetName();     }
      const char        *GetCleanName()          const { return GetCleanElectronsName(); }
      const char        *GetCleanMuonsName()     const { return fCleanMuonsName;         }
      const char        *GetGoodElectronsName()  const { return fGoodElectronsName;      }
      const char        *GetOutputName()         const { return GetCleanElectronsName(); }
      void               SetCleanElectronsName(const char *name) { fCleanElectrons->SetName(name);  }
      void               SetCleanName(const char *name)          { SetCleanElectronsName(name); }
      void               SetCleanMuonsName(const char *name)     { fCleanMuonsName     = name;  }
      void               SetGoodElectronsName(const char *name)  { fGoodElectronsName  = name;  }
      void               SetOutputName(const char *name)         { SetCleanElectronsName(name); }

    protected:
      void               Process();
      void               SlaveBegin();
      void               SlaveEnd();

      TString            fGoodElectronsName;  //name of good electrons (input)
      TString            fCleanMuonsName;     //name of clean muons (input)
      ElectronOArr*      fCleanElectrons; //clean electrons (output)
    
      ClassDef(ElectronCleaningMod, 2) // Electron cleaning module
  };
}
#endif
