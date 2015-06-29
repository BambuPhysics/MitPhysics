//--------------------------------------------------------------------------------------------------
// $Id: PhotonCleaningMod.h,v 1.5 2009/03/23 14:23:06 loizides Exp $
//
// PhotonCleaningMod
//
// This Module performs cleaning of jets, ie it removes jets which point 
// in the same direction as a clean isolated electrons.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_PHOTONCLEANINGMOD_H
#define MITPHYSICS_MODS_PHOTONCLEANINGMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/PhotonCol.h"

namespace mithep 
{
  class PhotonCleaningMod : public BaseMod
  {
    public:
      PhotonCleaningMod(const char *name="PhotonCleaningMod", 
                        const char *title="Photon cleaning module");

      const char      *GetCleanElectronsName()   const { return fCleanElectronsName;   }
      const char      *GetCleanName()            const { return GetCleanPhotonsName(); }
      const char      *GetCleanPhotonsName()     const { return fCleanPhotons->GetName();     }
      const char      *GetGoodPhotonsName()      const { return fGoodPhotonsName;      } 
      Double_t         GetMinDeltaRToElectron()  const { return fMinDeltaRToElectron;  }
      const char      *GetOutputName()           const { return GetCleanPhotonsName(); }
      void             SetCleanElectronsName(const char *name)  { fCleanElectronsName = name;}
      void             SetCleanName(const char *name)           { SetCleanPhotonsName(name);     }
      void             SetCleanPhotonsName(const char *name)    { fCleanPhotons->SetName(name); }
      void             SetGoodPhotonsName(const char *name)     { fGoodPhotonsName       = name; } 
      void             SetMinDeltaRToElectron(Double_t dr)      { fMinDeltaRToElectron   = dr;   }
      void             SetOutputName(const char *name)          { SetCleanPhotonsName(name);     }

    protected:
      void             Process();
      void             SlaveBegin();
      void             SlaveEnd();

      TString          fCleanElectronsName;   //name of clean electrons (input)
      TString          fGoodPhotonsName;      //name of good jets (input)
      PhotonOArr*      fCleanPhotons;     //clean photons (output)
      Double_t         fMinDeltaRToElectron;  //delta R threshold for separating electrons/photons
   
    ClassDef(PhotonCleaningMod, 1) // Photon cleaning module
  };
}
#endif
