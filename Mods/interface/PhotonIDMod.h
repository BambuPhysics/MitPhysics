//--------------------------------------------------------------------------------------------------
// PhotonIDMod
//
// Authors: Y.Iiyama
//--------------------------------------------------------------------------------------------------


#ifndef MITPHYSICS_MODS_PHOTONIDMOD_H
#define MITPHYSICS_MODS_PHOTONIDMOD_H

#include "MitPhysics/Mods/interface/IDMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"

namespace mithep {

  class PhotonIDMod : public IDMod {
  public:
    PhotonIDMod(char const* name = "PhotonIDMod", char const* title = "Photon Identification");
    ~PhotonIDMod();

  protected:
    void Process() override;
    void IdBegin() override;
    void IdTerminate() override;

    ClassDef(PhotonIDMod, 0)
  };

}

#endif
