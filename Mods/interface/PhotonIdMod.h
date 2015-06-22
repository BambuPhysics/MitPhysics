//--------------------------------------------------------------------------------------------------
// PhotonIdMod
//
// Authors: Y.Iiyama
//--------------------------------------------------------------------------------------------------


#ifndef MITPHYSICS_MODS_PHOTONIDMOD_H
#define MITPHYSICS_MODS_PHOTONIDMOD_H

#include "MitPhysics/Mods/interface/IdMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"

namespace mithep {

  class PhotonIdMod : public IdMod {
  public:
    PhotonIdMod(char const* name = "PhotonIdMod", char const* title = "Photon Identification");
    ~PhotonIdMod();

  protected:
    void Process() override;
    void IdBegin() override;
    void IdTerminate() override;

    ClassDef(PhotonIdMod, 0)
  };

}

#endif
