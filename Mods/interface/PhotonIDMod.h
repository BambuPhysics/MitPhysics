#ifndef MITPHYSICS_MODS_PHOTONIDMOD_H
#define MITPHYSICS_MODS_PHOTONIDMOD_H

#include "MitPhysics/Mods/interface/IDMod.h"
#include "MitAna/DataTree/interface/PhotonCol.h"

namespace mithep {

  class PhotonIDMod : public IDMod {
  public:
    PhotonIDMod(char const* name = "PhotonIDMod", char const* title = "Photon Identification");
    ~PhotonIDMod();

    void SetOutputName(char const* n) override { static_cast<PhotonOArr*>(fOutput)->SetName(n); }

  protected:
    void Process() override;
    void SlaveBegin() override;
    void SlaveTerminate() override;

    ClassDef(PhotonIDMod, 0)
  };

}

#endif
