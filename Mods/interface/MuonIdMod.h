//--------------------------------------------------------------------------------------------------
// MuonIdMod
//
// Authors:
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_MUONIDMOD_H
#define MITPHYSICS_MODS_MUONIDMOD_H

#include "MitPhysics/Mods/interface/IdMod.h"
#include "MitAna/DataTree/interface/MuonCol.h"

namespace mithep {

  class MuonIdMod : public IdMod {
  public:
    MuonIdMod(char const* name = "MuonIdMod", char const* title = "Muon Identification");
    ~MuonIdMod();

  protected:
    void Process() override;
    void IdBegin() override;
    void IdTerminate() override;

    ClassDef(MuonIdMod, 0)
  };

}

#endif
