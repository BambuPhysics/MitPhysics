//--------------------------------------------------------------------------------------------------
// MuonIDMod
//
// Authors:
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_MUONIDMOD_H
#define MITPHYSICS_MODS_MUONIDMOD_H

#include "MitPhysics/Mods/interface/IDMod.h"
#include "MitAna/DataTree/interface/MuonCol.h"

namespace mithep {

  class MuonIDMod : public IDMod {
  public:
    MuonIDMod(char const* name = "MuonIDMod", char const* title = "Muon Identification");
    ~MuonIDMod();

  protected:
    void Process() override;
    void IdBegin() override;
    void IdTerminate() override;

    ClassDef(MuonIDMod, 0)
  };

}

#endif
