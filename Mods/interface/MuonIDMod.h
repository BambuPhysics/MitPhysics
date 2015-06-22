//--------------------------------------------------------------------------------------------------
// MuonIDMod
//
// Authors: Y.Iiyama
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

    void SetOutputName(char const* n) override { static_cast<MuonOArr*>(fOutput)->SetName(n); }

  protected:
    void Process() override;
    void SlaveBegin() override;
    void SlaveTerminate() override;

    ClassDef(MuonIDMod, 0)
  };

}

#endif
