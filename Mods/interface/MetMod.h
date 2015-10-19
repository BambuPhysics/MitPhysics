//--------------------------------------------------------------------------------------------------
// MetMod
//
// Calculates MET out of the given input collection.
//
// Authors: Y.Iiyama
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_METMOD_H
#define MITPHYSICS_MODS_METMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/MetFwd.h"
#include "MitAna/DataTree/interface/ObjTypes.h"

#include "TString.h"

namespace mithep {

  class MetMod : public BaseMod {
  public:
    MetMod(char const* name = "MetMod", char const* title = "Met module");
    ~MetMod() {}

    void SetInputName(char const* name) { fInputName = name; }
    void SetOutputName(char const* name) { fOutputName = name; }
    void SetOutputType(UInt_t type) { fOutputType = type; }

    char const* GetOutputName() const { return fOutputName; }

  private:
    void SlaveBegin() override;
    void SlaveTerminate() override;
    void Process() override;

    TString fInputName{};
    TString fOutputName{};
    UInt_t fOutputType{mithep::kMet};

    mithep::MetOArr fOutput{};

    ClassDef(MetMod, 0)
  };

}

#endif
