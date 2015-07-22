//--------------------------------------------------------------------------------------------------
// EventCategoryMod
//
// Generic event category filter. Takes a NamedFastArrayBasic<Bool_t> and passes events if the
// specified index is true.
//
// Authors: Y.Iiyama
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_MODS_EVENTCATEGORYMOD_H
#define MITPHYSICS_MODS_EVENTCATEGORYMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"

#include "TString.h"

namespace mithep {

  class EventCategoryMod : public BaseMod {
  public:
    EventCategoryMod(char const* name = "EventCategoryMod", char const* title = "event category filter") : BaseMod(name, title) {}
    ~EventCategoryMod() {}

    char const* GetInputName() const { return fInputName; }
    UInt_t      GetCategory() const { return fCategory; }

    void SetInputName(char const* n) { fInputName = n; }
    void SetCategory(UInt_t c) { fCategory = c; }

  protected:
    void Process() override;

    TString fInputName{};
    UInt_t fCategory{0xffffffff};

    ClassDef(EventCategoryMod, 0)
  };

}

#endif
