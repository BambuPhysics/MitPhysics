#include "MitPhysics/Mods/interface/IDMod.h"

ClassImp(mithep::IDMod)

mithep::IDMod::IDMod(char const* name, char const* title) :
  BaseMod(name, title)
{
  fFlags.SetName(TString(name) + "Flags");
}

mithep::IDMod::~IDMod()
{
  delete fOutput;
}

Bool_t
mithep::IDMod::PublishOutput()
{
  if (!PublishObj(fOutput))
    return kFALSE;

  if (!fIsFilterMode) {
    if (!PublishObj(&fFlags))
      return kFALSE;
  }

  return kTRUE;
}

void
mithep::IDMod::RetractOutput()
{
  RetractObj(GetOutputName());
  if (!fIsFilterMode)
    RetractObj(GetFlagsName());
}
