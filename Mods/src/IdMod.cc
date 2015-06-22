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

void
mithep::IDMod::SlaveBegin()
{
  if (fIsFilterMode) {
    if (!PublishObj(fOutput))
      SendError(kAbortAnalysis, "SlaveBegin", "Cannot publish output");
  }
  else {
    if (!PublishObj(&fFlags))
      SendError(kAbortAnalysis, "SlaveBegin", "Cannot publish output");
  }

  AddTH1(fCutFlow, TString(GetName()) + "CutFlow", "Identification cut flow", 1, 0., 1.);

  IdBegin();
}

void
mithep::IDMod::SlaveTerminate()
{
  RetractObj(GetOutputName());

  IdTerminate();
}
