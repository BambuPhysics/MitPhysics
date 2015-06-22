#include "MitPhysics/Mods/interface/IdMod.h"

ClassImp(mithep::IdMod)

mithep::IdMod::IdMod(char const* name, char const* title) :
  BaseMod(name, title)
{
  fFlags.SetName(TString(name) + "Flags");
}

mithep::IdMod::~IdMod()
{
  delete fOutput;
}

void
mithep::IdMod::SlaveBegin()
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
mithep::IdMod::SlaveTerminate()
{
  RetractObj(GetOutputName());

  IdTerminate();
}
