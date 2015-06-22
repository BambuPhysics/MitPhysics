#include "MitPhysics/Mods/interface/MuonIDMod.h"

ClassImp(mithep::MuonIDMod)

mithep::MuonIDMod::MuonIDMod(char const* name/* = "MuonIDMod"*/, char const* title/* = "Muon Identification"*/) :
  IDMod(name, title)
{
  fOutput = new MuonOArr(32, TString(name) + "Output");
}

mithep::MuonIDMod::~MuonIDMod()
{
}

void
mithep::MuonIDMod::Process()
{
}

void
mithep::MuonIDMod::IdBegin()
{
}

void
mithep::MuonIDMod::IdTerminate()
{
}
