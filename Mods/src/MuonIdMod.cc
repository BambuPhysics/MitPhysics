#include "MitPhysics/Mods/interface/MuonIdMod.h"

ClassImp(mithep::MuonIdMod)

mithep::MuonIdMod::MuonIdMod(char const* name/* = "MuonIdMod"*/, char const* title/* = "Muon Identification"*/) :
  IdMod(name, title)
{
  fOutput = new MuonOArr(32, TString(name) + "Output");
}

mithep::MuonIdMod::~MuonIdMod()
{
}

void
mithep::MuonIdMod::Process()
{
}

void
mithep::MuonIdMod::IdBegin()
{
}

void
mithep::MuonIdMod::IdTerminate()
{
}
