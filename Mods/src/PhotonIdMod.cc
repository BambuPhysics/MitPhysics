#include "MitPhysics/Mods/interface/PhotonIdMod.h"

ClassImp(mithep::PhotonIdMod)

mithep::PhotonIdMod::PhotonIdMod(char const* name/* = "PhotonIdMod"*/, char const* title/* = "Photon Identification"*/) :
  IdMod(name, title)
{
  fOutput = new PhotonOArr(32, TString(name) + "Output");
}

mithep::PhotonIdMod::~PhotonIdMod()
{
}

void
mithep::PhotonIdMod::Process()
{
}

void
mithep::PhotonIdMod::IdBegin()
{
}

void
mithep::PhotonIdMod::IdTerminate()
{
}
