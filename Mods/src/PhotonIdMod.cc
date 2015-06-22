#include "MitPhysics/Mods/interface/PhotonIDMod.h"

ClassImp(mithep::PhotonIDMod)

mithep::PhotonIDMod::PhotonIDMod(char const* name/* = "PhotonIDMod"*/, char const* title/* = "Photon Identification"*/) :
  IDMod(name, title)
{
  fOutput = new PhotonOArr(32, TString(name) + "Output");
}

mithep::PhotonIDMod::~PhotonIDMod()
{
}

void
mithep::PhotonIDMod::Process()
{
}

void
mithep::PhotonIDMod::IdBegin()
{
}

void
mithep::PhotonIDMod::IdTerminate()
{
}
