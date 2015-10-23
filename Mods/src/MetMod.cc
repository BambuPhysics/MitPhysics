#include "MitPhysics/Mods/interface/MetMod.h"

#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/CaloMet.h"
#include "MitAna/DataTree/interface/ParticleCol.h"

ClassImp(mithep::MetMod)

mithep::MetMod::MetMod(char const* name/* = "MetMod"*/, char const* title/* = "Met module"*/) :
  BaseMod(name, title)
{
  fOutput.SetOwner(true);
}

void
mithep::MetMod::SlaveBegin()
{
  mithep::Met* outMet = 0;
  switch (fOutputType) {
  case mithep::kMet:
    outMet = new Met(0., 0.);
    break;
  case mithep::kCaloMet:
    outMet = new CaloMet(0., 0.);
    break;
  case mithep::kPFMet:
    outMet = new PFMet(0., 0.);
    break;
  default:
    SendError(kAbortModule, "SlaveBegin", "Invalid output type");
    return;
  }
  fOutput.AddOwned(outMet);

  fOutput.SetName(fOutputName);
  PublishObj(&fOutput);
}

void
mithep::MetMod::SlaveTerminate()
{
  RetractObj(fOutput.GetName());
}

void
mithep::MetMod::Process()
{
  auto* input = GetObject<mithep::ParticleCol>(fInputName);
  if (!input) {
    SendError(kAbortModule, "Process", "Input collection is null");
    return;
  }

  double mex = 0.;
  double mey = 0.;
  double mez = 0.;
  double sumPt = 0.;

  for (unsigned iP(0); iP != input->GetEntries(); ++iP) {
    auto& part(*input->At(iP));
    mex -= part.Px();
    mey -= part.Py();
    mez -= part.Pz();
    sumPt += part.Pt();
  }

  auto* outMet = fOutput.At(0);
  outMet->SetMex(mex);
  outMet->SetMey(mey);
  outMet->SetSumEt(sumPt);
  outMet->SetElongitudinal(mez);
}
