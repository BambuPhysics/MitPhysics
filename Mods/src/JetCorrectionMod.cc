#include "MitPhysics/Mods/interface/JetCorrectionMod.h"

#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"

using namespace mithep;

ClassImp(mithep::JetCorrectionMod)

//--------------------------------------------------------------------------------------------------
JetCorrectionMod::JetCorrectionMod(const char *name, const char *title) : 
  BaseMod(name, title)
{
  fCorrectedJets.SetOwner(true);
}

void
mithep::JetCorrectionMod::SetRhoType(RhoUtilities::RhoType type)
{
  // DEPRECATED FUNCTION
  // Use SetRhoAlgo instead

  switch(type) {
  case RhoUtilities::MIT_RHO_VORONOI_LOW_ETA:
    fRhoAlgo = mithep::PileupEnergyDensity::kLowEta;
    break;
  case RhoUtilities::MIT_RHO_VORONOI_HIGH_ETA:
    fRhoAlgo = mithep::PileupEnergyDensity::kHighEta;
    break;
  case RhoUtilities::MIT_RHO_RANDOM_LOW_ETA:
    fRhoAlgo = mithep::PileupEnergyDensity::kRandomLowEta;
    break;
  case RhoUtilities::MIT_RHO_RANDOM_HIGH_ETA:
    fRhoAlgo = mithep::PileupEnergyDensity::kRandom;
    break;
  case RhoUtilities::CMS_RHO_RHOKT6PFJETS:
    fRhoAlgo = mithep::PileupEnergyDensity::kKt6PFJets;
    break;
  default:
    fRhoAlgo = mithep::PileupEnergyDensity::kHighEta;
    break;
  }
}   

//--------------------------------------------------------------------------------------------------
void
JetCorrectionMod::SlaveBegin()
{
  // fCorrector = 0 is alllowed but should not happen here - the module is useless in this case
  MakeCorrector();

  PublishObj(&fCorrectedJets);

  if (fCorrector->HasSmearing()) {
    if (fGenJetsName == "")
      Warning("SlaveBegin", "Applying JER smearing without gen jets.");

    fCorrector->SetSmearingSeed(fSmearingSeed);
  }
}

void JetCorrectionMod::SlaveTerminate()
{
  RetractObj(fCorrectedJets.GetName());

  if (fOwnCorrector) {
    delete fCorrector;
    fCorrector = 0;
  }
}

//--------------------------------------------------------------------------------------------------
void 
JetCorrectionMod::Process()
{
  // Process entries of the tree. 
  auto* inJets = GetObject<JetCol>(fJetsName);
  if (!inJets) {
    SendError(kAbortModule, "Process", 
              "Pointer to input jet collection %s is null.",
              fJetsName.Data());
    return;
  }

  // get the energy density from the event
  auto* rhocol = GetObject<PileupEnergyDensityCol>(fRhoBranchName);
  double rho = 0.;
  if (rhocol) {
    rho = rhocol->At(0)->Rho(fRhoAlgo);
  }
  else if (fCorrector->IsEnabled(Jet::L1)) {
    SendError(kAbortModule, "Process", 
              "Pointer to input rho collection %s is null.",
              fRhoBranchName.Data());
    return;
  }

  GenJetCol const* genJets = 0;
  if (fCorrector->HasSmearing() && fGenJetsName != "")
    genJets = GetObject<GenJetCol>(fGenJetsName);

  fCorrectedJets.Reset();

  // loop over jets
  for (unsigned iJ = 0; iJ != inJets->GetEntries(); ++iJ) {
    //copy input jet, using special function to copy full derived class
    Jet* jet = inJets->At(iJ)->MakeCopy();

    fCorrector->Correct(*jet, rho);

    if (fCorrector->HasSmearing())
      fCorrector->Smear(*jet, rho, genJets);

    // add corrected jet to collection
    fCorrectedJets.AddOwned(jet);             
  }

  // sort according to ptrootcint forward declaration data members
  fCorrectedJets.Sort();
}

//--------------------------------------------------------------------------------------------------
void
JetCorrectionMod::AddCorrectionFromFile(char const* fileName, mithep::JetCorrector::Corrector::FactorType factorType/* = nFactortypes*/)
{
  MakeCorrector();
  fCorrector->AddParameterFile(fileName, factorType);
}

//--------------------------------------------------------------------------------------------------
void
JetCorrectionMod::SetCorrector(JetCorrector* corr)
{
  if (fCorrector && fOwnCorrector) {
    SendError(kAbortAnalysis, "SetCorrector", "Corrector already created");
  }
  else {
    fCorrector = corr;
  }
}

//--------------------------------------------------------------------------------------------------
void
JetCorrectionMod::MakeCorrector()
{
  if (!fCorrector) {
    fCorrector = new JetCorrector;
    fCorrector->SetUncertaintySigma(fSigma);
    fOwnCorrector = true;
  }
}
