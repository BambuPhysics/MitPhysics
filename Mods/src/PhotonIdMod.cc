#include "MitPhysics/Mods/interface/PhotonIdMod.h"
#include "MitAna/DataTree/interface/Names.h"

ClassImp(mithep::PhotonIdMod)

mithep::PhotonIdMod::PhotonIdMod(char const* name/* = "PhotonIdMod"*/, char const* title/* = "Photon Identification"*/) :
  IdMod<mithep::Photon>(name, title)
{
  fInputName = mithep::Names::gkPhotonBrn;
}

mithep::PhotonIdMod::~PhotonIdMod()
{
}

Bool_t
mithep::PhotonIdMod::IsGood(mithep::Photon const& photon)
{
  fCutFlow->Fill(cAll);

  if (photon.Pt() < fPtMin)
    return false;

  fCutFlow->Fill(cPt);

  if (photon.AbsEta() > fEtaMax)
    return false;

  fCutFlow->Fill(cEta);

  if (GetApplyPixelVeto() && photon.HasPixelSeed())
    return false;

  if (GetApplyElectronVeto() && !PhotonTools::PassElectronVeto(&photon, GetElectrons()))
    return false;

  if (GetApplyCSafeElectronVeto() && !PhotonTools::PassElectronVetoConvRecovery(&photon, GetElectrons(), GetConversions(), GetBeamSpot()->At(0)))
    return false;

  fCutFlow->Fill(cElectronVeto);

  //apply trigger matching
  if (fApplyTriggerMatching &&
      !PhotonTools::PassTriggerMatching(&photon, GetTrigObjects()))
    return false;

  fCutFlow->Fill(cTriggerMatching);

  //apply id cut
  if (!PassIdCut(photon))
    return false;

  fCutFlow->Fill(cId);

  //apply Isolation Cut
  if (!PassIsolationCut(photon))
    return false;

  fCutFlow->Fill(cIsolation);

  return true;
}

void
mithep::PhotonIdMod::IdBegin()
{
  fCutFlow->SetBins(nCuts, 0., double(nCuts));
  TAxis* xaxis = fCutFlow->GetXaxis();
  xaxis->SetBinLabel(cAll + 1, "All");
  xaxis->SetBinLabel(cPt + 1,  "Pt");
  xaxis->SetBinLabel(cEta + 1, "Eta");
  xaxis->SetBinLabel(cElectronVeto + 1, "ElectronVeto");
  xaxis->SetBinLabel(cTriggerMatching + 1, "TriggerMatching");
  xaxis->SetBinLabel(cId + 1,  "Id");
  xaxis->SetBinLabel(cIsolation + 1, "Isolation");
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::PhotonIdMod::PassIdCut(Photon const& pho)
{
  switch (fIdType) {
  case PhotonTools::kPhys14Tight:
  case PhotonTools::kPhys14Medium:
  case PhotonTools::kPhys14Loose:
  case PhotonTools::kHighPtV2:
  case PhotonTools::kSpring15Tight50ns:
  case PhotonTools::kSpring15Medium50ns:
  case PhotonTools::kSpring15Loose50ns:
  case PhotonTools::kSpring15Tight:
  case PhotonTools::kSpring15Medium:
  case PhotonTools::kSpring15Loose:
    return PhotonTools::PassID(&pho, PhotonTools::EPhIdType(fIdType));

  case PhotonTools::kNoId:
    return true;

  default:
    return false;
  }
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::PhotonIdMod::PassIsolationCut(Photon const& pho)
{
  switch (fIsoType) {
  case PhotonTools::kPhys14LooseIso:
  case PhotonTools::kPhys14MediumIso:
  case PhotonTools::kPhys14TightIso:
  case PhotonTools::kHighPtV2Iso:
  case PhotonTools::kSpring15Loose50nsIso:
  case PhotonTools::kSpring15Medium50nsIso:
  case PhotonTools::kSpring15Tight50nsIso:
  case PhotonTools::kSpring15LooseIso:
  case PhotonTools::kSpring15MediumIso:
  case PhotonTools::kSpring15TightIso:
    return PhotonTools::PassIsoFootprintRhoCorr(&pho, PhotonTools::EPhIsoType(fIsoType),
                                                GetVertices()->At(0), GetPFCandidates(),
                                                GetPileupEnergyDensity()->At(0)->Rho(fRhoAlgo));

  case PhotonTools::kNoIso:
    return true;

  default:
    return false;
  }
}
