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
  case PhotonTools::kSummer15Tight:
  case PhotonTools::kSummer15Medium:
  case PhotonTools::kSummer15Loose:
    return PhotonTools::PassID(&pho, PhotonTools::EPhIdType(fIdType));

  default:
    return false;
  }
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::PhotonIdMod::PassIsolationCut(Photon const& pho)
{
  switch (fIsoType) {
  case PhotonTools::kSummer15LooseIso:
  case PhotonTools::kSummer15MediumIso:
  case PhotonTools::kSummer15TightIso:
    return PhotonTools::PassIsoFootprintRhoCorr(&pho, PhotonTools::EPhIsoType(fIsoType),
                                                GetVertices()->At(0), GetPFCandidates(),
                                                GetPileupEnergyDensity()->At(0)->Rho(fRhoAlgo));

  case PhotonTools::kNoIso:
    return true;

  default:
    return false;
  }
}
