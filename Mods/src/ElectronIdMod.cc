#include "MitPhysics/Mods/interface/ElectronIdMod.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Init/interface/Constants.h"

ClassImp(mithep::ElectronIdMod)

//--------------------------------------------------------------------------------------------------
mithep::ElectronIdMod::ElectronIdMod(const char *name, const char *title) :
  IdMod<mithep::Electron>(name, title)
{
}

//--------------------------------------------------------------------------------------------------
void
mithep::ElectronIdMod::IdBegin()
{
  // If we use MVA Id, need to load MVA weights
  switch (fIdType) {
  case ElectronTools::kLikelihood:
    if (!fLH)
      SendError(kAbortAnalysis, "SlaveBegin", "Electron likelihood is not initialized.");
    break;

  default:
    break;
  }

  switch (fIsoType) {
  case ElectronTools::kCustomIso:
    SendError(kAbortAnalysis, "SlaveBegin", "Custom electron isolation is not yet implemented.");
    break;

  default:
    break;
  }

  fCutFlow->SetBins(nCuts, 0., double(nCuts));
  TAxis* xaxis = fCutFlow->GetXaxis();
  xaxis->SetBinLabel(cAll + 1, "All");
  xaxis->SetBinLabel(cIsEcalDriven + 1, "IsEcalDriven");
  xaxis->SetBinLabel(cPt + 1, "Pt");
  xaxis->SetBinLabel(cEta + 1, "Eta");
  xaxis->SetBinLabel(cEt + 1, "Et");
  xaxis->SetBinLabel(cFiducial + 1, "Fiducial");
  xaxis->SetBinLabel(cSpikeRemoval + 1, "SpikeRemoval");
  xaxis->SetBinLabel(cChargeFilter + 1, "ChargeFilter");
  xaxis->SetBinLabel(cTriggerMatching + 1, "TriggerMatching");
  xaxis->SetBinLabel(cConvFilterType1 + 1, "ConvFilterType1");
  xaxis->SetBinLabel(cConvFilterType2 + 1, "ConvFilterType2");
  xaxis->SetBinLabel(cNExpectedHits + 1, "NExpectedHits");
  xaxis->SetBinLabel(cD0 + 1, "D0");
  xaxis->SetBinLabel(cDZ + 1, "DZ");
  xaxis->SetBinLabel(cId + 1, "Id");
  xaxis->SetBinLabel(cIsolation + 1, "Isolation");
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::ElectronIdMod::IsGood(mithep::Electron const& electron)
{
  fCutFlow->Fill(cAll);

  if (electron.SCluster() == 0)
    return false;

  if (fApplyEcalSeeded && !electron.IsEcalDriven())
    return false;

  fCutFlow->Fill(cIsEcalDriven);

  if (electron.Pt() < fPtMin)
    return false;

  fCutFlow->Fill(cPt);

  if (electron.AbsEta() > fEtaMax)
    return false;

  fCutFlow->Fill(cEta);

  if (electron.SCluster()->Et() < fElectronEtMin)
    return false;

  fCutFlow->Fill(cEt);

  double scAbsEta = electron.SCluster()->AbsEta();

  if (fApplyEcalFiducial &&
      ((scAbsEta > gkEleEBEtaMax && scAbsEta < gkEleEEEtaMin) || scAbsEta > gkEleEEEtaMax))
    return false;

  fCutFlow->Fill(cFiducial);

  //apply ECAL spike removal
  if (fApplySpikeRemoval && !ElectronTools::PassSpikeRemovalFilter(&electron))
    return false;

  fCutFlow->Fill(cSpikeRemoval);

  // apply charge filter
  if (fChargeFilter && !ElectronTools::PassChargeFilter(&electron))
    return false;

  fCutFlow->Fill(cChargeFilter);

  //apply trigger matching
  if (fApplyTriggerMatching &&
      !ElectronTools::PassTriggerMatching(&electron, GetTrigObjects()))
    return false;

  fCutFlow->Fill(cTriggerMatching);

  // apply conversion filters
  if (fApplyConvFilterType1 &&
      !ElectronTools::PassConversionFilter(&electron,
                                           GetConversions(),
                                           GetBeamSpot()->At(0),
                                           0, 1e-6, 2.0, true, false))
    return false;

  fCutFlow->Fill(cConvFilterType1);

  if (fApplyConvFilterType2 && TMath::Abs(electron.ConvPartnerDCotTheta()) < 0.02 && TMath::Abs(electron.ConvPartnerDist()) < 0.02)
    return false;

  fCutFlow->Fill(cConvFilterType2);

  if (fApplyNExpectedHitsInnerCut && !ElectronTools::PassNExpectedHits(&electron, ElectronTools::EElIdType(fIdType), fInvertNExpectedHitsInnerCut))
    return false;

  fCutFlow->Fill(cNExpectedHits);

  // apply d0 cut
  if (fApplyD0Cut) {
    if (fWhichVertex >= -1) {
      if (!ElectronTools::PassD0Cut(&electron, GetVertices(), ElectronTools::EElIdType(fIdType), fWhichVertex))
        return false;
    }
    else if (!ElectronTools::PassD0Cut(&electron, GetBeamSpot(), ElectronTools::EElIdType(fIdType)))
      return false;
  }

  fCutFlow->Fill(cD0);

  // apply dz cut
  if (fApplyDZCut && !ElectronTools::PassDZCut(&electron, GetVertices(), ElectronTools::EElIdType(fIdType), fWhichVertex))
    return false;

  fCutFlow->Fill(cDZ);

  //apply id cut
  if (!PassIdCut(electron))
    return false;

  fCutFlow->Fill(cId);

  //apply Isolation Cut
  if (!PassIsolationCut(electron))
    return false;

  fCutFlow->Fill(cIsolation);

  return true;
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::ElectronIdMod::PassLikelihoodId(Electron const& ele)
{
  Double_t LikelihoodValue = ElectronTools::Likelihood(fLH, &ele);

  double likCut = fIdLikelihoodCut;
  if (likCut > -900) {
    if (ele.Pt() > 20) {
      if (ele.SCluster()->AbsEta() < gkEleEBEtaMax) {
        if (ele.NumberOfClusters() == 1)
          likCut = 3.5;
        else
          likCut = 4.0;
      }
      else  {
        if (ele.NumberOfClusters() == 1)
          likCut = 4.0;
        else
          likCut = 4.0;
      }
    }
    else {
      if (ele.SCluster()->AbsEta() < gkEleEBEtaMax) {
        if (ele.NumberOfClusters() == 1)
          likCut =  4.0;
        else
          likCut =  4.5;
      }
      else {
        if (ele.NumberOfClusters() == 1)
          likCut =  4.0;
        else
          likCut =  4.0;
      }
    }
  }

  return LikelihoodValue > likCut;
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::ElectronIdMod::PassIdCut(Electron const& ele)
{
  switch (fIdType) {
  case ElectronTools::kTight:
    return ele.PassTightID();

  case ElectronTools::kLoose:
    return ele.PassLooseID();

  case ElectronTools::kLikelihood:
    return ElectronTools::PassCustomID(&ele, ElectronTools::kVBTFWorkingPointFakeableId) && PassLikelihoodId(ele);

  case ElectronTools::kNoId:
    return true;

  case ElectronTools::kCustomIdLoose:
  case ElectronTools::kCustomIdTight:
  case ElectronTools::kVBTFWorkingPointFakeableId:
  case ElectronTools::kVBTFWorkingPoint95Id:
  case ElectronTools::kVBTFWorkingPoint90Id:
  case ElectronTools::kVBTFWorkingPoint85Id:
  case ElectronTools::kVBTFWorkingPoint80Id:
  case ElectronTools::kVBTFWorkingPointLowPtId:
  case ElectronTools::kVBTFWorkingPoint70Id:
    return ElectronTools::PassCustomID(&ele, ElectronTools::EElIdType(fIdType));

  case ElectronTools::kPhys14Veto:
  case ElectronTools::kPhys14Loose:
  case ElectronTools::kPhys14Medium:
  case ElectronTools::kPhys14Tight:
    return ElectronTools::PassID(&ele, ElectronTools::EElIdType(fIdType));

  default:
    return false;
  }
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::ElectronIdMod::PassIsolationCut(Electron const& ele)
{
  switch (fIsoType) {
  case ElectronTools::kPFIso:
  case ElectronTools::kPFIsoNoL:
    if (fIsoType == ElectronTools::kPFIso) {
      return ElectronTools::PassPFIso(&ele, ElectronTools::EElIsoType(fIsoType),
                                      GetPFCandidates(),
                                      GetVertices()->At(0));
    }
    else {
      return ElectronTools::PassPFIso(&ele, ElectronTools::EElIsoType(fIsoType),
                                      GetPFCandidates(),
                                      GetVertices()->At(0),
                                      GetNonIsolatedMuons(),
                                      GetNonIsolatedElectrons());
    }

  case ElectronTools::kVBTFWorkingPoint95IndividualIso:
  case ElectronTools::kVBTFWorkingPoint90IndividualIso:
  case ElectronTools::kVBTFWorkingPoint85IndividualIso:
  case ElectronTools::kVBTFWorkingPoint70IndividualIso:
  case ElectronTools::kVBTFWorkingPoint95CombinedIso:
  case ElectronTools::kVBTFWorkingPoint90CombinedIso:
  case ElectronTools::kVBTFWorkingPoint85CombinedIso:
  case ElectronTools::kVBTFWorkingPoint70CombinedIso:
    return ElectronTools::PassCustomIso(&ele, ElectronTools::EElIsoType(fIsoType));

  case ElectronTools::kPhys14VetoIso:
  case ElectronTools::kPhys14LooseIso:
  case ElectronTools::kPhys14MediumIso:
  case ElectronTools::kPhys14TightIso:
    return ElectronTools::PassIsoRhoCorr(&ele, ElectronTools::EElIsoType(fIsoType),
                                         GetPileupEnergyDensity()->At(0)->Rho(fRhoAlgo));

  case ElectronTools::kNoIso:
    return true;

  default:
    return false;
  }
}
