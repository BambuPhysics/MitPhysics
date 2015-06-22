#include "MitPhysics/Mods/interface/ElectronIdMod.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Init/interface/Constants.h"

#include <algorithm>
#include <limits>

using namespace mithep;

ClassImp(mithep::ElectronIdMod)

template<class T>
void
ElectronIdMod::GetAuxInput(ElectronIdMod::AuxInput inputCol, TObject const** aux)
{
  aux[inputCol] = GetObject<T>(fAuxInputNames[inputCol], true);
  if (!aux[inputCol])
    SendError(kAbortAnalysis, "GetAuxInput", "Could not retrieve auxiliary input " + fAuxInputNames[inputCol]);
}

//--------------------------------------------------------------------------------------------------
ElectronIdMod::ElectronIdMod(const char *name, const char *title) :
  IdMod(name, title)
{
  fOutput = new ElectronOArr(32, TString(name) + "Output");

  fAuxInputNames[kConversions] = Names::gkMvfConversionBrn;
  fAuxInputNames[kVertices] = ModNames::gkGoodVertexesName;
  fAuxInputNames[kBeamSpot] = Names::gkBeamSpotBrn;
  fAuxInputNames[kPFCandidates] = Names::gkPFCandidatesBrn;
  fAuxInputNames[kPileupEnergyDensity] = Names::gkPileupEnergyDensityBrn;
  fAuxInputNames[kNonIsolatedMuons] = "random";
  fAuxInputNames[kNonIsolatedElectrons] = "random";
}

//--------------------------------------------------------------------------------------------------
void
ElectronIdMod::IdBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the electron collection branch.

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
void
ElectronIdMod::Process()
{
  // Process entries of the tree.
  auto* electrons = GetObject<mithep::ElectronCol>(fInputName, true);
  if (!electrons)
    SendError(kAbortAnalysis, "Process", "Electrons not found");

  TObject const* aux[nAuxInputs] = {};

  //get trigger object collection if trigger matching is enabled
  if (fApplyTriggerMatching) {
    aux[kTrigObjects] = GetHLTObjects(fAuxInputNames[kTrigObjects]);
    if (!aux[kTrigObjects])
      SendError(kAbortAnalysis, "Process", "TrigObjects not found");
  }

  if (fApplyConvFilterType1)
    GetAuxInput<mithep::DecayParticleCol>(kConversions, aux);

  if (fApplyConvFilterType1 || (fApplyD0Cut && fWhichVertex < -1))
    GetAuxInput<mithep::BeamSpotCol>(kBeamSpot, aux);

  if ((fApplyD0Cut && fWhichVertex >= -1) || fApplyDZCut)
    GetAuxInput<mithep::VertexCol>(kVertices, aux);

  mithep::ElectronOArr* goodElectrons = 0;
  if (fIsFilterMode) {
    goodElectrons = static_cast<mithep::ElectronOArr*>(fOutput);
    goodElectrons->Reset();
  }
  else {
    fFlags.Resize(electrons->GetEntries());
    for (unsigned iE = 0; iE != electrons->GetEntries(); ++iE)
      fFlags.At(iE) = false;
  }

  for (UInt_t iE = 0; iE != electrons->GetEntries(); ++iE) {
    Electron const& electron = *electrons->At(iE);

    fCutFlow->Fill(cAll);

    if (electron.SCluster() == 0)
      continue;

    if (fApplyEcalSeeded && !electron.IsEcalDriven())
      continue;

    fCutFlow->Fill(cIsEcalDriven);

    if (electron.Pt() < fPtMin)
      continue;

    fCutFlow->Fill(cPt);

    if (electron.AbsEta() > fEtaMax)
      continue;

    fCutFlow->Fill(cEta);

    if (electron.SCluster()->Et() < fElectronEtMin)
      continue;

    fCutFlow->Fill(cEt);

    double scAbsEta = electron.SCluster()->AbsEta();

    if (fApplyEcalFiducial &&
        ((scAbsEta > gkEleEBEtaMax && scAbsEta < gkEleEEEtaMin) || scAbsEta > gkEleEEEtaMax))
      continue;

    fCutFlow->Fill(cFiducial);

    //apply ECAL spike removal
    if (fApplySpikeRemoval && !ElectronTools::PassSpikeRemovalFilter(&electron))
      continue;

    fCutFlow->Fill(cSpikeRemoval);

    // apply charge filter
    if (fChargeFilter && !ElectronTools::PassChargeFilter(&electron))
      continue;

    fCutFlow->Fill(cChargeFilter);

    //apply trigger matching
    if (fApplyTriggerMatching &&
        !ElectronTools::PassTriggerMatching(&electron, static_cast<mithep::TriggerObjectCol const*>(aux[kTrigObjects])))
      continue;

    fCutFlow->Fill(cTriggerMatching);

    // apply conversion filters
    if (fApplyConvFilterType1 &&
        !ElectronTools::PassConversionFilter(&electron,
                                             static_cast<mithep::DecayParticleCol const*>(aux[kConversions]),
                                             static_cast<mithep::BeamSpotCol const*>(aux[kBeamSpot])->At(0),
                                             0, 1e-6, 2.0, true, false))
      continue;

    fCutFlow->Fill(cConvFilterType1);

    if (fApplyConvFilterType2 && TMath::Abs(electron.ConvPartnerDCotTheta()) < 0.02 && TMath::Abs(electron.ConvPartnerDist()) < 0.02)
      continue;

    fCutFlow->Fill(cConvFilterType2);

    if (fApplyNExpectedHitsInnerCut && !ElectronTools::PassNExpectedHits(&electron, ElectronTools::EElIdType(fIdType), fInvertNExpectedHitsInnerCut))
      continue;

    fCutFlow->Fill(cNExpectedHits);

    // apply d0 cut
    if (fApplyD0Cut) {
      if (fWhichVertex >= -1) {
        if (!ElectronTools::PassD0Cut(&electron, static_cast<mithep::VertexCol const*>(aux[kVertices]), ElectronTools::EElIdType(fIdType), fWhichVertex))
          continue;
      }
      else if (!ElectronTools::PassD0Cut(&electron, static_cast<mithep::BeamSpotCol const*>(aux[kBeamSpot]), ElectronTools::EElIdType(fIdType)))
        continue;
    }

    fCutFlow->Fill(cD0);

    // apply dz cut
    if (fApplyDZCut && !ElectronTools::PassDZCut(&electron, static_cast<mithep::VertexCol const*>(aux[kVertices]), ElectronTools::EElIdType(fIdType), fWhichVertex))
      continue;

    fCutFlow->Fill(cDZ);

    //apply id cut
    if (!PassIdCut(electron, aux))
      continue;

    fCutFlow->Fill(cId);

    //apply Isolation Cut
    if (!PassIsolationCut(electron, aux))
      continue;

    fCutFlow->Fill(cIsolation);

    if (fIsFilterMode)
      goodElectrons->Add(&electron);
    else
      fFlags.At(iE) = true;
  }

  if (fIsFilterMode) {
    // sort according to pt
    goodElectrons->Sort();
  }
}

//--------------------------------------------------------------------------------------------------
Bool_t
ElectronIdMod::PassLikelihoodId(Electron const& ele)
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
mithep::ElectronIdMod::PassIdCut(Electron const& ele, TObject const**)
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
mithep::ElectronIdMod::PassIsolationCut(Electron const& ele, TObject const** aux)
{
  switch (fIsoType) {
  case ElectronTools::kPFIso:
  case ElectronTools::kPFIsoNoL:
    if (!aux[kPFCandidates])
      GetAuxInput<mithep::PFCandidateCol>(kPFCandidates, aux);
    if (!aux[kVertices])
      GetAuxInput<mithep::VertexCol>(kVertices, aux);

    if (fIsoType == ElectronTools::kPFIso) {
      return ElectronTools::PassPFIso(&ele, ElectronTools::EElIsoType(fIsoType),
                                      static_cast<mithep::PFCandidateCol const*>(aux[kPFCandidates]),
                                      static_cast<mithep::VertexCol const*>(aux[kVertices])->At(0));
    }
    else {
      if (!aux[kNonIsolatedMuons])
        GetAuxInput<mithep::MuonCol>(kNonIsolatedMuons, aux);
      if (!aux[kNonIsolatedElectrons])
        GetAuxInput<mithep::ElectronCol>(kNonIsolatedElectrons, aux);

      return ElectronTools::PassPFIso(&ele, ElectronTools::EElIsoType(fIsoType),
                                      static_cast<mithep::PFCandidateCol const*>(aux[kPFCandidates]),
                                      static_cast<mithep::VertexCol const*>(aux[kVertices])->At(0),
                                      static_cast<mithep::MuonCol const*>(aux[kNonIsolatedMuons]),
                                      static_cast<mithep::ElectronCol const*>(aux[kNonIsolatedElectrons]));
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
    if (!aux[kPileupEnergyDensity])
      GetAuxInput<mithep::PileupEnergyDensityCol>(kPileupEnergyDensity, aux);

    return ElectronTools::PassIsoRhoCorr(&ele, ElectronTools::EElIsoType(fIsoType),
                                         static_cast<mithep::PileupEnergyDensityCol const*>(aux[kPileupEnergyDensity])->At(0)->Rho(fRhoAlgo));

  case ElectronTools::kNoIso:
    return true;

  default:
    return false;
  }
}
