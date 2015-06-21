#include "MitPhysics/Mods/interface/ElectronIDMod.h"
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

ClassImp(mithep::ElectronIDMod)

//--------------------------------------------------------------------------------------------------
ElectronIDMod::ElectronIDMod(const char *name, const char *title) :
  IDMod(name, title)
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
void ElectronIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the electron collection branch.

  // If we use MVA ID, need to load MVA weights
  switch (fIDType) {
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

  if (!PublishOutput())
    SendError(kAbortAnalysis, "SlaveBegin", "Failed to publish output objects.");
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMod::Process()
{
  // Process entries of the tree.
  auto* electrons = GetObject<mithep::ElectronCol>(fInputName);
  if (!electrons)
    SendError(kAbortAnalysis, "Process", "Electrons not found");

  TObject const* aux[nAuxInputs] = {};

  //get trigger object collection if trigger matching is enabled
  if (fApplyTriggerMatching)
    aux[kTrigObjects] = GetHLTObjects(fAuxInputNames[kTrigObjects]);

  if (fApplyConvFilterType1)
    aux[kConversions] = GetObject<mithep::DecayParticleCol>(fAuxInputNames[kConversions], true);

  if (fApplyConvFilterType1 || (fApplyD0Cut && fWhichVertex < -1))
    aux[kBeamSpot] = GetObject<mithep::BeamSpotCol>(fAuxInputNames[kBeamSpot], true);

  if ((fApplyD0Cut && fWhichVertex >= -1) || fApplyDZCut)
    aux[kVertices] = GetObject<mithep::VertexCol>(fAuxInputNames[kVertices], true);

  auto* goodElectrons = static_cast<mithep::ElectronOArr*>(fOutput);
  goodElectrons->Clear();

  for (UInt_t iE = 0; iE != electrons->GetEntries(); ++iE) {
    Electron const& electron = *electrons->At(iE);

    if (electron.SCluster() == 0)
      continue;

    if (fApplyEcalSeeded && !electron.IsEcalDriven())
      continue;

    if (electron.Pt() < fPtMin)
      continue;

    if (electron.AbsEta() > fEtaMax)
      continue;

    if (electron.SCluster()->Et() < fElectronEtMin)
      continue;

    double scAbsEta = electron.SCluster()->AbsEta();

    if (fApplyEcalFiducial &&
        ((scAbsEta > gkEleEBEtaMax && scAbsEta < gkEleEEEtaMin) || scAbsEta > gkEleEEEtaMax))
      continue;

    //apply ECAL spike removal
    if (fApplySpikeRemoval && !ElectronTools::PassSpikeRemovalFilter(&electron))
      continue;

    // apply charge filter
    if (fChargeFilter && !ElectronTools::PassChargeFilter(&electron))
      continue;

    //apply trigger matching
    if (fApplyTriggerMatching &&
        !ElectronTools::PassTriggerMatching(&electron, static_cast<mithep::TriggerObjectCol const*>(aux[kTrigObjects])))
      continue;

    // apply conversion filters
    if (fApplyConvFilterType1 &&
        !ElectronTools::PassConversionFilter(&electron,
                                             static_cast<mithep::DecayParticleCol const*>(aux[kConversions]),
                                             static_cast<mithep::BeamSpotCol const*>(aux[kBeamSpot])->At(0),
                                             0, 1e-6, 2.0, true, false))
      continue;

    if (fApplyConvFilterType2 && TMath::Abs(electron.ConvPartnerDCotTheta()) < 0.02 && TMath::Abs(electron.ConvPartnerDist()) < 0.02)
      continue;

    if (fApplyNExpectedHitsInnerCut && !ElectronTools::PassNExpectedHits(&electron, ElectronTools::EElIdType(fIDType), fInvertNExpectedHitsInnerCut))
      continue;

    // apply d0 cut
    if (fApplyD0Cut) {
      if (fWhichVertex >= -1) {
        if (!ElectronTools::PassD0Cut(&electron, static_cast<mithep::VertexCol const*>(aux[kVertices]), ElectronTools::EElIdType(fIDType), fWhichVertex))
          continue;
      }
      else if (!ElectronTools::PassD0Cut(&electron, static_cast<mithep::BeamSpotCol const*>(aux[kBeamSpot]), ElectronTools::EElIdType(fIDType)))
        continue;
    }

    // apply dz cut
    if (fApplyDZCut && !ElectronTools::PassDZCut(&electron, static_cast<mithep::VertexCol const*>(aux[kVertices]), ElectronTools::EElIdType(fIDType), fWhichVertex))
      continue;

    //apply id cut
    if(!PassIDCut(electron, aux));
      continue;

    //apply Isolation Cut
    if (!PassIsolationCut(electron, aux))
      continue;
  }

  // sort according to pt
  goodElectrons->Sort();
}

//--------------------------------------------------------------------------------------------------
void
mithep::ElectronIDMod::SlaveTerminate()
{
  RetractOutput();
}

//--------------------------------------------------------------------------------------------------
Bool_t
ElectronIDMod::PassLikelihoodID(Electron const& ele)
{
  Double_t LikelihoodValue = ElectronTools::Likelihood(fLH, &ele);

  double likCut = fIDLikelihoodCut;
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
mithep::ElectronIDMod::PassIDCut(Electron const& ele, TObject const**)
{
  switch (fIDType) {
    case ElectronTools::kTight:
      return ele.PassTightID();

    case ElectronTools::kLoose:
      return ele.PassLooseID();

    case ElectronTools::kLikelihood:
      return ElectronTools::PassCustomID(&ele, ElectronTools::kVBTFWorkingPointFakeableId) && PassLikelihoodID(ele);

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
      return ElectronTools::PassCustomID(&ele, ElectronTools::EElIdType(fIDType));

    case ElectronTools::kPhys14Veto:
      return ElectronTools::PassCustomID(&ele, ElectronTools::kPhys14Veto);

    default:
      return false;
  }
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::ElectronIDMod::PassIsolationCut(Electron const& ele, TObject const** aux)
{
  switch (fIsoType) {
  case ElectronTools::kPFIso:
  case ElectronTools::kPFIsoNoL:
    if (!aux[kPFCandidates])
      aux[kPFCandidates] = GetObject<mithep::PFCandidateCol>(fAuxInputNames[kPFCandidates], true);
    if (!aux[kVertices])
      aux[kVertices] = GetObject<mithep::VertexCol>(fAuxInputNames[kVertices], true);

    if (fIsoType == ElectronTools::kPFIso) {
      return ElectronTools::PassPFIso(&ele, ElectronTools::EElIsoType(fIsoType),
                                      static_cast<mithep::PFCandidateCol const*>(aux[kPFCandidates]),
                                      static_cast<mithep::VertexCol const*>(aux[kVertices])->At(0));
    }
    else {
      if (!aux[kNonIsolatedMuons])
        aux[kNonIsolatedMuons] = GetObject<mithep::MuonCol>(fAuxInputNames[kNonIsolatedMuons], true);
      if (!aux[kNonIsolatedElectrons])
        aux[kNonIsolatedElectrons] = GetObject<mithep::ElectronCol>(fAuxInputNames[kNonIsolatedMuons], true);
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
    if (!aux[kPileupEnergyDensity])
      aux[kPileupEnergyDensity] = GetObject<mithep::PileupEnergyDensityCol>(fAuxInputNames[kPileupEnergyDensity], true);

    return ElectronTools::PassIsoRhoCorr(&ele, ElectronTools::EElIsoType(fIsoType),
                                         static_cast<mithep::PileupEnergyDensityCol const*>(aux[kPileupEnergyDensity])->At(0)->Rho(fRhoAlgo));

  case ElectronTools::kNoIso:
    return true;

  default:
    return false;
  }
}
