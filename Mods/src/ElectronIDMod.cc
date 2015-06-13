#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
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
  BaseMod(name,title),
  fGoodElectronsName(ModNames::gkGoodElectronsName),
  fElectronBranchName(Names::gkElectronBrn),
  fConversionBranchName(Names::gkMvfConversionBrn),
  fVertexName(ModNames::gkGoodVertexesName),
  fBeamSpotName(Names::gkBeamSpotBrn),
  fTrackName(Names::gkTrackBrn),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPileupEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPFNoPileUpName("PFNoPileUp"),
  fTrigObjectsName("HLTModTrigObjs"),
  fNonIsolatedMuonsName("random"),
  fNonIsolatedElectronsName("random"),

  fElectrons(0),
  fConversions(0),
  fVertices(0),
  fBeamSpot(0),
  fPFCandidates(0),
  fPileupEnergyDensity(0),
  fPFNoPileUpCands(0),
  fNonIsolatedMuons(0),
  fNonIsolatedElectrons(0),

  fElIdType(ElectronTools::kCustomIdTight),
  fElIsoType(ElectronTools::kPFIso),
  fApplyConvFilterType1(true),
  fApplyConvFilterType2(false),
  fApplyNExpectedHitsInnerCut(true),
  fInvertNExpectedHitsInnerCut(false),
  fCombinedIdCut(false),
  fApplySpikeRemoval(true),
  fApplyD0Cut(true),
  fApplyDZCut(true),
  fChargeFilter(true),
  fApplyTriggerMatching(false),
  fApplyEcalSeeded(false),
  fApplyEcalFiducial(false),
  fElectronsFromBranch(true),
  fRhoAlgo(mithep::PileupEnergyDensity::kHighEta),
  fWhichVertex(-1),
  fElectronPtMin(10),
  fElectronEtMin(0.0),
  fElectronEtaMax(2.5),
  fIDLikelihoodCut(-999.0),
  fIntRadius(0.0),

  fLH(0),
  fElectronIDMVA(0)
{
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we just request the electron collection branch.

  if (fCombinedIdCut) {
    fElIdType = ElectronTools::kNoId;
    fElIsoType = ElectronTools::kNoIso;
    fApplyConvFilterType1 = false;
    fApplyConvFilterType2 = false;
    fApplyD0Cut           = false;
    fApplyDZCut           = false;
  }

  ReqEventObject(fElectronBranchName, fElectrons, fElectronsFromBranch);

  if (fElIdType == ElectronTools::kHggLeptonTagId2012HCP)
    ReqEventObject(fVertexName, fVertices, true);
  else if (fApplyD0Cut || fApplyDZCut || fCombinedIdCut ||
           fElIdType == ElectronTools::kMVAID_BDTG_NoIPInfo ||
           fElIdType == ElectronTools::kMVAID_BDTG_WithIPInfo ||
           fElIdType == ElectronTools::kMVAID_BDTG_IDIsoCombined ||
           fElIdType == ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0 ||
           fElIdType == ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4 ||
           fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012HCP ||
           fElIsoType == ElectronTools::kPFIso ||
           fElIsoType == ElectronTools::kPFIsoNoL ||
           fElIsoType == ElectronTools::kPFIso_HWW2012TrigV0 ||
           fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012 ||
           fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012HCP)
    ReqEventObject(fVertexName, fVertices);

  if (fApplyConvFilterType1 || (fApplyD0Cut && fWhichVertex < -1))
    ReqEventObject(fBeamSpotName, fBeamSpot, true);

  if (fApplyConvFilterType1 || fCombinedIdCut)
    ReqEventObject(fConversionBranchName, fConversions, true);

  if (fElIdType == ElectronTools::kHggLeptonTagId2012HCP)
    ReqEventObject(fPFCandidatesName, fPFCandidates, true);

  if (fElIdType == ElectronTools::kMVAID_BDTG_IDIsoCombined ||
      fElIdType == ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0 ||
      fElIdType == ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4 ||
      fElIsoType == ElectronTools::kTrackJuraSliding ||
      fElIsoType == ElectronTools::kCombinedRelativeConeAreaCorrected ||
      fElIsoType == ElectronTools::kPFIso_HWW2012TrigV0 ||
      fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012 ||
      fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012HCP ||
      fElIsoType == ElectronTools::kPhys14VetoIso) {
    ReqEventObject(fPileupEnergyDensityName, fPileupEnergyDensity, true);
  }

  // Validate options

  // If we use MVA ID, need to load MVA weights
  switch (fElIdType) {
  case ElectronTools::kLikelihood:
    if (!fLH)
      SendError(kAbortAnalysis, "SlaveBegin", "Electron likelihood is not initialized.");
    break;

  case ElectronTools::kMVAID_BDTG_NoIPInfo:
  case ElectronTools::kMVAID_BDTG_WithIPInfo:
  case ElectronTools::kMVAID_BDTG_IDIsoCombined:
  case ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0:
  case ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4:
  case ElectronTools::kHggLeptonTagId2012HCP:
    if (!fElectronIDMVA)
      SendError(kAbortAnalysis, "SlaveBegin", "Electron ID MVA is not initialized.");
    break;

  default:
    break;
  }

  switch (fElIsoType) {
  case ElectronTools::kCustomIso:
    SendError(kAbortAnalysis, "SlaveBegin", "Custom electron isolation is not yet implemented.");
    break;

  case ElectronTools::kMVAIso_BDTG_IDIsoCombinedHWW2012TrigV4:
  case ElectronTools::kMVAIso_BDTG_IDIsoCombined:
    if (!fElectronIDMVA)
      SendError(kAbortAnalysis, "SlaveBegin", "Electron ID MVA is not initialized.");
    break;

  default:
    break;
  }

  SetCutValues();
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMod::Process()
{
  // Process entries of the tree.
  LoadEventObject(fElectronBranchName, fElectrons);
  if (!fElectrons)
    SendError(kAbortAnalysis, "Process", "Electrons not found");

  if (fApplyD0Cut || fApplyDZCut || fCombinedIdCut ||
      fElIdType == ElectronTools::kMVAID_BDTG_NoIPInfo ||
      fElIdType == ElectronTools::kMVAID_BDTG_WithIPInfo ||
      fElIdType == ElectronTools::kMVAID_BDTG_IDIsoCombined ||
      fElIdType == ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0 ||
      fElIdType == ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4 ||
      fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012HCP ||
      fElIsoType == ElectronTools::kPFIso ||
      fElIsoType == ElectronTools::kPFIsoNoL ||
      fElIsoType == ElectronTools::kPFIso_HWW2012TrigV0 ||
      fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012 ||
      fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012HCP) {
    LoadEventObject(fVertexName, fVertices);
    if (!fVertices)
      SendError(kAbortAnalysis, "Process", "Vertices not found");
  }

  if (fElIsoType == ElectronTools::kPFIsoNoL) {
    fNonIsolatedMuons = GetObjThisEvt<MuonCol>(fNonIsolatedMuonsName);
    fNonIsolatedElectrons = GetObjThisEvt<ElectronCol>(fNonIsolatedElectronsName);
    if (!fNonIsolatedMuons || !fNonIsolatedElectrons)
      SendError(kAbortAnalysis, "Process", "Nonisolated leptons not found");
  }

  if (fApplyConvFilterType1 || (fApplyD0Cut && fWhichVertex < -1)) {
    LoadEventObject(fBeamSpotName, fBeamSpot);
    if (!fBeamSpot)
      SendError(kAbortAnalysis, "Process", "Beam spot not found");
  }

  if (fApplyConvFilterType1 && fCombinedIdCut) {
    LoadEventObject(fConversionBranchName, fConversions);
    if (!fConversions)
      SendError(kAbortAnalysis, "Process", "Conversions not found");
  }

  if (fElIdType == ElectronTools::kHggLeptonTagId2012HCP) {
    LoadEventObject(fPFCandidatesName, fPFCandidates);
    if (!fPFCandidates)
      SendError(kAbortAnalysis, "Process", "PF candidates not found");
  }

  if (fElIsoType == ElectronTools::kTrackJuraSliding ||
      fElIsoType == ElectronTools::kCombinedRelativeConeAreaCorrected ||
      fElIsoType == ElectronTools::kMVAIso_BDTG_IDIsoCombined ||
      fElIdType  == ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0  ||
      fElIsoType == ElectronTools::kMVAIso_BDTG_IDIsoCombinedHWW2012TrigV4 ||
      fElIsoType == ElectronTools::kPFIso_HWW2012TrigV0      ||
      fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012 ||
      fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012HCP) {
    LoadEventObject(fPileupEnergyDensityName, fPileupEnergyDensity);
    if (!fPileupEnergyDensity)
      SendError(kAbortAnalysis, "Process", "Pileup energy density not found");
  }

  // Name is hardcoded, can be changed if someone feels to do it
  if (fElIsoType == ElectronTools::kPFIso_HWW2012TrigV0 || fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012 || fElIsoType == ElectronTools::kPFIso_HggLeptonTag2012HCP)
    fPFNoPileUpCands = GetObjThisEvt<PFCandidateCol>(fPFNoPileUpName);

  //get trigger object collection if trigger matching is enabled
  const TriggerObjectCol *trigObjs = 0;
  if (fApplyTriggerMatching)
    trigObjs = GetHLTObjects(fTrigObjectsName);

  ElectronOArr *GoodElectrons = new ElectronOArr;
  GoodElectrons->SetName(fGoodElectronsName);

  int NPass = 0;

  std::vector<std::pair<double, int> > elemvaidxs;

  for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *e = fElectrons->At(i);

    if (e->SCluster() == 0)
      continue;

    if (e->Pt() < fElectronPtMin)
      continue;

    if (e->SCluster()->Et() < fElectronEtMin)
      continue;

    if (e->AbsEta() > fElectronEtaMax)
      continue;

    if (fApplyEcalFiducial && ((e->SCluster()->AbsEta() > gkEleEBEtaMax && e->SCluster()->AbsEta() < gkEleEEEtaMin) || e->SCluster()->AbsEta() > gkEleEEEtaMax))
      continue;

    if (fApplyEcalSeeded && !e->IsEcalDriven())
      continue;

    //apply trigger matching
    if (fApplyTriggerMatching && !ElectronTools::PassTriggerMatching(e,trigObjs))
      continue;

    //apply ECAL spike removal
    if (fApplySpikeRemoval && !ElectronTools::PassSpikeRemovalFilter(e))
      continue;

    //apply Isolation Cut

    if (!PassIsolationCut(e))
      continue;

    // apply conversion filters
    if (fApplyConvFilterType1 && !ElectronTools::PassConversionFilter(e, fConversions, fBeamSpot->At(0), 0, 1e-6, 2.0, true, false))
      continue;

    if (fApplyConvFilterType2 && TMath::Abs(e->ConvPartnerDCotTheta()) < 0.02 && TMath::Abs(e->ConvPartnerDist()) < 0.02)
      continue;

    if (fApplyNExpectedHitsInnerCut && !ElectronTools::PassNExpectedHits(e, fElIdType, fInvertNExpectedHitsInnerCut))
      continue;

    // apply d0 cut
    if (fApplyD0Cut) {
      if (fWhichVertex >= -1) {
        if (!ElectronTools::PassD0Cut(e, fVertices, fElIdType, fWhichVertex))
          continue;
      }
      else if (!ElectronTools::PassD0Cut(e, fBeamSpot, fElIdType))
        continue;
    }

    // apply dz cut
    if (fApplyDZCut && !ElectronTools::PassDZCut(e, fVertices, fElIdType, fWhichVertex))
      continue;

    //apply id cut
    if(!PassIDCut(e))
      continue;

    if (fElIdType == ElectronTools::kHggLeptonTagId2012HCP) {
      double MVAValue = EvaluateMVAID(e);
      //fill temp vector of indexes and mva values for sorting
      elemvaidxs.push_back(std::make_pair(MVAValue, i));
    }

    // apply charge filter
    if (fChargeFilter && !ElectronTools::PassChargeFilter(e))
      continue;

    // apply full combined id, using Tight cuts
    if (fCombinedIdCut && ElectronTools::PassTightId(e, fVertices, fConversions, 2) != 15)
      continue;

    // add good electron
    if (fElIdType != ElectronTools::kHggLeptonTagId2012HCP) {
      // make sure to mark the selected electron (to be able to pickup by skimmer)
      e->Mark();
      // add it to the good electron collection (re-ordering follows below)
      GoodElectrons->Add(e);
    }

    NPass = NPass + 1;
  }

  // sort according to pt
  GoodElectrons->Sort();

  if (fElIdType == ElectronTools::kHggLeptonTagId2012HCP) {
    // sort by mva value (in descending order)
    std::sort(elemvaidxs.begin(), elemvaidxs.end(), std::greater<std::pair<double, int> >());

    // fill final list of electrons
    for (auto&& idx : elemvaidxs)
      GoodElectrons->Add(fElectrons->At(idx.second));
  }

  // add to event for other modules to use
  AddObjThisEvt(GoodElectrons);
}

//--------------------------------------------------------------------------------------------------
void
mithep::ElectronIDMod::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void
mithep::ElectronIDMod::SetCutValues()
{
  switch (fElIdType) {
  case ElectronTools::kPhys14Veto:
    fApplyD0Cut = true;
    fApplyDZCut = true;
    fInvertNExpectedHitsInnerCut = false;
    break;
  default:
    break;
  }
}

//--------------------------------------------------------------------------------------------------
Bool_t
ElectronIDMod::PassLikelihoodID(const Electron *ele) const
{
  Double_t LikelihoodValue = ElectronTools::Likelihood(fLH, ele);

  double likCut = fIDLikelihoodCut;
  if (likCut > -900) {
    if (ele->Pt() > 20) {
      if (ele->SCluster()->AbsEta() < gkEleEBEtaMax) {
        if (ele->NumberOfClusters() == 1)
          likCut = 3.5;
        else
          likCut = 4.0;
      }
      else  {
        if (ele->NumberOfClusters() == 1)
          likCut = 4.0;
        else
          likCut = 4.0;
      }
    }
    else {
      if (ele->SCluster()->AbsEta() < gkEleEBEtaMax) {
        if (ele->NumberOfClusters() == 1)
          likCut =  4.0;
        else
          likCut =  4.5;
      }
      else {
        if (ele->NumberOfClusters() == 1)
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
ElectronIDMod::PassMVAID(const Electron *el) const
{
  Double_t MVAValue = EvaluateMVAID(el);

  Double_t eta = el->SCluster()->AbsEta();
  Int_t MVABin = -1;
  if (fElIdType == ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0 ||
      fElIdType == ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4 ||
      fElIdType == ElectronTools::kHggLeptonTagId2012HCP) {
    if (el->Pt() < 20) {
      if (eta <  0.800)
        MVABin = 0;
      else if (eta < gkEleEBEtaMax)
        MVABin = 1;
      else
        MVABin = 2;
    }
    else {
      if (eta <  0.800)
        MVABin = 3;
      else if (eta < gkEleEBEtaMax)
        MVABin = 4;
      else
        MVABin = 5;
    }
  }
  else {
    if (el->Pt() < 20) {
      if(eta <  1.000)
        MVABin = 0;
      else if (eta < gkEleEBEtaMax)
        MVABin = 1;
      else
        MVABin = 2;
    }
    else {
      if (eta <  1.000)
        MVABin = 3;
      else if(eta < gkEleEBEtaMax)
        MVABin = 4;
      else
        MVABin = 5;
    }
  }

  std::vector<double> mvaCuts;
  switch (fElIdType) {
  case ElectronTools::kMVAID_BDTG_NoIPInfo:
    mvaCuts = {0.133,  0.465,  0.518,  0.942,  0.947,  0.878};
    break;
  case ElectronTools::kMVAID_BDTG_WithIPInfo:
    mvaCuts = {0.139,  0.525,  0.543,  0.947,  0.950,  0.884};
    break;
  case ElectronTools::kMVAID_BDTG_IDIsoCombined:
    mvaCuts = {0.4202, 0.6206, 0.6190, 0.9590, 0.9586, 0.9278};
    break;
  case ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0:
    mvaCuts = {0.000,  0.100,  0.620,  0.940,  0.850,  0.920};
    break;
  case ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4:
    mvaCuts = {0.123,  0.219,  0.509,  0.935,  0.889,  0.871};
    break;
  case ElectronTools::kHggLeptonTagId2012HCP:
    mvaCuts = {0.9,    0.9,    0.9,    0.9,    0.9,    0.9};
    break;
  default:
    return false;
  }

  return MVAValue > mvaCuts[MVABin];
}

Double_t
ElectronIDMod::EvaluateMVAID(const Electron *el) const
{
  Vertex const* vertex = fVertices->At(0);

  switch (fElIdType) {
  case ElectronTools::kMVAID_BDTG_IDIsoCombined:
    return fElectronIDMVA->MVAValue(el, vertex, fPFCandidates, fPileupEnergyDensity, fIntRadius);

  case ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0:
    return fElectronIDMVA->MVAValue(el, vertex, fPFCandidates, fPileupEnergyDensity, ElectronTools::kEleEANoCorr);

  case ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4:
    return fElectronIDMVA->MVAValue(el, vertex, fVertices, fPFCandidates, fPileupEnergyDensity, ElectronTools::kEleEANoCorr);

  case ElectronTools::kHggLeptonTagId2012HCP:
  case ElectronTools::kMVAID_BDTG_NoIPInfo:
  case ElectronTools::kMVAID_BDTG_WithIPInfo:
    return fElectronIDMVA->MVAValue(el, vertex);

  default:
    throw std::runtime_error("EvaluateMVAID called for non-MVA ID");
  }
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::ElectronIDMod::PassIDCut(Electron const* ele) const
{
  switch (fElIdType) {
    case ElectronTools::kTight:
      return ele->PassTightID();

    case ElectronTools::kLoose:
      return ele->PassLooseID();

    case ElectronTools::kLikelihood:
      return ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointFakeableId) && PassLikelihoodID(ele);

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
      return ElectronTools::PassCustomID(ele, fElIdType);

    case ElectronTools::kHggLeptonTagId:
      return ElectronTools::PassHggLeptonTagID(ele);

    case ElectronTools::kHggLeptonTagId2012:
      return ElectronTools::PassHggLeptonTagID2012(ele);

    case ElectronTools::kHggLeptonTagId2012HCP:
      return PassMVAID(ele);

    case ElectronTools::kMVAID_BDTG_NoIPInfo:
    case ElectronTools::kMVAID_BDTG_WithIPInfo:
    case ElectronTools::kMVAID_BDTG_IDIsoCombined:
    case ElectronTools::kMVAID_BDTG_IDHWW2012TrigV0:
    case ElectronTools::kMVAID_BDTG_IDIsoCombinedHWW2012TrigV4:
      return ElectronTools::PassCustomID(ele, ElectronTools::kVBTFWorkingPointFakeableId) && PassMVAID(ele);

    case ElectronTools::kPhys14Veto:
      return ElectronTools::PassCustomID(ele, ElectronTools::kPhys14Veto);

    default:
      return false;
  }
}

//--------------------------------------------------------------------------------------------------
Bool_t
mithep::ElectronIDMod::PassIsolationCut(Electron const* ele) const
{
  Double_t rho = fPileupEnergyDensity->At(0)->Rho(fRhoAlgo);
  Vertex const* vertex = fVertices->At(0);

  switch (fElIsoType) {
  case ElectronTools::kPFIso:
    return ElectronTools::PassPFIso(ele, fElIsoType, fPFCandidates, vertex);

  case ElectronTools::kPFIsoNoL:
    return ElectronTools::PassPFIso(ele, fElIsoType, fPFCandidates, vertex, fNonIsolatedMuons, fNonIsolatedElectrons);

  case ElectronTools::kVBTFWorkingPoint95IndividualIso:
  case ElectronTools::kVBTFWorkingPoint90IndividualIso:
  case ElectronTools::kVBTFWorkingPoint85IndividualIso:
  case ElectronTools::kVBTFWorkingPoint70IndividualIso:
  case ElectronTools::kVBTFWorkingPoint95CombinedIso:
  case ElectronTools::kVBTFWorkingPoint90CombinedIso:
  case ElectronTools::kVBTFWorkingPoint85CombinedIso:
  case ElectronTools::kVBTFWorkingPoint70CombinedIso:
    return ElectronTools::PassCustomIso(ele, fElIsoType);

  case ElectronTools::kPFIso_HWW2012TrigV0:
    return ElectronTools::PassIsoRhoCorr(ele, fElIsoType, rho, fPFNoPileUpCands, vertex);

  case ElectronTools::kPFIso_HggLeptonTag2012HCP:
    {
      Double_t distVtx = 999.0;
      for (UInt_t nv = 0; nv != fVertices->GetEntries(); ++nv) {
        double dz = TMath::Abs(ele->GsfTrk()->DzCorrected(*fVertices->At(nv)));
        if (dz > distVtx)
          continue;
        distVtx = dz;
        vertex = fVertices->At(nv);
      }
    }
    //fallthrough
  case ElectronTools::kPFIso_HggLeptonTag2012:
    return ElectronTools::PassIsoRhoCorr(ele, fElIsoType, rho, fPFNoPileUpCands, vertex);

  case ElectronTools::kMVAIso_BDTG_IDIsoCombinedHWW2012TrigV4:
    return ElectronTools::PassIso(ele, fElIsoType);

  case ElectronTools::kPhys14VetoIso:
    return ElectronTools::PassIsoRhoCorr(ele, fElIsoType, rho);

  case ElectronTools::kNoIso:
    return true;

  default:
    return false;
  }
}

void
mithep::ElectronIDMod::SetRhoType(RhoUtilities::RhoType type)
{
  SendError(kWarning, "SetRhoType", "ElectronIDMod::SetRhoType is deprecaed. Use SetRhoAlgo with mithep::PileupEnergyDensity::Algo instead.");

  switch (type) {
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
  }
}

void
mithep::ElectronIDMod::SetNExpectedHitsInnerCut(Int_t)
{
  SendError(kWarning, "SetNExpectedHitsInnerCut", "ElectronIDMod is not responsible for setting the cut value on the expected inner hits. Implement an appropriate cut for the given electron id in ElectronTools::PassNExpectedHits.");
}

void
mithep::ElectronIDMod::SetApplyD0Cut(Bool_t b)
{
  if (b)
    SendError(kWarning, "SetApplyD0Cut", "D0 cut value is not stored in the ElectronIDMod any more. Make sure that the appropriate cuts are impemented in the ElectronTools::PassD0Cut.");

  fApplyD0Cut = b;
}

void
mithep::ElectronIDMod::SetApplyDZCut(Bool_t b)
{
  if (b)
    SendError(kWarning, "SetApplyDZCut", "DZ cut value is not stored in the ElectronIDMod any more. Make sure that the appropriate cuts are impemented in the ElectronTools::PassDZCut.");

  fApplyDZCut = b;
}
