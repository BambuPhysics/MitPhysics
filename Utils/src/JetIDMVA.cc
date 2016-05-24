#include "MitPhysics/Utils/interface/JetIDMVA.h"
#include "MitPhysics/Utils/interface/JetTools.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include <stdexcept>

ClassImp(mithep::JetIDMVA)

using namespace mithep;

BitMask8 JetIDMVA::fgCorrectionMask(Long64_t(0x7)); // L1+L2+L3

JetIDMVA::JetIDMVA() :
  fVarNames{   // don't really need to have this member variable, but is convenient
    "nvtx",
    "jetPt",
    "jetEta",
    "jetPhi",
    "d0",
    "dZ",
    "beta",
    "betaStar",
    "nCharged",
    "nNeutrals",
    "dRMean",
    "ptD",
    "frac01",
    "frac02",
    "frac03",
    "frac04",
    "frac05",
    "dR2Mean",
    "DRweighted",
    "rho",
    "nParticles",
    "nCh",
    "majW",
    "minW",
    "fRing0",
    "fRing1",
    "fRing2",
    "fRing3",
    "pull",
    "min(pull,0.1)",
    "jetR",
    "jetRchg",
    "nTrueInt",
    "dRMatch"
  }
{
  assert(fVarNames[nVariables - 1] != ""); // checking for missing lines
}

JetIDMVA::~JetIDMVA()
{
  for (auto* reader : fReaders)
    delete reader;
}

void
JetIDMVA::Initialize(JetIDMVA::CutType cutType, JetIDMVA::MVAType mvaType,
                     TString const& weightsConfig, TString const& cutConfig)
{
  Initialize(cutType, mvaType, std::vector<TString>(1, weightsConfig), cutConfig);
}

void
JetIDMVA::Initialize(JetIDMVA::CutType cutType, JetIDMVA::MVAType mvaType,
                     std::vector<TString> const& weightsConfig, TString const& cutConfig)
{
  // need to update this interface to allow eta-binned weights

  if (cutConfig.Length() == 0)
    throw std::runtime_error("JetIDMVA missing necessary input files");

  if (fIsInitialized)
    throw std::runtime_error("Attempting to initialize JetIDMVA twice");

  fType = mvaType;

  TString lCutId;
  switch (fType) {
  case kBaseline:
    lCutId = "Baseline";
    break;
  case k42:
    lCutId = "PuJetIdOptMVA_wp";
    break;
  case k52:
    lCutId = "full_5x_wp";
    break;
  case kCut:
    lCutId = "PuJetIdCutBased_wp";
    break;
  case k53:
    lCutId = "full_53x_wp";
    break;
  case k53CHS:
  case k53BDTCHSFullPlusRMS:
    lCutId = "full_53x_chs_wp";
    break;
  case k53MET:
    lCutId = "met_53x_wp";
    break;
  case k53METFull:
    lCutId = "metfull_53x_wp";
    break;
  case k74CHS:
    lCutId = "full_74x_chs_wp";
    break;
  default:
     lCutId = "default";
    break;
  }

  TString lCutType;
  switch (cutType) {
  case kTight:
    lCutType = "Tight";
    break;
  case kMedium:
    lCutType = "Medium";
    break;
  case kLoose:
    lCutType = "Loose";
    break;
  case kMET:
    lCutType = "MET";
    break;
  default:
    break;
  }

  if (!InitializeCuts(cutConfig, lCutId, lCutType))
    throw std::runtime_error("JetIDMVA cut values not set!");

  if (fType == kCut) {
    fIsInitialized = kTRUE;
    return;
  }

  std::vector<std::vector<unsigned>> variables;
  std::vector<unsigned> spectators;
  fEtaBinLowEdges.assign(1, 0.);
  std::vector<unsigned> varIndexForBin(1, 0);
    
  switch (fType) {
  case kBaseline:
    variables.push_back({kNvtx, kJetPt, kJetEta, kJetPhi, kDZ, kD0, kBeta, kBetaStar,
          kNCharged, kNNeutrals, kDRMean, kFrac01, kFrac02, kFrac03, kFrac04, kFrac05});
    break;
  case k42:
    variables.push_back({kFrac01, kFrac02, kFrac03, kFrac04, kFrac05, kNvtx, kNNeutrals,
          kBeta, kBetaStar, kDZ, kNCharged});
    spectators = {kJetPt, kJetEta};
    break;
  case k52:
    variables.push_back({kFrac01, kFrac02, kFrac03, kFrac04, kFrac05, kDR2Mean, kNvtx, kNNeutrals,
          kBeta, kBetaStar, kDZ, kNCharged});
    spectators = {kJetPt, kJetEta};
    break;
  case k53:
  case k53CHS:
    variables.push_back({kNvtx, kDZ, kBeta, kBetaStar, kNCharged, kNNeutrals, kDR2Mean, kPtD,
          kFrac01, kFrac02, kFrac03, kFrac04, kFrac05});
    spectators = {kJetPt, kJetEta, kJetPhi};
    break;
  case k53BDTCHSFullPlusRMS:
    variables.push_back({kFrac01, kFrac02, kFrac03, kFrac04, kFrac05, kDR2Mean, kNvtx, kNNeutrals,
          kBeta, kBetaStar, kDZ, kNCharged});
    spectators = {kJetPt, kJetEta};
    break;
  case k74CHS:
    variables.push_back({kDR2Mean, kRho, kNParticles, kNCharged, kAxisMajor, kAxisMinor,
          kFrac01, kFrac02, kFrac03, kFrac04, kPtD, kBeta, kBetaStar, kPull, kJetR, kJetRchg});
    variables.push_back({kDR2Mean, kRho, kNParticles, kAxisMajor, kAxisMinor,
          kFrac01, kFrac02, kFrac03, kFrac04, kPtD, kPull, kJetR});
    spectators = {kJetPt, kJetEta, kNTrueInt, kDRMatch};
    fEtaBinLowEdges.push_back(2.);
    varIndexForBin.push_back(0); // variable set for [2-2.5]
    fEtaBinLowEdges.push_back(2.5);
    varIndexForBin.push_back(0); // variable set for [2.5-3]
    fEtaBinLowEdges.push_back(3.);
    varIndexForBin.push_back(1); // variable set for [3-5]
    break;
  case k53MET:
    variables.push_back({kNvtx, kJetPt, kJetEta, kJetPhi, kDZ, kBeta, kBetaStar, kNCharged, kNNeutrals,
          kDR2Mean, kPtD, kFrac01, kFrac02, kFrac03, kFrac04, kFrac05});
    break;
  case kQGP:
    variables.push_back({kNvtx, kJetPt, kJetEta, kJetPhi, kBeta, kBetaStar, kNCharged, kNNeutrals,
          kDR2Mean, kPtD, kFrac01, kFrac02, kFrac03, kFrac04, kFrac05});
    break;
  default:
    break;
  }
  fEtaBinLowEdges.push_back(5.);

  if (weightsConfig.size() != fEtaBinLowEdges.size() - 1)
    throw std::runtime_error("JetIDMVA wrong number of weight files!");

  for (bool& used : fVariableUsed)
    used = false;

  for (unsigned iBin(0); iBin != weightsConfig.size(); ++iBin) {
    auto* reader = new TMVA::Reader("!Color:!Silent");

    for (auto iV : variables[varIndexForBin[iBin]]) {
      reader->AddVariable(fVarNames[iV], fVariables + iV);
      fVariableUsed[iV] = true;
    }
    for (auto iS : spectators) {
      reader->AddSpectator(fVarNames[iS], fVariables + iS);
      fVariableUsed[iS] = true;
    }

    reader->BookMVA(fMethodName, weightsConfig[iBin]);

    fReaders.push_back(reader);
  }

  fIsInitialized = kTRUE;
}

//--------------------------------------------------------------------------------------------------
void
JetIDMVA::Initialize(JetIDMVA::CutType iCutType,
                     TString const& iLowPtWeights,
                     TString const& iHighPtWeights,
                     JetIDMVA::MVAType iType,
                     TString const& iCutFileName)
{
  Initialize(iCutType, iType, iHighPtWeights, iCutFileName);
}

//--------------------------------------------------------------------------------------------------
Bool_t
JetIDMVA::passCut(const PFJet *iJet, const Vertex *iVertex, const VertexCol *iVertices)
{
  MVAType currentType = fType;
  fType = kCut;
  Bool_t res = pass(iJet, iVertex, iVertices);
  fType = currentType;
  return res;
}

//--------------------------------------------------------------------------------------------------
Bool_t
JetIDMVA::pass(const PFJet *iJet, const Vertex *iVertex, const VertexCol *iVertices, Double_t rho)
{
  // A PF Jet with L1+L2+L3 corrections is expected.
  if (iJet->Corrections() != fgCorrectionMask)
    throw std::runtime_error("JetIDMVA works only with L1+L2+L3 corrected jets");

  double lEta = iJet->AbsEta();
  if (lEta > 4.99)
    return false;

  if(!JetTools::passPFId(iJet, JetTools::kPFLoose))
    return false;

  double lPt = iJet->Pt(); // use corrected Pt

  int lPtId = 3;
  if (lPt <= 10.)
    lPtId = 0;
  else if (lPt <= 20.)
    lPtId = 1;
  else if (lPt <= 30.)
    lPtId = 2;

  int lEtaId = 3;
  if (lEta <= 2.5)
    lEtaId = 0;
  else if (lEta <= 2.75)
    lEtaId = 1;
  else if (lEta <= 3.)
    lEtaId = 2;

  if (fType == kCut) {
    float betaStarModified = JetTools::betaStarClassic(iJet,iVertex,iVertices)/log(iVertices ->GetEntries()-0.64);
    float dR2Mean          = JetTools::dR2Mean(iJet,-1);

    if(betaStarModified < fBetaStarCut[lPtId][lEtaId] &&
       dR2Mean          < fRMSCut[lPtId][lEtaId])
      return true;
  }
  else {
    double lMVA = MVAValue(iJet, iVertex, iVertices, rho);
    double lMVACut = fMVACut[lPtId][lEtaId];
    if (lMVA > lMVACut)
      return true;
  }

  return false;
}

//--------------------------------------------------------------------------------------------------
Double_t
JetIDMVA::MVAValue(const PFJet *iJet, const Vertex *iVertex, //Vertex here is the PV
                   const VertexCol *iVertices, Double_t rho,
                   Bool_t printDebug)
{
  if (!fIsInitialized)
    throw std::runtime_error("Error: JetIDMVA not properly initialized.");

  // A PF Jet with L1+L2+L3 corrections is expected.
  if (iJet->Corrections() != fgCorrectionMask)
    throw std::runtime_error("JetIDMVA works only with L1+L2+L3 corrected jets");

  if (!JetTools::passPFId(iJet, JetTools::kPFLoose))
    return -2.;

  auto&& covariance(JetTools::W(iJet));
  double sumPt = JetTools::sumPt(iJet);

  //set all input variables
  if (fVariableUsed[kNvtx])
    fVariables[kNvtx]  = iVertices->GetEntries();
  if (fVariableUsed[kRho])
    fVariables[kRho] = rho;
  if (fVariableUsed[kJetPt])
    fVariables[kJetPt] = iJet->Pt();
  if (fVariableUsed[kJetEta])
    fVariables[kJetEta] = iJet->Eta();
  if (fVariableUsed[kJetPhi])
    fVariables[kJetPhi]    = iJet->Phi();
  if (fVariableUsed[kD0])
    fVariables[kD0]        = JetTools::impactParameter(iJet, iVertex);
  if (fVariableUsed[kDZ])
    fVariables[kDZ]        = JetTools::impactParameter(iJet, iVertex, true);
  if (fDZCut > 0.) {
    if (fVariableUsed[kBeta])
      fVariables[kBeta]      = JetTools::Beta(iJet, iVertex, fDZCut);
    if (fVariableUsed[kBetaStar])
      fVariables[kBetaStar]  = JetTools::betaStar(iJet, iVertex, iVertices, fDZCut);
  }
  else {
    if (fVariableUsed[kBeta])
      fVariables[kBeta]      = JetTools::BetaClassic(iJet, iVertex);
    if (fVariableUsed[kBetaStar])
      fVariables[kBetaStar]  = JetTools::betaStarClassic(iJet, iVertex, iVertices);
  }
  if (fVariableUsed[kNCharged] || fVariableUsed[kNCh])
    fVariables[kNCharged] = fVariables[kNCh] = iJet->ChargedMultiplicity();
  if (fVariableUsed[kNNeutrals])
    fVariables[kNNeutrals] = iJet->NeutralMultiplicity();
  if (fVariableUsed[kPtD])
    fVariables[kPtD] = covariance.ptD;
  if (fVariableUsed[kPull] || fVariableUsed[kMinPull01]) {
    fVariables[kPull] = JetTools::pull(iJet, fReproducePullBug);
    fVariables[kMinPull01] = std::min(double(fVariables[kPull]), 0.1);
  }
  if (fVariableUsed[kDRMean])
    fVariables[kDRMean] = JetTools::dRMean(iJet);
  if (fVariableUsed[kDR2Mean] || fVariableUsed[kDRWeighted])
    fVariables[kDR2Mean] = fVariables[kDRWeighted] = JetTools::dR2Mean(iJet);
  if (fVariableUsed[kFrac01] || fVariableUsed[kFRing0])
    fVariables[kFrac01] = fVariables[kFRing0] = JetTools::frac(iJet, 0.1, 0.,  -1);
  if (fVariableUsed[kFrac02] || fVariableUsed[kFRing1])
    fVariables[kFrac02] = fVariables[kFRing1] = JetTools::frac(iJet, 0.2, 0.1);
  if (fVariableUsed[kFrac03] || fVariableUsed[kFRing2])
    fVariables[kFrac03] = fVariables[kFRing2] = JetTools::frac(iJet, 0.3, 0.2);
  if (fVariableUsed[kFrac04] || fVariableUsed[kFRing3])
    fVariables[kFrac04] = fVariables[kFRing3] = JetTools::frac(iJet, 0.4, 0.3);
  if (fVariableUsed[kFrac05])
    fVariables[kFrac05] = JetTools::frac(iJet, 0.5, 0.4);
  if (fVariableUsed[kNParticles])
    fVariables[kNParticles] = iJet->NConstituents();
  if (fVariableUsed[kAxisMajor])
    fVariables[kAxisMajor] = covariance.majW;
  if (fVariableUsed[kAxisMinor])
    fVariables[kAxisMinor] = covariance.minW;
  if (fVariableUsed[kJetR]) {
    auto* leadCand = JetTools::leadCand(iJet, false);
    if (leadCand) // has to be nonnull
      fVariables[kJetR] = leadCand->Pt() / sumPt;
  }
  if (fVariableUsed[kJetRchg]) {
    auto* leadEmCand = JetTools::leadCand(iJet, false, PFCandidate::eGamma);
    if (!leadEmCand)
      leadEmCand = JetTools::trailCand(iJet);
    fVariables[kJetRchg] = leadEmCand->Pt() / sumPt;
  }
  if (fVariableUsed[kNTrueInt]) // never used!
    fVariables[kNTrueInt] = 0.;
  if (fVariableUsed[kDRMatch])
    fVariables[kDRMatch] = JetTools::dRMin(iJet);

  double absEta = iJet->AbsEta();
  if (absEta >= fEtaBinLowEdges.back())
    absEta = fEtaBinLowEdges.back() - 0.01;

  int iBin = -1;
  while (fEtaBinLowEdges[iBin + 1] < absEta)
    ++iBin;

  double lMVA = fReaders[iBin]->EvaluateMVA(fMethodName);

  if (printDebug) {
    std::cout << "Debug Jet MVA:" << std::endl;
    for (unsigned iV = 0; iV != nVariables; ++iV) {
      if (fVariableUsed[iV])
        std::cout << fVarNames[iV] << " = " << fVariables[iV] << std::endl;
    }
    std::cout << "=== : === " << std::endl;
    std::cout << "MVA value " << lMVA << std::endl;
  }

  return lMVA;
}

//--------------------------------------------------------------------------------------------------
Double_t*
JetIDMVA::QGValue(const PFJet *iJet, const Vertex *iVertex, //Vertex here is the PV
                  const VertexCol *iVertices,
                  const PileupEnergyDensityCol *iPileupEnergyDensity,
                  Bool_t printDebug)
{
  Double_t *lId = new double[3];
  lId[0] = -2;
  lId[1] = -2;
  lId[2] = -2;
  if (!fIsInitialized) {
    std::cout << "Error: JetIDMVA not properly initialized.\n";
    return lId;
  }
  if (!JetTools::passPFId(iJet, JetTools::kPFLoose))
    return lId;

  fVariables[kJetPt]       = iJet->Pt();
  if (fVariables[kJetPt] < 20)
    return lId;

  //set all input variables
  fVariables[kNvtx]      = iVertices->GetEntries();
  fVariables[kJetEta]    = iJet->RawMom().Eta();
  fVariables[kJetPhi]    = iJet->RawMom().Phi();
  fVariables[kD0]        = JetTools::impactParameter(iJet,iVertex);
  fVariables[kDZ]        = JetTools::impactParameter(iJet,iVertex,true);
  if (fDZCut > 0.) {
    fVariables[kBeta]      = JetTools::Beta(iJet, iVertex, fDZCut);
    fVariables[kBetaStar]  = JetTools::betaStar(iJet, iVertex, iVertices, fDZCut);
  }
  else {
    fVariables[kBeta]      = JetTools::BetaClassic(iJet, iVertex);
    fVariables[kBetaStar]  = JetTools::betaStarClassic(iJet, iVertex, iVertices);
  }
  fVariables[kNCharged]  = iJet->ChargedMultiplicity();
  fVariables[kNNeutrals] = iJet->NeutralMultiplicity();
  fVariables[kPtD]       = JetTools::W(iJet,-1).ptD;
  fVariables[kDRMean]    = JetTools::dRMean(iJet,-1);
  fVariables[kDR2Mean]   = JetTools::dR2Mean(iJet,-1);
  fVariables[kFrac01]    = JetTools::frac(iJet,0.1,0.);
  fVariables[kFrac02]    = JetTools::frac(iJet,0.2,0.1);
  fVariables[kFrac03]    = JetTools::frac(iJet,0.3,0.2);
  fVariables[kFrac04]    = JetTools::frac(iJet,0.4,0.3);
  fVariables[kFrac05]    = JetTools::frac(iJet,0.5,0.4);

  double absEta = iJet->AbsEta();
  if (absEta >= 5.)
    absEta = 4.99;

  int iBin = -1;
  while (fEtaBinLowEdges[iBin + 1] < absEta)
    ++iBin;

  lId[0] = fReaders[iBin]->EvaluateMulticlass(fMethodName)[0];
  lId[1] = fReaders[iBin]->EvaluateMulticlass(fMethodName)[1];
  lId[2] = fReaders[iBin]->EvaluateMulticlass(fMethodName)[2];

  return lId;
}

Bool_t
JetIDMVA::InitializeCuts(TString const& fileName, TString const& cutId, TString const& cutType)
{
  //Load Cut Matrix
  std::ifstream cutFile(fileName);

  //flag cuts being set. for fType == kCut we need two sets of cuts
  bool cutsSet[2] = {false, fType != kCut};

  while (true) {
    std::string line;
    std::getline(cutFile, line);
    if (!cutFile.good())
      break;

    // restream the line into words
    std::stringstream ss;
    ss.str(line);
    std::string word;

    // first word: cut ID (e.g. full_53x_wp)
    ss >> word;
    if (word != cutId)
      continue;

    // second word: cut type (e.g. Tight)
    ss >> word;
    if (word.find(cutType) == std::string::npos)
      continue;

    // A line with cutId and cutType is found.
    // Now set the target cut array, depending on the block type
    typedef float (*ArraysOfFour)[4];
    ArraysOfFour cutArray = 0;
    bool* flag = 0;
    if (fType == kCut) {
      if (word.find("BetaStar") == 0) {
        cutArray = fBetaStarCut;
        flag = cutsSet;
      }
      else if (word.find("RMS") == 0) {
        cutArray = fRMSCut;
        flag = cutsSet + 1;
      }
      else
        continue;
    }
    else {
      cutArray = fMVACut;
      flag = cutsSet;
    }

    // the next four lines in the config define cut values
    unsigned iL = 0;
    for (; iL != 4; ++iL) {
      std::getline(cutFile, line);
      if (!cutFile.good()) // input format not correct
        break;

      ss.clear();
      ss.str(line);
      // each line has four numbers = four bins in eta
      unsigned iE = 0;
      for (; iE != 4; ++iE) {
        ss >> cutArray[iL][iE];
        if (iE != 3 && !ss.good())
          break;
      }
      if (iE != 4) // input format not correct
        break;
    }
    if (iL != 4) {
      // input format not correct
      throw std::runtime_error(("JetIDMVA could not parse " + fileName).Data());
    }

    // config correctly read. set flag.
    *flag = true;

    if (cutsSet[0] && cutsSet[1]) // nothing more to read
      break;
  }

  cutFile.close();

  return cutsSet[0] && cutsSet[1];
}
