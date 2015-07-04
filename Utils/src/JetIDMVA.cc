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
    "dR2Mean"
  }
{
  assert(fVarNames[nVariables - 1] != ""); // checking for missing lines
}

JetIDMVA::~JetIDMVA()
{
  delete fReader;
}

void
JetIDMVA::Initialize(JetIDMVA::CutType cutType, JetIDMVA::MVAType mvaType,
                     TString const& weightsConfig, TString const& cutConfig)
{
  if (weightsConfig.Length() == 0 || cutConfig.Length() == 0)
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
    lCutId = "full_53x_chs_wp";
    break;
  case k53MET:
    lCutId = "met_53x_wp";
    break;
  case k53METFull:
    lCutId = "metfull_53x_wp";
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

  fReader        = new TMVA::Reader("!Color:!Silent:Error" );

  std::vector<unsigned> variables;
  std::vector<unsigned> spectators;
    
  switch (fType) {
  case kBaseline:
    variables = {kNvtx, kJetPt, kJetEta, kJetPhi, kDZ, kD0, kBeta, kBetaStar,
                 kNCharged, kNNeutrals, kDRMean, kFrac01, kFrac02, kFrac03, kFrac04, kFrac05};
    break;
  case k42:
    variables = {kFrac01, kFrac02, kFrac03, kFrac04, kFrac05, kNvtx, kNNeutrals,
                 kBeta, kBetaStar, kDZ, kNCharged};
    spectators = {kJetPt, kJetEta};
    break;
  case k52:
    variables = {kFrac01, kFrac02, kFrac03, kFrac04, kFrac05, kDR2Mean, kNvtx, kNNeutrals,
                 kBeta, kBetaStar, kDZ, kNCharged};
    spectators = {kJetPt, kJetEta};
    break;
  case k53:
  case k53CHS:
    variables = {kNvtx, kDZ, kBeta, kBetaStar, kNCharged, kNNeutrals, kDR2Mean, kPtD,
                 kFrac01, kFrac02, kFrac03, kFrac04, kFrac05};
    spectators = {kJetPt, kJetEta, kJetPhi};
    break;
  case k53MET:
    variables = {kNvtx, kJetPt, kJetEta, kJetPhi, kDZ, kBeta, kBetaStar, kNCharged, kNNeutrals,
                 kDR2Mean, kPtD, kFrac01, kFrac02, kFrac03, kFrac04, kFrac05};
    break;
  case kQGP:
    variables = {kNvtx, kJetPt, kJetEta, kJetPhi, kBeta, kBetaStar, kNCharged, kNNeutrals,
                 kDR2Mean, kPtD, kFrac01, kFrac02, kFrac03, kFrac04, kFrac05};
    break;
  default:
    break;
  }

  for (auto iV : variables)
    fReader->AddVariable(fVarNames[iV], fVariables + iV);
  for (auto iS : spectators)
    fReader->AddSpectator(fVarNames[iS], fVariables + iS);

  fReader->BookMVA(fMethodName, weightsConfig);

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
JetIDMVA::pass(const PFJet *iJet, const Vertex *iVertex, const VertexCol *iVertices)
{
  // A PF Jet with L1+L2+L3 corrections is expected.
  if (iJet->Corrections() != fgCorrectionMask)
    throw std::runtime_error("JetIDMVA works only with L1+L2+L3 corrected jets");

  double lEta = iJet->AbsEta();
  if (lEta > 4.99)
    return false;

  if(!JetTools::passPFLooseId(iJet))
    return false;

  double lPt = iJet->Pt(); // use corrected Pt
  int lPtId = 0;
  if (lPt > 10. && lPt < 20.)
    lPtId = 1;
  else if (lPt < 30.)
    lPtId = 2;
  else
    lPtId = 3;

  int lEtaId = 0;
  if (lEta > 2.5 && lEta < 2.75)
    lEtaId = 1;
  else if (lEta < 3.)
    lEtaId = 2;
  else
    lEtaId = 3;

  if (fType == kCut) {
    float betaStarModified = JetTools::betaStarClassic(iJet,iVertex,iVertices)/log(iVertices ->GetEntries()-0.64);
    float dR2Mean          = JetTools::dR2Mean(iJet,-1);

    if(betaStarModified < fBetaStarCut[lPtId][lEtaId] &&
       dR2Mean          < fRMSCut[lPtId][lEtaId])
      return true;
  }
  else {
    double lMVA = MVAValue(iJet, iVertex, iVertices);
    double lMVACut = fMVACut[lPtId][lEtaId];
    if (lMVA > lMVACut)
      return true;
  }

  return false;
}

//--------------------------------------------------------------------------------------------------
Double_t
JetIDMVA::MVAValue(const PFJet *iJet, const Vertex *iVertex, //Vertex here is the PV
                   const VertexCol *iVertices,
                   Bool_t printDebug)
{
  if (!fIsInitialized)
    throw std::runtime_error("Error: JetIDMVA not properly initialized.");

  // A PF Jet with L1+L2+L3 corrections is expected.
  if (iJet->Corrections() != fgCorrectionMask)
    throw std::runtime_error("JetIDMVA works only with L1+L2+L3 corrected jets");

  if (!JetTools::passPFLooseId(iJet))
    return -2.;

  //set all input variables
  fVariables[kNvtx]      = iVertices->GetEntries();
  fVariables[kJetPt]     = iJet->Pt();
  fVariables[kJetEta]    = iJet->Eta();
  fVariables[kJetPhi]    = iJet->Phi();
  fVariables[kD0]        = JetTools::impactParameter(iJet, iVertex);
  fVariables[kDZ]        = JetTools::impactParameter(iJet, iVertex, true);
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
  fVariables[kPtD]       = JetTools::W(iJet, -1, 0);
  fVariables[kDRMean]    = JetTools::dRMean(iJet, -1);
  fVariables[kDR2Mean]   = JetTools::dR2Mean(iJet, -1);
  fVariables[kFrac01]    = JetTools::frac(iJet, 0.1, 0.,  -1);
  fVariables[kFrac02]    = JetTools::frac(iJet, 0.2, 0.1, -1);
  fVariables[kFrac03]    = JetTools::frac(iJet, 0.3, 0.2, -1);
  fVariables[kFrac04]    = JetTools::frac(iJet, 0.4, 0.3, -1);
  fVariables[kFrac05]    = JetTools::frac(iJet, 0.5, 0.4, -1);

  double lMVA = 0;
  lMVA = fReader->EvaluateMVA(fMethodName);

  if (printDebug) {
    std::cout << "Debug Jet MVA:" << std::endl;
    for (unsigned iV = 0; iV != nVariables; ++iV)
      std::cout << fVarNames[iV] << " = " << fVariables[iV] << std::endl;
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
  if (!JetTools::passPFLooseId(iJet))
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
  fVariables[kPtD]       = JetTools::W(iJet,-1,0);
  fVariables[kDRMean]    = JetTools::dRMean(iJet,-1);
  fVariables[kDR2Mean]   = JetTools::dR2Mean(iJet,-1);
  fVariables[kFrac01]    = JetTools::frac(iJet,0.1,0., -1);
  fVariables[kFrac02]    = JetTools::frac(iJet,0.2,0.1,-1);
  fVariables[kFrac03]    = JetTools::frac(iJet,0.3,0.2,-1);
  fVariables[kFrac04]    = JetTools::frac(iJet,0.4,0.3,-1);
  fVariables[kFrac05]    = JetTools::frac(iJet,0.5,0.4,-1);

  double lMVA = 0;
  lId[0] = fReader->EvaluateMulticlass(fMethodName)[0];
  lId[1] = fReader->EvaluateMulticlass(fMethodName)[1];
  lId[2] = fReader->EvaluateMulticlass(fMethodName)[2];
  if (printDebug) {
    std::cout << "Debug Jet MVA:" << std::endl;
    for (unsigned iV = 0; iV != nVariables; ++iV)
      std::cout << fVarNames[iV] << " = " << fVariables[iV] << std::endl;
    std::cout << "=== : === " << std::endl;
    std::cout << "MVA value " << lMVA << std::endl;
  }
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
