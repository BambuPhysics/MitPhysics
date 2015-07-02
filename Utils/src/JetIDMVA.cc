#include "MitPhysics/Utils/interface/JetIDMVA.h"
#include "MitPhysics/Utils/interface/JetTools.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include <stdexcept>

ClassImp(mithep::JetIDMVA)

using namespace mithep;

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
  
  fIsInitialized = kTRUE;
  fType = mvaType;

  std::string lCutId;
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

  std::string lCutType;
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

  //Load Cut Matrix
  std::ifstream cutFile(cutConfig);

  bool cutsSet[2] = {false, false};

  while (true) {
    std::string line;
    std::getline(cutFile, line);
    if (!cutFile.good())
      break;

    // restream the line into words
    std::stringstream ss;
    ss.str(line);
    std::string word;
    ss >> word;
    if (word != lCutId)
      continue;
    ss >> word;
    if (word.find(lCutType) == std::string::npos)
      continue;

    // A line with lCutId and lCutType is found.
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

    // The following four lines in the config define cut values

    unsigned iL = 0;
    for (; iL != 4; ++iL) {
      std::getline(cutFile, line);
      if (!cutFile.good()) // input format not correct
        break;
      
      ss.str(line);
      // four bins in eta
      unsigned iE = 0;
      for (; iE != 4; ++iE) {
        ss >> cutArray[iL][iE];
        if (!ss.good())
          break;
      }
      if (iE != 4) // input format not correct
        break;
    }
    if (iL != 4) {
      // input format not correct
      throw std::runtime_error(("JetIDMVA could not parse " + cutConfig).Data());
    }

    // config correctly read. set flag.
    *flag = true;

    if (fType != kCut || (cutsSet[0] && cutsSet[1])) // nothing more to read
      break;
  }

  cutFile.close();

  if (fType == kCut) {
    if (!cutsSet[0] || !cutsSet[1])
      throw std::runtime_error("JetIDMVA cut values not set!");

    return;
  }

  if (!cutsSet[0])
    throw std::runtime_error("JetIDMVA cut values not set!");
  
  fReader        = new TMVA::Reader("!Color:!Silent:Error" );
  fReader->AddVariable("nvtx",      fVariables + kNVtx); 
  fReader->AddVariable("jetPt",     fVariables + kJPt1);  
  fReader->AddVariable("jetEta",    fVariables + kJEta1);
  fReader->AddVariable("dZ",        fVariables + kJDZ1);
  fReader->AddVariable("beta",      fVariables + kBeta);
  fReader->AddVariable("betaStar",  fVariables + kBetaStar);
  fReader->AddVariable("nCharged",  fVariables + kNCharged);
  fReader->AddVariable("nNeutrals", fVariables + kNNeutrals);
  fReader->AddVariable("frac01",    fVariables + kFrac01);
  fReader->AddVariable("frac02",    fVariables + kFrac02);
  fReader->AddVariable("frac03",    fVariables + kFrac03);
  fReader->AddVariable("frac04",    fVariables + kFrac04);
  fReader->AddVariable("frac05",    fVariables + kFrac05);

  if (fType == kBaseline) {
    fReader->AddVariable("jetPhi", fVariables + kJPhi1);             
    fReader->AddVariable("d0",     fVariables + kJD01);
    fReader->AddVariable("dRMean", fVariables + kDRMean);
  }
  if (fType == k52) {
    fReader->AddVariable("dR2Mean", fVariables + kDR2Mean);
  }
  if (fType == k53) {
    fReader->AddVariable("dR2Mean", fVariables + kDR2Mean);
    fReader->AddVariable("ptD",     fVariables + kPtD);
    fReader->AddSpectator("jetPhi", fVariables + kJPhi1);  
  } 
  if (fType == k53MET) {
    fReader->AddVariable("jetPhi",  fVariables + kJPhi1);  
    fReader->AddVariable("dR2Mean", fVariables + kDR2Mean);
    fReader->AddVariable("ptD",     fVariables + kPtD);
  } 
  if (fType == kQGP) {
    fReader->AddVariable("jetPhi", fVariables + kJPhi1);             
    fReader->AddVariable("dRMean", fVariables + kDRMean);
    fReader->AddVariable("ptD",    fVariables + kPtD);
  }

  fReader->BookMVA(fMethodName, weightsConfig);
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
  if(!JetTools::passPFLooseId(iJet))
    return false;
  if(iJet->Pt() < fJetPtMin)
    return false; 
  if(iJet->AbsEta() > 4.99)
    return false;

  double lPt = iJet->Pt();  
  int lPtId = 0; 
  if (lPt > 10. && lPt < 20.)
    lPtId = 1;
  else if (lPt < 30.)
    lPtId = 2;
  else
    lPtId = 3;
  
  double lEta = iJet->AbsEta();
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
       dR2Mean          < fRMSCut     [lPtId][lEtaId])
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

  if (!JetTools::passPFLooseId(iJet))
    return -2.;

  //set all input variables
  fVariables[kNVtx]      = iVertices->GetEntries();
  fVariables[kJPt1]      = iJet->Pt();
  fVariables[kJEta1]     = iJet->RawMom().Eta();
  fVariables[kJPhi1]     = iJet->RawMom().Phi();
  fVariables[kJD01]      = JetTools::impactParameter(iJet,iVertex);  
  fVariables[kJDZ1]      = JetTools::impactParameter(iJet,iVertex,true);
  fVariables[kBeta]      = JetTools::Beta(iJet,iVertex,fDZCut);
  fVariables[kBetaStar]  = JetTools::betaStar(iJet,iVertex,iVertices,fDZCut);
  fVariables[kNCharged]  = iJet->ChargedMultiplicity();
  fVariables[kNNeutrals] = iJet->NeutralMultiplicity();
  fVariables[kPtD]       = JetTools::W(iJet,-1,0);  
  fVariables[kDRMean]    = JetTools::dRMean (iJet,-1);
  fVariables[kDR2Mean]   = JetTools::dR2Mean(iJet,-1);
  fVariables[kFrac01]    = JetTools::frac(iJet,0.1,0., -1);
  fVariables[kFrac02]    = JetTools::frac(iJet,0.2,0.1,-1);
  fVariables[kFrac03]    = JetTools::frac(iJet,0.3,0.2,-1);
  fVariables[kFrac04]    = JetTools::frac(iJet,0.4,0.3,-1);
  fVariables[kFrac05]    = JetTools::frac(iJet,0.5,0.4,-1);

  double lMVA = 0;
  lMVA = fReader->EvaluateMVA(fMethodName);
  
  if (printDebug == kTRUE) {
    std::cout << "Debug Jet MVA: "
	      << fVariables[kNVtx]      << " "
	      << fVariables[kJPt1]      << " "
	      << fVariables[kJEta1]     << " "
	      << fVariables[kJPhi1]     << " "
	      << fVariables[kJD01]      << " "
	      << fVariables[kJDZ1]      << " "
	      << fVariables[kBeta]      << " "
	      << fVariables[kBetaStar]  << " "
	      << fVariables[kNCharged]  << " "
	      << fVariables[kNNeutrals] << " "
	      << fVariables[kDRMean]    << " "
	      << fVariables[kPtD]       << " "
	      << fVariables[kFrac01]    << " "
	      << fVariables[kFrac02]    << " "
	      << fVariables[kFrac03]    << " "
	      << fVariables[kFrac04]    << " "
	      << fVariables[kFrac05]    << " "
	      << fVariables[kDR2Mean]    
              << " === : === "
              << lMVA << " "    
              << std::endl;
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

  fVariables[kJPt1]       = iJet->Pt();
  if (fVariables[kJPt1] < 20)
    return lId;

  //set all input variables
  fVariables[kNVtx]       = iVertices->GetEntries();
  fVariables[kJEta1]      = iJet->RawMom().Eta();
  fVariables[kJPhi1]      = iJet->RawMom().Phi();
  fVariables[kJD01]       = JetTools::impactParameter(iJet,iVertex);  
  fVariables[kJDZ1]       = JetTools::impactParameter(iJet,iVertex,true);
  fVariables[kBeta]       = JetTools::Beta(iJet,iVertex,fDZCut);
  fVariables[kBetaStar]   = JetTools::betaStar(iJet,iVertex,iVertices,fDZCut);
  fVariables[kNCharged]   = iJet->ChargedMultiplicity();
  fVariables[kNNeutrals]  = iJet->NeutralMultiplicity();
  fVariables[kNParticles] = fVariables[kNCharged] + fVariables[kNNeutrals];
  fVariables[kPtD]        = JetTools::W(iJet,-1,0);  
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
  if (printDebug == kTRUE) {
    std::cout << "Debug Jet MVA: "
	      << fVariables[kNVtx]      << " "
	      << fVariables[kJPt1]      << " "
	      << fVariables[kJEta1]     << " "
	      << fVariables[kJPhi1]     << " "
	      << fVariables[kJD01]      << " "
	      << fVariables[kJDZ1]      << " "
	      << fVariables[kBeta]      << " "
	      << fVariables[kBetaStar]  << " "
	      << fVariables[kNCharged]  << " "
	      << fVariables[kNNeutrals] << " "
	      << fVariables[kDRMean]    << " "
	      << fVariables[kFrac01]    << " "
	      << fVariables[kFrac02]    << " "
	      << fVariables[kFrac03]    << " "
	      << fVariables[kFrac04]    << " "
	      << fVariables[kFrac05]    << " "
	      << fVariables[kDRMean]    
              << " === : === "
              << lMVA << " "    
              << std::endl;
  }
  return lId;
}
