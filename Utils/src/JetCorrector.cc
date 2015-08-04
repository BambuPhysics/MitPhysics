#include "MitPhysics/Utils/interface/JetCorrector.h"
#include "MitAna/DataTree/interface/CaloJet.h"
#include "MitAna/DataTree/interface/ObjTypes.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include <stdexcept>
#include <sys/stat.h>

ClassImp(mithep::JetCorrector)

mithep::JetCorrector::~JetCorrector()
{
  delete fCorrector;
}

void
mithep::JetCorrector::AddParameterFile(char const* fileName)
{
  mithep::Jet::ECorr lastLevel = fLevels.back();
  if (lastLevel == mithep::Jet::Custom) { // = Uncertainty
    std::cerr << "JetCorrector cannot add a correction level on top of Uncertainty." << std::endl;
    throw std::runtime_error("Configuration error");
  }

  struct stat buffer;
  if (stat(fileName, &buffer) != 0) {
    std::cerr << "File " << fileName << " does not exist." << std::endl;
    throw std::runtime_error("Configuration error");
  }

  try {
    fParameters.emplace_back(std::string(fileName));
  }
  catch (std::exception& ex) {
    std::cerr << "Exception in JetCorrector::AddParameterFile(" << fileName << "):" << std::endl;
    std::cerr << ex.what() << std::endl;
    throw;
  }

  auto it = fParameters.rbegin();
  mithep::Jet::ECorr level = TranslateLevel(it->definitions().level().c_str());
  if (level == mithep::Jet::Custom && lastLevel != mithep::Jet::L2 && lastLevel != mithep::Jet::L3) {
    // L2 for case of data L2L3 residual correction
    std::cerr << "Uncertainty is currently only available on top of L3 (+L2L3Residual) correction." << std::endl;
    throw std::runtime_error("Configuration error");
  }

  fLevels.push_back(level);

// Not enforcing ordered corrections because L2L3Residual is seen as L2 by the JetCorrectionParameters
//   // check correction ordering
//   lastLevel = level;
//   if (fParameters.size() > 1) {
//     ++it;
//     if (TranslateLevel(it->definitions().level().c_str()) >= lastLevel) {
//       std::cerr << "Exception in JetCorrector::AddParameterFile(" << fileName << "):" << std::endl;
//       std::cerr << "Correction parameters must be added in ascending order of correction levels" << std::endl;
//       throw std::runtime_error("Configuration error");
//     }
//   }
}

void
mithep::JetCorrector::ClearParameters()
{
  fParameters.clear();
  fLevels.clear();
}

void
mithep::JetCorrector::Initialize()
{
  delete fCorrector;
  fCorrector = new FactorizedJetCorrector(fParameters);
}

std::vector<Float_t>
mithep::JetCorrector::CorrectionFactors(mithep::Jet& jet, Double_t rho/* = 0.*/) const
{
  if (!IsInitialized())
    throw std::runtime_error("JetCorrector not initialized");

  auto&& rawMom = jet.RawMom();

  //compute correction factors
  fCorrector->setJetEta(rawMom.Eta());
  fCorrector->setJetPt(rawMom.Pt());
  fCorrector->setJetPhi(rawMom.Phi());
  fCorrector->setJetE(rawMom.E());

  fCorrector->setRho(rho);
  fCorrector->setJetA(jet.JetArea());
    
  //emf only valid for CaloJets
  if (jet.ObjType() == mithep::kCaloJet)
    fCorrector->setJetEMF(static_cast<mithep::CaloJet&>(jet).EnergyFractionEm());
  else
    fCorrector->setJetEMF(-99.0);

  std::vector<float>&& corrections(fCorrector->getSubCorrections());

  // if MaxCorrLevel is specified, downsize the corrections array
  // for data L1+L2+L3+L2L3Residual, the last (residual) appears as L2
  // therefore must check for level <= fMaxCorrLevel
  // also if uncertainty (using Custom) is required that is the last level
  if (fMaxCorrLevel != Jet::nECorrs) {
    for (unsigned iL = fLevels.size(); iL != 0; --iL) {
      if (fLevels[iL - 1] <= fMaxCorrLevel || fLevels[iL - 1] == mithep::Jet::Custom) {
        corrections.resize(iL);
        break;
      }
    }
  }

  return corrections;
}

Float_t
mithep::JetCorrector::CorrectionFactor(mithep::Jet& jet, Double_t rho/* = 0.*/) const
{
  auto&& factors(CorrectionFactors(jet, rho));
  if (factors.size() != 0)
    return factors.back();
  else
    return 1.;
}

void
mithep::JetCorrector::Correct(mithep::Jet& jet, Double_t rho/* = 0.*/) const
{
  //set and enable correction factors in the output jet

  auto&& corrections(CorrectionFactors(jet, rho));

  auto lastLevel = mithep::Jet::nECorrs;

  for (unsigned iC = 0; iC != corrections.size(); ++iC) {
    auto currentLevel = fLevels.at(iC);
    float currentCorrection = 1.;
    if (lastLevel == mithep::Jet::L3 && currentLevel == mithep::Jet::L2) {
      // special case for L2L3Residual
      // store L3*L2L3 as L3 correction
      if (iC > 1)
        currentCorrection = corrections.at(iC) / corrections.at(iC - 2);
      else
        currentCorrection = corrections.at(iC);

      currentLevel = mithep::Jet::L3;
    }
    else if (currentLevel == mithep::Jet::Custom) {
      currentLevel = lastLevel;
      // last cumulative (up to L3 or L2L3) + fSigma * uncertainty
      currentCorrection = corrections.at(iC - 1) + fSigma * corrections.at(iC);
    }
    else {
      if (iC > 0)
        currentCorrection = corrections.at(iC) / corrections.at(iC - 1);
      else
        currentCorrection = corrections.at(iC);
    }

    //set correction and enable
    jet.SetCorrectionScale(currentCorrection, currentLevel);
    jet.EnableCorrection(currentLevel);

    lastLevel = currentLevel;
  }
}

/*static*/
mithep::Jet::ECorr
mithep::JetCorrector::TranslateLevel(char const* levelName)
{
  std::string name(levelName);
  if (name == "L1Offset" || name == "L1FastJet" || name == "L1JPTOffset")
    return mithep::Jet::L1;
  else if (name == "L2Relative")
    return mithep::Jet::L2;
  else if (name == "L3Absolute")
    return mithep::Jet::L3;
  else if (name == "L4EMF")
    return mithep::Jet::L4;
  else if (name == "L5Flavor")
    return mithep::Jet::L5;
  else if (name == "L6SLB")
    return mithep::Jet::L6;
  else if (name == "L7Parton")
    return mithep::Jet::L7;
  else if (name == "Uncertainty")
    return mithep::Jet::Custom;
  else {
    std::cerr << "Exception in JetCorrector::TranslateLevel(): Unknown correction level " << name << std::endl;
    throw std::runtime_error("Unknown correction");
  }
}
