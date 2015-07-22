#include "MitPhysics/Utils/interface/JetCorrector.h"
#include "MitAna/DataTree/interface/CaloJet.h"
#include "MitAna/DataTree/interface/ObjTypes.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include <stdexcept>

ClassImp(mithep::JetCorrector)

mithep::JetCorrector::~JetCorrector()
{
  delete fCorrector;
}

void
mithep::JetCorrector::AddParameterFile(char const* fileName)
{
  try {
    fParameters.emplace_back(std::string(fileName));
  }
  catch (std::exception& ex) {
    std::cerr << "Exception in JetCorrector::AddParameterFile(" << fileName << "):" << std::endl;
    std::cerr << ex.what() << std::endl;
    throw;
  }

  auto it = fParameters.rbegin();
  mithep::Jet::ECorr lastLevel = TranslateLevel(it->definitions().level().c_str());
  fLevels.push_back(lastLevel);

// Not enforcing ordered corrections because L2L3Residual is seen as L2 by the JetCorrectionParameters
//   // check correction ordering
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
mithep::JetCorrector::Initialize()
{
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

  std::vector<float>&& corrections = fCorrector->getSubCorrections();
  if (fMaxCorrLevel != Jet::nECorrs) {
    for (unsigned iL = 1; iL <= fLevels.size(); ++iL) {
      if (fLevels[iL - 1] == fMaxCorrLevel) {
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
  auto&& corrections(CorrectionFactors(jet, rho));

  //set and enable correction factors in the output jet
  double cumulativeCorrection = 1.0;

  for (unsigned iC = 0; iC != corrections.size(); ++iC) {
    float currentCorrection = corrections.at(iC) / cumulativeCorrection;
    cumulativeCorrection = corrections.at(iC);

    //set correction and enable
    jet.SetCorrectionScale(currentCorrection, fLevels.at(iC));
    jet.EnableCorrection(fLevels.at(iC));
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
  else {
    std::cerr << "Exception in JetCorrector::TranslateLevel(): Unknown correction level " << name << std::endl;
    throw std::runtime_error("Unknown correction");
  }
}
