#include "MitPhysics/Utils/interface/JetCorrector.h"
#include "MitAna/DataTree/interface/CaloJet.h"
#include "MitAna/DataTree/interface/ObjTypes.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include <stdexcept>
#include <sys/stat.h>
#include <fstream>

#include "TROOT.h"
#include "TVector2.h"

ClassImp(mithep::JetCorrector)

mithep::JetCorrector::Corrector::Record::Record(unsigned nBinVars, unsigned nFormVars, std::string const& inputLine)
{
  // Record line formats
  //  [JEC]
  //   {binning (n = 2*nBinVars)} {# remaining words} {variable range (n = 2*nFormVars)} {parameters}
  //  [JEC uncertainty]
  //   {binning (n = 2*nBinVars)} {# remaining words} ({variable min} {uncertainty up} {uncertainty down}) x {# / 3}
  //  [JER scale factor]
  //   {binning (n = 2*nBinVars)} {# remaining words} {sf nominal} {sf down} {sf up}
  //  [JER resolution]
  //   {binning (n = 2*nBinVars)} {# remaining words} {variable range (n = 2*nFormVars)} {parameters}

  auto* words(TString(inputLine.c_str()).Tokenize(" "));

  auto readWord([words](unsigned idx)->double {
      return TString(words->At(idx)->GetName()).Atof();
    });

  unsigned iW(0);

  for (unsigned iV(0); iV != nBinVars; ++iV) {
    fBinVarLimits.emplace_back(readWord(iW), readWord(iW + 1));
    iW += 2;
  }

  int nWords(iW + readWord(iW));
  if (nWords != words->GetEntries()) {
    fBinVarLimits.pop_back();
    delete words;
    throw std::exception();
  }

  ++iW;

  for (unsigned iV(0); iV != nFormVars; ++iV) {
    fFormVarLimits.emplace_back(readWord(iW), readWord(iW + 1));
    iW += 2;
  }

  for (unsigned iP(0); iP != nWords - nBinVars * 2; ++iP)
    fParameters.push_back(readWord(iW++));

  delete words;
}

mithep::JetCorrector::Corrector::Corrector(char const* fileName, FactorType factorType/* = nFactorTypes*/)
{
  struct stat buffer;
  if (stat(fileName, &buffer) != 0) {
    std::cerr << "File " << fileName << " does not exist." << std::endl;
    throw std::runtime_error("Configuration error");
  }

  std::ifstream source(fileName);
  std::string linebuf;
  unsigned iL(0);

  while (std::getline(source, linebuf)) {
    ++iL;

    // read in the definition line
    std::string::size_type open(linebuf.find('{'));
    std::string::size_type close(linebuf.find('}'));
    if (open != std::string::npos && close != std::string::npos && open < close) {
      auto* words(TString(linebuf.substr(open + 1, close - open - 1).c_str()).Tokenize(" "));
      if (words->GetEntries() < 6) {
        std::cerr << "File " << fileName << " does not contain a proper JEC definition." << std::endl;
        delete words;
        throw std::runtime_error("Configuration error");
      }

      unsigned iW(0);

      unsigned nBinVar(TString(words->At(iW++)->GetName()).Atoi());
      if (nBinVar == 0) {
        std::cerr << "File " << fileName << " does not contain a proper JEC definition." << std::endl;
        delete words;
        throw std::runtime_error("Configuration error");
      }

      for (unsigned iV(0); iV != nBinVar; ++iV)
        fBinningVariables.push_back(ParseVariableName(words->At(iW++)->GetName()));

      unsigned nFormVar(TString(words->At(iW++)->GetName()).Atoi());

      for (unsigned iV(0); iV != nFormVar; ++iV)
        fFormulaVariables.push_back(ParseVariableName(words->At(iW++)->GetName()));

      TString formula(words->At(iW++)->GetName());
      if (formula != "" && formula != "None") {
        // TFormula has a global namespace - keep one copy from each input file (which may be used multiple times in a job)
        auto* func(gROOT->GetListOfFunctions()->FindObject(fileName));
        if (func)
          fFormula = static_cast<TFormula*>(func);
        else
          fFormula = new TFormula(fileName, formula);

        if (fFormula->GetNdim() != int(fFormulaVariables.size())) {
          std::cerr << "File " << fileName << " does not contain a proper JEC formula." << std::endl;
          delete words;
          throw std::runtime_error("Configuration error");
        }
      }

      fType = ParseTypeName(words->At(iW++)->GetName());

      if (fType == kCorrection || fType == kResponse) {
        fLevel = ParseLevelName(words->At(iW)->GetName());

        if (fLevel != mithep::Jet::Custom && !fFormula) {
          std::cerr << "Cannot initialize JEC without a formula." << std::endl;
          delete words;
          throw std::runtime_error("Configuration error");
        }
        if (fLevel == mithep::Jet::Custom && fFormulaVariables.size() != 1) {
          std::cerr << "JES uncertainty must be given as a grid of two variables." << std::endl;
          delete words;
          throw std::runtime_error("Configuration error");
        }
      }
      else if (fType == kResolution) {
        if (factorType == nFactorTypes) {
          std::cerr << "Resolution type (JetCorrector::Corrector::kPtResolution or kPhiResolution) must be set for resolution inputs." << std::endl;
          // otherwise the code cannot know if the returned value is a relative or an absolute resolution
          delete words;
          throw std::runtime_error("Configuration error");
        }

        fType = factorType;
      }

      delete words;

      break;
    }
  }

  while (std::getline(source, linebuf)) {
    ++iL;

    try {
      if (fType == kScaleFactor || (fType == kCorrection && fLevel == mithep::Jet::Custom))
        fRecords.emplace_back(fBinningVariables.size(), 0, linebuf);
      else
        fRecords.emplace_back(fBinningVariables.size(), fFormulaVariables.size(), linebuf);
    }
    catch (std::exception& ex) {
      std::cerr << "File " << fileName << " appears corrupt at line " << iL << std::endl;
      throw std::runtime_error("Configuration error");
    }
  }
}

Double_t
mithep::JetCorrector::Corrector::Eval(Double_t variables[nVarTypes]) const
{
  switch (fType) {
  case kCorrection:
    if (fLevel == mithep::Jet::Custom)
      return EvalUncertainty(variables);
    else
      return EvalCorrection(variables);

  case kResponse:
    std::cerr << "Response inversion not implemented yet. You are welcome to write it for us! (ref. CondFormats/JetMETObjects/src/SimpleJetCorrector.cc" << std::endl;
    throw std::runtime_error("NotImplemented");

  case kPtResolution:
    return EvalCorrection(variables) * variables[kJetPt];

  case kPhiResolution:
    return EvalCorrection(variables);

  case kScaleFactor:
    return EvalScaleFactor(variables);

  default:
    std::cerr << "Corrector type is not specified" << std::endl;
    throw std::runtime_error("NotImplemented");    
  }
}

Double_t
mithep::JetCorrector::Corrector::EvalCorrection(Double_t variables[nVarTypes]) const
{
  unsigned iRecord(FindRecord(variables));

  if (iRecord == fRecords.size()) // no record found
    return 1.;
      
  if (fInterpolate) {
    std::cerr << "JEC interpolation not implemented yet. You are welcome to write it for us! (ref. CondFormats/JetMETObjects/src/SimpleJetCorrector.cc" << std::endl;
    throw std::runtime_error("NotImplemented");
  }
  else {
    auto& record(fRecords[iRecord]);

    double x[4];
    std::vector<double> par;

    for (unsigned iX(0); iX != fFormulaVariables.size(); ++iX) {
      double val(variables[fFormulaVariables[iX]]);
      auto& limits(record.fFormVarLimits[iX]);
      if (val < limits.first)
        x[iX] = limits.first;
      else if (val >= limits.second)
        x[iX] = limits.second;
      else
        x[iX] = val;
    }

    for (unsigned iP(0); iP != record.fParameters.size(); ++iP)
      par.push_back(record.fParameters[iP]);

    return fFormula->EvalPar(x, par.data());
  }
}

Double_t
mithep::JetCorrector::Corrector::EvalUncertainty(Double_t variables[nVarTypes]) const
{
  unsigned iRecord(FindRecord(variables));

  if (iRecord == fRecords.size()) {
    // no record found
    std::cerr << "JetCorrector uncertainty out of range" << std::endl;
    return -999.;
  }

  auto& record(fRecords[iRecord]);

  double var(variables[fFormulaVariables[0]]);

  unsigned iP(0);
  unsigned iX(0);
  while (iP < record.fParameters.size() && var < record.fParameters[iP]) {
    iP += 3; // parameter array: pt0, uncertUp0, uncertDown0, pt1, uncertUp1, ...
    ++iX;
  }

  unsigned pOffset(fUncertaintyV >= 0 ? 1 : 2);

  if (iP == 0)
    return record.fParameters.at(pOffset);
  else if (iP >= record.fParameters.size())
    return record.fParameters.at(iX * 3 - 3 + pOffset);
  else {
    // linear interpolation
    double low(record.fParameters.at(iX * 3 - 3));
    double high(record.fParameters.at(iX * 3));
    double interval(high - low);
    return (record.fParameters.at(iX * 3 + pOffset) * (var - low) + record.fParameters.at(iX * 3 - 3 + pOffset) * (high - var)) / interval;
  }
}

Double_t
mithep::JetCorrector::Corrector::EvalScaleFactor(Double_t variables[nVarTypes]) const
{
  unsigned iRecord(FindRecord(variables));

  if (iRecord == fRecords.size()) {
    // no record found
    std::cerr << "JetCorrector uncertainty out of range" << std::endl;
    return 1.;
  }

  auto& record(fRecords[iRecord]);
  
  switch (fUncertaintyV) {
  case 0:
    return record.fParameters.at(0);
  case -1:
    return record.fParameters.at(1);
  case 1:
    return record.fParameters.at(2);
  default:
    return record.fParameters.at(0);
  };
}

UInt_t
mithep::JetCorrector::Corrector::FindRecord(Double_t variables[nVarTypes]) const
{
  unsigned iRecord(0);
  for (; iRecord != fRecords.size(); ++iRecord) {
    auto& record(fRecords[iRecord]);

    unsigned iVar(0);
    for (; iVar != fBinningVariables.size(); ++iVar) {
      double value(variables[fBinningVariables[iVar]]);
      
      if (record.fBinVarLimits[iVar].first > value || record.fBinVarLimits[iVar].second <= value)
        break;
    }
    if (iVar == fBinningVariables.size()) // all bin variables are in the range
      break;
  }

  return iRecord;
}

/*static*/
mithep::Jet::ECorr
mithep::JetCorrector::Corrector::ParseLevelName(char const* levelName)
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
    std::cerr << "Exception in JetCorrector::ParseLevelName(): Unknown correction level " << name << std::endl;
    throw std::runtime_error("Unknown correction");
  }
}

/*static*/
mithep::JetCorrector::Corrector::VarType
mithep::JetCorrector::Corrector::ParseVariableName(char const* varName)
{
  std::string name(varName);
  if (name == "JetEta")
    return mithep::JetCorrector::Corrector::kJetEta;
  else if (name == "NPV")
    return mithep::JetCorrector::Corrector::kNPV;
  else if (name == "JetPt")
    return mithep::JetCorrector::Corrector::kJetPt;
  else if (name == "JetPhi")
    return mithep::JetCorrector::Corrector::kJetPhi;
  else if (name == "JetE")
    return mithep::JetCorrector::Corrector::kJetE;
  else if (name == "JetEMF")
    return mithep::JetCorrector::Corrector::kJetEMF;
  else if (name == "JetA")
    return mithep::JetCorrector::Corrector::kJetA;
  else if (name == "Rho")
    return mithep::JetCorrector::Corrector::kRho;
  else if (name == "RelLepPt")
    return mithep::JetCorrector::Corrector::kRelLepPt;
  else if (name == "PtRel")
    return mithep::JetCorrector::Corrector::kPtRel;
  else {
    std::cerr << "Exception in JetCorrector::ParseVariableName(): Unknown variable " << name << std::endl;
    throw std::runtime_error("Unknown correction");
  }
}

/*static*/
mithep::JetCorrector::Corrector::FactorType
mithep::JetCorrector::Corrector::ParseTypeName(char const* typeName)
{
  TString name(typeName);
  name.ToLower();

  if (name == "correction")
    return mithep::JetCorrector::Corrector::kCorrection;
  else if (name == "response")
    return mithep::JetCorrector::Corrector::kResponse;
  else if (name == "resolution")
    return mithep::JetCorrector::Corrector::kResolution;
  else if (name == "scalefactor")
    return mithep::JetCorrector::Corrector::kScaleFactor;
  else {
    std::cerr << "Exception in JetCorrector::ParseTypeName(): Unknown factor type " << name << std::endl;
    throw std::runtime_error("Unknown correction");
  }
}

void
mithep::JetCorrector::AddParameterFile(char const* fileName, Corrector::FactorType factorType/* = Corrector::nFactorTypes*/)
{
  try {
    fCorrectors.emplace_back(fileName, factorType);
  }
  catch (std::exception& ex) {
    std::cerr << "Exception in JetCorrector::AddParameterFile(" << fileName << "):" << std::endl;
    std::cerr << ex.what() << std::endl;
    throw;
  }

  auto type = fCorrectors.back().GetType();
  if (type == Corrector::kPtResolution)
    fHasSmearing = true;
}

mithep::JetCorrector::JetCorrector()
{
  // Reserve enough space for all correctors. Otherwise emplace_back will
  // attempt a re-construction of the correctors at some point, which leads to mysterious segfaults.
  fCorrectors.reserve(mithep::Jet::nECorrs);
}

void
mithep::JetCorrector::ClearParameters()
{
  fCorrectors.clear();
  fMaxCorrLevel = mithep::Jet::nECorrs;
  fSigma = 0.;
  fHasSmearing = false;
}

std::vector<Float_t>
mithep::JetCorrector::CorrectionFactors(mithep::Jet const& jet, Double_t rho/* = 0.*/) const
{
  auto&& rawMom(jet.RawMom());

  //compute correction factors
  double vars[Corrector::nVarTypes];
  vars[Corrector::kJetEta] = rawMom.Eta();
  vars[Corrector::kJetPt] = rawMom.Pt();
  vars[Corrector::kJetPhi] = rawMom.Phi();
  vars[Corrector::kJetE] = rawMom.E();
  vars[Corrector::kRho] = rho;
  vars[Corrector::kJetA] = jet.JetArea();
    
  //emf only valid for CaloJets
  if (jet.ObjType() == mithep::kCaloJet)
    vars[Corrector::kJetEMF] = static_cast<mithep::CaloJet const&>(jet).EnergyFractionEm();
  else
    vars[Corrector::kJetEMF] = -99.0;

  std::vector<float> corrections;
  double factor(1.);
  for (auto& corr : fCorrectors) {
    if (corr.GetLevel() > fMaxCorrLevel)
      break;

    if (corr.GetEnabled()) {
      double scale(corr.Eval(vars));
      vars[Corrector::kJetPt] *= scale;
      vars[Corrector::kJetE] *= scale;
      factor *= scale;
    }

    corrections.push_back(factor);
  }

  return corrections;
}

Float_t
mithep::JetCorrector::CorrectionFactor(mithep::Jet const& jet, Double_t rho/* = 0.*/) const
{
  auto&& factors(CorrectionFactors(jet, rho));
  if (factors.size() != 0)
    return factors.back();
  else
    return 1.;
}

Float_t
mithep::JetCorrector::UncertaintyFactor(mithep::Jet const& jet) const
{
  Corrector const* uncertainty(0);
  for (auto& corr : fCorrectors) {
    if (corr.GetLevel() == mithep::Jet::Custom) {
      uncertainty = &corr;
      break;
    }
  }

  if (!uncertainty)
    return 1.;

  // uncertainty is calculated on fully corrected momentum
  double vars[Corrector::nVarTypes];
  vars[Corrector::kJetEta] = jet.Eta();
  vars[Corrector::kJetPt] = jet.Pt();

  // last cumulative (up to L3 or L2L3) + fSigma * uncertainty
  // getUncertainty(true): Upside uncertainty.
  // Usually downside is not used; probably the uncertainty is symmetric anyway..
  return 1. + fSigma * uncertainty->Eval(vars);
}

Bool_t
mithep::JetCorrector::IsEnabled(mithep::Jet::ECorr l) const
{
  if (l >= fMaxCorrLevel)
    return false;

  for (auto& corr : fCorrectors) {
    if (corr.GetLevel() == l && corr.GetEnabled())
      return true;
  }

  return false;
}

void
mithep::JetCorrector::Correct(mithep::Jet& jet, Double_t rho/* = 0.*/) const
{
  //set and enable correction factors in the output jet

  auto&& corrections(CorrectionFactors(jet, rho));

  auto lastLevel = mithep::Jet::nECorrs;

  for (unsigned iC(0); iC != corrections.size(); ++iC) {
    auto currentLevel(fCorrectors[iC].GetLevel());

    float currentCorrection(1.);
    if (lastLevel == mithep::Jet::L3 && currentLevel == mithep::Jet::L2) {
      // special case for L2L3Residual
      // store L3*L2L3 as L3 correction
      if (iC > 1)
        currentCorrection = corrections.at(iC) / corrections.at(iC - 2);
      else
        currentCorrection = corrections.at(iC);

      currentLevel = mithep::Jet::L3;
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

  if (fSigma != 0.) {
    jet.SetCorrectionScale(UncertaintyFactor(jet), mithep::Jet::Custom);
    jet.EnableCorrection(mithep::Jet::Custom);
  }
}

void
mithep::JetCorrector::Smear(mithep::Jet& jet, Double_t rho, mithep::GenJetCol const* genJets/* = 0*/)
{
  Corrector const* ptResolution = 0;
  Corrector const* phiResolution = 0;
  Corrector* scaleFactor = 0;
  for (auto& corr : fCorrectors) {
    if (corr.GetType() ==  Corrector::kPtResolution)
      ptResolution = &corr;
    if (corr.GetType() ==  Corrector::kPhiResolution)
      phiResolution = &corr;
    if (corr.GetType() ==  Corrector::kScaleFactor)
      scaleFactor = &corr;
  }
  if (ptResolution == 0 || scaleFactor == 0) {
    std::cerr << "Jet smearing requires pt resolution and scale factor input" << std::endl;
    throw std::runtime_error("Configuration");
  }

  auto&& rawMom(jet.RawMom());

  //compute correction factors
  double vars[Corrector::nVarTypes];
  vars[Corrector::kJetEta] = rawMom.Eta();
  vars[Corrector::kJetPt] = rawMom.Pt();
  vars[Corrector::kRho] = rho;

  double ptres = ptResolution->Eval(vars);

  scaleFactor->SetUncertainty(0);
  double sfNominal = scaleFactor->Eval(vars);
  double sf = 1.;
  if (fSigma == 0.) {
    sf = sfNominal;
  }
  else {
    if (fSigma > 0.)
      scaleFactor->SetUncertainty(1);
    else
      scaleFactor->SetUncertainty(-1);

    sf = sfNominal + (scaleFactor->Eval(vars) - sfNominal) * std::abs(fSigma);
  }

  double newPt = jet.Pt();
  double newPhi = jet.Phi();

  mithep::GenJet const* genJet = 0;

  if (genJets) {
    for (unsigned iJ = 0; iJ != genJets->GetEntries(); ++iJ) {
      genJet = genJets->At(iJ);
      if (mithep::MathUtils::DeltaR(*genJet, jet) < fGenJetMatchRadius && std::abs(genJet->Pt() - jet.Pt()) < 3. * ptres)
        break;
    }
  }

  if (genJet) {
    // just enhance the gen/reco discrepancy
    newPt = std::max(0., genJet->Pt() + (jet.Pt() - genJet->Pt()) * sf);
    newPhi = TVector2::Phi_mpi_pi(genJet->Phi() + (jet.Phi() - genJet->Phi()) * sf);
  }
  else {
    // apply additional smearing
    newPt = std::max(0., fRandom.Gaus(jet.Pt(), std::sqrt(sf * sf - 1.) * ptres));
    if (phiResolution) {
      double phires = phiResolution->Eval(vars);
      newPhi = TVector2::Phi_mpi_pi(fRandom.Gaus(jet.Phi(), std::sqrt(sf * sf - 1.) * phires));
    }
  }

  jet.SetRawPtEtaPhiM(newPt, jet.Eta(), newPhi, jet.Mass());
}
