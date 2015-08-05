//--------------------------------------------------------------------------------------------------
// JetCorrector
//
// A FactorizedJetCorrector wrapper. Used by JetCorrectionMod and MetCorrectionMod (and possibly
// by other modules)
//
// Authors: Y.Iiyama
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_JETCORRECTOR_H
#define MITPHYSICS_UTILS_JETCORRECTOR_H

#include "MitAna/DataTree/interface/Jet.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

#include <vector>
#include <algorithm>

class FactorizedJetCorrector;
class JetCorrectionUncertainty;

namespace mithep {

  class JetCorrector {
  public:
    JetCorrector() {}
    virtual ~JetCorrector();

    void AddParameterFile(char const* fileName);
    void ClearParameters();
    void SetMaxCorrLevel(mithep::Jet::ECorr m) { fMaxCorrLevel = m; }
    void SetUncertaintySigma(Double_t s) { fSigma = s; }
    void Initialize();

    std::vector<mithep::Jet::ECorr> const& GetLevels() const { return fLevels; }
    mithep::Jet::ECorr GetMaxCorrLevel() const { return fMaxCorrLevel; }
    // Returns a vector of cumulative corrections, e.g. {L1, L1*L2, L1*L2*L3, ...}
    std::vector<Float_t> CorrectionFactors(mithep::Jet&, Double_t rho = 0.) const;
    // Returns the full correction (i.e. the last element of CorrectionFactors)
    Float_t CorrectionFactor(mithep::Jet&, Double_t rho = 0.) const;
    Float_t UncertaintyFactor(mithep::Jet&) const;
    Bool_t IsInitialized() const;
    Bool_t IsEnabled(mithep::Jet::ECorr l) const
    { return l < fMaxCorrLevel && std::find(fLevels.begin(), fLevels.end(), l) != fLevels.end(); }

    void Correct(mithep::Jet&, Double_t rho = 0.) const;

    static mithep::Jet::ECorr TranslateLevel(char const* levelName);

  private:
    std::vector<JetCorrectorParameters> fParameters{};
    std::vector<mithep::Jet::ECorr> fLevels{};
    mithep::Jet::ECorr fMaxCorrLevel = mithep::Jet::nECorrs;
    FactorizedJetCorrector* fCorrector = 0;
    JetCorrectionUncertainty* fUncertainty = 0;
    Double_t fSigma = 0.;

    ClassDef(JetCorrector, 0)
  };

}

#endif
