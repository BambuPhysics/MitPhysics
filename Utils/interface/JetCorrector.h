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

#include "TFormula.h"

#include <vector>
#include <algorithm>

namespace mithep {

  class JetCorrector {
  public:
    class Corrector {
      // copied from CondFormats/JetMETObjects/interface/JetCorrectorParameters.h (7_4_X)
    public:
      enum VarType {
        kJetEta,
        kNPV,
        kJetPt,
        kJetPhi,
        kJetE,
        kJetEMF,
        kJetA,
        kRho,
        kRelLepPt,
        kPtRel,
        nVarTypes
      };

      struct Record {
        Record(unsigned nBinVars, unsigned nFormVars, unsigned nParameters, std::string const& inputLine);

        std::vector<std::pair<float, float>> fBinVarLimits;
        std::vector<std::pair<float, float>> fFormVarLimits;
        std::vector<float> fParameters;
      };
      
      Corrector(char const* fileName);
      ~Corrector() {}

      Jet::ECorr GetLevel() const { return fLevel; }
      Double_t Eval(Double_t [nVarTypes]) const;
      Bool_t GetEnabled() const { return fEnabled; }
      Bool_t GetInterpolate() const { return fInterpolate; }

      void SetEnabled(Bool_t b = kTRUE) { fEnabled = b; }
      void SetInterpolate(Bool_t b = kTRUE) { fInterpolate = b; }
      void SetUncertaintyUp(Bool_t b = kTRUE) { fUncertaintyUp = b; }
      
      static mithep::Jet::ECorr ParseLevelName(char const*);
      static VarType ParseVariableName(char const*);

    private:
      Double_t EvalCorrection(Double_t [nVarTypes]) const;
      Double_t EvalUncertainty(Double_t [nVarTypes]) const;
      UInt_t FindRecord(Double_t [nVarTypes]) const;

      Jet::ECorr fLevel{Jet::nECorrs};
      std::vector<VarType> fBinningVariables{};
      std::vector<VarType> fFormulaVariables{};
      std::vector<Record> fRecords{};
      Bool_t fIsResponse{kFALSE};
      TFormula* fFormula{0};
      Bool_t fEnabled{kTRUE};
      Bool_t fInterpolate{kFALSE}; // turned off in CMSSW
      Bool_t fUncertaintyUp{kTRUE};
    };
    
    JetCorrector();
    virtual ~JetCorrector() {}

    void AddParameterFile(char const* fileName);
    void ClearParameters();
    void SetMaxCorrLevel(mithep::Jet::ECorr m) { fMaxCorrLevel = m; }
    void SetUncertaintySigma(Double_t s) { fSigma = s; }
    void Initialize();

    mithep::Jet::ECorr GetMaxCorrLevel() const { return fMaxCorrLevel; }
    // Returns a vector of cumulative corrections, e.g. {L1, L1*L2, L1*L2*L3, ...}
    std::vector<Float_t> CorrectionFactors(mithep::Jet const&, Double_t rho = 0.) const;
    // Returns the full correction (i.e. the last element of CorrectionFactors)
    Float_t CorrectionFactor(mithep::Jet const&, Double_t rho = 0.) const;
    Float_t UncertaintyFactor(mithep::Jet const&) const;
    Bool_t IsEnabled(mithep::Jet::ECorr l) const;

    void Correct(mithep::Jet&, Double_t rho = 0.) const;

  private:
    std::vector<Corrector> fCorrectors{};
    mithep::Jet::ECorr fMaxCorrLevel = mithep::Jet::nECorrs;
    Double_t fSigma = 0.; // uncertainty

    ClassDef(JetCorrector, 0)
  };

}

#endif
