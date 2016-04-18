//--------------------------------------------------------------------------------------------------
// JetCorrector
//
// An original implementation of CMSSW FactorizedJetCorrector equivalent. Used by JetCorrectionMod
// and MetCorrectionMod (and possibly by other modules)
//
// Authors: Y.Iiyama
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_JETCORRECTOR_H
#define MITPHYSICS_UTILS_JETCORRECTOR_H

#include "MitAna/DataTree/interface/Jet.h"
#include "MitAna/DataTree/interface/GenJetCol.h"

#include "TFormula.h"
#include "TRandom3.h"

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

      enum FactorType {
        kCorrection,
        kResponse,
        kResolution, // dummy type not used
        kPtResolution,
        kPhiResolution,
        kScaleFactor,
        nFactorTypes
      };

      struct Record {
        Record(unsigned nBinVars, unsigned nFormVars, std::string const& inputLine);

        std::vector<std::pair<float, float>> fBinVarLimits;
        std::vector<std::pair<float, float>> fFormVarLimits;
        std::vector<float> fParameters;
      };
      
      Corrector(char const* fileName, FactorType factorType = nFactorTypes);
      ~Corrector() {}

      Jet::ECorr GetLevel() const { return fLevel; }
      Double_t Eval(Double_t [nVarTypes]) const;
      Bool_t GetEnabled() const { return fEnabled; }
      Bool_t GetInterpolate() const { return fInterpolate; }
      FactorType GetType() const { return fType; }

      void SetEnabled(Bool_t b = kTRUE) { fEnabled = b; }
      void SetInterpolate(Bool_t b = kTRUE) { fInterpolate = b; }
      void SetUncertainty(Int_t v) { fUncertaintyV = v; } // 0->nominal, -1->down, 1->up
      
      static mithep::Jet::ECorr ParseLevelName(char const*);
      static VarType ParseVariableName(char const*);
      static FactorType ParseTypeName(char const*);

    private:
      Double_t EvalCorrection(Double_t [nVarTypes]) const;
      Double_t EvalUncertainty(Double_t [nVarTypes]) const;
      Double_t EvalScaleFactor(Double_t [nVarTypes]) const;
      UInt_t FindRecord(Double_t [nVarTypes]) const;

      Jet::ECorr fLevel{Jet::nECorrs};
      std::vector<VarType> fBinningVariables{};
      std::vector<VarType> fFormulaVariables{};
      std::vector<Record> fRecords{};
      FactorType fType{nFactorTypes};
      TFormula* fFormula{0};
      Bool_t fEnabled{kTRUE};
      Bool_t fInterpolate{kFALSE}; // turned off in CMSSW
      Int_t fUncertaintyV{0}; // used in JEC uncertainty and JER scale factor
    };
    
    JetCorrector();
    virtual ~JetCorrector() {}

    // Need to specify the factor type for PtResolution vs PhiResolution
    // Former is returned as a relative resolution while the latter is absolute.
    void AddParameterFile(char const* fileName, Corrector::FactorType factorType = Corrector::nFactorTypes);
    void ClearParameters();
    void SetMaxCorrLevel(mithep::Jet::ECorr m) { fMaxCorrLevel = m; }
    void SetUncertaintySigma(Double_t s) { fSigma = s; }
    void SetGenJetMatchRadius(Double_t r) { fGenJetMatchRadius = r; }
    void SetSmearingSeed(Int_t s) { fRandom.SetSeed(s); }

    mithep::Jet::ECorr GetMaxCorrLevel() const { return fMaxCorrLevel; }
    // Returns a vector of cumulative corrections, e.g. {L1, L1*L2, L1*L2*L3, ...}
    std::vector<Float_t> CorrectionFactors(mithep::Jet const&, Double_t rho = 0.) const;
    // Returns the full correction (i.e. the last element of CorrectionFactors)
    Float_t CorrectionFactor(mithep::Jet const&, Double_t rho = 0.) const;
    Float_t UncertaintyFactor(mithep::Jet const&) const;
    Bool_t IsEnabled(mithep::Jet::ECorr l) const;
    Bool_t HasSmearing() const { return fHasSmearing; }

    // Correct the given jet with JEC parameters.
    void Correct(mithep::Jet&, Double_t rho = 0.) const;
    // Smear the given jet with JER parameters. If gen jet collection is passed, try match & scale first.
    void Smear(mithep::Jet&, Double_t rho, mithep::GenJetCol const* = 0); // cannot be const because of TRandom

  private:
    std::vector<Corrector> fCorrectors{};
    mithep::Jet::ECorr fMaxCorrLevel = mithep::Jet::nECorrs;
    Bool_t fHasSmearing = kFALSE;
    Double_t fSigma = 0.; // uncertainty
    Double_t fGenJetMatchRadius = 0.2; // should be 1/2 jet cone
    TRandom3 fRandom;

    ClassDef(JetCorrector, 0)
  };

}

#endif
