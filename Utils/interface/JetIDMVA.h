//--------------------------------------------------------------------------------------------------
// JetIDMVA
//
// Helper Class for Jet Id MVA
//
// Authors: P. Harris, Y.Iiyama
//--------------------------------------------------------------------------------------------------
#ifndef MITPHYSICS_UTILS_JetIDMVA_H
#define MITPHYSICS_UTILS_JetIDMVA_H

#include "MitAna/DataTree/interface/PFJetFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataCont/interface/Types.h"

namespace TMVA {
  class Reader;
}

namespace mithep {

  class JetIDMVA {
  public:
    enum MVAType {
      kBaseline,
      k42,
      k52,
      kCut,
      kQGP,
      k53,
      k53CHS,
      k53BDTCHSFullPlusRMS,
      k53MET,
      k53METFull,
      k74CHS,
      k76CHS,
      k80CHS,
      k81CHS,
      nMVATypes
    };

    enum CutType {
      kTight,
      kMedium,
      kLoose,
      kMET,
      nCutTypes
    };

    enum Variable {
      kNvtx,
      kJetPt,
      kJetEta,
      kJetPhi,
      kD0,
      kDZ,
      kBeta,
      kBetaStar,
      kNCharged,
      kNNeutrals,
      kDRMean,
      kPtD,
      kFrac01,
      kFrac02,
      kFrac03,
      kFrac04,
      kFrac05,
      kDR2Mean,
      kDRWeighted,
      kRho,
      kNParticles,
      kNCh,
      kAxisMajor,
      kAxisMinor,
      kFRing0,
      kFRing1,
      kFRing2,
      kFRing3,
      kPull,
      kMinPull01,
      kJetR,
      kJetRchg,
      kNTrueInt,
      kDRMatch,
      nVariables
    };

    JetIDMVA();
    virtual ~JetIDMVA();

    static BitMask8 fgCorrectionMask;

    void Initialize(JetIDMVA::CutType, JetIDMVA::MVAType, TString const& weightsConfig, TString const& cutConfig);
    void Initialize(JetIDMVA::CutType, JetIDMVA::MVAType, std::vector<TString> const& weightsConfig, TString const& cutConfig);

    // obsolete
    void Initialize(JetIDMVA::CutType,
                    TString const& iLowPtWeights = "",
                    TString const& iHighPtWeights = "",
                    JetIDMVA::MVAType = kBaseline,
                    TString const& iCutFileName = "");

    Bool_t IsInitialized() const { return fIsInitialized; }

    void SetReproducePullBug(Bool_t b) { fReproducePullBug = b; }
    void SetReproduceCovarianceBug(Bool_t b) { fReproduceCovarianceBug = b; }

    //Cut Based
    Bool_t passCut(PFJet const*, Vertex const*, VertexCol const*);

    //Corrected Jets
    Bool_t pass(PFJet const*, Vertex const*, VertexCol const*, Double_t rho = 0.);

    //What is this function? (Y.I. 2015.07.01)
    Double_t* QGValue(const PFJet *iJet,const Vertex *iVertex,const VertexCol *iVertices, //Vertex here is the PV
                      const PileupEnergyDensityCol *iPileupEnergyDensity,
                      Bool_t printDebug);

    //Corrected Jets
    Double_t MVAValue(const PFJet *iJet,const Vertex *iVertex,const VertexCol *iVertices, Double_t rho,
                      Bool_t printDebug=false);

    Float_t fDZCut = 0.2;             // dZ cut used in beta and beta* calculation to define association to PV
                                      // fDZCut <= 0. switches the beta and beta* definitions to the "classic"
                                      // versions, where vertex-track association in terms of fit weights is used.

  protected:
    Bool_t InitializeCuts(TString const& fileName, TString const& cutId, TString const& cutType);

    std::vector<double> fEtaBinLowEdges; // bin low edges (size fReaders + 1)
    std::vector<TMVA::Reader*> fReaders;
    TString       fMethodName = "JetIDMVAHighPt";
    MVAType       fType = nMVATypes;
    Bool_t        fIsInitialized = kFALSE;
    Float_t       fMVACut[4][4]{}; //Fix the cut array
    Float_t       fRMSCut[4][4]{};
    Float_t       fBetaStarCut[4][4]{};

    Float_t       fVariables[nVariables]{};
    Bool_t        fVariableUsed[nVariables]{};
    TString       fVarNames[nVariables]; // initialized in the Ctor

    Bool_t        fReproducePullBug = kFALSE; // CMSSW <= 7_6 had a bug in pull computation
    Bool_t        fReproduceCovarianceBug = kFALSE; // CMSSW <= 8_0 had a bug in covariance computation

    ClassDef(JetIDMVA,0)
  };
}


#endif
