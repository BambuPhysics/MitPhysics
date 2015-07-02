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
      k53MET,
      k53METFull,
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
      kNVtx,
      kJPt1,
      kJEta1,
      kJPhi1,
      kJD01,
      kJDZ1,
      kBeta,
      kBetaStar,
      kNCharged,
      kNNeutrals,
      kNParticles,
      kDRMean,
      kPtD,
      kFrac01,
      kFrac02,
      kFrac03,
      kFrac04,
      kFrac05,
      kDR2Mean,
      nVariables
    };

    JetIDMVA() {}
    virtual ~JetIDMVA();

    void Initialize(JetIDMVA::CutType, JetIDMVA::MVAType, TString const& weightsConfig, TString const& cutConfig);

    // obsolete
    void Initialize(JetIDMVA::CutType,
                    TString const& iLowPtWeights = "",
                    TString const& iHighPtWeights = "",
                    JetIDMVA::MVAType = kBaseline,
                    TString const& iCutFileName = "");

    Bool_t IsInitialized() const { return fIsInitialized; }

    //Cut Based
    Bool_t passCut(PFJet const*, Vertex const*, VertexCol const*);

    //Corrected Jets
    Bool_t pass(PFJet const*, Vertex const*, VertexCol const*);

    //What is this function? (Y.I. 2015.07.01)
    Double_t* QGValue(const PFJet *iJet,const Vertex *iVertex,const VertexCol *iVertices, //Vertex here is the PV
                      const PileupEnergyDensityCol *iPileupEnergyDensity,
                      Bool_t printDebug);

    //Corrected Jets
    Double_t MVAValue(const PFJet *iJet,const Vertex *iVertex,const VertexCol *iVertices,
                      Bool_t printDebug=false);

    Float_t fJetPtMin = 0.;
    Float_t fDZCut = 0.2;

  protected:
    TMVA::Reader* fReader = 0;
    TString       fMethodName = "JetIDMVAHighPt";
    MVAType       fType = nMVATypes;
    Bool_t        fIsInitialized = kFALSE;
    Float_t       fMVACut[4][4] = {}; //Fix the cut array
    Float_t       fRMSCut[4][4] = {};
    Float_t       fBetaStarCut[4][4] = {};

    Float_t      fVariables[nVariables] = {};

    ClassDef(JetIDMVA,0)
  };
}


#endif
