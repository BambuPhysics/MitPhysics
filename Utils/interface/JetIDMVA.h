//--------------------------------------------------------------------------------------------------
// JetIDMVA
//
// Helper Class for Jet Id MVA
//
// Authors: P. Harris
//--------------------------------------------------------------------------------------------------
#ifndef MITPHYSICS_UTILS_JetIDMVA_H
#define MITPHYSICS_UTILS_JetIDMVA_H

#include "MitAna/DataTree/interface/PFJetFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

// for Rho definitons
#include "MitPhysics/Utils/interface/RhoUtilities.h"

class FactorizedJetCorrector;

class TRandom3;
namespace TMVA {
  class Reader;
}

namespace mithep {
  class JetIDMVA {
    public:
      JetIDMVA();
      virtual ~JetIDMVA(); 

      enum MVAType {
        kBaseline  = 0,
	k42        = 1,
	k52        = 2,
	kCut       = 3,
	kQGP       = 4,
	k53        = 5,
	k53MET     = 6,
	k53METFull = 7
      };

      enum CutType {
        kTight     = 0,
        kMedium    = 1,
        kLoose     = 2,
        kMET       = 3
      };

      void     Initialize(JetIDMVA::CutType iCutType,
			  TString           iLowPtWeights ="$MIT_DATA/mva_JetID_lowpt.weights.xml",
			  TString           iHighPtWeights="$MIT_DATA/mva_JetID_highpt.weights.xml",
			  JetIDMVA::MVAType iType=kBaseline,
			  TString           iCutFileName  ="$CMSSW_BASE/src/MitPhysics/Utils/python/JetIdParams_cfi.py",bool i42=false);
    
      Bool_t   IsInitialized() const { return fIsInitialized; }
      Double_t MVAValue(    
			Float_t iNPV    ,
			Float_t iJPt1   ,
			Float_t iJEta1  ,
			Float_t iJPhi1  ,
			Float_t iJD01   ,
			Float_t iJDZ1   ,
			Float_t iBeta   ,
			Float_t iBetaStar,
			Float_t iNCharged,
			Float_t iNNeutrals,
			Float_t iDRMean  ,
			Float_t iFrac01  ,
			Float_t iFrac02  ,
			Float_t iFrac03  ,
			Float_t iFrac04  ,
			Float_t iFrac05  ,
			Float_t iDR2Mean ,
			Float_t iPtD
			);
      //Cut Based
      Bool_t passCut(const PFJet *iJet,const Vertex *iVertex,const VertexCol *iVertices);

      Bool_t passPt(const PFJet *iJet, FactorizedJetCorrector *iJetCorrector=0,
		    const PileupEnergyDensityCol *iPileupEnergyDensity=0,
		    RhoUtilities::RhoType type=RhoUtilities::DEFAULT);
      //UNcorrected Jets
      Bool_t   pass(const PFJet *iJet,const Vertex *iVertex,const VertexCol *iVertices,
		    FactorizedJetCorrector *iJetCorrector,
		    const PileupEnergyDensityCol *iPileupEnergyDensity,
		    RhoUtilities::RhoType type = RhoUtilities::DEFAULT);

      //Corrected Jets
      Bool_t   pass(const PFJet *iJet,const Vertex *iVertex,const VertexCol *iVertices);
		    			
      //Uncorrected Jets
      Double_t MVAValue(const PFJet *iJet,const Vertex *iVertex,const VertexCol *iVertices,
			FactorizedJetCorrector *iJetCorrector,
			const PileupEnergyDensityCol *iPileupEnergyDensity,
			Bool_t printDebug=false);

      Double_t* QGValue(const PFJet *iJet,const Vertex *iVertex,const VertexCol *iVertices, //Vertex here is the PV
			FactorizedJetCorrector *iJetCorrector,
			const PileupEnergyDensityCol *iPileupEnergyDensity,
			Bool_t printDebug);

      //Corrected Jets
      Double_t MVAValue(const PFJet *iJet,const Vertex *iVertex,const VertexCol *iVertices,
			Bool_t printDebug=false);


      double  correctedPt(const PFJet *iJet, FactorizedJetCorrector *iJetCorrector,
			  const PileupEnergyDensityCol *iPUEnergyDensity,
		          RhoUtilities::RhoType type = RhoUtilities::DEFAULT,int iId=-1);

      Float_t                  fJetPtMin;
      Float_t                  fDZCut;
      Bool_t                   f42;

    protected:      
      TMVA::Reader            *fReader;
      TMVA::Reader            *fLowPtReader;
      TString                  fLowPtMethodName;
      TString                  fHighPtMethodName;
      MVAType                  fType;
      CutType                  fCutType;
      Bool_t                   fIsInitialized;
      Float_t                  fMVACut[4][4]; //Fix the cut array
      Float_t                  fRMSCut[4][4];
      Float_t                  fBetaStarCut[4][4];

      Float_t fNVtx     ;
      Float_t fJPt1     ;
      Float_t fJEta1    ;
      Float_t fJPhi1    ;
      Float_t fJD01     ;
      Float_t fJDZ1     ;
      Float_t fBeta     ;
      Float_t fBetaStar ;
      Float_t fNCharged ;
      Float_t fNNeutrals;
      Float_t fNParticles;
      Float_t fDRMean   ;
      Float_t fPtD      ;
      Float_t fFrac01   ;
      Float_t fFrac02   ;
      Float_t fFrac03   ;
      Float_t fFrac04   ;
      Float_t fFrac05   ;
      Float_t fDR2Mean  ;

      ClassDef(JetIDMVA,0)
  };
}


#endif
