//--------------------------------------------------------------------------------------------------
// Recoil Tools
//
// Helper Class for Recoil Tools
//
// Authors: P. Harris
//--------------------------------------------------------------------------------------------------
#ifndef MITPHYSICS_UTILS_RecoilTools_H
#define MITPHYSICS_UTILS_RecoilTools_H

#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/PFJetFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitPhysics/Utils/interface/JetIDMVA.h"

namespace mithep {
  class RecoilTools {
    public:
    RecoilTools(TString iJetLowPtMVAFile = "",
		TString iJetHighPtMVAFile = "",
		TString iCutFile = "",
		bool i42=false,JetIDMVA::MVAType iType=JetIDMVA::kBaseline,bool iUseRho=true);
    virtual ~RecoilTools();
    JetIDMVA *fJetIDMVA;

    Met pfRecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,const PFCandidateCol *iCands);

    //Candidate filtered
    Met pfRecoil(double iPhi1,double iEta1,double iPhi2,double iEta2,const PFCandidateCol *iCands);
    Met pfCone  (double iPhi1,double iEta1,const PFCandidateCol *iCands,const Vertex *iVertex,
		 bool iCharge=false,Double_t iDZCut=0.3);

    Met trackMet(const PFCandidateCol *iCands,const Vertex *iVertex,Double_t iDZCut=0.1);

    //Candidate filtered
    Met trackRecoil(double iPhi1,double iEta1,double iPhi2,double iEta2,
		    const PFCandidateCol *iCands,const Vertex *iVertex,Double_t iDZCut=0.1);

    Met trackRecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
			     const PFCandidateCol *iCands,const Vertex *iVertex,double iDZCut=0.1);

    bool filter (const PFJet *iJet,Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2);

    bool filter(const PFCandidate *iCand,Double_t iPhi1,Double_t iEta1,Double_t iPhi2,Double_t iEta2);

    //Corrected Jets
    void addNeut(const PFJet *iJet,FourVectorM &iVec,Double_t &iSumEt,Double_t iRho,
		 int iSign=1);
    
    //Corrected Jets
    Met NoPUMet( const PFJetCol       *iJets,
		 const PFCandidateCol *iCands,const Vertex *iVertex,const VertexCol *iVertices,Double_t iRho,
		 Double_t iPhi1=1000,Double_t iEta1=1000,Double_t iPhi2=1000,Double_t iEta2=1000,
		 Double_t iDZCut=0.1);

    //Corrrected Jets
    Met NoPURecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,   
		   const PFJetCol       *iJets,
		   const PFCandidateCol *iCands,const Vertex *iVertex,const VertexCol *iVertices,Double_t iRho,
		   Double_t iPhi1=1000,Double_t iEta1=1000,Double_t iPhi2=1000,Double_t iEta2=1000,
		   Double_t iDZCut=0.1);

    //Corrected Jets
    Met PUCMet( const PFJetCol       *iJets,
		const PFCandidateCol *iCands,const Vertex *iVertex,const VertexCol *iVertices,Double_t iRho,bool iUseType1=false,
		Double_t iPhi1=1000,Double_t iEta1=1000,Double_t iPhi2=1000,Double_t iEta2=1000,
		Double_t iDZCut=0.1);

    //Corrected Jets
    Met PUCRecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
		  const PFJetCol       *iJets,
		  const PFCandidateCol *iCands,const Vertex *iVertex,const VertexCol *iVertices,Double_t iRho,bool iUseType1=false,
		  Double_t iPhi1=1000,Double_t iEta1=1000,Double_t iPhi2=1000,Double_t iEta2=1000,
		  Double_t iDZCut=0.1);

    //Corrected Jets
    Met PUMet( const PFJetCol       *iJets,
	       const PFCandidateCol *iCands,const Vertex *iVertex,const VertexCol *iVertices,Double_t iRho,
	       Double_t iPhi1=1000,Double_t iEta1=1000,Double_t iPhi2=1000,Double_t iEta2=1000,
	       Double_t iDZCut=0.2);

    //Corrected Jets
    Met PURecoil(Double_t iVisPt,Double_t iVisPhi,Double_t iVisSumEt,
		 const PFJetCol       *iJets,
		 const PFCandidateCol *iCands,const Vertex *iVertex,const VertexCol *iVertices,Double_t iRho,
		 Double_t iPhi1=1000,Double_t iEta1=1000,Double_t iPhi2=1000,Double_t iEta2=1000,
		 Double_t iDZCut=0.1);
    
    bool f42;
    bool fUseRho;
    ClassDef(RecoilTools, 1) // Recoil tools
  };
}
#endif
