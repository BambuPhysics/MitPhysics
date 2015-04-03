//--------------------------------------------------------------------------------------------------
// $Id $
//
// VertexTools
//
// Helper Class for photon Identification decisions.
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_VERTEXTOOLS_H
#define MITPHYSICS_UTILS_VERTEXTOOLS_H

#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/Electron.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/BaseVertex.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/VertexMVA.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"

#include "TMVA/Reader.h"

namespace TMVA {//MVA
  class Reader;
}

namespace mithep {
  typedef std::vector<double> VertexZarray; 
  typedef std::vector<const Track*> TrackArray;

  class VertexTools {
  public:

    VertexTools();
    virtual ~VertexTools() {}
    
    static double NewMass(const Photon* ph1, const Photon* ph2, const BaseVertex* vert);

    static VertexZarray ExtractZarray(const VertexCol* vcol, float zmin=-30, float zmax = 30,
				      const BaseVertex  *fBeamSpot = NULL);
    static VertexZarray ExtractZarray(float zmin=-30, float zmax=30, float step=0.05);
    
    static const Vertex* BestVtx(const PFCandidateCol *fPFJets, const VertexCol *c,
				 const BaseVertex  *fBeamSpot, FourVector diboso);


    static double BestVtx(const PFCandidateCol *fPFJets, VertexZarray zcol, 
			  const BaseVertex  *fBeamSpot, FourVector diboso);

    static double Prob(const PFCandidateCol *fPFJets, double zpos,  
		       const BaseVertex  *fBeamSpot, FourVector diboso);
    
    static double VertexWidth(const Vertex*,  const BaseVertex* );

    static void BanThisTrack(const Track*);
    static void Reset();
    
    static VertexTools* instance(const char* str){
      if(meobject == NULL){
	meobject = new VertexTools();
        meobject->InitM(str);
      }
      return meobject;
    }        

    // ----------------------------------------------------------
    // Methods (added by Fabian) on the EPS BaseLine Analysis
    const Vertex* findVtxBasicRanking(const Photon*           ph1, 
				      const Photon*           ph2, 
				      const BaseVertex*       bsp,
				      const VertexCol*        vtcs,
				      const DecayParticleCol* conv, Bool_t useMva, Double_t &vtxProb,
				      std::vector<int>   * debugInds = NULL ,
				      std::vector<double>* debugVals = NULL ,
				      std::vector<int>   * debugConv = NULL );
    // ----------------------------------------------------------

    
    static std::pair<double,double> VtxZFromConversion(const Photon *p, const DecayParticle *c, const BaseVertex *bsp);

    void InitM(const char* str);
    void InitP(int version = 1);
    
    Bool_t IsInitMvaM() const { return fIsInitMvaM; }
    Bool_t IsInitMvaP() const { return fIsInitMvaP; }
    
    static Double_t DeltaMassVtx(Double_t xp1, Double_t yp1, Double_t zp1,
            Double_t xp2, Double_t yp2, Double_t zp2,
	    Double_t xv,  Double_t yv,  Double_t zv,
            Double_t dz);
    
  private:
    
    double VtxMvaP(float ptbal, float ptasym, float logsumpt2, float limPullToConv, float nConv) const;
    
    static VertexTools *meobject;
    

    TString relname;

    TrackArray excluded;
    
    Bool_t fIsInitMvaM;
    Bool_t fIsInitMvaP;

    Float_t tmvar1, tmvar2, tmvar3, tmvar4, tmvar5, tmvar6;
    TMVA::Reader* reader;
    
    
    TMVA::Reader *readervtx;
    TMVA::Reader *readerevt;
    mutable Float_t fMvaPVars[5];
    mutable Float_t fMvaPEvtVars[8];
    
    TString fVtxProbMvaName;
    
    ClassDef(VertexTools, 0) // Muon tools
      };
}

#endif
