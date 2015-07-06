//--------------------------------------------------------------------------------------------------
// MVATools
//
// Helper Class for photon Identification decisions.
//
// Authors: M.Yang 2011/10/12
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_MVATOOLS_H
#define MITPHYSICS_UTILS_MVATOOLS_H

#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/DecayParticleFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/BeamSpotFwd.h"
#include "MitAna/DataTree/interface/TriggerObjectFwd.h"
#include "MitAna/DataTree/interface/MCParticleFwd.h"
#include "MitAna/DataTree/interface/SuperClusterFwd.h"
#include "MitAna/DataTree/interface/PFCandidateFwd.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

class TRandom3;
namespace TMVA {//MVA
  class Reader;
}

namespace mithep {
  class MVATools {
  public:
    MVATools();
    virtual ~MVATools() {}
    
    // MVA Typles
    enum IdMVAType {
      kNone      = 0,
      k2011IdMVA,
      k2012IdMVA_globe,
      k2012IdMVA,        // same as above, but more logical name...
      k2013FinalIdMVA_8TeV, 
      k2013FinalIdMVA_7TeV, 
      k2011IdMVA_HZg
    };
    
    //--------------------------
    //MVA
    //--------------------------
    
    //void InitializeMVA(int VariableType, TString EndcapWeights,TString BarrelWeights);
    void InitializeMVA(IdMVAType type);

    // removed this (fab): if a mod want to cuty on BDT values, compute them in the mod using GetMVAbdtValue(..) and make the cut in the mod
    //Bool_t PassMVASelection   (const Photon* p,const Vertex* vtx,const TrackCol* trackCol,const VertexCol* vtxCol,Double_t _tRho,Float_t bdtCutBarrel, Float_t bdtCutEndcap, const ElectronCol* els=0, Bool_t applyElectronVeto=kTRUE);
    // removed this (fab): no needed naywhere...
    //Int_t  PassElectronVetoInt(const Photon* p, const ElectronCol* els);
    
    Double_t GetMVAbdtValue   (const Photon* p,const Vertex* vtx,const TrackCol* trackCol,const VertexCol* vtxCol,Double_t _tRho, const PFCandidateCol *fPFCands=NULL ,const ElectronCol* els=0, Bool_t applyElectronVeto=kTRUE);
    
    // these we can remove at some point
/*     Float_t GetMVAbdtValue_2011      (const Photon* p,const Vertex* vtx,const TrackCol* trackCol,const VertexCol* vtxCol,Double_t _tRho,                                      const ElectronCol* els=0, Bool_t applyElectronVeto=kTRUE); */
/*     Float_t GetMVAbdtValue_2012_globe(const Photon* p,const Vertex* vtx,const TrackCol* trackCol,const VertexCol* vtxCol,Double_t _tRho, const PFCandidateCol *fPFCands=NULL ,const ElectronCol* els=0, Bool_t applyElectronVeto=kTRUE); */
    
    
  private:
    
    TMVA::Reader *fReaderEndcap;
    TMVA::Reader *fReaderBarrel;    
    
    IdMVAType fMVAType;
    std::vector<float> mvaVars;
    std::vector<std::string> varNames;
    std::map<std::string,unsigned int> mvaVarMapEB;
    std::map<std::string,unsigned int> mvaVarMapEE;
    
    // -------------------------------------------------
    // fab: these guys should all go away...
    //MVA Variables
/*     float HoE; */
/*     float covIEtaIEta; */
/*     float tIso1abs; */
/*     float tIso3abs; */
/*     float tIso2abs; */
/*     float R9; */
    
/*     float absIsoEcal; */
/*     float absIsoHcal; */
/*     float RelEMax; */
/*     float RelETop; */
/*     float RelEBottom; */
/*     float RelELeft; */
/*     float RelERight; */
/*     float RelE2x5Max; */
/*     float RelE2x5Top; */
/*     float RelE2x5Bottom; */
/*     float RelE2x5Left; */
/*     float RelE2x5Right; */
/*     float RelE5x5; */
    
/*     float EtaWidth; */
/*     float PhiWidth; */
/*     float CoviEtaiPhi; */
/*     float CoviPhiiPhi; */
    
/*     float NVertexes; */
/*     float RelPreshowerEnergy; */
    
/*     //variable for v2 and v1 */
/*     float RelIsoEcal; */
/*     float RelIsoHcal; */
/*     float tIso1; */
/*     float tIso3; */
/*     float tIso2; */
    
/*     float ScEta; */
    
/*     //variables used to compute mva variables */
    
/*     Bool_t PassElecVeto;   */
    
/*     double ecalIso3; */
/*     double ecalIso4; */
/*     double hcalIso4;  */
    
/*     unsigned int wVtxInd; */
    
/*     double trackIso1; */
    
/*     // track iso only */
/*     double trackIso3; */
    
/*     // track iso worst vtx */
/*     double trackIso2;  */
    
/*     double combIso1; */
/*     double combIso2; */
    
/*     double RawEnergy; */
    
/*     double dRTrack; */
    
/*     //spectator variables */
/*     double Pt_MVA; */
/*     double ScEta_MVA; */
    
/*     Bool_t isbarrel; */
    
/*     // check which category it is ... */
/*     int _tCat; */
    
/*     //1201 variable */
/*     float myphoton_pfchargedisogood03; */
/*     float myphoton_pfchargedisobad03;  */
/*     float myphoton_pfphotoniso03; */
/*     float myphoton_sieie;  */
/*     float myphoton_sieip; */
/*     float myphoton_etawidth;  */
/*     float myphoton_phiwidth; */
/*     float myphoton_r9;  */
/*     float myphoton_s4ratio; */
/*     float myphoton_SCeta;  */
/*     float event_rho; */
/*     float myphoton_ESEffSigmaRR; */
    
    //MVA  
/*     Bool_t PassMVA; */
/*     //Float_t bdt;   -> removed this, we should not have this a memeber variable !!! (fab) */
/*     Int_t PassElecVetoInt; */


    TMVA::Reader *reader;    
  };
}

#endif
