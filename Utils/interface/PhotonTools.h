//--------------------------------------------------------------------------------------------------
// $Id $
//
// PhotonTools
//
// Helper Class for photon Identification decisions.
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_PHOTONTOOLS_H
#define MITPHYSICS_UTILS_PHOTONTOOLS_H

#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"
#include "MitAna/DataTree/interface/Electron.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/BaseVertex.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/SuperCluster.h"
#include "MitAna/DataTree/interface/SuperClusterCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include "MitAna/DataTree/interface/PFCandidateCol.h"

class TRandom3;

namespace mithep {
  class PhotonTools {
    public:
      PhotonTools();
  
     enum DiphotonR9EtaCats {
        kCat1 = 0,       //barrel-barrel highr9/highr9
        kCat2,             //barrel-barrel highr9/lowr9+lowr9/lowr9
        kCat3,             //barrel-endcap+endcap/endcap highr9/highr9
        kCat4        //barrel-endcap+endcap-endcap highr9/lowr9 + lowr9-lowr9
      };
      
     enum DiphotonR9EtaPtCats {
       kOctCat0,
       kOctCat1,
       kOctCat2,
       kOctCat3,
       kOctCat4,
       kOctCat5,
       kOctCat6,
       kOctCat7
     };
      
     enum DiphotonR9EtaConversionCats {
        kNewCat1 = 0,       //barrel-barrel highr9/highr9
        kNewCat2,             //barrel-barrel highr9/lowr9+lowr9/lowr9 one/two conversion
        kNewCat3,             //barrel-barrel highr9/lowr9+lowr9/lowr9 no conversion        
        kNewCat4,             //barrel-endcap+endcap/endcap highr9/highr9
        kNewCat5,        //barrel-endcap+endcap-endcap highr9/lowr9 + lowr9-lowr9 one/two conversion
        kNewCat6        //barrel-endcap+endcap-endcap highr9/lowr9 + lowr9-lowr9 no conversion
      };      
      
     enum CiCBaseLineCats {
       kCiCNoCat = 0,
       kCiCCat1,
       kCiCCat2,
       kCiCCat3,
       kCiCCat4
     };     
     
     enum eScaleCats {
       kEBlowEtaGold = 0,
       kEBlowEtaBad,
       kEBhighEtaGold,
       kEBhighEtaBad,
       kEElowEtaGold,
       kEElowEtaBad,
       kEEhighEtaGold,
       kEEhighEtaBad
     };
       
     enum ShowerShapeScales {
       kNoShowerShapeScaling = 0,
       k2011ShowerShape,
       k2012ShowerShape
     };

    static eScaleCats EScaleCat(const Photon *p);

    // Methods for scaling/smearing Photons
    static void ScalePhoton(Photon* p, Double_t scale);
    static void SmearPhoton(Photon* p, TRandom3* rng, Double_t width, UInt_t iSeed);
    static void SmearPhotonError(Photon* p, Double_t width);
    static void ScalePhotonR9(Photon *p, Double_t scale);
    static void ScalePhotonError(Photon *p, Double_t scale);

    static void ScalePhotonShowerShapes(Photon *p, ShowerShapeScales scale);


    static Bool_t       PassSinglePhotonPresel(const Photon *p,const ElectronCol *els, const DecayParticleCol *conversions, const BaseVertex *bs, const TrackCol* trackCol, const Vertex *vtx, double rho, Bool_t applyElectronVeto = kTRUE, Bool_t invertElectronVeto = kFALSE);
    static Bool_t       PassSinglePhotonPreselPFISO(const Photon *p,const ElectronCol *els, const DecayParticleCol *conversions, const BaseVertex *bs, const TrackCol* trackCol,const Vertex *vtx, double rho, const PFCandidateCol *fPFCands, Bool_t applyElectronVeto = kTRUE, Bool_t invertElectronVeto = kFALSE);
    static Bool_t       PassSinglePhotonPreselPFISONoEcal(const Photon *p,const ElectronCol *els, const DecayParticleCol *conversions, const BaseVertex *bs, const TrackCol* trackCol,const Vertex *vtx, double rho, const PFCandidateCol *fPFCands, Bool_t applyElectronVeto = kTRUE, Bool_t invertElectronVeto = kFALSE);
    static Bool_t       PassSinglePhotonPreselPFISONoEcalNoPFChargedIso(const Photon *p,const ElectronCol *els, const DecayParticleCol *conversions, const BaseVertex *bs, const TrackCol* trackCol,const Vertex *vtx, double rho, const PFCandidateCol *fPFCands, Bool_t applyElectronVeto = kTRUE, Bool_t invertElectronVeto = kFALSE);
    static Bool_t       PassSinglePhotonPreselPFISO_NoTrigger(const Photon *p,const ElectronCol *els, const DecayParticleCol *conversions, const BaseVertex *bs, const TrackCol* trackCol,const Vertex *vtx, double rho, const PFCandidateCol *fPFCands, Bool_t applyElectronVeto = kTRUE, Bool_t invertElectronVeto = kFALSE);
    static Bool_t       PassConversionId(const Photon *p, const DecayParticle *c);
    static Bool_t       PassElectronVeto(const Photon *p, const ElectronCol *els);
    static Double_t     ElectronVetoCiC(const Photon *p, const ElectronCol *els);
    static Bool_t       PassElectronVetoConvRecovery(const Photon *p, const ElectronCol *els, const DecayParticleCol *conversions, const BaseVertex *v);
    static Bool_t       PassTriggerMatching(const Photon *p, const TriggerObjectCol *trigobjs);
    static const DecayParticle *MatchedConversion(const Photon *p, const DecayParticleCol *conversions, 
						  const BaseVertex *vtx, Int_t nWrongHitsMax=1, Double_t probMin=1e-6,
						  Double_t lxyMin = 2.0, Double_t dRMin = 0.1);
    static const DecayParticle *MatchedConversion(const SuperCluster *sc, const DecayParticleCol *conversions, 
                                                  const BaseVertex *vtx, Int_t nWrongHitsMax=1, Double_t probMin=1e-6,
                                                  Double_t lxyMin = 2.0, Double_t dRMin = 0.1);                                                   
    static const DecayParticle *MatchedConversion(const Track *t, const DecayParticleCol *conversions, 
						  const BaseVertex *vtx, Int_t nWrongHitsMax=1, Double_t probMin=1e-6,
						  Double_t lxyMin = 2.0);                                               
    static DiphotonR9EtaCats DiphotonR9EtaCat(const Photon *p1, const Photon *p2);
    static DiphotonR9EtaPtCats DiphotonR9EtaPtCat(const Photon *p1, const Photon *p2);
    static DiphotonR9EtaConversionCats DiphotonR9EtaConversionCat(const Photon *p1, const Photon *p2, const DecayParticleCol *conversions, const BaseVertex *v);
    static CiCBaseLineCats CiCBaseLineCat(const Photon *p);
    
    static const DecayParticle *MatchedCiCConversion(const Photon *p, const DecayParticleCol *conversions, 
						     Double_t dPhiMin=0.1, Double_t dEtaMin=0.1,Double_t dRMin=0.1, 
						     bool print   = false,
						     int* numLegs = NULL, int* convIdx = NULL);  // for debugging

                                                     
    static const Electron *MatchedElectron(const Photon *p, const ElectronCol *els);
    static const Photon *MatchedPhoton(const Electron *e, const PhotonCol *phs);
    static const SuperCluster *MatchedSC(const SuperCluster *psc, const SuperClusterCol *scs, Double_t drMin=0.3);

    static const SuperCluster *MatchedPFSC(const SuperCluster *psc, const PhotonCol *pfphos, const ElectronCol *eles, Double_t drMin=0.1);
    
    static bool PassCiCSelection(const Photon* ph, 
				 const Vertex* vtx, 
				 const TrackCol*    trackCol,
				 const ElectronCol* eleCol,
				 const VertexCol*   vtxCol,
				 double rho, double ptmin,
				 bool applyEleVeto = true,
				 bool print = false, float* kin=NULL);

    static bool PassCiCPFIsoSelection(const Photon* ph, 
				      const Vertex* vtx, 
				      const PFCandidateCol*    pfCol,
				      const VertexCol*   vtxCol,
				      double rho, double ptmin,bool dor9rescale, double p0b, double p1b,double p0e, double p1e, 
				      std::vector<double>* kin = NULL); 

    static bool PassCiCPFIsoSelectionWithEleVeto(const Photon* ph, 
						 const ElectronCol *els,
						 const DecayParticleCol *conversions, const BaseVertex *bs,
						 const Vertex* vtx, 
						 const PFCandidateCol*    pfCol,
						 const VertexCol*   vtxCol,
						 double rho, double ptmin,
						 Bool_t applyElectronVeto, Bool_t invertElectronVeto,
						 std::vector<double>* kin= NULL  // store variables for debugging...
						 );// add for mono photon                                
    
    static bool PassEgammaMediumSelectionWithEleVeto(const Photon* ph,
						     const ElectronCol *els,
						     const DecayParticleCol *conversions, const BaseVertex *bs,
						     const Vertex* vtx,
						     const PFCandidateCol* pfCol,
						     const VertexCol* vtxCol,
						     double rho, double ptmin,
						     Bool_t applyElectronVeto, Bool_t invertElectronVeto,
						     std::vector<double>* kin= NULL // store variables for debugging...
						     );// add for mono photon 

    static bool PassVgamma2011Selection(const Photon* ph, double rho);

    static const MCParticle *MatchMC(const Particle *ph, const MCParticleCol *c, Bool_t matchElectrons = kFALSE);
    ClassDef(PhotonTools, 0) // Muon tools
      };
}

#endif
