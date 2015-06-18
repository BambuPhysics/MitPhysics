#ifndef MITPHYSICS_UTILS_METLEPTONTOOLS_H
#define MITPHYSICS_UTILS_METLEPTONTOOLS_H

#include <TMatrixD.h>
#include "MitAna/DataTree/interface/PFJetFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/MuonFwd.h"
#include "MitAna/DataTree/interface/PFTauFwd.h"
#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/PhotonFwd.h"
#include "MitAna/DataTree/interface/PFMetFwd.h"
#include "MitAna/DataTree/interface/PFJetFwd.h"
#include "MitAna/DataTree/interface/PFCandidateFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityFwd.h"
#include "MitPhysics/Utils/interface/TauIsoMVA.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include <TVector3.h>
#include <TLorentzVector.h>

namespace mithep {
  class MetLeptonTools {
  public:
    MetLeptonTools();
    virtual ~MetLeptonTools() {}
    TauIsoMVA     *fTauIsoMVA; 
    bool           looseTauId(const PFTau *iTau,const PileupEnergyDensityCol* iPUEnergyDensity);
    static bool    looseEleId(const Electron *iElectron,const PileupEnergyDensityCol* iPUEnergyDensity,
			      const PFCandidateCol *iCands,const Vertex *iPV,const VertexCol *iVertices);
    static bool    looseMuId(const Muon *iMu,const PFCandidateCol *iCands,const Vertex *iPV,const VertexCol *iVertices);
    static bool    loosePhotonId(const Photon *iPhoton);
    static double  vis(const PFTau *iTau);
    static Float_t PFIsolation(const ChargedParticle *iLep,const PFCandidateCol *iCands);  
    static Float_t PFIsolationNoGamma(const ChargedParticle *iLep,const PFCandidateCol *iCands);  
    static Float_t isoPV(const ChargedParticle *iLep,const PFCandidateCol *iCands,
                         const Vertex *iPV,const VertexCol *iVertices,bool iEle=false);
    static Float_t chargedFracInCone(const Photon *iPhoton,const PFCandidateCol *iCands,
                                     const Vertex *iPV);
  };
}
#endif
