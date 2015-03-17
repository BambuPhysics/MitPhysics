#ifndef MITPHYSICS_UTILS_MVAVBF_H
#define MITPHYSICS_UTILS_MVAVBF_H

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
namespace TMVA {//MVA
  class Reader;
}

namespace mithep {
  class MVAVBF {
    
  public:
    
    MVAVBF();
    virtual ~MVAVBF() {}
    void InitializeMVA();
    float GetMVAbdtValue(float jet1pt, float jet2pt, float deltajeteta, float dijetmass, float zeppenfeld, float dphidijetgg, float diphoptOverdiphomass, float pho1ptOverdiphomass, float pho2ptOverdiphomass);
    
  private:
    
    TMVA::Reader *fReader;

    float _jet1pt;
    float _jet2pt;
    float _deltajeteta;
    float _dijetmass;
    float _zeppenfeld;
    float _dphidijetgg;
    float _diphoptOverdiphomass;
    float _pho1ptOverdiphomass;
    float _pho2ptOverdiphomass;

    TMVA::Reader *reader;    

    ClassDef(MVAVBF, 0) 
      };
}

#endif
