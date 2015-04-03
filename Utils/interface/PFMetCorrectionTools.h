#ifndef MITPHYSICS_UTILS_PFMetCorrectionTOOLS_H
#define MITPHYSICS_UTILS_PFMetCorrectionTOOLS_H

#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/Particle.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/CaloTowerCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include <TVector3.h>
#include <TLorentzVector.h>

namespace mithep {
  class PFMetCorrectionTools {
  public:
    PFMetCorrectionTools();
    virtual ~PFMetCorrectionTools() {}
    static Double_t ErrEt( Double_t Et, Double_t Eta);
    static void correctMet(Met *met, const Photon *phHard,const Photon *phSoft, Bool_t smearing, Bool_t scale, const PFJetCol *fPFJet, const GenJetCol *fGenJet, const JetCol *fcorrJet, UInt_t evt);
    static void shiftMet(Met *uncormet, Bool_t fIsData, Double_t spfMet);
    ClassDef(PFMetCorrectionTools, 0)
  };
}

#endif
