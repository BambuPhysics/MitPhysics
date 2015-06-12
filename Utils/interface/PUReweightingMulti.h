//--------------------------------------------------------------------------------------------------
// $Id: 
//
// PUReweightingMulti
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_PUREWEIGHTINGMULTI_H
#define MITPHYSICS_UTILS_PUREWEIGHTINGMULTI_H

#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include <THnSparse.h>
#include <TH3D.h>

namespace mithep {

  class PUReweightingMulti {

  public:

    PUReweightingMulti(const THnSparse *sparse, int maxn = 100, double minprob = 1e-12);

    virtual ~PUReweightingMulti();

    double TargetPdf(int *npus, int ndim) const;
    TH3D *Target3D();
    TH3D *Weights3D(const TH3D *source3d);
    
  protected:
    void Setup();
    
    
    const THnSparse *fSparse;
    int fNsBins;
    int fNDim;
    int fMaxN;
    double fMinProb;
    int fNCacheBins;
    Double_t *fWeights;    //[fNsBins]
    UShort_t *fIndices;    //[fNDim*fNsBins]
    Double_t *fCache;      //[fMaxN*fNCacheBins]
    Int_t    *fOffsets;    //[fNDim]
    Int_t    *fIndOffsets; //[fNDim]
    Double_t *fFactors;    //[fNDim]

    ClassDef(PUReweightingMulti, 1) // PUReweighting
      };


}

#endif
