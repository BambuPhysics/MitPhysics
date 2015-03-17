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
    Double_t *fWeights;
    UShort_t *fIndices;
    Double_t *fCache;
    Int_t    *fOffsets;
    Int_t    *fIndOffsets;
    Double_t *fFactors;

    ClassDef(PUReweightingMulti, 0) // PUReweighting
      };


}

#endif
