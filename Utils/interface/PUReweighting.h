//--------------------------------------------------------------------------------------------------
// $Id: 
//
// PUReweighting
//
// Authors: M. Zanetti
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_PUREWEIGHTING_H
#define MITPHYSICS_UTILS_PUREWEIGHTING_H

#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include <TH1F.h>

namespace mithep {

  class PUReweighting {

  public:

    PUReweighting(const char* fileName, const char* histoName);

    virtual ~PUReweighting() {}

    double reweightOOT(int inTimePUMultiplicity, int outOfTimePUMultiplicity);
			     
  private:
    TH1F * referenceHisto;


    ClassDef(PUReweighting, 1) // PUReweighting
      };


}

#endif
