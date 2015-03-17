//--------------------------------------------------------------------------------------------------
// $Id: MatchingTools.h,v 1.4 2009/06/28 08:03:09 loizides Exp $
//
// MatchingTools
//
// This class implements a couple of simple geometrical matching functions.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_MATCHINGTOOLS_H
#define MITPHYSICS_UTILS_MATCHINGTOOLS_H

#include <TObjArray.h>
#include <TMath.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataCont/interface/Collection.h"

namespace mithep
{
  class MatchingTools {
    public:
      template<class V1, class V2> 
      static const V2 *Closest(const V1 *v1, const Collection<V2> &col, Double_t maxR, 
                               Bool_t self=kFALSE);
      template<class V1, class V2> 
      static const V2 *Closest(const V1 *v1, const Collection<V2> &col, Double_t maxR,
                               Double_t maxPtDiff, Bool_t self=kFALSE);
      template<class V1, class V2> 
      static TObjArray *Closests(const V1 *v1, const Collection<V2> &col, Double_t maxR,
                                 Bool_t self=kFALSE);
  };

  //------------------------------------------------------------------------------------------------
  template<class V1, class V2>
  const V2 *MatchingTools::Closest(const V1 *v1, const Collection<V2> &col, Double_t maxR, 
                                   Bool_t self)
  {
    // Return closest geometrical neighbor in eta-phi within maximum DeltaR of maxR.
    // If self is kFALSE make sure that identical pointers are excluded.
    
    if (!v1)
      return 0;

    const V2 *res = 0;
    const UInt_t ents = col.GetEntries();
    if (ents>0) {
      Double_t dR = 1e30;
      for (UInt_t i = 0; i<ents; ++i) {
        const V2 *v2 = col.At(i);
        if (!v2) 
          continue;
        if (!self && (void*)v1==(void*)v2) 
          continue;
        Double_t diff = MathUtils::DeltaR(*v1,*v2);
        if ((diff<maxR) && (diff<dR)) {
          res = v2;
          dR  = diff;
        }
      }
    }
    return res;
  }

  //------------------------------------------------------------------------------------------------
  template<class V1, class V2>
  const V2 *MatchingTools::Closest(const V1 *v1, const Collection<V2> &col, Double_t maxR,
                                   Double_t maxPtDiff, Bool_t self)
  {
    // Return closest geometrical neighbor in eta-phi within maximum DeltaR of maxR.
    // Exclude particles where the relative pt differs by more than maxPtDiff.
    // If self is kFALSE make sure that identical pointers are excluded.
    
    if (!v1)
      return 0;

    const V2 *res = 0;
    const UInt_t ents = col.GetEntries();
    if (ents>0) {
      Double_t dR = 1e30;
      for (UInt_t i = 0; i<ents; ++i) {
        const V2 *v2 = col.At(i);
        if (!v2) 
          continue;
        if (!self && (void*)v1==(void*)v2) 
          continue;
        Double_t ptd = TMath::Abs(1-v2->Pt()/v1->Pt());
        if (ptd>maxPtDiff)
          continue;
        Double_t diff = MathUtils::DeltaR(*v1,*v2);
        if ((diff<maxR) && (diff<dR)) {
          res = v2;
          dR  = diff;
        }
      }
    }
    return res;
  }

  //------------------------------------------------------------------------------------------------
  template<class V1, class V2>
  TObjArray *MatchingTools::Closests(const V1 *v1, const Collection<V2> &col, Double_t maxR,
                                     Bool_t self)
  {
    // Return closest geometrical neighbors in eta-phi within maximum DeltaR of maxR.
    // If self is kFALSE make sure that identical pointers are excluded.

    if (!v1)
      return 0;

    TObjArray *res = new TObjArray;
    res->SetOwner(kFALSE);
    const UInt_t ents = col.GetEntries();
    for (UInt_t i = 0; i<ents; ++i) {
      const V2 *v2 = col.At(i);
      if (!v2)
        continue;
      if (!self && (void*)v1==(void*)v2) 
        continue;
      Double_t diff = MathUtils::DeltaR(*v1,*v2);
      if (diff<maxR)
        res->Add(const_cast<V2*>(v2));
    }
    return res;
  }
}
#endif
