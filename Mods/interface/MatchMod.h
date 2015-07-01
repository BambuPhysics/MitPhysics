//--------------------------------------------------------------------------------------------------
// MatchMod
// 
// This module compares two input collections geomatrically, and creates a new 
// collection of Pair objects that match  within a given radius in eta-phi.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

// Shouldn't this be in MitAna/PhysicsMods ?? (Y.Iiyama)

#ifndef MITPHYSICS_MODS_MATCHMOD_H
#define MITPHYSICS_MODS_MATCHMOD_H

#include "MitAna/DataCont/interface/Collection.h" 
#include "MitAna/DataCont/interface/ObjArray.h" 
#include "MitAna/DataTree/interface/Pair.h" 
#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitPhysics/Utils/interface/MatchingTools.h" 
#include <TH1D.h>

namespace mithep 
{
  template<class ClA, class ClB=ClA>
  class MatchMod : public BaseMod
  {
    public:
      MatchMod(const char *name="MatchMod", 
               const char *title="Generic matching module");

      Bool_t                   GetAddPairs(Bool_t b)     const { return fAddPairs;     }
      const char              *GetColNameA()             const { return fColNameA;     }
      const char              *GetColNameB()             const { return fColNameB;     }
      const char              *GetColNameC()             const { return fColNameC;     }
      const char              *GetOutputName()           const { return GetColNameC(); }
      const char              *GetInputNameA()           const { return GetColNameA(); }
      const char              *GetInputNameB()           const { return GetColNameB(); }
      Double_t                 GetMatchRadius(Double_t m)      { return fMaxR;         }
      void                     SetColNameA(const char *n)      { fColNameA = n;        }
      void                     SetColNameB(const char *n)      { fColNameB = n;        }
      void                     SetColNameC(const char *n)      { fColNameC = n;        }
      void                     SetInputNameA(const char *n)    { SetColNameA(n);       }
      void                     SetInputNameB(const char *n)    { SetColNameB(n);       }
      void                     SetOutputName(const char *n)    { SetColNameC(n);       }
      void                     SetMatchRadius(Double_t m)      { fMaxR     = m;        }
      void                     SetAddPairs(Bool_t b)           { fAddPairs = b;        }

    protected:
      void                     Process();
      void                     SlaveBegin();

      TString                  fColNameA;     //name of input  collection A
      TString                  fColNameB;     //name of input  collection B
      TString                  fColNameC;     //name of output collection C
      Double_t                 fMaxR;         //maximum radius for matching (def=0.1)
      Bool_t                   fAddPairs;     //=true then add Pair into output collection (def=1)

      ClassDef(MatchMod,1) // Generic matching module
  };
}

//--------------------------------------------------------------------------------------------------
template<class ClA, class ClB>
mithep::MatchMod<ClA,ClB>::MatchMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fMaxR(0.1),
  fAddPairs(kTRUE)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
template<class ClA, class ClB>
void mithep::MatchMod<ClA,ClB>::Process()
{
  // Process the two collection and find matching pairs.

  using namespace mithep;

  auto* colA = GetObject<Collection<ClA> >(GetColNameA());
  if (!colA) {
    this->SendError(kAbortModule, 
                    "Process", 
                    "Could not obtain collection with name %s!", GetColNameA());
    return;
  }
  
  auto* colB = GetObject<Collection<ClB> >(GetColNameB());
  if (!colB) {
    this->SendError(kAbortModule, 
                    "Process", 
                    "Could not obtain collection with name %s!", GetColNameB());
    return;
  }

  if (fAddPairs) {
    ObjArray<Pair<ClA,ClB> > *out = new ObjArray<Pair<ClA,ClB> >(0,GetOutputName());
    out->SetOwner(kTRUE);
    const UInt_t ents = colA->GetEntries();
    for(UInt_t i=0;i<ents;++i) {
      const ClA *pA = colA->At(i);
      const ClB *pB = MatchingTools::Closest(pA,*colB,fMaxR,kFALSE);
      if (!pB) 
        continue;
      out->AddOwned(new Pair<ClA,ClB>(pA,pB));
    }
    AddObjThisEvt(out);
  } else {
    ObjArray<ClA> *out = new ObjArray<ClA>(0,GetOutputName());
    out->SetOwner(kFALSE);
    const UInt_t ents = colA->GetEntries();
    for(UInt_t i=0;i<ents;++i) {
      const ClA *pA = colA->At(i);
      const ClB *pB = MatchingTools::Closest(pA,*colB,fMaxR,kFALSE);
      if (!pB) 
        continue;
      out->Add(pA);
    }
    AddObjThisEvt(out);
  }
}

//--------------------------------------------------------------------------------------------------
template<class ClA, class ClB>
void mithep::MatchMod<ClA,ClB>::SlaveBegin()
{
  // Run starting code on the client to do some consistency checks.


  if (fColNameA.IsNull()) {
    SendError(kAbortModule, "SlaveBegin", "No name given for input collection A.");
    return;
  }

  if (fColNameB.IsNull()) {
    SendError(kAbortModule, "SlaveBegin", "No name given for input collection B.");
    return;
  }

  if (fColNameC.IsNull()) {
    fColNameC="Matched";
    fColNameC+=fColNameA;
    fColNameC+=fColNameB;
    SendError(kWarning, "SlaveBegin", "No name given for output collection,"
              " will use \"%s\"", fColNameC.Data());
  }
}
#endif
