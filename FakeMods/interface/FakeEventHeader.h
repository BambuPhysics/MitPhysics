//--------------------------------------------------------------------------------------------------
// $Id: FakeEventHeader.h,v 1.3 2009/07/20 19:05:04 loizides Exp $
//
// FakeEventHeader
//
// Class to hold information specific for events where some fakable object has been promoted 
// to a (faked) lepton. 
//
// Authors: S. Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_FAKEMODS_FAKEEVENTHEADER_H
#define MITPHYSICS_FAKEMODS_FAKEEVENTHEADER_H
 
#include "MitAna/DataTree/interface/DataObject.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitPhysics/FakeMods/interface/FakeObjectFwd.h"

namespace mithep 
{
  class FakeObject;

  class FakeEventHeader : public DataObject
  {
    public:
      FakeEventHeader() : fWeight(1.0), fWeightLowError(0.0), fWeightHighError(0.0) {} 
      ~FakeEventHeader() {}

      Double_t                  Weight()              const { return  fWeight;                   }
      Double_t                  WeightLowError()      const { return  fWeightLowError;           }
      Double_t                  WeightHighError()     const { return  fWeightHighError;          }
      const FakeObject         *FakeObj(UInt_t i)     const { return  fFakeObjects.At(i);        }
      UInt_t                    FakeObjsSize()        const { return  fFakeObjects.GetEntries(); }
      const Jet                *UnfakedJet(UInt_t i)  const { return  fJets.At(i);               }
      UInt_t                    NJets()               const { return  fJets.GetEntries();        }
      void                      SetWeight(Double_t w)                { fWeight = w;              }
      void                      SetWeightLowError(Double_t error)    { fWeightLowError = error;  }
      void                      SetWeightHighError(Double_t error)   { fWeightHighError = error; }
      void                      AddJets(const Jet *jet)              { fJets.Add(jet);           }
      void                      AddFakeObject(const Particle *p, EObjType faketype, 
                                              Bool_t fakeTag, Bool_t mcTag);
      void                      AddFakeObject(const mithep::FakeObject *fo);

    protected:

      Double_t                  fWeight;          //fake event weight
      Double_t                  fWeightLowError;  //fake event weight low error
      Double_t                  fWeightHighError; //fake event weight high error
      FakeObjectArr             fFakeObjects;     //fake objects
      JetOArr                   fJets;            //collection of jets after some have been 

    ClassDef(FakeEventHeader, 1) // Fake event header class
  };
}

#endif
