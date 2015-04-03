// $Id: $

#include "MitPhysics/FakeMods/interface/FakeEventHeader.h"
#include "MitPhysics/FakeMods/interface/FakeObject.h"

ClassImp(mithep::FakeEventHeader)

void
mithep::FakeEventHeader::AddFakeObject(const Particle *p, EObjType faketype, 
                                       Bool_t fakeTag, Bool_t mcTag)  
{
  // Add new fake object
  mithep::FakeObject *newFake = fFakeObjects.AddNew();
  newFake->SetParticle(p);
  newFake->SetFakeType(faketype);
  newFake->SetFakeTag(fakeTag);
  newFake->SetMCTag(mcTag);  
}

void
mithep::FakeEventHeader::AddFakeObject(const mithep::FakeObject *fo)  
{
  // Add new fake object with default parameters taken from the original object.
  mithep::FakeObject *newFake = fFakeObjects.AddNew();
  newFake->SetParticle(fo->FakeParticle());
  newFake->SetFakeType(fo->ObjType());
  newFake->SetFakeTag(fo->FakeTag());
  newFake->SetMCTag(fo->MCTag());  
}
