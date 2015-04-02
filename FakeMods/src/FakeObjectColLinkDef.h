// $Id: ParticleColLinkDef.h,v 1.1 2009/06/15 15:00:15 loizides Exp $

#ifndef MITANA_DATATREE_FAKEOBJECTCOLLINKDEF_H
#define MITANA_DATATREE_FAKEOBJECTCOLLINKDEF_H
#include "MitPhysics/FakeMods/interface/FakeObjectCol.h"
#endif

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::FakeObject+;
#pragma link C++ class mithep::Collection<mithep::FakeObject>+;
#pragma link C++ class mithep::Array<mithep::FakeObject>+;
#pragma link C++ class mithep::ObjArray<mithep::FakeObject>+;
#pragma link C++ typedef mithep::FakeObjectCol;
#pragma link C++ typedef mithep::FakeObjectArr;
#pragma link C++ typedef mithep::FakeObjectOArr;
#endif
