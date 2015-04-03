// $Id: ParticleColLinkDef.h,v 1.1 2009/06/15 15:00:15 loizides Exp $

#ifndef MITANA_DATATREE_FAKEEVENTHEADERCOLLINKDEF_H
#define MITANA_DATATREE_FAKEEVENTHEADERCOLLINKDEF_H
#include "MitAna/DataCont/interface/Collection.h"
#include "MitAna/DataCont/interface/Array.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitPhysics/FakeMods/interface/FakeEventHeaderFwd.h"
#include "MitPhysics/FakeMods/interface/FakeEventHeader.h"
#include "MitPhysics/FakeMods/interface/FakeObjectFwd.h"
#include "MitPhysics/FakeMods/interface/FakeObject.h"
#endif

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::FakeEventHeader+;
#pragma link C++ class mithep::Collection<mithep::FakeEventHeader>+;
#pragma link C++ class mithep::Array<mithep::FakeEventHeader>+;
#pragma link C++ class mithep::ObjArray<mithep::FakeEventHeader>+;
#pragma link C++ typedef mithep::FakeEventHeaderCol;
#pragma link C++ typedef mithep::FakeEventHeaderOArr;

#pragma link C++ class mithep::FakeObject+;
#pragma link C++ class mithep::Collection<mithep::FakeObject>+;
#pragma link C++ class mithep::Array<mithep::FakeObject>+;
#pragma link C++ class mithep::ObjArray<mithep::FakeObject>+;
#pragma link C++ typedef mithep::FakeObjectCol;
#pragma link C++ typedef mithep::FakeObjectArr;
#pragma link C++ typedef mithep::FakeObjectOArr;
#endif
