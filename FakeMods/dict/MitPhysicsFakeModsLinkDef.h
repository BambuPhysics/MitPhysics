// $Id: MitPhysicsModsLinkDef.h,v 1.14 2009/04/30 08:09:32 loizides Exp $

#ifndef MITPHYSICS_FAKEMODS_LINKDEF_H
#define MITPHYSICS_FAKEMODS_LINKDEF_H
#include "MitPhysics/FakeMods/interface/FakeLeptonExampleAnaMod.h"
#include "MitPhysics/FakeMods/interface/FakeRate.h"
#include "MitPhysics/FakeMods/interface/GenFakesMod.h"
#include "MitPhysics/FakeMods/interface/GenFakeableObjsMod.h"
#endif

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::FakeLeptonExampleAnaMod+;
#pragma link C++ class mithep::GenFakesMod+;
#pragma link C++ class mithep::GenFakeableObjsMod+;
#pragma link C++ class mithep::FakeRate+;
#endif
