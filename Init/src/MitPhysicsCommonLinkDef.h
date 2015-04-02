// $Id: MitPhysicsModsLinkDef.h,v 1.2 2008/11/19 15:44:50 loizides Exp $

#ifndef MITPHYSICS_INIT_LINKDEF_H
#define MITPHYSICS_INIT_LINKDEF_H
#include "MitPhysics/Init/interface/ModNames.h"
#endif

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::ModNames+;
#endif
