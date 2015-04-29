// $Id: MitPhysicsSkimLinkDef.h,v 1.2 2012/06/03 20:25:56 paus Exp $

#ifndef MITPHYSICS_SKIM_LINKDEF_H
#define MITPHYSICS_SKIM_LINKDEF_H
#include "MitPhysics/Skim/interface/BaseH4lSkim.h"
#include "MitPhysics/Skim/interface/H4lSkim.h"
#include "MitPhysics/Skim/interface/H4lEleTagProbeSkim.h"
#include "MitPhysics/Skim/interface/H4lMuTagProbeSkim.h"
#include "MitPhysics/Skim/interface/H4lZPlusFakeSkim.h"
#include "MitPhysics/Skim/interface/H4lLightFlavorSkim.h"
#endif

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;
#pragma link C++ class mithep::BaseH4lSkim+;
#pragma link C++ class mithep::H4lSkim+;
#pragma link C++ class mithep::H4lEleTagProbeSkim+;
#pragma link C++ class mithep::H4lMuTagProbeSkim+;
#pragma link C++ class mithep::H4lZPlusFakeSkim+;
#pragma link C++ class mithep::H4lLightFlavorSkim+;
#endif
