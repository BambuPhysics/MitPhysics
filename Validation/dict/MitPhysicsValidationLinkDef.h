// $Id: MitPhysicsValidationLinkDef.h,v 1.1 2008/10/15 06:05:02 loizides Exp $

#ifndef MITPHYSICS_VALIDATION_LINKDEF_H
#define MITPHYSICS_VALIDATION_LINKDEF_H
#include "MitPhysics/Validation/interface/MCParticlesValMod.h"
#include "MitPhysics/Validation/interface/TracksValMod.h"
#include "MitPhysics/Validation/interface/JetValidationMod.h"
#endif

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::JetValidationMod+;
#pragma link C++ class mithep::MCParticlesValMod+;
#pragma link C++ class mithep::TracksValMod+;
#endif
