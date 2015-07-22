#ifndef MITPHYSICS_SELMODS_LINKDEF_H
#define MITPHYSICS_SELMODS_LINKDEF_H

#include "MitAna/DataCont/interface/Collection.h"
#include "MitAna/DataTree/interface/CompositeParticle.h"
#include "MitAna/DataTree/interface/Conversion.h"
#include "MitAna/DataTree/interface/Electron.h"
#include "MitAna/DataTree/interface/GenJet.h"
#include "MitAna/DataTree/interface/Jet.h"
#include "MitAna/DataTree/interface/Met.h"
#include "MitAna/DataTree/interface/Muon.h"
#include "MitAna/DataTree/interface/Particle.h"
#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/Track.h"
#include "MitAna/DataTree/interface/TriggerObject.h"
#include "MitPhysics/SelMods/interface/DilepSelMod.h"
#include "MitPhysics/SelMods/interface/GenericSelMod.h"
#include "MitPhysics/SelMods/interface/JetPlusIsoTrackSelMod.h"
#include "MitPhysics/SelMods/interface/LeptonPlusIsoTrackSelMod.h"
#include "MitPhysics/SelMods/interface/PhotonPlusIsoTrackSelMod.h"
#include "MitPhysics/SelMods/interface/HwwExampleAnalysisMod.h"
#include "MitPhysics/SelMods/interface/WBFExampleAnalysisMod.h"
#include "MitPhysics/SelMods/interface/EventCategoryMod.h"
#endif
 
#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::DilepSelMod+;
#pragma link C++ class mithep::GenericSelMod<mithep::Particle>+;
#pragma link C++ class mithep::GenericSelMod<mithep::CompositeParticle>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Conversion>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Electron>+;
#pragma link C++ class mithep::GenericSelMod<mithep::GenJet>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Jet>+;
#pragma link C++ class mithep::GenericSelMod<mithep::MCParticle>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Met>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Muon>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Photon>+;
#pragma link C++ class mithep::GenericSelMod<mithep::Track>+;
#pragma link C++ class mithep::GenericSelMod<mithep::TriggerObject>+;
#pragma link C++ class mithep::JetPlusIsoTrackSelMod+;
#pragma link C++ class mithep::LeptonPlusIsoTrackSelMod+;
#pragma link C++ class mithep::PhotonPlusIsoTrackSelMod+;
#pragma link C++ class mithep::HwwExampleAnalysisMod+;
#pragma link C++ class mithep::WBFExampleAnalysisMod+;
#pragma link C++ class mithep::EventCategoryMod+;
#endif
