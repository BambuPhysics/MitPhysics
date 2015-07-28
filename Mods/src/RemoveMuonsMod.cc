#include <iostream>

#include "MitPhysics/Mods/interface/RemoveMuonsMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/ParticleMapper.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"

using namespace mithep;

ClassImp(mithep::RemoveMuonsMod)

//--------------------------------------------------------------------------------------------------
RemoveMuonsMod::RemoveMuonsMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fTriggerObjects(0),
  fPFCandidates(0),
  fTriggerMatchName(""),
  fTriggerObjectsName(Names::gkHltObjBrn),
  fInputName(Names::gkPFCandidatesBrn),
  fOutputName("PFCandidateMinusMuons"),
  fInputFromBranch(true),
  fDeltaR(0.1)
{ }

//--------------------------------------------------------------------------------------------------
RemoveMuonsMod::~RemoveMuonsMod() {}

//--------------------------------------------------------------------------------------------------
void RemoveMuonsMod::SlaveBegin()
{
  // ===== load branches ====
  if(fInputFromBranch)
    ReqBranch(fInputName, fPFCandidates);
}

//--------------------------------------------------------------------------------------------------
void RemoveMuonsMod::SlaveTerminate() {}

//--------------------------------------------------------------------------------------------------
void RemoveMuonsMod::Process()
{
  // Process entries of the tree. 
  if (fInputFromBranch) 
    LoadBranch(fInputName);
  else
    LoadEventObject(fInputName, fPFCandidates);

  PFCandidateArr *PFParticlesMinusMuons = new PFCandidateArr(16);
  PFParticlesMinusMuons->SetName(fOutputName);

  LoadEventObject(fTriggerObjectsName, fTriggerObjects);


}
