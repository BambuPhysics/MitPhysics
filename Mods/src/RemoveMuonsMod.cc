#include <iostream>

#include "MitPhysics/Mods/interface/RemoveMuonsMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
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
  fPFCandidateName(Names::gkPFCandidatesBrn),
  fOutputName("PFCandidateMinusMuons"),
  fPFCandidateFromBranch(true),
  fDeltaR(0.1),
  fGettingTriggerMatch(false)
{ }

//--------------------------------------------------------------------------------------------------
RemoveMuonsMod::~RemoveMuonsMod() {}

//--------------------------------------------------------------------------------------------------
void RemoveMuonsMod::SlaveBegin()
{
  // ===== load branches ====
  if (fPFCandidateFromBranch)
    ReqBranch(fPFCandidateName, fPFCandidates);
  if (fMuonFromBranch)
    ReqBranch(fMuonName, fMuons);
}

//--------------------------------------------------------------------------------------------------
void RemoveMuonsMod::SlaveTerminate() {}

//--------------------------------------------------------------------------------------------------
void RemoveMuonsMod::Process()
{

  // Process entries of the tree. 
  if (fPFCandidateFromBranch) 
    LoadBranch(fPFCandidateName);
  else
    LoadEventObject(fPFCandidateName, fPFCandidates);

  if (fMuonFromBranch) 
    LoadBranch(fMuonName);
  else
    LoadEventObject(fMuonName, fMuons);

  ParticleMapper *mapper = new ParticleMapper();
  mapper->Initialize(*fPFCandidates,fDeltaR,fDeltaR);

  PFCandidateArr *PFParticlesMinusMuons = new PFCandidateArr(16);
  PFParticlesMinusMuons->SetName(fOutputName);

  LoadEventObject(fTriggerObjectsName, fTriggerObjects);

  std::vector<Int_t> removeThese;

  for (UInt_t i0 = 0; i0 < fTriggerObjects->GetEntries(); i0++) {
    if (fGettingTriggerMatch) {
      std::cout << fTriggerObjects->At(i0)->ModuleName() << std::endl;
      fGettingTriggerMatch = false;
    }
    if (fTriggerMatchName == fTriggerObjects->At(i0)->ModuleName()) {
      std::vector<Int_t> compareThese = mapper->GetNearEtaPhi(fTriggerObjects->At(i0)->Eta(),
                                                              fTriggerObjects->At(i0)->Phi());
      for (UInt_t i1 = 0; i1 < compareThese.size(); i1++){
        if (fPFCandidates->At(compareThese[i1])->PFType() == 3 &&
            MathUtils::DeltaR(fPFCandidates->At(compareThese[i1]),fTriggerObjects->At(i0)) < fDeltaR) {
          Bool_t found = false;
          for (UInt_t i2 = 0; i2 < removeThese.size(); i2++) {
            if (compareThese[i1] == removeThese[i2]) {
              found = true;
              break;
            }
          }
          if (!found)
            removeThese.push_back(compareThese[i1]);
        }
      }
    }
  }
  std::cout << removeThese.size() << std::endl;
}
