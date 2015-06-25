#include "MitAna/DataTree/interface/JetCol.h"
#include "MitPhysics/Mods/interface/PuppiMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/DataFormats/interface/Vect3.h"
#include "MitAna/DataTree/interface/VertexCol.h"

#include "TLorentzVector.h"

using namespace mithep;

ClassImp(mithep::PuppiMod)

//--------------------------------------------------------------------------------------------------
PuppiMod::PuppiMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fVertexesName(Names::gkPVBrn),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPuppiParticlesName("PuppiParticles"),
  fPFCandidates(0),
  fRMin(0.02),
  fR0(0.3),
  fBeta(1.0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
PuppiMod::~PuppiMod() {} 

//--------------------------------------------------------------------------------------------------
void PuppiMod::SlaveBegin()
{

  // ===== load branches ====
  ReqBranch(fVertexesName, fVertexes);
  ReqBranch(fPFCandidatesName, fPFCandidates);

}

//--------------------------------------------------------------------------------------------------
void PuppiMod::SlaveTerminate()
{
  // ===== deallocate memory ====
}

//--------------------------------------------------------------------------------------------------
void PuppiMod::Process()
{
  // Process entries of the tree. 

  LoadBranch(fVertexesName);
  LoadBranch(fPFCandidatesName);

  ParticleMapper *Mapper = new ParticleMapper();
  Mapper->Initialize(*fPFCandidates,fR0,fR0);

  const Int_t numCandidates = fPFCandidates->GetEntries();
  Double_t alphaF[numCandidates];
  Double_t alphaC[numCandidates];

  for(Int_t i0 = 0; i0 < numCandidates; i0++){
    const PFCandidate *iCandidate = fPFCandidates->At(i0);
    alphaF[i0] = 0;
    alphaC[i0] = 0;
    std::vector<Int_t> nearList = Mapper->GetSurrounding(i0);
    for(UInt_t i1 = 0; i1 < nearList.size(); i1++){
      Int_t i2 = nearList[i1];
      if(i0 == i2) continue;
      const PFCandidate *jCandidate = fPFCandidates->At(i2);
      Double_t dRTemp = MathUtils::DeltaR(iCandidate,jCandidate);
      if(dRTemp > fRMin && dRTemp < fR0){
        alphaF[i0] = alphaF[i0] + jCandidate->Pt()/(pow(jCandidate->Pt(),fBeta));
        if((abs(jCandidate->SourceVertex().Z() - fVertexes->At(0)->Position().Z()) < fDZCut) && 
           (MathUtils::AddInQuadrature(jCandidate->SourceVertex().X() - fVertexes->At(0)->Position().X(),
                                       jCandidate->SourceVertex().Y() - fVertexes->At(0)->Position().Y()) < fD0Cut)){
          alphaC[i0] = alphaC[i0] + jCandidate->Pt()/(pow(jCandidate->Pt(),fBeta));
        }
      }
    }
  }

  PFCandidate *PuppiParticle = fPFCandidates->At(0)->MakeCopy();

  // prepare the storage array for the corrected MET
  PFCandidateOArr *PuppiParticleCol = new PFCandidateOArr;
  PuppiParticleCol->SetOwner(kTRUE);
  PuppiParticleCol->SetName(fPuppiParticlesName);

  // add corrected met to collection
  PuppiParticleCol->AddOwned(PuppiParticle);
  
  // sort according to ptrootcint forward declaration data members
  PuppiParticleCol->Sort();
  
  // add to event for other modules to use
  AddObjThisEvt(PuppiParticleCol);
}
