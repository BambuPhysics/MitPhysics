#include <iostream>

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
  fBeta(1.0),
  fD0Cut(0.03),
  fDZCut(0.1)
{
  // Constructor.
  std::cout << "Hello, this is dog" << std::endl;
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

  // This mapper will return particles that are only close in eta, phi space
  ParticleMapper *Mapper = new ParticleMapper();
  Mapper->Initialize(*fPFCandidates,fR0,fR0);

  const Int_t numCandidates = fPFCandidates->GetEntries();
  Double_t alphaF[numCandidates];                             // This is alpha(F) of all particles for particle i
  Double_t alphaC[numCandidates];                             // This is alpha(C) of charged PV for particle i
  Int_t numCHPU = 0;                                          // Track the number of particles that are charged PU
  Double_t alphaFCHPU[numCandidates];                         // This is alphaF for charged PU particles
  Double_t alphaCCHPU[numCandidates];                         // This is alphaC for charged PU particles
  Int_t IndicesF[numCandidates];                              // A bunch of indices used for sorting
  Int_t IndicesC[numCandidates];
  Int_t numFCHPUis0 = 0;
  Int_t numCCHPUis0 = 0;
  Int_t isChargedPU = 0;                                      // Bool, but I made it an int in case we want more options

  for(Int_t i0 = 0; i0 < numCandidates; i0++){
    const PFCandidate *iCandidate = fPFCandidates->At(i0);
    // Determine if the PFCandidate is charged PU. This will help characterize neutral PU.
    if(iCandidate->PFType() == 1 && 
       ((fabs(iCandidate->SourceVertex().Z() - fVertexes->At(0)->Position().Z()) > fDZCut) ||
        (MathUtils::AddInQuadrature(iCandidate->SourceVertex().X() - fVertexes->At(0)->Position().X(),
                                    iCandidate->SourceVertex().Y() - fVertexes->At(0)->Position().Y()) > fD0Cut))){
      isChargedPU = 1;
    }
    else isChargedPU = 0;
      
    // Initialize alphas and indices for sorting
    alphaF[i0] = 0;
    alphaC[i0] = 0;
    IndicesF[i0] = i0;
    IndicesC[i0] = i0;
    // Mapper is getting nearby PFCandidates
    std::vector<Int_t> nearList = Mapper->GetSurrounding(i0);
    for(UInt_t i1 = 0; i1 < nearList.size(); i1++){
      Int_t i2 = nearList[i1];                                      // Getting index from mapper results
      if(i0 == i2) continue;                                        // Don't bother comparing a particle to itself
      const PFCandidate *jCandidate = fPFCandidates->At(i2);
      Double_t dRTemp = MathUtils::DeltaR(iCandidate,jCandidate);
      if(dRTemp > fRMin && dRTemp < fR0){                                          // Only look at other particles in this range
        alphaF[i0] = alphaF[i0] + jCandidate->Pt()/(pow(jCandidate->Pt(),fBeta));  // First do the sum inside the log
        if((jCandidate->PFType() == 1) && 
           (fabs(jCandidate->SourceVertex().Z() - fVertexes->At(0)->Position().Z()) < fDZCut) && 
           (MathUtils::AddInQuadrature(jCandidate->SourceVertex().X() - fVertexes->At(0)->Position().X(),
                                       jCandidate->SourceVertex().Y() - fVertexes->At(0)->Position().Y()) < fD0Cut)){
          alphaC[i0] = alphaC[i0] + jCandidate->Pt()/(pow(jCandidate->Pt(),fBeta));
        }
      }
    }
    if(alphaF[i0] == 0) alphaF[i0] = -100;                          // Take the logs and ignore sum == 0 particles
    else alphaF[i0] = TMath::Log(alphaF[i0]);
    if(alphaC[i0] == 0) alphaC[i0] = -100;
    else alphaC[i0] = TMath::Log(alphaC[i0]);
    if(isChargedPU == 1){
      if(alphaF[i0] == -100) numFCHPUis0++;                         // Count particles to ignore when taking the median
      if(alphaC[i0] == -100) numCCHPUis0++;
      alphaFCHPU[numCHPU] = alphaF[i0];
      alphaCCHPU[numCHPU] = alphaC[i0];
      numCHPU++;                                                    // Count the total number of charged PU particles
    }
  }

  // Set the indices for sorting for the charged PU alphas
  Int_t IndicesFCHPU[numCHPU];
  Int_t IndicesCCHPU[numCHPU];

  for(Int_t i0 = 0; i0 < numCHPU; i0++){
    IndicesFCHPU[i0] = i0;
    IndicesCCHPU[i0] = i0;
  }

  // Sort everything to find the median and drop particles faster
  // This gives back the Indices arrays in ordered form
  TMath::Sort(numCHPU,alphaFCHPU,IndicesFCHPU,0);
  TMath::Sort(numCHPU,alphaCCHPU,IndicesCCHPU,0);
  TMath::Sort(numCandidates,alphaF,IndicesF,0);
  TMath::Sort(numCandidates,alphaC,IndicesC,0);

  // Now we'll find the median and sigma (left-handed RMS) squared for each event
  Double_t alphaFMed = 0;
  Double_t alphaCMed = 0;
  Double_t sigma2F = 0;
  Double_t sigma2C = 0;

  Int_t medIndexF = (numCHPU + numFCHPUis0)/2;                       // These are sort of meta, but they are the index of the indices array that refers to the median
  Int_t medIndexC = (numCHPU + numCCHPUis0)/2;                       // Just watch how they are used...
  if(numCHPU == numFCHPUis0) alphaFMed = 0;
  else if((numCHPU - numFCHPUis0) % 2 == 0) alphaFMed = (alphaFCHPU[IndicesFCHPU[medIndexF - 1]] + alphaFCHPU[IndicesFCHPU[medIndexF]])/2;
  else alphaFMed = alphaFCHPU[IndicesFCHPU[medIndexF]];
  if(numCHPU == numCCHPUis0) alphaCMed = 0;
  else if((numCHPU - numCCHPUis0) % 2 == 0) alphaCMed = (alphaCCHPU[IndicesCCHPU[medIndexC - 1]] + alphaCCHPU[IndicesCCHPU[medIndexC]])/2;
  else alphaCMed = alphaCCHPU[IndicesCCHPU[medIndexC]];

  // Now compute the sigma2s
  for(Int_t i0 = numFCHPUis0; i0 < medIndexF; i0++) sigma2F = sigma2F + pow((alphaFMed-alphaFCHPU[IndicesFCHPU[i0]]),2);
  sigma2F = sigma2F/(medIndexF - numFCHPUis0);
  for(Int_t i0 = numCCHPUis0; i0 < medIndexC; i0++) sigma2C = sigma2C + pow((alphaCMed-alphaCCHPU[IndicesCCHPU[i0]]),2);
  sigma2C = sigma2C/(medIndexC - numCCHPUis0);

  // Time to compute the chi2 values for each particle and assign the weights
  Double_t chi2F[numCandidates];                                     // This is alpha(F) of all particles for particle i
  Double_t chi2C[numCandidates];                                     // This is alpha(C) of charged PV for particle i
  Double_t weightsF[numCandidates];
  Double_t weightsC[numCandidates];
  for(Int_t i0 = 0; i0 < numCandidates; i0++){
    if(alphaF[i0] < alphaFMed) chi2F[i0] = 0;
    else chi2F[i0] = pow((alphaF[i0] - alphaFMed),2)/sigma2F;
    if(alphaC[i0] < alphaCMed) chi2C[i0] = 0;
    else chi2C[i0] = pow((alphaC[i0] - alphaCMed),2)/sigma2C;

    weightsF[i0] = TMath::Prob(chi2F[i0],1);
    weightsC[i0] = TMath::Prob(chi2C[i0],1);
  }  

  // prepare the storage array for the PuppiParticles
  PFCandidateOArr *PuppiParticleCol = new PFCandidateOArr;
  PuppiParticleCol->SetOwner(kTRUE);
  PuppiParticleCol->SetName(fPuppiParticlesName);

  for(Int_t i0 = 0; i0 < numCandidates; i0++){
    PFCandidate *PuppiParticle = fPFCandidates->At(i0)->MakeCopy();
    if(PuppiParticle->PFType() == 6 || PuppiParticle->PFType() == 7)
      PuppiParticle->SetPtEtaPhiM(PuppiParticle->Pt()*(weightsF[i0]),PuppiParticle->Eta(),PuppiParticle->Phi(),PuppiParticle->Mass()*(weightsF[i0]));
    else PuppiParticle->SetPtEtaPhiM(PuppiParticle->Pt()*(weightsC[i0]),PuppiParticle->Eta(),PuppiParticle->Phi(),PuppiParticle->Mass()*(weightsC[i0]));

    // add PuppiParticle to the collection
    PuppiParticleCol->AddOwned(PuppiParticle);
  }
  
  // sort according to ptrootcint forward declaration data members ??
  PuppiParticleCol->Sort();

  // add to event for other modules to use
  AddObjThisEvt(PuppiParticleCol);
}
