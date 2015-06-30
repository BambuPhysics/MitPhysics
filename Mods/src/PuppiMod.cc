#include <iostream>

#include "MitAna/DataTree/interface/JetCol.h"
#include "MitPhysics/Mods/interface/PuppiMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/DataFormats/interface/Vect3.h"
#include "MitAna/DataTree/interface/VertexCol.h"

#include "TLorentzVector.h"
#include "Math/ProbFunc.h"

using namespace mithep;

ClassImp(mithep::PuppiMod)

//--------------------------------------------------------------------------------------------------
PuppiMod::PuppiMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fVertexesName(Names::gkPVBrn),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPuppiParticlesName("PuppiParticles"),
  fPFCandidates(0),
  fPuppiParticles(0),
  fRMin(0.02),
  fR0(0.3),
  fAlpha(1.0),
  fBeta(1.0),
  fD0Cut(0.03),
  fDZCut(0.1),
  fMinWeightCut(0.01),
  fMinNeutralPt(1.0),
  fMinNeutralPtSlope(0.005),
  fKeepPileup(kFALSE),
  fInvert(kFALSE),
  fApplyCHS(kTRUE)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
PuppiMod::~PuppiMod() {} 

//--------------------------------------------------------------------------------------------------
Int_t PuppiMod::GetParticleType(const PFCandidate *cand){
  if(cand->PFType() == 1){
    if((fabs(cand->SourceVertex().Z() - fVertexes->At(0)->Position().Z()) < fDZCut) &&
       (MathUtils::AddInQuadrature(cand->SourceVertex().X() - fVertexes->At(0)->Position().X(),
                                   cand->SourceVertex().Y() - fVertexes->At(0)->Position().Y()) < fD0Cut))
      return 1;                               // This is charged PV particle
    else return 2;                            // This is charged PU particle
  }
  else if(cand->PFType() != 6 && cand->PFType() != 7)
    return 3;                                 // This is neutral particle in the center
  else return 4;                              // This is neutral particle in the forward region
  return 0;                                   // This shouldn't really happen
}

//--------------------------------------------------------------------------------------------------
Double_t PuppiMod::Chi2fromDZ(Double_t dz){
  Double_t expSig = 1.0;                      // This is uncertainty of track, consider making it a parameter
  Double_t probPV = ROOT::Math::normal_cdf_c(fabs(dz),expSig) * 2.0;
  Double_t probPU = 1 - probPV;
  if(probPU == 1) probPU = 1 - 1e-16;         // This is a weird thing due to the next function I call
  Double_t chi = TMath::ChisquareQuantile(probPU,1);
  return pow(chi,2);
}

//--------------------------------------------------------------------------------------------------
void PuppiMod::SlaveBegin()
{
  // ===== load branches ====
  ReqBranch(fVertexesName, fVertexes);
  ReqBranch(fPFCandidatesName, fPFCandidates);

  // prepare the storage array for the PuppiParticles
  fPuppiParticles = new PFCandidateArr(16);
  fPuppiParticles->SetName(fPuppiParticlesName);
  PublishObj(fPuppiParticles);
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

  for(Int_t i0 = 0; i0 < numCandidates; i0++){
    const PFCandidate *iCandidate = fPFCandidates->At(i0);
    // Determine if the PFCandidate is charged PU. This will help characterize neutral PU.
    Int_t iParticleType = GetParticleType(iCandidate);
      
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
      if(dRTemp > fRMin && dRTemp < fR0){                                                      // Only look at other particles in this range
        Int_t jParticleType = GetParticleType(jCandidate);
        Double_t theAddition = (pow(jCandidate->Pt(),fAlpha))/(pow(dRTemp,fBeta));             // This is the thing we have to add inside the log
        alphaF[i0] = alphaF[i0] + theAddition;                                                 // First do the sum inside the log (alphaF)
        if(jParticleType == 1)                                                                 // If the particle is charged PV
          alphaC[i0] = alphaC[i0] + theAddition;                                               //   add to alphaC
      }
    }
    if(alphaF[i0] == 0) alphaF[i0] = -100;                          // Take the logs and ignore sum == 0 particles
    else alphaF[i0] = TMath::Log(alphaF[i0]);
    if(alphaC[i0] == 0) alphaC[i0] = -100;
    else alphaC[i0] = TMath::Log(alphaC[i0]);
    if(iParticleType == 2){                                         // If charged PU
      if(alphaF[i0] == -100) numFCHPUis0++;                         // Count particles to ignore when taking the median
      if(alphaC[i0] == -100) numCCHPUis0++;
      alphaFCHPU[numCHPU] = alphaF[i0];                             // Only intializing particles up to the number of charged PU
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

  fPuppiParticles->Delete();

  for(Int_t i0 = 0; i0 < numCandidates; i0++){
    // Now we are going to assign the weights
    Double_t chi2 = 0;
    Double_t weight = 0;
    Int_t CandidateType = GetParticleType(fPFCandidates->At(i0));
    if(fApplyCHS && CandidateType == 1) weight = 1;                                                 // If charged PV with CHS, the weight is 1
    else if(CandidateType == 3 || (not fApplyCHS && (CandidateType == 1 || CandidateType == 2))){   // If neutral central or not CHS, get central chi2
      if(alphaC[i0] > alphaCMed) chi2 = pow((alphaC[i0] - alphaCMed),2)/sigma2C;
    }
    else if(CandidateType == 4){                                                                    // If forward particle, get forward chi2
      if(alphaF[i0] > alphaFMed) chi2 = pow((alphaF[i0] - alphaFMed),2)/sigma2F;
    }

    if(chi2 > 0){                                                    // If chi2 value was assigned
      weight = ROOT::Math::chisquared_cdf(chi2,1);                   //   then make the weight
      if(weight < fMinWeightCut) weight = 0;                         // If less than the minimum cut, set weight back to zero
    }

    if((CandidateType == 3 || CandidateType == 4) &&                 // If neutral Pt is less than expected for given NPV
       (fPFCandidates->At(i0)->Pt()*(weight) < fMinNeutralPt + fMinNeutralPtSlope * (fVertexes->GetEntries())))
      weight = 0;                                                    //   set weight to zero
    if(fInvert) weight = 1.0 - weight;                               // Invert the weight here if asked for
    if(weight == 0 && not fKeepPileup) continue;                     // Throw out if we're not keeping it

    // add PuppiParticle to the collection
    PFCandidate *PuppiParticle = fPuppiParticles->Allocate();
    new (PuppiParticle) PFCandidate(*fPFCandidates->At(i0));
    if(weight < 1)                                                   // Weight the particle if required
      PuppiParticle->SetPtEtaPhiM(PuppiParticle->Pt()*(weight),PuppiParticle->Eta(),PuppiParticle->Phi(),PuppiParticle->Mass()*(weight));
  }
  fPuppiParticles->Trim();
}
