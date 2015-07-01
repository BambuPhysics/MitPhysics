#include <iostream>
#include <fstream>

#include "MitPhysics/Mods/interface/PuppiMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"

#include "TSystem.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/ProbFunc.h"

using namespace mithep;

ClassImp(mithep::PuppiMod)

//--------------------------------------------------------------------------------------------------
PuppiMod::PuppiMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fEtaConfigName(TString(gSystem->Getenv("CMSSW_BASE")) + "/src/MitPhysics/data/PuppiEta.cfg"),
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
  fRMSScaleFactor(1.0),
  fTrackUncertainty(1.0),
  fKeepPileup(kFALSE),
  fInvert(kFALSE),
  fApplyCHS(kTRUE),
  fApplyLowPUCorr(kTRUE),
  fUseEtaForAlgo(kFALSE),
  fEtaForAlgo(2.5)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
PuppiMod::~PuppiMod() {} 

//--------------------------------------------------------------------------------------------------
Int_t PuppiMod::GetParticleType(const PFCandidate *cand){
  if(fUseEtaForAlgo){
    Double_t checkEta = fabs(cand->Eta());
    if(checkEta > fEtaForAlgo){
      return 4;
    }
    else if(cand->PFType() != 1) return 3;
  }
  if(cand->PFType() == 1){
    if((fabs(cand->SourceVertex().Z() - fVertexes->At(0)->Position().Z()) < fDZCut) &&
       (MathUtils::AddInQuadrature(cand->SourceVertex().X() - fVertexes->At(0)->Position().X(),
				   cand->SourceVertex().Y() - fVertexes->At(0)->Position().Y()) < fD0Cut))
      return 1;                             // This is charged PV particle
    else return 2;                          // This is charged PU particle
  }
  else if(cand->PFType() != 6 && cand->PFType() != 7)
    return 3;                               // This is neutral particle in the center
  else return 4;                            // This is neutral particle in the forward region
}

//--------------------------------------------------------------------------------------------------
Int_t PuppiMod::GetEtaBin(const PFCandidate *cand){
  Int_t etaBin = -1;
  Double_t checkEta = fabs(cand->Eta());
  for(Int_t i0 = 0; i0 < fNumEtaBins; i0++){
    if(checkEta < fMaxEtas[i0]){
      etaBin = i0;                            // This relies on putting the eta bins
      break;                                  //   in increasing order
    }
  }
  return etaBin;                              // -1 means we don't care about the particle
}

//--------------------------------------------------------------------------------------------------
Double_t PuppiMod::Chi2fromDZ(Double_t dz){
  Double_t probPV = ROOT::Math::normal_cdf_c(fabs(dz),fTrackUncertainty) * 2.0;
  Double_t probPU = 1 - probPV;
  Double_t chi = 0;
  if(probPU == 1) chi = 100;                  // This doesn't exactly match CMSSW, but I like this better
  else chi = TMath::ChisquareQuantile(probPU,1);
  return pow(chi,2);
}

//--------------------------------------------------------------------------------------------------
void PuppiMod::SlaveBegin()
{
  // ===== load branches ====
  ReqBranch(fVertexesName, fVertexes);
  ReqBranch(fPFCandidatesName, fPFCandidates);

  // Prepare the storage array for the PuppiParticles
  fPuppiParticles = new PFCandidateArr(16);
  fPuppiParticles->SetName(fPuppiParticlesName);
  PublishObj(fPuppiParticles);

  // Read the configuration file
  std::cout << fEtaConfigName << std::endl;
  std::ifstream configFile;
  configFile.open(fEtaConfigName.Data());
  TString tempMaxEta;
  TString tempMinPt;
  TString tempMinNeutralPt;
  TString tempMinNeutralPtSlope;
  TString tempRMSEtaSF;
  TString tempMedEtaSF;
  TString tempEtaMaxExtrap;
  while(!configFile.eof()){
    configFile >> tempMaxEta >> tempMinPt >> tempMinNeutralPt >> tempMinNeutralPtSlope;
    configFile >> tempRMSEtaSF >> tempMedEtaSF >> tempEtaMaxExtrap;
    // If not blank or first line, store the bins
    if(tempEtaMaxExtrap != "" && tempMaxEta != "MaxEta"){
      fMaxEtas           .push_back(tempMaxEta           .Atof());
      fMinPts            .push_back(tempMinPt            .Atof());
      fMinNeutralPts     .push_back(tempMinNeutralPt     .Atof());
      fMinNeutralPtSlopes.push_back(tempMinNeutralPtSlope.Atof());
      fRMSEtaSFs         .push_back(tempRMSEtaSF         .Atof());
      fMedEtaSFs         .push_back(tempMedEtaSF         .Atof());
      fEtaMaxExtraps     .push_back(tempEtaMaxExtrap     .Atof());
    }
  }
  fNumEtaBins = fMaxEtas.size();
}

//--------------------------------------------------------------------------------------------------
void PuppiMod::SlaveTerminate()
{
  // ===== deallocate memory ====
  // Or you know, don't
}

//--------------------------------------------------------------------------------------------------
void PuppiMod::Process()
{
  // Process entries of the tree. 

  LoadBranch(fVertexesName);
  LoadBranch(fPFCandidatesName);

  // This mapper will return particles that are only close in eta, phi space
  ParticleMapper *Mapper = new ParticleMapper();
  Mapper->Initialize(*fPFCandidates,fR0,fR0,fMaxEtas[fNumEtaBins-1]);

  const Int_t numCandidates = fPFCandidates->GetEntries();
  Double_t alphaF[numCandidates];                             // This is alpha(F) of all particles for particle i
  Double_t alphaC[numCandidates];                             // This is alpha(C) of charged PV for particle i
  Int_t IndicesF[numCandidates];                              // A bunch of indices used for sorting
  Int_t IndicesC[numCandidates];
  Int_t numCHPU[fNumEtaBins];                                 // Track the number of particles that are charged PU for eta bin
  Int_t numFCHPUis0[fNumEtaBins];                             // Number of alphas to ignore when finding the median for each
  Int_t numCCHPUis0[fNumEtaBins];                             //   eta bin
  Double_t alphaFCHPU[fNumEtaBins][numCandidates];            // This is alphaF for charged PU particles for eta bin
  Double_t alphaCCHPU[fNumEtaBins][numCandidates];            // This is alphaC for charged PU particles for eta bin
  Int_t IndicesFCHPU[fNumEtaBins][numCandidates];
  Int_t IndicesCCHPU[fNumEtaBins][numCandidates];

  Int_t numCHPV[fNumEtaBins];                                 // Track the number of particles that are charged PV for eta bin
  Double_t alphaFCHPV[fNumEtaBins][numCandidates];            // This is alphaF for charged PV particles for eta bin
  Double_t alphaCCHPV[fNumEtaBins][numCandidates];            // This is alphaC for charged PV particles for eta bin

  for(Int_t i0 = 0; i0 < fNumEtaBins; i0++){                  // Initialize these for counting
    numCHPU[i0] = 0;
    numFCHPUis0[i0] = 0;
    numCCHPUis0[i0] = 0;
    numCHPV[i0] = 0;
  }

  for(Int_t i0 = 0; i0 < numCandidates; i0++){
    const PFCandidate *iCandidate = fPFCandidates->At(i0);
      
    // Initialize alphas and indices for sorting
    alphaF[i0] = 0;
    alphaC[i0] = 0;
    IndicesF[i0] = i0;
    IndicesC[i0] = i0;

    Int_t etaBin = GetEtaBin(iCandidate);
    if(etaBin < 0) continue;
    if(iCandidate->Pt() < fMinPts[etaBin]) continue;                // If not meeting Pt cut, say it's PU

    // Determine if the PFCandidate is charged PU. This will help characterize neutral PU.
    Int_t iParticleType = GetParticleType(iCandidate);
    // Mapper is getting nearby PFCandidates
    std::vector<Int_t> nearList = Mapper->GetSurrounding(i0);
    for(UInt_t i1 = 0; i1 < nearList.size(); i1++){
      Int_t i2 = nearList[i1];                                      // Getting index from mapper results
      if(i0 == i2) continue;                                        // Don't bother comparing a particle to itself
      const PFCandidate *jCandidate = fPFCandidates->At(i2);
      Double_t dRTemp = MathUtils::DeltaR(iCandidate,jCandidate);
      if(dRTemp > fRMin && dRTemp < fR0){                                           // Only look at other particles in this range
        Int_t jParticleType = GetParticleType(jCandidate);
        Double_t theAddition = (pow(jCandidate->Pt(),fAlpha))/(pow(dRTemp,fBeta));  // This is the thing we have to add inside the log
        alphaF[i0] = alphaF[i0] + theAddition;                                      // First do the sum inside the log (alphaF)
        if(jParticleType == 1)                                                      // If the particle is charged PV
          alphaC[i0] = alphaC[i0] + theAddition;                                    //   add to alphaC
      }
    }
    if(alphaF[i0] == 0) alphaF[i0] = -100;                          // Take the logs and ignore sum == 0 particles
    else alphaF[i0] = TMath::Log(alphaF[i0]);
    if(alphaC[i0] == 0) alphaC[i0] = -100;
    else alphaC[i0] = TMath::Log(alphaC[i0]);
    if(iParticleType == 2){                                         // If charged PU, we might store it in the proper eta bin
      Double_t checkEta = fabs(iCandidate->Eta());
      for(Int_t i1 = 0; i1 < fNumEtaBins; i1++){
        if(checkEta > fEtaMaxExtraps[i1]) continue;                 // If outside the binning that we are interested in, don't use CHPU
        if(alphaF[i0] == -100) numFCHPUis0[i1]++;                   // Count particles to ignore when taking the median
        if(alphaC[i0] == -100) numCCHPUis0[i1]++;
        alphaFCHPU[i1][numCHPU[i1]] = alphaF[i0];                   // Only intializing particles up to the number of charged PU
        alphaCCHPU[i1][numCHPU[i1]] = alphaC[i0];
        numCHPU[i1]++;                                              // Count the total number of charged PU particles
      }
    }
    else if(iParticleType == 1 && fApplyLowPUCorr){                 // If charged PV, and trying to correct
      Double_t checkEta = fabs(iCandidate->Eta());
      for(Int_t i1 = 0; i1 < fNumEtaBins; i1++){
        if(checkEta > fEtaMaxExtraps[i1]) continue;                 // If outside the binning that we are interested in, don't use CHPV
        alphaFCHPV[i1][numCHPV[i1]] = alphaF[i0];                   // Only intializing particles up to the number of charged PV
        alphaCCHPV[i1][numCHPV[i1]] = alphaC[i0];
        numCHPV[i1]++;                                              // Count the total number of charged PV particles
      }
    }
  }

  // Set the indices for sorting for the charged PU alphas
  for(Int_t i0 = 0; i0 < fNumEtaBins; i0++){
    for(Int_t i1 = 0; i1 < numCandidates; i1++){
      IndicesFCHPU[i0][i1] = i1;
      IndicesCCHPU[i0][i1] = i1;
    }
  }

  // Sort everything to find the median and drop particles faster
  // This gives back the Indices arrays in ordered form
  TMath::Sort(numCandidates,alphaF,IndicesF,0);
  TMath::Sort(numCandidates,alphaC,IndicesC,0);
  for(Int_t i0 = 0; i0 < fNumEtaBins; i0++){
    TMath::Sort(numCHPU[i0],alphaFCHPU[i0],IndicesFCHPU[i0],0);
    TMath::Sort(numCHPU[i0],alphaCCHPU[i0],IndicesCCHPU[i0],0);
  }

  // Now we'll find the median and sigma (left-handed RMS) squared for each event and eta bin
  Double_t alphaFMed[fNumEtaBins];
  Double_t alphaCMed[fNumEtaBins];
  Double_t sigma2F[fNumEtaBins];
  Double_t sigma2C[fNumEtaBins];

  for(Int_t i0 = 0; i0 < fNumEtaBins; i0++){
    Int_t medIndexF = (numCHPU[i0] + numFCHPUis0[i0])/2;              // These are sort of meta
    Int_t medIndexC = (numCHPU[i0] + numCCHPUis0[i0])/2;              // Just watch how they are used...
    if(numCHPU[i0] == numFCHPUis0[i0]) alphaFMed[i0] = 0;
    else if((numCHPU[i0] - numFCHPUis0[i0]) % 2 == 0) 
      alphaFMed[i0] = (alphaFCHPU[i0][IndicesFCHPU[i0][medIndexF - 1]] + alphaFCHPU[i0][IndicesFCHPU[i0][medIndexF]])/2;
    else alphaFMed[i0] = alphaFCHPU[i0][IndicesFCHPU[i0][medIndexF]];
    if(numCHPU[i0] == numCCHPUis0[i0]) alphaCMed[i0] = 0;
    else if((numCHPU[i0] - numCCHPUis0[i0]) % 2 == 0) 
      alphaCMed[i0] = (alphaCCHPU[i0][IndicesCCHPU[i0][medIndexC - 1]] + alphaCCHPU[i0][IndicesCCHPU[i0][medIndexC]])/2;
    else alphaCMed[i0] = alphaCCHPU[i0][IndicesCCHPU[i0][medIndexC]];
    
    // Now compute the sigma2s
    sigma2F[i0] = 0;
    sigma2C[i0] = 0;
    for(Int_t i1 = numFCHPUis0[i0]; i1 < medIndexF; i1++) 
      sigma2F[i0] = sigma2F[i0] + pow((alphaFMed[i0]-alphaFCHPU[i0][IndicesFCHPU[i0][i1]]),2);
    sigma2F[i0] = sigma2F[i0]/(medIndexF - numFCHPUis0[i0]);
    for(Int_t i1 = numCCHPUis0[i0]; i1 < medIndexC; i1++) 
      sigma2C[i0] = sigma2C[i0] + pow((alphaCMed[i0]-alphaCCHPU[i0][IndicesCCHPU[i0][i1]]),2);
    sigma2C[i0] = sigma2C[i0]/(medIndexC - numCCHPUis0[i0]);
    
    alphaFMed[i0] = alphaFMed[i0] * (fMedEtaSFs[i0]);                 // Scale the medians
    alphaCMed[i0] = alphaCMed[i0] * (fMedEtaSFs[i0]);

    sigma2F[i0] = sigma2F[i0] * (fRMSScaleFactor) * (fRMSEtaSFs[i0]); // Scale the sigmas
    sigma2C[i0] = sigma2C[i0] * (fRMSScaleFactor) * (fRMSEtaSFs[i0]);

    if(fApplyLowPUCorr){                                              // If we're applying LowPUCorr
      Int_t NumCorrF = 0;
      Int_t NumCorrC = 0;
      for(Int_t i1 = 0; i1 < numCHPV[i0]; i1++){                      //   count PV that are less than median
        if(alphaFCHPV[i0][i1] < alphaFMed[i0]) NumCorrF++;
        if(alphaCCHPV[i0][i1] < alphaCMed[i0]) NumCorrC++;
      }
      alphaFMed[i0] = alphaFMed[i0] 
        - sqrt(ROOT::Math::chisquared_quantile(double(NumCorrF)/
                                               double(numCHPV[i0] + numCHPU[i0] - numFCHPUis0[i0]),
                                               1.)
               *sigma2F[i0]);
      alphaCMed[i0] = alphaCMed[i0] 
        - sqrt(ROOT::Math::chisquared_quantile(double(NumCorrC)/
                                               double(numCHPV[i0] + numCHPU[i0] - numCCHPUis0[i0]),
                                               1.)
               *sigma2C[i0]);
    }
  }

  fPuppiParticles->Delete();

  for(Int_t i0 = 0; i0 < numCandidates; i0++){
    // Now we are going to assign the weights
    Double_t chi2 = 0;
    Double_t weight = 0;
    Int_t CandidateType = GetParticleType(fPFCandidates->At(i0));
    Int_t etaBin = GetEtaBin(fPFCandidates->At(i0));

    // If charged PV with CHS, the weight is 1
    if(fApplyCHS && CandidateType == 1) weight = 1;
    // If neutral central or not CHS, get central chi2
    else if(CandidateType == 3 || (not fApplyCHS && (CandidateType == 1 || CandidateType == 2)))
      {if(alphaC[i0] > alphaCMed[etaBin]) chi2 = pow((alphaC[i0] - alphaCMed[etaBin]),2)/sigma2C[etaBin];}
    // If forward particle, get forward chi2
    else if(CandidateType == 4)
      {if(alphaF[i0] > alphaFMed[etaBin]) chi2 = pow((alphaF[i0] - alphaFMed[etaBin]),2)/sigma2F[etaBin];}

    if(chi2 > 0){                                                    // If chi2 value was assigned
      weight = ROOT::Math::chisquared_cdf(chi2,1);                   //   then make the weight
      if(weight < fMinWeightCut) weight = 0;                         // If less than the minimum cut, set weight back to zero
    }

    if((CandidateType == 3 || CandidateType == 4) &&                 // If neutral Pt is less than expected for given NPV
       (fPFCandidates->At(i0)->Pt()*(weight) < fMinNeutralPts[etaBin] + fMinNeutralPtSlopes[etaBin] * (fVertexes->GetEntries())))
      weight = 0;                                                    //   set weight to zero
    if(fInvert) weight = 1.0 - weight;                               // Invert the weight here if asked for
    if(weight == 0 && not fKeepPileup) continue;                     // Throw out if we're not keeping it

    // add PuppiParticle to the collection
    PFCandidate *PuppiParticle = fPuppiParticles->Allocate();
    new (PuppiParticle) PFCandidate(*fPFCandidates->At(i0));
    if(weight < 1)                                                   // Weight the particle if required
      PuppiParticle->SetPtEtaPhiM(PuppiParticle->Pt()*(weight),PuppiParticle->Eta(),
                                  PuppiParticle->Phi(),PuppiParticle->Mass()*(weight));
  }
  fPuppiParticles->Trim();
}
