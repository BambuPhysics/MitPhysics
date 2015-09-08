#include <iostream>
#include <fstream>

#include "MitPhysics/Mods/interface/PuppiMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"

#include "TSystem.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/ProbFunc.h"

using namespace mithep;

ClassImp(mithep::PuppiMod)

//--------------------------------------------------------------------------------------------------
PuppiMod::PuppiMod(const char *name, const char *title) :
  BaseMod(name,title),
  fEtaConfigName(""),
  fVertexesName(Names::gkPVBrn),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPuppiParticlesName("PuppiParticles"),
  fInvertedParticlesName("PUPuppiParticles"),
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
  fNoLepton(kTRUE),
  fKeepPileup(kFALSE),
  fInvert(kFALSE),
  fApplyCHS(kTRUE),
  fApplyLowPUCorr(kTRUE),
  fUseEtaForAlgo(kTRUE),
  fEtaForAlgo(2.5),
  fBothPVandPU(kFALSE),
  fDumpingPuppi(kFALSE),
  fMapper(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
PuppiMod::~PuppiMod()
{
}

//--------------------------------------------------------------------------------------------------
PuppiMod::ParticleType
PuppiMod::GetParticleType(const PFCandidate *cand, Vertex const* pv) const
{
  auto pfType = cand->PFType();
  Bool_t charged = (pfType == PFCandidate::eHadron || pfType == PFCandidate::eElectron || pfType == PFCandidate::eMuon);

  if (fUseEtaForAlgo) {
    Double_t checkEta = cand->AbsEta();
    if (checkEta > fEtaForAlgo)
      return kNeutralForward;
    else if (not charged)
      return kNeutralCentral;
  }
  if (charged) {
    Double_t checkDZ = 0.;
    Double_t checkD0 = 0.;
    if (cand->HasTrackerTrk()) {
      checkDZ = cand->TrackerTrk()->DzCorrected(*pv);
      checkD0 = cand->TrackerTrk()->D0Corrected(*pv);
    }
    else if (cand->HasGsfTrk()) {
      checkDZ = cand->GsfTrk()->DzCorrected(*pv);
      checkD0 = cand->GsfTrk()->D0Corrected(*pv);
    }
    if (fabs(checkDZ) < fDZCut && checkD0 < fD0Cut)
      return kChargedPrimary;    // This is charged PV particle
    else
      return kChargedPU;         // This is charged PU particle
  }
  else if (pfType != PFCandidate::eHadronHF && pfType != PFCandidate::eEGammaHF)
    return kNeutralCentral;      // This is neutral particle in the center
  else
    return kNeutralForward;      // This is neutral particle in the forward region
}

//--------------------------------------------------------------------------------------------------
Int_t
PuppiMod::GetEtaBin(const PFCandidate *cand) const
{
  Int_t etaBin = -1;
  Double_t checkEta = cand->AbsEta();
  for (UInt_t iEtaBin = 0; iEtaBin < fNumEtaBins; iEtaBin++) {
    if (checkEta < fMaxEtas[iEtaBin]) {
      etaBin = iEtaBin;                       // This relies on putting the eta bins
      break;                                  //   in increasing order
    }
  }
  return etaBin;                              // -1 means we don't care about the particle
}

//--------------------------------------------------------------------------------------------------
Double_t
PuppiMod::Chi2fromDZ(Double_t dz) const
{
  Double_t probPV = ROOT::Math::normal_cdf_c(fabs(dz),fTrackUncertainty) * 2.0;
  Double_t probPU = 1 - probPV;
  Double_t chi = 0;
  if (probPU == 1)
    chi = 100;                                // This doesn't exactly match CMSSW, but I like this better
  else
    chi = TMath::ChisquareQuantile(probPU,1);
  return chi * chi;
}

//--------------------------------------------------------------------------------------------------
void
PuppiMod::SlaveBegin()
{
  // Read the configuration file
  if (fEtaConfigName == "") {
    SendError(kAbortAnalysis, "SlaveBegin", "Config name for Puppi not given.");
    return;
  }

  // Prepare the storage array for the PuppiParticles
  fPuppiParticles = new PFCandidateArr(16);
  fPuppiParticles->SetName(fPuppiParticlesName);
  PublishObj(fPuppiParticles);

  if (fBothPVandPU) {
    fInvertedParticles = new PFCandidateArr(16);
    fInvertedParticles->SetName(fInvertedParticlesName);
    PublishObj(fInvertedParticles);
  }

  Info("SlaveBegin", "Reading PUPPI config file " + fEtaConfigName);

  std::ifstream configFile;
  configFile.open(fEtaConfigName.Data());
  if (!configFile.is_open()) {
    SendError(kAbortAnalysis, "SlaveBegin", "Config file %s could not be opened.", fEtaConfigName.Data());
    return;
  }

  TString tempMaxEta;
  TString tempMinPt;
  TString tempMinNeutralPt;
  TString tempMinNeutralPtSlope;
  TString tempRMSEtaSF;
  TString tempMedEtaSF;
  TString tempEtaMaxExtrap;
  while(!configFile.eof()) {
    configFile >> tempMaxEta >> tempMinPt >> tempMinNeutralPt >> tempMinNeutralPtSlope;
    configFile >> tempRMSEtaSF >> tempMedEtaSF >> tempEtaMaxExtrap;
    // if not blank or first line, store the bins
    if (tempEtaMaxExtrap != "" && tempMaxEta != "MaxEta") {
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
  if (fDumpingPuppi)
    Info("SlaveBegin", "Number of Eta bins: %d", fNumEtaBins);

  fMapper = new ParticleMapper(fR0, fR0, fMaxEtas[fNumEtaBins-1]);
}

//--------------------------------------------------------------------------------------------------
void
PuppiMod::SlaveTerminate()
{
  // ===== deallocate memory ====
  RetractObj(fPuppiParticles->GetName());
  delete fPuppiParticles;
  delete fMapper;
  if (fBothPVandPU) {
    RetractObj(fInvertedParticles->GetName());
    delete fInvertedParticles;
  }
}

//--------------------------------------------------------------------------------------------------
void
PuppiMod::Process()
{
  // Process entries of the tree.
  auto* vertexes = GetObject<VertexCol>(fVertexesName);
  if (!vertexes)
    SendError(kAbortAnalysis, "Process", "PV collection not found");

  auto* pv = vertexes->At(0);

  auto* inPfCandidates = GetObject<PFCandidateCol>(fPFCandidatesName);
  PFCandidateCol* pfCandidates = NULL;
  PFCandidateOArr leptonCandidates;
  PFCandidateOArr tempCandidates;

  if (fNoLepton) {                                                                              // This is the default behavior of Puppi
    for (UInt_t iCand = 0; iCand < inPfCandidates->GetEntries(); iCand++) { 
      auto* cand = inPfCandidates->At(iCand);
      if (cand->PFType() != PFCandidate::eMuon && cand->PFType() != PFCandidate::eElectron)     // Take out leptons
        tempCandidates.Add(cand);
      else
        leptonCandidates.Add(cand);                                                             // Store them to add back at the end
    }
    pfCandidates = &tempCandidates;
  }
  else
    pfCandidates = inPfCandidates;

  if (fDumpingPuppi) {
    std::cout << "Primary Vertex location: " << pv->Position().x() << ", ";
    std::cout << pv->Position().y() << ", " << pv->Position().z() << std::endl;
  }

  // This mapper will return particles that are only close in eta, phi space
  fMapper->InitEvent(*pfCandidates);

  UInt_t numCandidates = pfCandidates->GetEntries();
  std::vector<Double_t> alphaF(numCandidates, 0.);                       // This is alpha(F) of all particles for particle i
  std::vector<Double_t> alphaC(numCandidates, 0.);                       // This is alpha(C) of charged PV for particle i
  std::vector<UInt_t> numFCHPUis0(fNumEtaBins, 0);                       // Number of alphas to ignore when finding the median
  std::vector<UInt_t> numCCHPUis0(fNumEtaBins, 0);                       //   for each eta bin
  std::vector<std::vector<Double_t>> alphaFCHPU(fNumEtaBins);            // This is alphaF for charged PU particles for eta bin
  std::vector<std::vector<Double_t>> alphaCCHPU(fNumEtaBins);            // This is alphaC for charged PU particles for eta bin
  std::vector<std::vector<Double_t>> alphaFCHPV(fNumEtaBins);            // This is alphaF for charged PV particles for eta bin
  std::vector<std::vector<Double_t>> alphaCCHPV(fNumEtaBins);            // This is alphaC for charged PV particles for eta bin

  for (UInt_t iCand = 0; iCand < numCandidates; iCand++) {
    const PFCandidate *iCandidate = pfCandidates->At(iCand);

    Int_t etaBin = GetEtaBin(iCandidate);
    if (etaBin < 0)
      continue;
    if (iCandidate->Pt() < fMinPts[etaBin])
      continue;                                                          // if not meeting Pt cut, say it's PU

    // Determine if the PFCandidate is charged PU. This will help characterize neutral PU.
    ParticleType iParticleType = GetParticleType(iCandidate, pv);

    // Mapper is getting nearby PFCandidates
    std::vector<UInt_t> nearList(fMapper->GetSurrounding(iCand));
    for (UInt_t iNeighbor : nearList) {
      if (iCand == iNeighbor)
        continue;                                                                    // Don't bother comparing a particle to itself

      const PFCandidate *jCandidate = pfCandidates->At(iNeighbor);
      Double_t dRTemp = MathUtils::DeltaR(iCandidate,jCandidate);
      if (dRTemp > fRMin && dRTemp < fR0) {                                          // Only look at other particles in this range
        ParticleType jParticleType = GetParticleType(jCandidate, pv);
        Double_t theAddition = (pow(jCandidate->Pt(),fAlpha))/(pow(dRTemp,fBeta));   // This is the thing we have to add inside the log
        alphaF[iCand] += theAddition;                                                // First do the sum inside the log (alphaF)
        if (jParticleType == kChargedPrimary)                                        // if the particle is charged PV
          alphaC[iCand] += theAddition;                                              //   add to alphaC
      }
    }
    
    if (alphaF[iCand] == 0)
      alphaF[iCand] = -100;                                                           // Take the logs and ignore sum == 0 particles
    else
      alphaF[iCand] = TMath::Log(alphaF[iCand]);

    if (alphaC[iCand] == 0)
      alphaC[iCand] = -100;
    else
      alphaC[iCand] = TMath::Log(alphaC[iCand]);

    if (iParticleType == kChargedPU) {                                  // if charged PU, we might store it in the proper eta bin
      Double_t checkEta = fabs(iCandidate->Eta());
      for (UInt_t iEtaBin = 0; iEtaBin < fNumEtaBins; iEtaBin++) {
        if (checkEta > fEtaMaxExtraps[iEtaBin])
          continue;                                                     // if outside the binning that we are interested in, don't use CHPU

        if (alphaF[iCand] <= 0)
          numFCHPUis0[iEtaBin]++;                                       // Count particles to ignore when taking the median

        if (alphaC[iCand] <= 0)
          numCCHPUis0[iEtaBin]++;

        alphaFCHPU[iEtaBin].push_back(alphaF[iCand]);                   // Only intializing particles up to the number of charged PU
        alphaCCHPU[iEtaBin].push_back(alphaC[iCand]);
      }
    }
    else if (iParticleType == kChargedPrimary && fApplyLowPUCorr) {     // if charged PV, and trying to correct
      Double_t checkEta = fabs(iCandidate->Eta());
      for (UInt_t iEtaBin = 0; iEtaBin < fNumEtaBins; iEtaBin++) {
        if (checkEta > fEtaMaxExtraps[iEtaBin])
          continue;                                                     // if outside the binning that we are interested in, don't use CHPV

        alphaFCHPV[iEtaBin].push_back(alphaF[iCand]);                   // Only intializing particles up to the number of charged PV
        alphaCCHPV[iEtaBin].push_back(alphaC[iCand]);
      }
    }
  }

  std::vector<UInt_t> numCHPU(fNumEtaBins);
  std::vector<UInt_t> numCHPV(fNumEtaBins);
  for (unsigned iEtaBin = 0; iEtaBin != fNumEtaBins; ++iEtaBin) {
    numCHPU[iEtaBin] = alphaFCHPU[iEtaBin].size();
    numCHPV[iEtaBin] = alphaFCHPV[iEtaBin].size();
  }

  // Sort alphas in ascending order to find the median and drop particles faster
  for (UInt_t iEtaBin = 0; iEtaBin < fNumEtaBins; iEtaBin++) {
    std::sort(alphaFCHPU[iEtaBin].begin(), alphaFCHPU[iEtaBin].end());
    std::sort(alphaCCHPU[iEtaBin].begin(), alphaCCHPU[iEtaBin].end());
    if (fDumpingPuppi) {
      std::cout << "alphaCs for median: " << std::endl;
      for (Double_t alC : alphaCCHPU[iEtaBin])
        std::cout << alC << std::endl;
      std::cout << "alphaFs for median: " << std::endl;
      for (Double_t alF : alphaFCHPU[iEtaBin])
        std::cout << alF << std::endl;
    }
  }

  // Now we'll find the median and sigma (left-handed RMS) squared for each event and eta bin
  std::vector<Double_t> alphaFMed(fNumEtaBins);
  std::vector<Double_t> alphaCMed(fNumEtaBins);
  std::vector<Double_t> sigma2F(fNumEtaBins, 0.);
  std::vector<Double_t> sigma2C(fNumEtaBins, 0.);

  for (UInt_t iEtaBin = 0; iEtaBin < fNumEtaBins; iEtaBin++) {
    UInt_t medIndexF = (numCHPU[iEtaBin] + numFCHPUis0[iEtaBin])/2;
    UInt_t medIndexC = (numCHPU[iEtaBin] + numCCHPUis0[iEtaBin])/2;
    if (fDumpingPuppi)
      std::cout << "In bin " << iEtaBin << " using " << medIndexC << "; " << medIndexF << std::endl;

    if (numCHPU[iEtaBin] == numFCHPUis0[iEtaBin])
      alphaFMed[iEtaBin] = 0;
    else if ((numCHPU[iEtaBin] - numFCHPUis0[iEtaBin]) % 2 == 0)
      alphaFMed[iEtaBin] = (alphaFCHPU[iEtaBin][medIndexF - 1] + alphaFCHPU[iEtaBin][medIndexF])/2;
    else
      alphaFMed[iEtaBin] = alphaFCHPU[iEtaBin][medIndexF];

    if (numCHPU[iEtaBin] == numCCHPUis0[iEtaBin])
      alphaCMed[iEtaBin] = 0;
    else if ((numCHPU[iEtaBin] - numCCHPUis0[iEtaBin]) % 2 == 0)
      alphaCMed[iEtaBin] = (alphaCCHPU[iEtaBin][medIndexC - 1] + alphaCCHPU[iEtaBin][medIndexC])/2;
    else
      alphaCMed[iEtaBin] = alphaCCHPU[iEtaBin][medIndexC];

    // Now compute the sigma2s
    for (UInt_t iAlpha = numFCHPUis0[iEtaBin]; iAlpha < medIndexF; iAlpha++)
      sigma2F[iEtaBin] += pow((alphaFMed[iEtaBin]-alphaFCHPU[iEtaBin][iAlpha]),2);

    if (medIndexF != numFCHPUis0[iEtaBin])
      sigma2F[iEtaBin] /= (medIndexF - numFCHPUis0[iEtaBin]);

    for (UInt_t iAlpha = numCCHPUis0[iEtaBin]; iAlpha < medIndexC; iAlpha++)
      sigma2C[iEtaBin] += pow((alphaCMed[iEtaBin]-alphaCCHPU[iEtaBin][iAlpha]),2);

    if (medIndexF != numFCHPUis0[iEtaBin])
      sigma2C[iEtaBin] /= (medIndexC - numCCHPUis0[iEtaBin]);

    alphaFMed[iEtaBin] *= fMedEtaSFs[iEtaBin];                       // Scale the medians
    alphaCMed[iEtaBin] *= fMedEtaSFs[iEtaBin];

    sigma2F[iEtaBin] *= fRMSScaleFactor * fRMSEtaSFs[iEtaBin];       // Scale the sigmas
    sigma2C[iEtaBin] *= fRMSScaleFactor * fRMSEtaSFs[iEtaBin];

    if (fApplyLowPUCorr) {                                           // if we're applying LowPUCorr
      Double_t NumCorrF = 0.;
      Double_t NumCorrC = 0.;
      for (UInt_t iCHPV = 0; iCHPV < numCHPV[iEtaBin]; iCHPV++) {    //   count PV that are less than median
        if (alphaFCHPV[iEtaBin][iCHPV] < alphaFMed[iEtaBin])
          NumCorrF += 1.;
        if (alphaCCHPV[iEtaBin][iCHPV] < alphaCMed[iEtaBin])
          NumCorrC += 1.;
      }

      Double_t quantF = ROOT::Math::chisquared_quantile(NumCorrF / double(numCHPV[iEtaBin] + numCHPU[iEtaBin] - numFCHPUis0[iEtaBin]), 1.);
      Double_t quantC = ROOT::Math::chisquared_quantile(NumCorrC / double(numCHPV[iEtaBin] + numCHPU[iEtaBin] - numFCHPUis0[iEtaBin]), 1.);
      alphaFMed[iEtaBin] -= sqrt(quantF * sigma2F[iEtaBin]);
      alphaCMed[iEtaBin] -= sqrt(quantC * sigma2C[iEtaBin]);
    }
  }

  fPuppiParticles->Delete();
  if (fBothPVandPU)
    fInvertedParticles->Delete();

  for (UInt_t iCand = 0; iCand < numCandidates; iCand++) {
    auto* cand = pfCandidates->At(iCand);

    // Now we are going to assign the weights
    Double_t chi2 = 0;
    Double_t weight = 0;
    ParticleType CandidateType = GetParticleType(cand, pv);
    Int_t etaBin = GetEtaBin(cand);

    // if charged PV with CHS, the weight is 1
    if (fApplyCHS && CandidateType == kChargedPrimary)
      weight = 1;
    // if neutral central or not CHS, get central chi2
    else if (CandidateType == kNeutralCentral || (!fApplyCHS && (CandidateType == kChargedPrimary || CandidateType == kChargedPU))) {
      if (alphaC[iCand] > alphaCMed[etaBin])
        chi2 = pow((alphaC[iCand] - alphaCMed[etaBin]),2)/sigma2C[etaBin];
    }
    // if for ward particle, get for ward chi2
    else if (CandidateType == kNeutralForward) {
      if (alphaF[iCand] > alphaFMed[etaBin])
        chi2 = pow((alphaF[iCand] - alphaFMed[etaBin]),2)/sigma2F[etaBin];
    }

    if (chi2 > 0) {                                                    // if chi2 value was assigned
      weight = ROOT::Math::chisquared_cdf(chi2,1);                     //   then make the weight
      if (weight < fMinWeightCut)
        weight = 0;                                                    // if less than the minimum cut, set weight back to zero
    }

    if ((CandidateType == kNeutralCentral || CandidateType == kNeutralForward) &&   // if neutral Pt is less than expected for given NPV
       (cand->Pt() * weight < fMinNeutralPts[etaBin] + fMinNeutralPtSlopes[etaBin] * (vertexes->GetEntries())))
      weight = 0;                                                                   //   set weight to zero

    if (fInvert && !fBothPVandPU)
      weight = 1.0 - weight;                                                        // Invert the weight here if asked for 

    if (weight != 0 || fKeepPileup) {
      // add PuppiParticle to the collection
      PFCandidate *PuppiParticle = fPuppiParticles->Allocate();
      new (PuppiParticle) PFCandidate(*cand);
      if (fDumpingPuppi) {
        std::cout << "=========================================================================================" << std::endl;
        std::cout << "PF Candidate Number: " << iCand << std::endl;
        std::cout << "PF Type: " << GetParticleType(PuppiParticle, pv) << " (" << PuppiParticle->PFType() << ")" << std::endl;
        if (PuppiParticle->HasTrackerTrk()) {
          std::cout << "Vertex distances: " << PuppiParticle->TrackerTrk()->DzCorrected(*pv) << "; ";
          std::cout << PuppiParticle->TrackerTrk()->D0Corrected(*pv) << std::endl;
        }
        else if (PuppiParticle->HasGsfTrk()) {
          std::cout << "Vertex distances: " << PuppiParticle->GsfTrk()->DzCorrected(*pv) << "; ";
          std::cout << PuppiParticle->GsfTrk()->D0Corrected(*pv) << std::endl;
        }
        std::vector<UInt_t> nearList(fMapper->GetSurrounding(iCand));
        std::cout << "Nearby PV particles: ";
        for (UInt_t iNeighbor : nearList) {
          if (iNeighbor == iCand)
            continue;
          
          if (GetParticleType(pfCandidates->At(iNeighbor), pv) == kChargedPrimary)
            std::cout << iNeighbor << " (" << MathUtils::DeltaR(PuppiParticle->Eta(),PuppiParticle->Phi(),
                                                                pfCandidates->At(iNeighbor)->Eta(),
                                                                pfCandidates->At(iNeighbor)->Phi()) <<  "), ";
        }
        std::cout << std::endl;
        std::cout << "Nearby other particles: ";
        for (UInt_t iNeighbor : nearList) {
          if (iNeighbor == iCand)
            continue;
          
          if (GetParticleType(pfCandidates->At(iNeighbor), pv) != kChargedPrimary)
            std::cout << iNeighbor << " (" << MathUtils::DeltaR(PuppiParticle->Eta(),PuppiParticle->Phi(),
                                                                pfCandidates->At(iNeighbor)->Eta(),
                                                                pfCandidates->At(iNeighbor)->Phi()) <<  "), ";
        }
        std::cout << std::endl;
        
        std::cout << "Pt: " << PuppiParticle->Pt() << "; Eta: " << PuppiParticle->Eta();
        std::cout << "; Phi: " << PuppiParticle->Phi() << "; Mass: " << PuppiParticle->Mass() << std::endl;
        std::cout << "Median Alphas:";
        std::cout << alphaCMed[GetEtaBin(PuppiParticle)] << "; " << alphaFMed[GetEtaBin(PuppiParticle)] << std::endl;
        std::cout << "Alpha RMS:";
        std::cout << sigma2C[GetEtaBin(PuppiParticle)] << "; " << sigma2F[GetEtaBin(PuppiParticle)] << std::endl;
        std::cout << "Weight: " << weight << "; " << alphaC[iCand] << "; " << alphaF[iCand] << std::endl;
      }
      
      if (weight < 1)                                                    // Weight the particle if required
        PuppiParticle->SetPtEtaPhiM(PuppiParticle->Pt() * weight,PuppiParticle->Eta(),
                                    PuppiParticle->Phi(),PuppiParticle->Mass() * weight);
    }
    if (fBothPVandPU) {
      weight = 1.0 - weight;
      if (weight == 0)
        continue;
      PFCandidate *PuppiParticle = fInvertedParticles->Allocate();
      new (PuppiParticle) PFCandidate(*pfCandidates->At(iCand));
      if (weight < 1)                                                    // Weight the particle if required
        PuppiParticle->SetPtEtaPhiM(PuppiParticle->Pt() * weight,PuppiParticle->Eta(),
                                    PuppiParticle->Phi(),PuppiParticle->Mass() * weight);
    }
  }

  if (fNoLepton) {                                                     // If no leptons were used, add them back in
    for (UInt_t iLepton = 0; iLepton < leptonCandidates.GetEntries(); iLepton++) {
      auto* cand = leptonCandidates.At(iLepton);
      if (GetParticleType(cand, pv) != kChargedPU) {                   // From PV, simply add with weight 1
        PFCandidate *PuppiParticle = fPuppiParticles->Allocate();
        new (PuppiParticle) PFCandidate(*cand);
      }
      else if (fBothPVandPU) {                                         // If PU, but using Inverted Puppi, add to inverted
        PFCandidate *PuppiParticle = fInvertedParticles->Allocate();
        new (PuppiParticle) PFCandidate(*cand);
      }
      else if (fKeepPileup) {                                          // If keeping pileup, just add with weight 0
        PFCandidate *PuppiParticle = fPuppiParticles->Allocate();
        new (PuppiParticle) PFCandidate(*cand);
        PuppiParticle->SetPtEtaPhiM(0,PuppiParticle->Eta(),
                                    PuppiParticle->Phi(),0);
      }
    }
  }

  if (fDumpingPuppi)
    fDumpingPuppi = false;                                             // You will not want to dump this more than once

  if (fBothPVandPU)
    fInvertedParticles->Trim();
  
  fPuppiParticles->Trim();
}
