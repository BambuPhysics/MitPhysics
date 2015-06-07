#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::JetCorrectionMod)

//--------------------------------------------------------------------------------------------------
JetCorrectionMod::JetCorrectionMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fJetsName(ModNames::gkPubJetsName),
  fCorrectedJetsName(ModNames::gkCorrectedJetsName),  
  fRhoBranchName("Rho"),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fEnabledL1Correction(kFALSE),
  rhoEtaMax(5.0),
  fJetCorrector(0),
  fEvtHdrName(Names::gkEvtHeaderBrn),
  fEventHeader(0),
  fRhoAlgo(mithep::PileupEnergyDensity::kHighEta)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
JetCorrectionMod::~JetCorrectionMod()
{
  if (fJetCorrector) {
    delete fJetCorrector;
    fJetCorrector = 0;
  }
}

void
mithep::JetCorrectionMod::SetRhoType(RhoUtilities::RhoType type)
{
  // DEPRECATED FUNCTION
  // Use SetRhoAlgo instead

  switch(type) {
  case RhoUtilities::MIT_RHO_VORONOI_LOW_ETA:
    fRhoAlgo = mithep::PileupEnergyDensity::kLowEta;
    break;
  case RhoUtilities::MIT_RHO_VORONOI_HIGH_ETA:
    fRhoAlgo = mithep::PileupEnergyDensity::kHighEta;
    break;
  case RhoUtilities::MIT_RHO_RANDOM_LOW_ETA:
    fRhoAlgo = mithep::PileupEnergyDensity::kRandomLowEta;
    break;
  case RhoUtilities::MIT_RHO_RANDOM_HIGH_ETA:
    fRhoAlgo = mithep::PileupEnergyDensity::kRandom;
    break;
  case RhoUtilities::CMS_RHO_RHOKT6PFJETS:
    fRhoAlgo = mithep::PileupEnergyDensity::kKt6PFJets;
    break;
  default:
    fRhoAlgo = mithep::PileupEnergyDensity::kHighEta;
    break;
  }
}   

//--------------------------------------------------------------------------------------------------
void JetCorrectionMod::SlaveBegin()
{
  //fill JetCorrectorParameters from files
  std::vector<JetCorrectorParameters> correctionParameters;
  for (std::vector<std::string>::const_iterator it = fCorrectionFiles.begin(); it!=fCorrectionFiles.end(); ++it) {
    correctionParameters.push_back(JetCorrectorParameters(*it));
  }
  
  //rho for L1 fastjet correction
  ReqBranch(fRhoBranchName, fRho);
  ReqBranch(fPFCandidatesName, fPFCandidates);
  ReqBranch(fEvtHdrName, fEventHeader);

  //initialize jet corrector class
  fJetCorrector = new FactorizedJetCorrector(correctionParameters);

  //keep track of which corrections are enabled
  for (std::vector<JetCorrectorParameters>::const_iterator it = correctionParameters.begin(); it != correctionParameters.end(); ++it) {
    std::string ss = it->definitions().level();
    if      (ss == "L1Offset")
      fEnabledCorrectionMask.SetBit(Jet::L1);    
    else if (ss == "L1FastJet")
      fEnabledCorrectionMask.SetBit(Jet::L1);    
    else if (ss == "L2Relative")
      fEnabledCorrectionMask.SetBit(Jet::L2);
    else if (ss == "L3Absolute")
      fEnabledCorrectionMask.SetBit(Jet::L3);
    else if (ss == "L4EMF")
      fEnabledCorrectionMask.SetBit(Jet::L4);
    else if (ss == "L5Flavor")
      fEnabledCorrectionMask.SetBit(Jet::L5);
    else if (ss == "L6SLB")
      fEnabledCorrectionMask.SetBit(Jet::L6);
    else if (ss == "L7Parton")
      fEnabledCorrectionMask.SetBit(Jet::L7);
  }

  for (UInt_t l=0; l<8; l=l+1) {
    if (fEnabledCorrectionMask.TestBit(l)) {
      fEnabledCorrections.push_back(Jet::ECorr(l));
    }
  }
}

//--------------------------------------------------------------------------------------------------
void JetCorrectionMod::Process()
{
  // Process entries of the tree. 

  const JetCol *inJets = GetObjThisEvt<JetCol>(fJetsName);
  if (!inJets) {
    SendError(kAbortModule, "Process", 
              "Pointer to input jet collection %s is null.",
              fJetsName.Data());
    return;
  }

  JetOArr *CorrectedJets = new JetOArr;
  CorrectedJets->SetOwner(kTRUE);
  CorrectedJets->SetName(fCorrectedJetsName);

  std::vector<float> corrections;

  // get the energy density from the event
  LoadBranch(fRhoBranchName);
  LoadBranch(fPFCandidatesName);
  LoadBranch(fEvtHdrName);

  // loop over jets
  for (UInt_t i=0; i<inJets->GetEntries(); ++i) {

    const Jet *inJet = inJets->At(i);

    //copy input jet, using special function to copy full derived class
    Jet *jet = inJet->MakeCopy();

    //cache uncorrected momentum
    const FourVectorM rawMom = jet->RawMom();

    //compute correction factors
    fJetCorrector->setJetEta(rawMom.Eta());
    fJetCorrector->setJetPt(rawMom.Pt());
    fJetCorrector->setJetPhi(rawMom.Phi());
    fJetCorrector->setJetE(rawMom.E());
    const PileupEnergyDensity *fR = fRho->At(0);

    Double_t theRho = fR->Rho(fRhoAlgo);

    fJetCorrector->setRho(theRho);
    fJetCorrector->setJetA(jet->JetArea());
    
    //emf only valid for CaloJets
    const CaloJet *caloJet = dynamic_cast<const CaloJet*>(jet);
    if (caloJet) {
      fJetCorrector->setJetEMF(caloJet->EnergyFractionEm());
    }
    else {
      fJetCorrector->setJetEMF(-99.0);
    }
    corrections = fJetCorrector->getSubCorrections();
    
    //set and enable correction factors in the output jet
    Double_t cumulativeCorrection = 1.0;

    for (UInt_t j=0; j<corrections.size(); ++j) {
      Double_t currentCorrection = corrections.at(j)/cumulativeCorrection;
      cumulativeCorrection = corrections.at(j);
      Jet::ECorr currentLevel = fEnabledCorrections.at(j);
      if (currentLevel==Jet::L1)
        jet->SetL1OffsetCorrectionScale(currentCorrection);
      else if (currentLevel==Jet::L2)
        jet->SetL2RelativeCorrectionScale(currentCorrection);
      else if (currentLevel==Jet::L3)
        jet->SetL3AbsoluteCorrectionScale(currentCorrection);
      else if (currentLevel==Jet::L4)
        jet->SetL4EMFCorrectionScale(currentCorrection);
      else if (currentLevel==Jet::L5)
        jet->SetL5FlavorCorrectionScale(currentCorrection);
      else if (currentLevel==Jet::L6)
        jet->SetL6LSBCorrectionScale(currentCorrection);
      else if (currentLevel==Jet::L7)
        jet->SetL7PartonCorrectionScale(currentCorrection);

      //enable corrections after setting them
      jet->EnableCorrection(currentLevel);
    }
    
    if (0) {
      printf("JetCorrectionMod(run/evt/lum/rhoEtaMax): %d %d %d %1f ",
	     fEventHeader->RunNum(),fEventHeader->EvtNum(),fEventHeader->LumiSec(),rhoEtaMax);
      printf("In Pt = %5f, Out Pt = %5f ",inJet->Pt(),jet->Pt());
      printf("In L1 = %5f, Out L1 = %5f ",
	     inJet->L1OffsetCorrectionScale(),jet->L1OffsetCorrectionScale());
      printf("In L2 = %5f, Out L2 = %5f ",
	     inJet->L2RelativeCorrectionScale(),jet->L2RelativeCorrectionScale());
      printf("In L3 = %5f, Out L3 = %5f ",
	     inJet->L3AbsoluteCorrectionScale(),jet->L3AbsoluteCorrectionScale());
      const PileupEnergyDensity *fR = fRho->At(0);
      printf("Rho(2.5) = %5f, Rho(5.0)  = %5f, Area: %5f\n",
	     fR->Rho(),fR->RhoHighEta(),jet->JetArea());
    }

    // add corrected jet to collection
    CorrectedJets->AddOwned(jet);             
  }

  // sort according to ptrootcint forward declaration data members
  CorrectedJets->Sort();
  
  // add to event for other modules to use
  AddObjThisEvt(CorrectedJets);
}


//--------------------------------------------------------------------------------------------------
void JetCorrectionMod::ApplyL1FastJetCorrection(float maxEta, bool useFixedGrid)
{
  fEnabledL1Correction = true;
  rhoEtaMax = maxEta;
  fUseFixedGrid = useFixedGrid;
}


//--------------------------------------------------------------------------------------------------
void JetCorrectionMod::ApplyL1FastJetCorrection(Jet *jet)
{
  double rho = 0;
  
  if(fUseFixedGrid) {
    // define 8 eta bins
    vector<Float_t> etabins;
    for (Int_t i=0;i<8;++i) etabins.push_back(-rhoEtaMax + 2*rhoEtaMax/7.0*i);
    //define 10 phi bins
    vector<Float_t> phibins;
    for (Int_t i=0;i<10;i++) phibins.push_back(-TMath::Pi()+(2*i+1)*TMath::TwoPi()/20.);
  
    Float_t etadist = etabins[1]-etabins[0];
    Float_t phidist = phibins[1]-phibins[0];
    Float_t etahalfdist = (etabins[1]-etabins[0])/2.;
    Float_t phihalfdist = (phibins[1]-phibins[0])/2.;

    vector<Float_t> sumPFNallSMDQ;
    sumPFNallSMDQ.reserve(80);
    for (UInt_t ieta=0;ieta<etabins.size();++ieta) {
      for (UInt_t iphi=0;iphi<phibins.size();++iphi) {
        Float_t pfniso_ieta_iphi = 0;
        assert(fPFCandidates);
        for(UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {      
          const PFCandidate *pfcand = fPFCandidates->At(i);
	  if (fabs(etabins[ieta] - pfcand->Eta()) > etahalfdist) continue;
	  if (fabs(MathUtils::DeltaPhi((Double_t)phibins[iphi],(Double_t)(pfcand->Phi()))) > phihalfdist) continue;
	  pfniso_ieta_iphi+=pfcand->Pt();
        }
        sumPFNallSMDQ.push_back(pfniso_ieta_iphi);
      }
    }

    Float_t evt_smdq = 0;
    sort(sumPFNallSMDQ.begin(),sumPFNallSMDQ.end());
  
    if(sumPFNallSMDQ.size()%2) evt_smdq = sumPFNallSMDQ[(sumPFNallSMDQ.size()-1)/2];
    else                       evt_smdq = (sumPFNallSMDQ[sumPFNallSMDQ.size()/2]+sumPFNallSMDQ[(sumPFNallSMDQ.size()-2)/2])/2.;
    rho = evt_smdq/(etadist*phidist);    
    
  } else {
    const PileupEnergyDensity *fR = fRho->At(0);

    rho = fR->Rho(fRhoAlgo);
  }
  
  Double_t l1Scale = (jet->Pt() - rho*jet->JetArea())/jet->Pt();
  l1Scale = (l1Scale>0) ? l1Scale : 0.0;
  
  jet->SetL1OffsetCorrectionScale(l1Scale);
}


//--------------------------------------------------------------------------------------------------
void JetCorrectionMod::AddCorrectionFromRelease(const std::string &path)
{
  edm::FileInPath file(path);
  AddCorrectionFromFile(file.fullPath());
}

//--------------------------------------------------------------------------------------------------
void JetCorrectionMod::AddCorrectionFromFile(const std::string &file)
{
  fCorrectionFiles.push_back(file);
}
