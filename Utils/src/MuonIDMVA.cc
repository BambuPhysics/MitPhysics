#include "MitPhysics/Utils/interface/MuonIDMVA.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>
#include <TRandom3.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

ClassImp(mithep::MuonIDMVA)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
MuonIDMVA::MuonIDMVA() :
  fMethodname("BDTG method"),
  fIsInitialized(kFALSE),
  fMVAType(MuonIDMVA::kUninitialized),
  fUseBinnedVersion(kTRUE),
  fNMVABins(0),
  fRhoAlgo(mithep::PileupEnergyDensity::kHighEta)
{
}


//--------------------------------------------------------------------------------------------------
MuonIDMVA::~MuonIDMVA()
{
  for(UInt_t i=0; i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i])
      delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void MuonIDMVA::Initialize(std::string const& methodName,
                           std::string const& weightsfile,
                           MuonIDMVA::MVAType type,
			   UInt_t algo)
{
  
  std::vector<std::string> tempWeightFileVector;
  tempWeightFileVector.push_back(weightsfile);
  Initialize(methodName,type,kFALSE,tempWeightFileVector,algo);
}

//--------------------------------------------------------------------------------------------------
void MuonIDMVA::Initialize(TString const& methodName,
                           TString const& Subdet0Pt10To20Weights , 
                           TString const& Subdet1Pt10To20Weights , 
                           TString const& Subdet2Pt10To20Weights,
                           TString const& Subdet0Pt20ToInfWeights,
                           TString const& Subdet1Pt20ToInfWeights, 
                           TString const& Subdet2Pt20ToInfWeights,
                           MuonIDMVA::MVAType type,
			   UInt_t algo) {

  std::vector<std::string> tempWeightFileVector;
  tempWeightFileVector.push_back(std::string(Subdet0Pt10To20Weights.Data()));
  tempWeightFileVector.push_back(std::string(Subdet1Pt10To20Weights.Data()));
  tempWeightFileVector.push_back(std::string(Subdet2Pt10To20Weights.Data()));
  tempWeightFileVector.push_back(std::string(Subdet0Pt20ToInfWeights.Data()));
  tempWeightFileVector.push_back(std::string(Subdet1Pt20ToInfWeights.Data()));
  tempWeightFileVector.push_back(std::string(Subdet2Pt20ToInfWeights.Data()));
  Initialize(std::string(methodName.Data()),type,kTRUE,tempWeightFileVector,algo);

}


//--------------------------------------------------------------------------------------------------
void MuonIDMVA::Initialize(std::string const& methodName,
                           MuonIDMVA::MVAType type,
                           Bool_t useBinnedVersion,
                           std::vector<std::string> const& weightsfiles,
			   UInt_t algo) {
  
  //clean up first
  for (uint i=0;i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
  fTMVAReader.clear();


  //initialize
  fIsInitialized = kTRUE;
  fMethodname = methodName;
  fMVAType = type;
  fUseBinnedVersion = useBinnedVersion;
  fRhoAlgo = algo;

  //Define expected number of bins
  UInt_t ExpectedNBins = 0;
  if (!fUseBinnedVersion) {
    ExpectedNBins = 1;
  } else if    (type == kV2 
                || type == kV3
                || type == kV8
                || type == kIDIsoCombinedDetIso
                || type == kIsoRingsV0
                || type == kIDV0 
                || type == kIDIsoCombinedIsoRingsV0
                ) {
    ExpectedNBins = 6;
  } else if (type == kIsoDeltaR){
    ExpectedNBins = 4;
  }

  fNMVABins = ExpectedNBins;

  //Check number of weight files given
  if (fNMVABins != weightsfiles.size() ) {
    std::cout << "Error: Expected Number of bins = " << fNMVABins << " does not equal to weightsfiles.size() = " 
              << weightsfiles.size() << std::endl;
    assert(fNMVABins == weightsfiles.size());
  }

  for(UInt_t i=0; i<fNMVABins; ++i) {
    TMVA::Reader *tmpTMVAReader = new TMVA::Reader( "!Color:!Silent:Error" );  
    tmpTMVAReader->SetVerbose(kTRUE);

    if (type == kV2) {
      tmpTMVAReader->AddVariable( "TkNchi2",              &fMVAVar_MuTkNchi2               );
      tmpTMVAReader->AddVariable( "GlobalNchi2",          &fMVAVar_MuGlobalNchi2           );
      tmpTMVAReader->AddVariable( "NValidHits",           &fMVAVar_MuNValidHits            );
      tmpTMVAReader->AddVariable( "NTrackerHits",         &fMVAVar_MuNTrackerHits          );
      tmpTMVAReader->AddVariable( "NPixelHits",           &fMVAVar_MuNPixelHits            );
      tmpTMVAReader->AddVariable( "NMatches",             &fMVAVar_MuNMatches              );
      tmpTMVAReader->AddVariable( "D0",                   &fMVAVar_MuD0                    );      
      tmpTMVAReader->AddVariable( "IP3d",                 &fMVAVar_MuIP3d                  );      
      tmpTMVAReader->AddVariable( "IP3dSig",              &fMVAVar_MuIP3dSig               );      
      tmpTMVAReader->AddVariable( "TrkKink",              &fMVAVar_MuTrkKink               );      
      tmpTMVAReader->AddVariable( "SegmentCompatibility", &fMVAVar_MuSegmentCompatibility  ); 
    }

    if (type == kV3) {
      tmpTMVAReader->AddVariable( "TkNchi2",              &fMVAVar_MuTkNchi2               );
      tmpTMVAReader->AddVariable( "GlobalNchi2",          &fMVAVar_MuGlobalNchi2           );
      tmpTMVAReader->AddVariable( "NValidHits",           &fMVAVar_MuNValidHits            );
      tmpTMVAReader->AddVariable( "NTrackerHits",         &fMVAVar_MuNTrackerHits          );
      tmpTMVAReader->AddVariable( "NPixelHits",           &fMVAVar_MuNPixelHits            );
      tmpTMVAReader->AddVariable( "NMatches",             &fMVAVar_MuNMatches              );
      tmpTMVAReader->AddVariable( "D0",                   &fMVAVar_MuD0                    );      
      tmpTMVAReader->AddVariable( "IP3d",                 &fMVAVar_MuIP3d                  );      
      tmpTMVAReader->AddVariable( "IP3dSig",              &fMVAVar_MuIP3dSig               );      
      tmpTMVAReader->AddVariable( "TrkKink",              &fMVAVar_MuTrkKink               );      
      tmpTMVAReader->AddVariable( "SegmentCompatibility", &fMVAVar_MuSegmentCompatibility  );      
      tmpTMVAReader->AddVariable( "CaloCompatibility",    &fMVAVar_MuCaloCompatibility     );      
      tmpTMVAReader->AddVariable( "HadEnergyOverPt",      &fMVAVar_MuHadEnergyOverPt       );      
      if (i==0 || i==2 || i==4) {
        tmpTMVAReader->AddVariable( "HoEnergyOverPt",     &fMVAVar_MuHoEnergyOverPt        );      
      }
      tmpTMVAReader->AddVariable( "EmEnergyOverPt",       &fMVAVar_MuEmEnergyOverPt        );      
      tmpTMVAReader->AddVariable( "HadS9EnergyOverPt",    &fMVAVar_MuHadS9EnergyOverPt     );      
      if (i==0 || i==2 || i==4) {
        tmpTMVAReader->AddVariable( "HoS9EnergyOverPt",   &fMVAVar_MuHoS9EnergyOverPt      );      
      }
      tmpTMVAReader->AddVariable( "EmS9EnergyOverPt",     &fMVAVar_MuEmS9EnergyOverPt      );      
    }

    if (type == kV8) {
      tmpTMVAReader->AddVariable( "TkNchi2",              &fMVAVar_MuTkNchi2               );
      tmpTMVAReader->AddVariable( "GlobalNchi2",          &fMVAVar_MuGlobalNchi2           );
      tmpTMVAReader->AddVariable( "NValidHits",           &fMVAVar_MuNValidHits            );
      tmpTMVAReader->AddVariable( "NTrackerHits",         &fMVAVar_MuNTrackerHits          );
      tmpTMVAReader->AddVariable( "NPixelHits",           &fMVAVar_MuNPixelHits            );
      tmpTMVAReader->AddVariable( "NMatches",             &fMVAVar_MuNMatches              );
      tmpTMVAReader->AddVariable( "D0",                   &fMVAVar_MuD0                    );      
      tmpTMVAReader->AddVariable( "IP3d",                 &fMVAVar_MuIP3d                  );      
      tmpTMVAReader->AddVariable( "IP3dSig",              &fMVAVar_MuIP3dSig               );      
      tmpTMVAReader->AddVariable( "TrkKink",              &fMVAVar_MuTrkKink               );      
      tmpTMVAReader->AddVariable( "SegmentCompatibility", &fMVAVar_MuSegmentCompatibility  );      
      tmpTMVAReader->AddVariable( "CaloCompatibility",    &fMVAVar_MuCaloCompatibility     );      
      tmpTMVAReader->AddVariable( "HadEnergyOverPt",      &fMVAVar_MuHadEnergyOverPt       );      
      tmpTMVAReader->AddVariable( "EmEnergyOverPt",       &fMVAVar_MuEmEnergyOverPt        );      
      tmpTMVAReader->AddVariable( "HadS9EnergyOverPt",    &fMVAVar_MuHadS9EnergyOverPt     );      
      tmpTMVAReader->AddVariable( "EmS9EnergyOverPt",     &fMVAVar_MuEmS9EnergyOverPt      );      
      tmpTMVAReader->AddVariable( "ChargedIso03OverPt",   &fMVAVar_MuChargedIso03OverPt    );
      tmpTMVAReader->AddVariable( "NeutralIso03OverPt",   &fMVAVar_MuNeutralIso03OverPt    );      
      tmpTMVAReader->AddVariable( "ChargedIso04OverPt",   &fMVAVar_MuChargedIso04OverPt    );
      tmpTMVAReader->AddVariable( "NeutralIso04OverPt",   &fMVAVar_MuNeutralIso04OverPt    );      
    }
    
    if (type == kIDIsoCombinedDetIso) {
      tmpTMVAReader->AddVariable( "TkNchi2",              &fMVAVar_MuTkNchi2               );
      tmpTMVAReader->AddVariable( "GlobalNchi2",          &fMVAVar_MuGlobalNchi2           );
      tmpTMVAReader->AddVariable( "NValidHits",           &fMVAVar_MuNValidHits            );
      tmpTMVAReader->AddVariable( "NTrackerHits",         &fMVAVar_MuNTrackerHits          );
      tmpTMVAReader->AddVariable( "NPixelHits",           &fMVAVar_MuNPixelHits            );
      tmpTMVAReader->AddVariable( "NMatches",             &fMVAVar_MuNMatches              );
      tmpTMVAReader->AddVariable( "D0",                   &fMVAVar_MuD0                    );      
      tmpTMVAReader->AddVariable( "IP3d",                 &fMVAVar_MuIP3d                  );      
      tmpTMVAReader->AddVariable( "IP3dSig",              &fMVAVar_MuIP3dSig               );      
      tmpTMVAReader->AddVariable( "TrkKink",              &fMVAVar_MuTrkKink               );      
      tmpTMVAReader->AddVariable( "SegmentCompatibility", &fMVAVar_MuSegmentCompatibility  );      
      tmpTMVAReader->AddVariable( "CaloCompatibility",    &fMVAVar_MuCaloCompatibility     );      
      tmpTMVAReader->AddVariable( "HadEnergyOverPt",      &fMVAVar_MuHadEnergyOverPt       );      
      tmpTMVAReader->AddVariable( "EmEnergyOverPt",       &fMVAVar_MuEmEnergyOverPt        );      
      tmpTMVAReader->AddVariable( "HadS9EnergyOverPt",    &fMVAVar_MuHadS9EnergyOverPt     );      
      tmpTMVAReader->AddVariable( "EmS9EnergyOverPt",     &fMVAVar_MuEmS9EnergyOverPt      );      
      tmpTMVAReader->AddVariable( "TrkIso03OverPt",       &fMVAVar_MuTrkIso03OverPt        );
      tmpTMVAReader->AddVariable( "EMIso03OverPt",        &fMVAVar_MuEMIso03OverPt         );
      tmpTMVAReader->AddVariable( "HadIso03OverPt",       &fMVAVar_MuHadIso03OverPt        );
      tmpTMVAReader->AddVariable( "TrkIso05OverPt",       &fMVAVar_MuTrkIso05OverPt        );
      tmpTMVAReader->AddVariable( "EMIso05OverPt",        &fMVAVar_MuEMIso05OverPt         );
      tmpTMVAReader->AddVariable( "HadIso05OverPt",       &fMVAVar_MuHadIso05OverPt        );
    }

    if (type == kIDV0) {
      tmpTMVAReader->AddVariable( "TkNchi2",              &fMVAVar_MuTkNchi2               );
      if (i!=4) tmpTMVAReader->AddVariable( "GlobalNchi2",&fMVAVar_MuGlobalNchi2           );
      if (i!=4) tmpTMVAReader->AddVariable( "NValidHits", &fMVAVar_MuNValidHits            );
      tmpTMVAReader->AddVariable( "NTrackerHits",         &fMVAVar_MuNTrackerHits          );
      tmpTMVAReader->AddVariable( "NPixelHits",           &fMVAVar_MuNPixelHits            );
      if (i!=5) tmpTMVAReader->AddVariable( "NMatches",   &fMVAVar_MuNMatches              );
      tmpTMVAReader->AddVariable( "TrkKink",              &fMVAVar_MuTrkKink               );      
      tmpTMVAReader->AddVariable( "SegmentCompatibility", &fMVAVar_MuSegmentCompatibility  );      
      tmpTMVAReader->AddVariable( "CaloCompatibility",    &fMVAVar_MuCaloCompatibility     );      
      tmpTMVAReader->AddVariable( "HadEnergy",            &fMVAVar_MuHadEnergy             );      
      tmpTMVAReader->AddVariable( "EmEnergy",             &fMVAVar_MuEmEnergy              );      
      tmpTMVAReader->AddVariable( "HadS9Energy",          &fMVAVar_MuHadS9Energy           );      
      tmpTMVAReader->AddVariable( "EmS9Energy",           &fMVAVar_MuEmS9Energy            );
    }

    if (type == kIsoRingsV0) {
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p0To0p1",         &fMVAVar_ChargedIso_DR0p0To0p1        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p1To0p2",         &fMVAVar_ChargedIso_DR0p1To0p2        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p2To0p3",         &fMVAVar_ChargedIso_DR0p2To0p3        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p3To0p4",         &fMVAVar_ChargedIso_DR0p3To0p4        );
      tmpTMVAReader->AddVariable( "ChargedIso_DR0p4To0p5",         &fMVAVar_ChargedIso_DR0p4To0p5        );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p0To0p1",           &fMVAVar_GammaIso_DR0p0To0p1          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p1To0p2",           &fMVAVar_GammaIso_DR0p1To0p2          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p2To0p3",           &fMVAVar_GammaIso_DR0p2To0p3          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p3To0p4",           &fMVAVar_GammaIso_DR0p3To0p4          );
      tmpTMVAReader->AddVariable( "GammaIso_DR0p4To0p5",           &fMVAVar_GammaIso_DR0p4To0p5          );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p0To0p1",   &fMVAVar_NeutralHadronIso_DR0p0To0p1  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p1To0p2",   &fMVAVar_NeutralHadronIso_DR0p1To0p2  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p2To0p3",   &fMVAVar_NeutralHadronIso_DR0p2To0p3  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p3To0p4",   &fMVAVar_NeutralHadronIso_DR0p3To0p4  );
      tmpTMVAReader->AddVariable( "NeutralHadronIso_DR0p4To0p5",   &fMVAVar_NeutralHadronIso_DR0p4To0p5  );
    }

    if (type == kIsoDeltaR) {
      tmpTMVAReader->AddVariable("PFCharged",                     &fMVAVar_MuRelIsoPFCharged       );
      tmpTMVAReader->AddVariable("PFNeutral",                     &fMVAVar_MuRelIsoPFNeutral       );
      tmpTMVAReader->AddVariable("PFPhotons",                     &fMVAVar_MuRelIsoPFPhotons       );
      tmpTMVAReader->AddVariable("SumDeltaR",                     &fMVAVar_MuDeltaRSum             );
      tmpTMVAReader->AddVariable("DeltaRMean",                    &fMVAVar_MuDeltaRMean            );
      tmpTMVAReader->AddVariable("Density",                       &fMVAVar_MuDensity               );
    }
    
    tmpTMVAReader->BookMVA(fMethodname , weightsfiles[i] );
    std::cout << "MVABin " << i << " : MethodName = " << fMethodname 
              << " , type == " << type << " , "
              << "Load weights file : " << weightsfiles[i] 
              << std::endl;
    fTMVAReader.push_back(tmpTMVAReader);

  }
  std::cout << "Muon ID MVA Completed\n";
}

//--------------------------------------------------------------------------------------------------
UInt_t MuonIDMVA::GetMVABin( double eta, double pt, 
                             Bool_t isGlobal, Bool_t isTrackerMuon) const {
  
  //Default is to return the first bin
  uint bin = 0;

  //return the first bin if not using binned version
  if (!fUseBinnedVersion) return 0;

  if (fMVAType == MuonIDMVA::kV2 
      || fMVAType == MuonIDMVA::kV3
      || fMVAType == MuonIDMVA::kV8
      || fMVAType == MuonIDMVA::kIDIsoCombinedDetIso) {
    if (pt < 14.5 && fabs(eta) < 1.5) bin = 0;
    if (pt < 14.5 && fabs(eta) >= 1.5) bin = 1;
    if (pt >= 14.5 && pt < 20 && fabs(eta) < 1.5) bin = 2;
    if (pt >= 14.5 && pt < 20 && fabs(eta) >= 1.5) bin = 3;
    if (pt >= 20 && fabs(eta) < 1.5) bin = 4;
    if (pt >= 20 && fabs(eta) >= 1.5) bin = 5;
  }

  if (fMVAType == MuonIDMVA::kIsoRingsV0 || fMVAType == MuonIDMVA::kIDV0
      || fMVAType == MuonIDMVA::kIDIsoCombinedIsoRingsV0) {
    if (isGlobal && isTrackerMuon) {
      if (pt < 10 && fabs(eta) < 1.479) bin = 0;
      if (pt >= 10 && fabs(eta) < 1.479) bin = 1;
      if (pt < 10 && fabs(eta) >= 1.479) bin = 2;
      if (pt >= 10 && fabs(eta) >= 1.479) bin = 3;
    } else if (!isGlobal && isTrackerMuon) {
      bin = 4;
    } else if (isGlobal && !isTrackerMuon) {
      bin = 5;
    } else {
      std::cout << "Warning: Muon is not a tracker muon. Such muons are not supported. \n";
      bin = 0;
    }
  }

  if (fMVAType == MuonIDMVA::kIsoDeltaR){
    if (pt <  20 && fabs(eta) <  1.479) bin = 0;
    if (pt <  20 && fabs(eta) >= 1.479) bin = 1;
    if (pt >= 20 && fabs(eta) <  1.479) bin = 2;
    if (pt >= 20 && fabs(eta) >= 1.479) bin = 3;
  }

  return bin;
}

//--------------------------------------------------------------------------------------------------
Double_t MuonIDMVA::MVAValue(Double_t MuPt , Double_t MuEta,
                             Double_t                   MuTkNchi2, 
                             Double_t                   MuGlobalNchi2, 
                             Double_t                   MuNValidHits, 
                             Double_t                   MuNTrackerHits, 
                             Double_t                   MuNPixelHits, 
                             Double_t                   MuNMatches, 
                             Double_t                   MuD0, 
                             Double_t                   MuIP3d, 
                             Double_t                   MuIP3dSig, 
                             Double_t                   MuTrkKink, 
                             Double_t                   MuSegmentCompatibility, 
                             Double_t                   MuCaloCompatibility, 
                             Double_t                   MuHadEnergyOverPt, 
                             Double_t                   MuHoEnergyOverPt, 
                             Double_t                   MuEmEnergyOverPt, 
                             Double_t                   MuHadS9EnergyOverPt, 
                             Double_t                   MuHoS9EnergyOverPt, 
                             Double_t                   MuEmS9EnergyOverPt,
                             Double_t                   MuChargedIso03OverPt,
                             Double_t                   MuNeutralIso03OverPt,
                             Double_t                   MuChargedIso04OverPt,
                             Double_t                   MuNeutralIso04OverPt, 
                             Bool_t                     printDebug                            
                             ) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: MuonIDMVA not properly initialized.\n"; 
    return -9999;
  }

  
  //set all input variables
  fMVAVar_MuTkNchi2              = MuTkNchi2; 
  fMVAVar_MuGlobalNchi2          = MuGlobalNchi2; 
  fMVAVar_MuNValidHits           = MuNValidHits; 
  fMVAVar_MuNTrackerHits         = MuNTrackerHits; 
  fMVAVar_MuNPixelHits           = MuNPixelHits;  
  fMVAVar_MuNMatches             = MuNMatches; 
  fMVAVar_MuD0                   = MuD0; 
  fMVAVar_MuIP3d                 = MuIP3d; 
  fMVAVar_MuIP3dSig              = MuIP3dSig; 
  fMVAVar_MuTrkKink              = MuTrkKink; 
  fMVAVar_MuSegmentCompatibility = MuSegmentCompatibility; 
  fMVAVar_MuCaloCompatibility    = MuCaloCompatibility; 
  fMVAVar_MuHadEnergyOverPt      = MuHadEnergyOverPt; 
  fMVAVar_MuHoEnergyOverPt       = MuHoEnergyOverPt; 
  fMVAVar_MuEmEnergyOverPt       = MuEmEnergyOverPt; 
  fMVAVar_MuHadS9EnergyOverPt    = MuHadS9EnergyOverPt; 
  fMVAVar_MuHoS9EnergyOverPt     = MuHoS9EnergyOverPt; 
  fMVAVar_MuEmS9EnergyOverPt     = MuEmS9EnergyOverPt; 
  fMVAVar_MuChargedIso03OverPt   = MuChargedIso03OverPt; 
  fMVAVar_MuNeutralIso03OverPt   = MuNeutralIso03OverPt; 
  fMVAVar_MuChargedIso04OverPt   = MuChargedIso04OverPt; 
  fMVAVar_MuNeutralIso04OverPt   = MuNeutralIso04OverPt; 

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  reader = fTMVAReader[GetMVABin(MuEta, MuPt, kTRUE, kTRUE )];
                                                
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug) {
    std::cout << "Debug Muon MVA: "
              << MuPt << " " << MuEta << " --> MVABin " << GetMVABin(MuEta, MuPt, kTRUE, kTRUE ) << " : "     
              << fMVAVar_MuTkNchi2              << " " 
              << fMVAVar_MuGlobalNchi2          << " " 
              << fMVAVar_MuNValidHits           << " " 
              << fMVAVar_MuNTrackerHits         << " " 
              << fMVAVar_MuNPixelHits           << " "  
              << fMVAVar_MuNMatches             << " " 
              << fMVAVar_MuD0                   << " " 
              << fMVAVar_MuIP3d                 << " " 
              << fMVAVar_MuIP3dSig              << " " 
              << fMVAVar_MuTrkKink              << " " 
              << fMVAVar_MuSegmentCompatibility << " " 
              << fMVAVar_MuCaloCompatibility    << " " 
              << fMVAVar_MuHadEnergyOverPt      << " " 
              << fMVAVar_MuHoEnergyOverPt       << " " 
              << fMVAVar_MuEmEnergyOverPt       << " " 
              << fMVAVar_MuHadS9EnergyOverPt    << " " 
              << fMVAVar_MuHoS9EnergyOverPt     << " " 
              << fMVAVar_MuEmS9EnergyOverPt     << " " 
              << fMVAVar_MuChargedIso03OverPt   << " " 
              << fMVAVar_MuNeutralIso03OverPt   << " " 
              << fMVAVar_MuChargedIso04OverPt   << " " 
              << fMVAVar_MuNeutralIso04OverPt   << " " 
              << " === : === "
              << mva 
              << std::endl;
  }

  return mva;
}


//--------------------------------------------------------------------------------------------------
Double_t MuonIDMVA::MVAValue(Double_t MuPt , Double_t MuEta,
                             Double_t                   MuTkNchi2, 
                             Double_t                   MuGlobalNchi2, 
                             Double_t                   MuNValidHits, 
                             Double_t                   MuNTrackerHits, 
                             Double_t                   MuNPixelHits, 
                             Double_t                   MuNMatches, 
                             Double_t                   MuD0, 
                             Double_t                   MuIP3d, 
                             Double_t                   MuIP3dSig, 
                             Double_t                   MuTrkKink, 
                             Double_t                   MuSegmentCompatibility, 
                             Double_t                   MuCaloCompatibility, 
                             Double_t                   MuHadEnergyOverPt, 
                             Double_t                   MuHoEnergyOverPt, 
                             Double_t                   MuEmEnergyOverPt, 
                             Double_t                   MuHadS9EnergyOverPt, 
                             Double_t                   MuHoS9EnergyOverPt, 
                             Double_t                   MuEmS9EnergyOverPt,
                             Double_t                   MuTrkIso03OverPt,
                             Double_t                   MuEMIso03OverPt,
                             Double_t                   MuHadIso03OverPt,
                             Double_t                   MuTrkIso05OverPt,
                             Double_t                   MuEMIso05OverPt,
                             Double_t                   MuHadIso05OverPt,
                             Bool_t                     printDebug                            
                             ) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: MuonIDMVA not properly initialized.\n"; 
    return -9999;
  }
  
  //set all input variables
  fMVAVar_MuTkNchi2              = MuTkNchi2; 
  fMVAVar_MuGlobalNchi2          = MuGlobalNchi2; 
  fMVAVar_MuNValidHits           = MuNValidHits; 
  fMVAVar_MuNTrackerHits         = MuNTrackerHits; 
  fMVAVar_MuNPixelHits           = MuNPixelHits;  
  fMVAVar_MuNMatches             = MuNMatches; 
  fMVAVar_MuD0                   = MuD0; 
  fMVAVar_MuIP3d                 = MuIP3d; 
  fMVAVar_MuIP3dSig              = MuIP3dSig; 
  fMVAVar_MuTrkKink              = MuTrkKink; 
  fMVAVar_MuSegmentCompatibility = MuSegmentCompatibility; 
  fMVAVar_MuCaloCompatibility    = MuCaloCompatibility; 
  fMVAVar_MuHadEnergyOverPt      = MuHadEnergyOverPt; 
  fMVAVar_MuHoEnergyOverPt       = MuHoEnergyOverPt; 
  fMVAVar_MuEmEnergyOverPt       = MuEmEnergyOverPt; 
  fMVAVar_MuHadS9EnergyOverPt    = MuHadS9EnergyOverPt; 
  fMVAVar_MuHoS9EnergyOverPt     = MuHoS9EnergyOverPt; 
  fMVAVar_MuEmS9EnergyOverPt     = MuEmS9EnergyOverPt; 
  fMVAVar_MuTrkIso03OverPt       = MuTrkIso03OverPt; 
  fMVAVar_MuEMIso03OverPt        = MuEMIso03OverPt; 
  fMVAVar_MuHadIso03OverPt       = MuHadIso03OverPt; 
  fMVAVar_MuTrkIso05OverPt       = MuTrkIso05OverPt; 
  fMVAVar_MuEMIso05OverPt        = MuEMIso05OverPt; 
  fMVAVar_MuHadIso05OverPt       = MuHadIso05OverPt; 

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  reader = fTMVAReader[GetMVABin(MuEta, MuPt, kTRUE, kTRUE )];
                                                
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug) {
    std::cout << "Debug Muon MVA: "
              << MuPt << " " << MuEta << " --> MVABin " << GetMVABin(MuEta, MuPt, kTRUE, kTRUE ) << " : "     
              << fMVAVar_MuTkNchi2              << " " 
              << fMVAVar_MuGlobalNchi2          << " " 
              << fMVAVar_MuNValidHits           << " " 
              << fMVAVar_MuNTrackerHits         << " " 
              << fMVAVar_MuNPixelHits           << " "  
              << fMVAVar_MuNMatches             << " " 
              << fMVAVar_MuD0                   << " " 
              << fMVAVar_MuIP3d                 << " " 
              << fMVAVar_MuIP3dSig              << " " 
              << fMVAVar_MuTrkKink              << " " 
              << fMVAVar_MuSegmentCompatibility << " " 
              << fMVAVar_MuCaloCompatibility    << " " 
              << fMVAVar_MuHadEnergyOverPt      << " " 
              << fMVAVar_MuHoEnergyOverPt       << " " 
              << fMVAVar_MuEmEnergyOverPt       << " " 
              << fMVAVar_MuHadS9EnergyOverPt    << " " 
              << fMVAVar_MuHoS9EnergyOverPt     << " " 
              << fMVAVar_MuEmS9EnergyOverPt     << " " 
              << fMVAVar_MuTrkIso03OverPt   << " " 
              << fMVAVar_MuEMIso03OverPt   << " " 
              << fMVAVar_MuHadIso03OverPt   << " " 
              << fMVAVar_MuTrkIso05OverPt   << " " 
              << fMVAVar_MuEMIso05OverPt   << " " 
              << fMVAVar_MuHadIso05OverPt   << " " 
              << " === : === "
              << mva 
              << std::endl;
  }

  return mva;
}


Double_t MuonIDMVA::MVAValue_IsoRings( Double_t MuPt,
                                       Double_t MuEta,
                                       Bool_t MuIsGlobal,
                                       Bool_t MuIsTracker,
                                       Double_t ChargedIso_DR0p0To0p1,
                                       Double_t ChargedIso_DR0p1To0p2,
                                       Double_t ChargedIso_DR0p2To0p3,
                                       Double_t ChargedIso_DR0p3To0p4,
                                       Double_t ChargedIso_DR0p4To0p5,
                                       Double_t GammaIso_DR0p0To0p1,
                                       Double_t GammaIso_DR0p1To0p2,
                                       Double_t GammaIso_DR0p2To0p3,
                                       Double_t GammaIso_DR0p3To0p4,
                                       Double_t GammaIso_DR0p4To0p5,
                                       Double_t NeutralHadronIso_DR0p0To0p1,
                                       Double_t NeutralHadronIso_DR0p1To0p2,
                                       Double_t NeutralHadronIso_DR0p2To0p3,
                                       Double_t NeutralHadronIso_DR0p3To0p4,
                                       Double_t NeutralHadronIso_DR0p4To0p5,
                                       Bool_t printDebug) {

  if (fMVAType != MuonIDMVA::kIsoRingsV0) {
    std::cout << "Error: This function is only supported for MVAType == kIsoRingsV0.\n" << std::endl;
    assert(kFALSE);
  }

  fMVAVar_MuPt = MuPt;
  fMVAVar_MuEta = MuEta;
  fMVAVar_ChargedIso_DR0p0To0p1 = ChargedIso_DR0p0To0p1;
  fMVAVar_ChargedIso_DR0p1To0p2 = ChargedIso_DR0p1To0p2;
  fMVAVar_ChargedIso_DR0p2To0p3 = ChargedIso_DR0p2To0p3;
  fMVAVar_ChargedIso_DR0p3To0p4 = ChargedIso_DR0p3To0p4;
  fMVAVar_ChargedIso_DR0p4To0p5 = ChargedIso_DR0p4To0p5;
  fMVAVar_GammaIso_DR0p0To0p1 = GammaIso_DR0p0To0p1;
  fMVAVar_GammaIso_DR0p1To0p2 = GammaIso_DR0p1To0p2;
  fMVAVar_GammaIso_DR0p2To0p3 = GammaIso_DR0p2To0p3;
  fMVAVar_GammaIso_DR0p3To0p4 = GammaIso_DR0p3To0p4;
  fMVAVar_GammaIso_DR0p4To0p5 = GammaIso_DR0p4To0p5;
  fMVAVar_NeutralHadronIso_DR0p0To0p1 = NeutralHadronIso_DR0p0To0p1;
  fMVAVar_NeutralHadronIso_DR0p1To0p2 = NeutralHadronIso_DR0p1To0p2;
  fMVAVar_NeutralHadronIso_DR0p2To0p3 = NeutralHadronIso_DR0p2To0p3;
  fMVAVar_NeutralHadronIso_DR0p3To0p4 = NeutralHadronIso_DR0p3To0p4;
  fMVAVar_NeutralHadronIso_DR0p4To0p5 = NeutralHadronIso_DR0p4To0p5;

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;

  if (printDebug == kTRUE) {
    std::cout <<" -> BIN: " << fMVAVar_MuEta << " " << fMVAVar_MuPt << " : " 
              << GetMVABin( fMVAVar_MuEta , fMVAVar_MuPt, MuIsGlobal, MuIsTracker) << std::endl;
  }
  reader = fTMVAReader[GetMVABin( fMVAVar_MuEta , fMVAVar_MuPt, MuIsGlobal, MuIsTracker)];                                              
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug == kTRUE) {

    std::cout << "Debug Muon MVA: \n";
    std::cout << fMVAVar_ChargedIso_DR0p0To0p1 << " "
              << fMVAVar_ChargedIso_DR0p1To0p2 << " "
              << fMVAVar_ChargedIso_DR0p2To0p3 << " "
              << fMVAVar_ChargedIso_DR0p3To0p4 << " "
              << fMVAVar_ChargedIso_DR0p4To0p5 << " "
              << fMVAVar_GammaIso_DR0p0To0p1 << " "
              << fMVAVar_GammaIso_DR0p1To0p2 << " "
              << fMVAVar_GammaIso_DR0p2To0p3 << " "
              << fMVAVar_GammaIso_DR0p3To0p4 << " "
              << fMVAVar_GammaIso_DR0p4To0p5 << " "
              << fMVAVar_NeutralHadronIso_DR0p0To0p1 << " "
              << fMVAVar_NeutralHadronIso_DR0p1To0p2 << " "
              << fMVAVar_NeutralHadronIso_DR0p2To0p3 << " "
              << fMVAVar_NeutralHadronIso_DR0p3To0p4 << " "
              << fMVAVar_NeutralHadronIso_DR0p4To0p5 << " "  
              << std::endl;
    std::cout << "MVA: " << mva << " "    
              << std::endl;    
  }  
  return mva;
}



Double_t MuonIDMVA::MVAValue_ID( Double_t MuPt,
                                 Double_t MuEta,
                                 Bool_t MuIsGlobal,
                                 Bool_t MuIsTracker,                                
                                 Double_t MuTkNchi2, 
                                 Double_t MuGlobalNchi2, 
                                 Double_t MuNValidHits, 
                                 Double_t MuNTrackerHits, 
                                 Double_t MuNPixelHits, 
                                 Double_t MuNMatches, 
                                 Double_t MuTrkKink, 
                                 Double_t MuSegmentCompatibility, 
                                 Double_t MuCaloCompatibility, 
                                 Double_t MuHadEnergy, 
                                 Double_t MuEmEnergy, 
                                 Double_t MuHadS9Energy, 
                                 Double_t MuEmS9Energy, 
                                 Bool_t printDebug) {

  if (fMVAType != MuonIDMVA::kIDV0) {
    std::cout << "Error: This function is only supported for MVAType == kIDV0.\n" << std::endl;
    assert(kFALSE);
  }

  fMVAVar_MuPt = MuPt;
  fMVAVar_MuEta = MuEta;

  fMVAVar_MuTkNchi2 = MuTkNchi2; 
  fMVAVar_MuGlobalNchi2 = MuGlobalNchi2; 
  fMVAVar_MuNValidHits = MuNValidHits; 
  fMVAVar_MuNTrackerHits = MuNTrackerHits; 
  fMVAVar_MuNPixelHits = MuNPixelHits; 
  fMVAVar_MuNMatches = MuNMatches; 
  fMVAVar_MuTrkKink = MuTrkKink; 
  fMVAVar_MuSegmentCompatibility = MuSegmentCompatibility; 
  fMVAVar_MuCaloCompatibility = MuCaloCompatibility; 
  fMVAVar_MuHadEnergy = MuHadEnergy; 
  fMVAVar_MuEmEnergy = MuEmEnergy; 
  fMVAVar_MuHadS9Energy = MuHadS9Energy; 
  fMVAVar_MuEmS9Energy = MuEmS9Energy; 

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;

  if (printDebug == kTRUE) {
    std::cout <<" -> BIN: " << fMVAVar_MuEta << " " << fMVAVar_MuPt << " : " 
              << GetMVABin( fMVAVar_MuEta , fMVAVar_MuPt, MuIsGlobal, MuIsTracker) << std::endl;
  }
  reader = fTMVAReader[GetMVABin( fMVAVar_MuEta , fMVAVar_MuPt, MuIsGlobal, MuIsTracker)];                                              
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug == kTRUE) {

    std::cout << "Debug Muon MVA: \n";
    std::cout << fMVAVar_MuTkNchi2              << " " 
              << fMVAVar_MuGlobalNchi2          << " " 
              << fMVAVar_MuNValidHits           << " " 
              << fMVAVar_MuNTrackerHits         << " " 
              << fMVAVar_MuNPixelHits           << " "  
              << fMVAVar_MuNMatches             << " "       
              << fMVAVar_MuTrkKink              << " " 
              << fMVAVar_MuSegmentCompatibility << " " 
              << fMVAVar_MuCaloCompatibility    << " " 
              << fMVAVar_MuHadEnergy            << " "  
              << fMVAVar_MuEmEnergy             << " " 
              << fMVAVar_MuHadS9Energy          << " "  
              << fMVAVar_MuEmS9Energy           << " "    
              << std::endl;
    std::cout << "MVA: " << mva << " "    
              << std::endl;    
  }
  return mva;
}



//--------------------------------------------------------------------------------------------------
Double_t MuonIDMVA::MVAValue(const Muon *mu, const Vertex *vertex, MuonTools *fMuonTools,
                             const PFCandidateCol *PFCands, 
                             const PileupEnergyDensityCol *PileupEnergyDensity, 
                             Bool_t printDebug) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: MuonIDMVA not properly initialized.\n"; 
    return -9999;
  }

  const Track *muTrk=0;
  if(mu->HasTrackerTrk())         { muTrk = mu->TrackerTrk();    }
  else if(mu->HasStandaloneTrk()) { muTrk = mu->StandaloneTrk(); } 
  
  Double_t muNchi2 = 0.0; 
  if(mu->HasGlobalTrk())          { muNchi2 = mu->GlobalTrk()->RChi2();     }
  else if(mu->HasStandaloneTrk()) { muNchi2 = mu->StandaloneTrk()->RChi2(); }
  else if(mu->HasTrackerTrk())    { muNchi2 = mu->TrackerTrk()->RChi2();    }

  Double_t ChargedIso03 = 0;
  Double_t NeutralIso03_05Threshold = 0;
  Double_t ChargedIso04 = 0;
  Double_t NeutralIso04_05Threshold = 0;
  ChargedIso03 = IsolationTools::PFMuonIsolation(mu, PFCands, vertex, 0.1, 99999, 0.3, 0.0, 0.0);
  NeutralIso03_05Threshold = IsolationTools::PFMuonIsolation(mu, PFCands, vertex, 0.0, 0.5, 0.3, 0.0, 0.0);
  ChargedIso04 = IsolationTools::PFMuonIsolation(mu, PFCands, vertex, 0.1, 99999, 0.4, 0.0, 0.0);
  NeutralIso04_05Threshold = IsolationTools::PFMuonIsolation(mu, PFCands, vertex, 0.0, 0.5, 0.4, 0.0, 0.0);
  
  Double_t Rho = PileupEnergyDensity->At(0)->Rho(fRhoAlgo);
 
  //set all input variables
  fMVAVar_MuTkNchi2              = muTrk->RChi2();
  fMVAVar_MuGlobalNchi2          = muNchi2;
  fMVAVar_MuNValidHits           = mu->NValidHits();
  fMVAVar_MuNTrackerHits         = muTrk->NHits();
  fMVAVar_MuNPixelHits           = muTrk->NPixelHits();
  fMVAVar_MuNMatches             = mu->NMatches();
  fMVAVar_MuD0                   = muTrk->D0Corrected(*vertex);
  fMVAVar_MuIP3d                 = mu->Ip3dPV();
  fMVAVar_MuIP3dSig              = mu->Ip3dPVSignificance();
  fMVAVar_MuTrkKink              = mu->TrkKink();
  fMVAVar_MuSegmentCompatibility = fMuonTools->GetSegmentCompatibility(mu);
  fMVAVar_MuCaloCompatibility    = fMuonTools->GetCaloCompatibility(mu, kTRUE, kTRUE);
  fMVAVar_MuHadEnergyOverPt      = (mu->HadEnergy() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuHadEnergy,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuHoEnergyOverPt       = (mu->HoEnergy() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuHoEnergy,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuEmEnergyOverPt       = (mu->EmEnergy() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuEmEnergy,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuHadS9EnergyOverPt    = (mu->HadS9Energy() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuHadS9Energy,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuHoS9EnergyOverPt     = (mu->HoS9Energy() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuHoS9Energy,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuEmS9EnergyOverPt     = (mu->EmS9Energy() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuEmS9Energy,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuChargedIso03OverPt   = (ChargedIso03 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuChargedIso03,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuNeutralIso03OverPt   = (NeutralIso03_05Threshold - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralIso03,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuChargedIso04OverPt   = (ChargedIso04 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuChargedIso04,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuNeutralIso04OverPt   = (NeutralIso04_05Threshold - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralIso04,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuTrkIso03OverPt       = (mu->IsoR03SumPt() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuTrkIso03,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuEMIso03OverPt        = (mu->IsoR03EmEt() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuEMIso03,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuHadIso03OverPt       = (mu->IsoR03HadEt() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuHadIso03,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuTrkIso05OverPt       = (mu->IsoR05SumPt() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuTrkIso05,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuEMIso05OverPt        = (mu->IsoR05EmEt() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuEMIso05,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuHadIso05OverPt       = (mu->IsoR05HadEt() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuHadIso05,muTrk->Eta()))/muTrk->Pt();

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  reader = fTMVAReader[GetMVABin(muTrk->Eta(), muTrk->Pt(), mu->IsGlobalMuon(), mu->IsTrackerMuon())];
                                                
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug) {
    std::cout << "Debug Muon MVA: "
              << mu->Pt() << " " << mu->Eta() << " " << mu->Phi() << " : "
              << muTrk->Pt() << " " << muTrk->Eta() << " --> MVABin " << GetMVABin(muTrk->Eta(), muTrk->Pt(), mu->IsGlobalMuon(), mu->IsTrackerMuon()) << " : "     
              << fMVAVar_MuTkNchi2              << " " 
              << fMVAVar_MuGlobalNchi2          << " " 
              << fMVAVar_MuNValidHits           << " " 
              << fMVAVar_MuNTrackerHits         << " " 
              << fMVAVar_MuNPixelHits           << " "  
              << fMVAVar_MuNMatches             << " " 
              << fMVAVar_MuD0                   << " " 
              << fMVAVar_MuIP3d                 << " " 
              << fMVAVar_MuIP3dSig              << " " 
              << fMVAVar_MuTrkKink              << " " 
              << fMVAVar_MuSegmentCompatibility << " " 
              << fMVAVar_MuCaloCompatibility    << " " 
              << fMVAVar_MuHadEnergyOverPt      << " " 
              << fMVAVar_MuHoEnergyOverPt       << " " 
              << fMVAVar_MuEmEnergyOverPt       << " " 
              << fMVAVar_MuHadS9EnergyOverPt    << " " 
              << fMVAVar_MuHoS9EnergyOverPt     << " " 
              << fMVAVar_MuEmS9EnergyOverPt     << " " 
              << fMVAVar_MuChargedIso03OverPt   << " " 
              << fMVAVar_MuNeutralIso03OverPt   << " " 
              << fMVAVar_MuChargedIso04OverPt   << " " 
              << fMVAVar_MuNeutralIso04OverPt   << " " 
              << fMVAVar_MuTrkIso03OverPt   << " " 
              << fMVAVar_MuEMIso03OverPt   << " " 
              << fMVAVar_MuHadIso03OverPt   << " " 
              << fMVAVar_MuTrkIso05OverPt   << " " 
              << fMVAVar_MuEMIso05OverPt   << " " 
              << fMVAVar_MuHadIso05OverPt   << " " 
              << " === : === "
              << mva 
              << std::endl;
  }

  return mva;
}


//--------------------------------------------------------------------------------------------------
Double_t MuonIDMVA::MVAValue(const Muon *mu, const Vertex *vertex, MuonTools *fMuonTools,
                             const PFCandidateCol *PFCands, 
                             const PileupEnergyDensityCol *PileupEnergyDensity, 
                             MuonTools::EMuonEffectiveAreaTarget EffectiveAreaTarget,
                             const ElectronCol *goodElectrons,
                             const MuonCol *goodMuons,
                             Bool_t printDebug) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: MuonIDMVA not properly initialized.\n"; 
    return -9999;
  }

  const Track *muTrk=0;
  if(mu->HasTrackerTrk())         { muTrk = mu->TrackerTrk();    }
  else if(mu->HasStandaloneTrk()) { muTrk = mu->StandaloneTrk(); } 
  
  Double_t muNchi2 = 0.0; 
  if(mu->HasGlobalTrk())          { muNchi2 = mu->GlobalTrk()->RChi2();     }
  else if(mu->HasStandaloneTrk()) { muNchi2 = mu->StandaloneTrk()->RChi2(); }
  else if(mu->HasTrackerTrk())    { muNchi2 = mu->TrackerTrk()->RChi2();    }

  Double_t ChargedIso03 = 0;
  Double_t NeutralIso03_05Threshold = 0;
  Double_t ChargedIso04 = 0;
  Double_t NeutralIso04_05Threshold = 0;
  ChargedIso03 = IsolationTools::PFMuonIsolation(mu, PFCands, vertex, 0.1, 99999, 0.3, 0.0, 0.0);
  NeutralIso03_05Threshold = IsolationTools::PFMuonIsolation(mu, PFCands, vertex, 0.0, 0.5, 0.3, 0.0, 0.0);
  ChargedIso04 = IsolationTools::PFMuonIsolation(mu, PFCands, vertex, 0.1, 99999, 0.4, 0.0, 0.0);
  NeutralIso04_05Threshold = IsolationTools::PFMuonIsolation(mu, PFCands, vertex, 0.0, 0.5, 0.4, 0.0, 0.0);
  
  Double_t Rho = PileupEnergyDensity->At(0)->Rho(fRhoAlgo);

  //set all input variables
  fMVAVar_MuPt                   = muTrk->Pt();
  fMVAVar_MuEta                  = muTrk->Eta();
  fMVAVar_MuTkNchi2              = muTrk->RChi2();
  fMVAVar_MuGlobalNchi2          = muNchi2;
  fMVAVar_MuNValidHits           = mu->NValidHits();
  fMVAVar_MuNTrackerHits         = muTrk->NHits();
  fMVAVar_MuNPixelHits           = muTrk->NPixelHits();
  fMVAVar_MuNMatches             = mu->NMatches();
  fMVAVar_MuD0                   = muTrk->D0Corrected(*vertex);
  fMVAVar_MuIP3d                 = mu->Ip3dPV();
  fMVAVar_MuIP3dSig              = mu->Ip3dPVSignificance();
  fMVAVar_MuTrkKink              = mu->TrkKink();
  fMVAVar_MuSegmentCompatibility = fMuonTools->GetSegmentCompatibility(mu);
  fMVAVar_MuCaloCompatibility    = fMuonTools->GetCaloCompatibility(mu, kTRUE, kTRUE);
  fMVAVar_MuHadEnergy            = mu->HadEnergy() ;
  fMVAVar_MuEmEnergy             = mu->EmEnergy();
  fMVAVar_MuHadS9Energy          = mu->HadS9Energy();
  fMVAVar_MuEmS9Energy           = mu->EmS9Energy();
  fMVAVar_MuChargedIso03OverPt   = (ChargedIso03 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuChargedIso03,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuNeutralIso03OverPt   = (NeutralIso03_05Threshold - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralIso03,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuChargedIso04OverPt   = (ChargedIso04 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuChargedIso04,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuNeutralIso04OverPt   = (NeutralIso04_05Threshold - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralIso04,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuTrkIso03OverPt       = (mu->IsoR03SumPt() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuTrkIso03,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuEMIso03OverPt        = (mu->IsoR03EmEt() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuEMIso03,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuHadIso03OverPt       = (mu->IsoR03HadEt() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuHadIso03,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuTrkIso05OverPt       = (mu->IsoR05SumPt() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuTrkIso05,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuEMIso05OverPt        = (mu->IsoR05EmEt() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuEMIso05,muTrk->Eta()))/muTrk->Pt();
  fMVAVar_MuHadIso05OverPt       = (mu->IsoR05HadEt() - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuHadIso05,muTrk->Eta()))/muTrk->Pt();


  Double_t tmpChargedIso_DR0p0To0p1  = 0;
  Double_t tmpChargedIso_DR0p1To0p2  = 0;
  Double_t tmpChargedIso_DR0p2To0p3  = 0;
  Double_t tmpChargedIso_DR0p3To0p4  = 0;
  Double_t tmpChargedIso_DR0p4To0p5  = 0;
  Double_t tmpGammaIso_DR0p0To0p1  = 0;
  Double_t tmpGammaIso_DR0p1To0p2  = 0;
  Double_t tmpGammaIso_DR0p2To0p3  = 0;
  Double_t tmpGammaIso_DR0p3To0p4  = 0;
  Double_t tmpGammaIso_DR0p4To0p5  = 0;
  Double_t tmpNeutralHadronIso_DR0p0To0p1  = 0;
  Double_t tmpNeutralHadronIso_DR0p1To0p2  = 0;
  Double_t tmpNeutralHadronIso_DR0p2To0p3  = 0;
  Double_t tmpNeutralHadronIso_DR0p3To0p4  = 0;
  Double_t tmpNeutralHadronIso_DR0p4To0p5  = 0;

  Double_t tmpMuDeltaRMean = 0;
  Double_t tmpMuDeltaRSum = 0;
  Double_t tmpMuDensity = 0;
  Double_t tmpMuNPFCand = 0;

  for (UInt_t p=0; p<PFCands->GetEntries();p++) {   
    const PFCandidate *pf = PFCands->At(p);
      
    //exclude the muon itself
    if(pf->TrackerTrk() && mu->TrackerTrk() &&
       pf->TrackerTrk() == mu->TrackerTrk()) continue;      

    //************************************************************
    // New Isolation Calculations
    //************************************************************
    Double_t dr = MathUtils::DeltaR(mu->Mom(), pf->Mom());

    if (dr < 0.5) {
      Bool_t IsLeptonFootprint = kFALSE;
      //************************************************************
      // Lepton Footprint Removal
      //************************************************************            
      if(goodElectrons) {
        for (UInt_t q=0; q < goodElectrons->GetEntries() ; ++q) {
	  //if pf candidate matches an electron passing ID cuts, then veto it
	  if(pf->GsfTrk() && goodElectrons->At(q)->GsfTrk() &&
	     pf->GsfTrk() == goodElectrons->At(q)->GsfTrk()) IsLeptonFootprint = kTRUE;
	  if(pf->TrackerTrk() && goodElectrons->At(q)->TrackerTrk() &&
	     pf->TrackerTrk() == goodElectrons->At(q)->TrackerTrk()) IsLeptonFootprint = kTRUE;
	  //if pf candidate lies in veto regions of electron passing ID cuts, then veto it
	  if(pf->BestTrk() && fabs(goodElectrons->At(q)->SCluster()->Eta()) >= 1.479 
             && MathUtils::DeltaR(goodElectrons->At(q)->Mom(), pf->Mom()) < 0.015) IsLeptonFootprint = kTRUE;
	  if(pf->PFType() == PFCandidate::eGamma && fabs(goodElectrons->At(q)->SCluster()->Eta()) >= 1.479 &&
	     MathUtils::DeltaR(goodElectrons->At(q)->Mom(), pf->Mom()) < 0.08) IsLeptonFootprint = kTRUE;
        }
      }
      if(goodMuons) {
        for (UInt_t q=0; q < goodMuons->GetEntries() ; ++q) {
	  //if pf candidate matches an muon passing ID cuts, then veto it
	  if(pf->TrackerTrk() && goodMuons->At(q)->TrackerTrk() &&
	     pf->TrackerTrk() == goodMuons->At(q)->TrackerTrk()) IsLeptonFootprint = kTRUE;
	  //if pf candidate lies in veto regions of muon passing ID cuts, then veto it
	  if(pf->BestTrk() && MathUtils::DeltaR(goodMuons->At(q)->Mom(), pf->Mom()) < 0.01) IsLeptonFootprint = kTRUE;
        }
      }

      if (!IsLeptonFootprint) {
	Bool_t passVeto = kTRUE;
	//Charged
        if(pf->BestTrk()) {	  	   
          if (!(fabs(pf->BestTrk()->DzCorrected(*vertex) - mu->BestTrk()->DzCorrected(*vertex)) < 0.2)) passVeto = kFALSE;
          //************************************************************
          // Veto any PFmuon, or PFEle
          if (pf->PFType() == PFCandidate::eElectron || pf->PFType() == PFCandidate::eMuon) passVeto = kFALSE;
          //************************************************************
          //************************************************************
          // Footprint Veto
          if (dr < 0.01) passVeto = kFALSE;
          //************************************************************
          if (passVeto) {
            if (dr < 0.1) tmpChargedIso_DR0p0To0p1 += pf->Pt();
            if (dr >= 0.1 && dr < 0.2) tmpChargedIso_DR0p1To0p2 += pf->Pt();
            if (dr >= 0.2 && dr < 0.3) tmpChargedIso_DR0p2To0p3 += pf->Pt();
            if (dr >= 0.3 && dr < 0.4) tmpChargedIso_DR0p3To0p4 += pf->Pt();
            if (dr >= 0.4 && dr < 0.5) tmpChargedIso_DR0p4To0p5 += pf->Pt();
          } //pass veto	   
        }
        //Gamma
        else if (pf->PFType() == PFCandidate::eGamma) {
          if (dr < 0.1) tmpGammaIso_DR0p0To0p1 += pf->Pt();
          if (dr >= 0.1 && dr < 0.2) tmpGammaIso_DR0p1To0p2 += pf->Pt();
          if (dr >= 0.2 && dr < 0.3) tmpGammaIso_DR0p2To0p3 += pf->Pt();
          if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += pf->Pt();
          if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += pf->Pt();
        }
        //NeutralHadron
        else {
          if (dr < 0.1) tmpNeutralHadronIso_DR0p0To0p1 += pf->Pt();
          if (dr >= 0.1 && dr < 0.2) tmpNeutralHadronIso_DR0p1To0p2 += pf->Pt();
          if (dr >= 0.2 && dr < 0.3) tmpNeutralHadronIso_DR0p2To0p3 += pf->Pt();
          if (dr >= 0.3 && dr < 0.4) tmpNeutralHadronIso_DR0p3To0p4 += pf->Pt();
          if (dr >= 0.4 && dr < 0.5) tmpNeutralHadronIso_DR0p4To0p5 += pf->Pt();
        }

        if (dr < 0.5) {
          tmpMuNPFCand++;
          tmpMuDeltaRMean += dr;
          tmpMuDeltaRSum  += dr;
          tmpMuDensity    += pf->Pt() / dr;
        }
      } //not lepton footprint
    } //in 1.0 dr cone
  } //loop over PF candidates
  
  //   Double_t fMVAVar_ChargedIso_DR0p0To0p1  = 0;
  //   Double_t fMVAVar_ChargedIso_DR0p1To0p2  = 0;
  //   Double_t fMVAVar_ChargedIso_DR0p2To0p3  = 0;
  //   Double_t fMVAVar_ChargedIso_DR0p3To0p4  = 0;
  //   Double_t fMVAVar_ChargedIso_DR0p4To0p5  = 0;
  //   Double_t fMVAVar_GammaIso_DR0p0To0p1  = 0;
  //   Double_t fMVAVar_GammaIso_DR0p1To0p2  = 0;
  //   Double_t fMVAVar_GammaIso_DR0p2To0p3  = 0;
  //   Double_t fMVAVar_GammaIso_DR0p3To0p4  = 0;
  //   Double_t fMVAVar_GammaIso_DR0p4To0p5  = 0;
  //   Double_t fMVAVar_NeutralHadronIso_DR0p0To0p1  = 0;
  //   Double_t fMVAVar_NeutralHadronIso_DR0p1To0p2  = 0;
  //   Double_t fMVAVar_NeutralHadronIso_DR0p2To0p3  = 0;
  //   Double_t fMVAVar_NeutralHadronIso_DR0p3To0p4  = 0;
  //   Double_t fMVAVar_NeutralHadronIso_DR0p4To0p5  = 0;

  fMVAVar_ChargedIso_DR0p0To0p1 = TMath::Min((tmpChargedIso_DR0p0To0p1)/mu->Pt(), 2.5);
  fMVAVar_ChargedIso_DR0p1To0p2 = TMath::Min((tmpChargedIso_DR0p1To0p2)/mu->Pt(), 2.5);
  fMVAVar_ChargedIso_DR0p2To0p3 = TMath::Min((tmpChargedIso_DR0p2To0p3)/mu->Pt(), 2.5);
  fMVAVar_ChargedIso_DR0p3To0p4 = TMath::Min((tmpChargedIso_DR0p3To0p4)/mu->Pt(), 2.5);
  fMVAVar_ChargedIso_DR0p4To0p5 = TMath::Min((tmpChargedIso_DR0p4To0p5)/mu->Pt(), 2.5); 
  fMVAVar_GammaIso_DR0p0To0p1 = TMath::Max(TMath::Min((tmpGammaIso_DR0p0To0p1 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuGammaIsoDR0p0To0p1, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_GammaIso_DR0p1To0p2 = TMath::Max(TMath::Min((tmpGammaIso_DR0p1To0p2 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuGammaIsoDR0p1To0p2, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_GammaIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpGammaIso_DR0p2To0p3 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuGammaIsoDR0p2To0p3, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_GammaIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpGammaIso_DR0p3To0p4 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuGammaIsoDR0p3To0p4, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_GammaIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpGammaIso_DR0p4To0p5 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuGammaIsoDR0p4To0p5, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p0To0p1 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p0To0p1 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralHadronIsoDR0p0To0p1, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p1To0p2 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p1To0p2 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralHadronIsoDR0p1To0p2, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p2To0p3 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralHadronIsoDR0p2To0p3, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p3To0p4 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralHadronIsoDR0p3To0p4, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p4To0p5 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralHadronIsoDR0p4To0p5, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);

  // Variables for dR MVA 
  fMVAVar_MuRelIsoPFCharged  = 0.;
  fMVAVar_MuRelIsoPFCharged += TMath::Min((tmpChargedIso_DR0p0To0p1)/mu->Pt(), 2.5);
  fMVAVar_MuRelIsoPFCharged += TMath::Min((tmpChargedIso_DR0p1To0p2)/mu->Pt(), 2.5);
  fMVAVar_MuRelIsoPFCharged += TMath::Min((tmpChargedIso_DR0p2To0p3)/mu->Pt(), 2.5);
  fMVAVar_MuRelIsoPFCharged += TMath::Min((tmpChargedIso_DR0p3To0p4)/mu->Pt(), 2.5);
  fMVAVar_MuRelIsoPFCharged += TMath::Min((tmpChargedIso_DR0p4To0p5)/mu->Pt(), 2.5);
  
  fMVAVar_MuRelIsoPFNeutral  = 0.;
  fMVAVar_MuRelIsoPFNeutral += TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p0To0p1 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralHadronIsoDR0p0To0p1, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_MuRelIsoPFNeutral += TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p1To0p2 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralHadronIsoDR0p1To0p2, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_MuRelIsoPFNeutral += TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p2To0p3 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralHadronIsoDR0p2To0p3, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_MuRelIsoPFNeutral += TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p3To0p4 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralHadronIsoDR0p3To0p4, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_MuRelIsoPFNeutral += TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p4To0p5 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuNeutralHadronIsoDR0p4To0p5, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);

  fMVAVar_MuRelIsoPFPhotons  = 0.;
  fMVAVar_MuRelIsoPFPhotons += TMath::Max(TMath::Min((tmpGammaIso_DR0p0To0p1 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuGammaIsoDR0p0To0p1, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_MuRelIsoPFPhotons += TMath::Max(TMath::Min((tmpGammaIso_DR0p1To0p2 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuGammaIsoDR0p1To0p2, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_MuRelIsoPFPhotons += TMath::Max(TMath::Min((tmpGammaIso_DR0p2To0p3 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuGammaIsoDR0p2To0p3, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_MuRelIsoPFPhotons += TMath::Max(TMath::Min((tmpGammaIso_DR0p3To0p4 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuGammaIsoDR0p3To0p4, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);
  fMVAVar_MuRelIsoPFPhotons += TMath::Max(TMath::Min((tmpGammaIso_DR0p4To0p5 - Rho*MuonTools::MuonEffectiveArea(MuonTools::kMuGammaIsoDR0p4To0p5, mu->Eta(), EffectiveAreaTarget))/mu->Pt(), 2.5), 0.0);

  fMVAVar_MuDeltaRMean      = tmpMuDeltaRMean/TMath::Max(1.0,tmpMuNPFCand);
  fMVAVar_MuDeltaRSum       = tmpMuDeltaRSum;
  fMVAVar_MuDensity         = tmpMuDensity;

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
 
  if (printDebug) {
    std::cout <<" -> BIN: " << fMVAVar_MuEta << " " << fMVAVar_MuPt << " : " 
              << GetMVABin(muTrk->Eta(), muTrk->Pt(), mu->IsGlobalMuon(), mu->IsTrackerMuon() )
              << std::endl;
  }

  reader = fTMVAReader[GetMVABin(muTrk->Eta(), muTrk->Pt(), mu->IsGlobalMuon(), mu->IsTrackerMuon() )];
                              
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug) {
    std::cout << "Debug Muon MVA: ";
    std::cout << " MuTkNchi2 " << fMVAVar_MuTkNchi2              
              << " MuGlobalNchi2 " << fMVAVar_MuGlobalNchi2          
              << " MuNValidHits " << fMVAVar_MuNValidHits           
              << " MuNTrackerHits " << fMVAVar_MuNTrackerHits         
              << " MuNPixelHits " << fMVAVar_MuNPixelHits           
              << " MuNMatches " << fMVAVar_MuNMatches             
              << " MuD0 " << fMVAVar_MuD0                
              << " MuIP3d " << fMVAVar_MuIP3d               
              << " MuIP3dSig " << fMVAVar_MuIP3dSig            
              << " MuTrkKink " << fMVAVar_MuTrkKink              
              << " MuSegmentCompatibility " << fMVAVar_MuSegmentCompatibility 
              << " MuCaloCompatibility " << fMVAVar_MuCaloCompatibility    
              << " MuHadEnergy " << fMVAVar_MuHadEnergy      
              << " MuEmEnergy " << fMVAVar_MuEmEnergy       
              << " MuHadS9Energy " << fMVAVar_MuHadS9Energy    
              << " MuEmS9Energy " << fMVAVar_MuEmS9Energy     
              << " eta " << fMVAVar_MuEta  
              << " pt " << fMVAVar_MuPt
              << " isoInfo: ";
    std::cout << fMVAVar_ChargedIso_DR0p0To0p1 << " "
              << fMVAVar_ChargedIso_DR0p1To0p2 << " "
              << fMVAVar_ChargedIso_DR0p2To0p3 << " "
              << fMVAVar_ChargedIso_DR0p3To0p4 << " "
              << fMVAVar_ChargedIso_DR0p4To0p5 << " "
              << fMVAVar_GammaIso_DR0p0To0p1 << " "
              << fMVAVar_GammaIso_DR0p1To0p2 << " "
              << fMVAVar_GammaIso_DR0p2To0p3 << " "
              << fMVAVar_GammaIso_DR0p3To0p4 << " "
              << fMVAVar_GammaIso_DR0p4To0p5 << " "
              << fMVAVar_NeutralHadronIso_DR0p0To0p1 << " "
              << fMVAVar_NeutralHadronIso_DR0p1To0p2 << " "
              << fMVAVar_NeutralHadronIso_DR0p2To0p3 << " "
              << fMVAVar_NeutralHadronIso_DR0p3To0p4 << " "
              << fMVAVar_NeutralHadronIso_DR0p4To0p5;
    std::cout << " MuRelIsoPFCharged: " << fMVAVar_MuRelIsoPFCharged
	      << " MuRelIsoPFNeutral: " << fMVAVar_MuRelIsoPFNeutral
	      << " MuRelIsoPFPhotons: " << fMVAVar_MuRelIsoPFPhotons
	      << " MuDeltaRMean: "      << fMVAVar_MuDeltaRMean
	      << " MuDeltaRMean: "      << fMVAVar_MuDeltaRMean
	      << " MuDensity: "         << fMVAVar_MuDensity;	      
    std::cout << " MVA: " << mva 
              << std::endl;
  }

  return mva;
}
