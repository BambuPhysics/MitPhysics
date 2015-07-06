#include "MitPhysics/Utils/interface/ElectronIDMVA.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>
#include <TRandom3.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"


ClassImp(mithep::ElectronIDMVA)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
ElectronIDMVA::ElectronIDMVA() :
fMethodname("BDTG method"),
fIsInitialized(kFALSE),
fMVAType(ElectronIDMVA::kUninitialized),
fUseBinnedVersion(kTRUE),
fNMVABins(0),
fRhoAlgo(mithep::PileupEnergyDensity::kHighEta)
{
  // Constructor.
}


//--------------------------------------------------------------------------------------------------
ElectronIDMVA::~ElectronIDMVA()
{
  for(UInt_t i=0; i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMVA::Initialize(std::string const& methodName,
                               std::string const& weightsfile,
                               ElectronIDMVA::MVAType type,
			       UInt_t algo/* = kHighEta*/)
{
  
  std::vector<std::string> tempWeightFileVector;
  tempWeightFileVector.push_back(weightsfile);
  Initialize(methodName,type,kFALSE,tempWeightFileVector,algo);
}

//--------------------------------------------------------------------------------------------------
void ElectronIDMVA::Initialize(TString const& methodName,
                               TString const& Subdet0Pt10To20Weights , 
                               TString const& Subdet1Pt10To20Weights , 
                               TString const& Subdet2Pt10To20Weights,
                               TString const& Subdet0Pt20ToInfWeights,
                               TString const& Subdet1Pt20ToInfWeights, 
                               TString const& Subdet2Pt20ToInfWeights,
                               ElectronIDMVA::MVAType type,
			       UInt_t algo/* = kHighEta*/)
{
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
void ElectronIDMVA::Initialize(std::string const& methodName,
                               ElectronIDMVA::MVAType type,
                               Bool_t useBinnedVersion,
                               std::vector<std::string> const& weightsfiles,
			       UInt_t algo/* = kHighEta*/)
{

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
  } else if (type == kBaseline 
             ||type == kNoIPInfo
             ||type == kWithIPInfo
             ||type == kIDIsoCombined) {
    ExpectedNBins = 6;
  } else if (type == kIDEGamma2012TrigV0 || 
             type == kIDEGamma2012NonTrigV0 || 
             type == kIDEGamma2012NonTrigV1 || 
             type == kIDHWW2012TrigV0 ||
             type == kIDIsoCombinedHWW2012TrigV4
             
    ) {
    ExpectedNBins = 6;
  } else if (type == kIsoRingsV0) {
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

    if (type == kBaseline) {
      tmpTMVAReader->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta            );
      tmpTMVAReader->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                   );
      tmpTMVAReader->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                   );
      tmpTMVAReader->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                    );
      tmpTMVAReader->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi            );
      tmpTMVAReader->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                    );
      tmpTMVAReader->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP    );      
    }    
    if (type == kNoIPInfo) {
      tmpTMVAReader->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta            );
      tmpTMVAReader->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                   );
      tmpTMVAReader->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                   );
      tmpTMVAReader->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                    );
      tmpTMVAReader->AddVariable( "EOverP",                &fMVAVar_EleEOverP                   );
      tmpTMVAReader->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout     );
      tmpTMVAReader->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi            );
      tmpTMVAReader->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                    );
      tmpTMVAReader->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP    );      
      tmpTMVAReader->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn      );
    }
    if (type == kWithIPInfo) {
      tmpTMVAReader->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta            );
      tmpTMVAReader->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                   );
      tmpTMVAReader->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                   );
      tmpTMVAReader->AddVariable( "D0",                    &fMVAVar_EleD0                       );
      tmpTMVAReader->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                    );
      tmpTMVAReader->AddVariable( "EOverP",                &fMVAVar_EleEOverP                   );
      tmpTMVAReader->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout     );
      tmpTMVAReader->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi            );
      tmpTMVAReader->AddVariable( "NBrem",                 &fMVAVar_EleNBrem                    );
      tmpTMVAReader->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP    );      
      tmpTMVAReader->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn      );
      tmpTMVAReader->AddVariable( "IP3d",                  &fMVAVar_EleIP3d                     );
      tmpTMVAReader->AddVariable( "IP3dSig",               &fMVAVar_EleIP3dSig                  );
    }
    if (type == kIDIsoCombined) {
      tmpTMVAReader->AddVariable( "SigmaIEtaIEta",         &fMVAVar_EleSigmaIEtaIEta            );
      tmpTMVAReader->AddVariable( "DEtaIn",                &fMVAVar_EleDEtaIn                   );
      tmpTMVAReader->AddVariable( "DPhiIn",                &fMVAVar_EleDPhiIn                   );
      tmpTMVAReader->AddVariable( "D0",                    &fMVAVar_EleD0                       );
      tmpTMVAReader->AddVariable( "FBrem",                 &fMVAVar_EleFBrem                    );
      tmpTMVAReader->AddVariable( "EOverP",                &fMVAVar_EleEOverP                   );
      tmpTMVAReader->AddVariable( "ESeedClusterOverPout",  &fMVAVar_EleESeedClusterOverPout     );
      tmpTMVAReader->AddVariable( "SigmaIPhiIPhi",         &fMVAVar_EleSigmaIPhiIPhi            );
      tmpTMVAReader->AddVariable( "OneOverEMinusOneOverP", &fMVAVar_EleOneOverEMinusOneOverP    );      
      tmpTMVAReader->AddVariable( "ESeedClusterOverPIn",   &fMVAVar_EleESeedClusterOverPIn      );
      tmpTMVAReader->AddVariable( "IP3d",                  &fMVAVar_EleIP3d                     );
      tmpTMVAReader->AddVariable( "IP3dSig",               &fMVAVar_EleIP3dSig                  );

      tmpTMVAReader->AddVariable( "GsfTrackChi2OverNdof",  &fMVAVar_EleGsfTrackChi2OverNdof     );
      tmpTMVAReader->AddVariable( "dEtaCalo",              &fMVAVar_EledEtaCalo                 );
      tmpTMVAReader->AddVariable( "dPhiCalo",              &fMVAVar_EledPhiCalo                 );
      tmpTMVAReader->AddVariable( "R9",                    &fMVAVar_EleR9                       );
      tmpTMVAReader->AddVariable( "SCEtaWidth",            &fMVAVar_EleSCEtaWidth               );
      tmpTMVAReader->AddVariable( "SCPhiWidth",            &fMVAVar_EleSCPhiWidth               );
      tmpTMVAReader->AddVariable( "CovIEtaIPhi",           &fMVAVar_EleCovIEtaIPhi              );
      if (i == 2 || i == 5) {
        tmpTMVAReader->AddVariable( "PreShowerOverRaw",      &fMVAVar_ElePreShowerOverRaw       );
      }
      tmpTMVAReader->AddVariable( "ChargedIso03",          &fMVAVar_EleChargedIso03OverPt       );
      tmpTMVAReader->AddVariable( "NeutralHadronIso03",    &fMVAVar_EleNeutralHadronIso03OverPt );
      tmpTMVAReader->AddVariable( "GammaIso03",            &fMVAVar_EleGammaIso03OverPt         );
      tmpTMVAReader->AddVariable( "ChargedIso04",          &fMVAVar_EleChargedIso04OverPt       );
      tmpTMVAReader->AddVariable( "NeutralHadronIso04",    &fMVAVar_EleNeutralHadronIso04OverPt );
      tmpTMVAReader->AddVariable( "GammaIso04",            &fMVAVar_EleGammaIso04OverPt         );

    }

    if (type == kIDEGamma2012TrigV0 || type == kIDHWW2012TrigV0) {
      // Pure tracking variables
      tmpTMVAReader->AddVariable("fbrem",           &fMVAVar_EleFBrem);
      tmpTMVAReader->AddVariable("kfchi2",          &fMVAVar_EleKFTrkChiSqr);
      tmpTMVAReader->AddVariable("kfhits",          &fMVAVar_EleKFTrkNLayers);  //Don't have this in (BAMBU <= 025)
      if(type == kIDEGamma2012TrigV0) 
         tmpTMVAReader->AddVariable("kfhitsall",       &fMVAVar_EleKFTrkNHits);
      tmpTMVAReader->AddVariable("gsfchi2",         &fMVAVar_EleGsfTrackChi2OverNdof);
      tmpTMVAReader->AddVariable("deta",            &fMVAVar_EleDEtaIn);
      tmpTMVAReader->AddVariable("dphi",            &fMVAVar_EleDPhiIn);
      tmpTMVAReader->AddVariable("detacalo",        &fMVAVar_EledEtaCalo);
      tmpTMVAReader->AddVariable("see",             &fMVAVar_EleSigmaIEtaIEta);
      tmpTMVAReader->AddVariable("spp",             &fMVAVar_EleSigmaIPhiIPhi);
      tmpTMVAReader->AddVariable("etawidth",        &fMVAVar_EleSCEtaWidth);
      tmpTMVAReader->AddVariable("phiwidth",        &fMVAVar_EleSCPhiWidth);
      tmpTMVAReader->AddVariable("e1x5e5x5",        &fMVAVar_EleE1x5OverE5x5);
      tmpTMVAReader->AddVariable("R9",              &fMVAVar_EleR9);
      tmpTMVAReader->AddVariable("HoE",             &fMVAVar_EleHoverE);
      tmpTMVAReader->AddVariable("EoP",             &fMVAVar_EleEOverP); 
      tmpTMVAReader->AddVariable("IoEmIoP",         &fMVAVar_EleOneOverEMinusOneOverP);
      tmpTMVAReader->AddVariable("eleEoPout",       &fMVAVar_EleEEleClusterOverPout); //Don't have this in (BAMBU <= 025)
      if(type == kIDEGamma2012TrigV0) 
        tmpTMVAReader->AddVariable("EoPout",          &fMVAVar_EleESeedClusterOverPout); 
      if (i == 2 || i == 5) {
        tmpTMVAReader->AddVariable( "PreShowerOverRaw",      &fMVAVar_ElePreShowerOverRaw       );
      }
      tmpTMVAReader->AddVariable( "d0",             &fMVAVar_EleD0);
      tmpTMVAReader->AddVariable( "ip3d",           &fMVAVar_EleIP3d);
    
      tmpTMVAReader->AddSpectator("eta",            &fMVAVar_EleEta);
      tmpTMVAReader->AddSpectator("pt",             &fMVAVar_ElePt);
    }

    if (type == kIDEGamma2012NonTrigV0 ) {
          // Pure tracking variables
      tmpTMVAReader->AddVariable("fbrem",           &fMVAVar_EleFBrem);
      tmpTMVAReader->AddVariable("kfchi2",          &fMVAVar_EleKFTrkChiSqr);
      tmpTMVAReader->AddVariable("kfhitsall",       &fMVAVar_EleKFTrkNHits);
      tmpTMVAReader->AddVariable("gsfchi2",         &fMVAVar_EleGsfTrackChi2OverNdof);
      tmpTMVAReader->AddVariable("deta",            &fMVAVar_EleDEtaIn);
      tmpTMVAReader->AddVariable("dphi",            &fMVAVar_EleDPhiIn);
      tmpTMVAReader->AddVariable("detacalo",        &fMVAVar_EledEtaCalo);
      tmpTMVAReader->AddVariable("see",             &fMVAVar_EleSigmaIEtaIEta);
      tmpTMVAReader->AddVariable("spp",             &fMVAVar_EleSigmaIPhiIPhi);
      tmpTMVAReader->AddVariable("etawidth",        &fMVAVar_EleSCEtaWidth);
      tmpTMVAReader->AddVariable("phiwidth",        &fMVAVar_EleSCPhiWidth);
      tmpTMVAReader->AddVariable("e1x5e5x5",        &fMVAVar_EleE1x5OverE5x5);
      tmpTMVAReader->AddVariable("R9",              &fMVAVar_EleR9);
      tmpTMVAReader->AddVariable("HoE",             &fMVAVar_EleHoverE);
      tmpTMVAReader->AddVariable("EoP",             &fMVAVar_EleEOverP); 
      tmpTMVAReader->AddVariable("IoEmIoP",         &fMVAVar_EleOneOverEMinusOneOverP);
      tmpTMVAReader->AddVariable("EoPout",          &fMVAVar_EleESeedClusterOverPout); 
      if (i==2 || i==5) {
        tmpTMVAReader->AddVariable("PreShowerOverRaw",&fMVAVar_ElePreShowerOverRaw);
      }
      tmpTMVAReader->AddSpectator("eta",            &fMVAVar_EleEta);
      tmpTMVAReader->AddSpectator("pt",             &fMVAVar_ElePt);
    }

    if (type == kIDEGamma2012NonTrigV1 ) {
          // Pure tracking variables
      tmpTMVAReader->AddVariable("fbrem",           &fMVAVar_EleFBrem);
      tmpTMVAReader->AddVariable("kfchi2",          &fMVAVar_EleKFTrkChiSqr);
      tmpTMVAReader->AddVariable("kfhits",          &fMVAVar_EleKFTrkNLayers);
      tmpTMVAReader->AddVariable("gsfchi2",         &fMVAVar_EleGsfTrackChi2OverNdof);
      tmpTMVAReader->AddVariable("deta",            &fMVAVar_EleDEtaIn);
      tmpTMVAReader->AddVariable("dphi",            &fMVAVar_EleDPhiIn);
      tmpTMVAReader->AddVariable("detacalo",        &fMVAVar_EledEtaCalo);
      tmpTMVAReader->AddVariable("see",             &fMVAVar_EleSigmaIEtaIEta);
      tmpTMVAReader->AddVariable("spp",             &fMVAVar_EleSigmaIPhiIPhi);
      tmpTMVAReader->AddVariable("etawidth",        &fMVAVar_EleSCEtaWidth);
      tmpTMVAReader->AddVariable("phiwidth",        &fMVAVar_EleSCPhiWidth);
      tmpTMVAReader->AddVariable("e1x5e5x5",        &fMVAVar_EleE1x5OverE5x5);
      tmpTMVAReader->AddVariable("R9",              &fMVAVar_EleR9);
      tmpTMVAReader->AddVariable("HoE",             &fMVAVar_EleHoverE);
      tmpTMVAReader->AddVariable("EoP",             &fMVAVar_EleEOverP); 
      tmpTMVAReader->AddVariable("IoEmIoP",         &fMVAVar_EleOneOverEMinusOneOverP);
      tmpTMVAReader->AddVariable("eleEoPout",       &fMVAVar_EleEEleClusterOverPout); 
      if (i==2 || i==5) {
        tmpTMVAReader->AddVariable("PreShowerOverRaw",&fMVAVar_ElePreShowerOverRaw);
      }
      tmpTMVAReader->AddSpectator("eta",            &fMVAVar_EleEta);
      tmpTMVAReader->AddSpectator("pt",             &fMVAVar_ElePt);
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
      tmpTMVAReader->AddSpectator( "eta",   &fMVAVar_EleEta );
      tmpTMVAReader->AddSpectator( "pt" ,   &fMVAVar_ElePt  );
    }

    if (type == kIDIsoCombinedHWW2012TrigV4) {

            // Pure tracking variables
      tmpTMVAReader->AddVariable("fbrem",                      &fMVAVar_EleFBrem);
      tmpTMVAReader->AddVariable("kfchi2",                     &fMVAVar_EleKFTrkChiSqr);
      tmpTMVAReader->AddVariable("kflayers",                   &fMVAVar_EleKFTrkNLayers);
      tmpTMVAReader->AddVariable("gsfchi2",                    &fMVAVar_EleGsfTrackChi2OverNdof);

      // Geometrical matchings
      tmpTMVAReader->AddVariable("deta",                       &fMVAVar_EleDEtaIn);
      tmpTMVAReader->AddVariable("dphi",                       &fMVAVar_EleDPhiIn);
      tmpTMVAReader->AddVariable("detacalo",                   &fMVAVar_EledEtaCalo);
    
      // Pure ECAL -> shower shapes
      tmpTMVAReader->AddVariable("see",                        &fMVAVar_EleSigmaIEtaIEta);
      tmpTMVAReader->AddVariable("spp",                        &fMVAVar_EleSigmaIPhiIPhi);
      tmpTMVAReader->AddVariable("etawidth",                   &fMVAVar_EleSCEtaWidth);
      tmpTMVAReader->AddVariable("phiwidth",                   &fMVAVar_EleSCPhiWidth);
      tmpTMVAReader->AddVariable("OneMinusSeedE1x5OverE5x5",   &fMVAVar_EleOneMinusE1x5OverE5x5);
      tmpTMVAReader->AddVariable("R9",                         &fMVAVar_EleR9);

      // Energy matching
      tmpTMVAReader->AddVariable("HoE",                        &fMVAVar_EleHoverE);
      tmpTMVAReader->AddVariable("EoP",                        &fMVAVar_EleEOverP); 
      tmpTMVAReader->AddVariable("IoEmIoP",                    &fMVAVar_EleOneOverEMinusOneOverP);
      tmpTMVAReader->AddVariable("EEleoPout",                  &fMVAVar_EleEEleClusterOverPout);
      if(i == 2 || i == 5) {
	tmpTMVAReader->AddVariable("PreShowerOverRaw",&fMVAVar_ElePreShowerOverRaw);
      }

      // IP
      tmpTMVAReader->AddVariable("d0",              &fMVAVar_EleD0);
      tmpTMVAReader->AddVariable("ip3d",            &fMVAVar_EleIP3d);

      //isolation variables
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
      tmpTMVAReader->AddVariable( "rho",                           &fMVAVar_Rho );

      //spectators
      tmpTMVAReader->AddSpectator("eta",            &fMVAVar_EleEta);
      tmpTMVAReader->AddSpectator("pt",             &fMVAVar_ElePt);

    }


    tmpTMVAReader->BookMVA(fMethodname , weightsfiles[i] );
    std::cout << "MVABin " << i << " : MethodName = " << fMethodname 
              << " , type == " << type << " , "
              << "Load weights file : " << weightsfiles[i] 
              << std::endl;
    fTMVAReader.push_back(tmpTMVAReader);

  }
  std::cout << "Electron ID MVA Completed\n";
}


//--------------------------------------------------------------------------------------------------
UInt_t ElectronIDMVA::GetMVABin( double eta, double pt) const {
  
    //Default is to return the first bin
    uint bin = 0;

    //return the first bin if not using binned version
    if (!fUseBinnedVersion) return 0;

    if (fMVAType == ElectronIDMVA::kBaseline 
        ||fMVAType == ElectronIDMVA::kNoIPInfo
        ||fMVAType == ElectronIDMVA::kWithIPInfo
        ||fMVAType == ElectronIDMVA::kIDIsoCombined) {
      if (pt < 20 && fabs(eta) < 1.0) bin = 0;
      if (pt < 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) bin = 1;
      if (pt < 20 && fabs(eta) >= 1.479) bin = 2;
      if (pt >= 20 && fabs(eta) < 1.0) bin = 3;
      if (pt >= 20 && fabs(eta) >= 1.0 && fabs(eta) < 1.479) bin = 4;
      if (pt >= 20 && fabs(eta) >= 1.479) bin = 5;
    }

    if (fMVAType == ElectronIDMVA::kIsoRingsV0) {
      if (pt < 10 && fabs(eta) < 1.479) bin = 0;
      if (pt < 10 && fabs(eta) >= 1.479) bin = 1;
      if (pt >= 10 && fabs(eta) < 1.479) bin = 2;
      if (pt >= 10 && fabs(eta) >= 1.479) bin = 3;
    }

    if (fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV0 ||
	fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1) {
      bin = 0;
      if (pt < 10 && fabs(eta) < 0.8) bin = 0;
      if (pt < 10 && fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) bin = 1;
      if (pt < 10 && fabs(eta) >= 1.479) bin = 2;
      if (pt >= 10 && fabs(eta) < 0.8) bin = 3;
      if (pt >= 10 && fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) bin = 4;
      if (pt >= 10 && fabs(eta) >= 1.479) bin = 5;
    }

    if (fMVAType == ElectronIDMVA::kIDEGamma2012TrigV0 || 
	fMVAType == ElectronIDMVA::kIDHWW2012TrigV0 ||
        fMVAType == ElectronIDMVA::kIDIsoCombinedHWW2012TrigV4
      ) {
      bin = 0;
      if (pt < 20 && fabs(eta) < 0.8) bin = 0;
      if (pt < 20 && fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) bin = 1;
      if (pt < 20 && fabs(eta) >= 1.479) bin = 2;
      if (pt >= 20 && fabs(eta) < 0.8) bin = 3;
      if (pt >= 20 && fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) bin = 4;
      if (pt >= 20 && fabs(eta) >= 1.479) bin = 5;
    }

    return bin;
}



//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(Double_t ElePt , Double_t EleEta,
                                 Double_t EleSigmaIEtaIEta,
                                 Double_t EleDEtaIn,
                                 Double_t EleDPhiIn,
                                 Double_t EleHoverE,
                                 Double_t EleD0,
                                 Double_t EleDZ,
                                 Double_t EleFBrem,
                                 Double_t EleEOverP,
                                 Double_t EleESeedClusterOverPout,
                                 Double_t EleSigmaIPhiIPhi,
                                 Double_t EleNBrem,
                                 Double_t EleOneOverEMinusOneOverP,
                                 Double_t EleESeedClusterOverPIn,
                                 Double_t EleIP3d,
                                 Double_t EleIP3dSig
  ) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  //set all input variables
  fMVAVar_EleSigmaIEtaIEta = EleSigmaIEtaIEta;
  fMVAVar_EleDEtaIn = EleDEtaIn;
  fMVAVar_EleDPhiIn = EleDPhiIn;
  fMVAVar_EleHoverE = EleHoverE;
  fMVAVar_EleD0 = EleD0;
  fMVAVar_EleDZ = EleDZ;
  fMVAVar_EleFBrem = EleFBrem;
  fMVAVar_EleEOverP = EleEOverP;
  fMVAVar_EleESeedClusterOverPout = EleESeedClusterOverPout;
  fMVAVar_EleSigmaIPhiIPhi = EleSigmaIPhiIPhi;
  fMVAVar_EleNBrem = EleNBrem;
  fMVAVar_EleOneOverEMinusOneOverP = EleOneOverEMinusOneOverP;
  fMVAVar_EleESeedClusterOverPIn = EleESeedClusterOverPIn;
  fMVAVar_EleIP3d = EleIP3d;
  fMVAVar_EleIP3dSig = EleIP3dSig;

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  reader = fTMVAReader[GetMVABin( EleEta, ElePt)];
                                                
  mva = reader->EvaluateMVA( fMethodname );

  return mva;
}

//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(Double_t ElePt , Double_t EleEta, Double_t PileupEnergyDensity,
                                 Double_t EleSigmaIEtaIEta,
                                 Double_t EleDEtaIn,
                                 Double_t EleDPhiIn,
                                 Double_t EleHoverE,
                                 Double_t EleD0,
                                 Double_t EleDZ,
                                 Double_t EleFBrem,
                                 Double_t EleEOverP,
                                 Double_t EleESeedClusterOverPout,
                                 Double_t EleSigmaIPhiIPhi,
                                 Double_t EleNBrem,
                                 Double_t EleOneOverEMinusOneOverP,
                                 Double_t EleESeedClusterOverPIn,
                                 Double_t EleIP3d,
                                 Double_t EleIP3dSig,
                                 Double_t EleGsfTrackChi2OverNdof,
                                 Double_t EledEtaCalo,
                                 Double_t EledPhiCalo,
                                 Double_t EleR9,
                                 Double_t EleSCEtaWidth,
                                 Double_t EleSCPhiWidth,
                                 Double_t EleCovIEtaIPhi,
                                 Double_t ElePreShowerOverRaw,
                                 Double_t EleChargedIso03,
                                 Double_t EleNeutralHadronIso03,
                                 Double_t EleGammaIso03,
                                 Double_t EleChargedIso04,
                                 Double_t EleNeutralHadronIso04,
                                 Double_t EleGammaIso04,
                                 Bool_t printDebug
  ) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Double_t Rho = 0;
  if (!(TMath::IsNaN(PileupEnergyDensity) || isinf(PileupEnergyDensity))) Rho = PileupEnergyDensity;

  //set all input variables
  fMVAVar_EleSigmaIEtaIEta = EleSigmaIEtaIEta;
  fMVAVar_EleDEtaIn = EleDEtaIn;
  fMVAVar_EleDPhiIn = EleDPhiIn;
  fMVAVar_EleHoverE = EleHoverE;
  fMVAVar_EleD0 = EleD0;
  fMVAVar_EleDZ = EleDZ;
  fMVAVar_EleFBrem = EleFBrem;
  fMVAVar_EleEOverP = EleEOverP;
  fMVAVar_EleESeedClusterOverPout = EleESeedClusterOverPout;
  fMVAVar_EleSigmaIPhiIPhi = EleSigmaIPhiIPhi;
  fMVAVar_EleNBrem = EleNBrem;
  fMVAVar_EleOneOverEMinusOneOverP = EleOneOverEMinusOneOverP;
  fMVAVar_EleESeedClusterOverPIn = EleESeedClusterOverPIn;
  fMVAVar_EleIP3d = EleIP3d;
  fMVAVar_EleIP3dSig = EleIP3dSig;
  fMVAVar_EleGsfTrackChi2OverNdof = EleGsfTrackChi2OverNdof;
  fMVAVar_EledEtaCalo = EledEtaCalo;
  fMVAVar_EledPhiCalo = EledPhiCalo;
  fMVAVar_EleR9 = EleR9;
  fMVAVar_EleSCEtaWidth = EleSCEtaWidth;
  fMVAVar_EleSCPhiWidth = EleSCPhiWidth;
  fMVAVar_EleCovIEtaIPhi = EleCovIEtaIPhi;
  fMVAVar_ElePreShowerOverRaw = ElePreShowerOverRaw;
  fMVAVar_EleChargedIso03OverPt 
    = (EleChargedIso03 
       - Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleChargedIso03, EleEta)) / ElePt;
  fMVAVar_EleNeutralHadronIso03OverPt 
    = (EleNeutralHadronIso03
       - Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso03, EleEta) 
       + Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso007,EleEta)) / ElePt;
  fMVAVar_EleGammaIso03OverPt 
    = (EleGammaIso03 
       - Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIso03, EleEta) 
       + Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIsoVetoEtaStrip03,EleEta))/ElePt;
  fMVAVar_EleChargedIso04OverPt 
    = (EleChargedIso04 
       - Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleChargedIso04, EleEta))/ElePt;
  fMVAVar_EleNeutralHadronIso04OverPt
    = (EleNeutralHadronIso04 
       - Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso04, EleEta) 
       + Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso007,EleEta))/ElePt;
  fMVAVar_EleGammaIso04OverPt 
    = (EleGammaIso04 
       - Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIso04, EleEta) 
       + Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIsoVetoEtaStrip04,EleEta))/ElePt;




  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  reader = fTMVAReader[GetMVABin( EleEta, ElePt)];
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug == kTRUE) {
    std::cout << "Debug Electron MVA: "
	 << ElePt << " " << EleEta << " " << " --> MVABin " << GetMVABin( EleEta, ElePt) << " : "     
	 << fMVAVar_EleSigmaIEtaIEta << " " 
	 << fMVAVar_EleDEtaIn << " " 
	 << fMVAVar_EleDPhiIn << " " 
	 << fMVAVar_EleHoverE << " " 
	 << fMVAVar_EleD0 << " " 
	 << fMVAVar_EleDZ << " " 
	 << fMVAVar_EleFBrem << " " 
	 << fMVAVar_EleEOverP << " " 
	 << fMVAVar_EleESeedClusterOverPout << " " 
	 << fMVAVar_EleSigmaIPhiIPhi << " " 
	 << fMVAVar_EleNBrem << " " 
	 << fMVAVar_EleOneOverEMinusOneOverP << " " 
	 << fMVAVar_EleESeedClusterOverPIn << " " 
	 << fMVAVar_EleIP3d << " " 
	 << fMVAVar_EleIP3dSig << " " 
	 << fMVAVar_EleGsfTrackChi2OverNdof << " "
	 << fMVAVar_EledEtaCalo << " "
	 << fMVAVar_EledPhiCalo << " "
	 << fMVAVar_EleR9 << " "
	 << fMVAVar_EleSCEtaWidth << " "
	 << fMVAVar_EleSCPhiWidth << " "
	 << fMVAVar_EleCovIEtaIPhi << " "
	 << fMVAVar_ElePreShowerOverRaw << " "
	 << fMVAVar_EleChargedIso03OverPt  << " "
	 << fMVAVar_EleNeutralHadronIso03OverPt  << " "
	 << fMVAVar_EleGammaIso03OverPt  << " "
	 << fMVAVar_EleChargedIso04OverPt  << " "
	 << fMVAVar_EleNeutralHadronIso04OverPt  << " "
	 << fMVAVar_EleGammaIso04OverPt  << " "
	 << " === : === "
	 << mva 
	 << std::endl;
  }

  return mva;
}

Double_t ElectronIDMVA::MVAValue_IsoRings( Double_t ElePt,
                            Double_t EleSCEta,
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

  if (fMVAType != ElectronIDMVA::kIsoRingsV0) {
    std::cout << "Error: This function is only supported for MVAType == kIsoRingsV0.\n" << std::endl;
    assert(kFALSE);
  }

  fMVAVar_ElePt = ElePt;
  fMVAVar_EleEta = EleSCEta;
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
    std::cout <<" -> BIN: " << fMVAVar_EleEta << " " << fMVAVar_ElePt << " : " << GetMVABin( fMVAVar_EleEta , fMVAVar_ElePt) << std::endl;
  }
  reader = fTMVAReader[GetMVABin( fMVAVar_EleEta , fMVAVar_ElePt)];                                              
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug == kTRUE) {

    std::cout << "Debug Electron MVA-ISO: ";
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

Double_t ElectronIDMVA::MVAValue_IDNonTrig( Double_t ElePt, 
                             Double_t EleSCEta, 
                             Double_t EleFBrem, 
                             Double_t EleKFTrkChiSqr,
                             Double_t EleKFTrkNHits,
                             Double_t EleGsfTrackChi2OverNdof,
                             Double_t EleDEtaIn, 
                             Double_t EleDPhiIn, 
                             Double_t EledEtaCalo,
                             Double_t EleSigmaIEtaIEta, 
                             Double_t EleSigmaIPhiIPhi, 
                             Double_t EleSCEtaWidth,
                             Double_t EleSCPhiWidth,
                             Double_t EleE1x5OverE5x5,
                             Double_t EleR9,
                             Double_t EleHoverE, 
                             Double_t EleEOverP, 
                             Double_t EleOneOverEMinusOneOverP, 
                             Double_t EleESeedClusterOverPout, 
                             Double_t ElePreShowerOverRaw,
                             Bool_t printDebug) {

  if (fMVAType != ElectronIDMVA::kIDEGamma2012NonTrigV0) {
    std::cout << "Error: This function is only supported for MVAType == kIDEGamma2012NonTrigV0.\n" << std::endl;
    assert(kFALSE);
  }

  fMVAVar_ElePt = ElePt; 
  fMVAVar_EleEta = EleSCEta; 
  fMVAVar_EleFBrem = EleFBrem; 
  fMVAVar_EleKFTrkChiSqr = EleKFTrkChiSqr;
  fMVAVar_EleKFTrkNHits = EleKFTrkNHits;
  fMVAVar_EleGsfTrackChi2OverNdof = EleGsfTrackChi2OverNdof;
  fMVAVar_EleDEtaIn = EleDEtaIn; 
  fMVAVar_EleDPhiIn = EleDPhiIn; 
  fMVAVar_EledEtaCalo = EledEtaCalo;
  fMVAVar_EleSigmaIEtaIEta = EleSigmaIEtaIEta; 
  fMVAVar_EleSigmaIPhiIPhi = EleSigmaIPhiIPhi; 
  fMVAVar_EleSCEtaWidth = EleSCEtaWidth;
  fMVAVar_EleSCPhiWidth = EleSCPhiWidth;
  fMVAVar_EleE1x5OverE5x5 = EleE1x5OverE5x5;
  fMVAVar_EleR9 = EleR9;
  fMVAVar_EleHoverE = EleHoverE; 
  fMVAVar_EleEOverP = EleEOverP; 
  fMVAVar_EleOneOverEMinusOneOverP = EleOneOverEMinusOneOverP; 
  fMVAVar_EleESeedClusterOverPout = EleESeedClusterOverPout; 
  fMVAVar_ElePreShowerOverRaw = ElePreShowerOverRaw;

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;

  if (printDebug == kTRUE) {
    std::cout <<" -> BIN: " << fMVAVar_EleEta << " " << fMVAVar_ElePt << " : " << GetMVABin( fMVAVar_EleEta , fMVAVar_ElePt) << std::endl;
  }
  reader = fTMVAReader[GetMVABin( fMVAVar_EleEta , fMVAVar_ElePt)];                                              
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug == kTRUE) {
    std::cout << "Debug Electron MVA: ";
    std::cout << " fbrem " <<  fMVAVar_EleFBrem  
              << " kfchi2 " << fMVAVar_EleKFTrkChiSqr  
              << " kfhits " << fMVAVar_EleKFTrkNLayers  
              << " kfhitsall " <<  fMVAVar_EleKFTrkNHits 
              << " gsfchi2 " << fMVAVar_EleGsfTrackChi2OverNdof  
              << " deta " <<  fMVAVar_EleDEtaIn  
              << " dphi " << fMVAVar_EleDPhiIn  
              << " detacalo " << fMVAVar_EledEtaCalo  
              << " see " << fMVAVar_EleSigmaIEtaIEta  
              << " spp " << fMVAVar_EleSigmaIPhiIPhi  
              << " etawidth " << fMVAVar_EleSCEtaWidth  
              << " phiwidth " << fMVAVar_EleSCPhiWidth  
              << " e1x5e5x5 " << fMVAVar_EleE1x5OverE5x5  
              << " R9 " << fMVAVar_EleR9  
              << " HoE " << fMVAVar_EleHoverE  
              << " EoP " << fMVAVar_EleEOverP  
              << " IoEmIoP " << fMVAVar_EleOneOverEMinusOneOverP  
              << " eleEoPout " << fMVAVar_EleESeedClusterOverPout  
              << " EoPout " << fMVAVar_EleESeedClusterOverPout  
              << " d0 " << fMVAVar_EleD0  
              << " ip3d " << fMVAVar_EleIP3d  
              << " eta " << fMVAVar_EleEta  
              << " pt " << fMVAVar_ElePt << std::endl;
    std::cout << "MVA: " << mva << " "    
              << std::endl;    
  }
  return mva;
}

//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(const Electron *ele, const Vertex *vertex, 
                                 const PFCandidateCol *PFCands, 
                                 const PileupEnergyDensityCol *PileupEnergyDensity,
                                 Double_t intRadius,
                                 Bool_t printDebug) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Double_t Rho = PileupEnergyDensity->At(0)->Rho(fRhoAlgo);

  //set all input variables
  fMVAVar_EleSigmaIEtaIEta = ele->CoviEtaiEta() ; 
  fMVAVar_EleDEtaIn = ele->DeltaEtaSuperClusterTrackAtVtx(); 
  fMVAVar_EleDPhiIn = ele->DeltaPhiSuperClusterTrackAtVtx(); 
  fMVAVar_EleHoverE = ele->HadronicOverEm(); 
  fMVAVar_EleD0 = ele->BestTrk()->D0Corrected(*vertex); 
  fMVAVar_EleDZ = ele->BestTrk()->DzCorrected(*vertex); 
  fMVAVar_EleFBrem = ele->FBrem(); 
  fMVAVar_EleEOverP = ele->ESuperClusterOverP(); 
  fMVAVar_EleESeedClusterOverPout = ele->ESeedClusterOverPout(); 
  if (!TMath::IsNaN(ele->SCluster()->Seed()->CoviPhiiPhi())) fMVAVar_EleSigmaIPhiIPhi = TMath::Sqrt(ele->SCluster()->Seed()->CoviPhiiPhi()); 
  else fMVAVar_EleSigmaIPhiIPhi = ele->CoviEtaiEta();
  fMVAVar_EleNBrem = ele->NumberOfClusters() - 1; 
  fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->EcalEnergy())) - 1.0 / ele->BestTrk()->P(); 
  fMVAVar_EleESeedClusterOverPIn = ele->ESeedClusterOverPIn(); 
  fMVAVar_EleIP3d = ele->Ip3dPV(); 
  fMVAVar_EleIP3dSig = ele->Ip3dPVSignificance(); 
  fMVAVar_EleGsfTrackChi2OverNdof = ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof();
  fMVAVar_EledEtaCalo =  ele->DeltaEtaSeedClusterTrackAtCalo();
  fMVAVar_EledPhiCalo = ele->DeltaPhiSeedClusterTrackAtCalo();
  fMVAVar_EleR9 = ele->SCluster()->R9();
  fMVAVar_EleSCEtaWidth = ele->SCluster()->EtaWidth();
  fMVAVar_EleSCPhiWidth = ele->SCluster()->PhiWidth();
  fMVAVar_EleCovIEtaIPhi = ele->SCluster()->Seed()->CoviEtaiPhi();
  fMVAVar_ElePreShowerOverRaw = ele->SCluster()->PreshowerEnergy() / ele->SCluster()->RawEnergy();
  fMVAVar_EleChargedIso03OverPt 
    = (IsolationTools::PFElectronIsolation(ele, PFCands, vertex, 0.1, 99999, 0.3, intRadius) 
       - Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleChargedIso03, ele->SCluster()->Eta())) / ele->Pt();
  fMVAVar_EleNeutralHadronIso03OverPt 
    = (IsolationTools::PFElectronIsolation(ele, PFCands, vertex, 0.1, 0.5, 0.3, intRadius, PFCandidate::eNeutralHadron) 
       - Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso03, ele->SCluster()->Eta()) 
       + Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso007,ele->SCluster()->Eta())) / ele->Pt();
  fMVAVar_EleGammaIso03OverPt 
    = (IsolationTools::PFElectronIsolation(ele, PFCands, vertex, 0.1, 0.5, 0.3, intRadius, PFCandidate::eGamma) 
       - Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIso03, ele->SCluster()->Eta()) 
       + Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIsoVetoEtaStrip03,ele->SCluster()->Eta())) / ele->Pt();
  fMVAVar_EleChargedIso04OverPt 
    = (IsolationTools::PFElectronIsolation(ele, PFCands, vertex, 0.1, 99999, 0.4, intRadius) 
       - Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleChargedIso04, ele->SCluster()->Eta())) / ele->Pt();
  fMVAVar_EleNeutralHadronIso04OverPt 
    = (IsolationTools::PFElectronIsolation(ele, PFCands, vertex, 0.1, 0.5, 0.4, intRadius, PFCandidate::eNeutralHadron) 
       - Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso04, ele->SCluster()->Eta()) 
       + Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIso007,ele->SCluster()->Eta())) / ele->Pt() ;
  fMVAVar_EleGammaIso04OverPt 
    = (IsolationTools::PFElectronIsolation(ele, PFCands, vertex, 0.1, 0.5, 0.4, intRadius, PFCandidate::eGamma) 
       - Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIso04, ele->SCluster()->Eta()) 
       + Rho * ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIsoVetoEtaStrip04,ele->SCluster()->Eta())) / ele->Pt();
  
  //Additional vars
  fMVAVar_EleEEleClusterOverPout = ele->EEleClusterOverPout();
  if (ele->TrackerTrk()) {
    fMVAVar_EleKFTrkChiSqr = ele->TrackerTrk()->RChi2();
    fMVAVar_EleKFTrkNHits = ele->TrackerTrk()->NHits();
  } else {
    fMVAVar_EleKFTrkChiSqr = -1;
    fMVAVar_EleKFTrkNHits = 0;
  }
  fMVAVar_EleE1x5OverE5x5 = ele->SCluster()->Seed()->E1x5() / ele->SCluster()->Seed()->E5x5();


  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  reader = fTMVAReader[GetMVABin( ele->SCluster()->Eta(), ele->Pt())];
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug == kTRUE) {
    std::cout << "Debug Electron MVA: "
              << ele->Pt() << " " << ele->Eta() << " " << ele->Phi() << " : "
              << ele->Pt() << " " << ele->SCluster()->AbsEta() << " --> MVABin " << GetMVABin( ele->SCluster()->Eta(), ele->Pt()) << " : "     
              << fMVAVar_EleSigmaIEtaIEta << " " 
              << fMVAVar_EleDEtaIn << " " 
              << fMVAVar_EleDPhiIn << " " 
              << fMVAVar_EleHoverE << " " 
              << fMVAVar_EleD0 << " " 
              << fMVAVar_EleDZ << " " 
              << fMVAVar_EleFBrem << " " 
              << fMVAVar_EleEOverP << " " 
              << fMVAVar_EleESeedClusterOverPout << " " 
              << fMVAVar_EleSigmaIPhiIPhi << " " 
              << fMVAVar_EleNBrem << " " 
              << fMVAVar_EleOneOverEMinusOneOverP << " " 
              << fMVAVar_EleESeedClusterOverPIn << " " 
              << fMVAVar_EleIP3d << " " 
              << fMVAVar_EleIP3dSig << " " 
              << fMVAVar_EleGsfTrackChi2OverNdof << " "
              << fMVAVar_EledEtaCalo << " "
              << fMVAVar_EledPhiCalo << " "
              << fMVAVar_EleR9 << " "
              << fMVAVar_EleSCEtaWidth << " "
              << fMVAVar_EleSCPhiWidth << " "
              << fMVAVar_EleCovIEtaIPhi << " "
              << fMVAVar_ElePreShowerOverRaw << " "
              << fMVAVar_EleKFTrkChiSqr  << " "
              << fMVAVar_EleKFTrkNHits  << " "
              << fMVAVar_EleE1x5OverE5x5  << " "
              << " ::: "
      
              << " === : === "
              << mva << " "    
              << std::endl;
    
  }

  return mva;
}

//--------------------------------------------------------------------------------------------------
Double_t ElectronIDMVA::MVAValue(const Electron *ele, const Vertex *vertex,
                                 Bool_t printDebug) {
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  //set all input variables
  fMVAVar_ElePt = ele->Pt();
  fMVAVar_EleEta = ele->SCluster()->Eta();
  fMVAVar_EleSigmaIEtaIEta = ele->CoviEtaiEta() ;

  if (fMVAType == ElectronIDMVA::kIDEGamma2012TrigV0 || 
      fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV0 ||
      fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1 ||
      fMVAType == ElectronIDMVA::kIDHWW2012TrigV0 
    ) {
    fMVAVar_EleDEtaIn = TMath::Min(fabs(double(ele->DeltaEtaSuperClusterTrackAtVtx())),0.06); ; 
    fMVAVar_EleDPhiIn = TMath::Min(fabs(double(ele->DeltaPhiSuperClusterTrackAtVtx())),0.6); 
    fMVAVar_EleFBrem = TMath::Max(double(ele->FBrem()),-1.0); 
    fMVAVar_EleEOverP = TMath::Min(double(ele->ESuperClusterOverP()), 20.0); 
    fMVAVar_EleESeedClusterOverPout = TMath::Min(double(ele->ESeedClusterOverPout()),20.0); 
    fMVAVar_EleEEleClusterOverPout = TMath::Min(double(ele->EEleClusterOverPout()),20.0); 
    fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->EcalEnergy())) - 1.0 / ele->P(); 
    fMVAVar_EleGsfTrackChi2OverNdof = TMath::Min(double( ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof()),200.0);
    fMVAVar_EledEtaCalo =  TMath::Min(fabs(double(ele->DeltaEtaSeedClusterTrackAtCalo())),0.2);
    fMVAVar_EleR9 = TMath::Min(double(ele->SCluster()->R9()), 5.0);   
  }
  else if (fMVAType == ElectronIDMVA::kIDIsoCombinedHWW2012TrigV4) {
    fMVAVar_EleDEtaIn = TMath::Min(fabs(double(ele->DeltaEtaSuperClusterTrackAtVtx())),0.06); ; 
    fMVAVar_EleDPhiIn = ele->DeltaPhiSuperClusterTrackAtVtx();
    fMVAVar_EleFBrem = TMath::Max(double(ele->FBrem()),-1.0); 
    fMVAVar_EleEOverP = TMath::Min(double(ele->ESuperClusterOverP()), 20.0); 
    fMVAVar_EleESeedClusterOverPout = TMath::Min(double(ele->ESeedClusterOverPout()),20.0); 
    fMVAVar_EleEEleClusterOverPout = TMath::Min(double(ele->EEleClusterOverPout()),20.0); 
    fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->EcalEnergy())) - 1.0 / ele->PIn(); 
    fMVAVar_EleGsfTrackChi2OverNdof = TMath::Min(double( ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof()),200.0);
    fMVAVar_EledEtaCalo =  ele->DeltaEtaSeedClusterTrackAtCalo();
    fMVAVar_EleR9 = TMath::Min(double(ele->SCluster()->R9()), 5.0);   
  }
  else {
    fMVAVar_EleDEtaIn = ele->DeltaEtaSuperClusterTrackAtVtx();  
    fMVAVar_EleDPhiIn = ele->DeltaPhiSuperClusterTrackAtVtx(); 
    fMVAVar_EleFBrem = ele->FBrem(); 
    fMVAVar_EleEOverP = ele->ESuperClusterOverP(); 
    fMVAVar_EleESeedClusterOverPout = ele->ESeedClusterOverPout(); 
    fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->EcalEnergy())) - 1.0 / ele->BestTrk()->P();
    fMVAVar_EleGsfTrackChi2OverNdof = ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof();
    fMVAVar_EledEtaCalo =  ele->DeltaEtaSeedClusterTrackAtCalo();
    fMVAVar_EleR9 = ele->SCluster()->R9();   
  }

  fMVAVar_EleHoverE = ele->HadronicOverEm(); 
  fMVAVar_EleD0 = ele->BestTrk()->D0Corrected(*vertex); 
  fMVAVar_EleDZ = ele->BestTrk()->DzCorrected(*vertex); 
  if (!TMath::IsNaN(ele->SCluster()->Seed()->CoviPhiiPhi())) fMVAVar_EleSigmaIPhiIPhi = TMath::Sqrt(ele->SCluster()->Seed()->CoviPhiiPhi()); 
  else {
    if (fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1 ) {
      fMVAVar_EleSigmaIPhiIPhi = 0;
    } else {
      fMVAVar_EleSigmaIPhiIPhi = ele->CoviEtaiEta();
    }
  }

  fMVAVar_EleNBrem = ele->NumberOfClusters() - 1; 
  fMVAVar_EleESeedClusterOverPIn = ele->ESeedClusterOverPIn(); 
  fMVAVar_EleIP3d = ele->Ip3dPV(); 
  fMVAVar_EleIP3dSig = ele->Ip3dPVSignificance(); 
  fMVAVar_EledPhiCalo = ele->DeltaPhiSeedClusterTrackAtCalo();
  fMVAVar_EleSCEtaWidth = ele->SCluster()->EtaWidth();
  fMVAVar_EleSCPhiWidth = ele->SCluster()->PhiWidth();
  fMVAVar_EleCovIEtaIPhi = ele->SCluster()->Seed()->CoviEtaiPhi();
  fMVAVar_ElePreShowerOverRaw = ele->SCluster()->PreshowerEnergy() / ele->SCluster()->RawEnergy();

  //Additional vars
  if (ele->TrackerTrk()) {
    if (fMVAType == ElectronIDMVA::kIDEGamma2012TrigV0 || 
        fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV0 ||
        fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1 ||
	fMVAType == ElectronIDMVA::kIDHWW2012TrigV0 ||
        fMVAType == ElectronIDMVA::kIDIsoCombinedHWW2012TrigV4  
      ) {
      fMVAVar_EleKFTrkChiSqr = TMath::Min(double(ele->TrackerTrk()->RChi2()),10.0);
    } else {
      fMVAVar_EleKFTrkChiSqr = ele->TrackerTrk()->RChi2();
    }
    fMVAVar_EleKFTrkNHits = ele->TrackerTrk()->NHits();
    fMVAVar_EleKFTrkNLayers = ele->CTFTrkNLayersWithMeasurement();
  } else {
    fMVAVar_EleKFTrkChiSqr = 0;
    fMVAVar_EleKFTrkNHits = -1;
    fMVAVar_EleKFTrkNLayers = -1;
  }
  
  if( ele->SCluster()->Seed()->E5x5() > 0.0 ) {
    if (fMVAType == ElectronIDMVA::kIDEGamma2012TrigV0 || 
        fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV0 ||
        fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1 ||
	fMVAType == ElectronIDMVA::kIDHWW2012TrigV0 ||
        fMVAType == ElectronIDMVA::kIDIsoCombinedHWW2012TrigV4  
      ) {
      fMVAVar_EleE1x5OverE5x5 = TMath::Min(TMath::Max(1 - double(ele->SCluster()->Seed()->E1x5()/ele->SCluster()->Seed()->E5x5()) , -1.0),2.0);
      fMVAVar_EleOneMinusE1x5OverE5x5 = TMath::Min(TMath::Max(1 - double(ele->SCluster()->Seed()->E1x5()/ele->SCluster()->Seed()->E5x5()) , -1.0),2.0);
    } else {
      fMVAVar_EleE1x5OverE5x5 = ele->SCluster()->Seed()->E1x5()/ele->SCluster()->Seed()->E5x5();
    }
  } else {
    fMVAVar_EleE1x5OverE5x5 = -1.0;
  }

  //Do Binding of MVA input variables
  if (   fMVAType == ElectronIDMVA::kIDEGamma2012TrigV0 
      || fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV0 
      || fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1 
      || fMVAType == ElectronIDMVA::kIsoRingsV0
      || fMVAType == ElectronIDMVA::kIDHWW2012TrigV0
    ) {
    bindVariables();
  }


  //***********************************************************************
  // Si [Oct 26 , 2012]
  // Code Below was used for the 2011 HWW analysis
  // But we deprecate it now 
  //***********************************************************************
//   fMVAVar_ElePt = ele->Pt();
//   fMVAVar_EleEta = ele->Eta();

//   //set all input variables
//   fMVAVar_EleSigmaIEtaIEta = ele->CoviEtaiEta() ; 
//   fMVAVar_EleDEtaIn = ele->DeltaEtaSuperClusterTrackAtVtx(); 
//   fMVAVar_EleDPhiIn = ele->DeltaPhiSuperClusterTrackAtVtx(); 
//   fMVAVar_EleHoverE = ele->HadronicOverEm(); 
//   fMVAVar_EleD0 = ele->BestTrk()->D0Corrected(*vertex); 
//   fMVAVar_EleDZ = ele->BestTrk()->DzCorrected(*vertex); 
//   fMVAVar_EleFBrem = ele->FBrem(); 
//   fMVAVar_EleEOverP = ele->ESuperClusterOverP(); 
//   fMVAVar_EleESeedClusterOverPout = ele->ESeedClusterOverPout(); 
//   if (!TMath::IsNaN(ele->SCluster()->Seed()->CoviPhiiPhi())) fMVAVar_EleSigmaIPhiIPhi = TMath::Sqrt(ele->SCluster()->Seed()->CoviPhiiPhi()); 
//   else fMVAVar_EleSigmaIPhiIPhi = 0;
//   fMVAVar_EleNBrem = ele->NumberOfClusters() - 1; 
//   fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->EcalEnergy())) - 1.0 / ele->BestTrk()->P(); 
//   fMVAVar_EleESeedClusterOverPIn = ele->ESeedClusterOverPIn(); 
//   fMVAVar_EleIP3d = ele->Ip3dPV(); 
//   fMVAVar_EleIP3dSig = ele->Ip3dPVSignificance(); 

//   fMVAVar_EleEEleClusterOverPout = ele->EEleClusterOverPout();
//   if (ele->TrackerTrk()) {
//     fMVAVar_EleKFTrkChiSqr = ele->TrackerTrk()->RChi2();
//     fMVAVar_EleKFTrkNHits = ele->TrackerTrk()->NHits();
//     fMVAVar_EleKFTrkNLayers = ele->CTFTrkNLayersWithMeasurement();
//   } else {
//     fMVAVar_EleKFTrkChiSqr = -1;
//     fMVAVar_EleKFTrkNHits = 0;
//     fMVAVar_EleKFTrkNLayers = 0;
//   }
//   fMVAVar_EleGsfTrackChi2OverNdof = ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof();
//   fMVAVar_EledEtaCalo =  ele->DeltaEtaSeedClusterTrackAtCalo();
//   fMVAVar_EleSCEtaWidth = ele->SCluster()->EtaWidth();
//   fMVAVar_EleSCPhiWidth = ele->SCluster()->PhiWidth();
//   fMVAVar_EleE1x5OverE5x5 = ele->SCluster()->Seed()->E1x5() / ele->SCluster()->Seed()->E5x5();
//   fMVAVar_EleR9 = ele->SCluster()->R9();
//   fMVAVar_EleHoverE = ele->HadronicOverEm(); 
//   fMVAVar_EleEOverP = ele->ESuperClusterOverP(); 
//   fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->EcalEnergy())) - 1.0 / ele->BestTrk()->P(); 
//   fMVAVar_EleR9 = ele->SCluster()->R9();
//   fMVAVar_ElePreShowerOverRaw = ele->SCluster()->PreshowerEnergy() / ele->SCluster()->RawEnergy();
    

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  reader = fTMVAReader[GetMVABin( ele->SCluster()->Eta(), ele->Pt())];                                              
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug == kTRUE) {
    std::cout << "Debug Electron MVA: "
              << ele->Pt() << " " << ele->Eta() << " " << ele->Phi() << " : "
              << ele->Pt() << " " << ele->SCluster()->AbsEta() << " --> MVABin " << GetMVABin( ele->SCluster()->Eta(), ele->Pt()) << " : "     
              << fMVAVar_EleSigmaIEtaIEta << " " 
              << fMVAVar_EleDEtaIn << " " 
              << fMVAVar_EleDPhiIn << " " 
              << fMVAVar_EleHoverE << " " 
              << fMVAVar_EleD0 << " " 
              << fMVAVar_EleDZ << " " 
              << fMVAVar_EleFBrem << " " 
              << fMVAVar_EleEOverP << " " 
              << fMVAVar_EleESeedClusterOverPout << " " 
              << fMVAVar_EleSigmaIPhiIPhi << " " 
              << fMVAVar_EleNBrem << " " 
              << fMVAVar_EleOneOverEMinusOneOverP << " " 
              << fMVAVar_EleESeedClusterOverPIn << " " 
              << fMVAVar_EleIP3d << " " 
              << fMVAVar_EleIP3dSig << " " 
              << fMVAVar_EleGsfTrackChi2OverNdof << " "
              << fMVAVar_EledEtaCalo << " "
              << fMVAVar_EledPhiCalo << " "
              << fMVAVar_EleR9 << " "
              << fMVAVar_EleSCEtaWidth << " "
              << fMVAVar_EleSCPhiWidth << " "
              << fMVAVar_EleCovIEtaIPhi << " "
              << fMVAVar_ElePreShowerOverRaw << " "
              << fMVAVar_EleKFTrkChiSqr  << " "
              << fMVAVar_EleKFTrkNHits  << " "
              << fMVAVar_EleE1x5OverE5x5  << " "
              << " === : === "
              << mva << " "    
              << std::endl;
    
  }



  return mva;
}


//--------------------------------------------------------------------------------------------------
//MVA Includes Isolation with removal of other leptons
//
Double_t ElectronIDMVA::MVAValue(const Electron *ele, const Vertex *vertex, 
                                 const PFCandidateCol *PFCands, 
                                 const PileupEnergyDensityCol *PileupEnergyDensity,
                                 ElectronTools::EElectronEffectiveAreaTarget EffectiveAreaTarget,
                                 const ElectronCol *goodElectrons/* = 0*/,
                                 const MuonCol *goodMuons/* = 0*/,
                                 Bool_t printDebug/* = kFALSE*/)
{
  
  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Double_t Rho = PileupEnergyDensity->At(0)->Rho(fRhoAlgo);

  //set all input variables
  fMVAVar_ElePt = ele->Pt();
  fMVAVar_EleEta = ele->SCluster()->Eta();
  fMVAVar_EleSigmaIEtaIEta = ele->CoviEtaiEta() ;

  if (fMVAType == ElectronIDMVA::kIDEGamma2012TrigV0 || 
      fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV0 ||
      fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1 ||
      fMVAType == ElectronIDMVA::kIDHWW2012TrigV0) {
    fMVAVar_EleDEtaIn = TMath::Min(fabs(double(ele->DeltaEtaSuperClusterTrackAtVtx())),0.06); ; 
    fMVAVar_EleDPhiIn = TMath::Min(fabs(double(ele->DeltaPhiSuperClusterTrackAtVtx())),0.6); 
    fMVAVar_EleFBrem = TMath::Max(double(ele->FBrem()),-1.0); 
    fMVAVar_EleEOverP = TMath::Min(double(ele->ESuperClusterOverP()), 20.0); 
    fMVAVar_EleESeedClusterOverPout = TMath::Min(double(ele->ESeedClusterOverPout()),20.0); 
    fMVAVar_EleEEleClusterOverPout = TMath::Min(double(ele->EEleClusterOverPout()),20.0); 
    fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->EcalEnergy())) - 1.0 / ele->P(); 
    fMVAVar_EleGsfTrackChi2OverNdof = TMath::Min(double( ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof()),200.0);
    fMVAVar_EledEtaCalo =  TMath::Min(fabs(double(ele->DeltaEtaSeedClusterTrackAtCalo())),0.2);
    fMVAVar_EleR9 = TMath::Min(double(ele->SCluster()->R9()), 5.0);   
  } else {
    fMVAVar_EleDEtaIn = ele->DeltaEtaSuperClusterTrackAtVtx();  
    fMVAVar_EleDPhiIn = ele->DeltaPhiSuperClusterTrackAtVtx(); 
    fMVAVar_EleFBrem = ele->FBrem(); 
    fMVAVar_EleEOverP = ele->ESuperClusterOverP(); 
    fMVAVar_EleESeedClusterOverPout = ele->ESeedClusterOverPout(); 
    fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->EcalEnergy())) - 1.0 / ele->BestTrk()->P();
    fMVAVar_EleGsfTrackChi2OverNdof = ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof();
    fMVAVar_EledEtaCalo =  ele->DeltaEtaSeedClusterTrackAtCalo();
    fMVAVar_EleR9 = ele->SCluster()->R9();   
  }

  fMVAVar_EleHoverE = ele->HadronicOverEm(); 
  fMVAVar_EleD0 = ele->BestTrk()->D0Corrected(*vertex); 
  fMVAVar_EleDZ = ele->BestTrk()->DzCorrected(*vertex); 
  if (!TMath::IsNaN(ele->SCluster()->Seed()->CoviPhiiPhi())) fMVAVar_EleSigmaIPhiIPhi = TMath::Sqrt(ele->SCluster()->Seed()->CoviPhiiPhi()); 
  else {
    if (fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1) {
      fMVAVar_EleSigmaIPhiIPhi = 0;
    } else {
      fMVAVar_EleSigmaIPhiIPhi = ele->CoviEtaiEta();
    }
  }

  fMVAVar_EleNBrem = ele->NumberOfClusters() - 1; 
  fMVAVar_EleESeedClusterOverPIn = ele->ESeedClusterOverPIn(); 
  fMVAVar_EleIP3d = ele->Ip3dPV(); 
  fMVAVar_EleIP3dSig = ele->Ip3dPVSignificance(); 
  fMVAVar_EledPhiCalo = ele->DeltaPhiSeedClusterTrackAtCalo();
  fMVAVar_EleSCEtaWidth = ele->SCluster()->EtaWidth();
  fMVAVar_EleSCPhiWidth = ele->SCluster()->PhiWidth();
  fMVAVar_EleCovIEtaIPhi = ele->SCluster()->Seed()->CoviEtaiPhi();
  fMVAVar_ElePreShowerOverRaw = ele->SCluster()->PreshowerEnergy() / ele->SCluster()->RawEnergy();

  //Additional vars
  if (ele->TrackerTrk()) {
    if (fMVAType == ElectronIDMVA::kIDEGamma2012TrigV0 || 
        fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV0 ||
        fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1 ||
	fMVAType == ElectronIDMVA::kIDHWW2012TrigV0 ) {
      fMVAVar_EleKFTrkChiSqr = TMath::Min(double(ele->TrackerTrk()->RChi2()),10.0);
    } else {
      fMVAVar_EleKFTrkChiSqr = ele->TrackerTrk()->RChi2();
    }
    fMVAVar_EleKFTrkNHits = ele->TrackerTrk()->NHits();
    fMVAVar_EleKFTrkNLayers = ele->CTFTrkNLayersWithMeasurement();
  } else {
    fMVAVar_EleKFTrkChiSqr = 0;
    fMVAVar_EleKFTrkNHits = -1;
    fMVAVar_EleKFTrkNLayers = -1;
  }
  
  if( ele->SCluster()->Seed()->E5x5() > 0.0 ) {
    if (fMVAType == ElectronIDMVA::kIDEGamma2012TrigV0 || 
        fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV0 ||
        fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1 ||
	fMVAType == ElectronIDMVA::kIDHWW2012TrigV0 ) {
      fMVAVar_EleE1x5OverE5x5 = TMath::Min(TMath::Max(1 - double(ele->SCluster()->Seed()->E1x5()/ele->SCluster()->Seed()->E5x5()) , -1.0),2.0);
    } else {
      fMVAVar_EleE1x5OverE5x5 = ele->SCluster()->Seed()->E1x5()/ele->SCluster()->Seed()->E5x5();
    }
  } else {
    fMVAVar_EleE1x5OverE5x5 = -1.0;
  }


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

  for (UInt_t p=0; p<PFCands->GetEntries();p++) {   
    const PFCandidate *pf = PFCands->At(p);
      
    //exclude the electron itself
    if(pf->GsfTrk() && ele->GsfTrk() &&
       pf->GsfTrk() == ele->GsfTrk()) continue;
    if(pf->TrackerTrk() && ele->TrackerTrk() &&
       pf->TrackerTrk() == ele->TrackerTrk()) continue;      

    //************************************************************
    // New Isolation Calculations
    //************************************************************
    Double_t dr = MathUtils::DeltaR(ele->Mom(), pf->Mom());

    if (dr < 1.0) {
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
	   if (!(fabs(pf->BestTrk()->DzCorrected(*vertex) - ele->BestTrk()->DzCorrected(*vertex)) < 0.2)) passVeto = kFALSE;
	   //************************************************************
	   // Veto any PFmuon, or PFEle
	   if (pf->PFType() == PFCandidate::eElectron || pf->PFType() == PFCandidate::eMuon) passVeto = kFALSE;
	   //************************************************************
	   //************************************************************
	   // Footprint Veto
	   if (fabs(ele->SCluster()->Eta()) >= 1.479 && dr < 0.015) passVeto = kFALSE;
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
	   //************************************************************
	   // Footprint Veto
	   if (fabs(ele->SCluster()->Eta()) >= 1.479) {
             if (dr < 0.08) passVeto = kFALSE;
	   }
	   //************************************************************
	   
	   if (passVeto) {
	     if (dr < 0.1) tmpGammaIso_DR0p0To0p1 += pf->Pt();
	     if (dr >= 0.1 && dr < 0.2) tmpGammaIso_DR0p1To0p2 += pf->Pt();
	     if (dr >= 0.2 && dr < 0.3) tmpGammaIso_DR0p2To0p3 += pf->Pt();
	     if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += pf->Pt();
	     if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += pf->Pt();
	   }
	 }
	 //NeutralHadron
	 else {
           if (dr < 0.1) tmpNeutralHadronIso_DR0p0To0p1 += pf->Pt();
           if (dr >= 0.1 && dr < 0.2) tmpNeutralHadronIso_DR0p1To0p2 += pf->Pt();
           if (dr >= 0.2 && dr < 0.3) tmpNeutralHadronIso_DR0p2To0p3 += pf->Pt();
           if (dr >= 0.3 && dr < 0.4) tmpNeutralHadronIso_DR0p3To0p4 += pf->Pt();
           if (dr >= 0.4 && dr < 0.5) tmpNeutralHadronIso_DR0p4To0p5 += pf->Pt();
	 }
      } //not lepton footprint
    } //in 1.0 dr cone
  } //loop over PF candidates

  fMVAVar_ChargedIso_DR0p0To0p1   = TMath::Min((tmpChargedIso_DR0p0To0p1)/ele->Pt(), 2.5);
  fMVAVar_ChargedIso_DR0p1To0p2   = TMath::Min((tmpChargedIso_DR0p1To0p2)/ele->Pt(), 2.5);
  fMVAVar_ChargedIso_DR0p2To0p3 = TMath::Min((tmpChargedIso_DR0p2To0p3)/ele->Pt(), 2.5);
  fMVAVar_ChargedIso_DR0p3To0p4 = TMath::Min((tmpChargedIso_DR0p3To0p4)/ele->Pt(), 2.5);
  fMVAVar_ChargedIso_DR0p4To0p5 = TMath::Min((tmpChargedIso_DR0p4To0p5)/ele->Pt(), 2.5); 
  fMVAVar_GammaIso_DR0p0To0p1 = TMath::Max(TMath::Min((tmpGammaIso_DR0p0To0p1 - Rho*ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIsoDR0p0To0p1, ele->SCluster()->Eta(), EffectiveAreaTarget))/ele->Pt(), 2.5), 0.0);
  fMVAVar_GammaIso_DR0p1To0p2 = TMath::Max(TMath::Min((tmpGammaIso_DR0p1To0p2 - Rho*ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIsoDR0p1To0p2, ele->SCluster()->Eta(), EffectiveAreaTarget))/ele->Pt(), 2.5), 0.0);
  fMVAVar_GammaIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpGammaIso_DR0p2To0p3 - Rho*ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIsoDR0p2To0p3, ele->SCluster()->Eta(), EffectiveAreaTarget))/ele->Pt(), 2.5), 0.0);
  fMVAVar_GammaIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpGammaIso_DR0p3To0p4 - Rho*ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIsoDR0p3To0p4, ele->SCluster()->Eta(), EffectiveAreaTarget))/ele->Pt(), 2.5), 0.0);
  fMVAVar_GammaIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpGammaIso_DR0p4To0p5 - Rho*ElectronTools::ElectronEffectiveArea(ElectronTools::kEleGammaIsoDR0p4To0p5, ele->SCluster()->Eta(), EffectiveAreaTarget))/ele->Pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p0To0p1 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p0To0p1 - Rho*ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIsoDR0p0To0p1, ele->SCluster()->Eta(), EffectiveAreaTarget))/ele->Pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p1To0p2 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p1To0p2 - Rho*ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIsoDR0p1To0p2, ele->SCluster()->Eta(), EffectiveAreaTarget))/ele->Pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p2To0p3 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p2To0p3 - Rho*ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIsoDR0p2To0p3, ele->SCluster()->Eta(), EffectiveAreaTarget))/ele->Pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p3To0p4 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p3To0p4 - Rho*ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIsoDR0p3To0p4, ele->SCluster()->Eta(), EffectiveAreaTarget))/ele->Pt(), 2.5), 0.0);
  fMVAVar_NeutralHadronIso_DR0p4To0p5 = TMath::Max(TMath::Min((tmpNeutralHadronIso_DR0p4To0p5 - Rho*ElectronTools::ElectronEffectiveArea(ElectronTools::kEleNeutralHadronIsoDR0p4To0p5, ele->SCluster()->Eta(), EffectiveAreaTarget))/ele->Pt(), 2.5), 0.0); 

  //Do Binding of MVA input variables
  if (   fMVAType == ElectronIDMVA::kIDEGamma2012TrigV0 
      || fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV0 
      || fMVAType == ElectronIDMVA::kIsoRingsV0
      || fMVAType == ElectronIDMVA::kIDHWW2012TrigV0) {
    bindVariables();
  }

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  
  reader = fTMVAReader[GetMVABin( fMVAVar_EleEta , fMVAVar_ElePt)];                                              
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug == kTRUE) {

    std::cout << "Debug Electron MVA-ID: "
              << fMVAVar_ElePt<< " " << fMVAVar_EleEta << " " 
	      << " --> MVABin " << GetMVABin( fMVAVar_EleEta , fMVAVar_ElePt) << " : "	  
              << " fbrem " <<  fMVAVar_EleFBrem  
              << " kfchi2 " << fMVAVar_EleKFTrkChiSqr  
              << " kfhits " << fMVAVar_EleKFTrkNLayers
              << " kfhitsall " << fMVAVar_EleKFTrkNHits
              << " gsfchi2 " << fMVAVar_EleGsfTrackChi2OverNdof  
              << " deta " <<  fMVAVar_EleDEtaIn  
              << " dphi " << fMVAVar_EleDPhiIn  
              << " detacalo " << fMVAVar_EledEtaCalo  
              << " see " << fMVAVar_EleSigmaIEtaIEta  
              << " spp " << fMVAVar_EleSigmaIPhiIPhi  
              << " etawidth " << fMVAVar_EleSCEtaWidth  
              << " phiwidth " << fMVAVar_EleSCPhiWidth  
              << " e1x5e5x5 " << fMVAVar_EleE1x5OverE5x5  
              << " R9 " << fMVAVar_EleR9  
              << " HoE " << fMVAVar_EleHoverE  
              << " EoP " << fMVAVar_EleEOverP  
              << " IoEmIoP " << fMVAVar_EleOneOverEMinusOneOverP  
              << " eleEoPout " << fMVAVar_EleEEleClusterOverPout  
              << " EoPout " << fMVAVar_EleESeedClusterOverPout  
	      << " PreShowerOverRaw" << fMVAVar_ElePreShowerOverRaw  
              << " d0 " << fMVAVar_EleD0  
              << " ip3d " << fMVAVar_EleIP3d  
              << " eta " << fMVAVar_EleEta  
              << " pt " << fMVAVar_ElePt
              << " === : === "
              << mva << " "    
              << std::endl;
    std::cout << "Debug Electron MVA-ISO: "
              << fMVAVar_ChargedIso_DR0p0To0p1 << " "
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
  }

  return mva;
}


//--------------------------------------------------------------------------------------------------
//MVA Includes Isolation 
//
Double_t ElectronIDMVA::MVAValue(const Electron *ele, const Vertex *vertex, 
                                 const VertexCol *primaryVertices,
                                 const PFCandidateCol *PFCands, 
                                 const PileupEnergyDensityCol *PileupEnergyDensity,
                                 ElectronTools::EElectronEffectiveAreaTarget EffectiveAreaTarget,
                                 Bool_t printDebug) {

  if (!fIsInitialized) { 
    std::cout << "Error: ElectronIDMVA not properly initialized.\n"; 
    return -9999;
  }

  Double_t Rho = PileupEnergyDensity->At(0)->Rho(fRhoAlgo);

  //set all input variables
  fMVAVar_ElePt = ele->Pt();
  fMVAVar_EleEta = ele->SCluster()->Eta();
  fMVAVar_EleSigmaIEtaIEta = ele->CoviEtaiEta() ;

  if (fMVAType == ElectronIDMVA::kIDEGamma2012TrigV0 || 
      fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV0 ||
      fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1 ||
      fMVAType == ElectronIDMVA::kIDHWW2012TrigV0 
    ) {
    fMVAVar_EleDEtaIn = TMath::Min(fabs(double(ele->DeltaEtaSuperClusterTrackAtVtx())),0.06); ; 
    fMVAVar_EleDPhiIn = TMath::Min(fabs(double(ele->DeltaPhiSuperClusterTrackAtVtx())),0.6); 
    fMVAVar_EleFBrem = TMath::Max(double(ele->FBrem()),-1.0); 
    fMVAVar_EleEOverP = TMath::Min(double(ele->ESuperClusterOverP()), 20.0); 
    fMVAVar_EleESeedClusterOverPout = TMath::Min(double(ele->ESeedClusterOverPout()),20.0); 
    fMVAVar_EleEEleClusterOverPout = TMath::Min(double(ele->EEleClusterOverPout()),20.0); 
    fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->EcalEnergy())) - 1.0 / ele->P(); 
    fMVAVar_EleGsfTrackChi2OverNdof = TMath::Min(double( ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof()),200.0);
    fMVAVar_EledEtaCalo =  TMath::Min(fabs(double(ele->DeltaEtaSeedClusterTrackAtCalo())),0.2);
    fMVAVar_EleR9 = TMath::Min(double(ele->SCluster()->R9()), 5.0);   
  }
  else if (fMVAType == ElectronIDMVA::kIDIsoCombinedHWW2012TrigV4) {
    fMVAVar_EleDEtaIn = TMath::Min(fabs(double(ele->DeltaEtaSuperClusterTrackAtVtx())),0.06); ; 
    fMVAVar_EleDPhiIn = ele->DeltaPhiSuperClusterTrackAtVtx();
    fMVAVar_EleFBrem = TMath::Max(double(ele->FBrem()),-1.0); 
    fMVAVar_EleEOverP = TMath::Min(double(ele->ESuperClusterOverP()), 20.0); 
    fMVAVar_EleESeedClusterOverPout = TMath::Min(double(ele->ESeedClusterOverPout()),20.0); 
    fMVAVar_EleEEleClusterOverPout = TMath::Min(double(ele->EEleClusterOverPout()),20.0); 
    fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->EcalEnergy())) - 1.0 / ele->PIn(); 
    fMVAVar_EleGsfTrackChi2OverNdof = TMath::Min(double( ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof()),200.0);
    fMVAVar_EledEtaCalo =  ele->DeltaEtaSeedClusterTrackAtCalo();
    fMVAVar_EleR9 = TMath::Min(double(ele->SCluster()->R9()), 5.0);   
  }
  else {
    fMVAVar_EleDEtaIn = ele->DeltaEtaSuperClusterTrackAtVtx();  
    fMVAVar_EleDPhiIn = ele->DeltaPhiSuperClusterTrackAtVtx(); 
    fMVAVar_EleFBrem = ele->FBrem(); 
    fMVAVar_EleEOverP = ele->ESuperClusterOverP(); 
    fMVAVar_EleESeedClusterOverPout = ele->ESeedClusterOverPout(); 
    fMVAVar_EleOneOverEMinusOneOverP = (1.0/(ele->EcalEnergy())) - 1.0 / ele->BestTrk()->P();
    fMVAVar_EleGsfTrackChi2OverNdof = ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof();
    fMVAVar_EledEtaCalo =  ele->DeltaEtaSeedClusterTrackAtCalo();
    fMVAVar_EleR9 = ele->SCluster()->R9();   
  }

  fMVAVar_EleHoverE = ele->HadronicOverEm(); 
  fMVAVar_EleD0 = ele->BestTrk()->D0Corrected(*vertex); 
  fMVAVar_EleDZ = ele->BestTrk()->DzCorrected(*vertex); 
  if (!TMath::IsNaN(ele->SCluster()->Seed()->CoviPhiiPhi())) fMVAVar_EleSigmaIPhiIPhi = TMath::Sqrt(ele->SCluster()->Seed()->CoviPhiiPhi()); 
  else {
    if (fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1 ) {
      fMVAVar_EleSigmaIPhiIPhi = 0;
    } else {
      fMVAVar_EleSigmaIPhiIPhi = ele->CoviEtaiEta();
    }
  }

  fMVAVar_EleNBrem = ele->NumberOfClusters() - 1; 
  fMVAVar_EleESeedClusterOverPIn = ele->ESeedClusterOverPIn(); 
  fMVAVar_EleIP3d = ele->Ip3dPV(); 
  fMVAVar_EleIP3dSig = ele->Ip3dPVSignificance(); 
  fMVAVar_EledPhiCalo = ele->DeltaPhiSeedClusterTrackAtCalo();
  fMVAVar_EleSCEtaWidth = ele->SCluster()->EtaWidth();
  fMVAVar_EleSCPhiWidth = ele->SCluster()->PhiWidth();
  fMVAVar_EleCovIEtaIPhi = ele->SCluster()->Seed()->CoviEtaiPhi();
  fMVAVar_ElePreShowerOverRaw = ele->SCluster()->PreshowerEnergy() / ele->SCluster()->RawEnergy();

  //Additional vars
  if (ele->TrackerTrk()) {
    if (fMVAType == ElectronIDMVA::kIDEGamma2012TrigV0 || 
        fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV0 ||
        fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1 ||
	fMVAType == ElectronIDMVA::kIDHWW2012TrigV0 ||
        fMVAType == ElectronIDMVA::kIDIsoCombinedHWW2012TrigV4  
      ) {
      fMVAVar_EleKFTrkChiSqr = TMath::Min(double(ele->TrackerTrk()->RChi2()),10.0);
    } else {
      fMVAVar_EleKFTrkChiSqr = ele->TrackerTrk()->RChi2();
    }
    fMVAVar_EleKFTrkNHits = ele->TrackerTrk()->NHits();
    fMVAVar_EleKFTrkNLayers = ele->CTFTrkNLayersWithMeasurement();
  } else {
    fMVAVar_EleKFTrkChiSqr = 0;
    fMVAVar_EleKFTrkNHits = -1;
    fMVAVar_EleKFTrkNLayers = -1;
  }
  
  if( ele->SCluster()->Seed()->E5x5() > 0.0 ) {
    if (fMVAType == ElectronIDMVA::kIDEGamma2012TrigV0 || 
        fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV0 ||
        fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV1 ||
	fMVAType == ElectronIDMVA::kIDHWW2012TrigV0 ||
        fMVAType == ElectronIDMVA::kIDIsoCombinedHWW2012TrigV4  
      ) {
      fMVAVar_EleE1x5OverE5x5 = TMath::Min(TMath::Max(1 - double(ele->SCluster()->Seed()->E1x5()/ele->SCluster()->Seed()->E5x5()) , -1.0),2.0);
      fMVAVar_EleOneMinusE1x5OverE5x5 = TMath::Min(TMath::Max(1 - double(ele->SCluster()->Seed()->E1x5()/ele->SCluster()->Seed()->E5x5()) , -1.0),2.0);
    } else {
      fMVAVar_EleE1x5OverE5x5 = ele->SCluster()->Seed()->E1x5()/ele->SCluster()->Seed()->E5x5();
    }
  } else {
    fMVAVar_EleE1x5OverE5x5 = -1.0;
  }


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

  for (UInt_t p=0; p<PFCands->GetEntries();p++) {   
    const PFCandidate *pf = PFCands->At(p);
      
    //************************************************************
    // New Isolation Calculations
    //************************************************************
    Double_t dr = MathUtils::DeltaR(ele->Mom(), pf->Mom());
    
    if (dr < 0.5) {
      
      Bool_t passVeto = kTRUE;
      //Charged
      if(pf->BestTrk()) {
        
        //*************************************************
        //Use only PFNoPU
        //*************************************************
        Bool_t isPFNoPU = kFALSE;
        if(pf->PFType() == PFCandidate::eHadron) {
          if(pf->HasTrackerTrk() && 
             primaryVertices->At(0)->HasTrack(pf->TrackerTrk()) &&
             primaryVertices->At(0)->TrackWeight(pf->TrackerTrk()) > 0) {
            isPFNoPU = kTRUE;
          } else { 
               
            Bool_t vertexFound = kFALSE;
            const Vertex *closestVtx = 0;
            Double_t dzmin = 10000;
               
            // loop over vertices
            for(UInt_t j = 0; j < primaryVertices->GetEntries(); j++) {
              const Vertex *vtx = primaryVertices->At(j);
              assert(vtx);
                 
              if(pf->HasTrackerTrk() && 
                 vtx->HasTrack(pf->TrackerTrk()) &&
                 vtx->TrackWeight(pf->TrackerTrk()) > 0) {
                vertexFound = kTRUE;
                closestVtx = vtx;
                break;
              }
              Double_t dz = fabs(pf->SourceVertex().Z() - vtx->Z());
              if(dz < dzmin) {
                closestVtx = vtx;
                dzmin = dz;
              }
            }
               
            Bool_t fCheckClosestZVertex = kTRUE; //we use option 1
            if(fCheckClosestZVertex) {
              // Fallback: if track is not associated with any vertex,
              // associate it with the vertex closest in z
              if(vertexFound || closestVtx != vertex) {
                isPFNoPU = kFALSE;
              } else {
                isPFNoPU = kTRUE;
              }
            } else {
              if(vertexFound && closestVtx != vertex) {
                isPFNoPU = kFALSE;
              } else {
                isPFNoPU = kTRUE;
              }
            }
          } //hadron & trk stuff
        } else { // hadron
          //
          isPFNoPU = kTRUE;
        }  
        if (!isPFNoPU) continue;

        //************************************************************
        // Veto any PFmuon, or PFEle
        if (pf->PFType() == PFCandidate::eElectron || pf->PFType() == PFCandidate::eMuon) passVeto = kFALSE;
        //************************************************************
        //************************************************************
        // Footprint Veto
        if (fabs(ele->SCluster()->Eta()) >= 1.479 && dr < 0.015) passVeto = kFALSE;
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
        //************************************************************
        // Footprint Veto
        if (fabs(ele->SCluster()->Eta()) > 1.479) {
          if (dr < 0.08) passVeto = kFALSE;
        }
        //************************************************************
	   
        if (passVeto) {
          if (dr < 0.1) tmpGammaIso_DR0p0To0p1 += pf->Pt();
          if (dr >= 0.1 && dr < 0.2) tmpGammaIso_DR0p1To0p2 += pf->Pt();
          if (dr >= 0.2 && dr < 0.3) tmpGammaIso_DR0p2To0p3 += pf->Pt();
          if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += pf->Pt();
          if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += pf->Pt();
        }
      }
      //NeutralHadron
      else {
        if (dr < 0.1) tmpNeutralHadronIso_DR0p0To0p1 += pf->Pt();
        if (dr >= 0.1 && dr < 0.2) tmpNeutralHadronIso_DR0p1To0p2 += pf->Pt();
        if (dr >= 0.2 && dr < 0.3) tmpNeutralHadronIso_DR0p2To0p3 += pf->Pt();
        if (dr >= 0.3 && dr < 0.4) tmpNeutralHadronIso_DR0p3To0p4 += pf->Pt();
        if (dr >= 0.4 && dr < 0.5) tmpNeutralHadronIso_DR0p4To0p5 += pf->Pt();
      }
    } //in 1.0 dr cone
  } //loop over PF candidates

  fMVAVar_ChargedIso_DR0p0To0p1   = tmpChargedIso_DR0p0To0p1/ele->Pt();
  fMVAVar_ChargedIso_DR0p1To0p2   = tmpChargedIso_DR0p1To0p2/ele->Pt();
  fMVAVar_ChargedIso_DR0p2To0p3 = tmpChargedIso_DR0p2To0p3/ele->Pt();
  fMVAVar_ChargedIso_DR0p3To0p4 = tmpChargedIso_DR0p3To0p4/ele->Pt();
  fMVAVar_ChargedIso_DR0p4To0p5 = tmpChargedIso_DR0p4To0p5/ele->Pt();
  fMVAVar_GammaIso_DR0p0To0p1 = tmpGammaIso_DR0p0To0p1/ele->Pt();
  fMVAVar_GammaIso_DR0p1To0p2 = tmpGammaIso_DR0p1To0p2/ele->Pt();
  fMVAVar_GammaIso_DR0p2To0p3 = tmpGammaIso_DR0p2To0p3/ele->Pt();
  fMVAVar_GammaIso_DR0p3To0p4 = tmpGammaIso_DR0p3To0p4/ele->Pt();
  fMVAVar_GammaIso_DR0p4To0p5 = tmpGammaIso_DR0p4To0p5/ele->Pt();
  fMVAVar_NeutralHadronIso_DR0p0To0p1 = tmpNeutralHadronIso_DR0p0To0p1/ele->Pt();
  fMVAVar_NeutralHadronIso_DR0p1To0p2 = tmpNeutralHadronIso_DR0p1To0p2/ele->Pt();
  fMVAVar_NeutralHadronIso_DR0p2To0p3 = tmpNeutralHadronIso_DR0p2To0p3/ele->Pt();
  fMVAVar_NeutralHadronIso_DR0p3To0p4 = tmpNeutralHadronIso_DR0p3To0p4/ele->Pt();
  fMVAVar_NeutralHadronIso_DR0p4To0p5 = tmpNeutralHadronIso_DR0p4To0p5/ele->Pt();
  fMVAVar_Rho = Rho;

  //Do Binding of MVA input variables
  if (   fMVAType == ElectronIDMVA::kIDEGamma2012TrigV0 
      || fMVAType == ElectronIDMVA::kIDEGamma2012NonTrigV0 
      || fMVAType == ElectronIDMVA::kIsoRingsV0
      || fMVAType == ElectronIDMVA::kIDHWW2012TrigV0
    ) {
    bindVariables();
  }

  Double_t mva = -9999;  
  TMVA::Reader *reader = 0;
  
  reader = fTMVAReader[GetMVABin( fMVAVar_EleEta , fMVAVar_ElePt)];                                              
  mva = reader->EvaluateMVA( fMethodname );

  if (printDebug == kTRUE) {

    std::cout << "Debug Electron MVA-ID: "
              << fMVAVar_ElePt<< " " << fMVAVar_EleEta << " " 
	      << " --> MVABin " << GetMVABin( fMVAVar_EleEta , fMVAVar_ElePt) << " : "	  
              << " fbrem " <<  fMVAVar_EleFBrem  
              << " kfchi2 " << fMVAVar_EleKFTrkChiSqr  
              << " kfhits " << fMVAVar_EleKFTrkNLayers
              << " kfhitsall " << fMVAVar_EleKFTrkNHits
              << " gsfchi2 " << fMVAVar_EleGsfTrackChi2OverNdof  
              << " deta " <<  fMVAVar_EleDEtaIn  
              << " dphi " << fMVAVar_EleDPhiIn  
              << " detacalo " << fMVAVar_EledEtaCalo  
              << " see " << fMVAVar_EleSigmaIEtaIEta  
              << " spp " << fMVAVar_EleSigmaIPhiIPhi  
              << " etawidth " << fMVAVar_EleSCEtaWidth  
              << " phiwidth " << fMVAVar_EleSCPhiWidth  
              << " e1x5e5x5 " << fMVAVar_EleOneMinusE1x5OverE5x5  
              << " R9 " << fMVAVar_EleR9  
              << " HoE " << fMVAVar_EleHoverE  
              << " EoP " << fMVAVar_EleEOverP  
              << " IoEmIoP " << fMVAVar_EleOneOverEMinusOneOverP  
              << " eleEoPout " << fMVAVar_EleEEleClusterOverPout  
              << " EoPout " << fMVAVar_EleESeedClusterOverPout  
	      << " PreShowerOverRaw" << fMVAVar_ElePreShowerOverRaw  
              << " d0 " << fMVAVar_EleD0  
              << " ip3d " << fMVAVar_EleIP3d  
              << " eta " << fMVAVar_EleEta  
              << " pt " << fMVAVar_ElePt
              << " === : === "
              << mva << " "    
              << std::endl;
    std::cout << "Debug Electron MVA-ISO: "
              << fMVAVar_ChargedIso_DR0p0To0p1 << " "
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
              << fMVAVar_Rho << " "
              << std::endl;
  }

  return mva;
}


void ElectronIDMVA::bindVariables() {

  // this binding is needed for variables that sometime diverge. 

  if(fMVAVar_EleFBrem < -1.)
    fMVAVar_EleFBrem = -1.;	
  
  fMVAVar_EleDEtaIn = fabs(fMVAVar_EleDEtaIn);
  if(fMVAVar_EleDEtaIn > 0.06)
    fMVAVar_EleDEtaIn = 0.06;
  
  
  fMVAVar_EleDPhiIn = fabs(fMVAVar_EleDPhiIn);
  if(fMVAVar_EleDPhiIn > 0.6)
    fMVAVar_EleDPhiIn = 0.6;
  
  
  if(fMVAVar_EleESeedClusterOverPout > 20.)
    fMVAVar_EleESeedClusterOverPout = 20.;
  
  if(fMVAVar_EleEOverP > 20.)
    fMVAVar_EleEOverP = 20.;
  
  if(fMVAVar_EleEEleClusterOverPout > 20.)
    fMVAVar_EleEEleClusterOverPout = 20.;
  
  
  fMVAVar_EledEtaCalo = fabs(fMVAVar_EledEtaCalo);
  if(fMVAVar_EledEtaCalo > 0.2)
    fMVAVar_EledEtaCalo = 0.2;
  
  
  if(fMVAVar_EleE1x5OverE5x5 < -1.)
    fMVAVar_EleE1x5OverE5x5 = -1;
  
  if(fMVAVar_EleE1x5OverE5x5 > 2.)
    fMVAVar_EleE1x5OverE5x5 = 2.; 
  
  
  
  if(fMVAVar_EleR9 > 5)
    fMVAVar_EleR9 = 5;
  
  if(fMVAVar_EleGsfTrackChi2OverNdof > 200.)
    fMVAVar_EleGsfTrackChi2OverNdof = 200;
  
  
  if(fMVAVar_EleKFTrkChiSqr > 10.)
    fMVAVar_EleKFTrkChiSqr = 10.;
  
  // Needed for a bug in CMSSW_420, fixed in more recent CMSSW versions
  if(std::isnan(fMVAVar_EleSigmaIPhiIPhi))
    fMVAVar_EleSigmaIPhiIPhi = 0.;	
  
  
  return;
}
