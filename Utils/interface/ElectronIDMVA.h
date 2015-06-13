//--------------------------------------------------------------------------------------------------
// ElectronIDMVA
//
// Helper Class for Electron Identification MVA
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_ElectronIDMVA_H
#define MITPHYSICS_UTILS_ElectronIDMVA_H

#include "MitAna/DataTree/interface/ElectronFwd.h"
#include "MitAna/DataTree/interface/VertexFwd.h"
#include "MitAna/DataTree/interface/TrackFwd.h"
#include "MitAna/DataTree/interface/Electron.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"

class TRandom3;
namespace TMVA {
  class Reader;
}

namespace mithep {
  class ElectronIDMVA {
    public:
      ElectronIDMVA();
      virtual ~ElectronIDMVA(); 

      enum MVAType {
        kUninitialized = 0,
        kBaseline,           // SigmaIEtaIEta, DEtaIn, DPhiIn, FBrem, SigmaIPhiIPhi, NBrem, OneOverEMinusOneOverP
        kNoIPInfo,                  // kBaseline + EOverP, ESeedClusterOverPout, ESeedClusterOverPIn
        kWithIPInfo,                // kBaseline + d0 , IP3d, IP3dSig
        kIDIsoCombined,             // ID variables , PFIso03 , PFIso04
        kIDEGamma2012TrigV0,        // EGamma certified (Spring 2012) ID-only MVA
        kIDEGamma2012NonTrigV0,     // EGamma certified (Spring 2012) ID-only MVA
        kIDEGamma2012NonTrigV1,     // EGamma certified (Spring 2012) ID-only MVA, "official" version
        kIsoRingsV0,                // Isolation MVA with IsoRings as input
        kIDHWW2012TrigV0,           // HWW certified (Spring 2012) ID-only MVA
        kIDIsoCombinedHWW2012TrigV4 // HWW ID+Iso Combined for trig ele
      };


      void     Initialize(std::string const& methodName,
                          std::string const& weightsfile,
                          ElectronIDMVA::MVAType,
			  UInt_t = mithep::PileupEnergyDensity::kHighEta);
      void     Initialize(std::string const& methodName,
                          ElectronIDMVA::MVAType,
                          Bool_t useBinnedVersion,
                          std::vector<std::string> const& weightsfiles,
			  UInt_t = mithep::PileupEnergyDensity::kHighEta);
      void     Initialize(TString const& methodName,
                          TString const& Subdet0Pt10To20Weights , 
                          TString const& Subdet1Pt10To20Weights , 
                          TString const& Subdet2Pt10To20Weights,
                          TString const& Subdet0Pt20ToInfWeights, 
                          TString const& Subdet1Pt20ToInfWeights, 
                          TString const& Subdet2Pt20ToInfWeights,
                          ElectronIDMVA::MVAType,
			  UInt_t = mithep::PileupEnergyDensity::kHighEta);
      
      Bool_t   IsInitialized() const { return fIsInitialized; }
      void     bindVariables();
      UInt_t   GetMVABin(double eta,double pt ) const;

      Double_t MVAValue(const Electron *ele, const Vertex *vertex, Bool_t printDebug = kFALSE);
      Double_t MVAValue(const Electron *ele, const Vertex *vertex, 
                        const PFCandidateCol *PFCands, 
                        const PileupEnergyDensityCol *PileupEnergyDensity,
                        Double_t intRadius,
                        Bool_t printDebug = kFALSE);
      Double_t MVAValue(const Electron *ele, const Vertex *vertex, 
                        const PFCandidateCol *PFCands, 
                        const PileupEnergyDensityCol *PileupEnergyDensity,
                        ElectronTools::EElectronEffectiveAreaTarget EffectiveAreaTarget,
                        const ElectronCol *goodElectrons = 0,
                        const MuonCol *goodMuons = 0,
                        Bool_t printDebug = kFALSE);
      Double_t MVAValue(const Electron *ele, const Vertex *vertex,
                        const VertexCol *primaryVertices,
                        const PFCandidateCol *PFCands,
                        const PileupEnergyDensityCol *PileupEnergyDensity,
                        ElectronTools::EElectronEffectiveAreaTarget EffectiveAreaTarget,
                        Bool_t printDebug = kFALSE);
      Double_t MVAValue(Double_t ElePt , Double_t EleEta,
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
                        Double_t EleIP3dSig);
      Double_t MVAValue(Double_t ElePt , Double_t EleEta,
                        Double_t PileupEnergyDensity,
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
                        Bool_t printDebug = kFALSE );
      Double_t MVAValue_IsoRings( Double_t ElePt,
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
                                  Bool_t printDebug = kFALSE);
      Double_t MVAValue_IDNonTrig( Double_t ElePt, 
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
                                   Bool_t printDebug = kFALSE);

    protected:      
      std::vector<TMVA::Reader*> fTMVAReader;
      TString                   fMethodname;
      Bool_t                    fIsInitialized;
      MVAType                   fMVAType;
      Bool_t                    fUseBinnedVersion;
      UInt_t                    fNMVABins;
      UInt_t                    fRhoAlgo;

      Float_t                   fMVAVar_ElePt; 
      Float_t                   fMVAVar_EleEta; 
      Float_t                   fMVAVar_EleSigmaIEtaIEta; 
      Float_t                   fMVAVar_EleDEtaIn; 
      Float_t                   fMVAVar_EleDPhiIn; 
      Float_t                   fMVAVar_EleHoverE; 
      Float_t                   fMVAVar_EleD0; 
      Float_t                   fMVAVar_EleDZ; 
      Float_t                   fMVAVar_EleFBrem; 
      Float_t                   fMVAVar_EleEOverP; 
      Float_t                   fMVAVar_EleESeedClusterOverPout; 
      Float_t                   fMVAVar_EleSigmaIPhiIPhi; 
      Float_t                   fMVAVar_EleNBrem; 
      Float_t                   fMVAVar_EleOneOverEMinusOneOverP; 
      Float_t                   fMVAVar_EleESeedClusterOverPIn; 
      Float_t                   fMVAVar_EleIP3d; 
      Float_t                   fMVAVar_EleIP3dSig; 
      Float_t                   fMVAVar_EleGsfTrackChi2OverNdof;
      Float_t                   fMVAVar_EledEtaCalo;
      Float_t                   fMVAVar_EledPhiCalo;
      Float_t                   fMVAVar_EleR9;
      Float_t                   fMVAVar_EleSCEtaWidth;
      Float_t                   fMVAVar_EleSCPhiWidth;
      Float_t                   fMVAVar_EleCovIEtaIPhi;
      Float_t                   fMVAVar_ElePreShowerOverRaw;
      Float_t                   fMVAVar_EleChargedIso03OverPt;
      Float_t                   fMVAVar_EleNeutralHadronIso03OverPt;
      Float_t                   fMVAVar_EleGammaIso03OverPt;
      Float_t                   fMVAVar_EleChargedIso04OverPt;
      Float_t                   fMVAVar_EleNeutralHadronIso04OverPt;
      Float_t                   fMVAVar_EleGammaIso04OverPt;

      Float_t                   fMVAVar_EleEEleClusterOverPout;
      Float_t                   fMVAVar_EleKFTrkChiSqr;
      Float_t                   fMVAVar_EleKFTrkNHits;
      Float_t                   fMVAVar_EleKFTrkNLayers;
      Float_t                   fMVAVar_EleE1x5OverE5x5;
      Float_t                   fMVAVar_EleOneMinusE1x5OverE5x5;

      Float_t                   fMVAVar_ChargedIso_DR0p0To0p1;
      Float_t                   fMVAVar_ChargedIso_DR0p1To0p2;
      Float_t                   fMVAVar_ChargedIso_DR0p2To0p3;
      Float_t                   fMVAVar_ChargedIso_DR0p3To0p4;
      Float_t                   fMVAVar_ChargedIso_DR0p4To0p5;
      Float_t                   fMVAVar_GammaIso_DR0p0To0p1;
      Float_t                   fMVAVar_GammaIso_DR0p1To0p2;
      Float_t                   fMVAVar_GammaIso_DR0p2To0p3;
      Float_t                   fMVAVar_GammaIso_DR0p3To0p4;
      Float_t                   fMVAVar_GammaIso_DR0p4To0p5;
      Float_t                   fMVAVar_NeutralHadronIso_DR0p0To0p1;
      Float_t                   fMVAVar_NeutralHadronIso_DR0p1To0p2;
      Float_t                   fMVAVar_NeutralHadronIso_DR0p2To0p3;
      Float_t                   fMVAVar_NeutralHadronIso_DR0p3To0p4;
      Float_t                   fMVAVar_NeutralHadronIso_DR0p4To0p5;

      Float_t                   fMVAVar_Rho;

    ClassDef(ElectronIDMVA, 1) // Muon tools
      };
}

#endif
