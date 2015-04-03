 // $Id $

#include "MitPhysics/SelMods/interface/HwwExampleAnalysisMod.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/MetTools.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "TFile.h"
#include "TTree.h"

using namespace mithep;
ClassImp(mithep::HwwExampleAnalysisMod)

//--------------------------------------------------------------------------------------------------
HwwExampleAnalysisMod::HwwExampleAnalysisMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fMuonBranchName(Names::gkMuonBrn),
  fMetName("NoDefaultNameSet"),
  fCleanJetsNoPtCutName("NoDefaultNameSet"),
  fVertexName(ModNames::gkGoodVertexesName),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fMuons(0),
  fMet(0),
  fVertices(0),
  fPFCandidates(0),
  fPFJetName0("AKt5PFJets"),
  fPFJet0(0),
  fNEventsSelected(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void HwwExampleAnalysisMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void HwwExampleAnalysisMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  // Load Branches
  ReqBranch(fMuonBranchName,   fMuons);
  ReqBranch(fPFCandidatesName, fPFCandidates);
  ReqBranch(fPFJetName0,       fPFJet0);

  //Create your histograms here

  //*************************************************************************************************
  // Selection Histograms
  //*************************************************************************************************
  AddTH1(fHWWSelection,"hHWWSelection", ";Cut Number;Number of Events",             17, -1.5, 15.5);
  AddTH1(fHWWToEESelection,"hHWWToEESelection", ";Cut Number;Number of Events",     17, -1.5, 15.5);
  AddTH1(fHWWToMuMuSelection,"hHWWToMuMuSelection", ";Cut Number;Number of Events", 17, -1.5, 15.5);
  AddTH1(fHWWToEMuSelection,"hHWWToEMuSelection", ";Cut Number;Number of Events",   17, -1.5, 15.5);
  AddTH1(fHWWToMuESelection,"hHWWToMuESelection", ";Cut Number;Number of Events",   17, -1.5, 15.5);

  //***********************************************************************************************
  // Histograms after preselection
  //***********************************************************************************************
  AddTH1(fLeptonEta          ,"hLeptonEta",";LeptonEta;Number of Events",100,-5.,5.0);
  AddTH1(fLeptonPtMax        ,"hLeptonPtMax",";Lepton P_t Max;Number of Events",150,0.,150.);
  AddTH1(fLeptonPtMin        ,"hLeptonPtMin",";Lepton P_t Min;Number of Events",150,0.,150.);
  AddTH1(fMetPtHist          ,"hMetPtHist",";Met;Number of Events",150,0.,300.);  
  AddTH1(fMetPhiHist         ,"hMetPhiHist",";#phi;Number of Events",28,-3.5,3.5);
  AddTH1(fUncorrMetPtHist    ,"hUncorrMetPtHist",";Met;Number of Events",150,0.,300.);  
  AddTH1(fUncorrMetPhiHist   ,"hUncorrMetPhiHist",";#phi;Number of Events",28,-3.5,3.5);
  AddTH1(fDeltaPhiLeptons    ,"hDeltaPhiLeptons",";#Delta#phi_{ll};Number of Events",90,0,180);
  AddTH1(fDeltaEtaLeptons    ,"hDeltaEtaLeptons",";#Delta#eta_{ll};Number of Events",100,-5.,5.0);
  AddTH1(fDileptonMass       ,"hDileptonMass",";Mass_{ll};Number of Events",150,0.,300.);
 
  //***********************************************************************************************
  // N-1 Histograms
  //***********************************************************************************************
  //All events
  AddTH1(fLeptonPtMax_NMinusOne            ,"hLeptonPtMax_NMinusOne",
                                            ";Lepton P_t Max;Number of Events",150,0.,150.);
  AddTH1(fLeptonPtMin_NMinusOne            ,"hLeptonPtMin_NMinusOne",
                                            ";Lepton P_t Min;Number of Events",150,0.,150.);
  AddTH1(fMetPtHist_NMinusOne              ,"hMetPtHist_NMinusOne",
                                            ";Met;Number of Events",150,0.,300.);  
  AddTH1(fMetPhiHist_NMinusOne             ,"hMetPhiHist_NMinusOne",
                                            ";#phi;Number of Events",28,-3.5,3.5);
  AddTH1(fMETdeltaPhilEtHist_NMinusOne     ,"hMETdeltaPhilEtHist_NMinusOne",
                                            ";METdeltaPhilEtHist;Number of Events",150,0.,300.);
  AddTH1(fNCentralJets_NMinusOne           ,"hNCentralJets_NMinusOne",
                                            ";Number of Central Jets;Number of Events",6,-0.5,5.5);
  AddTH1(fNSoftMuonsHist_NMinusOne        ,"hNSoftMuonsHist_NMinusOne",
                                            ";Number of Dirty Muons; Number of Events",6,-0.5,5.5);
  AddTH1(fDeltaPhiLeptons_NMinusOne        ,"hDeltaPhiLeptons_NMinusOne",
                                            ";#Delta#phi_{ll};Number of Events",90,0,180);
  AddTH1(fDeltaEtaLeptons_NMinusOne        ,"hDeltaEtaLeptons_NMinusOne",
                                            ";#Delta#eta_{ll};Number of Events",100,-5.0,5.0);
  AddTH1(fDileptonMass_NMinusOne           ,"hDileptonMass_NMinusOne",
                                            ";Mass_{ll};Number of Events",150,0.,300.);
  AddTH1(fMinDeltaPhiLeptonMet_NMinusOne   ,"hMinDeltaPhiLeptonMet_NMinusOne", 
                                            ";Min #Delta#phi_{l,Met};Number of Events",90,0.,180);

  //***********************************************************************************************
  // After all cuts Histograms
  //***********************************************************************************************
  AddTH1(fMinDeltaPhiLeptonMet_afterCuts    ,"hMinDeltaPhiLeptonMet_afterCuts", 
                                             ";Min #Delta#phi_{l,Met};Number of Events",90,0.,180);
  AddTH1(fMtLepton1_afterCuts               ,"hMtLepton1_afterCuts",
                                             ";M_t (Lepton1,Met);Number of Events",100,0.,200.);
  AddTH1(fMtLepton2_afterCuts               ,"hMtLepton2_afterCuts",
                                             ";M_t (Lepton2,Met);Number of Events",100,0.,200.);
  AddTH1(fMtHiggs_afterCuts                 ,"hMtHiggs_afterCuts",
                                             ";M_t (l1+l2+Met);Number of Events",150,0.,300.);
  AddTH1(fLeptonPtPlusMet_afterCuts         ,"hLeptonPtPlusMet_afterCuts",
                                             ";LeptonPtPlusMet;Number of Events",150,0., 300.);

}

//--------------------------------------------------------------------------------------------------
void HwwExampleAnalysisMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and  
  LoadBranch(fMuonBranchName);
  LoadBranch(fPFCandidatesName);
  LoadBranch(fPFJetName0);

  //Obtain all the good objects from the event cleaning module
  fVertices = GetObjThisEvt<VertexOArr>(fVertexName);
  ObjArray<Muon> *CleanMuons = dynamic_cast<ObjArray<Muon>* >(FindObjThisEvt(ModNames::gkCleanMuonsName));
  ObjArray<Electron> *CleanElectrons = dynamic_cast<ObjArray<Electron>* >(FindObjThisEvt(ModNames::gkCleanElectronsName));
  ParticleOArr *CleanLeptons = dynamic_cast<mithep::ParticleOArr*>
    (FindObjThisEvt(ModNames::gkMergedLeptonsName));
  ObjArray<Jet> *CleanJetsNoPtCut = dynamic_cast<ObjArray<Jet>* >
    (FindObjThisEvt(fCleanJetsNoPtCutName.Data()));
  TParameter<Double_t> *NNLOWeight = GetObjThisEvt<TParameter<Double_t> >("NNLOWeight");

  MetCol *met = dynamic_cast<ObjArray<Met>* >(FindObjThisEvt(fMetName));
  const Met *stdMet = 0;
  if (met) {
    stdMet = met->At(0);
  } else {
    std::cout << "Error: Met Collection " << fMetName << " could not be loaded." << std::endl;
    return;
  }

  //***********************************************************************************************
  //Kinematic PreSelection
  //***********************************************************************************************
  // At least two leptons in the event
  if (CleanLeptons->GetEntries() < 2) return;
  // Pt1 > 20 && Pt2 > 10
  if(CleanLeptons->At(0)->Pt() <= 20 || CleanLeptons->At(1)->Pt() <= 10) return;
  // opposite charge leptons
  if(CleanLeptons->At(0)->Charge() * CleanLeptons->At(1)->Charge() > 0) return;
    
  CompositeParticle *dilepton = new CompositeParticle();
  dilepton->AddDaughter(CleanLeptons->At(0));
  dilepton->AddDaughter(CleanLeptons->At(1));
   
  //***********************************************************************************************
  //Get Dirty Muons: Non-isolated Muons (exclude the clean muons)
  //***********************************************************************************************
  ObjArray<Muon> *SoftMuons = new ObjArray<Muon>;
  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i);
    if(!MuonTools::PassSoftMuonCut(mu, fVertices, 0.2)) continue;
    
    bool isCleanMuon = kFALSE;
    for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
      if(fMuons->At(i) == CleanMuons->At(j) &&
  	 CleanMuons->At(j)->Pt() > 10) isCleanMuon = kTRUE;
    }
    if(isCleanMuon == kFALSE) SoftMuons->Add(mu);
  }

  //***********************************************************************************************
  //|Z_vert-Z_l| maximum
  //***********************************************************************************************
  std::vector<double> leptonsDz;
  double zDiffMax = 0.0;
  if(fVertices->GetEntries() > 0) {
    for (UInt_t j=0; j<CleanMuons->GetEntries(); j++) {
      double pDz = CleanMuons->At(j)->BestTrk()->DzCorrected(*fVertices->At(0));
      leptonsDz.push_back(pDz);
    }
    for (UInt_t j=0; j<CleanElectrons->GetEntries(); j++) {   
      double pDz = CleanElectrons->At(j)->GsfTrk()->DzCorrected(*fVertices->At(0));
      leptonsDz.push_back(pDz);
    }
    for(UInt_t t=0; t<leptonsDz.size(); t++) {
      for(UInt_t i=t+1; i<leptonsDz.size(); i++) {
        if(TMath::Abs(leptonsDz[t]-leptonsDz[i]) > zDiffMax) zDiffMax = TMath::Abs(leptonsDz[t]-leptonsDz[i]);
      }
    }
    leptonsDz.clear();
  }

  //***********************************************************************************************
  //Define Event Variables
  //***********************************************************************************************
  //delta phi between the 2 leptons in degrees
  double deltaPhiLeptons = MathUtils::DeltaPhi(CleanLeptons->At(0)->Phi(), 
                                               CleanLeptons->At(1)->Phi())* 180.0 / TMath::Pi();

  double deltaEtaLeptons = CleanLeptons->At(0)->Eta() - CleanLeptons->At(1)->Eta();

  double deltaPhiDileptonMet = MathUtils::DeltaPhi(stdMet->Phi(), 
                                                   dilepton->Phi())*180.0 / TMath::Pi();

  double mtHiggs = TMath::Sqrt(2.0*dilepton->Pt() * stdMet->Pt()*
			       (1.0 - cos(deltaPhiDileptonMet * TMath::Pi() / 180.0)));

  //angle between MET and closest lepton
  double deltaPhiMetLepton[2] = {MathUtils::DeltaPhi(stdMet->Phi(), CleanLeptons->At(0)->Phi()),
                                 MathUtils::DeltaPhi(stdMet->Phi(), CleanLeptons->At(1)->Phi())};
  
  double mTW[2] = {TMath::Sqrt(2.0*CleanLeptons->At(0)->Pt()*stdMet->Pt()*
                               (1.0 - cos(deltaPhiMetLepton[0]))),
		   TMath::Sqrt(2.0*CleanLeptons->At(1)->Pt()*stdMet->Pt()*
                               (1.0 - cos(deltaPhiMetLepton[1])))};

  double minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
    deltaPhiMetLepton[0]:deltaPhiMetLepton[1];

  MetTools metTools(CleanMuons, CleanElectrons, fPFCandidates, fVertices->At(0), 0.1, 8.0, 5.0);
  double pMET[2] = {metTools.GetProjectedMet(CleanLeptons,stdMet),
  		    metTools.GetProjectedTrackMet(CleanLeptons)};

  double METdeltaPhilEt = TMath::Min(pMET[0],pMET[1]);

  //count the number of central Jets for vetoing and b-tagging
  vector<Jet*> sortedJetsAll;
  vector<Jet*> sortedJets;
  vector<Jet*> sortedJetsLowPt;
  for(UInt_t i=0; i<CleanJetsNoPtCut->GetEntries(); i++){
    if(CleanJetsNoPtCut->At(i)->RawMom().Pt() <= 7) continue;
    Jet* jet_a = new Jet(CleanJetsNoPtCut->At(i)->Px(),
   			 CleanJetsNoPtCut->At(i)->Py(),
   			 CleanJetsNoPtCut->At(i)->Pz(),
   			 CleanJetsNoPtCut->At(i)->E() );
    jet_a->SetMatchedMCFlavor(CleanJetsNoPtCut->At(i)->MatchedMCFlavor());
    jet_a->SetCombinedSecondaryVertexBJetTagsDisc(CleanJetsNoPtCut->At(i)->CombinedSecondaryVertexBJetTagsDisc());
    jet_a->SetCombinedSecondaryVertexMVABJetTagsDisc(CleanJetsNoPtCut->At(i)->CombinedSecondaryVertexMVABJetTagsDisc());
    jet_a->SetJetProbabilityBJetTagsDisc(CleanJetsNoPtCut->At(i)->JetProbabilityBJetTagsDisc());
    jet_a->SetJetBProbabilityBJetTagsDisc(CleanJetsNoPtCut->At(i)->JetBProbabilityBJetTagsDisc());
    jet_a->SetTrackCountingHighEffBJetTagsDisc(CleanJetsNoPtCut->At(i)->TrackCountingHighEffBJetTagsDisc());
    jet_a->SetTrackCountingHighPurBJetTagsDisc(CleanJetsNoPtCut->At(i)->TrackCountingHighPurBJetTagsDisc());
    jet_a->SetSimpleSecondaryVertexBJetTagsDisc(CleanJetsNoPtCut->At(i)->SimpleSecondaryVertexBJetTagsDisc());
    jet_a->SetSimpleSecondaryVertexHighEffBJetTagsDisc(CleanJetsNoPtCut->At(i)->SimpleSecondaryVertexHighEffBJetTagsDisc());
    jet_a->SetSimpleSecondaryVertexHighPurBJetTagsDisc(CleanJetsNoPtCut->At(i)->SimpleSecondaryVertexHighPurBJetTagsDisc());
    sortedJetsAll.push_back(jet_a);
  }

  for(UInt_t i=0; i<CleanJetsNoPtCut->GetEntries(); i++){
    if(TMath::Abs(CleanJetsNoPtCut->At(i)->Eta()) < 5.0 &&
       CleanJetsNoPtCut->At(i)->Pt() > 30.0){
      Jet* jet_b = new Jet(CleanJetsNoPtCut->At(i)->Px(),
     			   CleanJetsNoPtCut->At(i)->Py(),
   			   CleanJetsNoPtCut->At(i)->Pz(),
   			   CleanJetsNoPtCut->At(i)->E() );
      sortedJets.push_back(jet_b);
    }
  }

  for(UInt_t i=0; i<sortedJetsAll.size(); i++){
    bool overlap = kFALSE;
    for(UInt_t j=0; j<sortedJets.size(); j++){
      if(sortedJetsAll[i]->Pt() == sortedJets[j]->Pt() ||
        (sortedJetsAll[i]->CombinedSecondaryVertexBJetTagsDisc() == sortedJets[j]->CombinedSecondaryVertexBJetTagsDisc() &&
	 sortedJetsAll[i]->JetBProbabilityBJetTagsDisc()	 == sortedJets[j]->JetBProbabilityBJetTagsDisc() &&
	 sortedJetsAll[i]->TrackCountingHighPurBJetTagsDisc()	 == sortedJets[j]->TrackCountingHighPurBJetTagsDisc())
        ) {
        sortedJets[j]->SetMatchedMCFlavor(sortedJetsAll[i]->MatchedMCFlavor());
        sortedJets[j]->SetCombinedSecondaryVertexBJetTagsDisc(sortedJetsAll[i]->CombinedSecondaryVertexBJetTagsDisc());
        sortedJets[j]->SetCombinedSecondaryVertexMVABJetTagsDisc(sortedJetsAll[i]->CombinedSecondaryVertexMVABJetTagsDisc());
        sortedJets[j]->SetJetProbabilityBJetTagsDisc(sortedJetsAll[i]->JetProbabilityBJetTagsDisc());
        sortedJets[j]->SetJetBProbabilityBJetTagsDisc(sortedJetsAll[i]->JetBProbabilityBJetTagsDisc());
        sortedJets[j]->SetTrackCountingHighEffBJetTagsDisc(sortedJetsAll[i]->TrackCountingHighEffBJetTagsDisc());
        sortedJets[j]->SetTrackCountingHighPurBJetTagsDisc(sortedJetsAll[i]->TrackCountingHighPurBJetTagsDisc());
        sortedJets[j]->SetSimpleSecondaryVertexBJetTagsDisc(sortedJetsAll[i]->SimpleSecondaryVertexBJetTagsDisc());
        sortedJets[j]->SetSimpleSecondaryVertexHighEffBJetTagsDisc(sortedJetsAll[i]->SimpleSecondaryVertexHighEffBJetTagsDisc());
        sortedJets[j]->SetSimpleSecondaryVertexHighPurBJetTagsDisc(sortedJetsAll[i]->SimpleSecondaryVertexHighPurBJetTagsDisc());        
  	overlap = kTRUE;
        break;
      }
    }
    if(overlap == kFALSE){
      sortedJetsLowPt.push_back(sortedJetsAll[i]);
    }
  }
  double maxBtag = -99999.;
  for(UInt_t i=0; i<sortedJetsLowPt.size(); i++){
    if(sortedJetsLowPt[i]->TrackCountingHighEffBJetTagsDisc() > maxBtag){
      double dZAverageJetPt = 0.0;
      double sumJetPt = 0.0;
      double jetPt = 0.0;
      for(UInt_t iPF=0; iPF<fPFJet0->GetEntries(); iPF++){						  	      
        const PFJet *jet = fPFJet0->At(iPF);									
        if(MathUtils::DeltaR(jet->Mom(),sortedJetsLowPt[i]->Mom()) < 0.01){
          jetPt = jet->Pt();
          for (UInt_t npf=0; npf<jet->NPFCands();npf++) {
            const PFCandidate *pf = jet->PFCand(npf);
            if(pf->BestTrk()) {
              dZAverageJetPt = dZAverageJetPt + pf->Pt()*pf->Pt()*pf->BestTrk()->DzCorrected(*fVertices->At(0));
              sumJetPt = sumJetPt + pf->Pt()*pf->Pt();
            }
          }
          if(sumJetPt > 0) dZAverageJetPt = TMath::Abs(dZAverageJetPt)/sumJetPt;
          break;
        }
      } // loop over PF jets
      if(dZAverageJetPt < 2.0 && jetPt > 10){
        maxBtag  = sortedJetsLowPt[i]->TrackCountingHighEffBJetTagsDisc();
      }
    }
  }

  //Lepton Type
  int finalstateType = -1;
  if (CleanLeptons->At(0)->ObjType() == kMuon && CleanLeptons->At(1)->ObjType() == kMuon ){ // mumu
    finalstateType = 10;
  } else if(CleanLeptons->At(0)->ObjType() == kElectron && CleanLeptons->At(1)->ObjType() == kElectron ){ // ee
    finalstateType = 11;
  } else if(CleanLeptons->At(0)->ObjType() == kElectron && CleanLeptons->At(1)->ObjType() == kMuon) {
    finalstateType = 12;
  } else if(CleanLeptons->At(1)->ObjType() == kElectron && CleanLeptons->At(0)->ObjType() == kMuon) {
    finalstateType = 13;
  } else {
    std::cerr << "Error: finalstate lepton type not supported" << std::endl;
  }
                        
  double deltaPhiLLJet = 0.0;
  if(sortedJetsAll.size() > 0 && sortedJetsAll[0]->Pt() > 15.0 && (finalstateType == 10 || finalstateType == 11)){
    deltaPhiLLJet = MathUtils::DeltaPhi(dilepton->Phi(), sortedJetsAll[0]->Phi())*180.0/TMath::Pi();
  }

  //*********************************************************************************************
  //Define Cuts
  //*********************************************************************************************
  const int nCuts = 16;
  bool passCut[nCuts] = {false, false, false, false, false,
                         false, false, false, false, false,
			 false, false, false, false, false,
			 false};
  
  Bool_t PreselPtCut = kTRUE;
  if(CleanLeptons->At(0)->Pt() <= 20) PreselPtCut = kFALSE;
  if(CleanLeptons->At(1)->Pt() <= 10) PreselPtCut = kFALSE;
  //if(CleanLeptons->At(1)->ObjType() == kElectron && CleanLeptons->At(1)->Pt() <= 15) PreselPtCut = kFALSE;
  if(PreselPtCut == kTRUE)              passCut[0] = true;
  
  if(stdMet->Pt()    > 20.0)            passCut[1] = true;
  
  if(dilepton->Mass() > 12.0)           passCut[2] = true;
  
  if(sortedJets.size() < 1)             passCut[5] = true;

  if(deltaPhiLLJet < 165.0)		passCut[6] = true;

  if(SoftMuons->GetEntries() == 0)      passCut[7] = true;

  if(CleanLeptons->GetEntries() == 2)   passCut[8] = true;

  if(maxBtag < 2.1)                     passCut[9] = true;

  if (finalstateType == 10 || finalstateType == 11){ // mumu/ee
    if(fabs(dilepton->Mass()-91.1876)   > 15.0)   passCut[3] = true;
    if(METdeltaPhilEt > 37.0 + fVertices->GetEntries()/2.0) passCut[4] = true;
  }
  else { // mue/emu
    passCut[3] = true;
    if(METdeltaPhilEt > 20) passCut[4] = true;
  }

  if(CleanLeptons->At(0)->Pt() > 30)    passCut[10] = true;

  if(CleanLeptons->At(1)->Pt() > 25)    passCut[11] = true;

  if(dilepton->Mass() < 50)             passCut[12] = true;

  if(mtHiggs > 90.0 && mtHiggs < 160.0) passCut[13] = true;

  if(deltaPhiLeptons < 60.0)            passCut[14] = true;

  if(dilepton->Pt() > 45.0)             passCut[15] = true;

  //*********************************************************************************************
  //Make Selection Histograms. Number of events passing each level of cut
  //*********************************************************************************************  
  bool passAllCuts = true;
  for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
  if(passAllCuts) fNEventsSelected++;

  //Cut Selection Histograms
  fHWWSelection->Fill(-1,NNLOWeight->GetVal());
  if (finalstateType == 10 )
    fHWWToMuMuSelection->Fill(-1,NNLOWeight->GetVal());
  else if(finalstateType == 11 )
    fHWWToEESelection->Fill(-1,NNLOWeight->GetVal());
  else if(finalstateType == 12 )
    fHWWToEMuSelection->Fill(-1,NNLOWeight->GetVal());
  else if(finalstateType == 13 )
    fHWWToMuESelection->Fill(-1,NNLOWeight->GetVal());

  for (int k=0;k<nCuts;k++) {
    bool pass = true;
    bool passPreviousCut = true;
    for (int p=0;p<=k;p++) {
      pass = (pass && passCut[p]);
      if (p<k)
        passPreviousCut = (passPreviousCut&& passCut[p]);
    }
    
    if (pass) {
      fHWWSelection->Fill(k,NNLOWeight->GetVal());
      if (finalstateType == 10 )
        fHWWToMuMuSelection->Fill(k,NNLOWeight->GetVal());
      else if(finalstateType == 11)
        fHWWToEESelection->Fill(k,NNLOWeight->GetVal());
      else if(finalstateType == 12)
        fHWWToEMuSelection->Fill(k,NNLOWeight->GetVal());
      else if(finalstateType == 13)
        fHWWToMuESelection->Fill(k,NNLOWeight->GetVal());
    }
  }
  
  //*****************************************************************************************
  //Make Preselection Histograms  
  //*****************************************************************************************
  fLeptonEta->Fill(CleanLeptons->At(0)->Eta(),NNLOWeight->GetVal()); 
  fLeptonEta->Fill(CleanLeptons->At(1)->Eta(),NNLOWeight->GetVal());
  fLeptonPtMax->Fill(CleanLeptons->At(0)->Pt(),NNLOWeight->GetVal());
  fLeptonPtMin->Fill(CleanLeptons->At(1)->Pt(),NNLOWeight->GetVal());
  fMetPtHist->Fill(stdMet->Pt(),NNLOWeight->GetVal());                             
  fMetPhiHist->Fill(stdMet->Phi(),NNLOWeight->GetVal());                            
  fDeltaPhiLeptons->Fill(deltaPhiLeptons,NNLOWeight->GetVal());
  fDeltaEtaLeptons->Fill(deltaEtaLeptons,NNLOWeight->GetVal());
  fDileptonMass->Fill(dilepton->Mass(),NNLOWeight->GetVal());    

  //*********************************************************************************************
  // N-1 Histograms
  //*********************************************************************************************
  bool pass;;
  
  //N Jet Veto  
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 5) pass = (pass && passCut[k]);      
  }
  if (pass) {
    fNCentralJets_NMinusOne->Fill(sortedJets.size(),NNLOWeight->GetVal());
  }     
  
  // Final Met Cut
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 4) pass = (pass && passCut[k]);      
  }
  if (pass) {
    fMetPtHist_NMinusOne->Fill(stdMet->Pt(),NNLOWeight->GetVal());  
  }

  // dilepton mass
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 2 && k !=  3) pass = (pass && passCut[k]);    
  }
  if (pass) {
    fDileptonMass_NMinusOne->Fill(dilepton->Mass(),NNLOWeight->GetVal());
  }
  
  // Lepton Pt Max, Lepton Pt Min, DeltaPhiLeptons
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 0)
      pass = (pass && passCut[k]);
  }
  if (pass) {
    fLeptonPtMax_NMinusOne->Fill(CleanLeptons->At(0)->Pt(),NNLOWeight->GetVal());
    fLeptonPtMin_NMinusOne->Fill(CleanLeptons->At(1)->Pt(),NNLOWeight->GetVal());
    fDeltaPhiLeptons_NMinusOne->Fill(deltaPhiLeptons,NNLOWeight->GetVal()); 
  }
  
  // NSoftMuons
  pass = true;
  for (int k=0;k<nCuts;k++) {
    if (k != 7) pass = (pass && passCut[k]);    
  }
  if (pass) {
    fNSoftMuonsHist_NMinusOne->Fill(SoftMuons->GetEntries(),NNLOWeight->GetVal());
  }

  //*********************************************************************************************
  //Plots after all Cuts
  //*********************************************************************************************
  if (passAllCuts) {
    fMinDeltaPhiLeptonMet_afterCuts->Fill(minDeltaPhiMetLepton,NNLOWeight->GetVal());
    fMtLepton1_afterCuts->Fill(mTW[0],NNLOWeight->GetVal());
    fMtLepton2_afterCuts->Fill(mTW[1],NNLOWeight->GetVal());
    fMtHiggs_afterCuts->Fill(mtHiggs,NNLOWeight->GetVal());
    fLeptonPtPlusMet_afterCuts->Fill(CleanLeptons->At(0)->Pt()+CleanLeptons->At(1)->Pt()+stdMet->Pt(),NNLOWeight->GetVal());
  }
  
  delete dilepton;
  delete SoftMuons;
  for(UInt_t i=0; i<sortedJets.size();      i++) delete sortedJets[i];
  for(UInt_t i=0; i<sortedJetsAll.size();   i++) delete sortedJetsAll[i];
  return;
}

//--------------------------------------------------------------------------------------------------
void HwwExampleAnalysisMod::SlaveTerminate()
{
  
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.
  std::cout << "selected events on HwwExampleAnalysisMod: " << fNEventsSelected << std::endl;

} 
//--------------------------------------------------------------------------------------------------
void HwwExampleAnalysisMod::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.

}
