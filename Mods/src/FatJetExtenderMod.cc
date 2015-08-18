#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Mods/interface/FatJetExtenderMod.h"
#include "MitAna/DataTree/interface/PFJetCol.h"

#include "MitCommon/DataFormats/interface/Vect4M.h"
#include "MitCommon/DataFormats/interface/Vect3.h"
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/Utils/interface/Utils.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "QjetsPlugin.h"
#include "Qjets.h"

#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "MitAna/PhysicsUtils/interface/HEPTopTagger.h"

using namespace mithep;

ClassImp(mithep::FatJetExtenderMod)

//--------------------------------------------------------------------------------------------------
FatJetExtenderMod::FatJetExtenderMod(const char *name, const char *title) :
  BaseMod (name,title),
  fIsData (kFALSE),
  fQGTaggingActive (kTRUE),
  fQGTaggerCHS (kFALSE),
  fPublishOutput (kTRUE),
  fFatJetsName (""),
  fFatJetsFromBranch (kFALSE),
  fFatJets (0),
  fPFCandidatesName (Names::gkPFCandidatesBrn),
  fPFCandidates (0),
  fPileUpDenName(Names::gkPileupEnergyDensityBrn),
  fPileUpDenFromBranch(kTRUE),
  fPileUpDen(0),
  fVertexesName ("GoodVertexes"),
  fVertexesFromBranch(kFALSE),
  fVertexes(0),
  fXlFatJetsName ("XlFatJets"),
  fSoftDropZCut (0.1),
  fSoftDropR0 (.8),
  fPruneZCut (0.1),
  fPruneDistCut (0.5),
  fFilterN (3),
  fFilterRad (0.2),
  fTrimRad (0.2),
  fTrimPtFrac (0.05),
  fConeSize (0.8),
  fDeconstruct(0),
  fProcessNJets (4),
  fDoShowerDeconstruction(kFALSE),
  fBeVerbose(kFALSE),
  fDoECF(kFALSE),
  fNMaxMicrojets(5)
{
  // Constructor.

    // setup subjet branch names
    fXlSubJetsName[0] = fXlFatJetsName+"_SoftDropSubjets";
    fXlSubJetsName[1] = fXlFatJetsName+"_PrunedSubjets";
    fXlSubJetsName[2] = fXlFatJetsName+"_TrimmedSubjets";
    fXlSubJetsName[3] = fXlFatJetsName+"_CMSTTSubjets";
    fXlSubJetsName[4] = fXlFatJetsName+"_HEPTTSubjets";
    fXlSubJetsName[5] = fXlFatJetsName+"_NjettinessSubjets";


}

FatJetExtenderMod::~FatJetExtenderMod()
{
  // Destructor
  if (fXlSubJets){
    for(int i=0; i<XlSubJet::nSubJetTypes; ++i) {
      if (fSubJetFlags & (1<<i))
        delete fXlSubJets[i];
    }
  }

  if (fXlFatJets)
    delete fXlFatJets;

  if (fPruner)
    delete fPruner;
  if (fFilterer)
    delete fFilterer;
  if (fTrimmer)
    delete fTrimmer ;
  if (fCMSTopTagger)
    delete fCMSTopTagger;
  if (fCAJetDef)
    delete fCAJetDef;

  if (fActiveArea)
    delete fActiveArea;
  if (fAreaDefinition)
    delete fAreaDefinition;

  if (fQGTagger)
    delete fQGTagger;

  if (fDeconstruct)
    delete fDeconstruct;
  if (fParam)
    delete fParam;
  if (fSignal)
    delete fSignal;
  if (fBackground)
    delete fBackground;
  if (fISR)
    delete fISR;

  if (fStopwatch)
    delete fStopwatch;
}

//--------------------------------------------------------------------------------------------------
void FatJetExtenderMod::Process()
{
  if (fBeVerbose) fStopwatch->Start(kTRUE);
  // make sure the out collections are empty before starting
  fXlFatJets->Delete();
  for(unsigned int i=0; i<XlSubJet::nSubJetTypes; ++i) {
    if (fSubJetFlags & (1<<i))
      fXlSubJets[i]->Delete();
  }
  fFatJets = GetObject<JetCol>(fFatJetsName);
  fPFCandidates = GetObject<PFCandidateCol>(fPFCandidatesName);
  if (fQGTaggingActive){
    fPileUpDen = GetObject<PileupEnergyDensityCol>(fPileUpDenName);
    fVertexes = GetObject<VertexCol>(fVertexesName);
  }
if (fDebugFlag == 0) return;
  // Loop over PFCandidates and unmark them : necessary for skimming
  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i)
    fPFCandidates->At(i)->UnmarkMe();

  // Setup pileup density for QG computation
  if (fQGTaggingActive)
    fQGTagger->SetRhoIso(fPileUpDen->At(0)->RhoRandomLowEta());


  // Loop over jets
  for (UInt_t i=0; i<fFatJets->GetEntries(); ++i) {

    // consider only the first fProcessNJets jets
    if (i >= fProcessNJets)
      break;

    const FatJet *jet = dynamic_cast<const FatJet*>(fFatJets->At(i));
    if (! jet) {
      printf(" FatJetExtenderMod::Process() - ERROR - jets provided are not FatJets.");
      break;
    }

    // mark jet (and consequently its consituents) for further use in skim
    jet->Mark();
if (fDebugFlag == 1) continue;
    if (fBeVerbose) {
			fprintf(stderr,"Finished setup in %f seconds\n",fStopwatch->RealTime()); fStopwatch->Start();
		}
    FillXlFatJet(jet);

  }

  return;
}

//--------------------------------------------------------------------------------------------------
void FatJetExtenderMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis.
  // Initialize area caculation (done with ghost particles)

  // Create the new output collection
  fXlFatJets = new XlFatJetArr(16,fXlFatJetsName);
  for(unsigned int i = XlSubJet::kSoftDrop; i<XlSubJet::nSubJetTypes; ++i) {
    // only allocate memory for the subjets that are turned on
    if (fSubJetFlags & (1<<i))
      fXlSubJets[i] = new XlSubJetArr(16,fXlSubJetsName[i]);
  }
  // Publish collection for further usage in the analysis
  if (fPublishOutput) {
    PublishObj(fXlFatJets);
    for(unsigned int i = XlSubJet::kSoftDrop; i<XlSubJet::nSubJetTypes; ++i) {
      if (fSubJetFlags & (1<<i))
        PublishObj(fXlSubJets[i]);
    }
  }

  // Prepare pruner
  fPruner = new fastjet::Pruner(fastjet::cambridge_algorithm,fPruneZCut,fPruneDistCut);
  // Prepare filterer
  fFilterer = new fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm,fFilterRad),
                                  fastjet::SelectorNHardest(fFilterN));
  // Prepare trimmer
  fTrimmer = new fastjet::Filter(fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm,fTrimRad),
                                 fastjet::SelectorPtFractionMin(fTrimPtFrac)));

  // CA constructor (fConeSize = 0.6 for antiKt) - reproducing paper 1: http://arxiv.org/abs/1011.2268
  fCAJetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, fConeSize);
  fCMSTopTagger = new fastjet::CMSTopTagger();

  // Initialize area caculation (done with ghost particles)
  int activeAreaRepeats = 1;
  double ghostArea = 0.01;
  double ghostEtaMax = 7.0;
  fActiveArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
  fAreaDefinition = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*fActiveArea);

  // Initialize QGTagger class
  fQGTagger = new QGTagger(fQGTaggerCHS);

  // set up shower deconstruction stuff
  TString inputCard = Utils::GetEnv("MIT_DATA");
  inputCard += TString::Format("/SDAlgorithm/input_card_%i.dat",int(fConeSize*10));
  fParam = new AnalysisParameters(inputCard.Data());
  fSignal = new Deconstruction::TopGluonModel(*fParam);
  fBackground = new Deconstruction::BackgroundModel(*fParam);
  fISR = new Deconstruction::ISRModel(*fParam);
  fDeconstruct = new Deconstruction::Deconstruct(*fParam, *fSignal, *fBackground, *fISR);

  fStopwatch = new TStopwatch();

  return;
}

//--------------------------------------------------------------------------------------------------
void FatJetExtenderMod::SlaveTerminate()
{
}

//--------------------------------------------------------------------------------------------------
void FatJetExtenderMod::FillXlFatJet(const FatJet *fatJet)
{
  // Prepare and store in an array a new FatJet
  XlFatJet *xlFatJet = fXlFatJets->AddNew();
  // XlFatJet *xlFatJet = fXlFatJets->Allocate();
  // new (xlFatJet) XlFatJet(*fatJet);
if (fDebugFlag == 2) return;
  // Prepare and store QG tagging info
  float qgValue = -1.;
  if (fQGTaggingActive) {
    fQGTagger->CalculateVariables(fatJet, fVertexes);
    qgValue = fQGTagger->QGValue();
  }
  xlFatJet->SetQGTag(qgValue);
if (fDebugFlag == 3) return;
    if (fBeVerbose) {
			fprintf(stderr,"Finished QG-tagging in %f seconds\n",fStopwatch->RealTime()); fStopwatch->Start();
		}

  std::vector<fastjet::PseudoJet> fjParts;
  // Push all particle flow candidates of the input PFjet into fastjet particle collection
  for (UInt_t j=0; j<fatJet->NPFCands(); ++j) {
    const PFCandidate *pfCand = fatJet->PFCand(j);
    fjParts.push_back(fastjet::PseudoJet(pfCand->Px(),pfCand->Py(),pfCand->Pz(),pfCand->E()));
    fjParts.back().set_user_index(j);
  }

  // Setup the cluster for fastjet
  fastjet::ClusterSequenceArea *fjClustering =
    new fastjet::ClusterSequenceArea(fjParts,*fCAJetDef,*fAreaDefinition);
  // ---- Fastjet is ready ----


  // Produce a new set of jets based on the fastjet particle collection and the defined clustering
  // Cut off fat jets with pt < 10 GeV and consider only the hardest jet of the output collection
  std::vector<fastjet::PseudoJet> fjOutJets = sorted_by_pt(fjClustering->inclusive_jets(10.));

  // Check that the output collection size is non-null, otherwise nothing to be done further
  if (fjOutJets.size() < 1) {
    printf(" FatJetExtenderMod::FillXlFatJet() - WARNING - input FatJet produces null reclustering output!\n");
    if (fjOutJets.size()>0)
      fjClustering->delete_self_when_unused();
    delete fjClustering;

    return;
  }
if (fDebugFlag == 4) {
  if (fjOutJets.size()>0) fjClustering->delete_self_when_unused();
  delete fjClustering;
  return;
}

  fastjet::PseudoJet fjJet = fjOutJets[0];

    if (fBeVerbose) {
			fprintf(stderr,"Finished prepping fastjet in %f seconds\n",fStopwatch->RealTime()); fStopwatch->Start();
		}

  fastjet::contrib::Njettiness::AxesMode axisMode = fastjet::contrib::Njettiness::onepass_kt_axes;
  double beta = 1.0;
  double RNsub = fConeSize;
  double Rcutoff = 10000.;
  fastjet::contrib::Nsubjettiness  nSub1(1,axisMode,beta,RNsub,Rcutoff);
  fastjet::contrib::Nsubjettiness  nSub2(2,axisMode,beta,RNsub,Rcutoff);
  fastjet::contrib::Nsubjettiness  nSub3(3,axisMode,beta,RNsub,Rcutoff);
  fastjet::contrib::Nsubjettiness  nSub4(4,axisMode,beta,RNsub,Rcutoff);
  double tau1 = nSub1(fjJet);
  double tau2 = nSub2(fjJet);
  double tau3 = nSub3(fjJet);
  double tau4 = nSub4(fjJet);
  xlFatJet->SetTau1(tau1);
  xlFatJet->SetTau2(tau2);
  xlFatJet->SetTau3(tau3);
  xlFatJet->SetTau4(tau4);
if (fDebugFlag == 5) {
  if (fjOutJets.size()>0) fjClustering->delete_self_when_unused();
  delete fjClustering;
  return;
}

    if (fBeVerbose) {
			fprintf(stderr,"Finished njettiness calculation in %f seconds\n",fStopwatch->RealTime()); fStopwatch->Start();
		}

  if (fDoECF) {
    // Compute the energy correlation function ratios
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECR2b0  (2,0. ,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECR2b0p2(2,0.2,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECR2b0p5(2,0.5,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECR2b1  (2,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fastjet::contrib::EnergyCorrelatorDoubleRatio ECR2b2  (2,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
    double C2b0   = ECR2b0(fjJet);
    double C2b0p2 = ECR2b0p2(fjJet);
    double C2b0p5 = ECR2b0p5(fjJet);
    double C2b1   = ECR2b1(fjJet);
    double C2b2   = ECR2b2(fjJet);
    xlFatJet->SetC2b0(C2b0);
    xlFatJet->SetC2b0p2(C2b0p2);
    xlFatJet->SetC2b0p5(C2b0p5);
    xlFatJet->SetC2b1(C2b1);
    xlFatJet->SetC2b2(C2b2);

      if (fBeVerbose) {
  			fprintf(stderr,"Finished ECF calculation in %f seconds\n",fStopwatch->RealTime()); fStopwatch->Start();
  		}
  }

  // Compute Q-jets volatility
  std::vector<fastjet::PseudoJet> constits;
  GetJetConstituents(fjJet, constits, 0.01);
  double QJetVol = GetQjetVolatility(constits, 25, fCounter*25);
  fCounter++;
  constits.clear();
  xlFatJet->SetQJetVol(QJetVol);
if (fDebugFlag == 6) {
  if (fjOutJets.size()>0) fjClustering->delete_self_when_unused();
  delete fjClustering;
  return;
}

    if (fBeVerbose) {
			fprintf(stderr,"Finished Qjets in %f seconds\n",fStopwatch->RealTime()); fStopwatch->Start();
		}

  // do grooming and subjetting
  double thisJEC = xlFatJet->Pt()/xlFatJet->RawMom().Pt();
  fastjet::PseudoJet fjClusteredJets[XlSubJet::kTrimmed+1];
  fastjet::contrib::SoftDrop softDropSDb0(0.0, fSoftDropZCut, fSoftDropR0);
  fastjet::contrib::SoftDrop softDropSDb1(1.0, fSoftDropZCut, fSoftDropR0);
  fastjet::contrib::SoftDrop softDropSDb2(2.0, fSoftDropZCut, fSoftDropR0);
  fastjet::contrib::SoftDrop softDropSDbm1(-1.0, fSoftDropZCut, fSoftDropR0);
  softDropSDb0.set_tagging_mode();
  softDropSDb1.set_tagging_mode();
  softDropSDb2.set_tagging_mode();
  softDropSDbm1.set_tagging_mode();
  fjClusteredJets[XlSubJet::kSoftDrop] = softDropSDb0(fjJet);
  double MassSDb0 = fjClusteredJets[XlSubJet::kSoftDrop].m();
  double MassSDb1 = (softDropSDb1(fjJet)).m();
  double MassSDb2 = (softDropSDb2(fjJet)).m();
  double MassSDbm1 = (softDropSDbm1(fjJet)).m();
  xlFatJet->SetMassSDb0(MassSDb0*thisJEC);
  xlFatJet->SetMassSDb1(MassSDb1*thisJEC);
  xlFatJet->SetMassSDb2(MassSDb2*thisJEC);
  xlFatJet->SetMassSDbm1(MassSDbm1*thisJEC);

  fjClusteredJets[XlSubJet::kPruned] = (*fPruner)(fjJet);
  fjClusteredJets[XlSubJet::kTrimmed] = (*fTrimmer)(fjJet);
  double MassPruned = fjClusteredJets[XlSubJet::kPruned].m();
  double MassFiltered = ((*fFilterer)(fjJet)).m();
  double MassTrimmed = fjClusteredJets[XlSubJet::kTrimmed].m();
  xlFatJet->SetMassPruned(MassPruned*thisJEC);
  xlFatJet->SetMassFiltered(MassFiltered*thisJEC);
  xlFatJet->SetMassTrimmed(MassTrimmed*thisJEC);
if (fDebugFlag == 7) {
  if (fjOutJets.size()>0) fjClustering->delete_self_when_unused();
  delete fjClustering;
  return;
}

    if (fBeVerbose) {
			fprintf(stderr,"Finished grooming in %f seconds\n",fStopwatch->RealTime()); fStopwatch->Start();
		}

  // do CMS and HEP top tagging
  std::vector<fastjet::PseudoJet> lOutJets = sorted_by_pt(fjClustering->inclusive_jets(0.0));
  fastjet::PseudoJet iJet = lOutJets[0];
  HEPTopTagger hepTopJet = HEPTopTagger(*fjClustering,iJet);;
  fastjet::PseudoJet cmsTopJet = fCMSTopTagger->result(iJet);
if (fDebugFlag == 8) {
  if (fjOutJets.size()>0) fjClustering->delete_self_when_unused();
  delete fjClustering;
  return;
}

    if (fBeVerbose) {
			fprintf(stderr,"Finished CMSTT and HEPTT in %f seconds\n",fStopwatch->RealTime()); fStopwatch->Start();
		}

  // fill subjets
  Bool_t computedPullAngle = kFALSE;
  for (unsigned int iSJType = 0; iSJType!=XlSubJet::nSubJetTypes; ++iSJType) {
    if (!(fSubJetFlags & (1<<iSJType)))
      continue;
    // okay, we really do want to save this collection
    std::vector<fastjet::PseudoJet> fjSubjets;
    if (iSJType <= XlSubJet::kTrimmed) {
      // these have a common interface
      if (!fjClusteredJets[iSJType].has_constituents())
        // if subjet finding failed, skip this step
        continue;
      int nSubJets = std::min<unsigned int>(fjClusteredJets[iSJType].constituents().size(),3);
      if (nSubJets>0) {
        fjSubjets = fjClusteredJets[iSJType].associated_cluster_sequence()->exclusive_subjets(fjClusteredJets[iSJType],nSubJets);
        std::vector<fastjet::PseudoJet> fjSubJetsSorted = Sorted_by_pt_min_pt(fjSubjets,0.01);
        if (!computedPullAngle) {
          xlFatJet->SetPullAngle(GetPullAngle(fjSubJetsSorted,0.01));
          computedPullAngle = kTRUE;
        }
        FillXlSubJets(fjSubJetsSorted,xlFatJet,(ESubJetType)iSJType);
      }
    } else if (iSJType == XlSubJet::kCMSTT) {
      fjSubjets = cmsTopJet.pieces();
      FillXlSubJets(fjSubjets,xlFatJet,XlSubJet::kCMSTT);
    } else if (iSJType == XlSubJet::kHEPTT) {
      fjSubjets = hepTopJet.top_subjets();
      FillXlSubJets(fjSubjets,xlFatJet,XlSubJet::kHEPTT);
    }
  }
if (fDebugFlag == 9) {
  if (fjOutJets.size()>0) fjClustering->delete_self_when_unused();
  delete fjClustering;
  return;
}
    if (fBeVerbose) {
			fprintf(stderr,"Finished filling subjets in %f seconds\n",fStopwatch->RealTime()); fStopwatch->Start();
		}

  // take a shower
  if (fDoShowerDeconstruction) {
      // shower deconstruction
      double microconesize;
      if  (fMicrojetR0<0){
        // From Tobias:
        // 0..500   -> 0.3
        // 500..700 -> 0.2
        // 700..inf -> 0.
        if (fatJet->Pt() < 500)
         	 microconesize=0.3;
        else if (fatJet->Pt() < 700)
         	 microconesize=0.2;
        else
         	 microconesize=0.1;
      } else {
        microconesize = fMicrojetR0;
      }
      fastjet::JetDefinition reclustering(fastjet::JetAlgorithm::kt_algorithm, microconesize);
      fastjet::ClusterSequence * cs_micro = new fastjet::ClusterSequence(fjParts, reclustering);
      std::vector<fastjet::PseudoJet> microjets = fastjet::sorted_by_pt(cs_micro->inclusive_jets(10.));
      xlFatJet->SetNMicrojets(microjets.size());
      if (microjets.size()>fNMaxMicrojets)
        microjets.erase(microjets.begin()+fNMaxMicrojets,microjets.end());
      double Psignal = 0.0;
      double Pbackground = 0.0;
      try {
        double chi = fDeconstruct->deconstruct(microjets, Psignal, Pbackground);
        xlFatJet->SetChi(chi);
      } catch(Deconstruction::Exception &e) {
        std::cout << "Exception while running SD: " << e.what() << std::endl;
      }
    if (fBeVerbose) {
			fprintf(stderr,"Finished shower deconstruction in %f seconds\n",fStopwatch->RealTime()); fStopwatch->Start();
		}
    if (cs_micro->inclusive_jets(0.).size()>0)
      cs_micro->delete_self_when_unused();
    delete cs_micro;
  }

  // Store groomed 4-momenta, apply JEC
  fastjet::PseudoJet fj_tmp;
  if (fSubJetFlags & (1<<XlSubJet::kSoftDrop)){
    fj_tmp = fjClusteredJets[XlSubJet::kSoftDrop];
    xlFatJet->SetSoftDropP(GetCorrectedMomentum(fj_tmp,thisJEC));
  }
  if (fSubJetFlags & (1<<XlSubJet::kPruned)){
    fj_tmp = fjClusteredJets[XlSubJet::kPruned];
    xlFatJet->SetPrunedP(GetCorrectedMomentum(fj_tmp,thisJEC));
  }
  if (fSubJetFlags & (1<<XlSubJet::kTrimmed)){
    fj_tmp = fjClusteredJets[XlSubJet::kTrimmed];
    xlFatJet->SetTrimmedP(GetCorrectedMomentum(fj_tmp,thisJEC));
  }
  \
  // Store the color pull
  xlFatJet->SetPull(GetPull(fjJet,0.01).Mod());

  // Trim the output collections
  fXlFatJets->Trim();

  // Memory cleanup
  if (fjOutJets.size() > 0)
    fjClustering->delete_self_when_unused();
  delete fjClustering;
    if (fBeVerbose) {
			fprintf(stderr,"Finished filling and cleanup in %f seconds\n",fStopwatch->RealTime());
    }

  return;
}

//--------------------------------------------------------------------------------------------------

inline Vect4M FatJetExtenderMod::GetCorrectedMomentum(fastjet::PseudoJet fj_tmp, double thisJEC) {
  return Vect4M(thisJEC*fj_tmp.pt(),fj_tmp.eta(), fj_tmp.phi(),thisJEC*fj_tmp.m());
}

//--------------------------------------------------------------------------------------------------

void FatJetExtenderMod::FillXlSubJets(std::vector<fastjet::PseudoJet> &fjSubJets,
                                 XlFatJet *pFatJet, XlSubJet::ESubJetType subJetType)
{
  for (int iSJet=0; iSJet < (int) fjSubJets.size(); iSJet++) {
    XlSubJet *subJet = fXlSubJets[(unsigned int)subJetType]->AddNew();
    subJet->SetRawPtEtaPhiM(fjSubJets[iSJet].pt(),
                            fjSubJets[iSJet].eta(),
                            fjSubJets[iSJet].phi(),
                            fjSubJets[iSJet].m());

    // Store the QG tagging variable
    if (fQGTaggingActive)
      FillSubjetQGTagging(fjSubJets[iSJet], 0.01, subJet, pFatJet);

    // Store the subjet type value
    subJet->SetSubJetType(subJetType);

    // Add the subjet to the relative fatjet
    pFatJet->AddSubJet(subJet,subJetType);

  }

  return;
}

//--------------------------------------------------------------------------------------------------
std::vector <fastjet::PseudoJet>   FatJetExtenderMod::Sorted_by_pt_min_pt(std::vector <fastjet::PseudoJet> &jets,
                                                                     float jetPtMin)
{
  // First order collection by pt
  std::vector<fastjet::PseudoJet> sortedJets = sorted_by_pt(jets);

  // Loop on the sorted collection and erase jets below jetPtMin
  std::vector<fastjet::PseudoJet>::iterator it = sortedJets.begin();
  for ( ;  it != sortedJets.end(); ) {
    if (it->perp() < jetPtMin)
      it = sortedJets.erase(it);
    else
      it++;
  }

  // Return the reduced and sorted jet collection
  return sortedJets;
}

//--------------------------------------------------------------------------------------------------
void FatJetExtenderMod::GetJetConstituents(fastjet::PseudoJet &jet, std::vector <fastjet::PseudoJet> &constits,
                                      float constitsPtMin)
{
  // Loop on input jet constituents vector and discard very soft particles (ghosts)
  for (unsigned int iPart = 0; iPart < jet.constituents().size(); iPart++) {
    if (jet.constituents()[iPart].perp() < constitsPtMin)
      continue;
    constits.push_back(jet.constituents()[iPart]);
  }

  return;
}
//--------------------------------------------------------------------------------------------------
double FatJetExtenderMod::GetQjetVolatility(std::vector <fastjet::PseudoJet> &constits, int QJetsN, int seed)
{
  std::vector<float> qjetmasses;

  double zcut(0.1), dcut_fctr(0.5), exp_min(0.), exp_max(0.), rigidity(0.1), truncationFactor(0.0);

  QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity, truncationFactor);
  fastjet::JetDefinition qjet_def(&qjet_plugin);

  int nFailed = 0;
  for(unsigned int ii = 0 ; ii < (unsigned int) QJetsN ; ii++){
    qjet_plugin.SetRandSeed(seed+ii); // new feature in Qjets to set the random seed
    fastjet::ClusterSequence *qjet_seq =
      new fastjet::ClusterSequence(constits, qjet_def);

    if (!qjet_plugin.succeeded()) {
      delete qjet_seq; // inclusive_jets() not run yet, so no need to call delete_self_when_unused()
      return -(seed+ii);  // this will be the error value for when too many jets are left unmerged...needs more investigation
    }

    vector<fastjet::PseudoJet> inclusive_jets2 = sorted_by_pt(qjet_seq->inclusive_jets(10.0));
    // skip failed recombinations (with no output jets)
    if (inclusive_jets2.size() == 0) {
      nFailed++;
      continue;
    }
    else if (inclusive_jets2.size() > 1) {
      if (inclusive_jets2[1].pt() > inclusive_jets2[0].pt()*0.1) {
      // if more than one jet were found, probably don't trust the mass of the leading one
      // unless the subleading one is really small:
        nFailed++;
        continue;
      }
    }
    if (nFailed*5 > QJetsN) {
      if (inclusive_jets2.size() > 0)
        qjet_seq->delete_self_when_unused();
      delete qjet_seq;
      // if more than a fifth of the iterations fail, let's just give up
      return -1;
    }

    qjetmasses.push_back( inclusive_jets2[0].m() );
    // memory cleanup
    if (inclusive_jets2.size() > 0)
      qjet_seq->delete_self_when_unused();
    delete qjet_seq;
  }

  // find RMS of a vector
  float qjetsRMS = FindRMS( qjetmasses );
  // find mean of a vector
  float qjetsMean = FindMean( qjetmasses );
  float qjetsVolatility = qjetsRMS/qjetsMean;
  return qjetsVolatility;
}


//--------------------------------------------------------------------------------------------------
double FatJetExtenderMod::FindRMS(std::vector<float> qjetmasses)
{
  float mean = FindMean(qjetmasses);

  float totalsquared = 0.;
  float ctr = 0.;
  for (unsigned int i = 0; i < qjetmasses.size(); i++){
    totalsquared += (qjetmasses[i] - mean)*(qjetmasses[i] - mean) ;
    ++ctr;
  }
  float RMS = sqrt( totalsquared/ctr );
  return RMS;
}

//--------------------------------------------------------------------------------------------------
double FatJetExtenderMod::FindMean(std::vector<float> qjetmasses)
{
  float total = 0.;
  float ctr = 0.;
  for (unsigned int i = 0; i < qjetmasses.size(); i++){
      total = total + qjetmasses[i];
      ctr++;
  }
  return total/ctr;
}

//--------------------------------------------------------------------------------------------------
void FatJetExtenderMod::FillSubjetQGTagging(fastjet::PseudoJet &jet, float constitsPtMin,
                                       XlSubJet *pSubJet, XlFatJet *pFatJet)
{
  // Prepare a PFJet to compute the QGTagging
  PFJet pfJet;
  pfJet.SetRawPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),jet.m());

  // Loop on input jet constituents vector and discard very soft particles (ghosts)
  for (unsigned int iPart = 0; iPart < jet.constituents().size(); iPart++) {
    if (jet.constituents()[iPart].perp() < constitsPtMin)
      continue;
    int thisPFCandIndex = jet.constituents()[iPart].user_index();
    // Add the constituent to the PF subjet
    pfJet.AddPFCand(pFatJet->PFCand(thisPFCandIndex));
  }

  // Compute the subjet QGTagging
  if (fQGTaggingActive) {
    fQGTagger->CalculateVariables(&pfJet, fVertexes);
    pSubJet->SetQGTag(fQGTagger->QGValue());
    pSubJet->SetQGPtD(fQGTagger->GetPtD());
    pSubJet->SetQGAxis1(fQGTagger->GetAxis1());
    pSubJet->SetQGAxis2(fQGTagger->GetAxis2());
    pSubJet->SetQGMult(fQGTagger->GetMult());
  }

  return;
}

//--------------------------------------------------------------------------------------------------
TVector2 FatJetExtenderMod::GetPull(fastjet::PseudoJet &jet, float constitsPtMin)
{
  double dYSum   = 0;
  double dPhiSum = 0;
  // Loop on input jet constituents vector and discard very soft particles (ghosts)
  for (unsigned int iPart = 0; iPart < jet.constituents().size(); iPart++) {
    if (jet.constituents()[iPart].perp() < constitsPtMin)
      continue;
    double dY     = jet.constituents()[iPart].rapidity()-jet.rapidity();
    double dPhi   = MathUtils::DeltaPhi(jet.constituents()[iPart].phi(),jet.phi());
    double weight = jet.constituents()[iPart].pt()*sqrt(dY*dY + dPhi*dPhi);
    dYSum   += weight*dY;
    dPhiSum += weight*dPhi;
  }
  return TVector2(dYSum/jet.pt(), dPhiSum/jet.pt());
}

//--------------------------------------------------------------------------------------------------
double FatJetExtenderMod::GetPullAngle(std::vector<fastjet::PseudoJet> &fjSubJets, float constitsPtMin)
{
  // Subject collection already sorted by pt
  // Consider only the leading and the subleading for the pull angle computation
  // work in dy-dphi space of leading subjet

  // Exclude cases where there is no second subjet (input jet made by one particle)
  if (fjSubJets.size() < 2)
    return -20.;

  TVector2 lPull = GetPull(fjSubJets[0],constitsPtMin);
  TVector2 lJet(fjSubJets[1].rapidity()-fjSubJets[0].rapidity(),
                MathUtils::DeltaPhi(fjSubJets[1].phi(), fjSubJets[0].phi()));
  double lThetaP = lPull.DeltaPhi(lJet);
  return lThetaP;
}
