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

#include "fastjet/contrib/NjettinessDefinition.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "MitAna/PhysicsUtils/interface/HEPTopTagger.h"

#include <map>

using namespace mithep;

ClassImp(mithep::FatJetExtenderMod)

//--------------------------------------------------------------------------------------------------
FatJetExtenderMod::FatJetExtenderMod(const char *name, const char *title) :
BaseMod(name,title),
  fIsData(kFALSE),
  fQGTaggingActive(kTRUE),
  fQGTaggerCHS(kFALSE),
  fPublishOutput(kTRUE),
  fFatJetsName(""),
  fFatJetsFromBranch(kFALSE),
  fFatJets(0),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPFCandidates(0),
  fPileUpDenName(Names::gkPileupEnergyDensityBrn),
  fPileUpDenFromBranch(kTRUE),
  fPileUpDen(0),
  fVertexesName("GoodVertexes"),
  fVertexesFromBranch(kFALSE),
  fVertexes(0),
  fXlFatJetsName("XlFatJets"),
  fUseSoftDropLib(kTRUE),
  fSoftDropZCut(0.1),
  fSoftDropR0(.8),
  fPruneZCut(0.1),
  fPruneDistCut(0.5),
  fFilterN(3),
  fFilterRad(0.2),
  fTrimRad(0.2),
  fTrimPtFrac(0.05),
  fConeSize(0.8),
  fInputCard(""),
  fProcessNJets(4),
  fDoShowerDeconstruction(kFALSE),
  fBeVerbose(kFALSE),
  fDoECF(kFALSE),
  fDoQjets(kFALSE),
  fNMaxMicrojets(5),
  fNQjets(25),
  fJetAlgo(kAntiKt),
  fDoCMSandHTT(kFALSE),
  fCorrLevel(mNone),
  fCorrector(0)
{
  // Constructor.
  fJECFiles.clear();
}

FatJetExtenderMod::~FatJetExtenderMod()
{
}

//--------------------------------------------------------------------------------------------------
void FatJetExtenderMod::Process()
{
  if (fBeVerbose)
    fStopwatch->Start(kTRUE);

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
    if (!jet) {
      Error("Pocess", "Jets provided are not FatJets.");
      break;
    }

    // mark jet (and consequently its consituents) for further use in skim
    jet->Mark();
    if (fBeVerbose) {
      Info("Process", "Finished setup in %f seconds\n",fStopwatch->RealTime());
      fStopwatch->Start();
    }
    FillXlFatJet(jet);
  }
}

//--------------------------------------------------------------------------------------------------
void FatJetExtenderMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis.
  // setup subjet branch names
  fXlSubJetsName[0] = fXlFatJetsName+"_SoftDropSubjets";
  fXlSubJetsName[1] = fXlFatJetsName+"_PrunedSubjets";
  fXlSubJetsName[2] = fXlFatJetsName+"_TrimmedSubjets";
  fXlSubJetsName[3] = fXlFatJetsName+"_CMSTTSubjets";
  fXlSubJetsName[4] = fXlFatJetsName+"_HEPTTSubjets";
  fXlSubJetsName[5] = fXlFatJetsName+"_NjettinessSubjets";

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

  fNJettiness = new fastjet::contrib::Njettiness(fastjet::contrib::Njettiness::onepass_kt_axes, fastjet::contrib::NormalizedCutoffMeasure(1., fConeSize, 10000.));

  if (fDoECF) {
    fECR2b0 = new fastjet::contrib::EnergyCorrelatorDoubleRatio(2,0. ,fastjet::contrib::EnergyCorrelator::pt_R);
    fECR2b0p2 = new fastjet::contrib::EnergyCorrelatorDoubleRatio(2,0.2,fastjet::contrib::EnergyCorrelator::pt_R);
    fECR2b0p5 = new fastjet::contrib::EnergyCorrelatorDoubleRatio(2,0.5,fastjet::contrib::EnergyCorrelator::pt_R);
    fECR2b1 = new fastjet::contrib::EnergyCorrelatorDoubleRatio(2,1.0,fastjet::contrib::EnergyCorrelator::pt_R);
    fECR2b2 = new fastjet::contrib::EnergyCorrelatorDoubleRatio(2,2.0,fastjet::contrib::EnergyCorrelator::pt_R);
  }

  if (fUseSoftDropLib) {
    fSoftDrop[kSD0] = new fastjet::contrib::SoftDrop(0.0, fSoftDropZCut, fSoftDropR0);
    fSoftDrop[kSD1] = new fastjet::contrib::SoftDrop(1.0, fSoftDropZCut, fSoftDropR0);
    fSoftDrop[kSD2] = new fastjet::contrib::SoftDrop(2.0, fSoftDropZCut, fSoftDropR0);
    fSoftDrop[kSDm1] = new fastjet::contrib::SoftDrop(-1.0, fSoftDropZCut, fSoftDropR0);
    fSoftDrop[kSD0]->set_tagging_mode();
    fSoftDrop[kSD1]->set_tagging_mode();
    fSoftDrop[kSD2]->set_tagging_mode();
    fSoftDrop[kSDm1]->set_tagging_mode();
  }
  else {
    fSoftDropCalc = new SoftDropCalculator(fSoftDropZCut, fSoftDropR0);
    fSoftDropCalc->addBeta(0.);
    fSoftDropCalc->addBeta(1.);
    fSoftDropCalc->addBeta(2.);
    fSoftDropCalc->addBeta(-1.);
  }

  // Prepare pruner
  fPruner = new fastjet::Pruner(fastjet::cambridge_algorithm,fPruneZCut,fPruneDistCut);
  // Prepare filterer
  fFilterer = new fastjet::Filter(fastjet::JetDefinition(fastjet::cambridge_algorithm,fFilterRad),
                                  fastjet::SelectorNHardest(fFilterN));
  // Prepare trimmer
  fTrimmer = new fastjet::Filter(fastjet::JetDefinition(fastjet::kt_algorithm,fTrimRad),
                                 fastjet::SelectorPtFractionMin(fTrimPtFrac));

  if (fJetAlgo == kCambridgeAachen) {
    fJetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, fConeSize);
    fCAJetDef = fJetDef;
  } else if (fJetAlgo == kAntiKt) {
    fJetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, fConeSize);
    fCAJetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, fConeSize);
  } else if (fJetAlgo == kKt) {
    fJetDef = new fastjet::JetDefinition(fastjet::kt_algorithm, fConeSize);
    fCAJetDef = new fastjet::JetDefinition(fastjet::cambridge_algorithm, fConeSize);
  }

  if (fDoCMSandHTT)
    fCMSTopTagger = new fastjet::CMSTopTagger();

  // Initialize area caculation (done with ghost particles)
  int activeAreaRepeats = 1;
  double ghostArea = 0.01;
  double ghostEtaMax = 5.0;
  fActiveArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
  fAreaDefinition = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*fActiveArea);

  // Initialize QGTagger class
  if (fQGTaggingActive)
    fQGTagger = new QGTagger(fQGTaggerCHS);

  if (fDoShowerDeconstruction) {
    // set up shower deconstruction stuff
    if (fInputCard=="") {
      // default was never changed
      fInputCard = Utils::GetEnv("MIT_DATA");
      fInputCard += TString::Format("/SDAlgorithm/input_card_%i.dat",int(fConeSize*10));
    }
    fParam = new AnalysisParameters(fInputCard.Data());
    fSignal = new Deconstruction::TopGluonModel(*fParam);
    fBackground = new Deconstruction::BackgroundModel(*fParam);
    fISR = new Deconstruction::ISRModel(*fParam);
    fDeconstruct = new Deconstruction::Deconstruct(*fParam, *fSignal, *fBackground, *fISR);
  }

  if (fCorrLevel==mL2L3) {
    if (fJECFiles.size()==0)
      SendError(kAbortAnalysis,"SlaveBegin","L2/L3 corrections for groomed fatjets are requested but parameter files not provided.");
    // only case we need to redo corrections
    fCorrector = new JetCorrector();
    for (auto filePath : fJECFiles) {
      fCorrector->AddParameterFile(filePath.Data());
    }
  }

  fStopwatch = new TStopwatch();
}

//--------------------------------------------------------------------------------------------------
void FatJetExtenderMod::SlaveTerminate()
{
  RetractObj(fXlFatJets->GetName());
  for(int i=0; i<XlSubJet::nSubJetTypes; ++i) {
    if (fSubJetFlags & (1<<i))
      RetractObj(fXlSubJets[i]->GetName());
  }

  // Destructor
  if (fXlSubJets){
    for(int i=0; i<XlSubJet::nSubJetTypes; ++i) {
      if (fSubJetFlags & (1<<i))
        delete fXlSubJets[i];
    }
  }

  delete fXlFatJets;

  delete fNJettiness;

  delete fECR2b0;
  delete fECR2b0p2;
  delete fECR2b0p5;
  delete fECR2b1;
  delete fECR2b2;

  delete fPruner;
  delete fFilterer;
  delete fTrimmer ;

  for (auto* fsd : fSoftDrop)
    delete fsd;
  delete fSoftDropCalc;

  if (fCAJetDef!=fJetDef)
    delete fCAJetDef;
  delete fJetDef;

  delete fActiveArea;
  delete fAreaDefinition;
  delete fCMSTopTagger;

  delete fQGTagger;

  delete fDeconstruct;
  delete fParam;
  delete fSignal;
  delete fBackground;
  delete fISR;

  if (fCorrector)
    delete fCorrector;

  delete fStopwatch;
}

//--------------------------------------------------------------------------------------------------
void FatJetExtenderMod::FillXlFatJet(const FatJet *fatJet)
{
  // Prepare and store in an array a new FatJet
  XlFatJet *xlFatJet = fXlFatJets->Allocate();
  new (xlFatJet) XlFatJet(*fatJet);
  // Prepare and store QG tagging info
  float qgValue = -1.;
  if (fQGTaggingActive) {
    fQGTagger->CalculateVariables(fatJet, fVertexes);
    qgValue = fQGTagger->QGValue();
  }
  xlFatJet->SetQGTag(qgValue);

  if (fBeVerbose) {
    Info("FillXlFatJet", "Finished QG-tagging in %f seconds\n", fStopwatch->RealTime());
    fStopwatch->Start();
  }

  VPseudoJet fjParts;
  // Push all particle flow candidates of the input PFjet into fastjet particle collection
  for (UInt_t j=0; j<fatJet->NPFCands(); ++j) {
    const PFCandidate *pfCand = fatJet->PFCand(j);
    fjParts.emplace_back(pfCand->Px(),pfCand->Py(),pfCand->Pz(),pfCand->E());
    fjParts.back().set_user_index(j);
  }


  // Setup the cluster for fastjet
  fastjet::ClusterSequenceArea fjClustering(fjParts, *fJetDef, *fAreaDefinition);
  fastjet::ClusterSequenceArea fjCAClustering(fjParts, *fCAJetDef, *fAreaDefinition);

  // ---- Fastjet is ready ----

  VPseudoJet&& allJets(fjClustering.inclusive_jets(0.));
  VPseudoJet&& allCAJets(fjCAClustering.inclusive_jets(0.));

  // Consider only the hardest jet of the output collection
  // For nsubjettiness etc use the hardest above 10 GeV
  fastjet::PseudoJet* maxPtJet0 = 0;
  fastjet::PseudoJet* maxPtJet10 = 0;
  for (auto& jet : allJets) {
    double pt2 = jet.perp2();
    if (!maxPtJet0 || pt2 > maxPtJet0->perp2())
      maxPtJet0 = &jet;

    if (pt2 < 100.)
      continue;
    if (!maxPtJet10 || pt2 > maxPtJet10->perp2())
      maxPtJet10 = &jet;
  }
  fastjet::PseudoJet* maxCAPtJet0 = 0;
  fastjet::PseudoJet* maxCAPtJet10 = 0;
  if (fJetAlgo==kCambridgeAachen) {
    //  no need to recompute anything
    maxCAPtJet0 = maxPtJet0;
    maxCAPtJet10 = maxPtJet10;
  } else {
    for (auto& jet : allCAJets) {
      double pt2 = jet.perp2();
      if (!maxCAPtJet0 || pt2 > maxCAPtJet0->perp2())
        maxCAPtJet0 = &jet;

      if (pt2 < 100.)
        continue;
      if (!maxCAPtJet10 || pt2 > maxCAPtJet10->perp2())
        maxCAPtJet10 = &jet;
    }
  }


  // Check that the output collection size is non-null, otherwise nothing to be done further
  if (!maxPtJet10 || !maxCAPtJet10) {
    Warning("FillXlFatJet", "Input FatJet produces null reclustering output!\n");
    return;
  }

  if (fBeVerbose) {
    Info("FillXlFatJet", "Finished prepping fastjet in %f seconds",fStopwatch->RealTime());
    fStopwatch->Start();
  }

  VPseudoJet constituents(maxPtJet10->constituents());
  VPseudoJet CAconstituents(maxCAPtJet10->constituents());

  // we have now set up pseudo jets and jet definitions for CA and the XlFatJet algorithm (if different from CA)

  if (fDoECF) {
    // Compute the energy correlation function ratios
    // uses CA
    double C2b0   = (*fECR2b0)(*maxCAPtJet10);
    double C2b0p2 = (*fECR2b0p2)(*maxCAPtJet10);
    double C2b0p5 = (*fECR2b0p5)(*maxCAPtJet10);
    double C2b1   = (*fECR2b1)(*maxCAPtJet10);
    double C2b2   = (*fECR2b2)(*maxCAPtJet10);
    xlFatJet->SetC2b0(C2b0);
    xlFatJet->SetC2b0p2(C2b0p2);
    xlFatJet->SetC2b0p5(C2b0p5);
    xlFatJet->SetC2b1(C2b1);
    xlFatJet->SetC2b2(C2b2);

    if (fBeVerbose) {
      Info("FillXlFatJet", "Finished ECF calculation in %f seconds\n",fStopwatch->RealTime());
      fStopwatch->Start();
    }
  }

  if (fDoQjets) {
    // uses jet algo
    VPseudoJet constituentsNoGhost(FilterJetsByPt(constituents, 0.01));
    // Compute Q-jets volatility
    double QJetVol = GetQjetVolatility(constituentsNoGhost, fNQjets, fCounter*fNQjets);
    fCounter++;
    xlFatJet->SetQJetVol(QJetVol);

    if (fBeVerbose) {
      Info("FillXlFatJet", "Finished Qjets in %f seconds\n",fStopwatch->RealTime());
      fStopwatch->Start();
    }
  }

  // do grooming and subjetting
  double thisJEC = 1;
  if (fCorrLevel==mAll)
    thisJEC = xlFatJet->Pt()/xlFatJet->RawMom().Pt();
  else if (fCorrLevel==mL2L3) {
    Jet *l2l3Jet = fatJet->MakeCopy();
    fCorrector->Correct(*l2l3Jet);
    thisJEC = l2l3Jet->Pt()/l2l3Jet->RawMom().Pt();
  }
  fastjet::PseudoJet fjClusteredJets[XlSubJet::kTrimmed+1] = {
    fastjet::PseudoJet(),
    (*fPruner)(*maxPtJet10),
    (*fTrimmer)(*maxPtJet10)
  };

  double MassSDb0, MassSDb1, MassSDb2, MassSDbm1;
  if (fUseSoftDropLib) {
    if (fSoftDropR0 > 1.2) {
      fjClusteredJets[XlSubJet::kSoftDrop] = (*fSoftDrop[kSD1])(*maxPtJet10);
      MassSDb0 = ((*fSoftDrop[kSD0])(*maxPtJet10)).m();
      MassSDb1 = fjClusteredJets[XlSubJet::kSoftDrop].m();
    }
    else {
      fjClusteredJets[XlSubJet::kSoftDrop] = (*fSoftDrop[kSD0])(*maxPtJet10);
      MassSDb0 = fjClusteredJets[XlSubJet::kSoftDrop].m();
      MassSDb1 = ((*fSoftDrop[kSD1])(*maxPtJet10)).m();
    }

    MassSDb2 = ((*fSoftDrop[kSD2])(*maxPtJet10)).m();
    MassSDbm1 = ((*fSoftDrop[kSDm1])(*maxPtJet10)).m();
  }
  else {
    fSoftDropCalc->calculate(*maxPtJet10);

    if (fSoftDropR0 > 1.2) {
      fjClusteredJets[XlSubJet::kSoftDrop] = fSoftDropCalc->result()[1];
      MassSDb0 = fSoftDropCalc->result()[0].m();
      MassSDb1 = fjClusteredJets[XlSubJet::kSoftDrop].m();
    }
    else {
      fjClusteredJets[XlSubJet::kSoftDrop] = fSoftDropCalc->result()[0];
      MassSDb0 = fjClusteredJets[XlSubJet::kSoftDrop].m();
      MassSDb1 = fSoftDropCalc->result()[1].m();
    }

    MassSDb2 = fSoftDropCalc->result()[kSD2].m();
    MassSDbm1 = fSoftDropCalc->result()[kSDm1].m();
  }

  xlFatJet->SetMassSDb0(MassSDb0 * thisJEC);
  xlFatJet->SetMassSDb1(MassSDb1 * thisJEC);
  xlFatJet->SetMassSDb2(MassSDb2 * thisJEC);
  xlFatJet->SetMassSDbm1(MassSDbm1 * thisJEC);

  double MassPruned = fjClusteredJets[XlSubJet::kPruned].m();
  double MassFiltered = ((*fFilterer)(*maxPtJet10)).m();
  double MassTrimmed = fjClusteredJets[XlSubJet::kTrimmed].m();
  xlFatJet->SetMassPruned(MassPruned*thisJEC);
  xlFatJet->SetMassFiltered(MassFiltered*thisJEC);
  xlFatJet->SetMassTrimmed(MassTrimmed*thisJEC);

  if (fBeVerbose) {
    Info("FillXlFatJet", "Finished grooming in %f seconds\n",fStopwatch->RealTime());
    fStopwatch->Start();
  }

  // internal implementation of Nsubjettiness::result
  // getTau is getTauComponents().tau(), and there are bunch of unnecessary calculation done in the latter function. Can in principle go further here.
  // uses CA
  double tau1 = fNJettiness->getTau(1, CAconstituents);
  double tau2 = fNJettiness->getTau(2, CAconstituents);
  double tau3 = fNJettiness->getTau(3, CAconstituents);
  double tau4 = fNJettiness->getTau(4, CAconstituents);
  xlFatJet->SetTau1(tau1);
  xlFatJet->SetTau2(tau2);
  xlFatJet->SetTau3(tau3);
  xlFatJet->SetTau4(tau4);
  
  if (fBeVerbose) {
    Info("FillXlFatJet", "Finished njettiness calculation in %f seconds",fStopwatch->RealTime());
    fStopwatch->Start();
  }

  // fill subjets
  Bool_t computedPullAngle = kFALSE;
  for (unsigned int iSJType = 0; iSJType!=XlSubJet::nSubJetTypes; ++iSJType) {
    if (!(fSubJetFlags & (1<<iSJType)))
      continue;

    // okay, we really do want to save this collection
    if (iSJType <= XlSubJet::kTrimmed) {
      // these have a common interface
      if (!fjClusteredJets[iSJType].has_constituents())
        // if subjet finding failed, skip this step
        continue;

      int nSubJets = fjClusteredJets[iSJType].constituents().size();
      if (nSubJets > 3)
        nSubJets = 3;

      if (nSubJets != 0) {
        auto&& fjSubjets(fjClusteredJets[iSJType].associated_cluster_sequence()->exclusive_subjets(fjClusteredJets[iSJType], nSubJets));
        VPseudoJetPtr fjSubJetsSorted(Sorted_by_pt_min_pt(fjSubjets, 0.01));
        if (!computedPullAngle) {
          xlFatJet->SetPullAngle(GetPullAngle(fjSubJetsSorted, 0.01));
          computedPullAngle = kTRUE;
        }
        FillXlSubJets(fjSubJetsSorted, xlFatJet, (ESubJetType)iSJType);
      }
    }
  }

  if (fBeVerbose) {
    Info("FillXlFatJet", "Finished filling subjets in %f seconds\n",fStopwatch->RealTime());
    fStopwatch->Start();
  }

  // take a shower
  if (fDoShowerDeconstruction) {
    // shower deconstruction
    double microconesize;
    if  (fMicrojetConeSize<0){
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
      microconesize = fMicrojetConeSize;
    }
    fastjet::JetDefinition reclustering(fastjet::JetAlgorithm::kt_algorithm, microconesize);
    fastjet::ClusterSequence * cs_micro = new fastjet::ClusterSequence(fjParts, reclustering);
    VPseudoJet microjets = fastjet::sorted_by_pt(cs_micro->inclusive_jets(10.));
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
      Info("FillXlFatJet", "Finished shower deconstruction in %f seconds\n",fStopwatch->RealTime());
      fStopwatch->Start();
    }

    if (cs_micro->inclusive_jets(0.).size()>0)
      cs_micro->delete_self_when_unused();
    delete cs_micro;
  }

  // Store groomed 4-momenta, apply JEC
  fastjet::PseudoJet fj_tmp;
  if (fSubJetFlags & (1<<XlSubJet::kSoftDrop))
    xlFatJet->SetSoftDropP(GetCorrectedMomentum(fjClusteredJets[XlSubJet::kSoftDrop], thisJEC));
  if (fSubJetFlags & (1<<XlSubJet::kPruned))
    xlFatJet->SetPrunedP(GetCorrectedMomentum(fjClusteredJets[XlSubJet::kPruned], thisJEC));
  if (fSubJetFlags & (1<<XlSubJet::kTrimmed))
    xlFatJet->SetTrimmedP(GetCorrectedMomentum(fjClusteredJets[XlSubJet::kTrimmed], thisJEC));

  // Store the color pull
  xlFatJet->SetPull(GetPull(*maxPtJet10,0.01).Mod());

  // Trim the output collections
  fXlFatJets->Trim();

  // Memory cleanup
  if (fBeVerbose) {
    Info("Process","Finished filling and cleanup in %f seconds\n",fStopwatch->RealTime());
  }

  return;
}

//--------------------------------------------------------------------------------------------------

inline Vect4M FatJetExtenderMod::GetCorrectedMomentum(fastjet::PseudoJet const& fj_tmp, double thisJEC) {
  return Vect4M(thisJEC*fj_tmp.pt(),fj_tmp.eta(), fj_tmp.phi(),thisJEC*fj_tmp.m());
}

//--------------------------------------------------------------------------------------------------

void FatJetExtenderMod::FillXlSubJets(VPseudoJetPtr const& fjSubJets,
                                      XlFatJet* pFatJet, XlSubJet::ESubJetType subJetType)
{
  for (auto* fjSubJet : fjSubJets) {
    XlSubJet *subJet = fXlSubJets[subJetType]->AddNew();
    subJet->SetRawPtEtaPhiM(fjSubJet->pt(),
                            fjSubJet->eta(),
                            fjSubJet->phi(),
                            fjSubJet->m());

    // Store the QG tagging variable
    if (fQGTaggingActive)
      FillSubjetQGTagging(*fjSubJet, 0.01, subJet, pFatJet);

    // Store the subjet type value
    subJet->SetSubJetType(subJetType);

    // Add the subjet to the relative fatjet
    pFatJet->AddSubJet(subJet,subJetType);
  }

  return;
}

//--------------------------------------------------------------------------------------------------
FatJetExtenderMod::VPseudoJetPtr
FatJetExtenderMod::Sorted_by_pt_min_pt(VPseudoJet const& jets, float jetPtMin)
{
  std::map<double, fastjet::PseudoJet const*> sorter;
  for (auto&& jet : jets)
    sorter.emplace(jet.perp2(), &jet);

  double pt2Min = jetPtMin * jetPtMin;

  VPseudoJetPtr sortedJets;
  auto rEnd(sorter.rend());
  for (auto rItr(sorter.rbegin()); rItr != rEnd; ++rItr) {
    if (rItr->second->perp2() < pt2Min)
      break;
    sortedJets.push_back(rItr->second);
  }

  // Return the reduced and sorted jet collection
  return sortedJets;
}

//--------------------------------------------------------------------------------------------------
FatJetExtenderMod::VPseudoJet
FatJetExtenderMod::FilterJetsByPt(VPseudoJet const& jets, double ptMin)
{
  // Loop on input jet constituents vector and discard very soft particles (ghosts)
  VPseudoJet result;
  double ptMin2 = ptMin * ptMin;
  for (auto&& jet : jets) {
    if (jet.perp2() < ptMin2)
      continue;
    result.emplace_back(jet);
  }

  return result;
}
//--------------------------------------------------------------------------------------------------
double FatJetExtenderMod::GetQjetVolatility(VPseudoJet const& constits, int QJetsN, int seed)
{
  std::vector<float> qjetmasses;

  double zcut(0.1), dcut_fctr(0.5), exp_min(0.), exp_max(0.), rigidity(0.1), truncationFactor(0.0);

  QjetsPlugin qjet_plugin(zcut, dcut_fctr, exp_min, exp_max, rigidity, truncationFactor);
  fastjet::JetDefinition qjet_def(&qjet_plugin);

  int nFailed = 0;
  for(unsigned int ii = 0 ; ii < (unsigned int) QJetsN ; ii++){
    qjet_plugin.SetRandSeed(seed+ii); // new feature in Qjets to set the random seed
    fastjet::ClusterSequence qjet_seq(constits, qjet_def);

    if (!qjet_plugin.succeeded()) {
      return -(seed+ii);  // this will be the error value for when too many jets are left unmerged...needs more investigation
    }

    VPseudoJet inclusive_jets2 = fastjet::sorted_by_pt(qjet_seq.inclusive_jets(10.0));
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
      // if more than a fifth of the iterations fail, let's just give up
      return -1;
    }

    qjetmasses.push_back( inclusive_jets2[0].m() );
    // memory cleanup
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
void FatJetExtenderMod::FillSubjetQGTagging(fastjet::PseudoJet const& jet, float constitsPtMin,
                                            XlSubJet *pSubJet, XlFatJet const* pFatJet)
{
  // Prepare a PFJet to compute the QGTagging
  PFJet pfJet;
  pfJet.SetRawPtEtaPhiM(jet.pt(),jet.eta(),jet.phi(),jet.m());

  // Loop on input jet constituents vector and discard very soft particles (ghosts)
  VPseudoJet constituents(jet.constituents());
  double ptMin2 = constitsPtMin * constitsPtMin;

  for (auto&& cons : constituents) {
    if (cons.perp2() < ptMin2)
      continue;
    int thisPFCandIndex = cons.user_index();
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
}

//--------------------------------------------------------------------------------------------------
TVector2 FatJetExtenderMod::GetPull(fastjet::PseudoJet const& jet, float constitsPtMin)
{
  double dYSum   = 0;
  double dPhiSum = 0;
  // Loop on input jet constituents vector and discard very soft particles (ghosts)
  VPseudoJet constituents(jet.constituents());
  double ptMin2 = constitsPtMin * constitsPtMin;

  for (auto&& cons : constituents) {
    if (cons.perp2() < ptMin2)
      continue;

    double dY     = cons.rapidity()-jet.rapidity();
    double dPhi   = MathUtils::DeltaPhi(cons.phi(),jet.phi());
    double weight = cons.pt()*sqrt(dY*dY + dPhi*dPhi);
    dYSum   += weight*dY;
    dPhiSum += weight*dPhi;
  }

  return TVector2(dYSum/jet.pt(), dPhiSum/jet.pt());
}

//--------------------------------------------------------------------------------------------------
double FatJetExtenderMod::GetPullAngle(VPseudoJetPtr const& fjSubJets, float constitsPtMin)
{
  // Subject collection already sorted by pt
  // Consider only the leading and the subleading for the pull angle computation
  // work in dy-dphi space of leading subjet

  // Exclude cases where there is no second subjet (input jet made by one particle)
  if (fjSubJets.size() < 2)
    return -20.;

  TVector2 lPull = GetPull(*fjSubJets[0], constitsPtMin);
  TVector2 lJet(fjSubJets[1]->rapidity() - fjSubJets[0]->rapidity(),
                MathUtils::DeltaPhi(fjSubJets[1]->phi(), fjSubJets[0]->phi()));
  double lThetaP = lPull.DeltaPhi(lJet);
  return lThetaP;
}

void
FatJetExtenderMod::SoftDropCalculator::addBeta(double beta)
{
  betas_.push_back(beta);
  jetDone_.push_back(false);
  groomedJets_.push_back(fastjet::PseudoJet());
}

void
FatJetExtenderMod::SoftDropCalculator::calculate(PseudoJet const& inJet)
{
  // "parallelization" of RecursiveSymmetryCutBase (base class of SoftDrop) for
  // multiple beta values.
  // taking (mostly) the default setup for SoftDrop, i.e.
  //  _do_reclustering = true
  //  _recluster = 0 (uses C/A)
  //  _mu_cut = infinity
  //  _symmetry_measure = scalar_z
  //  _recursion_choice = larger_pt
  //  _subtractor = 0
  //  _grooming_mode = false (non-default)

  // IMPORTANT: at the moment we care only about the soft drop mass, and therefore
  // the groomed jets are not given the StructureType object.

  jetDone_.assign(jetDone_.size(), false);

  auto&& inJetAlgo(inJet.validated_cs()->jet_def().jet_algorithm());
  fastjet::PseudoJet jet(inJetAlgo == fastjet::cambridge_algorithm ? fastjet::contrib::Recluster(cambridge_algorithm, JetDefinition::max_allowable_R)(inJet) : inJet);

  // variables for tracking what will happen
  fastjet::PseudoJet piece1, piece2;
  
  // now recurse into the jet's structure
  while (jet.has_parents(piece1, piece2)) {
    
    // first sanity check: 
    // - zero or negative pts are not allowed for the input subjet
    // - zero or negative masses are not allowed for configurations
    //   in which the mass will effectively appear in a denominator
    //   (The masses will be checked later)
    if (jet.pt2() <= 0) {
      for (unsigned iJ = 0; iJ != jetDone_.size(); ++iJ) {
        if (!jetDone_[iJ])
          groomedJets_[iJ] = fastjet::PseudoJet();
      }

      return;
    }

    // determine the symmetry parameter
    // min(pt1, pt2)/(pt1+pt2), where the denominator is a scalar sum
    // of the two subjets
    double pt1 = piece1.pt();
    double pt2 = piece2.pt();
    // make sure denominator is non-zero
    if (pt1 + pt2 == 0.) {
      for (unsigned iJ = 0; iJ != jetDone_.size(); ++iJ) {
        if (!jetDone_[iJ])
          groomedJets_[iJ] = fastjet::PseudoJet();
      }

      return;
    }
    double sym = std::min(pt1, pt2) / (pt1 + pt2);

    double cutBase = piece1.squared_distance(piece2) / r0sq_;
    
    bool done = true;
    for (unsigned iJ = 0; iJ != betas_.size(); ++iJ) {
      if (jetDone_[iJ])
        continue;

      // make a tagging decision based on symmetry cut
      if (sym > symmetryCut_ * std::pow(cutBase, 0.5 * betas_[iJ])) {
        groomedJets_[iJ] = jet;
        jetDone_[iJ] = true;
        continue;
      }

      done = false;
    }

    if (done)
      break;
    
    // continue unclustering
    jet = pt1 > pt2 ? piece1 : piece2;
  }

  for (unsigned iJ = 0; iJ != jetDone_.size(); ++iJ) {
    if (!jetDone_[iJ]) {
      groomedJets_[iJ] = fastjet::PseudoJet();
      jetDone_[iJ] = true;
    }
  }
}
