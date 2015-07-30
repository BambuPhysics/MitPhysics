#include "MitAna/DataTree/interface/Names.h"
#include "MitPhysics/Mods/interface/FatJetExtenderMod.h"
#include "MitAna/DataTree/interface/PFJetCol.h"

#include "MitCommon/DataFormats/interface/Vect4M.h"
#include "MitCommon/DataFormats/interface/Vect3.h"
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/Utils/interface/Utils.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/PhysicsUtils/interface/CMSTopTagger.h"
#include "MitAna/PhysicsUtils/interface/HEPTopTagger.h"
#include "QjetsPlugin.h"
#include "Qjets.h"

#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/SoftDrop.hh"

using namespace mithep;

ClassImp(mithep::FatJetExtenderMod)

//--------------------------------------------------------------------------------------------------
FatJetExtenderMod::FatJetExtenderMod(const char *name, const char *title) :
  BaseMod (name,title),
  fIsData (kTRUE),
  fQGTaggingActive (kTRUE),
  fQGTaggerCHS (kFALSE),
  fPublishOutput (kTRUE),
  fFatJetsName (""),
  fFatJetsFromBranch (kTRUE),
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
  fSoftDropR0 (1.),
  fPruneZCut (0.1),
  fPruneDistCut (0.5),
  fFilterN (3),
  fFilterRad (0.2),
  fTrimRad (0.2),
  fTrimPtFrac (0.05),
  fConeSize (0.6),
  fDeconstruct(0),
  fDoSDSubJets(kFALSE),
  fDoPrunedSubJets(kFALSE),
  fDoTrimmedSubJets(kFALSE),
  fDoCMSTT(kFALSE),
  fDoHEPTT(kFALSE),
  fProcessNJets (3),
  fDoShowerDeconstruction(kTRUE)
{
  // Constructor.

    // setup subjet branch names
    fXlSubJetsName[0] = "SoftDropSubjets";
    fXlSubJetsName[1] = "PrunedSubjets";
    fXlSubJetsName[2] = "TrimmedSubjets";
    fXlSubJetsName[3] = "CMSTTSubjets";
    fXlSubJetsName[4] = "HEPTTSubjets";
    fXlSubJetsName[5] = "NjettinessSubjets";


}

FatJetExtenderMod::~FatJetExtenderMod()
{
  // Destructor
  if (fXlSubJets){
    for(int i=0; i<XlSubJet::nSubJetTypes; ++i)
      delete fXlSubJets[i];
  }

  if (fXlFatJets)
    delete fXlFatJets;

  delete fPruner;
  delete fFilterer;
  delete fTrimmer ;

  delete fCAJetDef;

  delete fActiveArea;
  delete fAreaDefinition;

  delete fQGTagger;

  delete fDeconstruct;
}

//--------------------------------------------------------------------------------------------------
void FatJetExtenderMod::Process()
{
  // make sure the out collections are empty before starting
  fXlFatJets->Delete();
  for(unsigned int i=0; i<XlSubJet::nSubJetTypes; ++i)
    fXlSubJets[i]->Delete();

  fFatJets = GetObject<FatJetCol>(fFatJetsName);
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

  // set up shower deconstruction stuff
  TString inputCard = Utils::GetEnv("CMSSW_BASE");
  inputCard += TString::Format("/src/MitAna/SDAlgorithm/config/input_card_%i.dat",int(fConeSize*10));
  AnalysisParameters param(inputCard.Data());
  Deconstruction::TopGluonModel *signal = new Deconstruction::TopGluonModel(param);
  Deconstruction::BackgroundModel *background = new Deconstruction::BackgroundModel(param);
  Deconstruction::ISRModel *isr = new Deconstruction::ISRModel(param);
  fDeconstruct = new Deconstruction::Deconstruct(param, *signal, *background, *isr);

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

    // perform Nsubjettiness analysis and fill the extended XlFatJet object
    // this method will also fill the SubJet collection
    FillXlFatJet(jet);

  }

  return;
}

//--------------------------------------------------------------------------------------------------
void FatJetExtenderMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis.
  ReqEventObject(fFatJetsName,fFatJets,fFatJetsFromBranch);
  ReqEventObject(fPfCandidatesName,fPfCandidates,fPfCandidatesFromBranch
  ReqEventObject(fPileUpDenName,fPileUpDen,fPileUpDenFromBranch);
  ReqEventObject(fVertexesName,fVertexes,fVertexesFromBranch);

  // Initialize area caculation (done with ghost particles)

  // Create the new output collection
  fXlFatJets = new XlFatJetArr(16,fXlFatJetsName);
  for(ESubJetType i = XlSubJet::kSoftDrop; i<XlSubJet::nSubJetTypes; ++i) {
    // only allocate memory for the subjets that are turned on
    if (fSubJetFlags & (1<<i))
      fXlSubJets[(unsigned int)i] = new XlSubJetArr(16,fXlSubJetsName[(unsigned int)i]);
  }
  // Publish collection for further usage in the analysis
  if (fPublishOutput) {
    PublishObj(fXlFatJets);
    for(ESubJetType i = XlSubJet::kSoftDrop; i<XlSubJet::nSubJetTypes; ++i) {
      if (fSubJetFlags & (1<<i))
        PublishObj(fXlSubJets[(unsigned int)i]);
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

  // Initialize area caculation (done with ghost particles)
  int activeAreaRepeats = 1;
  double ghostArea = 0.01;
  double ghostEtaMax = 7.0;
  fActiveArea = new fastjet::GhostedAreaSpec(ghostEtaMax,activeAreaRepeats,ghostArea);
  fAreaDefinition = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts,*fActiveArea);

  // Initialize QGTagger class
  fQGTagger = new QGTagger(fQGTaggerCHS);

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
  XlFatJet *xlFatJet = fXlFatJets->Allocate();
  new (xlFatJet) XlFatJet(*fatJet);

  // Prepare and store QG tagging info
  float qgValue = -1.;
  if (fQGTaggingActive) {
    fQGTagger->CalculateVariables(fatJet, fVertexes);
    qgValue = fQGTagger->QGValue();
  }
  xlFatJet->SetQGTag(qgValue);


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

    if (fjOutJets.size() > 0)
      fjClustering->delete_self_when_unused();
    delete fjClustering;

    return;
  }
  fastjet::PseudoJet fjJet = fjOutJets[0];

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

  // Compute Q-jets volatility
  std::vector<fastjet::PseudoJet> constits;
  GetJetConstituents(fjJet, constits, 0.01);
  double QJetVol = GetQjetVolatility(constits, 25, fCounter*25);
  fCounter++;
  constits.clear();

  // do grooming and subjetting
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

  fjClusteredJets[XlSubJet::kPruned] = (*fPruner)(fjJet);
  fjClusteredJets[XlSubJet::kTrimmed] = (*fTrimmer)(fjJet);
  double MassPruned = fjClusteredJets[XlSubJet::kPruned].m();
  double MassFiltered = ((*fFilterer)(fjJet)).m();
  double MassTrimmed = fjClusteredJets[XlSubJet::kTrimmed].m();

  // do CMS and HEP top tagging
  fastjet::CMSTopTagger* fCMSTopTagger = new fastjet::CMSTopTagger();
  std::vector<fastjet::PseudoJet> lOutJets = sorted_by_pt(fjClustering->inclusive_jets(0.0));
  fastjet::PseudoJet iJet = lOutJets[0];
  HEPTopTagger hepTopJet = HEPTopTagger(*fjClustering,iJet);;
  fastjet::PseudoJet cmsTopJet = fCMSTopTagger->result(iJet);


  // fill subjets
  Bool_t computedPullAngle = kFALSE;
  for (unsigned int iSJType = 0; iSJType!=XlSubJet::nSubJetTypes; ++iSJType) {
    if (fSubJetFlags & ~(1<<iSJType))
      continue;
    // okay, we really do want to save this collection
    std::vector<fastjet::PseudoJet> fjSubjets;
    if (iSJType <= XlSubJet::kTrimmed) {
      // these have a common interface
      int nSubJets = std::min<unsigned int>(fjClusteredJets[iSJType].constituents().size(),3);
      fjSubjets = fjClusteredJets[iSJType].associated_cluster_sequence()->exclusive_subjets(fjClusteredJets[iSJType],nSubJets);
      std::vector<fastjet::PseudoJet> fjSubJetsSorted = Sorted_by_pt_min_pt(fjTopSubJets,0.01);
      if (!computedPullAngle) {
        xlFatJet->SetPullAngle(GetPullAngle(fjSubJetsSorted,0.01));
      }
      FillXlSubJets(fjSubJetsSorted,xlFatJet,(ESubJetType)iSJType);
    } else if (iSJType == XlSubJet::kCMSTT) {
      fjSubjets = cmsTopJet.pieces();
      FillXlSubJets(fjSubjets,xlFatJet,XlSubJet::kCMSTT);
    } else if (iSJType == XlSubJet::kHEPTT) {
      fjSubjets = hepTopJet.top_subjets();
      FillXlSubJets(fjSubjets,xlFatJet,XlSubJet::kHEPTT);
    }
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
      if (microjets.size()>9)
        microjets.erase(microjets.begin()+9,microjets.end());
      double Psignal = 0.0;
      double Pbackground = 0.0;
      try {
        double chi = fDeconstruct->deconstruct(microjets, Psignal, Pbackground);
        xlFatJet->SetChi(chi);
      } catch(Deconstruction::Exception &e) {
        std::cout << "Exception while running SD: " << e.what() << std::endl;
      }
  }

  // ---- Fastjet is done ----
  double thisJEC = xlFatJet->Pt()/xlFatJet->RawMom().Pt();

  // Store the energy correlation values
  xlFatJet->SetC2b0(C2b0);
  xlFatJet->SetC2b0p2(C2b0p2);
  xlFatJet->SetC2b0p5(C2b0p5);
  xlFatJet->SetC2b1(C2b1);
  xlFatJet->SetC2b2(C2b2);

  // Store the groomed masses, apply JEC
  xlFatJet->SetMassSDb0(MassSDb0*thisJEC);
  xlFatJet->SetMassSDb1(MassSDb1*thisJEC);
  xlFatJet->SetMassSDb2(MassSDb2*thisJEC);
  xlFatJet->SetMassSDbm1(MassSDbm1*thisJEC);
  xlFatJet->SetMassPruned(MassPruned*thisJEC);
  xlFatJet->SetMassFiltered(MassFiltered*thisJEC);
  xlFatJet->SetMassTrimmed(MassTrimmed*thisJEC);

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

  return;
}

//--------------------------------------------------------------------------------------------------

inline Vect4M FatJetExtender::GetCorrectedMomentum(fastjet::PseudoJet fj_tmp, double thisJEC) {
  return Vect4M(thisJEC*fj_tmp.pt(),fj_tmp.eta(), fj_tmp.phi(),thisJEC*fj_tmp.m());
}

//--------------------------------------------------------------------------------------------------

void FatJetExtenderMod::FillXlSubJets(std::vector<fastjet::PseudoJet> &fjSubJets,
                                 XlFatJet *pFatJet, XlSubJet::ESubJetType subJetType)
{
  for (int iSJet=0; iSJet < (int) fjSubJets.size(); iSJet++) {
    XlSubJet *subJet = fXlSubJets[(unsigned int)subJetType]->Allocate();
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
    pFatJet->AddSubJet(subJet);

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
