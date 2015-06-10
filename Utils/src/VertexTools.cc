#include "MitPhysics/Utils/interface/VertexTools.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/PhotonTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include <TFile.h>
#include <TVector3.h>
#include <TSystem.h>

using namespace mithep;

VertexTools* VertexTools::meobject = NULL;

//--------------------------------------------------------------------------------------------------
VertexTools::VertexTools() :
  fIsInitMvaM(kFALSE),
  fIsInitMvaP(kFALSE),
  reader(0),
  readervtx(0),
  readerevt(0)
{
  
}

//--------------------------------------------------------------------------------------------------
void VertexTools::InitM(const char* str)  
{
 
  relname = str;
  reader = new TMVA::Reader( "!Color:!Silent" );    
  reader->AddVariable( "var1", &tmvar1 );
  reader->AddVariable( "var2", &tmvar2 );
  reader->AddVariable( "var3", &tmvar3 );
  reader->AddVariable( "var4", &tmvar4 );
  reader->AddVariable( "var5", &tmvar5 );
  reader->AddVariable( "var6", &tmvar6 ); 
  reader->BookMVA( "BDTG method",relname + TString("/src/MitPhysics/data/TMVAClassification_BDTG.weights.xml").Data());
  reader->BookMVA( "BDTD method",relname + TString("/src/MitPhysics/data/TMVAClassification_BDTD.weights.xml" ).Data());
  //reader->BookMVA( "CFMlpANN method", "/home/maxi/cms/root/TMVAClassification_CFMlpANN.weights.xml" );
  reader->BookMVA( "MLP method", relname + TString("/src/MitPhysics/data/TMVAClassification_MLP.weights.xml").Data());
  reader->BookMVA( "MLPBFGS method",relname + TString("/src/MitPhysics/data/TMVAClassification_MLPBFGS.weights.xml" ).Data());
  
  fIsInitMvaM = kTRUE;
  
}

//--------------------------------------------------------------------------------------------------
void VertexTools::InitP(int version)
{

  readervtx = new TMVA::Reader( "!Color:!Silent" );
  readerevt = new TMVA::Reader( "!Color:!Silent" );
  
  readervtx->AddVariable( "ptbal", &fMvaPVars[0] );
  readervtx->AddVariable( "ptasym", &fMvaPVars[1] );
  readervtx->AddVariable( "logsumpt2", &fMvaPVars[2] );
  readervtx->AddVariable( "limPullToConv", &fMvaPVars[3] );
  readervtx->AddVariable( "nConv", &fMvaPVars[4] );
  readervtx->BookMVA( "BDTCat", gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/TMVAClassification_BDTCat_conversions_tmva_407.weights.xml") );

  readerevt->AddVariable( "diphoPt0", &fMvaPEvtVars[0] );
  readerevt->AddVariable( "nVert", &fMvaPEvtVars[1] );
  readerevt->AddVariable( "MVA0", &fMvaPEvtVars[2] );
  readerevt->AddVariable( "MVA1", &fMvaPEvtVars[3] );
  readerevt->AddVariable( "dZ1", &fMvaPEvtVars[4] );
  readerevt->AddVariable( "MVA2", &fMvaPEvtVars[5] );
  readerevt->AddVariable( "dZ2", &fMvaPEvtVars[6] );
  readerevt->AddVariable( "nConv", &fMvaPEvtVars[7] );
  if (version==1) {
    readerevt->BookMVA( "BDTEvt", gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/TMVAClassification_evtBDTG_conversions_tmva_407.weights.xml") );  
  }
  else if (version==2) {
    readerevt->BookMVA( "BDTEvt", gSystem->Getenv("CMSSW_BASE") + TString("/src/MitPhysics/data/TMVAClassification_BDTvtxprob2012.weights.xml") );  
  }
  else assert(0);
  
  fIsInitMvaP = kTRUE;
}


double VertexTools::NewMass(const Photon* ph1, const Photon* ph2, const BaseVertex* vert)
{
  ThreeVector drv1 = (ThreeVector(ph1->CaloPos()) - vert->Position()).Unit();
  //ThreeVector drv1 = (ThreeVector(ph1->SCluster()->Point()) - vert->Position()).Unit();
  FourVector pho1c(drv1.X()*ph1->E(),drv1.Y()*ph1->E(),drv1.Z()*ph1->E(),ph1->E());
  ThreeVector drv2 = (ThreeVector(ph2->CaloPos()) - vert->Position()).Unit();
  FourVector pho2c(drv2.X()*ph2->E(),drv2.Y()*ph2->E(),drv2.Z()*ph2->E(),ph2->E());

  FourVector diboso = pho1c+pho2c;
  return diboso.M();
}

//--------------------------------------------------------------------------------------------------
VertexZarray VertexTools::ExtractZarray(const VertexCol* vcol, float zmin, float zmax,  const BaseVertex  *fBeamSpot)
{
  VertexZarray zs;
  if(vcol == NULL) return zs;

  for(unsigned vv = 0; vv < vcol->GetEntries(); vv++){    
    const Vertex* vert = vcol->At(vv);
    double zpos = vert->Z();
    if(fBeamSpot != NULL) 
      zpos = zpos - fBeamSpot->Z();
    if(zpos > zmin && zpos > zmin)
      zs.push_back(zpos);
  }
  return zs;
}

VertexZarray VertexTools::ExtractZarray(float zmin, float zmax, float step)
{
  VertexZarray zs;
  for(float zpos = zmin; zpos < zmax+step; zpos = zpos+step){
    zs.push_back(zpos);
  }
  return zs;
}

const Vertex* VertexTools::BestVtx( const PFCandidateCol *fPFJets, const VertexCol *c, 
				    const BaseVertex  *fBeamSpot, FourVector diboso) {
  
  if (!c || !c->GetEntries()) return NULL;
  
  double bestprob = -100.;
  const Vertex* bestvert = NULL;
  for(unsigned vv = 0; vv < c->GetEntries(); vv++){
    
    const Vertex* vert = c->At(vv);
    double zpos = vert->Z();     
    double prob = Prob(fPFJets, zpos, fBeamSpot, diboso);
    if(prob > bestprob){
      bestprob = prob;
      bestvert = vert;
    }
  }
  return bestvert;
}

double VertexTools::BestVtx( const PFCandidateCol *fPFJets, VertexZarray zcol,
			     const BaseVertex  *fBeamSpot, FourVector diboso){
  
  double bestprob = -100.;
  double bestz = -100.;
  for(unsigned vv = 0; vv < zcol.size(); vv++){
    double zpos = zcol[vv];     
    double prob = Prob(fPFJets, zpos, fBeamSpot, diboso);
    if(prob > bestprob){
      bestprob = prob;
      bestz = zpos;
    }
  }
  return bestz;
  
}

double VertexTools::Prob(const PFCandidateCol *fPFJets, double zpos,  
			 const BaseVertex  *fBeamSpot, FourVector diboso){
  
  double bosophi = diboso.Phi();
  double bosopt = diboso.Pt();
  
  Vertex* vert = new Vertex(0,0,zpos);
  double ZVC = vert->Z()-fBeamSpot->Z();
  VertexTools* vtool = VertexTools::instance("");
  
  Double_t sinsum = 0.;
  Double_t cossum = 0.;
  Double_t sumpt = 0.;
  Double_t ntrks = 0.;
  Double_t ntplus = 0.;
  Double_t ortminus = 0.;
  Double_t ortplus = 0.;
  Double_t bdplus = 0.;
  Double_t bdminus = 0.;
  Double_t zmean = 0.;
  Double_t zmeansq = 0.;
  Double_t ww = 0.;
  
  for(unsigned pfj = 0; pfj < fPFJets->GetEntries(); pfj++){
    const PFCandidate* pfca = fPFJets->At(pfj);
    if(! pfca->HasTrk()) continue;
    if(!(pfca->PFType() == PFCandidate::eX ||
	 pfca->PFType() == PFCandidate::eHadron || 
	 pfca->PFType() == PFCandidate::eElectron || 
	 pfca->PFType() == PFCandidate::eMuon) ) continue;
    const Track *t = pfca->Trk();
    if(fabs(t->DzCorrected(*fBeamSpot) - ZVC ) > 0.2) continue;
    
    if(pfca->Pt()<0.3 || pfca->Pt()>200) continue;

    std::vector<const Track*>::iterator itt;
    itt = find ((vtool->excluded).begin(), (vtool->excluded).end(), t);
    if(itt != (vtool->excluded).end()) continue;
    
    zmean = zmean + pfca->Pt()* t->DzCorrected(*fBeamSpot);
    zmeansq = zmeansq + pfca->Pt()*t->DzCorrected(*fBeamSpot)*t->DzCorrected(*fBeamSpot);
    ntrks++;
    sumpt  = sumpt + pfca->Pt()*pfca->Pt();
    ww = ww + pfca->Pt();
    
    sinsum = sinsum + pfca->Pt()*TMath::Sin(t->Phi());
    cossum = cossum + pfca->Pt()*TMath::Cos(t->Phi());
  }
  if(ntrks < 2 || !(sumpt > 0.) ) return 0;
  
  Double_t phim = TMath::ATan2(sinsum,cossum);
  zmean = zmean/ww;
  zmeansq = sqrt(zmeansq/ww - zmean*zmean);
    
  //--------------------
  Double_t bosoproj = 0.;
  Float_t ymean = 0.;
  Float_t ymsq = 0.;
  for(unsigned pfj = 0; pfj < fPFJets->GetEntries(); pfj++){
    const PFCandidate* pfca = fPFJets->At(pfj);
    if(! pfca->HasTrk()) continue;
    if(!(pfca->PFType() == PFCandidate::eX ||
	 pfca->PFType() == PFCandidate::eHadron ||
	 pfca->PFType() == PFCandidate::eElectron ||
	 pfca->PFType() == PFCandidate::eMuon) ) continue;
    
    const Track *t = pfca->Trk();
    if(fabs(t->DzCorrected(*fBeamSpot) - ZVC ) > 0.2) continue;
    //if(fabs(t->DzCorrected(*fBeamSpot) - ZVC ) > 3*zwidth) continue;
    if(pfca->Pt()<0.3 || pfca->Pt()>200) continue;

    std::vector<const Track*>::iterator itt;
    itt = find ((vtool->excluded).begin(), (vtool->excluded).end(), t);
    if(itt != (vtool->excluded).end()) continue;
    
    //Float_t phid = phim - t->Phi();
    Float_t phid = bosophi+3.14 - t->Phi();
    ymean = ymean + pfca->Pt()*TMath::Sin(phid);
    ymsq = ymsq + pow(pfca->Pt()*TMath::Sin(phid),2);
    bosoproj = bosoproj + pfca->Pt()*TMath::Cos(bosophi-t->Phi());
    if(TMath::Sin(phid) > 0.) ortplus = ortplus + pfca->Pt()*TMath::Sin(phid);
    else                      ortminus = ortminus + pfca->Pt()*fabs(TMath::Sin(phid));
    if(TMath::Cos(phid) > 0.) bdplus = bdplus + pfca->Pt()*TMath::Cos(phid);
    else                      bdminus = bdminus + pfca->Pt()*fabs(TMath::Cos(phid));
    if(TMath::Cos(phid) > 0.) ntplus = ntplus + TMath::Cos(phid);
  }
  
  Float_t phimsq = sqrt(ymsq/ntrks - pow(ymean/ntrks,2))*180/3.14;
  // not used... commented out by Fabian
  //Double_t A2n = (ortplus+ortminus) > 0.01 ? (bdplus+bdminus)/(ortplus+ortminus) : (bdplus+bdminus)/0.01;
  Double_t A1n = bdplus/bosopt;
  Double_t A3n = ntplus > 0 ? bdplus/ntplus : 0.;
  //-------------------
  Double_t angle = 180/3.14*TVector2::Phi_0_2pi(bosophi-phim);
  
  vtool->tmvar1 = ntrks;
  vtool->tmvar2 = sumpt;
  vtool->tmvar3 = A1n;
  vtool->tmvar4 = angle;
  vtool->tmvar5 = phimsq;
  vtool->tmvar6 = A3n;
  
  double p1 = vtool->reader->GetProba ( "BDTD method" );
  double p2 = vtool->reader->GetProba ( "BDTG method" );
  double p3 = vtool->reader->GetProba ( "MLP method");
  double p4 = vtool->reader->GetProba ( "MLPBFGS method");
  
  double retval = p1*p2*p3*p4;
  
  return retval;
}

double VertexTools::VertexWidth(const Vertex* vert,  const BaseVertex  *fBeamSpot){
  double width = 0.;
  double zmean = 0.;
  double zmeansq = 0.;
  double ww = 0.;
  if(vert == NULL) return width;
  for(unsigned i = 0; i < vert->NTracks(); i++){
    const Track *t = vert->Trk(i);
    if(t->Pt() < 0.3 ) continue;
    if(t->Pt() > 500.) continue;
    ww = ww + t->Pt();
    zmean = zmean + t->Pt()*t->DzCorrected(*fBeamSpot);
    zmeansq = zmeansq + t->Pt()*t->DzCorrected(*fBeamSpot)*t->DzCorrected(*fBeamSpot);
  }
  if( !(ww > 0.) ) return width;
  zmean = zmean/ww;

  width = sqrt(zmeansq/ww - zmean*zmean);
  return width;
  
}

void VertexTools::BanThisTrack(const Track* track){
  VertexTools* vtool = VertexTools::instance("");
  (vtool->excluded).push_back(track);
}
void VertexTools::Reset(){
  VertexTools* vtool = VertexTools::instance("");
  (vtool->excluded).clear();
}


//------------------------------------------------------------------------------------
// The below tools are from the H->2photons EPS2011 BaseLine (common) Selection
const Vertex* VertexTools::findVtxBasicRanking(const Photon*           ph1, 
					       const Photon*           ph2, 
					       const BaseVertex*       bsp,
					       const VertexCol*        vtcs,
					       const DecayParticleCol* conv, Bool_t useMva, Double_t &vtxProb,
					       std::vector<int>   * debugInds,
					       std::vector<double>* debugVals,
					       std::vector<int>   * debugConv) {
  
  //if (useMva) printf("using mva vertex selection\n");
  
  // check if all input is valid
  if( !ph1 || !ph2 || !bsp || !vtcs ) return NULL;
  // CAUTION: We allow for passing NULL for the Conversions, in that case only the simple Ranking is used.

  // here we will store the idx of the best Vtx
  unsigned int bestIdx = 0;
  UInt_t bestidxmva = 0;
  
  // using asd much as possible 'Globe' naming schemes...
  int*    ptbal_rank    = new int   [vtcs->GetEntries()];
  int*    ptasym_rank   = new int   [vtcs->GetEntries()];
  int*    total_rank    = new int   [vtcs->GetEntries()];
  double* ptbal         = new double[vtcs->GetEntries()];
  double* ptasym        = new double[vtcs->GetEntries()];
  double* sumpt2        = new double[vtcs->GetEntries()];  
  double* limPullToConv = new double[vtcs->GetEntries()];  
  double* mvaval        = new double[vtcs->GetEntries()];  

  
  unsigned int numVertices = vtcs->GetEntries();
  double       ptgg        = 0.;                     // stored for later in the conversion

  // loop over all the vertices...
  for(unsigned int iVtx = 0; iVtx < numVertices; ++iVtx) {
    
    const Vertex* tVtx = vtcs->At(iVtx);
    ptbal      [iVtx] = 0.0;
    ptasym     [iVtx] = 0.0;
    sumpt2     [iVtx] = 0.0; 
    limPullToConv [iVtx] = -1.0;
    ptbal_rank [iVtx] = 1;
    ptasym_rank[iVtx] = 1;
    
    // compute the photon momenta with respect to this Vtx
    FourVectorM newMomFst = ph1->MomVtx(tVtx->Position());
    FourVectorM newMomSec = ph2->MomVtx(tVtx->Position());
    
    FourVectorM higgsMom = newMomFst+newMomSec; 

    double ph1Eta = newMomFst.Eta();
    double ph2Eta = newMomSec.Eta();

    double ph1Phi = newMomFst.Phi();
    double ph2Phi = newMomSec.Phi();
    
    // loop over all tracks and computew variables for ranking...
    FourVectorM totTrkMom(0,0,0,0);
    for(unsigned int iTrk = 0; iTrk < tVtx->NTracks(); ++iTrk) {
      const Track* tTrk = tVtx->Trk(iTrk);
      
      sumpt2[iVtx] += tTrk->Pt()*tTrk->Pt();
      
      // compute distance between Trk and the Photons
      double tEta = tTrk->Eta();
      double tPhi = tTrk->Phi();
      double dEta1 = TMath::Abs(tEta-ph1Eta);
      double dEta2 = TMath::Abs(tEta-ph2Eta);
      double dPhi1 = TMath::Abs(tPhi-ph1Phi);
      double dPhi2 = TMath::Abs(tPhi-ph2Phi);
      if(dPhi1 > M_PI) dPhi1 = 2*M_PI - dPhi1;
      if(dPhi2 > M_PI) dPhi2 = 2*M_PI - dPhi2;

      double dR1 = TMath::Sqrt(dEta1*dEta1+dPhi1*dPhi1);
      double dR2 = TMath::Sqrt(dEta2*dEta2+dPhi2*dPhi2);
            
      if(dR1 < 0.05 || dR2 < 0.05) continue;
      totTrkMom = totTrkMom + tTrk->Mom4(0);
    }

    // compute the ranking variables for this Vtx
    double ptvtx   = totTrkMom.Pt();
    double pthiggs = higgsMom.Pt();
    if(iVtx == 0) ptgg = pthiggs;
    ptbal [iVtx]  = (totTrkMom.Px()*(newMomFst.Px()+newMomSec.Px()));
    ptbal [iVtx] += (totTrkMom.Py()*(newMomFst.Py()+newMomSec.Py()));
    ptbal [iVtx]  = -ptbal[iVtx]/pthiggs;
    ptasym[iVtx]  = (ptvtx - pthiggs)/(ptvtx + pthiggs);
    
    // do the little ranking acrobatics...
    for(unsigned int cVtx =0; cVtx < iVtx; ++cVtx) {
      if(ptbal [iVtx] > ptbal [cVtx])
	ptbal_rank[cVtx]++;
      else
	ptbal_rank[iVtx]++;
      if(ptasym [iVtx] > ptasym [cVtx])
	ptasym_rank[cVtx]++;
      else
	ptasym_rank[iVtx]++;
    }
  }
  
  // loop again over all Vertcices (*sigh*), compute total score and final rank
  // CAUTION: Total rank starts at '0', so the best ranked Vtx has RANK 0
  for(unsigned int iVtx = 0; iVtx < numVertices; ++iVtx) {
    ptasym_rank [iVtx] = ptbal_rank [iVtx]*ptasym_rank [iVtx]*(iVtx+1);
    total_rank  [iVtx] = 0;
    for(unsigned int cVtx =0; cVtx < iVtx; ++cVtx) {
      if(ptasym_rank [iVtx] > ptasym_rank [cVtx]) total_rank[iVtx]++;
      else if(ptasym_rank [iVtx] == ptasym_rank [cVtx]) {                   // use 'ptbal' as the tie-breaker
	if(ptbal_rank [iVtx] > ptbal_rank [cVtx]) total_rank[iVtx]++;
	else total_rank[cVtx]++;
      }
      else total_rank[cVtx]++;
    }
  }
  
  
  // find the best ranked Vertex so far....
  for(unsigned int iVtx = 0; iVtx < numVertices; ++iVtx) {
    if(total_rank[iVtx] == 0) bestIdx = iVtx;
  }
  
  // check if there's a conversion among the pre-selected photons
  // ...this will return NULL in case conv==NULL

  int numLegs1 = -1;
  int numLegs2 = -1;
  int convIdx1 = -1;
  int convIdx2 = -1;
  const DecayParticle* conv1 = PhotonTools::MatchedCiCConversion(ph1, conv, 0.1, 0.1, 0.1, false, &numLegs1, &convIdx1);
  const DecayParticle* conv2 = PhotonTools::MatchedCiCConversion(ph2, conv, 0.1, 0.1, 0.1, false, &numLegs2, &convIdx2);
  

  if ( debugConv ) {
    debugConv->push_back( conv1 ? numLegs1 : 0 );
    debugConv->push_back( conv2 ? numLegs2 : 0 );
    debugConv->push_back( conv1 ? convIdx1 : -1 );
    debugConv->push_back( conv2 ? convIdx2 : -1 );

    debugConv->push_back( conv->GetEntries() );
  } 
  
  double zconv  = 0.;
  double dzconv = 0.;
  int nConv = 0;
  
  if (conv1) nConv += 1;
  if (conv2) nConv += 1;

  double z1, dz1;
  double z2, dz2;
  
  if (conv1) {
    std::pair<double,double> zdz1 = VtxZFromConversion(ph1,conv1,bsp);
    z1 = zdz1.first;
    dz1 = zdz1.second;
    
    zconv = z1;
    dzconv = dz1;

    if( debugVals ) {
      debugVals->push_back(z1);
      debugVals->push_back(dz1);
      debugVals->push_back(conv1->Prob());
    }
  } else {
    if( debugVals ) {
      debugVals->push_back(-99.);
      debugVals->push_back(-99.);
      debugVals->push_back(-99.);
    }
  }

  if (conv2) {
    std::pair<double,double> zdz2 = VtxZFromConversion(ph2,conv2,bsp);
    z2 = zdz2.first;
    dz2 = zdz2.second;
    
    zconv = z2;
    dzconv = dz2;

    if( debugVals ) {
      debugVals->push_back(z2);
      debugVals->push_back(dz2);
      debugVals->push_back(conv2->Prob());
    }
  } else {
    if( debugVals ) {
      debugVals->push_back(-99.);
      debugVals->push_back(-99.);
      debugVals->push_back(-99.);
    }
  }

  if (conv1 && conv2) {
    zconv  = ( 1./(1./dz1/dz1 + 1./dz2/dz2 )*(z1/dz1/dz1 + z2/dz2/dz2) ) ;  // weighted average
    dzconv = TMath::Sqrt( 1./(1./dz1/dz1 + 1./dz2/dz2)) ;
  }
  



  //--------------------------------------------------------------------
  // start doing the Conversion acrobatics... 'copied' from the Globe...
  if(conv1 || conv2) {

    // loop over all ranked Vertices and choose the closest to the Conversion one
    int maxVertices = ( ptgg > 30 ? 3 : 5);
    double minDz = -1.; 
    
    for(unsigned int iVtx =0; iVtx < numVertices; ++iVtx) {

      limPullToConv[iVtx] = TMath::Abs(vtcs->At(iVtx)->Z()-zconv)/dzconv;

      if(total_rank[iVtx] < maxVertices) {
	const Vertex* tVtx = vtcs->At(iVtx);
	double tDz = TMath::Abs(zconv - tVtx->Z());

	if( (minDz < 0. || tDz < minDz) && ( tDz < dzconv ) ) {	  
	  minDz = tDz;
	  bestIdx = iVtx;
	}
      }    
    }
  }
  // END of Conversion Acrobatics
  //--------------------------------------------------------------------
  
  //final loop to compute mva values
  double mvamax = -1e6;
  for(unsigned int iVtx =0; iVtx < numVertices; ++iVtx) {    
    double mva = VtxMvaP(ptbal[iVtx],ptasym[iVtx],log(sumpt2[iVtx]),limPullToConv[iVtx],nConv);
    mvaval[iVtx] = mva;
    if (mva>mvamax) {
      mvamax = mva;
      bestidxmva = iVtx;
    }
    //printf("vtx %i: ptbal = %5f, ptasym = %5f, logsumpt2 = %5f, limpulltoconv = %5f, nconv = %i, mva = %5f\n",iVtx,ptbal[iVtx],ptasym[iVtx],log(sumpt2[iVtx]),limPullToConv[iVtx],nConv,mva);
  }

//   double mvahack = VtxMvaP(4.13519,-0.156296,3.17947,-1.0,0);
//   printf("mvahack = %5f\n",mvahack);
  
  //find second and third ranked vertices for event mva;
  UInt_t mvaidx1 = 0;
  mvamax = -1e6;
  for(unsigned int iVtx =0; iVtx < numVertices; ++iVtx) {    
    if (iVtx!=bestidxmva && mvaval[iVtx]>mvamax) {
      mvamax = mvaval[iVtx];
      mvaidx1 = iVtx;
    }
  }
  
  UInt_t mvaidx2 = 0;
  mvamax = -1e6;
  for(unsigned int iVtx =0; iVtx < numVertices; ++iVtx) {    
    if (iVtx!=bestidxmva && iVtx!=mvaidx1 && mvaval[iVtx]>mvamax) {
      mvamax = mvaval[iVtx];
      mvaidx2 = iVtx;
    }
  }  
  
  //compute per event mva output
  FourVectorM newMomFst = ph1->MomVtx(vtcs->At(bestidxmva)->Position());
  FourVectorM newMomSec = ph2->MomVtx(vtcs->At(bestidxmva)->Position());
  FourVectorM higgsMom = newMomFst+newMomSec;   
  
  fMvaPEvtVars[0] = higgsMom.Pt();
  fMvaPEvtVars[1] = numVertices;
  fMvaPEvtVars[2] = mvaval[bestidxmva];
  fMvaPEvtVars[3] = mvaval[mvaidx1];
  fMvaPEvtVars[4] = vtcs->At(mvaidx1)->Z() - vtcs->At(bestidxmva)->Z();
  fMvaPEvtVars[5] = mvaval[mvaidx2];
  fMvaPEvtVars[6] = vtcs->At(mvaidx2)->Z() - vtcs->At(bestidxmva)->Z();
  fMvaPEvtVars[7] = nConv;
  
  Double_t evtmva = readerevt->EvaluateMVA("BDTEvt");
  vtxProb = 1.-0.49*(evtmva+1.0);
  
//   printf("higgspt = %5f, numvert = %5f, mvabest = %5f, mva1 = %5f, dz1 = %5f, mva2 = %5f, dz2 = %5f, nconv = %5f",fMvaPEvtVars[0],fMvaPEvtVars[1],fMvaPEvtVars[2],fMvaPEvtVars[3],fMvaPEvtVars[4],fMvaPEvtVars[5],fMvaPEvtVars[6],fMvaPEvtVars[7]);
//   printf("vtxmva = %5f, vtxprob = %5f\n",evtmva,vtxProb);
//   
//   printf("e1 = %5f, sige1 = %5f\n",ph1->E(),ph1->EnergyErr());
//   printf("e2 = %5f, sige2 = %5f\n",ph2->E(),ph2->EnergyErr());
  
  // assign the debug value to vectors if present
  if ( debugInds && debugVals ) {
    debugInds->push_back(bestidxmva);
    debugInds->push_back(mvaidx1);		    
    debugInds->push_back(mvaidx2);

    debugVals->push_back(ptbal[bestidxmva]);
    debugVals->push_back(ptasym[bestidxmva]);
    debugVals->push_back(sumpt2[bestidxmva]);
    debugVals->push_back(limPullToConv[bestidxmva]);

    debugVals->push_back(vtcs->At(bestidxmva)->Z());
    debugVals->push_back(vtcs->At(mvaidx1)->Z());
    debugVals->push_back(vtcs->At(mvaidx2)->Z());

    debugVals->push_back(mvaval[bestidxmva]);
    debugVals->push_back(mvaval[mvaidx1]);
    debugVals->push_back(mvaval[mvaidx2]);
    


  }

  // delete the auxiliary dynamic arrays
  delete[] ptbal_rank  ;
  delete[] ptasym_rank ;
  delete[] ptbal       ;
  delete[] ptasym      ;
  delete[] sumpt2      ;
  delete[] limPullToConv;
  delete[] mvaval;

  delete[] total_rank  ;
  
  
  if (useMva) return vtcs->At(bestidxmva);
  else return vtcs->At(bestIdx);
}

//------------------------------------------------------------------------------------
double VertexTools::VtxMvaP(float ptbal, float ptasym, float logsumpt2, float limPullToConv, float nConv) const
{
  fMvaPVars[0] = ptbal;
  fMvaPVars[1] = ptasym;
  fMvaPVars[2] = logsumpt2;
  fMvaPVars[3] = limPullToConv;
  fMvaPVars[4] = nConv;
  
  return readervtx->EvaluateMVA("BDTCat");
  
}

//------------------------------------------------------------------------------------
//Compute contribution to relative uncertainty sigma_m/m from primary vertex location
//given ecal shower positions of two photons, plus the vtx z uncertainty (typically sqrt(2)*beamspot width)
//code originally from Y. Gershtein
Double_t VertexTools::DeltaMassVtx(Double_t xp1, Double_t yp1, Double_t zp1,
            Double_t xp2, Double_t yp2, Double_t zp2,
	    Double_t xv,  Double_t yv,  Double_t zv,
            Double_t dz)
{
  
      Double_t x1 = xp1 - xv;
      Double_t y1 = yp1 - yv;
      Double_t z1 = zp1 - zv;

      Double_t x2 = xp2 - xv;
      Double_t y2 = yp2 - yv;
      Double_t z2 = zp2 - zv;      
      
      Double_t r1 = sqrt(x1*x1+y1*y1+z1*z1);
      Double_t r2 = sqrt(x2*x2+y2*y2+z2*z2);
      Double_t phi1 = atan2(y1,x1);
      Double_t theta1 = atan2(sqrt(x1*x1+y1*y1),z1);
      Double_t phi2 = atan2(y2,x2);
      Double_t theta2 = atan2(sqrt(x2*x2+y2*y2),z2);

      Double_t sech1 = sin(theta1);
      Double_t tanh1 = cos(theta1);
      Double_t sech2 = sin(theta2);
      Double_t tanh2 = cos(theta2);
      Double_t cos12 = cos(phi1-phi2);

      Double_t rad1 = sech1*(sech1*tanh2-tanh1*sech2*cos12)/(1-tanh1*tanh2-sech1*sech2*cos12);
      Double_t rad2 = sech2*(sech2*tanh1-tanh2*sech1*cos12)/(1-tanh2*tanh1-sech2*sech1*cos12);

      return dz * 0.5*fabs(rad1/r1 + rad2/r2);
      
}

std::pair<double,double> VertexTools::VtxZFromConversion(const Photon *p, const DecayParticle *c, const BaseVertex *bsp) {
  
  //2012 ABC
  //const double dzpxb = 0.011;
  //const double dztib = 0.301;
  //const double dztob = 1.603;
  //const double dzpxf = 0.057;
  //const double dztid = 0.316;
  //const double dztec = 1.026;

  //2012 ABCD
  const double dzpxb = 0.011;
  const double dztib = 0.297;
  const double dztob = 1.530;
  const double dzpxf = 0.059;
  const double dztid = 0.315;
  const double dztec = 1.021;
 
  const double dzpxbsingle = 0.007;
  const double dztibsingle = 0.472;
  const double dztobsingle = 0.540;
  const double dzpxfsingle = 0.039;
  const double dztidsingle = 0.349;
  const double dztecsingle = 1.022;  
  
  
  double zconv = -99.;
  double dzconv = -99.;
  
  const mithep::ThreeVector caloPos(p->SCluster()->Point());
  double zconvsc   = c->Z0EcalVtxCiC(bsp->Position(), caloPos);
  double zconvtrk  = c->DzCorrected(bsp->Position()) + bsp->Z();
  
  if (c->NDaughters()==2) {
    if( p->SCluster()->AbsEta()<1.5 ) {
      double rho = c->Position().Rho();
      if     ( rho < 15. ) { dzconv = dzpxb; zconv = zconvtrk; }
      else if( rho < 60. ) { dzconv = dztib; zconv = zconvsc; }
      else                 { dzconv = dztob; zconv = zconvsc; }
    } else {
      double z = c->Position().Z();
      if     ( TMath::Abs(z) < 50. )   { dzconv = dzpxf; zconv = zconvtrk; }
      else if( TMath::Abs(z) < 100.)   { dzconv = dztid; zconv = zconvtrk; }
      else                             { dzconv = dztec; zconv = zconvsc;  }
    }  
  }
  else if (c->NDaughters()==1) {
    if( p->SCluster()->AbsEta()<1.5 ) {
      double rho = c->Position().Rho();
      if     ( rho < 15. ) { dzconv = dzpxbsingle; zconv = zconvtrk; }
      else if( rho < 60. ) { dzconv = dztibsingle; zconv = zconvsc; }
      else                 { dzconv = dztobsingle; zconv = zconvsc; }
    } else {
      double z = c->Position().Z();
      if     ( TMath::Abs(z) < 50. )   { dzconv = dzpxfsingle; zconv = zconvtrk; }
      else if( TMath::Abs(z) < 100.)   { dzconv = dztidsingle; zconv = zconvtrk; }
      else                             { dzconv = dztecsingle; zconv = zconvsc;  }
    }  
  }
  
  return std::pair<double,double>(zconv,dzconv);
  
  
}
