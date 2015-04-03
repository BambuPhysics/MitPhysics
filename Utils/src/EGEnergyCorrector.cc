// $Id: EGEnergyCorrector.cc,v 1.14 2013/11/13 20:51:33 bendavid Exp $

#include "MitPhysics/Utils/interface/EGEnergyCorrector.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooAbsPdf.h"
#include "RooConstVar.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooHybridBDTAutoPdf.h"
#include "HiggsAnalysis/GBRLikelihood/interface/RooDoubleCBFast.h"
#include "HiggsAnalysis/GBRLikelihood/interface/HybridGBRForest.h"
#include "HiggsAnalysis/GBRLikelihood/interface/HybridGBRForestD.h"
// #include "Cintex/Cintex.h"
#include <TFile.h>
#include <TRandom3.h>

ClassImp(mithep::EGEnergyCorrector)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
EGEnergyCorrector::EGEnergyCorrector() :
fReadereb(0),
fReaderebvariance(0),
fReaderee(0),
fReadereevariance(0),
fReaderebsemi(0),
fReadereesemi(0),
fReaderDebsemi(0),
fReaderDeesemi(0),
_mean(0),
_tgt(0),
_sigma(0),
_n1(0),
_n2(0),
_meanlim(0),
_sigmalim(0),
_n1lim(0),
_n2lim(0),
_pdf(0),
fMethodname("BDTG method"),
fIsInitialized(kFALSE),
fVals(0)
{
  // Constructor.
}


//--------------------------------------------------------------------------------------------------
EGEnergyCorrector::~EGEnergyCorrector()
{
  
  if (fVals) delete [] fVals;
  if (fReadereb) delete fReadereb;
  if (fReaderebvariance) delete fReaderebvariance;  
  if (fReaderee) delete fReaderee;
  if (fReadereevariance) delete fReadereevariance;
  if (fReaderebsemi) delete fReaderebsemi;
  if (fReadereesemi) delete fReadereesemi;
  if (fReaderDebsemi) delete fReaderDebsemi;
  if (fReaderDeesemi) delete fReaderDeesemi;  
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrector::Initialize(TString phfixstring, TString phfixfile, TString regweights, Int_t version=0) {
    if (fVals) delete [] fVals;
    if (fReadereb) delete fReadereb;
    if (fReaderebvariance) delete fReaderebvariance;  
    if (fReaderee) delete fReaderee;
    if (fReadereevariance) delete fReadereevariance;        
    if (fReaderebsemi) delete fReaderebsemi;
    if (fReadereesemi) delete fReadereesemi;
    if (fReaderDebsemi) delete fReaderDebsemi;
    if (fReaderDeesemi) delete fReaderDeesemi;
    
    if (version<=3) {
    
      fPhFix.initialise(std::string(phfixstring),std::string(phfixfile));


      
      fVals = new Float_t[73];

      // Is this necessary? (Y.I. Apr 2 2015)
      //      ROOT::Cintex::Cintex::Enable();   
      
      TFile *fgbr = new TFile(regweights,"READ");
      fReadereb = (GBRForest*)fgbr->Get("EBCorrection");
      fReaderebvariance = (GBRForest*)fgbr->Get("EBUncertainty");  
      fReaderee = (GBRForest*)fgbr->Get("EECorrection");
      fReadereevariance = (GBRForest*)fgbr->Get("EEUncertainty");      
      fgbr->Close();
      
      fIsInitialized = kTRUE;

    
    }
    else if (version==5) {
      fVals = new Float_t[37];
      
      TFile *fgbr = TFile::Open(regweights,"READ");
      fgbr->GetObject("EGRegressionForest_EB", fReaderebsemi);
      fgbr->GetObject("EGRegressionForest_EE", fReadereesemi);
      fgbr->Close();

      //recreate pdf with constraint transformations (can't load directly from file due to weird RooWorkspace IO features)
      
      _tgt = new RooRealVar("tgt","",1.);
      _mean = new RooRealVar("mean","",1.);
      _sigma = new RooRealVar("sigma","",1.);
      _n1 = new RooRealVar("n1","",2.);
      _n2 = new RooRealVar("n2","",2.);
      
      _sigmalim = new RooRealConstraint("sigmalim","",*_sigma,0.0002,0.5);
      _meanlim = new RooRealConstraint("meanlim","",*_mean,0.2,2.0);
      _n1lim = new RooRealConstraint("n1lim","",*_n1,1.01,110.);
      _n2lim = new RooRealConstraint("n2lim","",*_n2,1.01,110.);
      
      RooConstVar *cbmean = new RooConstVar("cbmean","",1.0);
      RooConstVar *alpha1 = new RooConstVar("alpha1","",2.0);
      RooConstVar *alpha2 = new RooConstVar("alpha2","",1.0);
      
      _pdf = new RooDoubleCBFast("sigpdf","",*_tgt,*cbmean,*_sigmalim,*alpha1,*_n1lim,*alpha2,*_n2lim);
      
      //add to RooArgList for proper garbage collection
      _args.addOwned(*_tgt);
      _args.addOwned(*_mean);
      _args.addOwned(*_sigma);
      _args.addOwned(*_n1);
      _args.addOwned(*_n2);
      _args.addOwned(*cbmean);
      _args.addOwned(*alpha1);
      _args.addOwned(*alpha2);
      _args.addOwned(*_sigmalim);
      _args.addOwned(*_meanlim);
      _args.addOwned(*_n1lim);
      _args.addOwned(*_n2lim);
      _args.addOwned(*_pdf);       
      
      //add target var to static normalization set to avoid memory churn/leak
      _normset.add(*_tgt);
      
      fIsInitialized = kTRUE;     
      
    }
    else if (version>=6 && version<=8) {
      fVals = new Float_t[37];
      
      TFile *fgbr = TFile::Open(regweights,"READ");
      fgbr->GetObject("EGRegressionForest_EB", fReaderDebsemi);
      fgbr->GetObject("EGRegressionForest_EE", fReaderDeesemi);
      fgbr->Close();      
      
      assert(fReaderDebsemi!=0);
      assert(fReaderDeesemi!=0);
      
      _tgt = new RooRealVar("tgt","",1.);
      _mean = new RooRealVar("mean","",1.);
      _sigma = new RooRealVar("sigma","",0.01);
      _n1 = new RooRealVar("n1","",2.);
      _n2 = new RooRealVar("n2","",2.);
      
      _sigmalim = new RooRealConstraint("sigmalim","",*_sigma,0.0002,0.5);
      _meanlim = new RooRealConstraint("meanlim","",*_mean,0.2,2.0);
      _n1lim = new RooRealConstraint("n1lim","",*_n1,1.01,5000.);
      _n2lim = new RooRealConstraint("n2lim","",*_n2,1.01,5000.);
      
      RooConstVar *alpha1 = new RooConstVar("alpha1","",2.0);
      RooConstVar *alpha2 = new RooConstVar("alpha2","",1.0);
      
      _pdf = new RooDoubleCBFast("sigpdf","",*_tgt,*_meanlim,*_sigmalim,*alpha1,*_n1lim,*alpha2,*_n2lim);
      
      //add to RooArgList for proper garbage collection
      _args.addOwned(*_tgt);
      _args.addOwned(*_mean);
      _args.addOwned(*_sigma);
      _args.addOwned(*_n1);
      _args.addOwned(*_n2);
      _args.addOwned(*alpha1);
      _args.addOwned(*alpha2);
      _args.addOwned(*_sigmalim);
      _args.addOwned(*_meanlim);
      _args.addOwned(*_n1lim);
      _args.addOwned(*_n2lim);
      _args.addOwned(*_pdf);      
      
      //add target var to static normalization set to avoid memory churn/leak
      _normset.add(*_tgt);      
      
      fIsInitialized = kTRUE;
      
    }

}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrector::CorrectEnergyWithError(Photon *p, const VertexCol *vtxs, Double_t rho, UInt_t version, Bool_t applyRescale) {
  
  
  if (version<=3) {
    std::pair<double,double> correction;
    if (version == 1) {
      correction = CorrectedEnergyWithError(p);
    }
    else if (version == 2) {
      correction = CorrectedEnergyWithErrorV2(p,vtxs);
    }
    else if (version == 3) {
      correction = CorrectedEnergyWithErrorV3(p,vtxs,rho,applyRescale);
    }

    //printf("photon: e = %5f, eerr = %5f, sceta = %5f\n",p->E(),p->EnergyErr(),p->SCluster()->Eta());
    //printf("gbr   : e = %5f, eerr = %5f\n",correction.first,correction.second);
    FourVectorM mom = p->Mom();
    double scale = correction.first/mom.E();
    p->SetMom(scale*mom.X(), scale*mom.Y(), scale*mom.Z(), scale*mom.E());
    p->SetEnergyErr(correction.second);
    p->SetEnergyErrSmeared(correction.second);
  }
  else {
    
    if (version == 5) {
      double ecor,mean,sigma,alpha1,n1,alpha2,n2,peakpdfval;
      CorrectedEnergyWithErrorV5(p,vtxs,rho,ecor,mean,sigma,alpha1,n1,alpha2,n2,peakpdfval);
      
      FourVectorM mom = p->Mom();
      double scale = ecor/mom.E();
      p->SetMom(scale*mom.X(), scale*mom.Y(), scale*mom.Z(), scale*mom.E());
      p->SetEnergyErr(sigma);
      p->SetEnergyErrSmeared(sigma);    
    }
    if (version>=6 && version<=8) {
      double ecor,sigEoverE,mean,sigma,alpha1,n1,alpha2,n2,peakpdfval;
      
      if (version==6) {
        CorrectedEnergyWithErrorV6(p,vtxs,rho,ecor,sigEoverE, mean,sigma,alpha1,n1,alpha2,n2,peakpdfval);
      }
      else if (version==7) {
        CorrectedEnergyWithErrorV7(p,vtxs,rho,ecor,sigEoverE, mean,sigma,alpha1,n1,alpha2,n2,peakpdfval);  
      }
      else if (version==8) {
        CorrectedEnergyWithErrorV8(p,vtxs,rho,ecor,sigEoverE, mean,sigma,alpha1,n1,alpha2,n2,peakpdfval);        
      }
      
      FourVectorM mom = p->Mom();
      double scale = ecor/mom.E();
      p->SetMom(scale*mom.X(), scale*mom.Y(), scale*mom.Z(), scale*mom.E());
      p->SetEnergyErr(sigEoverE);
      p->SetEnergyErrSmeared(sigEoverE);    
    }    
    
  }
  
  
  //printf("initial pt = %5f, regression pt = %5f\n",mom.Pt(),p->Pt());
  
  return;

}

//--------------------------------------------------------------------------------------------------
std::pair<double,double> EGEnergyCorrector::CorrectedEnergyWithError(const Photon *p) {
  
  const SuperCluster *s = p->SCluster();
 // const BasicCluster *b = s->Seed();  
  fPhFix.setup(p->E(),s->Eta(),s->Phi(),p->R9()); 
  
  double scpse = s->PreshowerEnergy();
  
  Bool_t isbarrel = (s->AbsEta()<1.5);

  if (isbarrel) {
    fVals[0]  = s->RawEnergy();
    fVals[1]  = p->R9();
    fVals[2]  = s->Eta();
    fVals[3]  = s->Phi();
    fVals[4]  = p->E55()/s->RawEnergy();
    fVals[5]  = fPhFix.etaC();
    fVals[6]  = fPhFix.etaS();
    fVals[7]  = fPhFix.etaM();
    fVals[8]  = fPhFix.phiC();
    fVals[9]  = fPhFix.phiS();
    fVals[10] = fPhFix.phiM();    
    fVals[11] = p->HadOverEm();
    fVals[12] = s->EtaWidth();
    fVals[13] = s->PhiWidth();
    fVals[14] = p->CoviEtaiEta();
  }
  else {
    fVals[0]  = s->RawEnergy();
    fVals[1]  = p->R9();
    fVals[2]  = s->Eta();
    fVals[3]  = s->Phi();
    fVals[4]  = p->E55()/s->RawEnergy();
    fVals[5]  = scpse/s->RawEnergy();
    fVals[6]  = fPhFix.xZ();
    fVals[7]  = fPhFix.xC();
    fVals[8]  = fPhFix.xS();
    fVals[9]  = fPhFix.xM();
    fVals[10] = fPhFix.yZ();
    fVals[11] = fPhFix.yC();
    fVals[12] = fPhFix.yS();
    fVals[13] = fPhFix.yM();
    fVals[14] = p->HadOverEm();
    fVals[15] = s->EtaWidth();
    fVals[16] = s->PhiWidth();
    fVals[17] = p->CoviEtaiEta();    
  }
    
  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *reader;
  const GBRForest *readervar;
  if (isbarrel) {
    den = s->RawEnergy();
    reader = fReadereb;
    readervar = fReaderebvariance;
  }
  else {
    den = s->RawEnergy() + scpse;
    reader = fReaderee;
    readervar = fReadereevariance;
  }
  
  Double_t ecor = reader->GetResponse(fVals)*den;
  Double_t ecorerr = readervar->GetResponse(fVals)*den*varscale;
  
  return std::pair<double,double>(ecor,ecorerr);
}

//--------------------------------------------------------------------------------------------------
std::pair<double,double> EGEnergyCorrector::CorrectedEnergyWithErrorV2(const Photon *p, const VertexCol *vtxs) {
  
  const SuperCluster *s = p->SCluster();
  const BasicCluster *b = s->Seed();
  
  const BasicCluster *b2 = 0;
  Double_t ebcmax = -99.;
  for (UInt_t i=0; i<s->ClusterSize(); ++i) {
    const BasicCluster *bc = s->Cluster(i);
    if (bc->Energy() > ebcmax && bc !=b) {
      b2 = bc;
      ebcmax = bc->Energy();
    }
  }

  const BasicCluster *bclast = 0;
  Double_t ebcmin = 1e6;
  for (UInt_t i=0; i<s->ClusterSize(); ++i) {
    const BasicCluster *bc = s->Cluster(i);
    if (bc->Energy() < ebcmin && bc !=b) {
      bclast = bc;
      ebcmin = bc->Energy();
    }
  }

  const BasicCluster *bclast2 = 0;
  ebcmin = 1e6;
  for (UInt_t i=0; i<s->ClusterSize(); ++i) {
    const BasicCluster *bc = s->Cluster(i);
    if (bc->Energy() < ebcmin && bc !=b && bc!=bclast) {
      bclast2 = bc;
      ebcmin = bc->Energy();
    }
  }

  
  double scpse = s->PreshowerEnergy();
  
  Bool_t isbarrel = (s->AbsEta()<1.5);
  Bool_t hasbc2 = b2 && b2->Energy()>0.;
  Bool_t hasbclast = bclast && bclast->Energy()>0.;
  Bool_t hasbclast2 = bclast2 && bclast2->Energy()>0.;
  
  
  if (isbarrel) {
    fVals[0]  = s->RawEnergy();
    fVals[1]  = p->R9();
    fVals[2]  = s->Eta();
    fVals[3]  = s->Phi();
    fVals[4]  = p->E55()/s->RawEnergy();   
    fVals[5] = p->HadOverEm();
    fVals[6] = s->EtaWidth();
    fVals[7] = s->PhiWidth();
    
    fVals[8] = b->Eta()-s->Eta();
    fVals[9] = MathUtils::DeltaPhi(b,s);
    fVals[10] = b->Energy()/s->RawEnergy();
    fVals[11] = b->E3x3()/b->Energy();
    fVals[12] = b->E5x5()/b->Energy();
    fVals[13] = TMath::Sqrt(b->CoviEtaiEta());
    fVals[14] = TMath::Sqrt(b->CoviPhiiPhi());
    fVals[15] = b->CoviEtaiPhi();
    fVals[16] = b->EMax()/b->Energy();
    fVals[17] = log(b->E2nd()/b->EMax());
    fVals[18] = log(b->ETop()/b->EMax());
    fVals[19] = log(b->EBottom()/b->EMax());
    fVals[20] = log(b->ELeft()/b->EMax());
    fVals[21] = log(b->ERight()/b->EMax());
    fVals[22] = (b->ETop()-b->EBottom())/(b->ETop()+b->EBottom());
    fVals[23] = (b->ELeft()-b->ERight())/(b->ELeft()+b->ERight());

    fVals[24] = hasbc2 ? (b2->Eta()-s->Eta()) : 0.;
    fVals[25] = hasbc2 ? MathUtils::DeltaPhi(b2,s) : 0.;
    fVals[26] = hasbc2 ? b2->Energy()/s->RawEnergy() : 0.;
    fVals[27] = hasbc2 ? b2->E3x3()/b2->Energy() : 0.;
    fVals[28] = hasbc2 ? b2->E5x5()/b2->Energy() : 0.;
    fVals[29] = hasbc2 ? TMath::Sqrt(b2->CoviEtaiEta()) : 0.;
    fVals[30] = hasbc2 ? TMath::Sqrt(b2->CoviPhiiPhi()) : 0.;
    //fVals[31] = hasbc2 ? b2->CoviEtaiPhi() : 0.;
    fVals[31] = hasbc2 ? b->CoviEtaiPhi() : 0.;
    fVals[32] = hasbc2 ? b2->EMax()/b2->Energy() : 0.;
    fVals[33] = hasbc2 ? log(b2->E2nd()/b2->EMax()) : 0.;
    fVals[34] = hasbc2 ? log(b2->ETop()/b2->EMax()) : 0.;
    fVals[35] = hasbc2 ? log(b2->EBottom()/b2->EMax()) : 0.;
    fVals[36] = hasbc2 ? log(b2->ELeft()/b2->EMax()) : 0.;
    fVals[37] = hasbc2 ? log(b2->ERight()/b2->EMax()) : 0.;
    fVals[38] = hasbc2 ? (b2->ETop()-b2->EBottom())/(b2->ETop()+b2->EBottom()) : 0.;
    fVals[39] = hasbc2 ? (b2->ELeft()-b2->ERight())/(b2->ELeft()+b2->ERight()) : 0.;

    fVals[40] = hasbclast ? (bclast->Eta()-s->Eta()) : 0.;
    fVals[41] = hasbclast ? MathUtils::DeltaPhi(bclast,s) : 0.;
    fVals[42] = hasbclast ? bclast->Energy()/s->RawEnergy() : 0.;
    fVals[43] = hasbclast ? bclast->E3x3()/bclast->Energy() : 0.;
    fVals[44] = hasbclast ? bclast->E5x5()/bclast->Energy() : 0.;
    fVals[45] = hasbclast ? TMath::Sqrt(bclast->CoviEtaiEta()) : 0.;
    fVals[46] = hasbclast ? TMath::Sqrt(bclast->CoviPhiiPhi()) : 0.;
    fVals[47] = hasbclast ? bclast->CoviEtaiPhi() : 0.;    

    fVals[48] = hasbclast2 ? (bclast2->Eta()-s->Eta()) : 0.;
    fVals[49] = hasbclast2 ? MathUtils::DeltaPhi(bclast2,s) : 0.;
    fVals[50] = hasbclast2 ? bclast2->Energy()/s->RawEnergy() : 0.;
    fVals[51] = hasbclast2 ? bclast2->E3x3()/bclast2->Energy() : 0.;
    fVals[52] = hasbclast2 ? bclast2->E5x5()/bclast2->Energy() : 0.;
    fVals[53] = hasbclast2 ? TMath::Sqrt(bclast2->CoviEtaiEta()) : 0.;
    fVals[54] = hasbclast2 ? TMath::Sqrt(bclast2->CoviPhiiPhi()) : 0.;
    fVals[55] = hasbclast2 ? bclast2->CoviEtaiPhi() : 0.;        
    
    fVals[56] = b->IEta();
    fVals[57] = b->IPhi();
    fVals[58] = b->IEta()%5;
    fVals[59] = b->IPhi()%2;
    fVals[60] = (TMath::Abs(b->IEta())<=25)*(b->IEta()%25) + (TMath::Abs(b->IEta())>25)*((b->IEta()-25*TMath::Abs(b->IEta())/b->IEta())%20);
    fVals[61] = b->IPhi()%20;
    fVals[62] = b->EtaCry();
    fVals[63] = b->PhiCry();

    fVals[64] = hasbc2 ? b2->IEta() : 0.;
    fVals[65] = hasbc2 ? b2->IPhi() : 0.;
    fVals[66] = hasbc2 ? b2->IEta()%5 : 0.;
    fVals[67] = hasbc2 ? b2->IPhi()%2 : 0.;
    fVals[68] = hasbc2 ? (TMath::Abs(b2->IEta())<=25)*(b2->IEta()%25) + (TMath::Abs(b2->IEta())>25)*((b2->IEta()-25*TMath::Abs(b2->IEta())/b2->IEta())%20) : 0.;
    fVals[69] = hasbc2 ? b2->IPhi()%20 : 0.;
    fVals[70] = hasbc2 ? b2->EtaCry() : 0.;
    fVals[71] = hasbc2 ? b2->PhiCry() : 0.;    
    
    fVals[72] = vtxs->GetEntries();
    
  }
  else {
    fVals[0]  = s->RawEnergy();
    fVals[1]  = p->R9();
    fVals[2]  = s->Eta();
    fVals[3]  = s->Phi();
    fVals[4]  = p->E55()/s->RawEnergy();
    fVals[5] = s->EtaWidth();
    fVals[6] = s->PhiWidth();
    fVals[7] = vtxs->GetEntries();
  }
    
  const Double_t varscale = 1.253;
  Double_t den;
  const GBRForest *reader;
  const GBRForest *readervar;
  if (isbarrel) {
    den = s->RawEnergy();
    reader = fReadereb;
    readervar = fReaderebvariance;
  }
  else {
    den = s->RawEnergy() + scpse;
    reader = fReaderee;
    readervar = fReadereevariance;
  }
  
  Double_t ecor = reader->GetResponse(fVals)*den;
  Double_t ecorerr = readervar->GetResponse(fVals)*den*varscale;
  
  return std::pair<double,double>(ecor,ecorerr);
}

//--------------------------------------------------------------------------------------------------
std::pair<double,double> EGEnergyCorrector::CorrectedEnergyWithErrorV3(const Photon *p, const VertexCol *vtxs, Double_t rho, Bool_t applyRescale) {
  
  const SuperCluster *s = p->SCluster();
  const BasicCluster *b = s->Seed();
  
  
  Bool_t isbarrel = (s->AbsEta()<1.5);
  
  //basic supercluster variables
  fVals[0]  = s->RawEnergy();
  fVals[1]  = s->Eta();
  fVals[2]  = s->Phi();
  fVals[3]  = p->R9();
  fVals[4]  = p->E55()/s->RawEnergy();  
  fVals[5] = s->EtaWidth();
  fVals[6] = s->PhiWidth();
  fVals[7] = s->NClusters();
  fVals[8] = p->HadOverEmTow();
  fVals[9] = rho;
  fVals[10] = vtxs->GetEntries();
  fVals[11] = b->Eta()-s->Eta();
  fVals[12] = atan2(sin(b->Phi()-s->Phi()), cos(b->Phi()-s->Phi()));
  fVals[13] = b->Energy()/s->RawEnergy();
  fVals[14] = b->E3x3()/b->Energy();
  fVals[15] = b->E5x5()/b->Energy();
  fVals[16] = TMath::Sqrt(b->CoviEtaiEta());
  fVals[17] = TMath::Sqrt(b->CoviPhiiPhi());
  fVals[18] = b->CoviEtaiPhi();  
  fVals[19] = b->EMax()/b->Energy();                       //crystal energy ratio gap variables   
  fVals[20] = b->E2nd()/b->Energy();
  fVals[21] = b->ETop()/b->Energy();
  fVals[22] = b->EBottom()/b->Energy();
  fVals[23] = b->ELeft()/b->Energy();
  fVals[24] = b->ERight()/b->Energy();
  fVals[25] = b->E2x5Max()/b->Energy();                       //crystal energy ratio gap variables   
  fVals[26] = b->E2x5Top()/b->Energy();
  fVals[27] = b->E2x5Bottom()/b->Energy();
  fVals[28] = b->E2x5Left()/b->Energy();
  fVals[29] = b->E2x5Right()/b->Energy();  
  
  if (isbarrel) {
    fVals[30] = b->IEta();
    fVals[31] = b->IPhi();
    fVals[32] = b->IEta()%5;
    fVals[33] = b->IPhi()%2;
    fVals[34] = (TMath::Abs(b->IEta())<=25)*(b->IEta()%25) + (TMath::Abs(b->IEta())>25)*((b->IEta()-25*TMath::Abs(b->IEta())/b->IEta())%20);
    fVals[35] = b->IPhi()%20;
    fVals[36] = b->EtaCry();
    fVals[37] = b->PhiCry();    
  }
  else {
    //preshower energy ratio (endcap only)
    fVals[30]  = s->PreshowerEnergy()/s->RawEnergy();
  }
    
  Double_t den;
  const GBRForest *reader;
  const GBRForest *readervar;
  if (isbarrel) {
    den = s->RawEnergy();
    reader = fReadereb;
    readervar = fReaderebvariance;
  }
  else {
    den = s->RawEnergy() + s->PreshowerEnergy();
    reader = fReaderee;
    readervar = fReadereevariance;
  }
  
  Double_t ecor = reader->GetResponse(fVals)*den;
  
  //apply shower shape rescaling - for Monte Carlo only, and only for calculation of energy uncertainty
  if (applyRescale) {
    if (isbarrel) {
      fVals[3] = 1.0045*p->R9() +0.001; //r9
      fVals[5] = 1.04302*s->EtaWidth() - 0.000618; //etawidth
      fVals[6] = 1.00002*s->PhiWidth() - 0.000371;  //phiwidth
      fVals[14] = fVals[3]*s->RawEnergy()/b->Energy();  //compute consistent e3x3/eseed after r9 rescaling
      if (fVals[15]<=1.0)  // rescale e5x5/eseed only if value is <=1.0, don't allow scaled values to exceed 1.0
        fVals[15] = TMath::Min(1.0,1.0022*b->E5x5()/b->Energy());

      fVals[4] = fVals[15]*b->Energy()/s->RawEnergy(); // compute consistent e5x5()/rawEnergy() after e5x5/eseed resacling  

      fVals[16] = 0.891832*TMath::Sqrt(b->CoviEtaiEta()) + 0.0009133; //sigietaieta
      fVals[17] = 0.993*TMath::Sqrt(b->CoviPhiiPhi());; //sigiphiiphi

      fVals[19] = 1.012*b->EMax()/b->Energy();                       //crystal energy ratio gap variables   
      fVals[20] = 1.0*b->E2nd()/b->Energy();
      fVals[21] = 0.94*b->ETop()/b->Energy();
      fVals[22] = 0.94*b->EBottom()/b->Energy();
      fVals[23] = 0.94*b->ELeft()/b->Energy();
      fVals[24] = 0.94*b->ERight()/b->Energy();
      fVals[25] = 1.006*b->E2x5Max()/b->Energy();                       //crystal energy ratio gap variables   
      fVals[26] = 1.09*b->E2x5Top()/b->Energy();
      fVals[27] = 1.09*b->E2x5Bottom()/b->Energy();
      fVals[28] = 1.09*b->E2x5Left()/b->Energy();
      fVals[29] = 1.09*b->E2x5Right()/b->Energy();

    }
    else {
      fVals[3] = 1.0086*p->R9() -0.0007;           //r9
      fVals[4] = TMath::Min(1.0,1.0022*b->E5x5()/s->RawEnergy());  //e5x5/rawenergy
      fVals[5] = 0.903254*s->EtaWidth() +0.001346;  //etawidth
      fVals[6] = 0.99992*s->PhiWidth() +  4.8e-07;   //phiwidth
      fVals[13] = TMath::Min(1.0,1.0022*b->Energy()/s->RawEnergy());  //eseed/rawenergy (practically equivalent to e5x5)


      fVals[14] = fVals[3]*s->RawEnergy()/b->Energy(); //compute consistent e3x3/eseed after r9 rescaling

      fVals[16] = 0.9947*TMath::Sqrt(b->CoviEtaiEta()) + 0.00003; //sigietaieta

      fVals[19] = 1.005*b->EMax()/b->Energy();                  //crystal energy ratio gap variables   
      fVals[20] = 1.02*b->E2nd()/b->Energy();
      fVals[21] = 0.96*b->ETop()/b->Energy();
      fVals[22] = 0.96*b->EBottom()/b->Energy();
      fVals[23] = 0.96*b->ELeft()/b->Energy();
      fVals[24] = 0.96*b->ERight()/b->Energy();
      fVals[25] = 1.0075*b->E2x5Max()/b->Energy();                          //crystal energy ratio gap variables   
      fVals[26] = 1.13*b->E2x5Top()/b->Energy();
      fVals[27] = 1.13*b->E2x5Bottom()/b->Energy();
      fVals[28] = 1.13*b->E2x5Left()/b->Energy();
      fVals[29] = 1.13*b->E2x5Right()/b->Energy();
    }

  }  
  
  Double_t ecorerr = readervar->GetResponse(fVals)*den;
  
  return std::pair<double,double>(ecor,ecorerr);
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrector::CorrectedEnergyWithErrorV5(const Photon *p, const VertexCol *vtxs, Double_t rho, Double_t &ecor, Double_t &mean, Double_t &sigma, Double_t &alpha1, Double_t &n1, Double_t &alpha2, Double_t &n2, Double_t &pdfpeakval) {
 
  const SuperCluster *s = p->SCluster();
  const BasicCluster *b = s->Seed();
  
  
  Bool_t isbarrel = (s->AbsEta()<1.5);
  
// //basic supercluster variables
  fVals[0] = s->RawEnergy();
  fVals[1] = s->Eta();
  fVals[2] = p->R9();
  fVals[3] = s->EtaWidth();
  fVals[4] = s->PhiWidth();
  fVals[5] = double(s->NClusters());
  fVals[6] = p->HadOverEmTow();
  fVals[7] = rho;
  fVals[8] = double(vtxs->GetEntries());

  //seed basic cluster variables
  fVals[9] = b->Eta()-s->Eta();
  fVals[10] = atan2(sin(b->Phi()-s->Phi()), cos(b->Phi()-s->Phi()));
  fVals[11] = b->Energy()/s->RawEnergy();
  fVals[12] = b->E3x3()/b->E5x5();
  fVals[13] = sqrt(b->CoviEtaiEta()); //sigietaieta
  fVals[14] = sqrt(b->CoviPhiiPhi()); //sigiphiiphi
  fVals[15] = b->CoviEtaiPhi();   //sigietaiphi
  fVals[16] = b->EMax()/b->E5x5(); //crystal energy ratio gap variables
  fVals[17] = b->E2nd()/b->E5x5();
  fVals[18] = b->ETop()/b->E5x5();
  fVals[19] = b->EBottom()/b->E5x5();
  fVals[20] = b->ELeft()/b->E5x5();
  fVals[21] = b->ERight()/b->E5x5();
  fVals[22] = b->E2x5Max()/b->E5x5(); //crystal energy ratio gap variables
  fVals[23] = b->E2x5Top()/b->E5x5();
  fVals[24] = b->E2x5Bottom()/b->E5x5();
  fVals[25] = b->E2x5Left()/b->E5x5();
  fVals[26] = b->E2x5Right()/b->E5x5();

  if (isbarrel) {
    //additional energy ratio (always ~1 for endcap, therefore only included for barrel)
    fVals[27] = b->E5x5()/b->Energy();
    
    fVals[28] = b->IEta();
    fVals[29] = b->IPhi()%18;  //should really be (iphi-1)%20 ie duplicate with the later variable
    fVals[30] = b->IEta()%5;
    fVals[31] = b->IPhi()%2;
    fVals[32] = (TMath::Abs(b->IEta())<=25)*(b->IEta()%25) + (TMath::Abs(b->IEta())>25)*((b->IEta()-25*TMath::Abs(b->IEta())/b->IEta())%20); //should be ieta->ieta-1 for %
    fVals[33] = b->IPhi()%20; //should really be (iphi-1)%20
    fVals[34] = b->EtaCry();
    fVals[35] = b->PhiCry();        
    

  }
  else {
    //preshower energy ratio (endcap only)
    fVals[27] = s->PreshowerEnergy()/s->RawEnergy();
  }
  
  double den;
  HybridGBRForest *forest;
  if (isbarrel) {
    den = s->RawEnergy();
    forest = fReaderebsemi;
  }
  else {
    den = s->RawEnergy() + s->PreshowerEnergy();
    forest = fReadereesemi;
  }
  
  _tgt->setVal(1.0); //evaluate pdf at peak position
  
  //set raw response variables from GBRForest
  _sigma->setVal(forest->GetResponse(&fVals[0],0));
  _mean->setVal(forest->GetResponse(&fVals[0],1));
  _n1->setVal(forest->GetResponse(&fVals[0],2));
  _n2->setVal(forest->GetResponse(&fVals[0],3));
  
  //retrieve final pdf parameter values from transformed forest outputs
  ecor = den/_meanlim->getVal();
  mean = _meanlim->getVal();
  sigma = _sigmalim->getVal();
  alpha1 = 2.0; //alpha hardcoded in this version of the regression
  n1 = _n1lim->getVal();
  alpha2 = 1.0; //alpha hardcoded in this version of the regression
  n2 = _n2lim->getVal();
  
  //possible memory leak
  //pdfpeakval = _pdf->getVal(&_normset);
  pdfpeakval = 0.;
  
  return;
  
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrector::CorrectedEnergyWithErrorV6(const Photon *p, const VertexCol *vtxs, Double_t rho, double &ecor, double &sigEoverE, double &cbmean, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval) {
 
  const SuperCluster *s = p->SCluster();
  const BasicCluster *b = s->Seed();
  
  
  Bool_t isbarrel = (s->AbsEta()<1.5);
  
// //basic supercluster variables
  fVals[0] = s->RawEnergy();
  fVals[1] = s->Eta();
  fVals[2] = s->Phi();
  fVals[3] = p->R9();
  fVals[4] = s->EtaWidth();
  fVals[5] = s->PhiWidth();
  fVals[6] = double(s->NClusters());
  fVals[7] = p->HadOverEmTow();
  fVals[8] = rho;
  fVals[9] = double(vtxs->GetEntries());

  fVals[10] = b->Eta()-s->Eta();
  fVals[11] = atan2(sin(b->Phi()-s->Phi()), cos(b->Phi()-s->Phi()));
  fVals[12] = b->Energy()/s->RawEnergy();
  fVals[13] = b->E3x3()/b->E5x5();
  fVals[14] = sqrt(b->CoviEtaiEta()); //sigietaieta
  fVals[15] = sqrt(b->CoviPhiiPhi()); //sigiphiiphi
  fVals[16] = b->CoviEtaiPhi();   //sigietaiphi
  fVals[17] = b->EMax()/b->E5x5(); //crystal energy ratio gap variables
  fVals[18] = b->E2nd()/b->E5x5();
  fVals[19] = b->ETop()/b->E5x5();
  fVals[20] = b->EBottom()/b->E5x5();
  fVals[21] = b->ELeft()/b->E5x5();
  fVals[22] = b->ERight()/b->E5x5();
  fVals[23] = b->E2x5Max()/b->E5x5(); //crystal energy ratio gap variables
  fVals[24] = b->E2x5Top()/b->E5x5();
  fVals[25] = b->E2x5Bottom()/b->E5x5();
  fVals[26] = b->E2x5Left()/b->E5x5();
  fVals[27] = b->E2x5Right()/b->E5x5();

  if (isbarrel) {
    //additional energy ratio (always ~1 for endcap, therefore only included for barrel)
    fVals[28] = b->E5x5()/b->Energy();
    
    //local coordinates and crystal indices (barrel only)
    
    //seed cluster    
    fVals[29] = b->IEta(); //crystal ieta
    fVals[30] = b->IPhi(); //crystal iphi
    fVals[31] = (b->IEta()-1*TMath::Abs(b->IEta())/b->IEta())%5; //submodule boundary eta symmetry
    fVals[32] = (b->IPhi()-1)%2; //submodule boundary phi symmetry
    fVals[33] = (TMath::Abs(b->IEta())<=25)*((b->IEta()-1*TMath::Abs(b->IEta())/b->IEta())%25) + (TMath::Abs(b->IEta())>25)*((b->IEta()-26*TMath::Abs(b->IEta())/b->IEta())%20); //module boundary eta approximate symmetry
    fVals[34] = (b->IPhi()-1)%20; //module boundary phi symmetry
    fVals[35] = b->EtaCry(); //local coordinates with respect to closest crystal center at nominal shower depth
    fVals[36] = b->PhiCry();
    
  }
  else {
    //preshower energy ratio (endcap only)
    fVals[28] = s->PreshowerEnergy()/s->RawEnergy();
  }
  
  double den;
  HybridGBRForestD *forest;
  if (isbarrel) {
    den = s->RawEnergy();
    forest = fReaderDebsemi;
  }
  else {
    den = s->RawEnergy() + s->PreshowerEnergy();
    forest = fReaderDeesemi;
  }
    
  //set raw response variables from GBRForest
  _sigma->setVal(forest->GetResponse(&fVals[0],0));
  _mean->setVal(forest->GetResponse(&fVals[0],1));
  _n1->setVal(forest->GetResponse(&fVals[0],2));
  _n2->setVal(forest->GetResponse(&fVals[0],3));
  
  //retrieve final pdf parameter values from transformed forest outputs
  cbmean = _meanlim->getVal();
  cbsigma = _sigmalim->getVal();
  cbalpha1 = 2.0; //alpha hardcoded in this version of the regression
  cbn1 = _n1lim->getVal();
  cbalpha2 = 1.0; //alpha hardcoded in this version of the regression
  cbn2 = _n2lim->getVal();
  
  _tgt->setVal(cbmean); //evaluate pdf at peak position
  pdfpeakval = _pdf->getVal(*_tgt);
    
  //set final energy and relative energy resolution
  ecor = den*cbmean;
  sigEoverE = cbsigma/cbmean;
    
  return;
  
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrector::CorrectedEnergyWithErrorV7(const Photon *p, const VertexCol *vtxs, Double_t rho, double &ecor, double &sigEoverE, double &cbmean, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval) {
 
  const SuperCluster *s = p->SCluster();
  const BasicCluster *b = s->Seed();
  
  
  Bool_t isbarrel = (s->AbsEta()<1.5);
  
// //basic supercluster variables
  fVals[0] = s->RawEnergy();
  fVals[1] = s->Eta();
  fVals[2] = p->R9();
  fVals[3] = s->EtaWidth();
  fVals[4] = s->PhiWidth();
  fVals[5] = double(s->NClusters());
  fVals[6] = p->HadOverEmTow();
  fVals[7] = rho;
  fVals[8] = double(vtxs->GetEntries());

  fVals[9]  = b->Eta()-s->Eta();
  fVals[10] = atan2(sin(b->Phi()-s->Phi()), cos(b->Phi()-s->Phi()));
  fVals[11] = b->Energy()/s->RawEnergy();
  fVals[12] = b->E3x3()/b->E5x5();
  fVals[13] = sqrt(b->CoviEtaiEta()); //sigietaieta
  fVals[14] = sqrt(b->CoviPhiiPhi()); //sigiphiiphi
  fVals[15] = b->CoviEtaiPhi();   //sigietaiphi
  fVals[16] = b->EMax()/b->E5x5(); //crystal energy ratio gap variables
  fVals[17] = b->E2nd()/b->E5x5();
  fVals[18] = b->ETop()/b->E5x5();
  fVals[19] = b->EBottom()/b->E5x5();
  fVals[20] = b->ELeft()/b->E5x5();
  fVals[21] = b->ERight()/b->E5x5();
  fVals[22] = b->E2x5Max()/b->E5x5(); //crystal energy ratio gap variables
  fVals[23] = b->E2x5Top()/b->E5x5();
  fVals[24] = b->E2x5Bottom()/b->E5x5();
  fVals[25] = b->E2x5Left()/b->E5x5();
  fVals[26] = b->E2x5Right()/b->E5x5();

  if (isbarrel) {
    //additional energy ratio (always ~1 for endcap, therefore only included for barrel)
    fVals[27] = b->E5x5()/b->Energy();
    
    //local coordinates and crystal indices (barrel only)
    
    //seed cluster    
    fVals[28] = b->IEta(); //crystal ieta
    fVals[29] = (b->IEta()-1*TMath::Abs(b->IEta())/b->IEta())%5; //submodule boundary eta symmetry
    fVals[30] = (b->IPhi()-1)%2; //submodule boundary phi symmetry
    fVals[31] = (TMath::Abs(b->IEta())<=25)*((b->IEta()-1*TMath::Abs(b->IEta())/b->IEta())%25) + (TMath::Abs(b->IEta())>25)*((b->IEta()-26*TMath::Abs(b->IEta())/b->IEta())%20); //module boundary eta approximate symmetry
    fVals[32] = (b->IPhi()-1)%20; //module boundary phi symmetry
    fVals[33] = b->EtaCry(); //local coordinates with respect to closest crystal center at nominal shower depth
    fVals[34] = b->PhiCry();
    
  }
  else {
    //preshower energy ratio (endcap only)
    fVals[27] = s->PreshowerEnergy()/s->RawEnergy();
  }
  
  double den;
  HybridGBRForestD *forest;
  if (isbarrel) {
    den = s->RawEnergy();
    forest = fReaderDebsemi;
  }
  else {
    den = s->RawEnergy() + s->PreshowerEnergy();
    forest = fReaderDeesemi;
  }
    
  //set raw response variables from GBRForest
  _sigma->setVal(forest->GetResponse(&fVals[0],0));
  _mean->setVal(forest->GetResponse(&fVals[0],1));
  _n1->setVal(forest->GetResponse(&fVals[0],2));
  _n2->setVal(forest->GetResponse(&fVals[0],3));
  
  //retrieve final pdf parameter values from transformed forest outputs
  cbmean = _meanlim->getVal();
  cbsigma = _sigmalim->getVal();
  cbalpha1 = 2.0; //alpha hardcoded in this version of the regression
  cbn1 = _n1lim->getVal();
  cbalpha2 = 1.0; //alpha hardcoded in this version of the regression
  cbn2 = _n2lim->getVal();
  
  _tgt->setVal(cbmean); //evaluate pdf at peak position
  pdfpeakval = _pdf->getVal(*_tgt);
    
  //set final energy and relative energy resolution
  ecor = den*cbmean;
  sigEoverE = cbsigma/cbmean;
    
  return;
  
}

//--------------------------------------------------------------------------------------------------
void EGEnergyCorrector::CorrectedEnergyWithErrorV8(const Photon *p, const VertexCol *vtxs, Double_t rho, double &ecor, double &sigEoverE, double &cbmean, double &cbsigma, double &cbalpha1, double &cbn1, double &cbalpha2, double &cbn2, double &pdfpeakval) {
 
  const SuperCluster *s = p->SCluster();
  
  Bool_t isbarrel = (s->AbsEta()<1.5);
  
  if (isbarrel) {
    CorrectedEnergyWithErrorV6(p,vtxs,rho,ecor,sigEoverE,cbmean,cbsigma,cbalpha1,cbn1,cbalpha2,cbn2,pdfpeakval);
  }
  else {
    CorrectedEnergyWithErrorV7(p,vtxs,rho,ecor,sigEoverE,cbmean,cbsigma,cbalpha1,cbn1,cbalpha2,cbn2,pdfpeakval);
  }
  
  return;
  
}
