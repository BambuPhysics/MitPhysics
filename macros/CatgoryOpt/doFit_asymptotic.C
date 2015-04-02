#if !defined(__CLING__) || defined(__ROOTCLING__)
#include <TROOT.h>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TF1.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TNtuple.h"
#include "TMath.h"
#include <sstream>

#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealIntegral.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooSimultaneous.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooFitResult.h"
#include "RooFFTConvPdf.h"
#include "RooPlot.h"
#include "RooMinimizer.h"
#include "RooNLLVar.h"
#include "RooCachedPdf.h"
#include "TTimeStamp.h"
#include "RooConstVar.h"
#include "RooExponential.h"
#include "RooBernstein.h" 
#endif

using namespace std;

void doFit_asymptotic(int which = 2,
		      int numCat  = 5,
		      /*double splitVal1 = 0.54,
		      double splitVal2 = 0.92,
		      double splitVal3 = 0.865,
		      double splitVal4 = 0,
		      double splitVal5 = 0.71,*/
                      /*double splitVal1 = 0.89,
		      double splitVal2 = 0.72,
		      double splitVal3 = 0.55,
		      double splitVal4 = 0.05,
		      double splitVal5 = -100,*/
                      double splitVal1 = 0.05,
		      double splitVal2 = 0.55,
		      double splitVal3 = 0.74,
		      double splitVal4 = 0.89,
		      double splitVal5 = -100,
		      //TString fileName = "appOutput_SM_Oct14.root",
		      //TString modelDir = "/afs/cern.ch/user/f/fabstoec/scratch0/optimizeMVACats/MVAmodel/"
		      //TString fileName = "Tmva_AppOutput_SMDipho_2012Jan16_wrongBR.root", 
		      //TString fileName = "Tmva_AppOutput_SMDipho_2012Jan16.root",
		      //TString fileName = "Tmva_AppOutput_SMDipho_2012Jan16_NewVBF.root", 
                      //TString fileName = "Tmva_AppOutput_SMDipho_2012Jan16_NewVBF_CorrScaleFactor.root", 
                      TString fileName = "Tmva_AppOutput_SMDipho_2012Jan16_NewVBF_120.root ", 
                      TString modelDir = "/afs/cern.ch/user/m/mingyang/optimizeMVACats/MVAmodel/"
		      ) {

  
  if(numCat > 6) return;//set maximum number of categories

  //===============1.remove the printout from roofit=========================
  RooMsgService::instance().setSilentMode(true);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Caching);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Minimization);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Plotting);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Fitting);
  RooMsgService::instance().getStream(0).removeTopic(RooFit::Eval);
  RooMsgService::instance().getStream(1).removeTopic(RooFit::Eval);
  
  //============2.define nonfit variables==========================
  //--------input file-------------------
  TString modelName = modelDir+fileName;
  //--------cross section-----
  double theXS   = 1.0;
  double maxXS   = 10.0;
  //--------the signal strength modifier range ------
  double muMin  = 0.;
  double muMax  = maxXS/theXS;
  double theMu  = 0.;//set the initial value of Mu to 0 (background only hypothesis)
  //--------parameters--------------------------------
  double nsig;
  double nbgini;
  double nbg;
  double exp;
  double nbg_update;
  double exp_update;
  //--------set up category bourndary variables--------
  int startCat = 0;
  int endCat   = numCat+1;

  double* catVals    = new double[numCat+1];//ming:array of mva category cut values

  double allVals[7];//ming:all the boundary values
  allVals[0] = splitVal1;
  allVals[1] = splitVal2;
  allVals[2] = splitVal3;
  allVals[3] = splitVal4;
  allVals[4] = splitVal5;
  //allVals[5] = 1.;
  //allVals[6] = -0.5;
  allVals[5] = 1.;
  allVals[6] = -1;

  double allValsSorted[7];

  for(int i=0; i<7; ++i) {//make allValsSorted increase with the order
    allValsSorted[i]=allVals[i];
    int tSpot = i;
    
    for(int j = i-1; j>=0; --j) {
      if(allValsSorted[j] > allValsSorted[tSpot]) {
	// switch
	double tempVal = allValsSorted[j];	  
	allValsSorted[j] = allValsSorted[tSpot];
	allValsSorted[tSpot] = tempVal;
	tSpot--;
      } else break;
    }
  }  
  
  for(int i=6-numCat; i<=6; ++i)
    catVals[i-6+numCat] = allValsSorted[i];//catVals from -1 to 1 increase with the order
  
  //==========3.read in the input tree==============  
  TFile* modFile = new TFile(modelName.Data());
  TNtuple* tup = (TNtuple*) modFile->FindObjectAny("MITMVAtuple");

  //============fit=============
  //-----------set the independent variable----------------
  RooRealVar* mass = new RooRealVar("mass","",100.,180.);
  mass->setRange(100.,180.);
  mass->setBins(10000,"cache");//bins for numrical integration
  mass->setBins(320);//bins for fitting and plotting
  //-----------set the variable to make dataset-------------
  RooRealVar* weight = new RooRealVar("weight","",0.,100.);
  RooRealVar* proc   = new RooRealVar("proc"  ,"",0.,10.);
  RooRealVar* mva    = new RooRealVar("mva"   ,"",-10.,10.);
  RooRealVar* vbftag = new RooRealVar("vbftag"   ,"",0.,1.); 
  //----------prepare the dataset---------------------------
  TString* dataNames      = new TString[numCat+1];
  TString* sigNames       = new TString[numCat+1];
  RooDataHist** sigCat = new RooDataHist*[numCat+1];
  RooDataSet** dataCat = new RooDataSet*[numCat+1];
  RooDataHist** dataHCat = new RooDataHist*[numCat+1];
  TString baseCut = "mass > 100 && mass < 180 ";
  TString baseBG  = baseCut + TString(" && proc > 3 && proc < 11");//MC background
  TString baseSIG = baseCut + TString(" && proc < 4 ");//MC signal
  TString baseObsData = baseCut + TString(" && proc==11");//obs data for count
  int*    numObsData= new int[numCat+1];
  //----------prepare the pdf-------------------------------
  //signal model
  RooAbsPdf**   sigpdfcat_part1 = new RooAbsPdf*[numCat+1];
  RooAbsPdf**   sigpdfcat_part2 = new RooAbsPdf*[numCat+1];
  RooAbsPdf**   sigpdfcat       = new RooAbsPdf*[numCat+1];
  RooRealVar**  sigMean_part1   = new RooRealVar*[numCat+1];
  RooRealVar**  sigWidth_part1  = new RooRealVar*[numCat+1];
  RooRealVar**  sigMean_part2   = new RooRealVar*[numCat+1];
  RooRealVar**  sigWidth_part2  = new RooRealVar*[numCat+1];
  RooRealVar**  sigRatio        = new RooRealVar*[numCat+1];
  double* nominalSignal   = new double[numCat+1];
  TH1D** sigH = new TH1D*[numCat+1];//ming:used to get nominalSignal and input to RooDataHist
  TH1D** ObsDataH = new TH1D*[numCat+1];//obs data for count
  //background model
  //exp
  RooExponential**  bgpdfcat_exp = new RooExponential*[numCat+1];
  RooRealVar**      bgpdfval_exp = new RooRealVar*[numCat+1];
  //bs
  RooBernstein**    bgpdfcat_bs =  new RooBernstein*[numCat+1];
  RooConstVar*      bgpdfval0_bs = new RooConstVar; 
  RooRealVar**      bgpdfval1_bs = new RooRealVar*[numCat+1];
  RooRealVar**      bgpdfval2_bs = new RooRealVar*[numCat+1];
  RooRealVar**      bgpdfval3_bs = new RooRealVar*[numCat+1];
  RooRealVar**      bgpdfval4_bs = new RooRealVar*[numCat+1];
  //total model
  RooAddPdf**      totpdfcat = new RooAddPdf*[numCat+1];
  RooRealVar** nsigcat = new RooRealVar*[numCat+1];
  RooRealVar** nbgcat  = new RooRealVar*[numCat+1];
  //----------nll for the mu given the dataset---------------
  RooNLLVar** totNllExt = new RooNLLVar*[numCat+1];//total Nll
  //--------draw option------------------
  /*bool drawSigModel = false;
    bool drawBkgModel = false;
    bool draw95CLlimitModel = false;
    bool drawcan = false;*/
  bool drawSigModel = true;
  bool drawBkgModel = true;
  bool draw95CLlimitModel = true;
  bool drawcan = false;
  RooPlot** frame = new RooPlot*[20];//define an array of pointer;frame[i] is a pointer
  TCanvas* thisCan = new TCanvas("thisCan","MyCanvas",1200,(numCat+1)*300); 
  TCanvas* extraCan = new TCanvas("extraCan","extraCan",800,600); 
  thisCan->Divide(3,numCat+1);
  gStyle->SetMarkerSize(0.5);
  //---------count option------------------
  bool enableCount = true;
  //============category boundaries===============
  std::stringstream mvaCut;
  std::stringstream mvaCutLabel;
  TString SIGstring;
  TString BGstring;
  TString ObsDatastring;
  for(int i=0; i<numCat+1; ++i) {
    if(i!= numCat) {
      mvaCut.str("");
      mvaCut << " && mva > "<<catVals[i]<<" && mva < "<<catVals[i+1];//cut values for cat i
      std::cout<<mvaCut.str().c_str()<<std::endl;  
      mvaCutLabel.str("");
      mvaCutLabel << " mva > "<<catVals[i]<<" && mva < "<<catVals[i+1];//cut values for cat i 
    }
    if(i== numCat) {
     mvaCutLabel.str("");
     mvaCutLabel << " Dijet-TAG ";
    }
    //-----obsData for count-----
    //read in the obsdata mass hist for category i
    std::stringstream pSS; 
    pSS.str("");
    pSS << "ObsDataH" <<i;    
    ObsDataH[i] = new TH1D(pSS.str().c_str(),"",320,100.,180.);//0.25 GeV per bin
    TString theDrawString0 = TString("mass>>")+TString(pSS.str().c_str());
    if(i!=numCat){
      ObsDatastring = baseObsData+TString(mvaCut.str().c_str())+TString(" && vbftag==0");
    }
    if(i==numCat){
      ObsDatastring = baseObsData+TString(" && vbftag==1 && mva>0.05");
    }
    extraCan->cd();
    tup->Draw(theDrawString0.Data(),ObsDatastring.Data());
    numObsData[i] = ObsDataH[i]->GetEntries();//number of events in the signal mass hist
     
    //-----signal model-----
    //read in the signal mass hist for category i
    pSS.str("");
    pSS << "sigH" <<i;    
    sigH[i] = new TH1D(pSS.str().c_str(),"",320,100.,180.);//0.25 GeV per bin
    TString theDrawString = TString("mass>>")+TString(pSS.str().c_str());
    //TString SIGstring = TString("weight*(")+baseSIG+TString(mvaCut.str().c_str())+TString(")");//signal cut along with cat cut
    if(i!=numCat){
      SIGstring = TString("weight*(")+baseSIG+TString(mvaCut.str().c_str())+TString(" && vbftag==0")+TString(")");
    }
    if(i==numCat){
      SIGstring = TString("weight*(")+baseSIG+TString(" && vbftag==1")+TString(" && mva>0.05")+TString(")");
    }
    extraCan->cd();
    tup->Draw(theDrawString.Data(),SIGstring.Data());
    //sigH[i] = (TH1D*)gPad->GetPrimitive(pSS.str().c_str());
    nominalSignal[i] = sigH[i]->GetSumOfWeights();//number of events in the signal mass hist
    //obtain the signal RooDataHist to be fit for category i
    pSS.str("");
    pSS << "sigdata" << i;    
    sigNames[i]    = TString(pSS.str().c_str());
    sigCat[i]      = new RooDataHist( sigNames[i].Data(),"",*mass,sigH[i]);//hist used to extract signal model for category i
    //set fit parameter
    pSS.str("");
    pSS << "nsigcat" << i;
    nsigcat[i] = new RooRealVar(pSS.str().c_str(),"",0.);//initialize number of signal events
    nsigcat[i]->setVal(nominalSignal[i]*theMu);//set nsigcat[i] to constant
    pSS.str("");
    pSS<<"sigMean_part1"<<i;
    sigMean_part1[i]  = new RooRealVar(pSS.str().c_str(),"",125.,100.,180.);
    sigMean_part1[i] ->removeRange();
    pSS.str("");
    pSS<<"sigWidth_part1"<<i;
    sigWidth_part1[i] = new RooRealVar(pSS.str().c_str(),"",1.5,-5.,5.);
    sigWidth_part1[i] ->removeRange();    
    pSS.str("");
    pSS<<"sigMean_part2"<<i;
    sigMean_part2[i]  = new RooRealVar(pSS.str().c_str(),"",125.,100.,180.);
    sigMean_part2[i] ->removeRange();
    pSS.str("");
    pSS<<"sigWidth_part2"<<i;
    sigWidth_part2[i] = new RooRealVar(pSS.str().c_str(),"",2.5,-5.,5.);
    sigWidth_part2[i] ->removeRange();    
    pSS.str("");
    pSS<<"sigRatio"<<i;
    sigRatio[i]  = new RooRealVar(pSS.str().c_str(),"",0.5,0.1,2.);
    sigRatio[i]->removeRange();
    //signal pdf
    pSS.str("");
    pSS  << "sigpdfcat_part1" << i;
    sigpdfcat_part1[i]    = new RooGaussian ( pSS.str().c_str(),"",*mass,*sigMean_part1[i],*sigWidth_part1[i]);
    pSS.str("");
    pSS  << "sigpdfcat_part2" << i;
    sigpdfcat_part2[i]    = new RooGaussian ( pSS.str().c_str(),"",*mass,*sigMean_part2[i],*sigWidth_part2[i]);
    pSS.str("");
    pSS << "sigpdfcat" << i;    
    sigpdfcat[i]   = new RooAddPdf(pSS.str().c_str(),"",*sigpdfcat_part1[i],*sigpdfcat_part2[i],*sigRatio[i]);//ratio*s1+(1-ratio)*s2
    //fit 
    sigpdfcat[i]->fitTo(*sigCat[i]);
    //fix signal model parameters
    sigMean_part1[i]  ->setConstant();
    sigWidth_part1[i] ->setConstant();
    sigMean_part2[i]  ->setConstant();
    sigWidth_part2[i] ->setConstant();
    sigRatio[i] ->setConstant();
    if( drawSigModel ) {
      frame[3*i+1]= mass->frame();
      pSS.str("");
      pSS << "125 GeV Higgs MC and Model (5.089 fb-1) ";
      frame[3*i+1]->SetTitle(pSS.str().c_str());
      sigCat[i]->plotOn(frame[3*i+1]);
      sigpdfcat[i]->plotOn(frame[3*i+1]);
      thisCan->cd(3*i+1);
      frame[3*i+1]->Draw();
      TLatex* text1=new TLatex(3.5,23.5,mvaCutLabel.str().c_str());
      text1->SetNDC();
      text1->SetTextAlign(13);
      text1->SetX(0.4);//(0.940);
      text1->SetY(0.8);
      text1->SetTextFont(42);
      text1->SetTextSize(0.065);// dflt=28
      text1->Draw();
      std::stringstream NSigLabel;
      NSigLabel.str("");
      NSigLabel <<"NSig = "<<nominalSignal[i];
      printf("nsig:%d:\n",nominalSignal[i]);
      TLatex* text2=new TLatex(3.5,23.5,NSigLabel.str().c_str());
      text2->SetNDC();
      text2->SetTextAlign(13);
      text2->SetX(0.4);//(0.940);
      text2->SetY(0.7);
      text2->SetTextFont(42);
      text2->SetTextSize(0.065);// dflt=28
      text2->Draw();
    }
    
    //-----background model-----
    //obtain the backgroud RooDataHist to be fit for category i
    pSS.str("");
    pSS << "data_cat" << i;        
    dataNames[i]   = TString(pSS.str().c_str());
    if(i!=numCat){
      BGstring  = baseBG +TString(mvaCut.str().c_str())+TString(" && vbftag==0");
    }
    if(i==numCat){
      BGstring  = baseBG +TString(" && vbftag==1");
    }
    dataCat[i]     = new RooDataSet (dataNames[i].Data(),"",tup,RooArgSet(*mass,*weight,*proc,*mva,*vbftag),BGstring.Data(), "weight");//MC bkg dataset used to extract background model for cat i
    //set fit parameter and pdf
    pSS.str("");
    pSS << "nbgcat" << i;
    nbgini=dataCat[i]->sumEntries();
    nbgcat[i] = new RooRealVar(pSS.str().c_str(),"",nbgini,0,50000);//number of background events;initial value total number of events in the mc background dataset
    nbgcat[i]->removeRange();
    switch(which){
    case 1:
      {
      pSS.str("");
      pSS << "e" <<i;
      bgpdfval_exp[i] = new RooRealVar(pSS.str().c_str(),"",0.,-100.,100.);//exponet parameter
      bgpdfval_exp[i]->removeRange();
      pSS.str("");
      pSS << "bgpdfcat_exp" << i;
      bgpdfcat_exp[i] = new RooExponential(pSS.str().c_str(),"",*mass,*bgpdfval_exp[i]);
      bgpdfcat_exp[i]->fitTo(*dataCat[i]);
      if( drawBkgModel ) {
	frame[3*i+2]= mass->frame();
	pSS.str("");
	//pSS << "MC Bkg Data And Model 4thBSPol (5.089 fb-1) "<<mvaCut.str().c_str();    
        pSS << "MC Bkg Data And Model exp";    
	frame[3*i+2]->SetTitle(pSS.str().c_str());
	dataCat[i]->plotOn(frame[3*i+2]);
	bgpdfcat_exp[i]->plotOn(frame[3*i+2]);
	thisCan->cd(3*i+2);
	frame[3*i+2]->Draw(); 
        std::stringstream NBkgLabel;
	NBkgLabel.str("");
	NBkgLabel <<"NBkg = "<<nbgini;
	TLatex* text3=new TLatex(3.5,23.5,NBkgLabel.str().c_str());
	text3->SetNDC();
	text3->SetTextAlign(13);
	text3->SetX(0.6);//(0.940);
	text3->SetY(0.7);
	text3->SetTextFont(42);
	text3->SetTextSize(0.065);// dflt=28
	text3->Draw();
      }
      //-----geneate asimov dataset from background model-----
      pSS.str("");
      pSS << "dataHCat" << i;        
      dataHCat[i] = bgpdfcat_exp[i]->generateBinned(*mass,dataCat[i]->sumEntries(),RooFit::Asimov(),RooFit::Name(pSS.str().c_str()));
      //---to get the total pdf---
      pSS.str("");
      pSS << "totpdfcat" <<i;    
      totpdfcat[i] = new RooAddPdf(pSS.str().c_str(),"",RooArgList(*sigpdfcat[i],*bgpdfcat_exp[i]),RooArgList(*nsigcat[i],*nbgcat[i]));//build the sig plus bkg model
      totpdfcat[i]->fitTo(*dataHCat[i],RooFit::Extended(true));
      break;
      }
    case 2:
      {
      pSS.str("");
      pSS << "c0";
      bgpdfval0_bs= new RooConstVar(pSS.str().c_str(),"",1.0);//exponet parameter
      pSS.str("");
      pSS << "c1" <<i;
      bgpdfval1_bs[i] = new RooRealVar(pSS.str().c_str(),"",0.,0,10000.);//exponet parameter
      //bgpdfval1_bs[i]->removeRange(); 
      pSS.str("");
      pSS << "c2" <<i;
      bgpdfval2_bs[i] = new RooRealVar(pSS.str().c_str(),"",0.,0,10000.);//exponet parameter
      //bgpdfval2_bs[i]->removeRange(); 
      pSS.str("");
      pSS << "c3" <<i;
      bgpdfval3_bs[i] = new RooRealVar(pSS.str().c_str(),"",0.,0.,10000.);//exponet parameter
      //bgpdfval3_bs[i]->removeRange(); 
      pSS.str("");
      pSS << "c4" <<i;
      bgpdfval4_bs[i] = new RooRealVar(pSS.str().c_str(),"",0.,0.,10000.);//exponet parameter
      //bgpdfval4_bs[i]->removeRange();  
      pSS.str("");
      pSS << "bgpdfcat_bs" << i;
      bgpdfcat_bs[i] = new RooBernstein(pSS.str().c_str(),"",*mass,RooArgList(*bgpdfval0_bs,*bgpdfval1_bs[i],*bgpdfval2_bs[i],*bgpdfval3_bs[i],*bgpdfval4_bs[i]));
      bgpdfcat_bs[i]->fitTo(*dataCat[i]);
      //draw the sig RooDataHist and sig model
      if( drawBkgModel ) {
	frame[3*i+2]= mass->frame();
	pSS.str("");
	//pSS << "MC Bkg Data And Model 4thBSPol (5.089 fb-1) "<<mvaCut.str().c_str();    
        pSS << "MC Bkg Data And Model 4thBSPol";    
	frame[3*i+2]->SetTitle(pSS.str().c_str());
	dataCat[i]->plotOn(frame[3*i+2]);
	bgpdfcat_bs[i]->plotOn(frame[3*i+2]);
	thisCan->cd(3*i+2);
	frame[3*i+2]->Draw(); 
	std::stringstream NBkgLabel;
	NBkgLabel.str("");
	NBkgLabel <<"NBkg = "<<nbgini;
	TLatex* text3=new TLatex(3.5,23.5,NBkgLabel.str().c_str());
	text3->SetNDC();
	text3->SetTextAlign(13);
	text3->SetX(0.6);//(0.940);
	text3->SetY(0.7);
	text3->SetTextFont(42);
	text3->SetTextSize(0.065);// dflt=28
	text3->Draw();
      }
      //-----geneate asimov dataset from background model-----
      pSS.str("");
      pSS << "dataHCat" << i;        
      dataHCat[i] = bgpdfcat_bs[i]->generateBinned(*mass,dataCat[i]->sumEntries(),RooFit::Asimov(),RooFit::Name(pSS.str().c_str()));
      //---to get the total pdf---
      pSS.str("");
      pSS << "totpdfcat" <<i;    
      totpdfcat[i] = new RooAddPdf(pSS.str().c_str(),"",RooArgList(*sigpdfcat[i],*bgpdfcat_bs[i]),RooArgList(*nsigcat[i],*nbgcat[i]));//build the sig plus bkg model
      totpdfcat[i]->fitTo(*dataHCat[i],RooFit::Extended(true));
      break;
      }
    }
    totNllExt[i] = (RooNLLVar*) totpdfcat[i]->createNLL(*dataHCat[i],RooFit::Extended(true));//calculate nll for each category
  }

  // ===========95% C.L.upper limit on theMu===========
  //-----NLL for background only model theMu=0-----
  double minNll    =  0.;
  for(int i=startCat; i<endCat; ++i) minNll += totNllExt[i]->getVal();//ming:the total Nll for theMu=0 (background only hypothesis)
  std::cout<<"  start   = "<<minNll<<std::endl;

  //-----Target NLL for 95% C.L.upper limit on theMu-----
  //double cl = 0.95;
  double targetNll = minNll + 1.92 ; //0.5*pow(ROOT::Math::normal_quantile(1.-0.5*(1.-cl),1.0), 2.);
  
  //-----set the tolerance for the 95% C.L.upper limit on theMu-----
  double tRelAcc = 0.001;
  double tAbsAcc = 0.0005;

  //-----initial value for the eErr and theMu-----
  double rErr  = 0.5*(muMax - muMin);
  theMu = 0.5*(muMax + muMin);
  
  //-----interate theMu till the rErr is less than either absolute value 0.0005 or the relative error 0.001 
  while ( rErr > std::max(tRelAcc * theMu, tAbsAcc) ) {
    double nextNll = 0.;
    for(int i=startCat; i<endCat; ++i) {
      nsigcat[i]->setVal(nominalSignal[i]*theMu);//update nsigcat with new theMu
      totpdfcat[i]->fitTo(*dataHCat[i],RooFit::Extended(true));//refit the totpdfcat with updated theMu
      nextNll += totNllExt[i]->getVal();//totNllExt update with the totpdfcat
      std::cout<<"  for mu="<<theMu<<"  nll = "<<nextNll<<"    ( "<<targetNll<<" )"<<std::endl;
    }
    if(nextNll < targetNll)
      muMin = theMu;
    else
      muMax = theMu;    
    rErr  = 0.5*(muMax - muMin);
    theMu = 0.5*(muMax + muMin);
  }
  
  double tR_XS = theMu   *theXS;//ming: the 95% C.L.limit on signal strength modifier   
  double tR_Er = rErr  *theXS;//ming: the error on the limit
  
  std::cout<<std::endl<<"95% C.L. Limit Mu="<<setprecision(5)<<"  [ RES ]  ( "<<tR_XS<<"  +-  "<<tR_Er<<" ) x sigma (FP/SM)"<<std::endl;

  std::stringstream label95Mu;
  label95Mu.str("");
  label95Mu << "95% C.L. Limit Mu="<<setprecision(5)<<tR_XS<<"  +-  "<<tR_Er; 

  if(draw95CLlimitModel){
    for(int i=startCat; i<endCat; ++i) {
      frame[3*i+3]= mass->frame();
      std::stringstream pSS;
      pSS.str("");
      pSS << "Asimov Data and 95CL limit Mu Model"; 
      frame[3*i+3]->SetTitle(pSS.str().c_str());
      dataHCat[i]->plotOn(frame[3*i+3]);
      totpdfcat[i]->plotOn(frame[3*i+3]);
      thisCan->cd(3*i+3);
      frame[3*i+3]->Draw();
      TLatex* text4=new TLatex(3.5,23.5,label95Mu.str().c_str());
      text4->SetNDC();
      text4->SetTextAlign(13);
      text4->SetX(0.2);//(0.940);
      text4->SetY(0.8);
      text4->SetTextFont(42);
      text4->SetTextSize(0.065);// dflt=28
      text4->Draw();
    }
  }
  
  if(drawcan){
    std::stringstream CanName;
    CanName.str("");
    CanName << "/afs/cern.ch/user/m/mingyang/optimizeMVACats/plot/"<<numCat<<"cat_95CLlimitMu_"<<setprecision(5)<<tR_XS<<".eps";  
    thisCan->SaveAs(CanName.str().c_str());
  }
  
  for(int i=0; i<numCat+1; ++i){
    if(enableCount){
      printf("cat:%d NumObsData:%d\n",i,numObsData[i]);
    }
  }
  return;
}

