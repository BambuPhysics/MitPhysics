#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include "LHAPDF/LHAPDF.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/PDFProducerMod.h"

using namespace mithep;

ClassImp(mithep::PDFProducerMod)

//--------------------------------------------------------------------------------------------------
PDFProducerMod::PDFProducerMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPrintDebug(kFALSE),
  fMCEvInfoName(Names::gkMCEvtInfoBrn),
  fPDFName("cteq65.LHgrid"),
  fRunPDF(kFALSE),
  fIsData(kFALSE),
  fMCEventInfo(0)
{
  // Constructor
}

//--------------------------------------------------------------------------------------------------
void PDFProducerMod::Process()
{
  // Process entries of the tree.

  // Array to be filled
  UInt_t nmembers = LHAPDF::numberPDF() + 1;
  FArrDouble *PDFArr = new FArrDouble(nmembers);

  if (fIsData == kFALSE) {
    LoadBranch(fMCEvInfoName);

    Double_t Q    = fMCEventInfo->Scale();
    Int_t    id1  = fMCEventInfo->Id1();
    Double_t x1   = fMCEventInfo->X1();
    Double_t pdf1 = fMCEventInfo->Pdf1();
    Int_t    id2  = fMCEventInfo->Id2();
    Double_t x2   = fMCEventInfo->X2();
    Double_t pdf2 = fMCEventInfo->Pdf2();

    if (GetFillHist()) {
      hDPDFHisto[0]->Fill(TMath::Min(Q,999.999));
      hDPDFHisto[1]->Fill(TMath::Min(pdf1,999.999));
      hDPDFHisto[2]->Fill(TMath::Min(pdf2,999.999));
      hDPDFHisto[3]->Fill(TMath::Min(x1,0.999));
      hDPDFHisto[4]->Fill(TMath::Min(x2,0.999));
      if (x1 + x2 > 0) {
	hDPDFHisto[5]->Fill(TMath::Min(x1/(x1+x2),0.999));
      }
    }

    if (fRunPDF == kTRUE) {
      ParticleOArr *leptons = GetObjThisEvt<ParticleOArr>(ModNames::gkMergedLeptonsName);
      if (leptons->GetEntries() >= 2) { // Nlep >= 2 to fill it
	if (fPrintDebug) 
      std::cout << "Start loop over PDF members:" << std::endl;

	for (UInt_t i=0; i<nmembers; ++i) {
	  LHAPDF::usePDFMember(i);
	  Double_t newpdf1 = LHAPDF::xfx(x1, Q, id1)/x1;
	  Double_t newpdf2 = LHAPDF::xfx(x2, Q, id2)/x2;
	  Double_t TheWeight = newpdf1/pdf1*newpdf2/pdf2;

	  if (fPrintDebug) {
        std::cout << i << " --> " << pdf1 << " "   << pdf2 << " | "
                  << x1   << " "   << x2   << " | "
                  << id1  << " "   << id2  << " | "
                  << Q    << " : " <<  TheWeight << std::endl;
	  }
	  if (GetFillHist()) {
            hDPDFHisto[6]->Fill(TMath::Min(TheWeight,99.999));
	    hDPDFHisto[7]->Fill(TheWeight);
	  }
	  PDFArr->Add(TheWeight);
	}
      } // Nlep >= 2 to fill it
      else {
	for (UInt_t i=0; i<nmembers; ++i) PDFArr->Add(1.0);
      }
    }
  }
  AddObjThisEvt(PDFArr, "PDFWeights");
}

//--------------------------------------------------------------------------------------------------
void PDFProducerMod::SlaveBegin()
{
  // Setup LHAPDF and book branch and histograms if wanted.

  gSystem->Setenv("LHAPATH","");
  LHAPDF::setVerbosity(LHAPDF::SILENT);
  LHAPDF::initPDFSet(fPDFName.Data());
  LHAPDF::getDescription();

  if (fIsData == kFALSE) {
    ReqBranch(fMCEvInfoName, fMCEventInfo);
  }

  if (GetFillHist()) {
    char sb[1024];
    sprintf(sb,"hDPDFHisto_%d", 0); hDPDFHisto[0] = new TH1D(sb,sb,500,0,1000); 
    sprintf(sb,"hDPDFHisto_%d", 1); hDPDFHisto[1] = new TH1D(sb,sb,500,0,1000); 
    sprintf(sb,"hDPDFHisto_%d", 2); hDPDFHisto[2] = new TH1D(sb,sb,500,0,1000); 
    sprintf(sb,"hDPDFHisto_%d", 3); hDPDFHisto[3] = new TH1D(sb,sb,100,0,1); 
    sprintf(sb,"hDPDFHisto_%d", 4); hDPDFHisto[4] = new TH1D(sb,sb,100,0,1); 
    sprintf(sb,"hDPDFHisto_%d", 5); hDPDFHisto[5] = new TH1D(sb,sb,100,0,1); 
    sprintf(sb,"hDPDFHisto_%d", 6); hDPDFHisto[6] = new TH1D(sb,sb,500,0,100); 
    sprintf(sb,"hDPDFHisto_%d", 7); hDPDFHisto[7] = new TH1D(sb,sb,1,0,1); 
    for (int i=0; i<8; i++) {
      AddOutput(hDPDFHisto[i]);
    }
  }
}
