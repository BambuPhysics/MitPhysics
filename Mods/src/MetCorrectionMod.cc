#include "MitAna/DataTree/interface/JetCol.h"
#include "MitPhysics/Mods/interface/MetCorrectionMod.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/DataFormats/interface/Vect3.h"
#include "MitAna/DataTree/interface/VertexCol.h"

#include "TLorentzVector.h"

using namespace mithep;

ClassImp(mithep::MetCorrectionMod)

//--------------------------------------------------------------------------------------------------
  MetCorrectionMod::MetCorrectionMod(const char *name, const char *title) : 
    BaseMod(name,title),
    fMetName("PFMet"),
    fCorrectedMetName("PFMetT0T1Shift"),
    fJetsName(ModNames::gkPubJetsName),
    fCorrectedJetsName(ModNames::gkCorrectedJetsName),  
    fPFCandidatesName(Names::gkPFCandidatesBrn),
    fVertexName(ModNames::gkGoodVertexesName),
    fMinDz(0.2),
    fApplyType0(kTRUE),
    fApplyType1(kTRUE),
    fApplyShift(kTRUE),
    fExprType0("-(-0.703151*x)*(1.0 + TMath::Erf(-0.0303531*TMath::Power(x, 0.909209)))"),
    fExprShiftDataPx("+4.83642e-02 + 2.48870e-01*x"),
    fExprShiftDataPy("-1.50135e-01 - 8.27917e-02*x"),
    fExprShiftMCPx  ("+1.62861e-01 - 2.38517e-02*x"),
    fExprShiftMCPy  ("+3.60860e-01 - 1.30335e-01*x"),
    fIsData(kTRUE),
    fPrint(kFALSE),
    fPFMet(0),      
    fPFCandidates(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
MetCorrectionMod::~MetCorrectionMod() {} 

//--------------------------------------------------------------------------------------------------
void MetCorrectionMod::SlaveBegin()
{
  // ===== initialize the formulae ====  
  //type0 formula: function of Pt of vectorial sum of PU particles
  fFormulaType0       = new TFormula("formulaType0",       fExprType0);
  //XY shift formula: function of nVtx
  fFormulaShiftDataPx = new TFormula("formulaShiftDataPx", fExprShiftDataPx);
  fFormulaShiftDataPy = new TFormula("formulaShiftDataPy", fExprShiftDataPy);
  fFormulaShiftMCPx   = new TFormula("formulaShiftMCPx",   fExprShiftMCPx);
  fFormulaShiftMCPy   = new TFormula("formulaShiftMCPy",   fExprShiftMCPy);

  // ===== load PFMET and PFCandidates branches ====
  ReqBranch(fMetName, fPFMet);
  ReqBranch(fPFCandidatesName, fPFCandidates);

}

//--------------------------------------------------------------------------------------------------
void MetCorrectionMod::SlaveTerminate()
{
  // ===== deallocate memory ====
  delete fFormulaType0;
  delete fFormulaShiftDataPx;
  delete fFormulaShiftDataPy;
  delete fFormulaShiftMCPx;
  delete fFormulaShiftMCPy;

}

//--------------------------------------------------------------------------------------------------
void MetCorrectionMod::Process()
{
  // Process entries of the tree. 

  LoadBranch(fMetName);
  const VertexCol *inVertices          = GetObjThisEvt<VertexCol>(fVertexName);
  if (!fPFMet) {
    SendError(kAbortModule, "Process", 
              "Pointer to input met %s is null.",
              fMetName.Data());
    return;
  }
  if (!inVertices) {
    SendError(kAbortModule, "Process", 
              "Pointer to input vertices %s is null.",
              fVertexName.Data());
    return;
  }    

  // prepare the storage array for the corrected MET
  MetOArr *CorrectedMetCol = new MetOArr;
  CorrectedMetCol->SetOwner(kTRUE);
  CorrectedMetCol->SetName(fCorrectedMetName);

  // initialize the corrected met to the uncorrected one
  Met *CorrectedMet = fPFMet->At(0)->MakeCopy();
  
  // ===== Type 0 corrections, to mitigate pileup ====
  // https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1473/1.html  
  if (fApplyType0) {
    LoadBranch(fPFCandidatesName);
    if (!fPFCandidates) {
      SendError(kAbortModule, "Process","Pointer to input PFCandidates %s is null.",
                fPFCandidatesName.Data());
      return;
    }
    // prepare the 4-mom sum of the charged PU candidates
    TLorentzVector sumPUMom(0,0,0,0);

    // get the Z position of the PV
    Double_t ZofPV = inVertices->At(0)->Z();

    for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {      

      const PFCandidate *pfcand = fPFCandidates->At(i);
      // exclude non PU candidates
      if (fabs(pfcand->SourceVertex().Z() - ZofPV) < fMinDz)
	continue;
      // consider only charged hadrons, electrons and muons
      if (pfcand->PFType() != PFCandidate::eHadron &&
          pfcand->PFType() != PFCandidate::eElectron &&
          pfcand->PFType() != PFCandidate::eMuon)
	continue;
      TLorentzVector thisMom(0,0,0,0);
      thisMom.SetPtEtaPhiE(pfcand->Mom().Pt(),pfcand->Mom().Eta(),pfcand->Mom().Phi(),
			   pfcand->Mom().E());
      sumPUMom += thisMom;

      // debug
      if (fPrint) {
        std::cout << "PFCand index " << i << " :: this Vtx dZ: "
                  << fabs(pfcand->SourceVertex().Z() - ZofPV) << std::endl;
        std::cout << "PFCand index " << i << " :: sumPUMom Pt: " << sumPUMom.Pt() << std::endl;
      }
            
    }
    
    // compute the MET Type 0 correction
    Double_t sumPUPt  = sumPUMom.Pt();
    Double_t sumPUPhi = sumPUMom.Phi();
    Double_t sumPUPtCorr = fFormulaType0->Eval(sumPUPt);
    Double_t sumPUPxCorr = TMath::Cos(sumPUPhi)*sumPUPtCorr;
    Double_t sumPUPyCorr = TMath::Sin(sumPUPhi)*sumPUPtCorr;

    // correct the MET
    CorrectedMet->SetMex(CorrectedMet->Mex() + sumPUPxCorr);
    CorrectedMet->SetMey(CorrectedMet->Mey() + sumPUPyCorr);

    // debug
    if (fPrint) {
      std::cout << "\n" << std::endl;
      std::cout << "Final sumPUMom Pt Corr: " << sumPUPtCorr << std::endl;
      std::cout << "Final sumPUMom Phi    : " << sumPUPhi << std::endl;
      std::cout << "Final sumPUMom Px Corr: " << sumPUPxCorr << std::endl;
      std::cout << "Final sumPUMom Py Corr: " << sumPUPyCorr << std::endl;
      std::cout << "raw Met Pt            : " << fPFMet->At(0)->Pt() << std::endl;
      std::cout << "cor Met Pt            : " << CorrectedMet->Pt() << std::endl;
      std::cout << "+++++++ End of type 0 correction scope +++++++\n\n" << std::endl;
    }
    
  } // end Type 0 correction scope

  // ===== Type 1 corrections, to propagate JEC to MET ====
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis#Type_I_Correction
  if (fApplyType1) {
    
    const JetCol *inJets     = GetObjThisEvt<JetCol>(fJetsName);
    const JetCol *inCorrJets = GetObjThisEvt<JetCol>(fCorrectedJetsName);
    if (!inJets) {
      SendError(kAbortModule, "Process", 
                "Pointer to input jet collection %s is null.",
                fJetsName.Data());
      return;
    }
    if (!inJets) {
      SendError(kAbortModule, "Process", 
                "Pointer to input corrected jet collection %s is null.",
                fCorrectedJetsName.Data());
      return;
    }
    if (inJets->GetEntries() != inCorrJets->GetEntries())  {
      SendError(kAbortModule, "Process", 
                "Input corrected and uncorrected jet collections have different size.");
      return;
    }

    // prepare the 4-mom for the Type1 correction
    TLorentzVector type1Mom(0,0,0,0);

    for (UInt_t i=0; i<inJets->GetEntries(); ++i) {
      const Jet *inJet = inJets->At(i);      
      const Jet *inCorrJet = inCorrJets->At(i);      

      // do not propagate JEC for soft jets
      if (inCorrJet->Mom().Pt() < 10.) continue;      
      TLorentzVector thisJetMom(0,0,0,0);
      thisJetMom.SetPtEtaPhiE(inJet->Mom().Pt(),inJet->Mom().Eta(),inJet->Mom().Phi(),
			      inJet->Mom().E());
      TLorentzVector thisCorrJetMom(0,0,0,0);
      thisCorrJetMom.SetPtEtaPhiE(inCorrJet->Mom().Pt(),inCorrJet->Mom().Eta(),
				  inCorrJet->Mom().Phi(),inCorrJet->Mom().E());

      // compute the MET Type 1 correction
      type1Mom = type1Mom + thisJetMom - thisCorrJetMom;
      
      // debug
      if (fPrint) {
        std::cout << "Jet index " << i << " :: raw jet Pt:   " << inJet->Mom().Pt() << std::endl;
        std::cout << "Jet index " << i << " :: cor jet Pt:   " << inCorrJet->Mom().Pt() << std::endl;
        std::cout << "Jet index " << i << " :: type1 cor Pt: " << type1Mom.Pt() << std::endl;
      }
    }
    
    // correct the MET
    CorrectedMet->SetMex(CorrectedMet->Mex() + type1Mom.Px());
    CorrectedMet->SetMey(CorrectedMet->Mey() + type1Mom.Py());

    // debug
    if (fPrint) {
      std::cout << "\n" << std::endl;
      std::cout << "Final type1 cor Pt: " << type1Mom.Pt() << std::endl;
      std::cout << "raw Met Pt        : " << fPFMet->At(0)->Pt() << std::endl;
      std::cout << "cor Met Pt        : " << CorrectedMet->Pt() << std::endl;
      std::cout << "+++++++ End of type 1 correction scope +++++++\n\n" << std::endl;
    }
  }

  // ===== XY Shift correction, to reduce the MET azimuthal modulation ====
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis#xy_Shift_Correction
  // NB: the correction in CMSSW is applied with a minus sign, not as noted on the twiki
  if (fApplyShift) {
    // prepare the correction containers
    Double_t xyShiftCorrX, xyShiftCorrY;
    // number of vertices in Double format: to be used in the correction formula
    Double_t nVtx = inVertices->GetEntries() * 1.;
    // compute the XY Shift correction
    if (fIsData) {
      xyShiftCorrX = fFormulaShiftDataPx->Eval(nVtx) * -1.;
      xyShiftCorrY = fFormulaShiftDataPy->Eval(nVtx) * -1.;
    }
    else {
      xyShiftCorrX = fFormulaShiftMCPx->Eval(nVtx) * -1.;
      xyShiftCorrY = fFormulaShiftMCPy->Eval(nVtx) * -1.;
    }

    // correct the MET
    CorrectedMet->SetMex(CorrectedMet->Mex() + xyShiftCorrX);
    CorrectedMet->SetMey(CorrectedMet->Mey() + xyShiftCorrY);

    // debug
    if (fPrint) {
      std::cout << "XY shift cor Pt: " << sqrt(xyShiftCorrX*xyShiftCorrX + xyShiftCorrY*xyShiftCorrY) << std::endl;
      std::cout << "raw Met Pt     : " << fPFMet->At(0)->Pt() << std::endl;
      std::cout << "cor Met Pt     : " << CorrectedMet->Pt() << std::endl;
      std::cout << "+++++++ End of XY shift correction scope +++++++\n\n" << std::endl;
    }
  }

  // add corrected met to collection
  CorrectedMetCol->AddOwned(CorrectedMet);
  
  // sort according to ptrootcint forward declaration data members
  CorrectedMetCol->Sort();
  
  // add to event for other modules to use
  AddObjThisEvt(CorrectedMetCol);
}
