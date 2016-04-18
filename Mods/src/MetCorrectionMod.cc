#include "MitPhysics/Mods/interface/MetCorrectionMod.h"

#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include "TVector2.h"
#include <cstring>

using namespace mithep;

ClassImp(mithep::MetCorrectionMod)

//--------------------------------------------------------------------------------------------------
void
MetCorrectionMod::SlaveBegin()
{
  switch (fOutputType) {
  case kMet:
    fOutput = new MetArr(1, fOutputName);
    break;
  case kCaloMet:
    fOutput = new CaloMetArr(1, fOutputName);
    break;
  case kPFMet:
    fOutput = new PFMetArr(1, fOutputName);
    break;
  default:
    SendError(kAbortAnalysis, "SlaveBegin",
              "Unknown output Met type %d", fOutputType);
    return;
  }

  PublishObj(fOutput);

  if (fApplyType0) {
    //type0 formula: function of Pt of vectorial sum of PU particles
    if (!fFormulaType0)
      MakeFormula(0);
  }

  if (fApplyType1 && fJetCorrector)
    fJetCorrector->SetUncertaintySigma(fJESUncertaintySigma);

  if (fApplyShift) {
    //XY shift formula: function of nVtx
    if (!fFormulaShiftPx)
      MakeFormula(1);
    if (!fFormulaShiftPy)
      MakeFormula(2);
  }
}

//--------------------------------------------------------------------------------------------------
void
MetCorrectionMod::SlaveTerminate()
{
  if (fOwnJetCorrector) {
    delete fJetCorrector;
    fJetCorrector = 0;
  }
  RetractObj(fOutput->GetName());
  delete fOutput;
  fOutput = 0;
}

//--------------------------------------------------------------------------------------------------
void
MetCorrectionMod::Process()
{
  // Process entries of the tree.

  std::vector<Muon const*> goodMuons;

  if (fMuonGeometricMatch) {
    auto* muonCol = GetObject<MuonCol>("Muons");
    for (unsigned iM = 0; iM != muonCol->GetEntries(); ++iM) {
      auto& muon(*muonCol->At(iM));
      if (muon.IsGlobalMuon() || muon.IsStandaloneMuon())
        goodMuons.push_back(&muon);
    }
  }

  auto* metCol = GetObject<MetCol>(fMetName);
  if (!metCol) {
    SendError(kAbortModule, "Process",
              "Pointer to input met %s is null.",
              fMetName.Data());
    return;
  }

  if (metCol->GetEntries() == 0) {
    SendError(kAbortModule, "Process",
              "Met collection %s is empty",
              fMetName.Data());
    return;
  }

  auto& inMet(*metCol->At(0));
  if ((fOutputType == kPFMet && inMet.ObjType() != kPFMet) ||
      (fOutputType == kCaloMet && inMet.ObjType() != kCaloMet)) {
    SendError(kAbortModule, "Process",
              "Input Met type does not match the output type.");
    return;
  }

  VertexCol* inVertices = 0;
  if (fApplyType0 || fApplyShift) {
    inVertices = GetObject<VertexCol>(fVertexName);
    if (!inVertices) {
      SendError(kAbortModule, "Process",
                "Pointer to input vertices %s is null.",
                fVertexName.Data());
      return;
    }
  }

  PFCandidateCol const* pfCandidates = 0;
  if (fApplyType0 || fApplyUnclustered) {
    pfCandidates = GetObject<PFCandidateCol>(fPFCandidatesName);
    if (!pfCandidates) {
      SendError(kAbortModule, "Process","Pointer to input PFCandidates %s is null.",
                fPFCandidatesName.Data());
      return;
    }
  }
  
  JetCol const* jets = 0;
  if (fApplyType1 || fApplyUnclustered) {
    jets = GetObject<JetCol>(fJetsName);
    if (!jets) {
      SendError(kAbortModule, "Process",
                "Pointer to input jet collection %s is null.",
                fJetsName.Data());
      return;
    }
  }

  // prepare the storage array for the corrected MET
  fOutput->Reset();

  // type0, type1, shift, unclustered
  enum CorrectionType {
    kType0,
    kType1,
    kShift,
    kUnclustered,
    nCorrectionTypes
  };

  TVector2 metCorrection[nCorrectionTypes] = {};
  double sumEtCorrection[nCorrectionTypes] = {};

  // ===== Type 0 corrections, to mitigate pileup ====
  // https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/1473/1.html
  // !!! Not checked prior to Run 2 - use at your own risk
  if (fApplyType0) {
    // prepare the 4-mom sum of the charged PU candidates
    double puX = 0.;
    double puY = 0.;
    double puSum = 0.;

    // get the Z position of the PV
    Double_t ZofPV = inVertices->At(0)->Z();

    for (UInt_t i=0; i<pfCandidates->GetEntries(); ++i) {

      const PFCandidate *pfcand = pfCandidates->At(i);
      // exclude non PU candidates
      if (fabs(pfcand->SourceVertex().Z() - ZofPV) < fMinDz)
	continue;
      // consider only charged hadrons, electrons and muons
      if (pfcand->PFType() != PFCandidate::eHadron &&
          pfcand->PFType() != PFCandidate::eElectron &&
          pfcand->PFType() != PFCandidate::eMuon)
	continue;

      puX += pfcand->Px();
      puY += pfcand->Py();
      puSum += pfcand->Pt();

      // debug
      if (fPrint) {
        std::cout << "PFCand index " << i << " :: this Vtx dZ: "
                  << fabs(pfcand->SourceVertex().Z() - ZofPV) << std::endl;
        std::cout << "PFCand index " << i << " :: sumPUMom Pt: " << std::sqrt(puX * puX + puY * puY) << std::endl;
      }

    }

    // compute the MET Type 0 correction
    Double_t sumPUPt  = std::sqrt(puX * puX + puY * puY);
    Double_t sumPUPtCorr = fFormulaType0->Eval(sumPUPt);
    Double_t sumPUPxCorr = puX * sumPUPtCorr / sumPUPt;
    Double_t sumPUPyCorr = puY * sumPUPtCorr / sumPUPt;
    Double_t sumPUSumEtCorr = puSum * sumPUPtCorr / sumPUPt;

    // correct the MET
    metCorrection[kType0].Set(sumPUPxCorr, sumPUPyCorr);
    sumEtCorrection[kType0] -= sumPUSumEtCorr;

    // debug
    if (fPrint) {
      std::cout << "\n" << std::endl;
      std::cout << "Final sumPUMom Pt Corr: " << sumPUPtCorr << std::endl;
      std::cout << "Final sumPUMom Px Corr: " << sumPUPxCorr << std::endl;
      std::cout << "Final sumPUMom Py Corr: " << sumPUPyCorr << std::endl;
      std::cout << "+++++++ End of type 0 correction scope +++++++\n\n" << std::endl;
    }

  } // end Type 0 correction scope

  // ===== Type 1 corrections, to propagate JEC to MET =====
  // ===== Also collecting unclustered energy as we loop through jets =====
  // https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMetAnalysis#Type_I_Correction
  if (fApplyType1 || fApplyUnclustered) {
    auto* inRho = GetObject<PileupEnergyDensityCol>(Names::gkPileupEnergyDensityBrn);
    if (!inRho) {
      SendError(kAbortModule, "Process",
                "Pointer to input rho collection is null.");
      return;
    }

    double rho = inRho->At(0)->Rho(fRhoAlgo);

    for (unsigned iJ = 0; iJ != jets->GetEntries(); ++iJ) {
      auto& inJet = *jets->At(iJ);
      auto&& inJetRawMom(inJet.RawMom());

      // filter on |eta|
      double absEta = inJet.AbsEta();

      if (absEta > fMaxJetEta)
        continue;

      // filter on jet properties
      if (fMaxEMFraction > 0.) {
        if (inJet.ObjType() == kPFJet) {
          auto& inPFJet = static_cast<PFJet const&>(inJet);
          if ((inPFJet.ChargedEmEnergy() + inPFJet.NeutralEmEnergy()) / inJetRawMom.E() > fMaxEMFraction)
            continue;
        }
        else
          SendError(kWarning, "Process", "MaxEMFraction is set but the input jet correction is not PF.");
      }

      if (fSkipMuons) {
        if (inJet.ObjType() == kPFJet) {
          auto& inPFJet = static_cast<PFJet const&>(inJet);
          for (unsigned iC = 0; iC != inPFJet.NConstituents(); ++iC) {
            auto* cand = inPFJet.PFCand(iC);
            bool isMuon = false;

            if (fMuonGeometricMatch) {
              for (auto* muon : goodMuons) {
                if (MathUtils::DeltaR(*cand, *muon) < 0.01) {
                  isMuon = true;
                  break;
                }
              }
            }
            else {
              isMuon = (cand->Mu() && (cand->Mu()->IsGlobalMuon() || cand->Mu()->IsStandaloneMuon()));
            }

            if (isMuon)
              inJetRawMom -= cand->Mom();
          }
        }
        else
          SendError(kWarning, "Process", "SkipMuons is set but the input jet correction is not PF.");
      }

      // filter on corrected pt

      double fullCorr = 1.; // L1L2L3
      double offsetCorr = 1.; // L1 only

      if (fJetCorrector) {
        // will correct the jet pt internally
        // save the current correction level first
        auto currentMax = fJetCorrector->GetMaxCorrLevel();
        fJetCorrector->SetMaxCorrLevel(Jet::L3);

        if (absEta < 9.9) {
          std::vector<float>&& corr(fJetCorrector->CorrectionFactors(inJet, rho));
          fullCorr = corr.back() * fJetCorrector->UncertaintyFactor(inJet);
          offsetCorr = corr.front();
        }
        else {
          auto modJet(inJet);
          modJet.SetRawPtEtaPhiM(inJetRawMom.Pt(), inJet.Eta() / absEta * 9.9, inJetRawMom.Phi(), inJetRawMom.M());
          std::vector<float>&& corr(fJetCorrector->CorrectionFactors(modJet, rho));
          fullCorr = corr.back() * fJetCorrector->UncertaintyFactor(modJet);
          offsetCorr = corr.front();
        }

        // recover the correction level
        fJetCorrector->SetMaxCorrLevel(currentMax);
      }
      else {
        // use the correction already stored in the jet
        // will cause potential discrepancy with the internal correction version for jets with |eta| > 9.9
        // also will not use JEC uncertainty
        for (unsigned iC : {Jet::L1, Jet::L2, Jet::L3})
          fullCorr *= inJet.CorrectionScale(iC);

        offsetCorr = inJet.L1OffsetCorrectionScale();
      }

      auto fullCorrMom = inJetRawMom * fullCorr;

      if (fullCorrMom.Pt() < fMinJetPt) continue;

      if (fApplyUnclustered) {
        // this jet is considered "clustered". See below for unclustered energy computation.
        metCorrection[kUnclustered] -= TVector2(inJetRawMom.Px(), inJetRawMom.Py());
        sumEtCorrection[kUnclustered] -= inJetRawMom.Pt();
      }

      if (fApplyType1) {
        auto offsetCorrMom = inJetRawMom * offsetCorr;

        // compute the MET Type 1 correction
        metCorrection[kType1] -= TVector2(fullCorrMom.Px() - offsetCorrMom.Px(), fullCorrMom.Py() - offsetCorrMom.Py());
        sumEtCorrection[kType1] += (fullCorrMom.Et() - offsetCorrMom.Et());

        // debug
        if (fPrint) {
          std::cout << "Jet index " << iJ << " :: raw jet Pt:   " << inJetRawMom.Pt() << std::endl;
          std::cout << "Jet index " << iJ << " :: cor jet Pt:   " << fullCorrMom.Pt() << std::endl;
          std::cout << "Jet index " << iJ << " :: offset cor Pt: " << offsetCorrMom.Pt() << std::endl;
        }
      }
    }

    // debug
    if (fApplyType1 && fPrint) {
      std::cout << "\n" << std::endl;
      std::cout << "Final type1 cor Pt: " << metCorrection[kType1].Mod() << std::endl;
      std::cout << "raw Met Pt        : " << inMet.Pt() << std::endl;
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
      xyShiftCorrX = fFormulaShiftPx->Eval(nVtx) * -1.;
      xyShiftCorrY = fFormulaShiftPy->Eval(nVtx) * -1.;

    metCorrection[kShift].Set(xyShiftCorrX, xyShiftCorrY);

    // debug
    if (fPrint) {
      std::cout << "XY shift cor Pt: " << metCorrection[kShift].Mod() << std::endl;
      std::cout << "raw Met Pt     : " << inMet.Pt() << std::endl;
      std::cout << "+++++++ End of XY shift correction scope +++++++\n\n" << std::endl;
    }
  }

  if (fApplyUnclustered) {
  // computation of unclustered energy:
  // vec{Met} = - vec{unclustered} - vec{jets}
  // SumEt = sumUnclustered + sumJets
  // => vec{unclustered} = - vec{Met} - vec{jets}
  //    SumUnclustered = SumEt - sumJets
    if (fPrint) {
      std::cout << "Unclustered cor for MET: (" << inMet.Mex() << ", " << inMet.Mey() << ")" << std::endl;
      std::cout << " total jet pT:           (" << -metCorrection[kUnclustered].X() << ", " << -metCorrection[kUnclustered].Y() << ")" << std::endl;
    }

    metCorrection[kUnclustered] -= TVector2(inMet.Mex(), inMet.Mey());
    sumEtCorrection[kUnclustered] += inMet.SumEt();

    if (fPrint)
      std::cout << " total unclustered:      (" << metCorrection[kUnclustered].X() << ", " << metCorrection[kUnclustered].Y() << ")" << std::endl;

    // multiply by the variation rate to arrive at the actual correction
    metCorrection[kUnclustered] *= fUnclusteredVariation;
    sumEtCorrection[kUnclustered] *= fUnclusteredVariation;

    if (fPrint)
      std::cout << " correction:             (" << metCorrection[kUnclustered].X() << ", " << metCorrection[kUnclustered].Y() << ")" << std::endl;
  }

  // initialize the corrected met to the uncorrected one
  Met* corrected = 0;
  switch (fOutputType) {
  case kMet:
    corrected = static_cast<MetArr*>(fOutput)->Allocate();
    new (corrected) Met(inMet);
    break;
  case kCaloMet:
    corrected = static_cast<CaloMetArr*>(fOutput)->Allocate();
    new (corrected) CaloMet(static_cast<CaloMet&>(inMet));
    break;
  case kPFMet:
    corrected = static_cast<PFMetArr*>(fOutput)->Allocate();
    new (corrected) PFMet(static_cast<PFMet&>(inMet));
    break;
  }
  double mex = inMet.Mex();
  double mey = inMet.Mey();
  double sumEt = inMet.SumEt();

  if (fPrint)
    std::cout << "Original MET: (" << mex << ", " << mey << ")" << std::endl;

  for (unsigned iC = 0; iC != nCorrectionTypes; ++iC) {
    if (fPrint)
      std::cout << " correction " << iC << ": (" << metCorrection[iC].X() << ", " << metCorrection[iC].Y() << ")" << std::endl;

    mex += metCorrection[iC].X();
    mey += metCorrection[iC].Y();
    sumEt += sumEtCorrection[iC];
  }

  if (fPrint)
    std::cout << "Corrected MET: (" << mex << ", " << mey << ")" << std::endl;

  corrected->SetMex(mex);
  corrected->SetMey(mey);
  corrected->SetSumEt(sumEt);
}

void
MetCorrectionMod::AddJetCorrectionFromFile(char const* file)
{
  MakeJetCorrector();
  fJetCorrector->AddParameterFile(file);
}

void
MetCorrectionMod::SetJetCorrector(JetCorrector* corr)
{
  if (fJetCorrector && fOwnJetCorrector) {
    SendError(kAbortAnalysis, "SetJetCorrector", "Corrector already created");
  }
  else {
    fJetCorrector = corr;
  }
}

void
MetCorrectionMod::MakeJetCorrector()
{
  if (!fJetCorrector) {
    fJetCorrector = new JetCorrector;
    fOwnJetCorrector = true;
  }
}

void
MetCorrectionMod::MakeFormula(UInt_t idx, char const* expr/* = ""*/)
{
  // idx is only used internally, so no enum is defined.
  switch (idx) {
  case 0:
    delete fFormulaType0;
    if (std::strlen(expr) == 0)
      expr = "-(-0.703151*x)*(1.0 + TMath::Erf(-0.0303531*TMath::Power(x, 0.909209)))";
    fFormulaType0 = new TFormula("formulaType0", expr);
    break;

  case 1:
    delete fFormulaShiftPx;
    if (std::strlen(expr) == 0) {
      if (fIsData)
        expr = "+4.83642e-02 + 2.48870e-01*x";
      else
        expr = "+1.62861e-01 - 2.38517e-02*x";
    }
    fFormulaShiftPx = new TFormula("formulaShiftPx", expr);
    break;

  case 2:
    delete fFormulaShiftPy;
    if (std::strlen(expr) == 0) {
      if (fIsData)
        expr = "-1.50135e-01 - 8.27917e-02*x";
      else
        expr = "+3.60860e-01 - 1.30335e-01*x";
    }
    fFormulaShiftPy = new TFormula("formulaShiftPy", expr);
    break;

  default:
    break;
  }
}
