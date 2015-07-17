#include "MitPhysics/Mods/interface/PFTauIdMod.h"
#include "MitAna/DataTree/interface/Names.h"

ClassImp(mithep::PFTauIdMod)

mithep::PFTauIdMod::PFTauIdMod(char const* name/* = "PFTauIdMod"*/, char const* title/* = "Tau identification module"*/) :
  IdMod<PFTau>(name, title)
{
  fInputName = mithep::Names::gkHPSTauBrn;
}

Bool_t
mithep::PFTauIdMod::IsGood(mithep::PFTau const& tau)
{
  fCutFlow->Fill(cAll);

  if (fIsHPSSel) {
    // basic selections
    if (tau.Pt() < fPtMin)
      return false;

    fCutFlow->Fill(cPt);

    if (tau.AbsEta() > fEtaMax)
      return false;

    fCutFlow->Fill(cEta);

    double cutFlow = cEta + 1.;
    for (auto& disc : fDiscriminators) {
      unsigned iD = disc.first;
      CutConfig& cut(disc.second);
      if (cut.second) {
        if (tau.PFTauDiscriminator(iD) < cut.first)
          return false;
      }
      else {
        if (tau.PFTauDiscriminator(iD) > cut.first)
          return false;
      }
      fCutFlow->Fill(cutFlow);
      cutFlow += 1.;
    }
  }
  else {
    if (TMath::Abs(tau.Charge()) - 1.0 > 0.0001)
      return false;

    fCutFlow->Fill(cCharge);

    if (tau.NSignalPFCands() == 0)
      return false;

    fCutFlow->Fill(cNPF);

    if (!tau.LeadChargedHadronPFCand())
      return false;

    fCutFlow->Fill(cLeadPF);

    if (tau.LeadChargedHadronPFCand()->Pt() <= fPtLeadChargedHadronPFCandMin)
      return false;

    fCutFlow->Fill(cLeadPt);

    if (tau.IsoChargedHadronPtSum() >= fIsoChargedHadronPtSumMax)
      return false;

    fCutFlow->Fill(cChargedPtSum);

    if (tau.IsoGammaEtSum() >= fIsoGammaEtSumMax)
      return false;

    fCutFlow->Fill(cGammaEtSum);

    CompositeParticle tauSystem;
    CompositeParticle tauChargedSystem;
    UInt_t nTrk = 0;
    for (UInt_t iPF = 0; iPF != tau.NSignalPFCands(); ++iPF) {
      tauSystem.AddDaughter(tau.SignalPFCand(iPF));
      if (tau.SignalPFCand(iPF) != 0 &&
          tau.SignalPFCand(iPF)->Charge() != 0){
        ++nTrk;
        tauChargedSystem.AddDaughter(tau.SignalPFCand(iPF));
      }
    }
    if (nTrk != 1 && nTrk != 3)
      return false;

    fCutFlow->Fill(cNTrk);

    if (tauSystem.Pt() <= fPtMin)
      return false;

    fCutFlow->Fill(cSystemPt);

    if (tauChargedSystem.Mass() <= fSignalMassMin || tauChargedSystem.Mass() >= fSignalMassMax)
      return false;

    fCutFlow->Fill(cSystemMass);
  }

  return true;
}

void
mithep::PFTauIdMod::IdBegin()
{
  TAxis* xaxis = fCutFlow->GetXaxis();

  if (fIsHPSSel) {
    fCutFlow->SetBins(nHPSCuts, 0., double(nHPSCuts));
    xaxis->SetBinLabel(cAll + 1, "All");
    xaxis->SetBinLabel(cPt + 1, "Pt");
    xaxis->SetBinLabel(cEta + 1, "Eta");
    int iX = fCutFlow->GetNbinsX() + 1;
    for (auto& disc : fDiscriminators) {
      fCutFlow->SetBins(iX, 0., iX);
      xaxis->SetBinLabel(iX, PFTau::PFTauDiscriminatorName(disc.first));
      ++iX;
    }
  }
  else {
    fCutFlow->SetBins(nNonHPSCuts, 0., double(nNonHPSCuts));
    xaxis->SetBinLabel(cAll + 1, "All");
    xaxis->SetBinLabel(cPt + 1, "Pt");
    xaxis->SetBinLabel(cEta + 1, "Eta");
    xaxis->SetBinLabel(cCharge + 1, "Charge");
    xaxis->SetBinLabel(cNPF + 1, "NPF");
    xaxis->SetBinLabel(cLeadPF + 1, "LeadPF");
    xaxis->SetBinLabel(cLeadPt + 1, "LeadPt");
    xaxis->SetBinLabel(cChargedPtSum + 1, "ChargedPtSum");
    xaxis->SetBinLabel(cGammaEtSum + 1, "GammaEtSum");
    xaxis->SetBinLabel(cNTrk + 1, "NTrk");
    xaxis->SetBinLabel(cSystemPt + 1, "SystemPt");
    xaxis->SetBinLabel(cSystemMass + 1, "SystemMass");
  }
}

Long64_t
mithep::PFTauIdMod::GetDiscriminatorMask() const
{
  Long64_t m = 0;
  for (auto& disc : fDiscriminators)
    m |= (1 << disc.first);

  return m;
}
