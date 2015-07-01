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

    for (unsigned iB = 0; iB != PFTau::nDiscriminators; ++iB) {
      if (((fDiscriminators >> iB) & 1) == 1) {
        if (tau.PFTauDiscriminator(iB) < 0.5)
          return false;

        fCutFlow->Fill(cDisc + iB);
      }
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
    for (unsigned iD = 0; iD != PFTau::nDiscriminators; ++iD)
      xaxis->SetBinLabel(cDisc + iD + 1, PFTau::PFTauDiscriminatorName(iD));
  }
  else {
    fCutFlow->SetBins(nNonHPSCuts, 0., double(nNonHPSCuts));
    xaxis->SetBinLabel(cAll + 1, "All");
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
