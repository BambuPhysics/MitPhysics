#include <TH1D.h>
#include <TH2D.h>
#include <TParameter.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataTree/interface/MetCol.h" 	 
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MCParticleFwd.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"

using namespace mithep;

ClassImp(mithep::GeneratorMod)

//--------------------------------------------------------------------------------------------------
GeneratorMod::GeneratorMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fIsData(kFALSE),
  fPrintDebug(kFALSE),
  fCopyArrays(kFALSE),
  fMCPartName(Names::gkMCPartBrn),
  fMCMETName(ModNames::gkMCMETName),
  fMCLeptonsName(ModNames::gkMCLeptonsName),
  fMCAllLeptonsName(ModNames::gkMCAllLeptonsName),
  fMCTausName(ModNames::gkMCTausName),
  fMCNeutrinosName(ModNames::gkMCNeutrinosName),
  fMCQuarksName(ModNames::gkMCQuarksName),
  fMCqqHsName(ModNames::gkMCqqHsName),
  fMCBosonsName(ModNames::gkMCBosonsName),
  fMCPhotonsName(ModNames::gkMCPhotonsName),
  fMCRadPhotonsName(ModNames::gkMCRadPhotonsName),
  fMCISRPhotonsName(ModNames::gkMCISRPhotonsName),
  fPtLeptonMin(0.0),
  fEtaLeptonMax(5.0),
  fPtPhotonMin(0.0),
  fEtaPhotonMax(5.0),
  fPtRadPhotonMin(0.0),
  fEtaRadPhotonMax(5.0),
  fPdgIdCut(0),
  fMassMinCut(-FLT_MAX),
  fMassMaxCut(FLT_MAX),
  fApplyISRFilter(kFALSE),
  fApplyVVFilter(kFALSE),
  fApplyVGFilter(kFALSE),
  fAllowWWEvents(kTRUE),
  fAllowWZEvents(kFALSE),
  fAllowZZEvents(kFALSE),
  fFilterBTEvents(kFALSE)
{
  // Constructor
  fGenLeptons = new MCParticleArr();
  fGenLeptons->SetName(TString("Pub") + ModNames::gkMCLeptonsName);
  fGenAllLeptons = new MCParticleArr();
  fGenAllLeptons->SetName(TString("Pub") + ModNames::gkMCAllLeptonsName);
  fGenTaus = new MCParticleArr();      
  fGenTaus->SetName(TString("Pub") + ModNames::gkMCTausName);
  fGenNeutrinos = new MCParticleArr(); 
  fGenNeutrinos->SetName(TString("Pub") + ModNames::gkMCNeutrinosName);
  fGenQuarks = new MCParticleArr();    
  fGenQuarks->SetName(TString("Pub") + ModNames::gkMCQuarksName);
  fGenqqHs = new MCParticleArr();      
  fGenqqHs->SetName(TString("Pub") + ModNames::gkMCqqHsName);
  fGenBosons = new MCParticleArr();    
  fGenBosons->SetName(TString("Pub") + ModNames::gkMCBosonsName);
  fGenPhotons = new MCParticleArr();   
  fGenPhotons->SetName(TString("Pub") + ModNames::gkMCPhotonsName);
  fGenRadPhotons = new MCParticleArr();
  fGenRadPhotons->SetName(TString("Pub") + ModNames::gkMCRadPhotonsName);
  fGenISRPhotons = new MCParticleArr();
  fGenISRPhotons->SetName(TString("Pub") + ModNames::gkMCISRPhotonsName);
}


//--------------------------------------------------------------------------------------------------
GeneratorMod::~GeneratorMod()
{
  // Destructor.

  delete fGenLeptons;
  delete fGenAllLeptons;
  delete fGenTaus;      
  delete fGenNeutrinos; 
  delete fGenQuarks;    
  delete fGenqqHs;      
  delete fGenBosons;    
  delete fGenPhotons;   
  delete fGenRadPhotons;
  delete fGenISRPhotons;
}

//--------------------------------------------------------------------------------------------------
void GeneratorMod::Process()
{
  // Process entries of the tree.

  // these arrays will be filled in the loop of particles
  MetOArr *GenMet               = new MetOArr;
  GenMet->SetName(fMCMETName);
  GenMet->SetOwner(kTRUE);
  MCParticleOArr *GenLeptons    = new MCParticleOArr;
  GenLeptons->SetName(fMCLeptonsName);
  MCParticleOArr *GenAllLeptons = new MCParticleOArr;
  GenAllLeptons->SetName(fMCAllLeptonsName);
  MCParticleOArr *GenTaus       = new MCParticleOArr; 
  GenTaus->SetName(fMCTausName);
  GenTaus->SetOwner(kTRUE);
  MCParticleOArr *GenNeutrinos  = new MCParticleOArr;
  GenNeutrinos->SetName(fMCNeutrinosName);
  MCParticleOArr *GenQuarks     = new MCParticleOArr;
  GenQuarks->SetName(fMCQuarksName);
  MCParticleOArr *GenqqHs       = new MCParticleOArr;
  GenqqHs->SetName(fMCqqHsName);
  MCParticleOArr *GenBosons     = new MCParticleOArr;
  GenBosons->SetName(fMCBosonsName);
  MCParticleOArr *GenPhotons    = new MCParticleOArr;
  GenPhotons->SetName(fMCPhotonsName);
  MCParticleOArr *GenRadPhotons = new MCParticleOArr;
  GenRadPhotons->SetName(fMCRadPhotonsName);
  MCParticleOArr *GenISRPhotons = new MCParticleOArr;
  GenISRPhotons->SetName(fMCISRPhotonsName);

  MCParticleOArr *GenTempMG0     = new MCParticleOArr;
  MCParticleOArr *GenTempLeptons = new MCParticleOArr;

  Bool_t isOld = kFALSE;
  Int_t sumV[2] = {0, 0}; // W, Z
  Int_t sumVVFlavor[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  Double_t totalMET[3] = {0.0, 0.0, 0.0};
  Bool_t isqqH = kFALSE;
  Double_t ptMin = 999999.;
  if (fIsData == kFALSE) {
    if (fPrintDebug) 
      printf("\n************ Next Event ************\n\n");

    // load MCParticle branch
    auto* particles = GetObject<MCParticleCol>(fMCPartName);
    
    unsigned int nGenWMuons = 0;
    unsigned int nGenWElectrons = 0; //for MonoJet generator-level information
    
    for (UInt_t i=0; i<particles->GetEntries(); ++i) {
      const MCParticle *p = particles->At(i);

      if (fPrintDebug) 
	p->Print("l");

      // rad photons, includes gamma from WWGamma vertex.
      if (p->Is(MCParticle::kGamma) && p->Status() == 1 && 
          p->Pt() > fPtRadPhotonMin && p->AbsEta() < fEtaRadPhotonMax && 
          p->DistinctMother() && p->DistinctMother()->Status() == 3 &&
	 (p->DistinctMother()->Is(MCParticle::kEl)  || p->DistinctMother()->Is(MCParticle::kMu) ||
          p->DistinctMother()->Is(MCParticle::kTau) || p->DistinctMother()->Is(MCParticle::kW))   ) {
	CompositeParticle *object = new CompositeParticle();
	object->AddDaughter(p);
	object->AddDaughter(p->DistinctMother());
	if (object->Mass() > 1.0 || p->DistinctMother()->Is(MCParticle::kW))
	  GenRadPhotons->Add(p);
	delete object;
      }

      // ISR photons
      if (p->Is(MCParticle::kGamma) && p->Status() == 1 && 
          p->Pt() > fPtRadPhotonMin && p->AbsEta() < fEtaRadPhotonMax &&
          p->DistinctMother() && p->DistinctMother()->IsParton()	) {
	GenISRPhotons->Add(p);
      }

      // MET computation at generation level
      if (p->Status() == 1 && !p->IsNeutrino()) {
	totalMET[0] = totalMET[0] + p->Px();
	totalMET[1] = totalMET[1] + p->Py();
	totalMET[2] = totalMET[2] + p->Pz();
      }

      if (!p->IsGenerated()) continue;

      if ((((p->Is(MCParticle::kEl) || p->Is(MCParticle::kMu)) && 
            !p->HasMother(MCParticle::kTau, kFALSE)) || p->Is(MCParticle::kTau)) && 
           p->Status() == 3) {
        GenTempLeptons->Add(p);
      }

	//For Monojet generator-level information:
	if (p->Is(MCParticle::kMu)  && p->Status() == 3 && !p->HasMother(MCParticle::kTau, kFALSE) && p->HasMother(MCParticle::kW, kFALSE)) {
		nGenWMuons++;		
	}
	if (p->Is(MCParticle::kEl)  && p->Status() == 3 && !p->HasMother(MCParticle::kTau, kFALSE) && p->HasMother(MCParticle::kW, kFALSE)) {
		nGenWElectrons++;	
	}
	
	
	
      if ((((p->Is(MCParticle::kEl) || p->Is(MCParticle::kMu)) && 
            !p->HasMother(MCParticle::kTau, kFALSE)) || p->Is(MCParticle::kTau)) && 
           p->Status() == 3 &&
	  (p->HasMother(MCParticle::kW, kFALSE) || p->HasMother(MCParticle::kZ, kFALSE))) {
        if (p->Pt() < ptMin) ptMin = p->Pt();
      }

      // all muons/electrons
      if ((p->Is(MCParticle::kEl) || p->Is(MCParticle::kMu)) && p->Status() == 1) {
	if (p->Pt() > fPtLeptonMin && p->AbsEta() < fEtaLeptonMax) {
          GenAllLeptons->Add(p);
	}
	Bool_t isGoodLepton = kFALSE;
	const MCParticle *pm = p;
	while (pm->HasMother() && isGoodLepton == kFALSE) {
          if (pm->PdgId() == 92) // string reached, terminate loop
            break;
          if (pm->Mother()->Is(MCParticle::kZ)  || pm->Mother()->Is(MCParticle::kW)  ||
              pm->Mother()->Is(MCParticle::kZp) || pm->Mother()->Is(MCParticle::kWp) ||
              pm->Mother()->Is(MCParticle::kH)  || pm->Mother()->Is(MCParticle::kH0) ||
	      pm->Mother()->Is(MCParticle::kA0) || pm->Mother()->Is(MCParticle::kHp)) {
            GenLeptons->Add(p);
            isGoodLepton = kTRUE;
            break;
          }
	  else if (pm->Mother()->Is(MCParticle::kPi0) || pm->Mother()->Is(MCParticle::kEta)) {
            // this is fake, but it is a trick to get rid of these cases and abort the loop
            break;
          } 
          pm = pm->Mother();
	}
      }

      // hadronic taus
      else if (p->Is(MCParticle::kTau) && p->Status() == 2) {
	if (!p->HasDaughter(MCParticle::kEl) && !p->HasDaughter(MCParticle::kMu)) {
          const MCParticle *tv = p->FindDaughter(MCParticle::kTauNu);
          if (tv) {
            MCParticle *pm_f = new MCParticle(*p);
            pm_f->SetMom(p->Px()-tv->Px(), p->Py()-tv->Py(),
                	 p->Pz()-tv->Pz(), p->E()-tv->E());
            GenTaus->AddOwned(pm_f);
          }
	  else {
            //SendError(kWarning, "Process", "Could not find a tau neutrino!");
          }
	}
      }

      // neutrinos
      else if (p->Status() == 1 && p->IsNeutrino()) {
	GenNeutrinos->Add(p);
      }

      // quarks from W/Z decays or top particles
      else if (p->IsQuark() && p->HasMother()) {
        GenQuarks->Add(p);
      }

      // qqH, information about the forward jets
      else if (isqqH == kFALSE && p->Is(MCParticle::kH)) {
	isqqH = kTRUE;
	const MCParticle *pq1 = particles->At(i-1);
	const MCParticle *pq2 = particles->At(i-2);
	if (!pq1 || !pq2) {
            SendError(kWarning, "Process", "Could not find quark pair!");
	} else if (pq1->IsQuark()   && pq2->IsQuark()   && 
                   pq1->HasMother() && pq2->HasMother() &&
                   pq1->Mother() == pq2->Mother()) {
          GenqqHs->Add(pq1);
          GenqqHs->Add(pq2);
	}

	if (p->Status() == 3)  
          GenBosons->Add(p); // take higgs boson in account here rather in next else if 
      }

      // information about bosons: W, Z, h, Z', W', H0, A0, H+
      else if ((p->Status() == 3 &&
               (p->Is(MCParticle::kZ)    || p->Is(MCParticle::kW)   || p->Is(MCParticle::kH) ||
        	p->Is(MCParticle::kZp)   || p->Is(MCParticle::kZpp) ||
        	p->Is(MCParticle::kH0)   || p->Is(MCParticle::kA0)  || p->Is(MCParticle::kHp))) ||
	       (p->Status() == 2 &&
	       (p->Is(MCParticle::kJPsi) || p->Is(MCParticle::kUpsilon)))) {
	GenBosons->Add(p);
	if     (p->Is(MCParticle::kW)) sumV[0]++;
	else if (p->Is(MCParticle::kZ)) sumV[1]++;
	if     (p->Is(MCParticle::kW) && p->HasDaughter(MCParticle::kMu)  && 
        	p->HasDaughter(MCParticle::kMuNu))
          sumVVFlavor[0]++;
	else if (p->Is(MCParticle::kW) && p->HasDaughter(MCParticle::kEl)  && 
        	p->HasDaughter(MCParticle::kElNu))
          sumVVFlavor[1]++;
	else if (p->Is(MCParticle::kW) && p->HasDaughter(MCParticle::kTau) && 
        	p->HasDaughter(MCParticle::kTauNu))
          sumVVFlavor[2]++;
	else if (p->Is(MCParticle::kZ) && p->HasDaughter(MCParticle::kMu,kTRUE) && 
        	p->HasDaughter(-1*MCParticle::kMu,kTRUE))
          sumVVFlavor[3]++;
	else if (p->Is(MCParticle::kZ) && p->HasDaughter(MCParticle::kEl,kTRUE) && 
        	p->HasDaughter(-1*MCParticle::kEl,kTRUE))
          sumVVFlavor[4]++;
	else if (p->Is(MCParticle::kZ) && p->HasDaughter(MCParticle::kTau,kTRUE) && 
        	p->HasDaughter(-1*MCParticle::kTau,kTRUE))
          sumVVFlavor[5]++;
	else if (p->Is(MCParticle::kZ) && p->HasDaughter(MCParticle::kMuNu,kTRUE) && 
        	p->HasDaughter(-1*MCParticle::kMuNu,kTRUE))
          sumVVFlavor[6]++;
	else if (p->Is(MCParticle::kZ) && p->HasDaughter(MCParticle::kElNu,kTRUE) && 
        	p->HasDaughter(-1*MCParticle::kElNu,kTRUE))
          sumVVFlavor[7]++;
	else if (p->Is(MCParticle::kZ) && p->HasDaughter(MCParticle::kTauNu,kTRUE) && 
        	p->HasDaughter(-1*MCParticle::kTauNu,kTRUE))
          sumVVFlavor[8]++;
      }

      // photons
      else if (p->Status() == 1 && p->Is(MCParticle::kGamma) &&
               p->Pt() > fPtPhotonMin && p->AbsEta() < fEtaPhotonMax) {
	GenPhotons->Add(p);
      }

      // W/Z -> lnu for Madgraph
      if (p->IsParton() && p->NDaughters() >= 2) {
	CompositeParticle diBoson;
	if (p->HasDaughter(MCParticle::kMu) && p->HasDaughter(MCParticle::kMuNu)) {
          isOld = kFALSE;
	  for (UInt_t nl = 0; nl < GenTempMG0->GetEntries(); nl++) {
	    if (p->FindDaughter(MCParticle::kMu) == GenTempMG0->At(nl)) {
	      isOld = kTRUE;
	      break;
	    }
	  }
	  if (isOld == kFALSE) {
	    GenTempMG0->Add(p->FindDaughter(MCParticle::kMu));
	    diBoson.AddDaughter(p->FindDaughter(MCParticle::kMu));
            diBoson.AddDaughter(p->FindDaughter(MCParticle::kMuNu));
	    sumV[0]++;
            sumVVFlavor[0]++;
            if (GetFillHist()) 
              hDVMass[0]->Fill(TMath::Min(diBoson.Mass(),199.999));
            const MCParticle *tmp_mu = p->FindDaughter(MCParticle::kMu);
            while (tmp_mu->HasDaughter(MCParticle::kMu) && 
          	   tmp_mu->FindDaughter(MCParticle::kMu)->IsGenerated())
              tmp_mu = tmp_mu->FindDaughter(MCParticle::kMu);	  

            GenLeptons->Add(tmp_mu);
	  }
	}
	if (p->HasDaughter(MCParticle::kEl) && p->HasDaughter(MCParticle::kElNu)) {
          isOld = kFALSE;
	  for (UInt_t nl = 0; nl < GenTempMG0->GetEntries(); nl++) {
	    if (p->FindDaughter(MCParticle::kEl) == GenTempMG0->At(nl)) {
	      isOld = kTRUE;
	      break;
	    }
	  }
	  if (isOld == kFALSE) {
	    GenTempMG0->Add(p->FindDaughter(MCParticle::kEl));
            diBoson.AddDaughter(p->FindDaughter(MCParticle::kEl));
            diBoson.AddDaughter(p->FindDaughter(MCParticle::kElNu));
	    sumV[0]++;
            sumVVFlavor[1]++;
            if (GetFillHist()) 
              hDVMass[1]->Fill(TMath::Min(diBoson.Mass(),199.999));
            const MCParticle *tmp_e = p->FindDaughter(MCParticle::kEl);
            while (tmp_e->HasDaughter(MCParticle::kEl) && 
        	   tmp_e->FindDaughter(MCParticle::kEl)->IsGenerated())
              tmp_e = tmp_e->FindDaughter(MCParticle::kEl);       
            GenLeptons->Add(tmp_e);
          }
	}
	if (p->HasDaughter(MCParticle::kTau) && p->HasDaughter(MCParticle::kTauNu)) {
          isOld = kFALSE;
	  for (UInt_t nl = 0; nl < GenTempMG0->GetEntries(); nl++) {
	    if (p->FindDaughter(MCParticle::kTau) == GenTempMG0->At(nl)) {
	      isOld = kTRUE;
	      break;
	    }
	  }
	  if (isOld == kFALSE) {
	    GenTempMG0->Add(p->FindDaughter(MCParticle::kTau));
            diBoson.AddDaughter(p->FindDaughter(MCParticle::kTau));
            diBoson.AddDaughter(p->FindDaughter(MCParticle::kTauNu));
	    sumV[0]++;
            sumVVFlavor[2]++;
            if (GetFillHist()) 
              hDVMass[2]->Fill(TMath::Min(diBoson.Mass(),199.999));
            const MCParticle *tau = p->FindDaughter(MCParticle::kTau);
            if (tau->HasDaughter(MCParticle::kMu)) 
              GenLeptons->Add(tau->FindDaughter(MCParticle::kMu));
            if (tau->HasDaughter(MCParticle::kEl)) 
              GenLeptons->Add(tau->FindDaughter(MCParticle::kEl));
            if (tau->HasDaughter(MCParticle::kTau)) {
              const MCParticle *tau_second = tau->FindDaughter(MCParticle::kTau);
              if (tau_second->HasDaughter(MCParticle::kMu)) 
        	GenLeptons->Add(tau_second->FindDaughter(MCParticle::kMu));
              if (tau_second->HasDaughter(MCParticle::kEl)) 
        	GenLeptons->Add(tau_second->FindDaughter(MCParticle::kEl));
            }
	  }
	}
	if (p->HasDaughter(MCParticle::kMu,kTRUE) && p->HasDaughter(-1*MCParticle::kMu,kTRUE)) {
          isOld = kFALSE;
	  for (UInt_t nl = 0; nl < GenTempMG0->GetEntries(); nl++) {
	    if (p->FindDaughter(MCParticle::kMu,kTRUE) == GenTempMG0->At(nl)) {
	      isOld = kTRUE;
	      break;
	    }
	  }
	  if (isOld == kFALSE) {
	    GenTempMG0->Add(p->FindDaughter(MCParticle::kMu,kTRUE));
            diBoson.AddDaughter(p->FindDaughter(MCParticle::kMu,kTRUE));
            diBoson.AddDaughter(p->FindDaughter(-1*MCParticle::kMu,kTRUE));
	    sumV[1]++;
            sumVVFlavor[3]++;
            if (GetFillHist()) 
              hDVMass[3]->Fill(TMath::Min(diBoson.Mass(),199.999));
            const MCParticle *tmp_mu0 = p->FindDaughter(MCParticle::kMu,kTRUE);
            while (tmp_mu0->HasDaughter(MCParticle::kMu) && 
          	   tmp_mu0->FindDaughter(MCParticle::kMu)->IsGenerated())
              tmp_mu0 = tmp_mu0->FindDaughter(MCParticle::kMu);	    
            const MCParticle *tmp_mu1 = p->FindDaughter(-1*MCParticle::kMu,kTRUE);
            while (tmp_mu1->HasDaughter(MCParticle::kMu) && 
          	   tmp_mu1->FindDaughter(MCParticle::kMu)->IsGenerated())
              tmp_mu1 = tmp_mu1->FindDaughter(MCParticle::kMu);	    
            GenLeptons->Add(tmp_mu0);
            GenLeptons->Add(tmp_mu1);
	  }
	}
	if (p->HasDaughter(MCParticle::kEl,kTRUE) && p->HasDaughter(-1*MCParticle::kEl,kTRUE)) {
          isOld = kFALSE;
	  for (UInt_t nl = 0; nl < GenTempMG0->GetEntries(); nl++) {
	    if (p->FindDaughter(MCParticle::kEl,kTRUE) == GenTempMG0->At(nl)) {
	      isOld = kTRUE;
	      break;
	    }
	  }
	  if (isOld == kFALSE) {
	    GenTempMG0->Add(p->FindDaughter(MCParticle::kEl,kTRUE));
            diBoson.AddDaughter(p->FindDaughter(MCParticle::kEl,kTRUE));
            diBoson.AddDaughter(p->FindDaughter(-1*MCParticle::kEl,kTRUE));
	    sumV[1]++;
            sumVVFlavor[4]++;
            if (GetFillHist()) 
              hDVMass[4]->Fill(TMath::Min(diBoson.Mass(),199.999));
            const MCParticle *tmp_e0 = p->Daughter(0);
            while (tmp_e0->HasDaughter(MCParticle::kEl) && 
          	   tmp_e0->FindDaughter(MCParticle::kEl)->IsGenerated())
              tmp_e0 = tmp_e0->FindDaughter(MCParticle::kEl);	  
            const MCParticle *tmp_e1 = p->Daughter(1);
            while (tmp_e1->HasDaughter(MCParticle::kEl) && 
          	   tmp_e1->FindDaughter(MCParticle::kEl)->IsGenerated())
              tmp_e1 = tmp_e1->FindDaughter(MCParticle::kEl);	  
            GenLeptons->Add(tmp_e0);
            GenLeptons->Add(tmp_e1);
	  }
	}
	if (p->HasDaughter(MCParticle::kTau,kTRUE) && p->HasDaughter(-1*MCParticle::kTau,kTRUE)) {
          isOld = kFALSE;
	  for (UInt_t nl = 0; nl < GenTempMG0->GetEntries(); nl++) {
	    if (p->FindDaughter(MCParticle::kTau,kTRUE) == GenTempMG0->At(nl)) {
	      isOld = kTRUE;
	      break;
	    }
	  }
	  if (isOld == kFALSE) {
	    GenTempMG0->Add(p->FindDaughter(MCParticle::kTau,kTRUE));
            diBoson.AddDaughter(p->FindDaughter(MCParticle::kTau,kTRUE));
            diBoson.AddDaughter(p->FindDaughter(-1*MCParticle::kTau,kTRUE));
	    sumV[1]++;
            sumVVFlavor[5]++;
            if (GetFillHist()) 
              hDVMass[5]->Fill(TMath::Min(diBoson.Mass(),199.999));
            const MCParticle *tau0 = p->Daughter(0);
            if (tau0->HasDaughter(MCParticle::kMu)) 
              GenLeptons->Add(tau0->FindDaughter(MCParticle::kMu));
            if (tau0->HasDaughter(MCParticle::kEl)) 
              GenLeptons->Add(tau0->FindDaughter(MCParticle::kEl));
            const MCParticle *tau1 = p->Daughter(1);
            if (tau1->HasDaughter(MCParticle::kMu)) 
              GenLeptons->Add(tau1->FindDaughter(MCParticle::kMu));
            if (tau1->HasDaughter(MCParticle::kEl)) 
              GenLeptons->Add(tau1->FindDaughter(MCParticle::kEl));
            if (tau0->HasDaughter(MCParticle::kTau)) {
              const MCParticle *tau0_second = tau0->FindDaughter(MCParticle::kTau);
              if (tau0_second->HasDaughter(MCParticle::kMu)) 
        	GenLeptons->Add(tau0_second->FindDaughter(MCParticle::kMu));
              if (tau0_second->HasDaughter(MCParticle::kEl)) 
        	GenLeptons->Add(tau0_second->FindDaughter(MCParticle::kEl));
            }
            if (tau1->HasDaughter(MCParticle::kTau)) {
              const MCParticle *tau1_second = tau1->FindDaughter(MCParticle::kTau);
              if (tau1_second->HasDaughter(MCParticle::kMu)) 
        	GenLeptons->Add(tau1_second->FindDaughter(MCParticle::kMu));
              if (tau1_second->HasDaughter(MCParticle::kEl)) 
        	GenLeptons->Add(tau1_second->FindDaughter(MCParticle::kEl));
            }
	  }
	}
	if (p->HasDaughter(MCParticle::kMuNu,kTRUE) && p->HasDaughter(-1*MCParticle::kMuNu,kTRUE)) {
          isOld = kFALSE;
	  for (UInt_t nl = 0; nl < GenTempMG0->GetEntries(); nl++) {
	    if (p->FindDaughter(MCParticle::kMuNu,kTRUE) == GenTempMG0->At(nl)) {
	      isOld = kTRUE;
	      break;
	    }
	  }
	  if (isOld == kFALSE) {
	    GenTempMG0->Add(p->FindDaughter(MCParticle::kMuNu,kTRUE));
            diBoson.AddDaughter(p->FindDaughter(MCParticle::kMuNu,kTRUE));
            diBoson.AddDaughter(p->FindDaughter(-1*MCParticle::kMuNu,kTRUE));
	    sumV[1]++;
            sumVVFlavor[6]++;
            if (GetFillHist()) 
              hDVMass[6]->Fill(TMath::Min(diBoson.Mass(),199.999));
	  }
	}
	if (p->HasDaughter(MCParticle::kElNu,kTRUE) && p->HasDaughter(-1*MCParticle::kElNu,kTRUE)) {
          isOld = kFALSE;
	  for (UInt_t nl = 0; nl < GenTempMG0->GetEntries(); nl++) {
	    if (p->FindDaughter(MCParticle::kElNu,kTRUE) == GenTempMG0->At(nl)) {
	      isOld = kTRUE;
	      break;
	    }
	  }
	  if (isOld == kFALSE) {
	    GenTempMG0->Add(p->FindDaughter(MCParticle::kElNu,kTRUE));
            diBoson.AddDaughter(p->FindDaughter(MCParticle::kElNu,kTRUE));
            diBoson.AddDaughter(p->FindDaughter(-1*MCParticle::kElNu,kTRUE));
	    sumV[1]++;
            sumVVFlavor[7]++;
            if (GetFillHist()) 
              hDVMass[7]->Fill(TMath::Min(diBoson.Mass(),199.999));
	  }
	}
	if (p->HasDaughter(MCParticle::kTauNu,kTRUE) && p->HasDaughter(-1*MCParticle::kTauNu,kTRUE)) {
          isOld = kFALSE;
	  for (UInt_t nl = 0; nl < GenTempMG0->GetEntries(); nl++) {
	    if (p->FindDaughter(MCParticle::kTauNu,kTRUE) == GenTempMG0->At(nl)) {
	      isOld = kTRUE;
	      break;
	    }
	  }
	  if (isOld == kFALSE) {
	    GenTempMG0->Add(p->FindDaughter(MCParticle::kTauNu,kTRUE));
            diBoson.AddDaughter(p->FindDaughter(MCParticle::kTauNu,kTRUE));
            diBoson.AddDaughter(p->FindDaughter(-1*MCParticle::kTauNu,kTRUE));
	    sumV[1]++;
            sumVVFlavor[8]++;
            if (GetFillHist()) 
              hDVMass[8]->Fill(TMath::Min(diBoson.Mass(),199.999));
	  }
	}
      }

      // t -> lnu for Madgraph
      if (p->Is(MCParticle::kTop)) {
	CompositeParticle diBoson;
	if (p->HasDaughter(MCParticle::kMu) && p->HasDaughter(MCParticle::kMuNu)) {
          diBoson.AddDaughter(p->FindDaughter(MCParticle::kMu));
          diBoson.AddDaughter(p->FindDaughter(MCParticle::kMuNu));
          if (GetFillHist()) 
            hDVMass[9]->Fill(TMath::Min(diBoson.Mass(),199.999));
          GenLeptons->Add(p->FindDaughter(MCParticle::kMu));
	}    
	else if (p->HasDaughter(MCParticle::kEl) && p->HasDaughter(MCParticle::kElNu)) {
          diBoson.AddDaughter(p->FindDaughter(MCParticle::kEl));
          diBoson.AddDaughter(p->FindDaughter(MCParticle::kElNu));
          if (GetFillHist()) 
            hDVMass[10]->Fill(TMath::Min(diBoson.Mass(),199.999));
          GenLeptons->Add(p->FindDaughter(MCParticle::kEl));
	}    
	else if (p->HasDaughter(MCParticle::kTau) && p->HasDaughter(MCParticle::kTauNu)) {
          diBoson.AddDaughter(p->FindDaughter(MCParticle::kTau));
          diBoson.AddDaughter(p->FindDaughter(MCParticle::kTauNu));
          if (GetFillHist()) 
            hDVMass[11]->Fill(TMath::Min(diBoson.Mass(),199.999));
          const MCParticle *tau = p->FindDaughter(MCParticle::kTau);
          if (tau->HasDaughter(MCParticle::kMu)) 
            GenLeptons->Add(tau->FindDaughter(MCParticle::kMu));
          if (tau->HasDaughter(MCParticle::kEl)) 
            GenLeptons->Add(tau->FindDaughter(MCParticle::kEl));
	}
	else if (!p->HasDaughter(MCParticle::kW)) {
          for (UInt_t nd=0; nd<p->NDaughters(); ++nd)
            if (p->Daughter(nd)->IsNot(MCParticle::kBottom) &&
        	p->Daughter(nd)->IsNot(MCParticle::kGamma)) 
              diBoson.AddDaughter(p->Daughter(nd));
          if (GetFillHist()) 
            hDVMass[12]->Fill(TMath::Min(diBoson.Mass(),199.999));
	}
      }

      // mass cut for given pid
      if (fPdgIdCut && p->Is(fPdgIdCut) && 
	 (p->Mass() < fMassMinCut || p->Mass() > fMassMaxCut)) {
	delete GenTempMG0;
	delete GenTempLeptons;
	SkipEvent();
	return;
      }   
    } // end loop of particles
    
    //Filling generated muon/electron histograms for MonoJet selection
    if (GetFillHist()) {
    	if (nGenWMuons > 0 && nGenWElectrons == 0) {
    	hDGenWMuons->Fill(nGenWMuons);
    	hDGenWElectrons->Fill(0);
    	}
    	else if (nGenWElectrons > 0 && nGenWMuons == 0) {
    	hDGenWElectrons->Fill(nGenWElectrons);
    	hDGenWMuons->Fill(nGenWMuons);
    	}
    	else{
    	hDGenWElectrons->Fill(0);
    	hDGenWMuons->Fill(0);
    	}
    }

    delete GenTempMG0;
    if (GetFillHist()) {
      hDGenAllLeptons[5]->Fill(GenTempLeptons->GetEntries());
      Double_t theHighMass = 0.0;
      Double_t theLowMass  = 1000.;
      if (GenTempLeptons->GetEntries() >= 3) {
        for (unsigned int i=0; i<GenTempLeptons->GetEntries(); i++) {
          const MCParticle *geni = GenTempLeptons->At(i);
          for (unsigned int j=i+1; j<GenTempLeptons->GetEntries(); j++) {
  	    const MCParticle *genj = GenTempLeptons->At(j);
  	    CompositeParticle dilepton;
  	    dilepton.AddDaughter(geni);
  	    dilepton.AddDaughter(genj);
            if (geni->AbsPdgId() == genj->AbsPdgId() && geni->PdgId() != genj->PdgId()) {
  	      if (dilepton.Mass() < theLowMass) {theLowMass = dilepton.Mass();}         
            } // same-flavor, opposite-charge
  	    if (dilepton.Mass() > theHighMass) {theHighMass = dilepton.Mass();}         
          } // loop j
        } // loop i
      } // at least three leptons
      hDGenAllLeptons[6]->Fill(TMath::Min(theLowMass,99.999));
      hDGenAllLeptons[7]->Fill(TMath::Min(theHighMass,99.999));
    }
    if (fApplyVGFilter == kTRUE) {
      Double_t theLowMass  = 1000.;
      for (unsigned int i=0; i<GenTempLeptons->GetEntries(); i++) {
        const MCParticle *geni = GenTempLeptons->At(i);
        for (unsigned int j=i+1; j<GenTempLeptons->GetEntries(); j++) {
          const MCParticle *genj = GenTempLeptons->At(j);
          CompositeParticle dilepton;
          dilepton.AddDaughter(geni);
          dilepton.AddDaughter(genj);
          if (geni->AbsPdgId() == genj->AbsPdgId() && geni->PdgId() != genj->PdgId()) {
            if (dilepton.Mass() < theLowMass) {theLowMass = dilepton.Mass();}	     
          } // same-flavor, opposite-charge
        } // loop j
      } // loop i
      if (theLowMass > 12) {delete GenTempLeptons; SkipEvent(); return;}
    }
    delete GenTempLeptons;

    // --------------------------------
    // Begin special study about VVjets
    // --------------------------------
    if (sumV[0] + 4*sumV[1] == 2 || sumV[0] + 4*sumV[1] == 5 || sumV[0] + 4*sumV[1] == 8) {
      MCParticleOArr *GenTempMG1    = new MCParticleOArr;
      Double_t diBosonMass[2] = {0., 0.};
      for (UInt_t i=0; i<particles->GetEntries(); ++i) {
	const MCParticle *p = particles->At(i);

	if (p->IsParton() && p->NDaughters() >= 2) {
	  CompositeParticle diBoson;
	  if (p->HasDaughter(MCParticle::kMu) && p->HasDaughter(MCParticle::kMuNu)) {
            isOld = kFALSE;
	    for (UInt_t nl = 0; nl < GenTempMG1->GetEntries(); nl++) {
	      if (p->FindDaughter(MCParticle::kMu) == GenTempMG1->At(nl)) {
		isOld = kTRUE;
		break;
	      }
	    }
	    if (isOld == kFALSE) {
	      GenTempMG1->Add(p->FindDaughter(MCParticle::kMu));
	      diBoson.AddDaughter(p->FindDaughter(MCParticle::kMu));
              diBoson.AddDaughter(p->FindDaughter(MCParticle::kMuNu));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 2)
        	hDVVMass[0]->Fill(TMath::Min(diBoson.Mass(),199.999));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 5)
        	hDVVMass[1]->Fill(TMath::Min(diBoson.Mass(),199.999));
	    }
	  }
	  if (p->HasDaughter(MCParticle::kEl) && p->HasDaughter(MCParticle::kElNu)) {
            isOld = kFALSE;
	    for (UInt_t nl = 0; nl < GenTempMG1->GetEntries(); nl++) {
	      if (p->FindDaughter(MCParticle::kEl) == GenTempMG1->At(nl)) {
		isOld = kTRUE;
		break;
	      }
	    }
	    if (isOld == kFALSE) {
	      GenTempMG1->Add(p->FindDaughter(MCParticle::kEl));
              diBoson.AddDaughter(p->FindDaughter(MCParticle::kEl));
              diBoson.AddDaughter(p->FindDaughter(MCParticle::kElNu));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 2)
        	hDVVMass[2]->Fill(TMath::Min(diBoson.Mass(),199.999));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 5)
        	hDVVMass[3]->Fill(TMath::Min(diBoson.Mass(),199.999));
            }
	  }
	  if (p->HasDaughter(MCParticle::kTau) && p->HasDaughter(MCParticle::kTauNu)) {
            isOld = kFALSE;
	    for (UInt_t nl = 0; nl < GenTempMG1->GetEntries(); nl++) {
	      if (p->FindDaughter(MCParticle::kTau) == GenTempMG1->At(nl)) {
		isOld = kTRUE;
		break;
	      }
	    }
	    if (isOld == kFALSE) {
	      GenTempMG1->Add(p->FindDaughter(MCParticle::kTau));
              diBoson.AddDaughter(p->FindDaughter(MCParticle::kTau));
              diBoson.AddDaughter(p->FindDaughter(MCParticle::kTauNu));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 2)
        	hDVVMass[4]->Fill(TMath::Min(diBoson.Mass(),199.999));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 5)
        	hDVVMass[5]->Fill(TMath::Min(diBoson.Mass(),199.999));
	    }
	  }
	  if (p->HasDaughter(MCParticle::kMu,kTRUE) && p->HasDaughter(-1*MCParticle::kMu,kTRUE)) {
            isOld = kFALSE;
	    for (UInt_t nl = 0; nl < GenTempMG1->GetEntries(); nl++) {
	      if (p->FindDaughter(MCParticle::kMu,kTRUE) == GenTempMG1->At(nl)) {
		isOld = kTRUE;
		break;
	      }
	    }
	    if (isOld == kFALSE) {
	      GenTempMG1->Add(p->FindDaughter(MCParticle::kMu,kTRUE));
              diBoson.AddDaughter(p->FindDaughter(MCParticle::kMu,kTRUE));
              diBoson.AddDaughter(p->FindDaughter(-1*MCParticle::kMu,kTRUE));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 5)
        	hDVVMass[6]->Fill(TMath::Min(diBoson.Mass(),199.999));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 8)
        	hDVVMass[7]->Fill(TMath::Min(diBoson.Mass(),199.999));
	    }
	  }
	  if (p->HasDaughter(MCParticle::kEl,kTRUE) && p->HasDaughter(-1*MCParticle::kEl,kTRUE)) {
            isOld = kFALSE;
	    for (UInt_t nl = 0; nl < GenTempMG1->GetEntries(); nl++) {
	      if (p->FindDaughter(MCParticle::kEl,kTRUE) == GenTempMG1->At(nl)) {
		isOld = kTRUE;
		break;
	      }
	    }
	    if (isOld == kFALSE) {
	      GenTempMG1->Add(p->FindDaughter(MCParticle::kEl,kTRUE));
              diBoson.AddDaughter(p->FindDaughter(MCParticle::kEl,kTRUE));
              diBoson.AddDaughter(p->FindDaughter(-1*MCParticle::kEl,kTRUE));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 5)
        	hDVVMass[8]->Fill(TMath::Min(diBoson.Mass(),199.999));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 8)
        	hDVVMass[9]->Fill(TMath::Min(diBoson.Mass(),199.999));
	    }
	  }
	  if (p->HasDaughter(MCParticle::kTau,kTRUE) && p->HasDaughter(-1*MCParticle::kTau,kTRUE)) {
            isOld = kFALSE;
	    for (UInt_t nl = 0; nl < GenTempMG1->GetEntries(); nl++) {
	      if (p->FindDaughter(MCParticle::kTau,kTRUE) == GenTempMG1->At(nl)) {
		isOld = kTRUE;
		break;
	      }
	    }
	    if (isOld == kFALSE) {
	      GenTempMG1->Add(p->FindDaughter(MCParticle::kTau,kTRUE));
              diBoson.AddDaughter(p->FindDaughter(MCParticle::kTau,kTRUE));
              diBoson.AddDaughter(p->FindDaughter(-1*MCParticle::kTau,kTRUE));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 5)
        	hDVVMass[10]->Fill(TMath::Min(diBoson.Mass(),199.999));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 8)
        	hDVVMass[11]->Fill(TMath::Min(diBoson.Mass(),199.999));
	    }
	  }
	  if (p->HasDaughter(MCParticle::kMuNu,kTRUE) && p->HasDaughter(-1*MCParticle::kMuNu,kTRUE)) {
            isOld = kFALSE;
	    for (UInt_t nl = 0; nl < GenTempMG1->GetEntries(); nl++) {
	      if (p->FindDaughter(MCParticle::kMuNu,kTRUE) == GenTempMG1->At(nl)) {
		isOld = kTRUE;
		break;
	      }
	    }
	    if (isOld == kFALSE) {
	      GenTempMG1->Add(p->FindDaughter(MCParticle::kMuNu,kTRUE));
              diBoson.AddDaughter(p->FindDaughter(MCParticle::kMuNu,kTRUE));
              diBoson.AddDaughter(p->FindDaughter(-1*MCParticle::kMuNu,kTRUE));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 5)
        	hDVVMass[12]->Fill(TMath::Min(diBoson.Mass(),199.999));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 8)
        	hDVVMass[13]->Fill(TMath::Min(diBoson.Mass(),199.999));
	    }
	  }
	  if (p->HasDaughter(MCParticle::kElNu,kTRUE) && p->HasDaughter(-1*MCParticle::kElNu,kTRUE)) {
            isOld = kFALSE;
	    for (UInt_t nl = 0; nl < GenTempMG1->GetEntries(); nl++) {
	      if (p->FindDaughter(MCParticle::kElNu,kTRUE) == GenTempMG1->At(nl)) {
		isOld = kTRUE;
		break;
	      }
	    }
	    if (isOld == kFALSE) {
	      GenTempMG1->Add(p->FindDaughter(MCParticle::kElNu,kTRUE));
              diBoson.AddDaughter(p->FindDaughter(MCParticle::kElNu,kTRUE));
              diBoson.AddDaughter(p->FindDaughter(-1*MCParticle::kElNu,kTRUE));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 5)
        	hDVVMass[14]->Fill(TMath::Min(diBoson.Mass(),199.999));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 8)
        	hDVVMass[15]->Fill(TMath::Min(diBoson.Mass(),199.999));
	    }
	  }
	  if (p->HasDaughter(MCParticle::kTauNu,kTRUE) && 
              p->HasDaughter(-1*MCParticle::kTauNu,kTRUE)) {
            isOld = kFALSE;
	    for (UInt_t nl = 0; nl < GenTempMG1->GetEntries(); nl++) {
	      if (p->FindDaughter(MCParticle::kTauNu,kTRUE) == GenTempMG1->At(nl)) {
		isOld = kTRUE;
		break;
	      }
	    }
	    if (isOld == kFALSE) {
	      GenTempMG1->Add(p->FindDaughter(MCParticle::kTauNu,kTRUE));
              diBoson.AddDaughter(p->FindDaughter(MCParticle::kTauNu,kTRUE));
              diBoson.AddDaughter(p->FindDaughter(-1*MCParticle::kTauNu,kTRUE));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 5)
        	hDVVMass[16]->Fill(TMath::Min(diBoson.Mass(),199.999));
              if (GetFillHist() && sumV[0] + 4*sumV[1] == 8)
        	hDVVMass[17]->Fill(TMath::Min(diBoson.Mass(),199.999));
	    }
	  }
	  if      (diBosonMass[0] <= 0) diBosonMass[0] = diBoson.Mass();
	  else if (diBosonMass[1] <= 0) diBosonMass[1] = diBoson.Mass();
	  if (diBoson.Mass() > 0) {
            for (UInt_t i=0; i<GenQuarks->GetEntries(); i++) {
	      if (GenQuarks->At(i)->AbsPdgId() != MCParticle::kBottom) continue;
	      CompositeParticle topCand;
	      topCand.AddDaughter(&diBoson);
	      topCand.AddDaughter(GenQuarks->At(i));
	      if (fFilterBTEvents == kTRUE && topCand.Mass() > 170.0 && topCand.Mass() < 180.0) SkipEvent();
	    }
	  }
	}
	else if (p->Status() == 3 && (p->Is(MCParticle::kZ)  || p->Is(MCParticle::kW))) {
	  if      (diBosonMass[0] <= 0) diBosonMass[0] = p->Mass();
	  else if (diBosonMass[1] <= 0) diBosonMass[1] = p->Mass();
          if (GetFillHist()) {
	    if     (sumV[0] + 4*sumV[1] == 2 && p->Is(MCParticle::kW) && 
                    p->HasDaughter(MCParticle::kMu)  && 
                    p->HasDaughter(MCParticle::kMuNu))
	      hDVVMass[0]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 2 && p->Is(MCParticle::kW) && 
                    p->HasDaughter(MCParticle::kEl)  && 
                    p->HasDaughter(MCParticle::kElNu))
	      hDVVMass[2]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 2 && p->Is(MCParticle::kW) && 
                    p->HasDaughter(MCParticle::kTau) && 
                    p->HasDaughter(MCParticle::kTauNu))
	      hDVVMass[4]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 5 && p->Is(MCParticle::kW) && 
                    p->HasDaughter(MCParticle::kMu)  && 
                    p->HasDaughter(MCParticle::kMuNu))
	      hDVVMass[1]->Fill(p->Mass());
	    else if (sumV[0] + 4*sumV[1] == 5 && p->Is(MCParticle::kW) && 
                    p->HasDaughter(MCParticle::kEl)  && 
                    p->HasDaughter(MCParticle::kElNu))
	      hDVVMass[3]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 5 && p->Is(MCParticle::kW) && 
                    p->HasDaughter(MCParticle::kTau) && 
                    p->HasDaughter(MCParticle::kTauNu))
	      hDVVMass[5]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 5 && p->Is(MCParticle::kZ) && 
                    p->HasDaughter(MCParticle::kMu,kTRUE) && 
                    p->HasDaughter(-1*MCParticle::kMu,kTRUE))
	      hDVVMass[6]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 5 && p->Is(MCParticle::kZ) && 
                    p->HasDaughter(MCParticle::kEl,kTRUE) && 
                    p->HasDaughter(-1*MCParticle::kEl,kTRUE))
	      hDVVMass[8]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 5 && p->Is(MCParticle::kZ) && 
                    p->HasDaughter(MCParticle::kTau,kTRUE) && 
                    p->HasDaughter(-1*MCParticle::kTau,kTRUE))
	      hDVVMass[10]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 5 && p->Is(MCParticle::kZ) && 
                    p->HasDaughter(MCParticle::kMuNu,kTRUE) && 
                    p->HasDaughter(-1*MCParticle::kMuNu,kTRUE))
	      hDVVMass[12]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 5 && p->Is(MCParticle::kZ) && 
                    p->HasDaughter(MCParticle::kElNu,kTRUE) && 
                    p->HasDaughter(-1*MCParticle::kElNu,kTRUE))
	      hDVVMass[14]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 5 && p->Is(MCParticle::kZ) && 
                    p->HasDaughter(MCParticle::kTauNu,kTRUE) && 
                    p->HasDaughter(-1*MCParticle::kTauNu,kTRUE))
	      hDVVMass[16]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 8 && p->Is(MCParticle::kZ) && 
                    p->HasDaughter(MCParticle::kMu,kTRUE) && 
                    p->HasDaughter(-1*MCParticle::kMu,kTRUE))
	      hDVVMass[7]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 8 && p->Is(MCParticle::kZ) && 
                    p->HasDaughter(MCParticle::kEl,kTRUE) && 
                    p->HasDaughter(-1*MCParticle::kEl,kTRUE))
	      hDVVMass[9]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 8 && p->Is(MCParticle::kZ) && 
                    p->HasDaughter(MCParticle::kTau,kTRUE) && 
                    p->HasDaughter(-1*MCParticle::kTau,kTRUE))
	      hDVVMass[11]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 8 && p->Is(MCParticle::kZ) && 
                    p->HasDaughter(MCParticle::kMuNu,kTRUE) && 
                    p->HasDaughter(-1*MCParticle::kMuNu,kTRUE))
	      hDVVMass[13]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 8 && p->Is(MCParticle::kZ) && 
                    p->HasDaughter(MCParticle::kElNu,kTRUE) && 
                    p->HasDaughter(-1*MCParticle::kElNu,kTRUE))
	      hDVVMass[15]->Fill(TMath::Min(p->Mass(),199.999));
	    else if (sumV[0] + 4*sumV[1] == 8 && p->Is(MCParticle::kZ) && 
                    p->HasDaughter(MCParticle::kTauNu,kTRUE) && 
                    p->HasDaughter(-1*MCParticle::kTauNu,kTRUE))
	      hDVVMass[17]->Fill(TMath::Min(p->Mass(),199.999));
	  }
	  if (p->Is(MCParticle::kW)) {
            for (UInt_t i=0; i<GenQuarks->GetEntries(); i++) {
	      if (GenQuarks->At(i)->AbsPdgId() != MCParticle::kBottom) continue;
	      CompositeParticle topCand;
	      topCand.AddDaughter(p);
	      topCand.AddDaughter(GenQuarks->At(i));
	      if (fFilterBTEvents == kTRUE && topCand.Mass() > 170.0 && topCand.Mass() < 180.0) SkipEvent();
	    }
	  }
	}
      } // end loop of particles
      if (GetFillHist()) {
	if (diBosonMass[0] > 70 && diBosonMass[0] < 110 && 
            diBosonMass[1] > 70 && diBosonMass[1] < 110) {
          if (sumV[0] + 4*sumV[1] == 2) {
            if      (sumVVFlavor[0] == 2)    		        hDVVMass[18]->Fill(0.);
            else if (sumVVFlavor[1] == 2)    		        hDVVMass[18]->Fill(1.);
            else if (sumVVFlavor[2] == 2)                        hDVVMass[18]->Fill(2.);
            else if (sumVVFlavor[0] == 1 && sumVVFlavor[1] == 1) hDVVMass[18]->Fill(3.);
            else if (sumVVFlavor[0] == 1 && sumVVFlavor[2] == 1) hDVVMass[18]->Fill(4.);
            else if (sumVVFlavor[1] == 1 && sumVVFlavor[2] == 1) hDVVMass[18]->Fill(5.);
            else                                                hDVVMass[18]->Fill(6.);
          }
          if (sumV[0] + 4*sumV[1] == 5) {
            if      (sumVVFlavor[3] == 1 && sumVVFlavor[0] == 1)  hDVVMass[19]->Fill(0.);
            else if (sumVVFlavor[3] == 1 && sumVVFlavor[1] == 1)  hDVVMass[19]->Fill(1.);
            else if (sumVVFlavor[3] == 1 && sumVVFlavor[2] == 1)  hDVVMass[19]->Fill(2.);
            else if (sumVVFlavor[4] == 1 && sumVVFlavor[0] == 1)  hDVVMass[19]->Fill(3.);
            else if (sumVVFlavor[4] == 1 && sumVVFlavor[1] == 1)  hDVVMass[19]->Fill(4.);
            else if (sumVVFlavor[4] == 1 && sumVVFlavor[2] == 1)  hDVVMass[19]->Fill(5.);
            else if (sumVVFlavor[5] == 1 && sumVVFlavor[0] == 1)  hDVVMass[19]->Fill(6.);
            else if (sumVVFlavor[5] == 1 && sumVVFlavor[1] == 1)  hDVVMass[19]->Fill(7.);
            else if (sumVVFlavor[5] == 1 && sumVVFlavor[2] == 1)  hDVVMass[19]->Fill(8.);
            else  {                         			 hDVVMass[19]->Fill(9.);
	    /*for (int i=0; i<9; i++) cout << sumVVFlavor[i] << " " ;cout << endl;*/}
          }
          if (sumV[0] + 4*sumV[1] == 8 &&
             sumVVFlavor[3] + sumVVFlavor[4] + sumVVFlavor[5] == 2) {
            if      (sumVVFlavor[3] == 2)  			 hDVVMass[20]->Fill(0.);
            else if (sumVVFlavor[4] == 2)  			 hDVVMass[20]->Fill(1.);
            else if (sumVVFlavor[5] == 2)  			 hDVVMass[20]->Fill(2.);
            else if (sumVVFlavor[3] == 1 && sumVVFlavor[4] == 1)  hDVVMass[20]->Fill(3.);
            else if (sumVVFlavor[3] == 1 && sumVVFlavor[5] == 1)  hDVVMass[20]->Fill(4.);
            else if (sumVVFlavor[4] == 1 && sumVVFlavor[5] == 1)  hDVVMass[20]->Fill(5.);
            else                                                 hDVVMass[20]->Fill(6.);
          }
          else if (sumV[0] + 4*sumV[1] == 8) {
            if      (sumVVFlavor[6] == 2)  			 hDVVMass[21]->Fill(0.);
            else if (sumVVFlavor[7] == 2)  			 hDVVMass[21]->Fill(1.);
            else if (sumVVFlavor[8] == 2)  			 hDVVMass[21]->Fill(2.);
            else if (sumVVFlavor[3] == 1 && sumVVFlavor[6] == 1)  hDVVMass[21]->Fill(3.);
            else if (sumVVFlavor[3] == 1 && sumVVFlavor[7] == 1)  hDVVMass[21]->Fill(4.);
            else if (sumVVFlavor[3] == 1 && sumVVFlavor[8] == 1)  hDVVMass[21]->Fill(5.);
            else if (sumVVFlavor[4] == 1 && sumVVFlavor[6] == 1)  hDVVMass[21]->Fill(6.);
            else if (sumVVFlavor[4] == 1 && sumVVFlavor[7] == 1)  hDVVMass[21]->Fill(7.);
            else if (sumVVFlavor[4] == 1 && sumVVFlavor[8] == 1)  hDVVMass[21]->Fill(8.);
            else if (sumVVFlavor[5] == 1 && sumVVFlavor[6] == 1)  hDVVMass[21]->Fill(9.);
            else if (sumVVFlavor[5] == 1 && sumVVFlavor[7] == 1)  hDVVMass[21]->Fill(10.);
            else if (sumVVFlavor[5] == 1 && sumVVFlavor[8] == 1)  hDVVMass[21]->Fill(11.);
            else if (sumVVFlavor[6] == 1 && sumVVFlavor[7] == 1)  hDVVMass[21]->Fill(12.);
            else if (sumVVFlavor[6] == 1 && sumVVFlavor[8] == 1)  hDVVMass[21]->Fill(13.);
            else if (sumVVFlavor[7] == 1 && sumVVFlavor[8] == 1)  hDVVMass[21]->Fill(14.);
            else                                                 hDVVMass[21]->Fill(15.);
          }
	} // 60<mV1/2<120
	if (sumV[0] + 4*sumV[1] == 2) 
          hDVVMass[22]->Fill(TMath::Min(TMath::Min(diBosonMass[0],diBosonMass[1]),199.999));
	if (sumV[0] + 4*sumV[1] == 2) 
          hDVVMass[23]->Fill(TMath::Min(TMath::Max(diBosonMass[0],diBosonMass[1]),199.999));
	if (sumV[0] + 4*sumV[1] == 5) 
          hDVVMass[24]->Fill(TMath::Min(TMath::Min(diBosonMass[0],diBosonMass[1]),199.999));
	if (sumV[0] + 4*sumV[1] == 5) 
          hDVVMass[25]->Fill(TMath::Min(TMath::Max(diBosonMass[0],diBosonMass[1]),199.999));
	if (sumV[0] + 4*sumV[1] == 8) 
          hDVVMass[26]->Fill(TMath::Min(TMath::Min(diBosonMass[0],diBosonMass[1]),199.999));
	if (sumV[0] + 4*sumV[1] == 8) 
          hDVVMass[27]->Fill(TMath::Min(TMath::Max(diBosonMass[0],diBosonMass[1]),199.999));
      }
      delete GenTempMG1;
    } // WW, WZ or ZZ
    // --------------------------------
    // End special study about VVjets
    // --------------------------------

    Met *theMET = new Met(-totalMET[0], -totalMET[1]);
    theMET->SetElongitudinal(-totalMET[2]);
    GenMet->AddOwned(theMET);

    // sort according to pt
    GenLeptons->Sort();
    GenAllLeptons->Sort();
    GenTaus->Sort();
    GenNeutrinos->Sort();
    GenQuarks->Sort();
    GenqqHs->Sort();
    GenBosons->Sort();
    GenPhotons->Sort();
    GenRadPhotons->Sort();
    GenISRPhotons->Sort();
  } // Only for Monte Carlo

  // add objects to this event for other modules to use
  AddObjThisEvt(GenMet);  
  AddObjThisEvt(GenLeptons);  
  AddObjThisEvt(GenAllLeptons);
  AddObjThisEvt(GenTaus);
  AddObjThisEvt(GenNeutrinos);
  AddObjThisEvt(GenQuarks);
  AddObjThisEvt(GenqqHs);
  AddObjThisEvt(GenBosons);
  AddObjThisEvt(GenPhotons);
  AddObjThisEvt(GenRadPhotons);
  AddObjThisEvt(GenISRPhotons);

  // --------------------------------
  // Copy these Collections into the Arrays for Publication for Output Module
  // --------------------------------
  fGenLeptons->Delete();
  fGenAllLeptons->Delete();
  fGenTaus->Delete();
  fGenNeutrinos->Delete();
  fGenQuarks->Delete();
  fGenqqHs->Delete();
  fGenBosons->Delete();
  fGenPhotons->Delete();
  fGenRadPhotons->Delete();
  fGenISRPhotons->Delete();

  for (UInt_t i=0; i < GenLeptons->GetEntries(); ++i) {
    mithep::MCParticle *genParticle = fGenLeptons->AddNew();
    genParticle->SetPtEtaPhiM(GenLeptons->At(i)->Pt(), GenLeptons->At(i)->Eta(), GenLeptons->At(i)->Phi(), GenLeptons->At(i)->Mass());
    genParticle->SetPdgId(GenLeptons->At(i)->PdgId());
    genParticle->SetStatus(GenLeptons->At(i)->Status());
  }
  for (UInt_t i=0; i < GenAllLeptons->GetEntries(); ++i) {
    mithep::MCParticle *genParticle = fGenAllLeptons->AddNew();
    genParticle->SetPtEtaPhiM(GenAllLeptons->At(i)->Pt(), GenAllLeptons->At(i)->Eta(), GenAllLeptons->At(i)->Phi(), GenAllLeptons->At(i)->Mass());
    genParticle->SetPdgId(GenAllLeptons->At(i)->PdgId());
    genParticle->SetStatus(GenAllLeptons->At(i)->Status());
  }
  for (UInt_t i=0; i < GenTaus->GetEntries(); ++i) {
    mithep::MCParticle *genParticle = fGenTaus->AddNew();
    genParticle->SetPtEtaPhiM(GenTaus->At(i)->Pt(), GenTaus->At(i)->Eta(), GenTaus->At(i)->Phi(), GenTaus->At(i)->Mass());
    genParticle->SetPdgId(GenTaus->At(i)->PdgId());
    genParticle->SetStatus(GenTaus->At(i)->Status());
  }
  for (UInt_t i=0; i < GenNeutrinos->GetEntries(); ++i) {
    mithep::MCParticle *genParticle = fGenNeutrinos->AddNew();
    genParticle->SetPtEtaPhiM(GenNeutrinos->At(i)->Pt(), GenNeutrinos->At(i)->Eta(), GenNeutrinos->At(i)->Phi(), GenNeutrinos->At(i)->Mass());
    genParticle->SetPdgId(GenNeutrinos->At(i)->PdgId());
    genParticle->SetStatus(GenNeutrinos->At(i)->Status());
  }
  for (UInt_t i=0; i < GenQuarks->GetEntries(); ++i) {
    mithep::MCParticle *genParticle = fGenQuarks->AddNew();
    genParticle->SetPtEtaPhiM(GenQuarks->At(i)->Pt(), GenQuarks->At(i)->Eta(), GenQuarks->At(i)->Phi(), GenQuarks->At(i)->Mass());
    genParticle->SetPdgId(GenQuarks->At(i)->PdgId());
    genParticle->SetStatus(GenQuarks->At(i)->Status());
  }
  for (UInt_t i=0; i < GenqqHs->GetEntries(); ++i) {
    mithep::MCParticle *genParticle = fGenqqHs->AddNew();
    genParticle->SetPtEtaPhiM(GenqqHs->At(i)->Pt(), GenqqHs->At(i)->Eta(), GenqqHs->At(i)->Phi(), GenqqHs->At(i)->Mass());
    genParticle->SetPdgId(GenqqHs->At(i)->PdgId());
    genParticle->SetStatus(GenqqHs->At(i)->Status());
  }

  if (fCopyArrays) {
    // --------------------------------
    // Copy these Collections into the Arrays for Publication for Output Module
    // --------------------------------
    fGenLeptons->Delete();
    fGenAllLeptons->Delete();
    fGenTaus->Delete();
    fGenNeutrinos->Delete();
    fGenQuarks->Delete();
    fGenqqHs->Delete();
    fGenBosons->Delete();
    fGenPhotons->Delete();
    fGenRadPhotons->Delete();
    fGenISRPhotons->Delete();

    for (UInt_t i=0; i < GenLeptons->GetEntries(); ++i) {
      mithep::MCParticle *genParticle = fGenLeptons->AddNew();
      genParticle->SetPtEtaPhiM(GenLeptons->At(i)->Pt(), GenLeptons->At(i)->Eta(), GenLeptons->At(i)->Phi(), GenLeptons->At(i)->Mass());
      genParticle->SetPdgId(GenLeptons->At(i)->PdgId());
      genParticle->SetStatus(GenLeptons->At(i)->Status());
    }
    for (UInt_t i=0; i < GenAllLeptons->GetEntries(); ++i) {
      mithep::MCParticle *genParticle = fGenAllLeptons->AddNew();
      genParticle->SetPtEtaPhiM(GenAllLeptons->At(i)->Pt(), GenAllLeptons->At(i)->Eta(), GenAllLeptons->At(i)->Phi(), GenAllLeptons->At(i)->Mass());
      genParticle->SetPdgId(GenAllLeptons->At(i)->PdgId());
      genParticle->SetStatus(GenAllLeptons->At(i)->Status());
    }
    for (UInt_t i=0; i < GenTaus->GetEntries(); ++i) {
      mithep::MCParticle *genParticle = fGenTaus->AddNew();
      genParticle->SetPtEtaPhiM(GenTaus->At(i)->Pt(), GenTaus->At(i)->Eta(), GenTaus->At(i)->Phi(), GenTaus->At(i)->Mass());
      genParticle->SetPdgId(GenTaus->At(i)->PdgId());
      genParticle->SetStatus(GenTaus->At(i)->Status());
    }
    for (UInt_t i=0; i < GenNeutrinos->GetEntries(); ++i) {
      mithep::MCParticle *genParticle = fGenNeutrinos->AddNew();
      genParticle->SetPtEtaPhiM(GenNeutrinos->At(i)->Pt(), GenNeutrinos->At(i)->Eta(), GenNeutrinos->At(i)->Phi(), GenNeutrinos->At(i)->Mass());
      genParticle->SetPdgId(GenNeutrinos->At(i)->PdgId());
      genParticle->SetStatus(GenNeutrinos->At(i)->Status());
    }
    for (UInt_t i=0; i < GenQuarks->GetEntries(); ++i) {
      mithep::MCParticle *genParticle = fGenQuarks->AddNew();
      genParticle->SetPtEtaPhiM(GenQuarks->At(i)->Pt(), GenQuarks->At(i)->Eta(), GenQuarks->At(i)->Phi(), GenQuarks->At(i)->Mass());
      genParticle->SetPdgId(GenQuarks->At(i)->PdgId());
      genParticle->SetStatus(GenQuarks->At(i)->Status());
    }
    for (UInt_t i=0; i < GenqqHs->GetEntries(); ++i) {
      mithep::MCParticle *genParticle = fGenqqHs->AddNew();
      genParticle->SetPtEtaPhiM(GenqqHs->At(i)->Pt(), GenqqHs->At(i)->Eta(), GenqqHs->At(i)->Phi(), GenqqHs->At(i)->Mass());
      genParticle->SetPdgId(GenqqHs->At(i)->PdgId());
      genParticle->SetStatus(GenqqHs->At(i)->Status());
    }
    for (UInt_t i=0; i < GenBosons->GetEntries(); ++i) {
      mithep::MCParticle *genParticle = fGenBosons->AddNew();
      genParticle->SetPtEtaPhiM(GenBosons->At(i)->Pt(), GenBosons->At(i)->Eta(), GenBosons->At(i)->Phi(), GenBosons->At(i)->Mass());
      genParticle->SetPdgId(GenBosons->At(i)->PdgId());
      genParticle->SetStatus(GenBosons->At(i)->Status());
    }
    for (UInt_t i=0; i < GenPhotons->GetEntries(); ++i) {
      mithep::MCParticle *genParticle = fGenPhotons->AddNew();
      genParticle->SetPtEtaPhiM(GenPhotons->At(i)->Pt(), GenPhotons->At(i)->Eta(), GenPhotons->At(i)->Phi(), GenPhotons->At(i)->Mass());
      genParticle->SetPdgId(GenPhotons->At(i)->PdgId());
      genParticle->SetStatus(GenPhotons->At(i)->Status());
    }
    for (UInt_t i=0; i < GenRadPhotons->GetEntries(); ++i) {
      mithep::MCParticle *genParticle = fGenRadPhotons->AddNew();
      genParticle->SetPtEtaPhiM(GenRadPhotons->At(i)->Pt(), GenRadPhotons->At(i)->Eta(), GenRadPhotons->At(i)->Phi(), GenRadPhotons->At(i)->Mass());
      genParticle->SetPdgId(GenRadPhotons->At(i)->PdgId());
      genParticle->SetStatus(GenRadPhotons->At(i)->Status());
    }
    for (UInt_t i=0; i < GenISRPhotons->GetEntries(); ++i) {
      mithep::MCParticle *genParticle = fGenISRPhotons->AddNew();
      genParticle->SetPtEtaPhiM(GenISRPhotons->At(i)->Pt(), GenISRPhotons->At(i)->Eta(), GenISRPhotons->At(i)->Phi(), GenISRPhotons->At(i)->Mass());
      genParticle->SetPdgId(GenISRPhotons->At(i)->PdgId());
      genParticle->SetStatus(GenISRPhotons->At(i)->Status());
    }
  }

  for (UInt_t i=0; i < GenBosons->GetEntries(); ++i) {
    mithep::MCParticle *genParticle = fGenBosons->AddNew();
    genParticle->SetPtEtaPhiM(GenBosons->At(i)->Pt(), GenBosons->At(i)->Eta(), GenBosons->At(i)->Phi(), GenBosons->At(i)->Mass());
    genParticle->SetPdgId(GenBosons->At(i)->PdgId());
    genParticle->SetStatus(GenBosons->At(i)->Status());
  }
  for (UInt_t i=0; i < GenPhotons->GetEntries(); ++i) {
    mithep::MCParticle *genParticle = fGenPhotons->AddNew();
    genParticle->SetPtEtaPhiM(GenPhotons->At(i)->Pt(), GenPhotons->At(i)->Eta(), GenPhotons->At(i)->Phi(), GenPhotons->At(i)->Mass());
    genParticle->SetPdgId(GenPhotons->At(i)->PdgId());
    genParticle->SetStatus(GenPhotons->At(i)->Status());
  }
  for (UInt_t i=0; i < GenRadPhotons->GetEntries(); ++i) {
    mithep::MCParticle *genParticle = fGenRadPhotons->AddNew();
    genParticle->SetPtEtaPhiM(GenRadPhotons->At(i)->Pt(), GenRadPhotons->At(i)->Eta(), GenRadPhotons->At(i)->Phi(), GenRadPhotons->At(i)->Mass());
    genParticle->SetPdgId(GenRadPhotons->At(i)->PdgId());
    genParticle->SetStatus(GenRadPhotons->At(i)->Status());
  }
  for (UInt_t i=0; i < GenISRPhotons->GetEntries(); ++i) {
    mithep::MCParticle *genParticle = fGenISRPhotons->AddNew();
    genParticle->SetPtEtaPhiM(GenISRPhotons->At(i)->Pt(), GenISRPhotons->At(i)->Eta(), GenISRPhotons->At(i)->Phi(), GenISRPhotons->At(i)->Mass());
    genParticle->SetPdgId(GenISRPhotons->At(i)->PdgId());
    genParticle->SetStatus(GenISRPhotons->At(i)->Status());
  }

  // Apply WW filter (without filling all histograms)
  Bool_t passVVFilter = kFALSE;
  if ((fAllowWWEvents == kTRUE && sumV[0] + 4*sumV[1] == 2 )  ||
      (fAllowWZEvents == kTRUE && sumV[0] + 4*sumV[1] == 6 )  ||
       (fAllowZZEvents == kTRUE && sumV[0] + 4*sumV[1] == 8 )   )
    passVVFilter = kTRUE;
  if (fApplyVVFilter && !passVVFilter) {
    SkipEvent();
    return;
  }

  // fill histograms if requested
  if (GetFillHist()) {
    // MET
    hDGenMet[0]->Fill(GenMet->At(0)->Pt());
    hDGenMet[1]->Fill(GenMet->At(0)->Px());
    hDGenMet[2]->Fill(GenMet->At(0)->Py());
    hDGenMet[3]->Fill(GenMet->At(0)->Elongitudinal());

    // pt min for leptons from W/Z
    hDGenPtMin->Fill(TMath::Min(ptMin, 199.999));

    //if (GenLeptons->GetEntries() != 3) {SkipEvent(); return;}
    //if (GenLeptons->At(0)->Pt()     <= 20.0 || GenLeptons->At(1)->Pt()     <= 10.0 || GenLeptons->At(2)->Pt()     <= 10.0 ||
    //   GenLeptons->At(0)->AbsEta() >=  2.5 || GenLeptons->At(1)->AbsEta() >=  2.5 || GenLeptons->At(2)->AbsEta() >=  2.5) {SkipEvent(); return;}
  
    // leptons
    hDGenLeptons[0]->Fill(GenLeptons->GetEntries());
    for (UInt_t i=0; i<GenLeptons->GetEntries(); i++) {
      hDGenLeptons[1]->Fill(GenLeptons->At(i)->Pt());
      hDGenLeptons[2]->Fill(TMath::Abs(GenLeptons->At(i)->Eta()));
      hDGenLeptons[3]->Fill(GenLeptons->At(i)->PhiDeg());
      for (UInt_t j=i+1; j<GenLeptons->GetEntries(); j++) {
        CompositeParticle *dilepton = new CompositeParticle();
        dilepton->AddDaughter(GenLeptons->At(i));
        dilepton->AddDaughter(GenLeptons->At(j));
        hDGenLeptons[4]->Fill(dilepton->Mass());
	delete dilepton;
      }
    }
    // looking at events with two leptons
    if (GenLeptons->GetEntries() == 2) {
      hDGenLeptons[5]->Fill(TMath::Min(TMath::Max(TMath::Abs(GenLeptons->At(0)->Eta()),
                                                  TMath::Abs(GenLeptons->At(1)->Eta())),
                                       4.999));
      hDGenLeptons[6]->Fill(TMath::Min(TMath::Min(TMath::Abs(GenLeptons->At(0)->Eta()),
                                                  TMath::Abs(GenLeptons->At(1)->Eta())),
                                       4.999));
      if (TMath::Abs(GenLeptons->At(0)->Eta()) < 2.5 &&
          TMath::Abs(GenLeptons->At(1)->Eta()) < 2.5) {
        hDGenLeptons[7]->Fill(TMath::Min(GenLeptons->At(0)->Pt(),199.999));
        if (GenLeptons->At(0)->Pt() > 20.0) {
          hDGenLeptons[8]->Fill(TMath::Min(GenLeptons->At(1)->Pt(),199.999));
          if (GenLeptons->At(1)->Pt() > 10.0) {
            CompositeParticle *dilepton = new CompositeParticle();
            dilepton->AddDaughter(GenLeptons->At(0));
            dilepton->AddDaughter(GenLeptons->At(1));
            hDGenLeptons[9]->Fill(TMath::Min(dilepton->Mass(),999.999));
	    hDGenLeptons[10]->Fill(MathUtils::DeltaPhi(GenLeptons->At(0)->Phi(),
	    					       GenLeptons->At(1)->Phi())
	        				       * 180./ TMath::Pi());
	    hDGenLeptons[11]->Fill(MathUtils::DeltaR(*GenLeptons->At(0), 
            					     *GenLeptons->At(1)));
	    if (dilepton->Mass() > 12.0) {
	      hDGenLeptons[12]->Fill(MathUtils::DeltaPhi(GenLeptons->At(0)->Phi(),
	                                                 GenLeptons->At(1)->Phi())
						         * 180./ TMath::Pi());
	      hDGenLeptons[13]->Fill(MathUtils::DeltaR(*GenLeptons->At(0), 
                                                       *GenLeptons->At(1)));
	    }
	    delete dilepton;
	  }
	}
      }
    }
    // looking at events with three leptons
    if (GenLeptons->GetEntries() == 3) {
      if (TMath::Abs(GenLeptons->At(0)->Eta()) < 2.5 &&
          TMath::Abs(GenLeptons->At(1)->Eta()) < 2.5 &&
          TMath::Abs(GenLeptons->At(2)->Eta()) < 2.5) {
        hDGenLeptons[14]->Fill(TMath::Min(GenLeptons->At(0)->Pt(),199.999));
        if (GenLeptons->At(0)->Pt() > 20.0) {
          hDGenLeptons[15]->Fill(TMath::Min(GenLeptons->At(1)->Pt(),199.999));
          hDGenLeptons[16]->Fill(TMath::Min(GenLeptons->At(2)->Pt(),199.999));
          if (GenLeptons->At(1)->Pt() > 10.0 && GenLeptons->At(2)->Pt() > 10.0) {
            CompositeParticle *dilepton01 = new CompositeParticle();
            dilepton01->AddDaughter(GenLeptons->At(0));
            dilepton01->AddDaughter(GenLeptons->At(1));
            CompositeParticle *dilepton02 = new CompositeParticle();
            dilepton02->AddDaughter(GenLeptons->At(0));
            dilepton02->AddDaughter(GenLeptons->At(2));
            CompositeParticle *dilepton12 = new CompositeParticle();
            dilepton12->AddDaughter(GenLeptons->At(1));
            dilepton12->AddDaughter(GenLeptons->At(2));
            hDGenLeptons[17]->Fill(TMath::Min(dilepton01->Mass(),999.999));
            hDGenLeptons[17]->Fill(TMath::Min(dilepton02->Mass(),999.999));
            hDGenLeptons[17]->Fill(TMath::Min(dilepton12->Mass(),999.999));
            CompositeParticle *trilepton = new CompositeParticle();
            trilepton->AddDaughter(GenLeptons->At(0));
            trilepton->AddDaughter(GenLeptons->At(1));
            trilepton->AddDaughter(GenLeptons->At(2));
            hDGenLeptons[18]->Fill(TMath::Min(trilepton->Mass(),999.999));
	    Double_t deltaR[3] = {MathUtils::DeltaR(*GenLeptons->At(0),
                                                    *GenLeptons->At(1)),
                                  MathUtils::DeltaR(*GenLeptons->At(0), 
                                                    *GenLeptons->At(2)), 
                                  MathUtils::DeltaR(*GenLeptons->At(1), 
                                                    *GenLeptons->At(2))};
	    Double_t deltaRMin = deltaR[0];
            for (Int_t i=1; i<3; i++) 
              if (deltaRMin > deltaR[i]) 
                deltaRMin = deltaR[i];
            hDGenLeptons[19]->Fill(deltaRMin);

	    delete dilepton01;
	    delete dilepton02;
	    delete dilepton12;
	    delete trilepton;
	  }
	}
      }
    }
    // looking at events with four leptons
    if (GenLeptons->GetEntries() == 4) {
      if (TMath::Abs(GenLeptons->At(0)->Eta()) < 2.5 &&
          TMath::Abs(GenLeptons->At(1)->Eta()) < 2.5 &&
          TMath::Abs(GenLeptons->At(2)->Eta()) < 2.5 &&
          TMath::Abs(GenLeptons->At(3)->Eta()) < 2.5) {
        hDGenLeptons[20]->Fill(TMath::Min(GenLeptons->At(0)->Pt(),199.999));
        if (GenLeptons->At(0)->Pt() > 20.0) {
          hDGenLeptons[21]->Fill(TMath::Min(GenLeptons->At(1)->Pt(),199.999));
          hDGenLeptons[22]->Fill(TMath::Min(GenLeptons->At(2)->Pt(),199.999));
          hDGenLeptons[23]->Fill(TMath::Min(GenLeptons->At(3)->Pt(),199.999));
          if (GenLeptons->At(1)->Pt() > 10.0 && GenLeptons->At(2)->Pt() > 10.0 &&
	      GenLeptons->At(3)->Pt() > 10.0) {
            CompositeParticle *dilepton01 = new CompositeParticle();
            dilepton01->AddDaughter(GenLeptons->At(0));
            dilepton01->AddDaughter(GenLeptons->At(1));
            CompositeParticle *dilepton02 = new CompositeParticle();
            dilepton02->AddDaughter(GenLeptons->At(0));
            dilepton02->AddDaughter(GenLeptons->At(2));
            CompositeParticle *dilepton03 = new CompositeParticle();
            dilepton03->AddDaughter(GenLeptons->At(0));
            dilepton03->AddDaughter(GenLeptons->At(3));
            CompositeParticle *dilepton12 = new CompositeParticle();
            dilepton12->AddDaughter(GenLeptons->At(1));
            dilepton12->AddDaughter(GenLeptons->At(2));
            CompositeParticle *dilepton13 = new CompositeParticle();
            dilepton13->AddDaughter(GenLeptons->At(1));
            dilepton13->AddDaughter(GenLeptons->At(3));
            CompositeParticle *dilepton23 = new CompositeParticle();
            dilepton23->AddDaughter(GenLeptons->At(2));
            dilepton23->AddDaughter(GenLeptons->At(3));
            hDGenLeptons[24]->Fill(TMath::Min(dilepton01->Mass(),999.999));
            hDGenLeptons[24]->Fill(TMath::Min(dilepton02->Mass(),999.999));
            hDGenLeptons[24]->Fill(TMath::Min(dilepton03->Mass(),999.999));
            hDGenLeptons[24]->Fill(TMath::Min(dilepton12->Mass(),999.999));
            hDGenLeptons[24]->Fill(TMath::Min(dilepton13->Mass(),999.999));
            hDGenLeptons[24]->Fill(TMath::Min(dilepton23->Mass(),999.999));
            CompositeParticle *fourlepton = new CompositeParticle();
            fourlepton->AddDaughter(GenLeptons->At(0));
            fourlepton->AddDaughter(GenLeptons->At(1));
            fourlepton->AddDaughter(GenLeptons->At(2));
            fourlepton->AddDaughter(GenLeptons->At(3));
            hDGenLeptons[25]->Fill(TMath::Min(fourlepton->Mass(),999.999));
	    Double_t deltaR[6] = {MathUtils::DeltaR(*GenLeptons->At(0),
                                                    *GenLeptons->At(1)),
                                  MathUtils::DeltaR(*GenLeptons->At(0), 
                                                    *GenLeptons->At(2)),
                                  MathUtils::DeltaR(*GenLeptons->At(0), 
                                                    *GenLeptons->At(3)), 
                                  MathUtils::DeltaR(*GenLeptons->At(1), 
                                                    *GenLeptons->At(2)),
                                  MathUtils::DeltaR(*GenLeptons->At(1), 
                                                    *GenLeptons->At(3)),
                                  MathUtils::DeltaR(*GenLeptons->At(2), 
                                                    *GenLeptons->At(3))};
	    Double_t deltaRMin = deltaR[0];
            for (Int_t i=1; i<6; i++) 
              if (deltaRMin > deltaR[i]) 
                deltaRMin = deltaR[i];
            hDGenLeptons[26]->Fill(deltaRMin);

	    delete dilepton01;
	    delete dilepton02;
	    delete dilepton03;
	    delete dilepton12;
	    delete dilepton13;
	    delete dilepton23;
	    delete fourlepton;
	  }
	}
      }
    }

    // all leptons
    hDGenAllLeptons[0]->Fill(GenAllLeptons->GetEntries());
    for (UInt_t i=0; i<GenAllLeptons->GetEntries(); i++) {
      hDGenAllLeptons[1]->Fill(GenAllLeptons->At(i)->Pt());
      hDGenAllLeptons[2]->Fill(GenAllLeptons->At(i)->AbsEta());
      hDGenAllLeptons[3]->Fill(GenAllLeptons->At(i)->PhiDeg());
    }
    if (GenAllLeptons->GetEntries() >= 2) hDGenAllLeptons[4]->Fill(GenAllLeptons->At(1)->Pt());

    // taus
    hDGenTaus[0]->Fill(GenTaus->GetEntries());
    for (UInt_t i=0; i<GenTaus->GetEntries(); i++) {
      hDGenTaus[1]->Fill(GenTaus->At(i)->Pt());
      hDGenTaus[2]->Fill(GenTaus->At(i)->AbsEta());
      hDGenTaus[3]->Fill(GenTaus->At(i)->PhiDeg());
    }

    // neutrinos
    hDGenNeutrinos[0]->Fill(GenNeutrinos->GetEntries());
    CompositeParticle *neutrinoTotal = new CompositeParticle();
    for (UInt_t i=0; i<GenNeutrinos->GetEntries(); i++) {
      if (GenNeutrinos->At(i)->HasMother())
        neutrinoTotal->AddDaughter(GenNeutrinos->At(i));
    }
    if (GenNeutrinos->GetEntries() > 0) {
      hDGenNeutrinos[1]->Fill(neutrinoTotal->Pt());
      hDGenNeutrinos[2]->Fill(neutrinoTotal->AbsEta());
      hDGenNeutrinos[3]->Fill(neutrinoTotal->PhiDeg());    
    }
    delete neutrinoTotal;
    
    // quarks
    int nBQuarks=0; int nTQuarks=0;
    hDGenQuarks[0]->Fill(GenQuarks->GetEntries());
    for (UInt_t i=0; i<GenQuarks->GetEntries(); i++) {
      for (UInt_t j=i+1; j<GenQuarks->GetEntries(); j++) {
        CompositeParticle *dijet = new CompositeParticle();
        dijet->AddDaughter(GenQuarks->At(i));
        dijet->AddDaughter(GenQuarks->At(j));
        hDGenQuarks[1]->Fill(dijet->Pt());
        hDGenQuarks[2]->Fill(dijet->Mass());
	if (TMath::Abs(GenQuarks->At(i)->Eta()) < 2.5 && 
            TMath::Abs(GenQuarks->At(j)->Eta()) < 2.5) {
          hDGenQuarks[3]->Fill(dijet->Pt());
          hDGenQuarks[4]->Fill(dijet->Mass());
	}
	delete dijet;
      }
      // b quark info
      if    (GenQuarks->At(i)->AbsPdgId() == MCParticle::kBottom && !GenQuarks->At(i)->FindMother(MCParticle::kTop) && 
                                                                    !GenQuarks->At(i)->FindMother(MCParticle::kBottom)) {
        nBQuarks++;
	hDGenQuarks[5]->Fill(GenQuarks->At(i)->Pt());
        hDGenQuarks[6]->Fill(GenQuarks->At(i)->Eta());      
        hDGenQuarks[7]->Fill(GenQuarks->At(i)->Phi());      
        if (GenLeptons->GetEntries() >= 2 && 
	    GenLeptons->At(0)->Pt() > 20  &&
	    GenLeptons->At(1)->Pt() > 15) {
          if (TMath::Abs(GenLeptons->At(0)->Eta()) < 2.5 &&
              TMath::Abs(GenLeptons->At(1)->Eta()) < 2.5) {
            hDGenQuarks[8]->Fill(GenQuarks->At(i)->Pt());	
            hDGenQuarks[9]->Fill(GenQuarks->At(i)->Eta());      
            hDGenQuarks[10]->Fill(GenQuarks->At(i)->Phi());	
	  }
	}
      }
      // t quark info
      else if (GenQuarks->At(i)->AbsPdgId() == MCParticle::kTop && !GenQuarks->At(i)->FindMother(MCParticle::kTop)) {
        nTQuarks++;
        hDGenQuarks[11]->Fill(GenQuarks->At(i)->Pt());
        hDGenQuarks[12]->Fill(GenQuarks->At(i)->Eta());      
        hDGenQuarks[13]->Fill(GenQuarks->At(i)->Phi());      
      }
      // light quark info
      else {
        hDGenQuarks[14]->Fill(GenQuarks->At(i)->Pt());
        hDGenQuarks[15]->Fill(GenQuarks->At(i)->Eta());      
        hDGenQuarks[16]->Fill(GenQuarks->At(i)->Phi());      
      }
    }

    // filter events with either b or top quarks
    if (fFilterBTEvents == kTRUE && nTQuarks > 0) SkipEvent();

    // wbf
    if (GenqqHs->GetEntries() == 2) {
      hDGenWBF[0]->Fill(MathUtils::DeltaPhi(GenqqHs->At(0)->Phi(),
    					    GenqqHs->At(1)->Phi()) * 180./ TMath::Pi());
      hDGenWBF[1]->Fill(TMath::Abs(GenqqHs->At(0)->Eta()-GenqqHs->At(1)->Eta()));
      hDGenWBF[2]->Fill(TMath::Max(GenqqHs->At(0)->Pt(),GenqqHs->At(1)->Pt()));
      hDGenWBF[3]->Fill(TMath::Min(GenqqHs->At(0)->Pt(),GenqqHs->At(1)->Pt()));
      CompositeParticle *diqq = new CompositeParticle();
      diqq->AddDaughter(GenqqHs->At(0));
      diqq->AddDaughter(GenqqHs->At(1));
      hDGenWBF[4]->Fill(diqq->Mass());
      delete diqq;
    }

    // bosons
    hDGenBosons[0]->Fill(GenBosons->GetEntries());
    for (UInt_t i=0; i<GenBosons->GetEntries(); i++) {
      hDGenBosons[1]->Fill(GenBosons->At(i)->Pt());
      hDGenBosons[2]->Fill(GenBosons->At(i)->Eta());
      hDGenBosons[3]->Fill(TMath::Min(GenBosons->At(i)->Mass(),1999.999));
      hDGenBosons[4]->Fill(TMath::Min(GenBosons->At(i)->Mass(),199.999));
      if (GenBosons->At(i)->Is(MCParticle::kW))
        hDGenBosons[5]->Fill(TMath::Min(GenBosons->At(i)->Mass(),199.999));
      if (GenBosons->At(i)->Is(MCParticle::kZ))
        hDGenBosons[6]->Fill(TMath::Min(GenBosons->At(i)->Mass(),199.999));
      hDGenBosons[7]->Fill(GenBosons->At(i)->Rapidity());
      if (GenBosons->At(i)->Mass() > 60 && GenBosons->At(i)->Mass() < 120) {
        hDGenBosons[8] ->Fill(GenBosons->At(i)->Pt());
        hDGenBosons[9] ->Fill(GenBosons->At(i)->Eta());
        hDGenBosons[10]->Fill(GenBosons->At(i)->Rapidity());        
      }
    }
    hDGenBosons[11]->Fill(TMath::Min((double)(sumV[0] + 4*sumV[1]),12.4999));

    // photons
    hDGenPhotons[0]->Fill(GenPhotons->GetEntries());
    for (UInt_t i=0; i<GenPhotons->GetEntries(); i++) {
      hDGenPhotons[1]->Fill(GenPhotons->At(i)->Pt());
      hDGenPhotons[2]->Fill(GenPhotons->At(i)->Eta());
    } 

    // Rad photons
    hDGenRadPhotons[0]->Fill(GenRadPhotons->GetEntries());
    for (UInt_t i=0; i<GenRadPhotons->GetEntries(); i++) {
      hDGenRadPhotons[1]->Fill(TMath::Min(GenRadPhotons->At(i)->Pt(),199.999));
      hDGenRadPhotons[2]->Fill(TMath::Min(GenRadPhotons->At(i)->AbsEta(),4.999));
      hDGenRadPhotons[3]->Fill(TMath::Min((double)GenRadPhotons->At(i)->DistinctMother()->Status(),
                                          19.499));
      hDGenRadPhotons[4]->Fill(GenRadPhotons->At(i)->IsGenerated() + 
                               2*GenRadPhotons->At(i)->IsSimulated());
      if (GenRadPhotons->At(i)->DistinctMother()) {
        hDGenRadPhotons[5]->Fill(TMath::Min(
                                 MathUtils::DeltaR(*GenRadPhotons->At(i),
                                                   *GenRadPhotons->At(i)->DistinctMother()),
                                 4.999));
        CompositeParticle *object = new CompositeParticle();
        object->AddDaughter(GenRadPhotons->At(i));
        object->AddDaughter(GenRadPhotons->At(i)->DistinctMother());
        hDGenRadPhotons[6]->Fill(TMath::Min(object->Mass(),99.999));
        delete object;
      }
      Int_t Mother = 0;
      if (GenRadPhotons->At(i)->DistinctMother() &&
         GenRadPhotons->At(i)->DistinctMother()->Is(MCParticle::kMu)) Mother = 1;
      hDGenRadPhotons[7]->Fill(Mother);
    }

    // ISR photons
    hDGenISRPhotons[0]->Fill(GenISRPhotons->GetEntries());
    for (UInt_t i=0; i<GenISRPhotons->GetEntries(); i++) {
      hDGenISRPhotons[1]->Fill(TMath::Min(GenISRPhotons->At(i)->Pt(),199.999));
      hDGenISRPhotons[2]->Fill(TMath::Min(GenISRPhotons->At(i)->AbsEta(),4.999));
      hDGenISRPhotons[3]->Fill(TMath::Min((Double_t)GenISRPhotons->At(i)->DistinctMother()->Status(),
                                          19.499));
      hDGenISRPhotons[4]->Fill(GenISRPhotons->At(i)->IsGenerated() + 
                               2*GenISRPhotons->At(i)->IsSimulated());
      if (GenISRPhotons->At(i)->DistinctMother()) {
      	hDGenISRPhotons[5]->Fill(TMath::Min(
      				 MathUtils::DeltaR(*GenISRPhotons->At(i),
      						   *GenISRPhotons->At(i)->DistinctMother()),4.999));
      	CompositeParticle *object = new CompositeParticle();
      	object->AddDaughter(GenISRPhotons->At(i));
      	object->AddDaughter(GenISRPhotons->At(i)->DistinctMother());
      	hDGenISRPhotons[6]->Fill(TMath::Min(object->Mass(),99.999));
      	delete object;
      }
    }
  } // Fill histograms

  // Apply ISR+Rad filter (but filling all histograms)
  if (fApplyISRFilter == kTRUE &&
     (GenISRPhotons->GetEntries() > 0 || GenRadPhotons->GetEntries() > 0)) {
    SkipEvent();
  }

}

//--------------------------------------------------------------------------------------------------
void GeneratorMod::SlaveBegin()
{
  // Book branch and histograms if wanted.

  // Publish Arrays For the Output Module
  PublishObj(fGenLeptons);  
  PublishObj(fGenAllLeptons);
  PublishObj(fGenTaus);
  PublishObj(fGenNeutrinos);
  PublishObj(fGenQuarks);
  PublishObj(fGenqqHs);
  PublishObj(fGenBosons);
  PublishObj(fGenPhotons);
  PublishObj(fGenRadPhotons);
  PublishObj(fGenISRPhotons);

  // fill histograms
  if (GetFillHist()) {
    // MET
    AddTH1(hDGenMet[0],"hDGenMet_0","Gen MET Pt;p_{t} [GeV];#",200,0,200);
    AddTH1(hDGenMet[1],"hDGenMet_1","Gen MET Px;p_{x} [GeV];#",400,-200,200); 
    AddTH1(hDGenMet[2],"hDGenMet_2","Gen MET Py;p_{y} [GeV];#",400,-200,200); 
    AddTH1(hDGenMet[3],"hDGenMet_3","Gen MET Pz;p_{z} [GeV];#",400,-1000,1000);

    // pt min for leptons from W/Z
    AddTH1(hDGenPtMin, "hDGenPtMin","Pt min leptons from W/Z;p_{t} [GeV];#",200,0.0,200.0); 

    // leptons from W/Z
    AddTH1(hDGenLeptons[0], "hDGenLeptons_0",
           "Number of leptons from W/Z;N_{leptons};#",10,-0.5,9.5); 
    AddTH1(hDGenLeptons[1], "hDGenLeptons_1","Pt leptons from W/Z;p_{t} [GeV];#",100,0.0,200.0); 
    AddTH1(hDGenLeptons[2], "hDGenLeptons_2","Eta leptons from W/Z;#eta;#",50,0.0,5.0); 
    AddTH1(hDGenLeptons[3], "hDGenLeptons_3","Phi leptons from W/Z;#phi;#",90,0.0,180.0); 
    AddTH1(hDGenLeptons[4], "hDGenLeptons_4","Dilepton mass from W/Z;m_{ll};#",1000,0.0,1000.0);
    AddTH1(hDGenLeptons[5], "hDGenLeptons_5","Eta Max for 2 lepton case;#eta;#",50,0.0,5.0); 
    AddTH1(hDGenLeptons[6], "hDGenLeptons_6","Eta Min for 2 lepton case;#eta;#",50,0.0,5.0); 
    AddTH1(hDGenLeptons[7], "hDGenLeptons_7",
           "Pt Max for 2 lepton case;p_{t} [GeV];#",100,0.0,200.0); 
    AddTH1(hDGenLeptons[8], "hDGenLeptons_8",
           "Pt Min for 2 lepton case;p_{t} [GeV];#",100,0.0,200.0); 
    AddTH1(hDGenLeptons[9], "hDGenLeptons_9",
           "DiLepton mass for 2 lepton case;p_{t} [GeV];#",1000,0.0,1000.0);
    AddTH1(hDGenLeptons[10],"hDGenLeptons_10",
           "Delta Phi ll for 2 lepton case;#Delta#phi_{ll};#",90,0.0,180.0); 
    AddTH1(hDGenLeptons[11],"hDGenLeptons_11","Delta R ll;#Delta R_{ll};#",100,0.0,5.0); 
    AddTH1(hDGenLeptons[12],"hDGenLeptons_12",
           "Delta Phi ll for 2 lepton case;#Delta#phi_{ll};#",90,0.0,180.0); 
    AddTH1(hDGenLeptons[13],"hDGenLeptons_13","Delta R ll;#Delta R_{ll};#",100,0.0,5.0); 
    AddTH1(hDGenLeptons[14],"hDGenLeptons_14",
           "Pt Max for 3 lepton case;p_{t} [GeV];#",100,0.0,200.0); 
    AddTH1(hDGenLeptons[15],"hDGenLeptons_15",
           "Pt 2nd for 3 lepton case;p_{t} [GeV];#",100,0.0,200.0); 
    AddTH1(hDGenLeptons[16],"hDGenLeptons_16",
           "Pt Min for 3 lepton case;p_{t} [GeV];#",100,0.0,200.0); 
    AddTH1(hDGenLeptons[17],"hDGenLeptons_17",
           "Dilepton mass for 3 lepton case;m_{ll};#",1000,0.0,1000.0); 
    AddTH1(hDGenLeptons[18],"hDGenLeptons_18",
           "Trilepton mass for 3 lepton case;m_{lll};#",1000,0.0,1000.0); 
    AddTH1(hDGenLeptons[19],"hDGenLeptons_19",
           "Delta R Minimum between leptons for 3 lepton case;#Delta R_{ll};#",100,0.0,5.0); 
    AddTH1(hDGenLeptons[20],"hDGenLeptons_20",
           "Pt Max for 4 lepton case;p_{t} [GeV];#",100,0.0,200.0); 
    AddTH1(hDGenLeptons[21],"hDGenLeptons_21",
           "Pt 2nd for 4 lepton case;p_{t} [GeV];#",100,0.0,200.0); 
    AddTH1(hDGenLeptons[22],"hDGenLeptons_22",
           "Pt 3rd for 4 lepton case;p_{t} [GeV];#",100,0.0,200.0); 
    AddTH1(hDGenLeptons[23],"hDGenLeptons_23",
           "Pt 4th for 4 lepton case;#",100,0.0,200.0); 
    AddTH1(hDGenLeptons[24],"hDGenLeptons_24",
           "Dilepton mass for 4 lepton case;m_{ll};#",1000,0.0,1000.0); 
    AddTH1(hDGenLeptons[25],"hDGenLeptons_25",
           "Fourlepton mass for 3 lepton case;m_{llll};#",1000,0.0,1000.0); 
    AddTH1(hDGenLeptons[26],"hDGenLeptons_26",
           "Delta R Minimum between leptons for 4 lepton case;#Delta R_{ll};#",100,0.0,5.0); 
           
           //leptons for monojet generator-level info
    AddTH1(hDGenWMuons,"hDGenWMuons","Number of Muons in Generated W Events;N_{muons};#",10,-0.5,9.5);
    AddTH1(hDGenWElectrons,"hDGenWElectrons","Number of Electrons in Generated W Events;N_{electrons};#",10,-0.5,9.5);
    
    // all leptons
    AddTH1(hDGenAllLeptons[0], "hDGenAllLeptons_0",
           "Number of all leptons;N_{leptons};#",10,-0.5,9.5); 
    AddTH1(hDGenAllLeptons[1], "hDGenAllLeptons_1","Pt all leptons;p_{t} [GeV];#",400,0.0,200.0); 
    AddTH1(hDGenAllLeptons[2], "hDGenAllLeptons_2","Eta all leptons;#eta;#",50,0.0,5.0); 
    AddTH1(hDGenAllLeptons[3], "hDGenAllLeptons_3","Phi all leptons;#phi;#",90,0.0,180.0); 
    AddTH1(hDGenAllLeptons[4], "hDGenAllLeptons_4","Pt second lepton;p_{t} [GeV];#",400,0.0,200.0); 
    AddTH1(hDGenAllLeptons[5], "hDGenAllLeptons_5",
           "Number of all leptons including taus;N_{leptons};#",10,-0.5,9.5); 
    AddTH1(hDGenAllLeptons[6], "hDGenAllLeptons_6","MinMass; [GeV];#",100,0.0,100.0); 
    AddTH1(hDGenAllLeptons[7], "hDGenAllLeptons_7","MaxMass; [GeV];#",100,0.0,100.0); 

    // taus
    AddTH1(hDGenTaus[0], "hDGenTaus_0","Number of taus;N_{tau};#",10,-0.5,9.5); 
    AddTH1(hDGenTaus[1], "hDGenTaus_1","Pt taus;p_{t} [GeV];#",100,0.0,200.0); 
    AddTH1(hDGenTaus[2], "hDGenTaus_2","Eta taus;#eta;#",50,0.0,5.0); 
    AddTH1(hDGenTaus[3], "hDGenTaus_3","Phi taus;#phi;#",90,0.0,180.0); 

    // neutrinos
    AddTH1(hDGenNeutrinos[0], "hDGenNeutrinos_0","Number of neutrinos;N_{#nu};#",10,-0.5,9.5); 
    AddTH1(hDGenNeutrinos[1], "hDGenNeutrinos_1","Pt neutrinos;p_{t} [GeV];#",100,0.0,200.0);
    AddTH1(hDGenNeutrinos[2], "hDGenNeutrinos_2","Eta neutrinos;#eta;#",50,0.0,5.0); 
    AddTH1(hDGenNeutrinos[3], "hDGenNeutrinos_3","Phi neutrinos;#phi;#",90,0.0,180.0); 

    // quarks
    AddTH1(hDGenQuarks[0], "hDGenQuarks_0", "Number of quarks;N_{quarks};#",10,-0.5,9.5); 
    AddTH1(hDGenQuarks[1], "hDGenQuarks_1", "dijet pt for quarks;p_{t} [GeV];#",200,0.0,400.); 
    AddTH1(hDGenQuarks[2], "hDGenQuarks_2", "dijet mass for quarks;m_{jj};#",2000,0.0,2000.);
    AddTH1(hDGenQuarks[3], "hDGenQuarks_3", 
           "dijet pt for quarks with |eta|<2.5;p_{t} [GeV];#",200,0.0,400.); 
    AddTH1(hDGenQuarks[4], "hDGenQuarks_4", 
           "dijet mass for quarks with |eta|<2.5;m_{jj};#",2000,0.0,2000.);
    AddTH1(hDGenQuarks[5], "hDGenQuarks_5", "Pt for b quarks;p_{t} [GeV];#",200,0.0,400.); 
    AddTH1(hDGenQuarks[6], "hDGenQuarks_6", "Eta for b quarks;#eta;#",200,-10.0,10.); 
    AddTH1(hDGenQuarks[7], "hDGenQuarks_7", 
           "Phi for b quarks;#phi;#",200,-TMath::Pi(),TMath::Pi());
    AddTH1(hDGenQuarks[8], "hDGenQuarks_8", 
           "Pt for b quarks with |eta|<2.5;p_{t} [GeV];#",200,0.0,400.); 
    AddTH1(hDGenQuarks[9], "hDGenQuarks_9", 
           "Eta for b quarks with |eta|<2.5;#eta;#",200,-10.0,10.); 
    AddTH1(hDGenQuarks[10],"hDGenQuarks_10",
           "Phi for b quarks with |eta|<2.5;#phi;#",200,-TMath::Pi(),TMath::Pi());
    AddTH1(hDGenQuarks[11],"hDGenQuarks_11","Pt for t quarks;p_{t} [GeV];#",200,0.0,400.); 
    AddTH1(hDGenQuarks[12],"hDGenQuarks_12","Eta for t quarks;#eta;#",200,-10.0,10.); 
    AddTH1(hDGenQuarks[13],"hDGenQuarks_13",
           "Phi for t quarks;#phi;#",200,-TMath::Pi(),TMath::Pi());
    AddTH1(hDGenQuarks[14],"hDGenQuarks_14","Pt for light quarks;p_{t} [GeV];#",200,0.0,400.); 
    AddTH1(hDGenQuarks[15],"hDGenQuarks_15","Eta for light quarks;#eta;#",200,-10.0,10.); 
    AddTH1(hDGenQuarks[16],"hDGenQuarks_16",
           "Phi for light quarks;#phi;#",200,-TMath::Pi(),TMath::Pi());

    // qqH
    AddTH1(hDGenWBF[0], "hDGenWBF_0",
           "Delta Phi jj for WBF quarks;#Delta Phi_{jj};#",90,0.0,180.);
    AddTH1(hDGenWBF[1], "hDGenWBF_1",
           "Delta Eta jj for WBF quarks;#Delta #eta_{jj};#",100,0.0,10.);  
    AddTH1(hDGenWBF[2], "hDGenWBF_2",
           "Pt max for WBF quarks;p_{t} [GeV];#",200,0.0,400.); 
    AddTH1(hDGenWBF[3], "hDGenWBF_3",
           "Pt min for WBF quarks;p_{t} [GeV];#",200,0.0,400.);
    AddTH1(hDGenWBF[4], "hDGenWBF_4",
           "dijet mass for WBF quarks;m_{jj};#",200,0.0,4000.);

    // bosons
    AddTH1(hDGenBosons[0], "hDGenBosons_0", "Number of bosons;N_{bosons};#",10,-0.5,9.5); 
    AddTH1(hDGenBosons[1], "hDGenBosons_1", "Pt of bosons;p_{t} [GeV];#",200,0.0,400.0); 
    AddTH1(hDGenBosons[2], "hDGenBosons_2", "Eta of bosons;#eta;#",100,-10.0,10.0); 
    AddTH1(hDGenBosons[3], "hDGenBosons_3", "Mass of bosons;Mass;#",2000,0.0,2000.0);
    AddTH1(hDGenBosons[4], "hDGenBosons_4", "Mass of bosons;m_{V};#",200,0.0,200.0);
    AddTH1(hDGenBosons[5], "hDGenBosons_5", "Mass of W bosons;m_{W};#",200,0.0,200.0);
    AddTH1(hDGenBosons[6], "hDGenBosons_6", "Mass of Z bosons;m_{Z};#",200,0.0,200.0);
    AddTH1(hDGenBosons[7], "hDGenBosons_7", "Rapidity of bosons;rapidity;#",100,-10.0,10.0); 
    AddTH1(hDGenBosons[8], "hDGenBosons_8", "Pt of bosons;p_{t} [GeV];#",200,0.0,400.0); 
    AddTH1(hDGenBosons[9], "hDGenBosons_9", "Eta of bosons;#eta;#",100,-10.0,10.0); 
    AddTH1(hDGenBosons[10],"hDGenBosons_10","Rapidity of bosons;rapidity;#",100,-10.0,10.0); 
    AddTH1(hDGenBosons[11],"hDGenBosons_11","Number of W bosons + 4 * Z bosons;Number;#",13,-0.5,12.5); 

    // photons
    AddTH1(hDGenPhotons[0], "hDGenPhotons_0", "Number of photons;N_{photons};#",10,-0.5,9.5); 
    AddTH1(hDGenPhotons[1], "hDGenPhotons_1", "Pt of photons;p_{t} [GeV];#",200,0.0,400.0); 
    AddTH1(hDGenPhotons[2], "hDGenPhotons_2", "Eta of photons;#eta;#",100,-5.0,5.0); 

    //  rad photons
    AddTH1(hDGenRadPhotons[0], "hDGenRadPhotons_0", 
           "Number of radiative photons;N_{photons};#",10,-0.5,9.5); 
    AddTH1(hDGenRadPhotons[1], "hDGenRadPhotons_1", 
           "Pt of radiative photons;p_{t} [GeV];#",400,0.0,200.0);
    AddTH1(hDGenRadPhotons[2], "hDGenRadPhotons_2", 
           "Eta of radiative photons;#eta;#",100,0.0,5.0); 
    AddTH1(hDGenRadPhotons[3], "hDGenRadPhotons_3", 
           "Status of mother of radiative photons;Status;#",20,-0.5,19.5); 
    AddTH1(hDGenRadPhotons[4], "hDGenRadPhotons_4", 
           "IsGenerated+2*IsSimulated of radiative photons;IsGenerated+2*IsSimulated;#",
           4,-0.5,3.5); 
    AddTH1(hDGenRadPhotons[5], "hDGenRadPhotons_5", 
           "Delta R between photon and mother of radiative photons;#Delta R;#",500,0.0,5.0); 
    AddTH1(hDGenRadPhotons[6], "hDGenRadPhotons_6", 
           "Mass photon-photon mother;Mass;#",500,0.0,100.0); 
    AddTH1(hDGenRadPhotons[7], "hDGenRadPhotons_7", 
           "Number of radiative photon with muon as a mother;Status;#",2,-0.5,1.5); 

    //  ISR photons
    AddTH1(hDGenISRPhotons[0], "hDGenISRPhotons_0", 
           "Number of ISR photons;N_{photons};#",10,-0.5,9.5); 
    AddTH1(hDGenISRPhotons[1], "hDGenISRPhotons_1", 
           "Pt of ISR photons;p_{t} [GeV];#",400,0.0,200.0);
    AddTH1(hDGenISRPhotons[2], "hDGenISRPhotons_2", 
           "Eta of ISR photons;#eta;#",100,0.0,5.0); 
    AddTH1(hDGenISRPhotons[3], "hDGenISRPhotons_3", 
           "Status of mother of radiative photons;#eta;#",20,-0.5,19.5);
    AddTH1(hDGenISRPhotons[4], "hDGenISRPhotons_4", 
           "IsGenerated+2*IsSimulated of radiative photons;IsGenerated+2*IsSimulated;#",
           4,-0.5,3.5); 
    AddTH1(hDGenISRPhotons[5], "hDGenISRPhotons_5", 
           "Delta R between photon and mother of ISR photons;#Delta R;#",500,0.0,5.0); 
    AddTH1(hDGenISRPhotons[6], "hDGenISRPhotons_6", 
           "Mass photon-photon mother;Mass;#",500,0.0,100.0); 

    // auxiliar for MG studies
    AddTH1(hDVMass[0], "hDVMass_0", "Mass of munu candidates  ;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVMass[1], "hDVMass_1", "Mass of elnu candidates  ;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVMass[2], "hDVMass_2", "Mass of taunu candidates ;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVMass[3], "hDVMass_3", "Mass of mumu candidates  ;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVMass[4], "hDVMass_4", "Mass of ee candidates    ;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVMass[5], "hDVMass_5", "Mass of tautau candidates;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVMass[6], "hDVMass_6", "Mass of numunumu candidates;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVMass[7], "hDVMass_7", "Mass of nuenue candidates;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVMass[8], "hDVMass_8", "Mass of nutaunutau candidates;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVMass[9], "hDVMass_9", 
           "Mass of munu candidates for t events ;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVMass[10],"hDVMass_10",
           "Mass of elnu candidates for t events ;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVMass[11],"hDVMass_11",
           "Mass of taunu candidates for t events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVMass[12],"hDVMass_12",
           "Mass of qq candidates for t events;Mass [GeV];#",200,0.,200.);

    // Special study about VVjets
    AddTH1(hDVVMass[0], "hDVVMass_0", "Mass of munu for WW events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[1], "hDVVMass_1", "Mass of munu WZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[2], "hDVVMass_2", "Mass of elnu WW events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[3], "hDVVMass_3", "Mass of elnu WZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[4], "hDVVMass_4", "Mass of taunu WW events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[5], "hDVVMass_5", "Mass of taunu WZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[6], "hDVVMass_6", "Mass of mumu WZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[7], "hDVVMass_7", "Mass of mumu ZZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[8], "hDVVMass_8", "Mass of ee WZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[9], "hDVVMass_9", "Mass of ee ZZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[10],"hDVVMass_10","Mass of tautau WZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[11],"hDVVMass_11","Mass of tautau ZZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[12],"hDVVMass_12","Mass of numunumu WZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[13],"hDVVMass_13","Mass of numunumu ZZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[14],"hDVVMass_14","Mass of nuenue WZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[15],"hDVVMass_15","Mass of nuenue ZZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[16],"hDVVMass_16","Mass of nutaunutau WZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[17],"hDVVMass_17","Mass of nutaunutau ZZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[18],"hDVVMass_18","Ratios for WW events;Type;#",7,-0.5,6.5); 
    AddTH1(hDVVMass[19],"hDVVMass_19","Ratios for WZ events;Type;#",10,-0.5,9.5); 
    AddTH1(hDVVMass[20],"hDVVMass_20","Ratios for ZZ2l events;Type;#",7,-0.5,6.5); 
    AddTH1(hDVVMass[21],"hDVVMass_21","Ratios for ZZ4l events;Type;#",16,-0.5,15.5);
    AddTH1(hDVVMass[22],"hDVVMass_22","Maximum mass for WW events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[23],"hDVVMass_23","Minimum mass for WW events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[24],"hDVVMass_24","Maximum mass for WZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[25],"hDVVMass_25","Minimum mass for WZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[26],"hDVVMass_26","Maximum mass for ZZ events;Mass [GeV];#",200,0.,200.);
    AddTH1(hDVVMass[27],"hDVVMass_27","Minimum mass for ZZ events;Mass [GeV];#",200,0.,200.);

  }
}
