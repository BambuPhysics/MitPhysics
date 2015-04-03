#include <TH1D.h>
#include <TH2D.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/PartonFlavorHistoryMod.h"

using namespace mithep;
using std::cout;
using std::endl;
using std::string;

ClassImp(mithep::PartonFlavorHistoryMod)

//--------------------------------------------------------------------------------------------------
PartonFlavorHistoryMod::PartonFlavorHistoryMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fMCPartName(Names::gkMCPartBrn),
  fMCSampleType("NotSet"),
  fApplyPartonFlavorFilter(kFALSE),
  fMCType(kMCTypeUndef),
  fParticles(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void PartonFlavorHistoryMod::Process()
{
  // Process entries of the tree.

  // load MCParticle branch
  LoadBranch(fMCPartName);

  //**************************************
  //Parton Flavor Classification For WJets
  //**************************************
  MCParticleOArr *ClassificationPartons = new MCParticleOArr;
  vector<string> FlavorSources;
  vector<MCParticle*> SisterPartons;
  for (UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);
    string FlavorSource = "NULL";  
    Bool_t FoundProgenitor = kFALSE;

    //look for status 2 quarks with daughters - denoted "the quark"
    if ((p->Status() == 2) && 
        (p->IsQuark() && p->AbsPdgId() != 6) && 
        (p->NDaughters() > 0)) {
      
      //make collection of all parents
      MCParticleOArr *parents = new MCParticleOArr;
      const MCParticle *temp = p;
      while (temp->HasMother()) {
        parents->Add(temp->Mother());
        temp = temp->Mother();
      }
      
      MDB(kAnalysis, 4) {
        cout << "quark: " << i << " " << p->PdgId() << " " 
             << p->Status() << " " 
             << p->Pt() << " " 
             << p->Eta() << " " 
             << p->Phi() << " " << endl;
      }
      
      // matrix element:    status 3 parent with precisely 2 grandparents outside of initial
      //                    section (0-5), with same ID as the quark.
      // flavor excitation: only one outgoing parton
      // gluon splitting:   parent is a quark of different flavor than the quark, or parent is a
      //                    gluon. can be ISR or FSR.
      // true decay:        from resonance like top, higgs, etc.  if the quark has no status 3
      //                    partons of same flavor then it is gluon splitting

      // loop over all parents of the quark.
      for (UInt_t j=0; j<parents->GetEntries(); ++j) {
        Int_t parentIndex = -1;
        for (UInt_t k=0; k<fParticles->GetEntries(); ++k) {
          if (fParticles->At(k) == parents->At(j)) {
            parentIndex = k;
            break;
          }
        }
        
        MDB(kAnalysis, 4)
          cout << "parents: " << j << " " << parents->At(j)->PdgId() << " " 
               << parents->At(j)->Status() << " " 
               << parents->At(j)->Pt() << " " 
               << parents->At(j)->Eta() << " " 
               << parents->At(j)->Phi() << " " 
               << parentIndex << endl;
        
        // matrix element          
        if ((parentIndex > 5) && 
            (parents->At(j)->Status() == 3) && 
            (parents->At(j)->PdgId() == p->PdgId())) {            
          if (FlavorSource == "NULL") 
            FlavorSource = "MatrixElement";  
          FoundProgenitor = kTRUE;
        }
	else if ((parents->At(j)->IsGluon()) || 
                   (parents->At(j)->IsQuark() && parents->At(j)->AbsPdgId() != 6)) {
          // gluon splitting
          if (FlavorSource == "NULL") 
            FlavorSource = "GluonSplitting";
          FoundProgenitor = kTRUE;
        }
	else if (parents->At(j)->AbsPdgId() == 2212) {
          if (FlavorSource == "NULL") 
            FlavorSource = "InitialState";
        }
	else if (!parents->At(j)->IsGluon()) {
          // parent is not a gluon or light quark
          if (FlavorSource == "NULL") 
            FlavorSource = "FlavorDecay";
          FoundProgenitor = kTRUE;
        }        
      }

      MDB(kAnalysis, 4)
        cout << "FlavorSource : " << FlavorSource << " " << FoundProgenitor << endl; 
  
      ClassificationPartons->Add(p);
      FlavorSources.push_back(FlavorSource);
      delete parents;
    }
  }

  for (UInt_t i=0; i<ClassificationPartons->GetEntries(); ++i) {
    MCParticle *sister = 0;
    for (UInt_t j=0; j<ClassificationPartons->GetEntries(); ++j) {
      if ((ClassificationPartons->At(i)->PdgId() == -1*ClassificationPartons->At(j)->PdgId()) &&
          (ClassificationPartons->At(i)->Status() == ClassificationPartons->At(j)->Status())) {
        sister = ClassificationPartons->At(j);
        break;
      }
    }
    SisterPartons.push_back(sister);
  }

  //************************************************
  //Count number of B and C quarks, and find sisters
  //************************************************
  Int_t NumberOfBQuarks                       = 0;
  Int_t NumberOfCQuarks                       = 0;
  Int_t NumberOfBQuarksFromMatrixElement      = 0;
  Int_t NumberOfBQuarksFromGluonSplitting     = 0;
  Int_t NumberOfCQuarksFromMatrixElement      = 0;
  Int_t NumberOfCQuarksFromGluonSplitting     = 0;
  Double_t maxDRForBSistersFromMatrixElement  = 0.0;
  Double_t minDRForBSistersFromMatrixElement  = -1.0;
  Double_t maxDRForCSistersFromMatrixElement  = 0.0;
  Double_t minDRForCSistersFromMatrixElement  = -1.0;
  Double_t maxDRForBSistersFromGluonSplitting = 0.0;
  Double_t minDRForBSistersFromGluonSplitting = -1.0;
  Double_t maxDRForCSistersFromGluonSplitting = 0.0;
  Double_t minDRForCSistersFromGluonSplitting = -1.0;

  assert(SisterPartons.size() == FlavorSources.size());
  assert(SisterPartons.size() == ClassificationPartons->GetEntries());

  //we change classification to flavor excitation if the quark had no sister
  for (UInt_t i=0; i<ClassificationPartons->GetEntries(); ++i) {
    if (FlavorSources[i] == "GluonSplitting" && !SisterPartons[i])
      FlavorSources[i] = "FlavorExcitation";
  }

  for (UInt_t i=0; i<ClassificationPartons->GetEntries(); ++i) {
    MDB(kAnalysis, 4) {
      cout << "Classification Partons: " 
           << ClassificationPartons->At(i)->PdgId() << " " 
           << ClassificationPartons->At(i)->Status() << " " 
           << ClassificationPartons->At(i)->Pt() << " " 
           << ClassificationPartons->At(i)->Eta() << " " 
           << ClassificationPartons->At(i)->Phi() << " " 
           << FlavorSources[i] << " ";
    }

    Double_t dr = -1.0;
    if (SisterPartons[i]) {
      MDB(kAnalysis, 4) {
        cout << " : Sister : " 
             << SisterPartons[i]->PdgId() << " "
             << SisterPartons[i]->Status() << " "
             << SisterPartons[i]->Pt() << " "
             << SisterPartons[i]->Eta() << " "
             << SisterPartons[i]->Phi() << " ";
      }
      dr = MathUtils::DeltaR(*ClassificationPartons->At(i),
                             *SisterPartons[i]);
    }

    MDB(kAnalysis, 4) { 
      cout << endl; 
    }

    if (ClassificationPartons->At(i)->AbsPdgId() == 5) {
      NumberOfBQuarks++;
      if (FlavorSources[i] == "GluonSplitting") {
        NumberOfBQuarksFromGluonSplitting++;
        if (dr > -1.0) {
          if (dr > maxDRForBSistersFromGluonSplitting) {
            maxDRForBSistersFromGluonSplitting = dr;
          }
          if ((dr < minDRForBSistersFromGluonSplitting) || 
              (minDRForBSistersFromGluonSplitting < 0.0)) {
            minDRForBSistersFromGluonSplitting = dr;
          }
        }
      } else if ((FlavorSources[i] == "MatrixElement") || 
                 (FlavorSources[i] == "FlavorExcitation")) {
        NumberOfBQuarksFromMatrixElement++;
        if (dr > -1.0) {
          if (dr > maxDRForBSistersFromMatrixElement) {
            maxDRForBSistersFromMatrixElement = dr;
          }
          if ((dr < minDRForBSistersFromMatrixElement) || 
              (minDRForBSistersFromMatrixElement < 0.0)) {
            minDRForBSistersFromMatrixElement = dr;
          }
        }
      }
    }

    if (ClassificationPartons->At(i)->AbsPdgId() == 4) { 
      NumberOfCQuarks++;
      if (FlavorSources[i] == "GluonSplitting") {
        NumberOfCQuarksFromGluonSplitting++;
        if (dr > -1.0) {
          if (dr > maxDRForCSistersFromGluonSplitting) {
            maxDRForCSistersFromGluonSplitting = dr;
          }
          if ((dr < minDRForCSistersFromGluonSplitting) || 
              (minDRForCSistersFromGluonSplitting < 0.0)) {
            minDRForCSistersFromGluonSplitting = dr;
          }
        }
      } else if ((FlavorSources[i] == "MatrixElement") || 
                 (FlavorSources[i] == "FlavorExcitation")) {
        NumberOfCQuarksFromMatrixElement++;
        if (dr > -1.0) {
          if (dr > maxDRForCSistersFromMatrixElement) {
            maxDRForCSistersFromMatrixElement = dr;
          }
          if ((dr < minDRForCSistersFromMatrixElement) || 
              (minDRForCSistersFromMatrixElement < 0.0)) {
            minDRForCSistersFromMatrixElement = dr;
          }
        }
      }      
    }
  }

  //*******************************
  //Determine Flavor Classification
  //*******************************
  Int_t FlavorClassification = -1;

  if (NumberOfBQuarks == 0 && NumberOfCQuarks == 0) {
    FlavorClassification = 11;
  } else if (NumberOfBQuarksFromMatrixElement >= 2 && minDRForBSistersFromMatrixElement > 0.5) {
    FlavorClassification = 1;
  } else if (NumberOfBQuarksFromMatrixElement == 1) {
    FlavorClassification = 2;
  } else if (NumberOfCQuarksFromMatrixElement >= 2 && minDRForCSistersFromMatrixElement > 0.5) {
    FlavorClassification = 3;
  } else if (NumberOfCQuarksFromMatrixElement == 1) {
    FlavorClassification = 4;
  } else if (NumberOfBQuarksFromGluonSplitting >= 2 && maxDRForBSistersFromGluonSplitting < 0.5) {
    FlavorClassification = 5;
  } else if (NumberOfCQuarksFromGluonSplitting >= 2 && maxDRForCSistersFromGluonSplitting < 0.5) {
    FlavorClassification = 6;
  } else {
    FlavorClassification = 0;
  }

  delete ClassificationPartons;
  fFlavorClassification->Fill(FlavorClassification);

  //**********************************
  //Filter Events based on sample type
  //**********************************
  if (fApplyPartonFlavorFilter) {
    //For WJets only accept 5,6,11. Light flavor + bb/cc from gluon splitting at small angle
    if (fMCType == kMCTypeVLightJets ) {
      if (!(FlavorClassification ==  5 || FlavorClassification == 6 
         || FlavorClassification == 11 || FlavorClassification == 0 
         || FlavorClassification ==  4)
        ) {
        SkipEvent();
        return;
      }
    }
    //For VQQ only accept 1,2,3,4. bb/cc from matrix element at large angle
    else if (fMCType == kMCTypeVQQ) {
      if (!(FlavorClassification == 1 || FlavorClassification == 2
         || FlavorClassification == 3 || FlavorClassification == 4)
        ) {
        SkipEvent();
        return;
      }
    }
    
    //For Wc only accept 4. b from matrix element at large angle
    else if (fMCType == kMCTypeWc) {
      if (!(FlavorClassification == 4))
      {      
        SkipEvent();
        return;
      }
    }
    
    else {
      cout << "Warning: The given MCType " << fMCType << " was not recognized.\n";
    }
  }

}

//--------------------------------------------------------------------------------------------------
void PartonFlavorHistoryMod::SlaveBegin()
{
  // Book branch and histograms if wanted.

  ReqBranch(fMCPartName, fParticles);

  AddTH1(fFlavorClassification ,"hFlavorClassification",";FlavorClassification;Number of Events",
         12,-0.5,11.5);
  
  if (string(fMCSampleType.Data()) == "NotSet") {
    cout << "Warning: MCSampleType is not set. Parton Flavor Filter will not be applied.\n";
  }

  if (string(fMCSampleType.Data()) == "kMCTypeVLightJets") {
    fMCType = kMCTypeVLightJets;
  } else if (string(fMCSampleType.Data()) == "kMCTypeWc") {
    fMCType = kMCTypeWc;
  } else if (string(fMCSampleType.Data()) == "kMCTypeVQQ") {
    fMCType = kMCTypeVQQ;
  } else {
    fMCType = kMCTypeUndef;
    fApplyPartonFlavorFilter = false;
  }
}
