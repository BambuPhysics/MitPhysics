#include "MitPhysics/Utils/interface/QGTagger.h"

#include "TSystem.h"

ClassImp(mithep::QGTagger)

using namespace mithep;
 
QGTagger::QGTagger(bool useCHS)
{
  TString dataDir(gSystem->Getenv("MIT_DATA"));
  if (dataDir.Length() == 0)
    throw std::runtime_error("MIT_DATA environment is not set.");

  // Constructor
  qgLikelihood = new QGLikelihoodCalculator(dataDir, useCHS);
}

QGTagger::~QGTagger() 
{
  // Destructor.
}

Float_t QGTagger::QGValue() {
 
  // check if all the needed quantities are there
 const char * variable_names[] = {"pt", "eta", "rhoIso", "mult", "ptD", "axis2"};

  for (uint i=0; i<sizeof(variable_names)/sizeof(variable_names[0]); i++) {
    if (!variables.count(variable_names[i])) {
      std::cout<<"[QGTagger] WARNING: "<<variable_names[i]<<" has not been filled, Q-G discrimanation not computed"<<std::endl;
      return -1;
    }
  }

  return qgLikelihood->QGvalue(variables);
}

void QGTagger::CalculateVariables(const PFJet *jet, const VertexCol *vertices){

  // quantities we need to fill: pt, eta, rhoIso, mult, ptD, axis2
  variables["pt"] = jet->Pt(); // no correction applied, assume we get already corrected/calibrated jets
  variables["eta"] = jet->Eta();  
  
  // counters
  UInt_t chargedMultiplicity = 0, neutralMultiplicity = 0;
  // global quantities
  Float_t sumWeight = 0., sumDeltaEta = 0., sumDeltaPhi = 0., sumDeltaEta2 = 0., sumDeltaPhi2 = 0., sumDeltaEtaDeltaPhi = 0., sumPt = 0.;


  // loop over constituents
  for(UInt_t i=0;i<jet->NPFCands();i++){

    // charged candidates
    const Track  * track = jet->PFCand(i)->TrackerTrk();
    bool isGoodForAxis = false;
    if (track) { 

      // get closest vertex
      Float_t deltaZ = 999;
      Int_t closestVertexIndex = -1;
      for (UInt_t j=0;j<vertices->GetEntries();j++){
        if (TMath::Abs(track->DzCorrected(*vertices->At(j))) < deltaZ) {
	  deltaZ = fabs(track->DzCorrected(*vertices->At(j)));
	  closestVertexIndex = j;
	}
      }
      
      // assume the first vertex is the good one 
      if (closestVertexIndex == 0 and track->Quality().Quality(TrackQuality::highPurity)) {
	deltaZ = fabs(track->DzCorrected(*vertices->At(0)));
	Float_t deltaZErr = sqrt(pow(track->DszErr(),2) + pow(vertices->At(0)->ZErr(),2));
	if (deltaZ/deltaZErr < 5) {
	  isGoodForAxis = true;
	  Float_t deltaXY = track->D0Corrected(*vertices->At(0));
	  Float_t deltaXYErr = sqrt(pow(track->DxyErr(),2) + pow(vertices->At(0)->XErr(),2) + pow(vertices->At(0)->YErr(),2));
	  if ( fabs(deltaXY)/deltaXYErr < 5) chargedMultiplicity++;
	}
      }
    }

    // neutral candidates
    else {
      if (jet->PFCand(i)->Pt() > 1.0) neutralMultiplicity++;
      isGoodForAxis = true;
    }

    // update global quantities
    Float_t deltaEta = jet->PFCand(i)->Eta() - jet->Eta();
    Float_t deltaPhi = 2*atan( tan( (jet->PFCand(i)->Phi() - jet->Phi())/2 ) );
    Float_t weight = pow(jet->PFCand(i)->Pt(),2);
    if (isGoodForAxis){
      sumWeight += weight;
      sumPt += jet->PFCand(i)->Pt();
      sumDeltaEta += deltaEta*weight;
      sumDeltaPhi += deltaPhi*weight;
      sumDeltaEta2 += deltaEta*deltaEta*weight;
      sumDeltaEtaDeltaPhi += deltaEta*deltaPhi*weight;
      sumDeltaPhi2 += deltaPhi*deltaPhi*weight;
    }
  }    

  // calculate axis and ptD
  Float_t a = 0., b = 0., c = 0.;
  Float_t averageDeltaEta = 0., averageDeltaPhi = 0., averageDeltaEta2 = 0., averageDeltaPhi2 = 0.;
  if (sumWeight > 0) {
    variables["ptD"] = sqrt(sumWeight)/sumPt;
    averageDeltaEta = sumDeltaEta/sumWeight;
    averageDeltaPhi = sumDeltaPhi/sumWeight;
    averageDeltaEta2 = sumDeltaEta2/sumWeight;
    averageDeltaPhi2 = sumDeltaPhi2/sumWeight;
    a = averageDeltaEta2 - averageDeltaEta*averageDeltaEta;                          
    b = averageDeltaPhi2 - averageDeltaPhi*averageDeltaPhi;                          
    c = -(sumDeltaEtaDeltaPhi/sumWeight - averageDeltaEta*averageDeltaPhi);                
  }
  else
    variables["ptD"] = 0;

  Float_t delta = sqrt(fabs((a-b)*(a-b)+4*c*c));
  if (a+b+delta > 0)
    variables["axis1"] = sqrt(0.5*(a+b+delta));
  else
    variables["axis1"] = 0.;
  if (a+b-delta > 0)
    variables["axis2"] = sqrt(0.5*(a+b-delta));
  else
    variables["axis2"] = 0.;

  variables["mult"] = (chargedMultiplicity + neutralMultiplicity);
}

Float_t QGTagger::GetPtD()
{
  return variables["ptD"];
}

Float_t QGTagger::GetAxis1()
{
  return variables["axis1"];
}

Float_t QGTagger::GetAxis2() {
  return variables["axis2"];
}

Float_t QGTagger::GetMult() {
  return variables["mult"];
}
