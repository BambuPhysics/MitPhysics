#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/SeparatePileUpMod.h"

using namespace mithep;

ClassImp(mithep::SeparatePileUpMod)

//------------------------------------------------------------------------------
SeparatePileUpMod::SeparatePileUpMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fPFCandidatesName(Names::gkPFCandidatesBrn),
  fPFPileUpName("PFPileUp"),
  fPFNoPileUpName("PFNoPileUp"),
  fAllVertexName(Names::gkPVBrn),
  fVertexName(ModNames::gkGoodVertexesName),
  fCheckClosestZVertex(kTRUE),
  fUseAllVertices(kTRUE)
{
  // Constructor.
}

//------------------------------------------------------------------------------
void SeparatePileUpMod::Process()
{
  // Process entries of the tree. 
  auto* pfCandidates = GetObject<PFCandidateCol>(fPFCandidatesName);

  PFCandidateOArr *pfPileUp = new PFCandidateOArr;
  pfPileUp->SetName(fPFPileUpName);
  
  PFCandidateOArr *pfNoPileUp = new PFCandidateOArr;
  pfNoPileUp->SetName(fPFNoPileUpName);

  auto* allVertices = GetObject<VertexCol>(fAllVertexName);
  
  VertexCol const* vertices = 0;
  if (fUseAllVertices == kTRUE)
    vertices = allVertices;
  else
    vertices = GetObject<VertexOArr>(fVertexName);

  for (UInt_t i = 0; i < pfCandidates->GetEntries(); i++) {
    const PFCandidate *pf = pfCandidates->At(i);
    assert(pf);

    if (pf->PFType() == PFCandidate::eHadron) {
      if (pf->HasTrackerTrk() && 
         vertices->At(0)->HasTrack(pf->TrackerTrk()) &&
         vertices->At(0)->TrackWeight(pf->TrackerTrk()) > 0)
      {
        pfNoPileUp->Add(pf);
      }
      else {
        Bool_t vertexFound = kFALSE;
        const Vertex *closestVtx = 0;
        Double_t dzmin = 10000;

	for (UInt_t j = 0; j < allVertices->GetEntries(); j++) {
	  const Vertex *vtx = allVertices->At(j);
	  assert(vtx);

	  if (pf->HasTrackerTrk() && 
	     vtx->HasTrack(pf->TrackerTrk()) &&
	     vtx->TrackWeight(pf->TrackerTrk()) > 0) {
	    vertexFound = kTRUE;
	    closestVtx = vtx;
	    break;
	  }
	  Double_t dz = fabs(pf->SourceVertex().Z() - vtx->Z());
	  if (dz < dzmin) {
	    closestVtx = vtx;
	    dzmin = dz;
	  }
	}

	if (fCheckClosestZVertex) {
	  // Fallback: if track is not associated with any vertex,
	  // associate it with the vertex closest in z
	  if (vertexFound || closestVtx != vertices->At(0))
	    pfPileUp->Add(pf);
	  else
	    pfNoPileUp->Add(pf);
	}
	else {
	  if (vertexFound && closestVtx != vertices->At(0))
	    pfPileUp->Add(pf);
	  else
	    pfNoPileUp->Add(pf); // Ridiculous but that's how it is
	}
      }
    }
    else {
      pfNoPileUp->Add(pf);
    }
  } // Loop over PF candidates

  // add to event for other modules to use
  AddObjThisEvt(pfPileUp);  
  AddObjThisEvt(pfNoPileUp);  
}
