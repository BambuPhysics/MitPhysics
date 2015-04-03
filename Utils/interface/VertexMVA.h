//--------------------------------------------------------------------------------------------------
// $Id $
//
// VertexMVA
//
// Helper Class for photon Identification decisions.
//
// Authors: J.Bendavid
//--------------------------------------------------------------------------------------------------

#ifndef MITPHYSICS_UTILS_VERTEXMVA_H
#define MITPHYSICS_UTILS_VERTEXMVA_H

#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/Electron.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/BaseVertex.h"
#include "MitAna/DataTree/interface/DecayParticleCol.h"
#include "MitAna/DataTree/interface/VertexCol.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/TriggerObjectCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

namespace mithep {
  class VertexMVA {
    public:
      VertexMVA();
      virtual ~VertexMVA() {}
  
      Double_t GetProb(const Vertex *v) const;
      

    ClassDef(VertexMVA, 0) // Vertex tools
  };
}

#endif
