#ifndef RHO_UTILITIES_H
#define RHO_UTILITIES_H

namespace mithep{

  struct RhoUtilities {
    enum RhoType {
      DEFAULT = 0,
      MIT_RHO_VORONOI_LOW_ETA,       // was RhoLowEta();
      MIT_RHO_VORONOI_HIGH_ETA,      // was Rho();
      MIT_RHO_RANDOM_LOW_ETA,        // was RhoRandomLowEta();
      MIT_RHO_RANDOM_HIGH_ETA,       // was RhoRandom();
      CMS_RHO_RHOKT6PFJETS,           // was RhoKt6PFJets();
      CMS_RHO_FIXEDGRIDFASTJETALL
    };
  };
}

#endif
