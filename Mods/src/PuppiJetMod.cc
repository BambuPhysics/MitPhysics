#include "MitPhysics/Mods/interface/PuppiJetMod.h"

using namespace mithep;

templateClassImp(PuppiJetMod)

template<> Bool_t PuppiJetMod<FatJet>::PassJet(fastjet::PseudoJet &p) {
  return (p.pt()>200);
}

template<> Bool_t PuppiJetMod<PFJet>::PassJet(fastjet::PseudoJet &p) {
  return (p.pt()>10);
}
