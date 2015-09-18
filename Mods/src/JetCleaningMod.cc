#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/Photon.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"

using namespace mithep;

ClassImp(mithep::JetCleaningMod)

//--------------------------------------------------------------------------------------------------
JetCleaningMod::JetCleaningMod(const char *name, const char *title) :
  BaseMod(name, title),
  fCleanElectronsName(ModNames::gkCleanElectronsName),
  fCleanMuonsName(ModNames::gkCleanMuonsName),
  fCleanPhotonsName(ModNames::gkCleanPhotonsName),
  fCleanTausName(ModNames::gkCleanTausName),
  fGoodJetsName(ModNames::gkGoodJetsName),
  fCleanJets(new JetOArr),
  fMinDeltaR{0.3, 0.3, 0.3, 0.3}
{
  fCleanJets->SetName(ModNames::gkCleanJetsName);
}

JetCleaningMod::~JetCleaningMod()
{
  delete fCleanJets;
}

//--------------------------------------------------------------------------------------------------
void
JetCleaningMod::Process()
{
  // Process entries of the tree.

  fCleanJets->Reset();

  // get input collections
  const JetCol* GoodJets = GetObject<JetCol>(fGoodJetsName);

  ParticleCol const* cleanCollections[4]{};
  if (!fCleanElectronsName.IsNull())
    cleanCollections[0] = GetObject<ParticleCol>(fCleanElectronsName);
  if (!fCleanMuonsName.IsNull())
    cleanCollections[1] = GetObject<ParticleCol>(fCleanMuonsName);
  if (!fCleanTausName.IsNull())
    cleanCollections[2] = GetObject<ParticleCol>(fCleanTausName);
  if (!fCleanPhotonsName.IsNull())
    cleanCollections[3] = GetObject<ParticleCol>(fCleanPhotonsName);

  // create output collection
  // remove any jet that overlaps in eta, phi with an isolated electron.
  for (UInt_t i=0; i<GoodJets->GetEntries(); ++i) {
    const Jet *jet = GoodJets->At(i);

    // check for overlap with clean objects
    unsigned iC = 0;
    for (; iC != 4; ++iC) {
      auto* collection = cleanCollections[iC];
      if (!collection)
        continue;

      unsigned iO = 0;
      for (; iO != collection->GetEntries(); ++iO) {
        auto* part = collection->At(iO);
        if (part->ObjType() == kPhoton) {
          if (MathUtils::DeltaR(static_cast<Photon const*>(part)->SCluster(), jet) < fMinDeltaR[iC])
            break;
        }
        else {
          if (MathUtils::DeltaR(part, jet) < fMinDeltaR[iC])
            break;
        }
      }

      if (iO != collection->GetEntries()) {
        // had a match
        break;
      }
    }

    if (iC != 4) {
      // one of the collections had a match
      continue;
    }

    fCleanJets->Add(jet);
  }

  // sort according to pt
  fCleanJets->Sort();
}

void
JetCleaningMod::SlaveBegin()
{
  PublishObj(fCleanJets);
}

void
JetCleaningMod::SlaveTerminate()
{
  RetractObj(fCleanJets->GetName());
}
