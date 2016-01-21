#include "MitPhysics/SelMods/interface/BadEventsFilterMod.h"

#include "MitAna/DataTree/interface/EvtSelData.h"

#include "TFile.h"

#include <algorithm>

ClassImp(mithep::BadEventsFilterMod)

void
mithep::BadEventsFilterMod::SlaveBegin()
{
  fBitMask = 0;

  auto* file = GetCurrentFile();
  if (!file)
    return;

  auto* nameTree = dynamic_cast<TTree*>(file->Get(fLabelTreeName));
  if (!nameTree) {
    SendError(kWarning, "SlaveBegin", "EvtSelData label names are not stored in the file.");
    return;
  }

  std::vector<std::string>* filterLabels(new std::vector<std::string>);
  TBranch* labelBranch(0);
  nameTree->SetBranchAddress(fLabelBranchName, &filterLabels, &labelBranch);
  if (!labelBranch) {
    SendError(kWarning, "SlaveBegin", "EvtSelData label names are not stored in the file.");
    return;
  }

  labelBranch->GetEntry(0);

  for (auto& filt : fEnabledFilters) {
    auto itr(std::find(filterLabels->begin(), filterLabels->end(), filt));
    if (itr == filterLabels->end()) {
      SendError(kWarning, "SlaveBegin", ("MET filter with label " + filt + " is not defined.").c_str());
      continue;
    }

    fBitMask |= (1 << (itr - filterLabels->begin()));
  }

  delete nameTree;
  delete filterLabels;

  if (GetFillHist()) {
    AddTH1(hCounter, "hMETFilterCounter", "Number of events flagged bad", fEnabledFilters.size(), 0., double(fEnabledFilters.size()));
    int iX = 1;
    for (unsigned iB = 0; iB != filterLabels->size(); ++iB) {
      if (((fBitMask >> iB) & 1) != 0)
        hCounter->GetXaxis()->SetBinLabel(iX++, filterLabels->at(iB).c_str());
    }
  }
}

void
mithep::BadEventsFilterMod::Process()
{
  auto* evtSelData = GetObject<mithep::EvtSelData>(fEvtSelDataName);
  if (!evtSelData) {
    SendError(kWarning, "Process", "EvtSelData " + fEvtSelDataName + " not found.");
    return;
  }

  Int_t word = evtSelData->metFiltersWord();

  if (GetFillHist()) {
    int iX = 1;
    for (unsigned iB = 0; iB != 8 * sizeof(Int_t); ++iB) {
      if (((fBitMask >> iB) & 1) == 0)
        continue;
      
      if (((word >> iB) & 1) == 0)
        hCounter->Fill(iX - 0.5);

      ++iX;
    }
  }

  if ((fBitMask & word) != fBitMask)
    SkipEvent();

  return;
}

void
mithep::BadEventsFilterMod::SetFilter(char const* name, Bool_t enable/* = kTRUE*/)
{
  auto itr(std::find(fEnabledFilters.begin(), fEnabledFilters.end(), name));
  if (enable && itr == fEnabledFilters.end())
    fEnabledFilters.emplace_back(name);
  else if (!enable && itr != fEnabledFilters.end())
    fEnabledFilters.erase(itr);
}
