#include "MitPhysics/SelMods/interface/BadEventsFilterMod.h"

#include "MitAna/DataTree/interface/EvtSelData.h"
#include "MitAna/DataTree/interface/EventHeader.h"

#include "TFile.h"
#include "TObjArray.h"

#include <algorithm>

ClassImp(mithep::BadEventsFilterMod)

mithep::BadEventsFilterMod::BadEventsFilterMod(char const* name, char const* title) :
  BaseMod(name, title),
  fFilterNames(),
  fTagResults(0, "BadEventsFilter")
{
  fFilterNames.SetName("BadEventsFilterNames");
  fFilterNames.SetOwner(true);
}

void
mithep::BadEventsFilterMod::AddEventList(char const* name, char const* fileName)
{
  // Input must be an ASCII file with format run:lumi:event, one entry per row.

  auto& list = fEventLists[name];

  std::ifstream input(fileName);
  std::string line;

  while (true) {
    std::getline(input, line);
    if (!input.good())
      break;

    unsigned run(std::atoi(line.substr(0, line.find(":")).c_str()));
    unsigned lumi(std::atoi(line.substr(line.find(":") + 1, line.rfind(":")).c_str()));
    unsigned event(std::atoi(line.substr(line.rfind(":") + 1).c_str()));

    list.emplace(run, lumi, event);
  }
}

void
mithep::BadEventsFilterMod::SlaveBegin()
{
  if (GetFillHist()) {
    AddTH1(hCounter, "hMETFilterCounter", "Number of events flagged bad", fEnabledFilters.size(), 0., double(fEnabledFilters.size()));
  }

  if (fTaggingMode) {
    PublishObj(&fFilterNames);
    PublishObj(&fTagResults);
  }
}

void
mithep::BadEventsFilterMod::BeginRun()
{
  if (fReload == 1) {
    fBitMask = 0;

    fFilterNames.Clear();
    fTagResults.Clear();

    auto* file = GetCurrentFile();
    if (!file)
      return;

    auto* nameTree = dynamic_cast<TTree*>(file->Get(fLabelTreeName));
    if (!nameTree) {
      // is bambu version <= 042
      SendError(kWarning, "BeginRun", "EvtSelData label names are not stored in the file.");

      enum StaticFilter {
        kHBHENoiseFilter,
        kECALDeadCellFilter,
        kTrackingFailureFilter,
        kEEBadScFilter,
        kECALaserCorrFilter,
        kManyStripClusFilter,
        kTooManyStripClusFilter,
        kLogErrorTooManyClustersFilter,
        kCSCTightHaloFilter,
        kCSCLooseHaloFilter,
        nStaticFilters
      };

      TString filterNames[nStaticFilters] = {
        "HBHENoiseFilter",
        "ECALDeadCellFilter",
        "TrackingFailureFilter",
        "EEBadScFilter",
        "ECALaserCorrFilter",
        "ManyStripClusFilter",
        "TooManyStripClusFilter",
        "LogErrorTooManyClustersFilter",
        "CSCTightHaloFilter",
        "CSCLooseHaloFilter"
      };

      for (auto& filt : fEnabledFilters) {
        unsigned iF = std::find(filterNames, filterNames + nStaticFilters, filt) - filterNames;
        if (iF != nStaticFilters)
          fBitMask |= (1 << iF);
        else
          SendError(kAbortAnalysis, "BeginRun", ("MET filter with label " + filt + " is not defined.").c_str());
      }

      if (hCounter) {
        int iX = 1;
        for (unsigned iB = 0; iB != nStaticFilters; ++iB) {
          if (((fBitMask >> iB) & 1) != 0)
            hCounter->GetXaxis()->SetBinLabel(iX++, filterNames[iB]);
        }
        for (auto& nameAndList : fEventLists)
          hCounter->GetXaxis()->SetBinLabel(iX++, nameAndList.first.c_str());
      }

      if (fTaggingMode) {
        for (unsigned iB = 0; iB != nStaticFilters; ++iB) {
          if (((fBitMask >> iB) & 1) != 0) {
            fFilterNames.Add(new TObjString(filterNames[iB]));
            fTagResults.Add(kFALSE);
          }
        }

        for (auto& nameAndList : fEventLists) {
          fFilterNames.Add(new TObjString(nameAndList.first.c_str()));
          fTagResults.Add(kFALSE);
        }
      }

      // 2 -> never try reload again
      fReload = 2;
      return;
    }

    fReload = 0;

    std::vector<std::string>* filterLabels(new std::vector<std::string>);
    TBranch* labelBranch(0);
    nameTree->SetBranchAddress(fLabelBranchName, &filterLabels, &labelBranch);
    if (!labelBranch)
      SendError(kAbortAnalysis, "BeginRun", "EvtSelData label names are not stored in the file.");

    labelBranch->GetEntry(0);

    for (auto& filt : fEnabledFilters) {
      auto itr(std::find(filterLabels->begin(), filterLabels->end(), filt));
      if (itr == filterLabels->end())
        SendError(kAbortAnalysis, "BeginRun", ("MET filter with label " + filt + " is not defined.").c_str());

      fBitMask |= (1 << (itr - filterLabels->begin()));
    }

    if (hCounter) {
      int iX = 1;
      for (unsigned iB = 0; iB != filterLabels->size(); ++iB) {
        if (((fBitMask >> iB) & 1) != 0)
          hCounter->GetXaxis()->SetBinLabel(iX++, filterLabels->at(iB).c_str());
      }
      for (auto& nameAndList : fEventLists)
        hCounter->GetXaxis()->SetBinLabel(iX++, nameAndList.first.c_str());
    }

    if (fTaggingMode) {
      for (unsigned iB = 0; iB != filterLabels->size(); ++iB) {
        if (((fBitMask >> iB) & 1) != 0) {
          fFilterNames.Add(new TObjString(filterLabels->at(iB).c_str()));
          fTagResults.Add(kFALSE);
        }
      }

      for (auto& nameAndList : fEventLists) {
        fFilterNames.Add(new TObjString(nameAndList.first.c_str()));
        fTagResults.Add(kFALSE);
      }
    }

    delete nameTree;
    delete filterLabels;
  }
}

Bool_t
mithep::BadEventsFilterMod::Notify()
{
  if (fReload == 0)
    fReload = 1;

  return true;
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

  unsigned iF = 0;

  if (fTaggingMode) {
    for (unsigned iB = 0; iB != 8 * sizeof(Int_t); ++iB) {
      if (((fBitMask >> iB) & 1) == 0)
        continue;
      
      // In tagging mode, flag is true (normal decision) if bit is off (bad event)
      bool bitFailed = ((word >> iB) & 1) == 0;
      if (fNormalDecision)
        fTagResults.At(iF) = bitFailed;
      else
        fTagResults.At(iF) = !bitFailed;

      ++iF;
    }
  }

  bool badEvent = (fBitMask & word) != fBitMask;

  if ((fTaggingMode || !badEvent) && fEventLists.size() != 0) {
    EventID id(GetEventHeader()->RunNum(), GetEventHeader()->LumiSec(), GetEventHeader()->EvtNum());

    for (auto& nameAndList : fEventLists) {
      auto& list(nameAndList.second);
      auto itr = list.find(id);
      if (itr != list.end()) { // event found in the list
        badEvent = true;
        fTagResults.At(iF) = fNormalDecision;
      }
      else
        fTagResults.At(iF) = !fNormalDecision;

      ++iF;
    }
  }

  if (!fTaggingMode) {
    if (fNormalDecision && badEvent)
      SkipEvent();
    else if (!fNormalDecision && !badEvent)
      SkipEvent();
  }
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
