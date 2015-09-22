import sys
import os
from MitAna.TreeMod.bambu import mithep, analysis

from MitPhysics.SelMods.BadEventsFilterMod import badEventsFilterMod
from MitPhysics.Mods.GoodPVFilterMod import goodPVFilterMod
from MitPhysics.Mods.JetIdMod import jetIdMod

mitdata = os.environ['MIT_DATA']

# "typedefs"
kMet = mithep.MonoXSkimMod.kMet
kDielectron = mithep.MonoXSkimMod.kDielectron
kDimuon = mithep.MonoXSkimMod.kDimuon
kSingleElectron = mithep.MonoXSkimMod.kSingleElectron
kSingleMuon = mithep.MonoXSkimMod.kSingleMuon
kPhoton = mithep.MonoXSkimMod.kPhoton

kMonoJet = mithep.MonoXSkimMod.kMonoJet
kMonoPhoton = mithep.MonoXSkimMod.kMonoPhoton

BExpr = mithep.BooleanMod.Expression

monojetTriggerNames = {}
monojetTriggerNames[kMet] = ["HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v*", "HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v*"]
monojetTriggerNames[kSingleMuon] = ["HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v*", "HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v*", "HLT_IsoMu27_v*"]
monojetTriggerNames[kDimuon] = ["HLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v*", "HLT_PFMETNoMu90_NoiseCleaned_PFMHTNoMu90_IDTight_v*", "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*"]
monojetTriggerNames[kSingleElectron] = ["HLT_Ele27_eta2p1_WPLoose_Gsf_v*"]
monojetTriggerNames[kDielectron] = ["HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"]
monojetTriggerNames[kPhoton] = ["HLT_Photon175_v*"]

monophotonTriggerNames = {}
monophotonTriggerNames[kMet] = ['HLT_Photon175_v*']
monophotonTriggerNames[kSingleMuon] = ["HLT_IsoMu27_v*"]
monophotonTriggerNames[kDimuon] = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*"]
monophotonTriggerNames[kSingleElectron] = ["HLT_Ele27_eta2p1_WPLoose_Gsf_v*"]
monophotonTriggerNames[kDielectron] = ["HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*"]
monophotonTriggerNames[kPhoton] = ["HLT_Photon175_v*"]

hltMod = mithep.HLTMod(
    BitsName = "HLTBits",
    TrigObjsName = "MonoJetTriggerObjects",
    AbortIfNotAccepted = analysis.isRealData
)
addedTriggers = []

for l in monojetTriggerNames.values():
    for t in l:
        if t not in addedTriggers:
            hltMod.AddTrigger(t)
            addedTriggers.append(t)

for l in monophotonTriggerNames.values():
    for t in l:
        if t not in addedTriggers:
            hltMod.AddTrigger(t)
            addedTriggers.append(t)

jetCorrection = mithep.JetCorrectionMod(
    InputName = 'AKt4PFJetsCHS',
    CorrectedJetsName = 'CorrectedJets',
    RhoAlgo = mithep.PileupEnergyDensity.kFixedGridFastjetAll
)

metCorrection = mithep.MetCorrectionMod('MetCorrection',
    InputName = 'PFMet',
    OutputName = 'PFType1CorrectedMet',
    JetsName = 'AKt4PFJetsCHS',
    RhoAlgo = mithep.PileupEnergyDensity.kFixedGridFastjetAll,
    MaxEMFraction = 0.9,
    SkipMuons = True
)
metCorrection.ApplyType0(False)
metCorrection.ApplyType1(True)
metCorrection.ApplyShift(False)
metCorrection.IsData(analysis.isRealData)

if analysis.isRealData:
    jecSources = [
        "Summer15_50nsV5_DATA_L1FastJet_AK4PFchs.txt",
        "Summer15_50nsV5_DATA_L2Relative_AK4PFchs.txt",
        "Summer15_50nsV5_DATA_L3Absolute_AK4PFchs.txt",
        "Summer15_50nsV5_DATA_L2L3Residual_AK4PFchs.txt"
    ]
else:
    jecSources = [
        "Summer15_50nsV5_MC_L1FastJet_AK4PFchs.txt",
        "Summer15_50nsV5_MC_L2Relative_AK4PFchs.txt",
        "Summer15_50nsV5_MC_L3Absolute_AK4PFchs.txt"
    ]

for jec in jecSources:
    jetCorrection.AddCorrectionFromFile(mitdata + '/JEC/' + jec)
    metCorrection.AddJetCorrectionFromFile(mitdata + '/JEC/' + jec)

muonId = mithep.MuonIdMod("MuonId",
    MuonClassType = mithep.MuonTools.kPFGlobalorTracker,
    IdType = mithep.MuonTools.kNoId,
    IsoType = mithep.MuonTools.kNoIso,
    ApplyD0Cut = True,
    ApplyDZCut = True,
    WhichVertex = 0,
    PtMin = 10.,
    EtaMax = 2.4,
    OutputName = "VetoMuons"
)

electronId = mithep.ElectronIdMod("ElectronId",
    PtMin = 10.,  
    EtaMax = 2.4,
    ApplyEcalFiducial = True,
    IdType = mithep.ElectronTools.kSummer15Veto,
    IsoType = mithep.ElectronTools.kSummer15VetoIso,
    ConversionsName = "Conversions",
    ApplyConversionFilterType1 = True,
    ApplyConversionFilterType2 = False,
    ApplyD0Cut = True,
    ApplyDZCut = True,
    WhichVertex = 0,
    OutputName = "VetoElectrons"
)

photonId = mithep.PhotonIdMod("VetoPhotonId",
    PtMin = 15.,
    OutputName = "VetoPhotons",
    IdType = mithep.PhotonTools.kSummer15Loose,
    IsoType = mithep.PhotonTools.kSummer15LooseIso,
    ApplyElectronVeto = True
)

jetId = jetIdMod.clone(
    MVATrainingSet = mithep.JetIDMVA.nMVATypes,
    PtMin = 15.
)

monojetPrefilter = mithep.MonoXSkimMod('MonojetPrefilter',
    MetName = metCorrection.GetOutputName(),
    AnalysisType = kMonoJet,
    GoodElectronsName = mithep.Names.gkElectronBrn,
    GoodMuonsName = mithep.Names.gkMuonBrn,
    GoodPhotonsName = mithep.Names.gkPhotonBrn,
    JetsName = jetCorrection.GetOutputName(),
    CategoryFlagsName = 'MonojetPrefilter',
    MinMetPt = 150.,
    MinMonoXPt = 90.
)
for cat, triggers in monojetTriggerNames.items():
    for t in triggers:
        monojetPrefilter.AddTriggerName(cat, t)

monophotonPrefilter = mithep.MonoXSkimMod('MonophotonPrefilter',
    MetName = metCorrection.GetOutputName(),
    AnalysisType = kMonoPhoton,
    GoodElectronsName = mithep.Names.gkElectronBrn,
    GoodMuonsName = mithep.Names.gkMuonBrn,
    GoodPhotonsName = mithep.Names.gkPhotonBrn,
    CategoryFlagsName = 'MonophotonPrefilter',
    MinMetPt = 150.,
    MinMonoXPt = 60.
)
for cat, triggers in monophotonTriggerNames.items():
    for t in triggers:
        monophotonPrefilter.AddTriggerName(cat, t)

prefilter = mithep.BooleanMod('Prefilter',
    Expression = BExpr(monojetPrefilter, monophotonPrefilter, BExpr.kOR)
)

monojetSkim = monojetPrefilter.clone('MonojetSkim',
    GoodElectronsName = electronId.GetOutputName(),
    GoodMuonsName = muonId.GetOutputName(),
    GoodPhotonsName = photonId.GetOutputName(),
    JetsName = jetId.GetOutputName(),
    CategoryFlagsName = 'MonojetSkim'
)
for cat, triggers in monojetTriggerNames.items():
    for t in triggers:
        monojetSkim.AddTriggerName(cat, t)

monophotonSkim = monophotonPrefilter.clone('MonophotonSkim',
    GoodElectronsName = electronId.GetOutputName(),
    GoodMuonsName = muonId.GetOutputName(),
    GoodPhotonsName = photonId.GetOutputName(),
    CategoryFlagsName = 'MonophotonSkim'
)
for cat, triggers in monophotonTriggerNames.items():
    for t in triggers:
        monophotonSkim.AddTriggerName(cat, t)

skim = mithep.BooleanMod('Skim',
    Expression = BExpr(monojetSkim, monophotonSkim, BExpr.kOR)
)

output = mithep.OutputMod('Output',
    MaxFileSize = 10 * 1024, # 10 GB - should never exceed
    FileName = 'monox',
    PathName = ".",
    CheckTamBr = False,
    KeepTamBr = False,
    CheckBrDep = True,
    UseBrDep = True
)
output.Keep("*")
output.Drop("L1TechBits*")
output.Drop("L1AlgoBits*")
output.Drop("MCVertexes")
output.Drop("PFEcal*SuperClusters")
output.Drop("*Tracks")
output.Drop("StandaloneMuonTracksWVtxConstraint")
output.Drop("PrimaryVertexesBeamSpot")
output.Drop("InclusiveSecondaryVertexes")
output.Drop("CosmicMuons")
output.Drop("MergedElectronsStable")
output.Drop("MergedConversions*")
output.Drop("AKT8GenJets")
output.Drop("AKt4PFJets")
output.Drop("DCASig")
output.AddNewBranch(metCorrection.GetOutputName())
output.AddNewBranch(monojetSkim.GetCategoryFlagsName())
output.AddNewBranch(monophotonSkim.GetCategoryFlagsName())

output.AddCondition(skim)

analysis.setSequence(
    badEventsFilterMod * 
    (
        (metCorrection * jetCorrection) +
        monojetPrefilter +
        monophotonPrefilter +
        (
            prefilter * goodPVFilterMod * electronId * muonId * photonId * jetId *
            (
                monojetSkim +
                monophotonSkim +
                skim
            )
        )
    ) +
    output
)
