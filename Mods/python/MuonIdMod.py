from MitAna.TreeMod.bambu import mithep

muonIdMod = mithep.MuonIdMod(
    OutputName = 'LooseMuons',
    MuonClassType = mithep.MuonTools.kAll,
    IdType = mithep.MuonTools.kNoId,
    IsoType = mithep.MuonTools.kNoIso,
    PFNoPileupCandidatesName = 'pfNoPU',
    PFPileupCandidatesName = 'pfPU',
    ApplyD0Cut = False,
    ApplyDZCut = False,
    PtMin = 3.,
    EtaMax = 2.4
)
