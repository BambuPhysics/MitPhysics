from MitAna.TreeMod.bambu import mithep

muonIdMod = mithep.MuonIdMod(
    OutputName = 'LooseMuons',
    IdType = mithep.MuonTools.kNoId,
    IsoType = mithep.MuonTools.kPFIsoBetaPUCorrected,
    PFNoPileupCandidatesName = 'pfNoPU',
    PFPileupCandidatesName = 'pfPU',
    PtMin = 10.,
    EtaMax = 2.4
)
