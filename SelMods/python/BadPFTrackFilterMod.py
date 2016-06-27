from MitAna.TreeMod.bambu import mithep

badPFTrackFilterMod = mithep.BadPFTrackFilterMod(
    MuonsName = mithep.Names.gkMuonBrn,
    PFCandidatesName = mithep.Names.gkPFCandidatesBrn,
    MinPt = 100.,
    MinTrackRelErr = 0.5,
    MinRelDPt = -0.5
)
