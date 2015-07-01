import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

import sys

options = VarParsing('analysis')
options.register('globalTag', default = '', mult = VarParsing.multiplicity.singleton, mytype = VarParsing.varType.string, info = 'global tag')
options.register('payload', default = 'AK4PFchs', mult = VarParsing.multiplicity.singleton, mytype = VarParsing.varType.string, info = 'DB payload name')

options.parseArguments()

if not options.globalTag:
    print 'Global tag not set'
    sys.exit(1)

process = cms.Process("JEC")

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = options.globalTag

process.saveCorrections = cms.EDAnalyzer('JetCorrectorDBReader', 
    payloadName    = cms.untracked.string(options.payload),
    printScreen    = cms.untracked.bool(True),
    createTextFile = cms.untracked.bool(True),
    globalTag      = cms.untracked.string(options.globalTag)
)

process.p = cms.Path(
    process.saveCorrections
)
