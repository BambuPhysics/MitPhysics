## getJECLocal
## Run this script as a CMSSW job to dump the JEC constants into text files from the conditions DB.

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

import sys

options = VarParsing('analysis')
options.register('globalTag', default = '', mult = VarParsing.multiplicity.singleton, mytype = VarParsing.varType.string, info = 'global tag')
options.register('jecTag', default = '', mult = VarParsing.multiplicity.singleton, mytype = VarParsing.varType.string, info = 'JEC tag (optional - only used for output file names)')
options.register('payload', default = 'AK4PFchs', mult = VarParsing.multiplicity.singleton, mytype = VarParsing.varType.string, info = 'DB payload name')

options.parseArguments()

if not options.globalTag:
    print 'Global tag not set'
    sys.exit(1)

if not options.jecTag:
    options.jecTag = options.globalTag

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
    globalTag      = cms.untracked.string(options.jecTag)
)

process.p = cms.Path(
    process.saveCorrections
)
