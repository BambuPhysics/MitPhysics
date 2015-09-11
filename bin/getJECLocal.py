## getJECLocal
## Run this script as a CMSSW job to dump the JEC constants into text files from sqlite files provided by the JetMET group.

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

import sys

options = VarParsing('analysis')
options.register('source', default = '', mult = VarParsing.multiplicity.singleton, mytype = VarParsing.varType.string, info = 'source sqlite file')
options.register('jecTag', default = 'Summer15_50nsV5_DATA_AK4PF', mult = VarParsing.multiplicity.singleton, mytype = VarParsing.varType.string, info = 'payload tag')
options.register('payload', default = 'AK4PF', mult = VarParsing.multiplicity.singleton, mytype = VarParsing.varType.string, info = 'DB payload name')

options.parseArguments()

if not options.source:
    print 'Source not given'
    sys.exit(1)

if '/' in options.source:
    print 'Source sqlite file must be in the current directory (cannot have / in path name).'
    sys.exit(1)

process = cms.Process("JEC")

process.load("CondCore.DBCommon.CondDBCommon_cfi")
from CondCore.DBCommon.CondDBSetup_cfi import *
process.jec = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(0)
    ),
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('JetCorrectionsRecord'),
            tag    = cms.string('JetCorrectorParametersCollection_' + options.jecTag),
            label  = cms.untracked.string(options.payload)
        ),
    ), 
    connect = cms.string('sqlite:' + options.source)
)

process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))

process.load('Configuration.StandardSequences.Services_cff')

process.saveCorrections = cms.EDAnalyzer('JetCorrectorDBReader', 
    payloadName    = cms.untracked.string(options.payload),
    printScreen    = cms.untracked.bool(True),
    createTextFile = cms.untracked.bool(True),
    globalTag      = cms.untracked.string(options.jecTag)
)

process.p = cms.Path(
    process.saveCorrections
)
