import FWCore.ParameterSet.Config as cms

process = cms.Process("HoughCheck")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START60_V4::All'

### standard includes
process.load('Configuration/StandardSequences/Services_cff')
process.load('Configuration.StandardSequences.GeometryPilot2_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_6_1_0-START61_V8/RelValSingleMuPt10/GEN-SIM-RECO/v1/00000/D2D2E77F-F94C-E211-86F6-003048FFD71A.root'
#    '/store/relval/CMSSW_6_1_0-START61_V8/RelValSingleMuPt1/GEN-SIM-RECO/v1/00000/E82A4907-664D-E211-9671-003048678F8A.root'
#    '/store/relval/CMSSW_6_1_0-PU_START61_V8/RelValTTbar/GEN-SIM-RECO/v1/00000/48D7BE40-DC4D-E211-A2E6-003048D373AE.root'
    )
#    skipevents = cms.untracked.uint32(0)
)

process.houghcheck = cms.EDAnalyzer('HoughCheckOnTracks',
        tracks = cms.untracked.InputTag('TrackRefitter'),
        TTRHBuilder = cms.string('WithAngleAndTemplate'),
        algoSel = cms.vuint32(),
        layerSel = cms.vuint32(),
        maxDoca = cms.untracked.double(100.),
        nBinsDoca = cms.untracked.int32(200),
        maxKappa = cms.untracked.double(0.002),
        nBinsKappa = cms.untracked.int32(200),
        nBinsPhi = cms.untracked.int32(200),
        maxAbsZ0 = cms.untracked.double(265.),
        nBinsZ0 = cms.untracked.int32(1),
        maxAbsLambda = cms.untracked.double(1.4),
        nBinsLambda = cms.untracked.int32(1),
        verbosity = cms.untracked.uint32(2)
)


process.p = cms.Path(process.TrackRefitter*process.houghcheck)
