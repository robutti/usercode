import FWCore.ParameterSet.Config as cms

process = cms.Process("HoughCheck")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger") 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')

### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

#process.load('MLoVetere.HTTrackSeeding.HoughTransformSeedLayersSingleStep_cfi')
process.load('ERobutti.HoughTransChecks.HoughTransformSeedLayersMatched_cfi')

process.clustToHits = cms.Sequence(
    process.siPixelRecHits*process.siStripMatchedRecHits
    )

process.beamSpot = cms.Sequence(
    process.offlineBeamSpot
    )

process.tracking = cms.Sequence(
    process.trackingGlobalReco
    )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/relval/CMSSW_6_2_0/RelValSingleMuPt1/GEN-SIM-RECO/PRE_ST62_V8-v3/00000/22B9BCD6-60EC-E211-924E-003048FEAEB0.root'
#        '/store/relval/CMSSW_6_2_0/RelValSingleMuPt10/GEN-SIM-RECO/PRE_ST62_V8-v3/00000/FEB7D35C-5CEC-E211-80AA-003048FEB8EE.root'
#        '/store/relval/CMSSW_6_2_0/RelValTTbar/GEN-SIM-RECO/PU_PRE_ST62_V8-v2/00000/06277119-A9EC-E211-B39F-003048F1BFE0.root'
        )
#    skipEvents = cms.untracked.uint32(1)
    )

process.houghcheck = cms.EDAnalyzer('HoughCheckStubs',
#                                    tracks = cms.InputTag('generalTracks'),
                                    tracks = cms.InputTag('TrackRefitter'),
                                    TTRHBuilder = cms.string('WithAngleAndTemplate'),
#                                    algoSel = cms.vuint32(9, 10, 11, 12, 13, 14),
                                    seedingLayers  =  cms.string('HoughTransformSeedLayersPixelsAndMatchedOneSet'),
#                                    minPar = cms.untracked.vdouble(-0.5, -0.005, -3.14159, -50., -2.5),
#                                    maxPar = cms.untracked.vdouble(0.5, 0.005, 3.14159, 50., 2.5),
                                    nBins = cms.untracked.vint32(20, 20, 20, 20, 20),
                                    binOverlap = cms.double(0.25),
                                    voteThr = cms.uint32(2),
                                    cleanupSeeds = cms.untracked.bool(False),
                                    outRootFile = cms.untracked.string('houghCheck_stubs.root'),
                                    verbosity = cms.untracked.uint32(1)
)

process.p = cms.Path(process.TrackRefitter*process.htSeedLayers*process.houghcheck)
#process.p = cms.Path(process.clustToHits*process.beamSpot*process.tracking*process.htSeedLayers*process.houghcheck)
#process.p = cms.Path(process.htSeedLayers*process.houghcheck)
