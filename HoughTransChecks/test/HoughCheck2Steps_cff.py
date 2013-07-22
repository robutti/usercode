import FWCore.ParameterSet.Config as cms

process = cms.Process("HoughCheck")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger") 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

### conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'START61_V8::All'

### standard includes
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.Geometry_cff')
process.load("Configuration.StandardSequences.RawToDigi_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

process.load("RecoTracker.TrackProducer.TrackRefitters_cff")

process.load('MLoVetere.HTTrackSeeding.HoughTransformSeedLayersMultiStep_cfi')

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
#        '/store/relval/CMSSW_6_1_0-START61_V8/RelValSingleMuPt10/GEN-SIM-RECO/v1/00000/D2D2E77F-F94C-E211-86F6-003048FFD71A.root'
#        '/store/relval/CMSSW_6_1_0-START61_V8/RelValSingleMuPt1/GEN-SIM-RECO/v1/00000/E82A4907-664D-E211-9671-003048678F8A.root'
        '/store/relval/CMSSW_6_1_0-PU_START61_V8/RelValTTbar/GEN-SIM-RECO/v1/00000/48D7BE40-DC4D-E211-A2E6-003048D373AE.root'
        )
#    skipEvents = cms.untracked.uint32(1)
    )

process.houghcheck = cms.EDAnalyzer('HoughCheck2Steps',
#                                    tracks = cms.InputTag('generalTracks'),
                                    tracks = cms.InputTag('TrackRefitter'),
                                    TTRHBuilder = cms.string('WithAngleAndTemplate'),
                                    algoSel = cms.vuint32(9, 10, 11, 12, 13, 14),
                                    seedingLayers  =  cms.string('HoughTransformSeedLayersAllHitsOneSet'),
                                    minPar = cms.untracked.vdouble(-80., -0.075, -3.14159, -450., -2.5),
                                    maxPar = cms.untracked.vdouble(80., 0.075, 3.14159, 450., 2.5),
                                    nBins = cms.untracked.vint32(4, 50, 200, 6, 50),
                                    phiBinOverlap = cms.double(0.25),
                                    etaBinOverlap = cms.double(0.25),
                                    xyVoteThr = cms.uint32(2),
                                    xyzVoteThr = cms.uint32(2),
                                    cleanupSeeds = cms.untracked.bool(False),
                                    outRootFile = cms.untracked.string(houghCheck_2steps.root),
                                    verbosity = cms.untracked.uint32(1)
)

process.p = cms.Path(process.TrackRefitter*process.htSeedLayers*process.houghcheck)
#process.p = cms.Path(process.clustToHits*process.beamSpot*process.tracking*process.htSeedLayers*process.houghcheck)
#process.p = cms.Path(process.htSeedLayers*process.houghcheck)
