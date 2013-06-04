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

process.houghcheck = cms.EDAnalyzer('HoughCheckAllHits',
#                                    tracks = cms.InputTag('generalTracks'),
                                    tracks = cms.InputTag('TrackRefitter'),
                                    algoSel = cms.vuint32(5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
                                    seedingLayers  =  cms.string('HoughTransformSeedLayersAllHitsOneSet'),
                                    rangeDoca = cms.untracked.vdouble(-80., 80.),
                                    nBinsDoca = cms.untracked.int32(4),
                                    rangeSqrtK = cms.untracked.vdouble(-0.075, 0.075),
                                    nBinsSqrtK = cms.untracked.int32(50),
#                                    rangePhi = cms.untracked.vdouble(0., 3.14159),
                                    nBinsPhi = cms.untracked.int32(200),
                                    rangeZ0 = cms.untracked.vdouble(-450., 450.),
                                    nBinsZ0 = cms.untracked.int32(6),
                                    rangeEta = cms.untracked.vdouble(-2.5, 2.5),
                                    nBinsEta = cms.untracked.int32(50),
                                    verbosity = cms.untracked.uint32(1)
)

process.p = cms.Path(process.TrackRefitter*process.htSeedLayers*process.houghcheck)
#process.p = cms.Path(process.clustToHits*process.beamSpot*process.tracking*process.htSeedLayers*process.houghcheck)
#process.p = cms.Path(process.htSeedLayers*process.houghcheck)
