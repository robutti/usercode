import FWCore.ParameterSet.Config as cms

process = cms.Process("HoughCheck")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger = cms.Service("MessageLogger") 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

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
        '/store/relval/CMSSW_6_2_0/RelValTTbar/GEN-SIM-RECO/PU_PRE_ST62_V8-v2/00000/06277119-A9EC-E211-B39F-003048F1BFE0.root'
        )
#    skipEvents = cms.untracked.uint32(1)
    )

process.houghcheck = cms.EDAnalyzer('HoughCheck2Steps',
#                                    tracks = cms.InputTag('generalTracks'),
                                    tracks = cms.InputTag('TrackRefitter'),
                                    TTRHBuilder = cms.string('WithAngleAndTemplate'),
#                                    algoSel = cms.vuint32(9, 10, 11, 12, 13, 14),
                                    seedingLayers  =  cms.string('HoughTransformSeedLayersAllHitsOneSet'),
                                    minPar = cms.untracked.vdouble(-80., -0.075, -3.14159, -450., -2.5),
                                    maxPar = cms.untracked.vdouble(80., 0.075, 3.14159, 450., 2.5),
                                    nBins = cms.untracked.vint32(2, 2, 2, 2, 2),
                                    phiBinOverlap = cms.double(0.25),
                                    etaBinOverlap = cms.double(0.25),
                                    xyVoteThr = cms.uint32(2),
                                    xyzVoteThr = cms.uint32(2),
                                    cleanupSeeds = cms.untracked.bool(False),
                                    outRootFile = cms.untracked.string('houghCheck_2steps.root'),
                                    verbosity = cms.untracked.uint32(3)
)

process.p = cms.Path(process.TrackRefitter*process.htSeedLayers*process.houghcheck)
#process.p = cms.Path(process.clustToHits*process.beamSpot*process.tracking*process.htSeedLayers*process.houghcheck)
#process.p = cms.Path(process.htSeedLayers*process.houghcheck)
