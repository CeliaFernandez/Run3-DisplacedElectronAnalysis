import FWCore.ParameterSet.Config as cms

ntuples = cms.EDAnalyzer('ntuplizer',
    nameOfOutput = cms.string('output.root'),
    isData = cms.bool(True),
    EventInfo = cms.InputTag("generator"),
    RunInfo = cms.InputTag("generator"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    reducedBarrelRecHitCollection = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
    reducedEndcapRecHitCollection = cms.InputTag("reducedEgamma", "reducedEERecHits"),
    photonCollection = cms.InputTag("slimmedPhotons"),
    lowPtElectronCollection = cms.InputTag("slimmedLowPtElectrons"),
    electronCollection = cms.InputTag("slimmedElectrons"),
    primaryVertexCollection  = cms.InputTag("offlineSlimmedPrimaryVertices")
)


