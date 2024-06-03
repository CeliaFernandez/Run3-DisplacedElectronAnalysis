import FWCore.ParameterSet.Config as cms
import os


process = cms.Process("demo")
#process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.Services_cff')

# Debug printout and summary.
process.load("FWCore.MessageService.MessageLogger_cfi")

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True),
  # Set up multi-threaded run. Must be consistent with config.JobType.numCores in crab_cfg.py.
  #numberOfThreads=cms.untracked.uint32(8)
)

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag

# Select number of events to be processed
nEvents = 1000
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents) )

# Read events
#listOfFiles = ['/store/data/Run2022F/EGamma/MINIAOD/PromptReco-v1/000/360/389/00000/02d37797-3eda-47c4-bba9-1f28b1402dd7.root']
listOfFiles = ['file:/eos/user/f/fernance/Run3-Analyses/SnT-Supercluster/CMSSW_12_4_11_patch3/src/EXO-Run3Summer22EEMiniAODv3-00772.root']
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring( listOfFiles ),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents = cms.untracked.uint32(0)
  )
process.GlobalTag = GlobalTag(process.GlobalTag, '124X_mcRun3_2022_realistic_v5')

## Define the process to run 
## 
process.load("Analysis.Run3-SuperclusterAnalysis.ntuples_cfi")

process.p = cms.EndPath(process.ntuples)

