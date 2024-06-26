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
nEvents = -1
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(nEvents) )

# Read events
listOfFiles = []
## Signal
#listOfFiles.append('/store/mc/Run3Summer22MiniAODv4/SquarkToNeutralinoTo2LNu-MSquark_350_MChi_150_ctau_25mm_TuneCP5_13p6TeV_pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2530000/57d6cf9e-007e-4d84-94b0-82b00c6a9314.root')
#listOfFiles.append('/store/mc/Run3Summer22MiniAODv4/SquarkToNeutralinoTo2LNu-MSquark_350_MChi_150_ctau_25mm_TuneCP5_13p6TeV_pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2530000/a67c5e70-431d-4e65-a6d9-5af6e7645ddb.root')
#listOfFiles.append('/store/mc/Run3Summer22MiniAODv4/SquarkToNeutralinoTo2LNu-MSquark_350_MChi_150_ctau_25mm_TuneCP5_13p6TeV_pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/30000/315bb0b7-3fac-4cec-a4ae-10d9edb896b7.root')
#listOfFiles.append('/store/mc/Run3Summer22MiniAODv4/SquarkToNeutralinoTo2LNu-MSquark_350_MChi_150_ctau_25mm_TuneCP5_13p6TeV_pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/30000/8e290069-3f1c-4d40-96d0-e29fd70dbd92.root')
#listOfFiles.append('/store/mc/Run3Summer22MiniAODv4/SquarkToNeutralinoTo2LNu-MSquark_350_MChi_150_ctau_25mm_TuneCP5_13p6TeV_pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/30000/96cf6b75-1b1e-48ab-9d5b-75daeceb0b95.root')
#listOfFiles.append('/store/mc/Run3Summer22MiniAODv4/SquarkToNeutralinoTo2LNu-MSquark_350_MChi_150_ctau_25mm_TuneCP5_13p6TeV_pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/30000/96fffa9c-c057-41f8-86a0-516f2e3162ca.root')
#listOfFiles.append('/store/mc/Run3Summer22MiniAODv4/SquarkToNeutralinoTo2LNu-MSquark_350_MChi_150_ctau_25mm_TuneCP5_13p6TeV_pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/40000/793356db-61d5-47e3-a7f0-2d73849055f9.root')
#listOfFiles.append('/store/mc/Run3Summer22MiniAODv4/SquarkToNeutralinoTo2LNu-MSquark_350_MChi_150_ctau_25mm_TuneCP5_13p6TeV_pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/40000/99fba257-a893-48ef-b71f-d9cbb9857377.root')
#listOfFiles.append('/store/mc/Run3Summer22MiniAODv4/SquarkToNeutralinoTo2LNu-MSquark_350_MChi_150_ctau_25mm_TuneCP5_13p6TeV_pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/40000/c2ff12c8-dec4-465f-a939-978755861aef.root')
#listOfFiles.append('/store/mc/Run3Summer22MiniAODv4/SquarkToNeutralinoTo2LNu-MSquark_350_MChi_150_ctau_25mm_TuneCP5_13p6TeV_pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/50000/68f92996-453d-412e-9141-2c75830cfff9.root')
#listOfFiles.append('/store/mc/Run3Summer22MiniAODv4/SquarkToNeutralinoTo2LNu-MSquark_350_MChi_150_ctau_25mm_TuneCP5_13p6TeV_pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/50000/f8176527-8861-4f76-8765-8c4b58f5c1f8.root')
## Drell-Yan
listOfFiles.append('/store/mc/Run3Summer22MiniAODv4/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/30000/38a57f5e-e24f-4b08-a1a4-b83501c9aa02.root')
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

