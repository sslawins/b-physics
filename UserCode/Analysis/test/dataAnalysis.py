#!/cvmfs/cms.cern.ch/el8_amd64_gcc12/cms/cmssw/CMSSW_14_1_1/external/el8_amd64_gcc12/bin/python3

import FWCore.ParameterSet.Config as cms
import os
import sys
import subprocess


process = cms.Process("MojaAnaliza")


# MessageLogger & co.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

#dataDir = '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/OMTF/TrackingVertexing/BsToMuMuGamma_14_0_15_patch1_25_09_2024/TSG-Run3Summer22EEGS_Run2022_BsToMuMuGamma_14_0_15_patch1_25_09_2024/BsToMuMuGamma_14_0_15_patch1_25_09_2024/240925_134701/0000/'
dataDir = 'root://cms-xrd-global.cern.ch//store/data/Run2022C/ParkingDoubleMuonLowMass0/AOD/10Dec2022-v2/2540000/'

lsCommand = 'ls -1 ' + dataDir + '| grep root'
#print('Command: ', lsCommand)



dir = subprocess.Popen(lsCommand, stdout=subprocess.PIPE, shell=True, text=True)
lsOutput = dir.communicate()[0]

files = []
for f in lsOutput.split():
    #print(dataDir + f)
    files.append(dataDir + f)  # Full path to the files with 'file:' prefix


files = ['root://cms-xrd-global.cern.ch//store/data/Run2022C/ParkingDoubleMuonLowMass0/AOD/10Dec2022-v2/2540000/019a2652-4483-484d-9892-0887e18413ec.root',
         'root://cms-xrd-global.cern.ch//store/data/Run2022C/ParkingDoubleMuonLowMass0/AOD/10Dec2022-v2/2540000/01da8b6d-946e-4331-9ef5-bb12e6d9a94f.root',
         'root://cms-xrd-global.cern.ch//store/data/Run2022C/ParkingDoubleMuonLowMass0/AOD/10Dec2022-v2/2540000/02cb4e7e-bf61-4038-8ecd-223583720404.root',
         'root://cms-xrd-global.cern.ch//store/data/Run2022C/ParkingDoubleMuonLowMass0/AOD/10Dec2022-v2/2540000/0383f1aa-aa59-4e47-abb0-3d7391213bd8.root',
         'root://cms-xrd-global.cern.ch//store/data/Run2022C/ParkingDoubleMuonLowMass0/AOD/10Dec2022-v2/2540000/056a1f3c-8c9e-4869-b3df-51f81cc517df.root']

print('Number of files: ', len(files))

# input files (up to 255 files accepted)

#process.source = cms.Source('PoolSource', fileNames =cms.untracked.vstring("file:data.root") )
process.source = cms.Source('PoolSource', fileNames =cms.untracked.vstring(files) )
process.source.skipEvents = cms.untracked.uint32(0)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000))

process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.Geometry.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("Configuration.StandardSequences.Reconstruction_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_data', '')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.MessageLogger.suppressWarning  = cms.untracked.vstring('Geometry','AfterSource','L1T')
process.options = cms.untracked.PSet( wantSummary=cms.untracked.bool(False))

process.analiza= cms.EDAnalyzer("DataAnalysis",
  muonSrc = cms.InputTag("slimmedMuons"),
  candidateSrc = cms.InputTag("packedPFCandidates"),
  displacedSrc = cms.InputTag("slimmedDisplacedMuons"),
  outHist = cms.string('histosData.root')
)

process.MyPath = cms.Path(process.analiza)
process.schedule = cms.Schedule(process.MyPath)


print("All files set for analysis.")

