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


files = ['root://cms-xrd-global.cern.ch//store/data/Run2022C/ParkingDoubleMuonLowMass0/AOD/10Dec2022-v2/2540000/07e71d4d-7c8e-400b-8500-53dde9a4a700.root']

print('Number of files: ', len(files))

# input files (up to 255 files accepted)

#process.source = cms.Source('PoolSource', fileNames =cms.untracked.vstring("file:data.root") )
process.source = cms.Source('PoolSource', fileNames =cms.untracked.vstring(files) )
process.source.skipEvents = cms.untracked.uint32(0)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))

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
  processName = cms.string("HLT"),
  outHist = cms.string('histosData.root')
)

process.MyPath = cms.Path(process.analiza)
process.schedule = cms.Schedule(process.MyPath)


print("All files set for analysis.")

