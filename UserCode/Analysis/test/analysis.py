import FWCore.ParameterSet.Config as cms
import os
import sys
import subprocess

process = cms.Process("MojaAnaliza")

path = '/eos/cms/store/group/dpg_trigger/comm_trigger/L1Trigger/OMTF/TrackingVertexing/BsToMuMuGamma_14_0_15_patch1_25_09_2024/TSG-Run3Summer22EEGS_Run2022_BsToMuMuGamma_14_0_15_patch1_25_09_2024/BsToMuMuGamma_14_0_15_patch1_25_09_2024/240925_134701/0000/'
flist=[]
for file in os.listdir(path):
  if len(file) > 4 and file[-4:]=='root':
    flist.append('file:'+path+file)

# print('flist', flist)

# MessageLogger & co.
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1)
process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False))

# input files (up to 255 files accepted)
# process.source = cms.Source('PoolSource', fileNames =cms.untracked.vstring(flist))
process.source = cms.Source("PoolSource",
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODAPVv2/BdToKPi_BMuonFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/MINIAODSIM/BPH_106X_mcRun2_asymptotic_preVFP_v11-v2/2550000/220F4B68-DEFE-334E-9FCD-ECD84A0737DC.root',
#'root://xrootd-cms.infn.it//store/data/Run2023D/ParkingDoubleMuonLowMass0/MINIAOD/22Sep2023_v1-v1/2550000/0419eec5-0ae4-4732-8f06-6d72dd25a149.root',

'/store/mc/Run3Winter23MiniAOD/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/GTv3Digi_GTv3_MiniGTv3_126X_mcRun3_2023_forPU65_v3-v2/2540000/27f6ecbd-6839-49f9-86e7-b3c957ae1f46.root',
#'root://cms-xrd-global.cern.ch//store/mc/Run3Winter23MiniAOD/BsToMuMuG_MuGFilter_SoftQCDnonD_TuneCP5_13p6TeV_pythia8-evtgen/MINIAODSIM/GTv3Digi_GTv3_MiniGTv3_126X_mcRun3_2023_forPU65_v3-v2/2540000/27f6ecbd-6839-49f9-86e7-b3c957ae1f46.root',
)
)
process.source.skipEvents = cms.untracked.uint32(0)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

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

process.analiza= cms.EDAnalyzer("Analysis",
  outHist = cms.string('histos_bs.root')
)

process.MyPath = cms.Path(process.analiza)
process.schedule = cms.Schedule(process.MyPath)
