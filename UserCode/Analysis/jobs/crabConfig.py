from WMCore.Configuration import Configuration

config = Configuration()
config.section_("General")
#config.General.requestName = setEra+setType 
#config.General.workArea = setType+setId+'_jobVer'+jobVer
config.General.workArea = 'phi_gamma_data'
config.General.requestName = 'pg_data2024'
config.General.transferLogs = True 
config.General.transferOutputs = True 

config.section_("Data")

#config.Data.inputDataset = '/ParkingDoubleMuonLowMass0/Run2023D-22Sep2023_v1-v1/MINIAOD'
config.Data.inputDataset = '/ParkingDoubleMuonLowMass0/Run2024B-PromptReco-v1/MINIAOD'

config.Data.lumiMask= '2024B_Golden.json'
#config.Data.lumiMask= 'Cert_Collisions2023_366442_370790_Golden.json'


#config.Data.lumiMask='Cert_Collisions2022_eraG_362433_362760_Golden.json'
#config.Data.lumiMask='Cert_Collisions2022_eraD_357538_357900_Golden.json'
#config.Data.lumiMask='Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
#config.Data.lumiMask='Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'

#config.Data.runRange = '362719-362760'
config.Data.runRange = '379000-379200'

config.Data.useParent = False 
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
#config.Data.splitting = 'Automatic'
#config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 100 #number of files per jobs
config.Data.totalUnits =  -1 #number of event
#config.Data.outLFNDirBase = '/store/user/konec/crabout/'
config.Data.outLFNDirBase = '/store/user/kszlezak/crabout/'
config.Data.outputDatasetTag = 'phi_gamma3'

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'dataAnalysis.py'
#config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles = ['dataTestJob.root']

config.section_("Site")
#config.Site.whitelist = ['T3_CH_CERNCAF']
#config.Site.whitelist = ['T2_CH_CERN']
#config.Site.storageSite = 'T2_PL_Swierk'
config.Site.storageSite = 'T3_CH_CERNBOX'
#config.Site.blacklist = ['T2_KR_*','T2_CN_*','T2_BR_*','T2_US_Florida','T2_US_UCSD']
