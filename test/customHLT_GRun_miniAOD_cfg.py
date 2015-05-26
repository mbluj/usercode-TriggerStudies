# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: customHLT --conditions MCRUN2_72_V1 -n 10 --step=HLT:User --eventcontent FEVTDEBUGHLT --datatier GEN-SIM-DIGI-RAW --no_exec --processName=TauHLT --filein /store/relval/CMSSW_7_2_1/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_72_V1-v1/00000/04949CA1-275D-E411-99B8-0025905B855C.root --fileout file:customHLT.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('TauHLT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
#process.load('HLTrigger.Configuration.HLT_User_cff') #Custom menu
#process.load('HLTrigger.Configuration.HLT_GRun_cff') #GRun from release
process.load('TriggerStudies.Tau.HLT_GRun_cff') #GRun from a custom dump
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/00A5809F-5581-E411-AFF4-0025905A6136.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/045F75B8-5781-E411-BFA6-0025905A60D0.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/0A1CECDC-5881-E411-9B6A-0025905B85F6.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/0ACE805D-5581-E411-AC6E-0025905A6092.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/1847904D-5681-E411-B92A-0025905B8572.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/1A44C3BD-5481-E411-98EB-003048FFCB6A.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/2A13600A-5581-E411-93FC-0025905A60B2.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/2E3EFE40-5681-E411-92F3-002618943924.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/36332F45-5681-E411-B7D3-0025905A60EE.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/3807A247-5681-E411-82FC-0025905A6076.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/4045DB40-5681-E411-84D1-0026189438B0.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/42A8AEA4-5581-E411-A26C-0025905A60E0.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/4867FBB3-5781-E411-81AF-0025905B8582.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/4EE9CD99-5381-E411-AA94-0025905964C4.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/524FA19A-5681-E411-91CF-0025905B8610.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/585BC3C4-5781-E411-9AD4-0025905A606A.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/66CD8D09-5581-E411-B531-0025905A60A0.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/6AB5F9D5-5381-E411-8CBB-0025905A48FC.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/6EB012A3-5581-E411-8813-0025905A611E.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/7801EEA3-5581-E411-9AC0-002590596498.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/80F522D9-5381-E411-8EA4-0025905A60B4.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/844CDEB9-5781-E411-BE80-0025905A605E.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/8EFFA798-5681-E411-B6AB-0025905A6092.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/98ACE28B-5381-E411-9BC5-002618943985.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/9A6DA044-5681-E411-8EEF-0025905B8596.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/9CB319F2-5581-E411-A754-0025905B857C.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/A0B8DC9E-5581-E411-BCC7-0025905A6092.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/A26504D4-5381-E411-8710-002618943904.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/B03D9457-5581-E411-BF9A-0025905A609A.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/B0A2C39E-5681-E411-99EC-0025905A6066.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/B610038F-5381-E411-9DB2-002618943933.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/BC75F047-5681-E411-A433-0025905B85EE.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/C4F39099-5681-E411-BCDE-0025905A6094.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/C681EB45-5681-E411-BAF3-0025905B860C.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/CC868C8D-5381-E411-8C39-002618B27F8A.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/CCC6959B-5681-E411-A8DD-0025905A612C.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/D8E2309A-5381-E411-A4B3-002590596486.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/DA3571F8-5581-E411-A2AF-0025905A6090.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/E0F355A0-5581-E411-9961-0025905B85D0.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/E2D59DB3-5781-E411-BF98-0025905B8576.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/E42F9C56-5581-E411-AF2E-0025905B860E.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/F2F6073F-5481-E411-8D52-0025905964B6.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/F453160D-5581-E411-A7CE-0025905964B4.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/F6FAEA42-5681-E411-ABD3-002618FDA204.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-DIGI-RAW-HLTDEBUG/PU25ns_MCRUN2_73_V7-v1/00000/FEE80348-5681-E411-8C5E-0025905A48F2.root"
    ),
    fileNames = cms.untracked.vstring(
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO//PU25ns_MCRUN2_73_V7-v1/00000/1EF54505-8081-E411-8C29-00261894393A.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO//PU25ns_MCRUN2_73_V7-v1/00000/2431F681-7581-E411-91D3-00261894397E.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO//PU25ns_MCRUN2_73_V7-v1/00000/4CBEE6E7-7881-E411-9E77-0025905A608A.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO//PU25ns_MCRUN2_73_V7-v1/00000/4E455D6D-7681-E411-8436-0025905B8590.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO//PU25ns_MCRUN2_73_V7-v1/00000/6E61F881-7B81-E411-9651-0025905A6076.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO//PU25ns_MCRUN2_73_V7-v1/00000/7EA06BE5-7381-E411-A18C-0025905A48B2.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO//PU25ns_MCRUN2_73_V7-v1/00000/AEA42B0C-7281-E411-97BA-0025905A613C.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO//PU25ns_MCRUN2_73_V7-v1/00000/B0448D2D-7581-E411-B208-0025905A608A.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO//PU25ns_MCRUN2_73_V7-v1/00000/C87041F1-7281-E411-B8D0-0025905A609A.root",
        "/store/relval/CMSSW_7_3_0/RelValQQH1352T_Tauola_13/GEN-SIM-RECO//PU25ns_MCRUN2_73_V7-v1/00000/D83DCA84-7481-E411-83F5-002618943956.root"
        )
)

##################################
# GlobalTag
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'MCRUN2_73_V7', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'PHYS14_25_V1', '')

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(False),
    wantSummary = cms.untracked.bool(True),                                                                                 
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.FwkReport.reportEvery = 1

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.0.pre0 $'),
    annotation = cms.untracked.string('miniAOD nevts:1'),
    name = cms.untracked.string('Applications')
)

##################################
# Add scheduled MiniAOD sequences
from TriggerStudies.Tau.miniAOD_my_PAT_cff import *

# Output definition
addMiniAODOutput(process, outFile='miniAOD.root')
process.MINIAODSIMoutput.outputCommands.append('keep *_TriggerResults_*_*')
process.MINIAODSIMoutput.outputCommands.append('keep *_hltL1extraParticles_*_*')

# PAT & MiniAOD scheduled sequences
addMiniAODSched(process)
if 'patTrigger' in process.__dict__:
    process.patTrigger.processName = cms.string('TauHLT')
    process.patAndMiniAODPath.remove(process.patTrigger)
    process.MINIAODSIMoutput_step.replace(
        process.MINIAODSIMoutput,
        process.patTrigger+process.MINIAODSIMoutput        
    )
if 'patTriggerEvent' in process.__dict__:
    process.patTriggerEvent.processName = cms.string('TauHLT')
    process.patAndMiniAODPath.remove(process.patTriggerEvent)
    process.MINIAODSIMoutput_step.replace(
        process.MINIAODSIMoutput,
        process.patTriggerEvent+process.MINIAODSIMoutput        
    )
if 'selectedPatTrigger' in process.__dict__:
    process.patAndMiniAODPath.remove(process.selectedPatTrigger)
    process.MINIAODSIMoutput_step.replace(
        process.MINIAODSIMoutput,
        process.selectedPatTrigger+process.MINIAODSIMoutput        
    )



# Other statements
from HLTrigger.Configuration.CustomConfigs import ProcessName
process = ProcessName(process)


##################################
# Schedule definition
process.schedule = cms.Schedule()
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.patAndMiniAODPath])
process.schedule.extend([process.endjob_step,process.MINIAODSIMoutput_step])

##################################
# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

#MB ->
# load PostLS1 customisation
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1
process = customisePostLS1(process)

# adapt HLT modules to the correct process name
if 'hltTrigReport' in process.__dict__:
    process.hltTrigReport.HLTriggerResults             = cms.InputTag( 'TriggerResults', '', 'TauHLT' )

if 'hltPreExpressCosmicsOutputSmart' in process.__dict__:
    process.hltPreExpressCosmicsOutputSmart.hltResults = cms.InputTag( 'TriggerResults', '', 'TauHLT' )

if 'hltPreExpressOutputSmart' in process.__dict__:
    process.hltPreExpressOutputSmart.hltResults        = cms.InputTag( 'TriggerResults', '', 'TauHLT' )

if 'hltPreDQMForHIOutputSmart' in process.__dict__:
    process.hltPreDQMForHIOutputSmart.hltResults       = cms.InputTag( 'TriggerResults', '', 'TauHLT' )

if 'hltPreDQMForPPOutputSmart' in process.__dict__:
    process.hltPreDQMForPPOutputSmart.hltResults       = cms.InputTag( 'TriggerResults', '', 'TauHLT' )

if 'hltPreHLTDQMResultsOutputSmart' in process.__dict__:
    process.hltPreHLTDQMResultsOutputSmart.hltResults  = cms.InputTag( 'TriggerResults', '', 'TauHLT' )

if 'hltPreHLTDQMOutputSmart' in process.__dict__:
    process.hltPreHLTDQMOutputSmart.hltResults         = cms.InputTag( 'TriggerResults', '', 'TauHLT' )

if 'hltPreHLTMONOutputSmart' in process.__dict__:
    process.hltPreHLTMONOutputSmart.hltResults         = cms.InputTag( 'TriggerResults', '', 'TauHLT' )

if 'hltDQMHLTScalers' in process.__dict__:
    process.hltDQMHLTScalers.triggerResults            = cms.InputTag( 'TriggerResults', '', 'TauHLT' )
    process.hltDQMHLTScalers.processname               = 'TauHLT'

if 'hltDQML1SeedLogicScalers' in process.__dict__:
    process.hltDQML1SeedLogicScalers.processname       = 'TauHLT'


#MB <-
#MB ->
# override the L1 menu from an Xml file
process.l1GtTriggerMenuXml = cms.ESProducer("L1GtTriggerMenuXmlProducer",
  TriggerMenuLuminosity = cms.string('startup'),
  DefXmlFile = cms.string('L1Menu_Collisions2015_25ns_v2_L1T_Scales_20141121_Imp0_0x1030.xml'),
  VmeXmlFile = cms.string('')
)
process.L1GtTriggerMenuRcdSource = cms.ESSource("EmptyESSource",
  recordName = cms.string('L1GtTriggerMenuRcd'),
  iovIsRunNotTime = cms.bool(True),
  firstValid = cms.vuint32(1)
)
process.es_prefer_l1GtParameters = cms.ESPrefer('L1GtTriggerMenuXmlProducer','l1GtTriggerMenuXml')

# customize the L1 emulator to run customiseL1EmulatorFromRaw with HLT to switchToSimStage1Digis
process.load( 'Configuration.StandardSequences.RawToDigi_cff' )
process.load( 'Configuration.StandardSequences.SimL1Emulator_cff' )
import L1Trigger.Configuration.L1Trigger_custom
#                    
 
# 2015 Run2 emulator
import L1Trigger.L1TCalorimeter.L1TCaloStage1_customForHLT
process = L1Trigger.L1TCalorimeter.L1TCaloStage1_customForHLT.customiseL1EmulatorFromRaw( process )
#
process = L1Trigger.Configuration.L1Trigger_custom.customiseResetPrescalesAndMasks( process )
# customize the HLT to use the emulated results
import HLTrigger.Configuration.customizeHLTforL1Emulator
process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToL1Emulator( process )
process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToSimStage1Digis( process )
#MB <-
#MB ->
#CSC unpacking
process.cscReEmulTriggerPrimitiveDigis.CSCComparatorDigiProducer = cms.InputTag("simMuonCSCDigis","MuonCSCComparatorDigi")
#MB <-

# End of customisation functions


