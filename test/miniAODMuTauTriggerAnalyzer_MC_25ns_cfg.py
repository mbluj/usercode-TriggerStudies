# example from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD

import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

isMC = True

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) 
    #input = cms.untracked.int32(10000) 
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/tmp/mbluj/miniAOD-prod_PAT_GRun_my.root',
        #'root://se.cis.gov.pl:1094//dpm/cis.gov.pl/home/cms/store/user/bluj/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/GGH125ToTauTau_MiniAOD_GRunV47_v3/c72d2c1790a6b87324592002758e580a/miniAOD_10_1_og1.root'
	#
        #'root://xrootd-cms.infn.it//store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/254/833/00000/FCB9CBAE-104B-E511-BF89-02163E011F65.root',
	#
	'root://cms-xrd-global.cern.ch///store/mc/RunIISpring15DR74/WJetsToLNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v1/00000/048FB1EE-33FD-E411-A2BA-0025905A6094.root'
    )
)

#myJSON = 'prod/Cert_246908-251252_13TeV_PromptReco_Collisions15_JSON.txt'
#import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = myJSON).getVLuminosityBlockRange()
#print "WARNING: JSON imported to source for running w/o Crab" 

process.evtCounter = cms.EDAnalyzer(
    "BasicGenEventInfoAnalyzer",
    isMC = cms.bool(isMC),
    genEvtInfo = cms.InputTag("generator")
)


process.muLooseTau = cms.EDAnalyzer(
    "MiniAODMuTauTriggerAnalyzer",
    bits = cms.InputTag("TriggerResults","","HLT"), # use correct process name
    prescales = cms.InputTag("patTrigger"), 
    objects = cms.InputTag("selectedPatTrigger"),
    muons = cms.InputTag("slimmedMuons"),
    #muons = cms.InputTag("preSelectedMuons"),
    #taus = cms.InputTag("slimmedTaus"),
    taus = cms.InputTag("preSelectedTaus"),
    l1CenJets = cms.InputTag("l1extraParticles:Central"), #FIXME: from RECO      
    l1TauJets = cms.InputTag("l1extraParticles:Tau"), #FIXME: from RECO      
    l1IsoTaus = cms.InputTag("l1extraParticles:IsoTau"), #FIXME: from RECO (missing in Phys14 and 50ns samples)
    met = cms.InputTag("slimmedMETs"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muonTriggers = cms.vstring(  #version number is ignored, so can be replaced by wildcard (*) or dropped
        "HLT_IsoMu24_eta2p1_v*", # unp'ed at 1.4e34
        "HLT_IsoMu20_eta2p1_v*", # unp'ed at 7e33
        "HLT_IsoMu17_eta2p1_v*", # unp'ed at 5e33
        "HLT_IsoMu16_eta2p1_CaloMET30_v*",
    ),
    muonFilters = cms.vstring(
        "hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09", # HLT_IsoMu24_eta2p1_v1
        "hltL3crIsoL1sMu16Eta2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p09", # HLT_IsoMu20_eta2p1_v1
        "hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09",  # HLT_IsoMu17_eta2p1_v1
        "hltL3crIsoL1sMu14erETM30L1f0L2f14QL3f10QL3L3trkIsoFiltered0p09", # HLT_IsoMu16_eta2p1_CaloMET30_v
    ),
    tauFilters = cms.vstring(
        #
        "hltOverlapFilterIsoMu17LooseIsoPFTau20", # HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1
        "hltL1sMu16erTauJet20er", # L1 of HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v1
        #
        "hltOverlapFilterSingleIsoMu17LooseIsoPFTau20", # HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v1
        #
        "hltOverlapFilterIsoMu17MediumIsoPFTau40Reg", # HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v1
        "hltL1sMu16erIsoTau36er", # L1 of HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v1
        "hltL2Tau35eta2p2", # L2 of HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v1
        "hltL2IsoTau35eta2p2", # L2.5 of HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v1
        "hltOverlapFilterIsoMu17L2IsoTau35", # L2.5 of HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v1
        #
        "hltOverlapFilterIsoMu24LooseIsoPFTau20", # HLT_IsoMu24_eta2p1_LooseIsoPFTau20_v1
        ## common loosePFTau20 filters (w/o checking overlap wrt lepton)
        "hltPFTau20Track",
        "hltPFTau20TrackLooseIso",
        "hltPFTau20TrackLooseIsoAgainstMuon",
        ## common mediumPFTau40 filters (w/o checking overlap wrt lepton)
        "hltPFTau40TrackPt1Reg",
        "hltPFTau40TrackPt1MediumIsolationReg",
        "hltPFTau40TrackPt1MediumIsolationL1HLTMatchedReg",
	## filters from mu+Tau50Trk30 
        "hltSingleL2Tau35eta2p2", # L2
        "hltPFTau50Track", # L3 
        "hltPFTau50TrackPt30", # L3 TrkPt>30
        "hltPFTau50TrackPt30LooseAbsOrRelIso", # L3 Iso
    ),
    minPtMu = cms.double(17.),
    maxEtaMu = cms.double(2.1),
    isoMu = cms.double(0.15),
    minPtTau = cms.double(18.),
    maxEtaTau = cms.double(2.3),
    tauIds = cms.vstring(
        "decayModeFindingNewDMs", 
        #"byLooseCombinedIsolationDeltaBetaCorr3Hits", #applied at presel level
        "againstMuonTight3", 
    ),
    tauIdsForTrees = cms.vstring(
        "decayModeFinding",
	"decayModeFindingNewDMs",
        "byCombinedIsolationDeltaBetaCorrRaw3Hits",
	"chargedIsoPtSum",
	"neutralIsoPtSum",
	"puCorrPtSum",
        "againstElectronVLooseMVA5",
	"againstElectronLooseMVA5",
	"againstElectronMediumMVA5",
        "againstMuonLoose3",
        "againstMuonTight3",
    ),
    muTightId = cms.bool(False),
    checkMCMatch = cms.bool(isMC),
    isMC = cms.bool(isMC),
    genEvtInfo = cms.InputTag("generator"),
)

import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
process.hltFilter = hlt.hltHighLevel.clone(
    TriggerResultsTag = process.muLooseTau.bits,
    HLTPaths = process.muLooseTau.muonTriggers,
    throw = False
)

process.preSelectedMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    #cut = cms.string("pt > 20. && abs(eta) < 2.1 && isPFMuon && isGlobalMuon && userIsolation(\'pat::PfChargedAllIso\') < 1.5") # PfChargedAllIso=13= \'pfChargedAll\'which is not acceptd due to some reason
    #cut = cms.string("pt > 20. && abs(eta) < 2.1 && isPFMuon && isGlobalMuon")
    cut = cms.string("pt > 17. && abs(eta) < 2.1") #keep low threshold for mu+MET events
)

process.preSelectedTaus = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    #cut = cms.string("pt > 18. && abs(eta) < 2.3 && tauID(\'decayModeFindingNewDMs\')> 0.5 && tauID(\'againstMuonLoose3\') > 0.5 && tauID(\'chargedIsoPtSum\') < 2.0" ) # Loose <2.0, Med <1.0, Tight<0.8
    cut = cms.string("pt > 18. && abs(eta) < 2.3 && tauID(\'decayModeFindingNewDMs\')> 0.5 && tauID(\'againstMuonLoose3\') > 0.5 && tauID(\'byCombinedIsolationDeltaBetaCorrRaw3Hits\') < 2.0" ) # Loose <2.0, Med <1.0, Tight<0.8
)

process.countPreSelMuons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("preSelectedMuons")
)

process.countPreSelTaus = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("preSelectedTaus")
)

process.TFileService = cms.Service(
    "TFileService", 
    fileName = cms.string("muTauTrgAna.root") 
)

process.p0 = cms.Path(
    process.evtCounter
)

process.p = cms.Path(
    process.hltFilter +
    process.preSelectedMuons + process.countPreSelMuons +
    process.preSelectedTaus + process.countPreSelTaus +
    process.muLooseTau
)
