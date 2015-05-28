# example from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD

import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
)
process.maxEvents = cms.untracked.PSet( 
    input = cms.untracked.int32(-1) 
    #input = cms.untracked.int32(100) 
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/tmp/mbluj/miniAOD-prod_PAT_GRun_my.root',
        #'root://se.cis.gov.pl:1094//dpm/cis.gov.pl/home/cms/store/user/bluj/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/GGH125ToTauTau_MiniAOD_GRunV47_v1/8511c3cd7b71dd75a559c65d182442f5/miniAOD_100_1_Xv6.root'
        #'root://se.cis.gov.pl:1094//dpm/cis.gov.pl/home/cms/store/user/bluj/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/GGH125ToTauTau_MiniAOD_GRunV47_v2/6b3acb073896b48a28b982ccc80b3330/miniAOD_100_1_nMQ.root'
        'root://se.cis.gov.pl:1094//dpm/cis.gov.pl/home/cms/store/user/bluj/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/GGH125ToTauTau_MiniAOD_GRunV47_v3/c72d2c1790a6b87324592002758e580a/miniAOD_10_1_og1.root'
    )
)

process.eLooseTau = cms.EDAnalyzer(
    "MiniAODETauTriggerAnalyzer",
    #bits = cms.InputTag("TriggerResults","","HLT"), # use correct process name
    bits = cms.InputTag("TriggerResults","","TauHLT"), # use correct process name
    prescales = cms.InputTag("patTrigger"), 
    objects = cms.InputTag("selectedPatTrigger"),
    #electrons = cms.InputTag("slimmedElectrons"),
    electrons = cms.InputTag("preSelectedElectrons"),
    #taus = cms.InputTag("slimmedTaus"),
    taus = cms.InputTag("preSelectedTaus"),
    #l1CenJets = cms.InputTag("l1extraParticles:Central"), #FIXME: from RECO
    #l1TauJets = cms.InputTag("l1extraParticles:Tau"), #FIXME: from RECO
    #l1IsoTaus = cms.InputTag("l1extraParticles:IsoTau"), #FIXME: from RECO (missing in Phys14 and 50ns samples)
    l1CenJets = cms.InputTag("hltL1extraParticles:Central"), #FIXME: from TauHLT (missing in v2)
    l1TauJets = cms.InputTag("hltL1extraParticles:Tau"), #FIXME: from TauHLT (missing in v2)
    l1IsoTaus = cms.InputTag("hltL1extraParticles:IsoTau"), #FIXME: from TauHLT (missing in v2 and 50ns samples)
    met = cms.InputTag("slimmedMETs"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    electronTriggers = cms.vstring(  #version number is ignored, so can be replaced by wildcard (*) or dropped
        "HLT_Ele27_eta2p1_WP75_Gsf_v*", # unp'ed at 7e33
        "HLT_Ele32_eta2p1_WP75_Gsf_v*", # unp'ed at 1.4e34
        #"HLT_Ele27_eta2p1_WP85_Gsf_v*", # Phys'14
        #"HLT_Ele32_eta2p1_WP85_Gsf_v*", # Phys'14
    ),
    electronFilters = cms.vstring(
        "hltEle27WP75GsfTrackIsoFilter", # HLT_Ele27_eta2p1_WP75_Gsf_v1
        "hltEle32WP75GsfTrackIsoFilter", # HLT_Ele32_eta2p1_WP75_Gsf_v1
        #"hltEle27WP85GsfTrackIsoFilter", # HLT_Ele27_eta2p1_WP85_Gsf_v1
        #"hltEle32WP85GsfTrackIsoFilter", # HLT_Ele32_eta2p1_WP85_Gsf_v1
    ),
    tauFilters = cms.vstring(
        #
        #"hltOverlapFilterIsoEle22WP85GsfLooseIsoPFTau20", # HLT_Ele22_eta2p1_WP85_Gsf_LooseIsoPFTau20_v1, Phys'14
        #"hltL1sL1IsoEG20erTauJet20er", # L1 of HLT_Ele22_eta2p1_WP85_Gsf_LooseIsoPFTau20_v1, Phys'14
        "hltOverlapFilterIsoEle22WP75GsfLooseIsoPFTau20", # HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v1
        "hltL1sL1IsoEG20erTauJet20er", # L1 of HLT_Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20_v1
        #        
        "hltOverlapFilterIsoEle27WP75GsfLooseIsoPFTau20", # HLT_Ele27_eta2p1_WP75_Gsf_LooseIsoPFTau20_v1
        #
        "hltOverlapFilterIsoEle32WP75GsfLooseIsoPFTau20", # HLT_Ele32_eta2p1_WP75_Gsf_LooseIsoPFTau20_v1
        ## common loosePFTau20 filters (w/o checking overlap wrt lepton)
        "hltPFTau20Track",
        "hltPFTau20TrackLooseIso",
    ),
    minPtEl = cms.double(30.),
    maxEtaEl = cms.double(2.1),
    isoEl = cms.double(0.15),
    elId = cms.string("NoImplemented"),
    minPtTau = cms.double(18.),
    maxEtaTau = cms.double(2.3),
    tauIds = cms.vstring(
        "decayModeFinding", 
        #"byLooseCombinedIsolationDeltaBetaCorr3Hits", #applied at presel level
        "againstElectronLoose", #FIXME better anti-e should be used when ready 
    ),
    checkMCMatch = cms.bool(True)
)

import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
process.hltFilter = hlt.hltHighLevel.clone(
    TriggerResultsTag = process.eLooseTau.bits,
    HLTPaths = process.eLooseTau.electronTriggers,
    throw = False
)

process.preSelectedElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt > 30. && abs(eta) < 2.1")
)

process.preSelectedTaus = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string("pt > 18. && abs(eta) < 2.3 && tauID(\'decayModeFinding\')> 0.5 && tauID(\'againstMuonLoose3\') > 0.5 && tauID(\'chargedIsoPtSum\') < 2.0" ) # Loose <2.0, Med <1.0, Tight<0.8
)

process.countPreSelElectrons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("preSelectedElectrons")
)

process.countPreSelTaus = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("preSelectedTaus")
)

process.TFileService = cms.Service(
    "TFileService", 
    fileName = cms.string("elTauTrgAna.root") 
)

process.p = cms.Path(
    process.hltFilter +
    process.preSelectedElectrons + process.countPreSelElectrons +
    process.preSelectedTaus + process.countPreSelTaus +
    process.eLooseTau
)



