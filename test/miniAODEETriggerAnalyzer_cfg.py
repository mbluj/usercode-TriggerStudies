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
        'root://se.cis.gov.pl:1094//dpm/cis.gov.pl/home/cms/store/user/bluj/DYToEE_M-50_Tune4C_13TeV-pythia8/DYToEE_MiniAOD_GRunV47_v2/6b3acb073896b48a28b982ccc80b3330/miniAOD_10_1_CBb.root',
        'root://se.cis.gov.pl:1094//dpm/cis.gov.pl/home/cms/store/user/bluj/DYToEE_M-50_Tune4C_13TeV-pythia8/DYToEE_MiniAOD_GRunV47_v2/6b3acb073896b48a28b982ccc80b3330/miniAOD_11_1_FEP.root',
    )
)

process.eE = cms.EDAnalyzer(
    "MiniAODEETriggerAnalyzer",
    #bits = cms.InputTag("TriggerResults","","HLT"), # use correct process name
    bits = cms.InputTag("TriggerResults","","TauHLT"), # use correct process name
    prescales = cms.InputTag("patTrigger"), 
    objects = cms.InputTag("selectedPatTrigger"),
    #electrons = cms.InputTag("slimmedElectrons"),
    electrons = cms.InputTag("preSelectedElectrons"),
    met = cms.InputTag("slimmedMETs"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    tagTriggers = cms.vstring(  #version number is ignored, so can be replaced by wildcard (*) or dropped
        "HLT_Ele27_eta2p1_WP75_Gsf_v*", # unp'ed at 7e33
        "HLT_Ele32_eta2p1_WP75_Gsf_v*", # unp'ed at 1.4e34
        #"HLT_Ele27_eta2p1_WP85_Gsf_v*", # Phys'14
        #"HLT_Ele32_eta2p1_WP85_Gsf_v*", # Phys'14  
    ),
    tagFilters = cms.vstring(
        "hltEle27WP75GsfTrackIsoFilter", # HLT_Ele27_eta2p1_WP75_Gsf_v1
        "hltEle32WP75GsfTrackIsoFilter", # HLT_Ele32_eta2p1_WP75_Gsf_v1
        #"hltEle27WP85GsfTrackIsoFilter", # HLT_Ele27_eta2p1_WP85_Gsf_v1
        #"hltEle32WP85GsfTrackIsoFilter", # HLT_Ele32_eta2p1_WP85_Gsf_v1
    ),
    probeFilters = cms.vstring(
        # control path HLT_DoubleEle24_22_eta2p1_WP75_Gsf_v1 
        # (for most of the following filters saveTags=False thus info not stored)
        "hltEle24Ele22leg2EtFilter",
        "hltEle24Ele22WP75leg2ClusterShapeFilter",
        "hltEle24Ele22WP75leg2HcEFilter",
        "hltEle24Ele22WP75leg2EcalIsoFilter",
        "hltEle24Ele22WP75leg2HcalIsoFilter",
        "hltEle24Ele22WP75Gsfleg2OneOESuperMinusOneOPFilter",
        "hltEle24Ele22WP75Gsfleg2Chi2Filter",
        "hltEle24Ele22WP75Gsfleg2DetaFilter",
        "hltEle24Ele22WP75Gsfleg2DphiFilter",
        "hltEle24Ele22WP75Gsfleg2TrackIsoFilter",
        # control path HLT_Ele22_eta2p1_WP75_Gsf_v1 (p'ed in data, for x-check with MC)
        # (for most of the following filters saveTags=False thus info not stored)
        "hltSingleEG22EtFilter",
        "hltSingleEle22WP75ClusterShapeFilter",
        "hltSingleEle22WP75HcEFilter",
        "hltSingleEle22WP75EcalIsoFilter",
        "hltSingleEle22WP75HcalIsoFilter",
        "hltSingleEle22WP75GsfOneOESuperMinusOneOPFilter",
        "hltSingleEle22WP75GsfChi2Filter",
        "hltSingleEle22WP75GsfDetaFilter",
        "hltSingleEle22WP75GsfDphiFilter",
        "hltSingleEle22WP75GsfTrackIsoFilter",
    ),
    l1NoIsoEGs = cms.InputTag("l1extraParticles:NonIsolated"), #FIXME: from RECO
    l1IsoEGs = cms.InputTag("l1extraParticles:Isolated"), #FIXME: from RECO
    #l1NoIsoEGs = cms.InputTag("hltL1extraParticles:NonIsolated"), #FIXME: from TauHLT (missing in v2)
    #l1IsoEGs = cms.InputTag("hltL1extraParticles:Isolated"), #FIXME: from TauHLT (missing in v2)
    minPtTag = cms.double(30.),
    maxEtaTag = cms.double(2.1),
    isoTag = cms.double(0.15),
    tagId = cms.string("NoImplemented"),
    minPtProbe = cms.double(20.),
    maxEtaProbe = cms.double(2.1),
    isoProbe = cms.double(0.1),
    probeId = cms.string("NoImplemented"),
    checkMCMatch = cms.bool(True)
)

import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
process.hltFilter = hlt.hltHighLevel.clone(
    TriggerResultsTag = process.eE.bits,
    HLTPaths = process.eE.tagTriggers,
    throw = False
)

process.preSelectedElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt > 20. && abs(eta) < 2.1")
)

process.countPreSelElectrons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("preSelectedElectrons")
)

process.TFileService = cms.Service(
    "TFileService", 
    fileName = cms.string("eETrgAna.root") 
)

process.p = cms.Path(
    process.hltFilter +
    process.preSelectedElectrons + process.countPreSelElectrons +
    process.eE
)



