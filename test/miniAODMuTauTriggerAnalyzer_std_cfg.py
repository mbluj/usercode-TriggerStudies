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
        'root://se.cis.gov.pl:1094//dpm/cis.gov.pl/home/cms/store/user/bluj/GluGluToHToTauTau_M-125_13TeV-powheg-pythia6/GGH125ToTauTau_MiniAOD_GRunV47_v1/8511c3cd7b71dd75a559c65d182442f5/miniAOD_100_1_Xv6.root'
    )
)

process.muLooseTau = cms.EDAnalyzer(
    "MiniAODMuTauTriggerAnalyzer",
    bits = cms.InputTag("TriggerResults","","HLT"), # use correct process name
    #bits = cms.InputTag("TriggerResults","","TauHLT"), # use correct process name
    prescales = cms.InputTag("patTrigger"), 
    objects = cms.InputTag("selectedPatTrigger"),
    #muons = cms.InputTag("slimmedMuons"),
    muons = cms.InputTag("preSelectedMuons"),
    #taus = cms.InputTag("slimmedTaus"),
    taus = cms.InputTag("preSelectedTaus"),
    met = cms.InputTag("slimmedMETs"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    muonTriggers = cms.vstring(  #version number is ignored, so can be replaced by wildcard (*) or dropped
        "HLT_IsoMu24_eta2p1_v*",
        "HLT_IsoMu24_eta2p1_IterTrk02_v*", # Phys'14
    ),
    muonFilters = cms.vstring(
        "hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09", # HLT_IsoMu24_eta2p1_v1
        "hltL3crIsoL1sMu20Eta2p1L1f0L2f20QL3f24QL3crIsoRhoFiltered0p15IterTrk02", # HLT_IsoMu24_eta2p1_IterTrk02_v1
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
    ),
    minPtMu = cms.double(25.),
    maxEtaMu = cms.double(2.1),
    isoMu = cms.double(0.15),
    minPtTau = cms.double(18.),
    maxEtaTau = cms.double(2.3),
    tauIds = cms.vstring(
        "decayModeFinding", 
        #"byLooseCombinedIsolationDeltaBetaCorr3Hits", #applied at presel level
        "againstMuonTight3", 
    ),
    checkMCMatch = cms.bool(True)
)

process.muMediumTau = process.muLooseTau.clone(
    taus = cms.InputTag("preSelectedMedTaus"),
    maxEtaTau = cms.double(2.1),
    tauIds = cms.vstring(
        "decayModeFinding",
        #"byMediumCombinedIsolationDeltaBetaCorr3Hits", #applied at presel level
        "againstMuonTight3",
    ),
    tauFilters = cms.vstring(
        "hltOverlapFilterIsoMu17MediumIsoPFTau40Reg", # HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v1
        "hltL1sMu16erIsoTau36er", # L1 of HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v1
        "hltL2Tau35eta2p2", # L2 of HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v1
        "hltL2IsoTau35eta2p2", # L2.5 of HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v1
        "hltOverlapFilterIsoMu17L2IsoTau35", # L2.5 of HLT_IsoMu17_eta2p1_MediumIsoPFTau40_Trk1_eta2p1_Reg_v1
    ),
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
    #cut = cms.string("pt > 25. && abs(eta) < 2.1 && isPFMuon && isGlobalMuon && userIsolation(\'pat::PfChargedAllIso\') < 1.5") # PfChargedAllIso=13= \'pfChargedAll\'which is not acceptd due to some reason
    #cut = cms.string("pt > 25. && abs(eta) < 2.1 && isPFMuon && isGlobalMuon")
    cut = cms.string("pt > 25. && abs(eta) < 2.1")
)

process.preSelectedTaus = cms.EDFilter(
    "PATTauSelector",
    src = cms.InputTag("slimmedTaus"),
    cut = cms.string("pt > 18. && abs(eta) < 2.3 && tauID(\'decayModeFinding\')> 0.5 && tauID(\'againstMuonLoose3\') > 0.5 && tauID(\'chargedIsoPtSum\') < 2.0" ) # Loose <2.0, Med <1.0, Tight<0.8
)

process.preSelectedMedTaus = process.preSelectedTaus.clone(
    cut = cms.string("pt > 18. && abs(eta) < 2.1 && tauID(\'decayModeFinding\')> 0.5 && tauID(\'againstMuonLoose3\') > 0.5 && tauID(\'chargedIsoPtSum\') < 1.0" ) # Loose <2.0, Med <1.0, Tight<0.8
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
process.countPreSelMedTaus = process.countPreSelTaus.clone(
    src = cms.InputTag("preSelectedMedTaus")
)

process.TFileService = cms.Service(
    "TFileService", 
    fileName = cms.string("muTauTrgAna.root") 
)

process.p = cms.Path(
    process.hltFilter +
    process.preSelectedMuons + process.countPreSelMuons +
    process.preSelectedTaus + process.countPreSelTaus +
    process.muLooseTau
)

'''
process.p2 = cms.Path(
    process.hltFilter +
    process.preSelectedMuons + process.countPreSelMuons +
    process.preSelectedMedTaus + process.countPreSelMedTaus +
    process.muMediumTau
)
'''
