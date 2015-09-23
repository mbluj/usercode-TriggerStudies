# example from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD

import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

isMC = False

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
        #"root://cms-xrd-global.cern.ch///store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/009D49A5-7314-E511-84EF-0025905A605E.root"
        "root://xrootd-cms.infn.it//store/data/Run2015C/SingleMuon/MINIAOD/PromptReco-v1/000/254/833/00000/FCB9CBAE-104B-E511-BF89-02163E011F65.root"
    )
)

process.evtCounter = cms.EDAnalyzer(
    "BasicGenEventInfoAnalyzer",
    isMC = cms.bool(isMC),
    genEvtInfo = cms.InputTag("generator")
)

process.muMu = cms.EDAnalyzer(
    "MiniAODMuMuTriggerAnalyzer",
    bits = cms.InputTag("TriggerResults","","HLT"), # use correct process name
    prescales = cms.InputTag("patTrigger"), 
    objects = cms.InputTag("selectedPatTrigger"),
    #muons = cms.InputTag("slimmedMuons"),
    muons = cms.InputTag("preSelectedMuons"),
    l1Muons = cms.InputTag("l1extraParticles"),
    met = cms.InputTag("slimmedMETs"),
    vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
    tagTriggers = cms.vstring(  #version number is ignored, so can be replaced by wildcard (*) or dropped
        "HLT_IsoMu24_eta2p1_v*", # unp'ed at 1.4e34 
        "HLT_IsoMu20_eta2p1_v*", # unp'ed at 7e33 
        "HLT_IsoMu17_eta2p1_v*", # unp'ed at 5e33
    ),
    tagFilters = cms.vstring(
        "hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09", # HLT_IsoMu24_eta2p1_v1
        "hltL3crIsoL1sMu16Eta2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p09", # HLT_IsoMu20_eta2p1_v1
        "hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09", # HLT_IsoMu17_eta2p1_v1
    ),
    probeFilters = cms.vstring(
        # control path HLT_DoubleIsoMu17_eta2p1_v1
        "hltL1sL1DoubleMu125", # L1Mu Pt>12,5
        "hltL1fL1sDoubleMu125L1Filtered16er", # L1Mu with Pt>16,16, |eta|<2.1
        "hltL2fL1sDoubleMu125L1f16erL2Filtered10Q", # L2Mu with Pt>10
        "hltL3fL1sDoubleMu125L1f16erL2f10QL3Filtered17Q", # L3Mu with Pt>17
        "hltL3DzL1sDoubleMu125L1f16erL2f10QL3f17QL3DzFiltered0p2", # Dz(mu,mu)<0.2 (sanity check)
        "hltL3crIsoL1sDoubleMu125L1f16erL2f10QL3f17QL3Dz0p2L3crIsoRhoFiltered0p15IterTrk02", # full isolation
        # control path HLT_IsoMu17_eta2p1_v1 (p'ed in data, for x-check with MC)
        "hltL1sSingleMu16er", #L1 seed
        "hltL1fL1sSingleMu16erL1Filtered0", # L1Mu with Pt>16
        "hltL2fL1sSingleMu16erL1f0L2Filtered10Q", # L2Mu with Pt>10 
        "hltL3fL1sSingleMu16erL1f0L2f10QL3Filtered17Q", # L3Mu with Pt>17
        "hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09", # full isolation
    ),
    minPtTag = cms.double(21.),
    maxEtaTag = cms.double(2.1),
    isoTag = cms.double(0.15),
    minPtProbe = cms.double(15.),
    maxEtaProbe = cms.double(2.1),
    isoProbe = cms.double(0.15),
    tagTightId = cms.bool(False),
    probeTightId = cms.bool(False),
    checkMCMatch = cms.bool(isMC),
    isMC = cms.bool(isMC),
    genEvtInfo = cms.InputTag("generator"),
)

import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt
process.hltFilter = hlt.hltHighLevel.clone(
    TriggerResultsTag = process.muMu.bits,
    HLTPaths = process.muMu.tagTriggers,
    throw = False
)

process.preSelectedMuons = cms.EDFilter(
    "PATMuonSelector",
    src = cms.InputTag("slimmedMuons"),
    #cut = cms.string("pt > 15. && abs(eta) < 2.1 && isPFMuon && isGlobalMuon && userIsolation(\'pat::PfChargedAllIso\') < 1.5") # PfChargedAllIso=13= \'pfChargedAll\'which is not acceptd due to some reason
    #cut = cms.string("pt > 15. && abs(eta) < 2.1 && isPFMuon && isGlobalMuon")
    cut = cms.string("pt > 15. && abs(eta) < 2.1")
)

process.countPreSelMuons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("preSelectedMuons")
)

process.TFileService = cms.Service(
    "TFileService", 
    fileName = cms.string("muMuTrgAna.root") 
)

process.p0 = cms.Path(
    process.evtCounter
)

process.p = cms.Path(
    process.hltFilter +
    process.preSelectedMuons + process.countPreSelMuons +
    process.muMu
)



