# Configuration by Michal Bluj based on auto generated configuration file
# cff
# Revision 1.0.pre0
# 29/01/2015

import FWCore.ParameterSet.Config as cms

###############################
def addMiniAODOutput(process, outFile='miniAOD-prod_PAT.root'):

    #####################
    # Output definition
    process.load('Configuration.EventContent.EventContent_cff')
    process.MINIAODSIMoutput = cms.OutputModule(
        "PoolOutputModule",
        compressionLevel = cms.untracked.int32(4),
        compressionAlgorithm = cms.untracked.string('LZMA'),
        eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
        outputCommands = process.MINIAODSIMEventContent.outputCommands,
        fileName = cms.untracked.string( outFile ),
        dataset = cms.untracked.PSet(
            filterName = cms.untracked.string(''),
            dataTier = cms.untracked.string('')
        ),
        dropMetaData = cms.untracked.string('ALL'),
        fastCloning = cms.untracked.bool(False),
        overrideInputFileSplitLevels = cms.untracked.bool(True)
    )
    # Additional output definition

    #####################
    # Path and EndPath definitions
    process.load('Configuration.StandardSequences.EndOfProcess_cff')
    process.endjob_step = cms.EndPath(process.endOfProcess)
    process.MINIAODSIMoutput_step = cms.EndPath(process.MINIAODSIMoutput)


def addMiniAODSched(process):
    #####################
    # PAT and MiniAOD
    process.load('Configuration.StandardSequences.PATMC_cff')

    # Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
    from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllMC 

    #call to customisation function miniAOD_customizeAllMC imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
    process = miniAOD_customizeAllMC(process)

    # Define scheduled PAT sequence
    #process.load('PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff') #MB

    process.prePatSequenceSched = cms.Sequence(
        process.patPFMetOriginalReserved +
        # reduced e/gamma information
        process.reducedEgamma +
        # pruned genParticles
        process.prunedGenParticlesWithStatusOne + process.prunedGenParticles +
        # e-id's
        process.egmGsfElectronIDSequence +
        # pfCandidates filtering
        process.particleFlowPtrs + 
        process.pfNoPileUpSequence +
        process.pfNoPileUpJMESequence +
        # top taggers
        process.cmsttRaw + process.caTopTagInfos +
        # groomed/trimmed/filtered/jets
        process.ak8PFJetsCHSPruned + process.ak8PFJetsCHSTrimmed+ process.ak8PFJetsCHSFiltered +
        process.ak8PFJetsCHSPrunedLinks+ process.ak8PFJetsCHSTrimmedLinks + process.ak8PFJetsCHSFilteredLinks +
        # Njetiness
        process.NjettinessAK8 +
        # puJetId
        process.pileupJetId +
        # b-tag
        process.impactParameterTagInfos + process.secondaryVertexTagInfos +
        # prePAT slimming
        process.offlineSlimmedPrimaryVertices +
        process.packedGenParticles +
        process.packedPFCandidates +
        process.lostTracks +
        process.slimmedGenJets
    )    

    #add AK8 jets to PAT seqneces
    process.makeJetsAK8 = cms.Sequence(
        process.jetTracksAssociatorAtVertexAK8 +
        process.patJetCorrFactorsAK8 +
        process.patJetChargeAK8 + 
        process.patJetPartonMatchAK8 + process.patJetGenJetMatchAK8 +
        process.patJetPartonsLegacy + process.patJetPartonAssociationLegacyAK8 + process.patJetFlavourAssociationLegacyAK8 +
        process.patJetPartons + process.patJetFlavourAssociationAK8 +
        process.patJetsAK8
    )
    process.patCandidates.replace(
        process.makePatJets,
        process.makePatJets + process.makeJetsAK8
    )
    process.selectedPatCandidates.replace(
        process.selectedPatJets,
        process.selectedPatJets + process.selectedPatJetsAK8
   )

   #add jets for MEt unc to PAT seqneces
    process.makeJetsAK4PFForMetUnc = cms.Sequence(
        process.patJetCorrFactorsAK4PFForMetUnc +
        process.patJetCharge + 
        process.patJetPartonMatchAK4PFForMetUnc + process.patJetGenJetMatchAK4PFForMetUnc +
        process.patJetPartonsLegacy + process.patJetPartonAssociationLegacyAK4PFForMetUnc + process.patJetFlavourAssociationLegacyAK4PFForMetUnc +
        process.patJetPartonsForMetUnc + process.patJetFlavourAssociationAK4PFForMetUnc +
        process.patJetsAK4PFForMetUnc
    )
    process.patCandidates.replace(
        process.makePatJets,
        process.makePatJets + process.makeJetsAK4PFForMetUnc
    )
    process.selectedPatCandidates.replace(
        process.selectedPatJets,
        process.selectedPatJets + process.selectedPatJetsAK4PFForMetUnc
    )

    # Scheduled PAT sequence
    process.patSequenceSched = cms.Sequence(
        process.patCandidates +
        process.selectedPatCandidates +
        # met uncertainty
        process.producePatPFMETCorrectionsUncOriginalReserved +
        process.pfType1MEtUncertaintySequence +
        # pat trigger
        process.patTrigger #+ process.patTriggerEvent
    )
    # Scheduled slimming sequence for MiniAOD
    process.miniAODSlimmingSequenceSched = cms.Sequence(
        process.selectedPatTrigger +  
        process.slimmedJets + process.slimmedJetsAK8 +
        process.slimmedElectrons +
        process.slimmedMuons + 
        process.slimmedPhotons +
        process.slimmedTaus+ 
        process.slimmedSecondaryVertices +
        process.slimmedMETs
    )
    # Scheduled PAT + MiniAOD slimming path
    process.patAndMiniAODPath = cms.Path(
        process.prePatSequenceSched +
        process.patSequenceSched +
        process.miniAODSlimmingSequenceSched
    )

# End of customisation functions

