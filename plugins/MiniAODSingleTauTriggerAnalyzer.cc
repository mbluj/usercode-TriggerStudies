//
// M. Bluj, National Centre for Nuclear Research, Poland
//

// system include files
#include <memory>
#include <cmath>

//
#include <vector>
#include <string>
#include <algorithm>
#include <set>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

class MiniAODSingleTauTriggerAnalyzer : public edm::EDAnalyzer {
public:
  explicit MiniAODSingleTauTriggerAnalyzer(const edm::ParameterSet&);
  ~MiniAODSingleTauTriggerAnalyzer() {}

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  TTree *initTree(edm::Service<TFileService> &fs, std::string name="tree");
  void bookVariable(TTree *t=0, std::string var="foo");
  void cleanFilterVars();
  std::string triggerNameWithoutVersion(const std::string &triggerName);

  bool isGoodTau(const pat::Tau &aTau);

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

  edm::EDGetTokenT<pat::TauCollection> taus_;
  edm::EDGetTokenT<l1extra::L1JetParticleCollection> l1CenJets_, l1TauJets_, l1IsoTaus_;
  edm::EDGetTokenT<pat::METCollection> mets_;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;
  edm::EDGetTokenT<GenEventInfoProduct> genEvtInfo_;

  std::vector<std::string> tauFilters_;
  
  std::map<std::string, float> treeVars_;
  
  double minPtTau_, maxEtaTau_;
  std::vector<std::string> tauIds_;
  std::vector<std::string> tauIdsForTrees_;

  bool checkMCMatch_;
  bool isMC_;

  TTree *tree_;
};

MiniAODSingleTauTriggerAnalyzer::MiniAODSingleTauTriggerAnalyzer(const edm::ParameterSet& iConfig):
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  taus_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  l1CenJets_(consumes<l1extra::L1JetParticleCollection>(iConfig.getParameter<edm::InputTag>("l1CenJets"))),
  l1TauJets_(consumes<l1extra::L1JetParticleCollection>(iConfig.getParameter<edm::InputTag>("l1TauJets"))),
  l1IsoTaus_(consumes<l1extra::L1JetParticleCollection>(iConfig.getParameter<edm::InputTag>("l1IsoTaus"))),
  mets_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("met"))),
  vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  genEvtInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEvtInfo"))),
  tauFilters_(iConfig.getParameter<std::vector<std::string> >("tauFilters")),
  minPtTau_(iConfig.getParameter<double>("minPtTau")),
  maxEtaTau_(iConfig.getParameter<double>("maxEtaTau")),
  tauIds_(iConfig.getParameter<std::vector<std::string> >("tauIds")),
  tauIdsForTrees_(iConfig.getParameter<std::vector<std::string> >("tauIdsForTrees")),
  checkMCMatch_(iConfig.getParameter<bool>("checkMCMatch")),
  isMC_(iConfig.getParameter<bool>("isMC"))
{
  edm::Service<TFileService> fs;
  tree_ = initTree(fs,"muTauTriggerTree");
  bookVariable(tree_,"tauPt");
  bookVariable(tree_,"tauEta");
  bookVariable(tree_,"tauPhi");
  bookVariable(tree_,"tauM");
  bookVariable(tree_,"tauCharge");
  bookVariable(tree_,"tauDM");
  bookVariable(tree_,"tauLeadChPt");
  bookVariable(tree_,"muGenMatch");
  bookVariable(tree_,"tauGenMatch");
  bookVariable(tree_,"MEt");
  bookVariable(tree_,"MEtPhi");
  bookVariable(tree_,"nVtx");
  bookVariable(tree_,"weight");

  for(unsigned int i=0; i<tauFilters_.size(); ++i)
    bookVariable(tree_,tauFilters_[i]);
  for(unsigned int i=0; i<tauIdsForTrees_.size(); ++i)
    bookVariable(tree_,tauIdsForTrees_[i]);
  bookVariable(tree_,"l1CenJet");
  bookVariable(tree_,"l1TauJet");
  bookVariable(tree_,"l1IsoTau");
}

TTree * MiniAODSingleTauTriggerAnalyzer::initTree(edm::Service<TFileService> &fs, std::string name)
{
  TTree *aTree = fs->make<TTree>( name.c_str(), name.c_str() );

  treeVars_["run"];
  aTree->Branch("run", &treeVars_["run"], "run/F");
  treeVars_["lumi"];
  aTree->Branch("lumi", &treeVars_["lumi"], "lumi/F");
  treeVars_["event"];
  aTree->Branch("event", &treeVars_["event"], "event/F");

  return aTree;
}

void MiniAODSingleTauTriggerAnalyzer::bookVariable(TTree *t, std::string var)
{
  if(!t) return;
  treeVars_[var];
  t->Branch(var.c_str(), &treeVars_[var], (var + "/F").c_str());

  return;
}

void MiniAODSingleTauTriggerAnalyzer::cleanFilterVars(){
  for(std::string &label : tauFilters_)    
    treeVars_[label] = 0.;

  treeVars_["l1CenJet"] = 0.;
  treeVars_["l1TauJet"] = 0.;
  treeVars_["l1IsoTau"] = 0.;
}

std::string MiniAODSingleTauTriggerAnalyzer::triggerNameWithoutVersion(const std::string &triggerName)
{
  std::string triggerNameWithoutVersion(triggerName);
  if( !triggerNameWithoutVersion.length() > 0 ) return triggerNameWithoutVersion; 
  if(triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] != '*') {
    while(triggerNameWithoutVersion.length() > 0
	  && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] >= '0'
	  && triggerNameWithoutVersion[triggerNameWithoutVersion.length()-1] <= '9') {
      triggerNameWithoutVersion.replace(triggerNameWithoutVersion.length()-1, 1, "");
    }
    //triggerNameWithoutVersion += "*";
  }
  else
    triggerNameWithoutVersion.replace(triggerNameWithoutVersion.length()-1, 1, "");
  return triggerNameWithoutVersion;
}

bool MiniAODSingleTauTriggerAnalyzer::isGoodTau(const pat::Tau &aTau)
{
  //std::cout<<"Checking tau: pt="<<aTau.pt()<<", eta="<<aTau.eta()<<", phi="<<aTau.phi()<<std::endl;
  // kinematics
  if(aTau.pt() < minPtTau_) return false;
  if(fabs( aTau.eta() ) > maxEtaTau_) return false;
  // Id
  for(unsigned int i=0; i<tauIds_.size(); ++i)
    if(aTau.tauID(tauIds_[i]) < 0.5) return false;

  return true;
}

void MiniAODSingleTauTriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(taus_, taus);
  edm::Handle<l1extra::L1JetParticleCollection> l1CenJets;
  iEvent.getByToken(l1CenJets_, l1CenJets);
  edm::Handle<l1extra::L1JetParticleCollection> l1TauJets;
  iEvent.getByToken(l1TauJets_, l1TauJets);
  edm::Handle<l1extra::L1JetParticleCollection> l1IsoTaus;
  iEvent.getByToken(l1IsoTaus_, l1IsoTaus);

  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(mets_, mets);

  edm::Handle<reco::VertexCollection> vtxs;
  iEvent.getByToken(vertices_, vtxs);

  treeVars_["run"]   = iEvent.id().run();
  treeVars_["lumi"]  = iEvent.id().luminosityBlock();
  treeVars_["event"] = iEvent.id().event();
  treeVars_["nVtx"]  = vtxs->size();
  if(isMC_){
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByToken(genEvtInfo_,genEvtInfo);
    //get event weight
    treeVars_["weight"] = genEvtInfo->weight();
  }
  else{
    treeVars_["weight"] = 1;
  }

  // bool isTriggered = false;
  // const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  // //std::cout << "\n === TRIGGER PATHS === " << std::endl;
  // for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
  //   if(muonTriggers_.find( triggerNameWithoutVersion( names.triggerName(i) ) ) != muonTriggers_.end() ) {
  //     /*
  //     std::cout << "Trigger " << names.triggerName(i) << 
  // 	", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
  // 	": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") 
  // 		<< std::endl;
  //     */
  //     //treeVars_[ triggerNameWithoutVersion( names.triggerName(i) )] = (triggerBits->accept(i) ? 1. : 0.)*triggerPrescales->getPrescaleForIndex(i);
  //     treeVars_[ triggerNameWithoutVersion( names.triggerName(i) )] = (triggerBits->accept(i) ? 1. : 0.);//FIXME: removed prescales for mismatched menu and prescale table (as in samples with rerun HLT) 
  //     isTriggered |= triggerBits->accept(i);
  //   }
  // }
  // if(!isTriggered) return;

  //At least one vertex
  if(vtxs->size() == 0) return;
  treeVars_["Vx"] = vtxs->at(0).x();
  treeVars_["Vy"] = vtxs->at(0).y();
  treeVars_["Vz"] = vtxs->at(0).z();

  if(mets->size() > 0){
    treeVars_["MEt"] = mets->at(0).pt();
    treeVars_["MEtPhi"] = mets->at(0).phi();
  }
    
  for(const pat::Tau &tau : *taus) {
    if( !isGoodTau(tau) ) continue;
    // Vertex
    if( fabs(tau.vertex().z() - vtxs->at(0).z())>0.0001 ) continue;
    pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
    if( fabs( packedLeadTauCand->dz(vtxs->at(0).position()) ) > 0.2) continue; 
    //std::cout<<"Tau: pt="<<tau.pt()<<", eta="<<tau.eta()<<", phi="<<tau.phi()<<std::endl;
    treeVars_["tauPt"] = tau.pt();
    treeVars_["tauEta"] = tau.eta();
    treeVars_["tauPhi"] = tau.phi();
    treeVars_["tauM"] = tau.mass();
    treeVars_["tauCharge"] = tau.charge();
    treeVars_["tauDM"] = tau.decayMode();
    treeVars_["tauLeadChPt"] = packedLeadTauCand->pt();
    for(std::string &label : tauIdsForTrees_){
      if(tau.isTauIDAvailable(label))
	treeVars_[label] = tau.tauID(label);
      else
	treeVars_[label] = -1.;
    }
    // gen match
    if(!checkMCMatch_ || tau.genJet() ) 
      treeVars_["tauGenMatch"] = 1;
    else
      treeVars_["tauGenMatch"] = 0;
      
    //std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
    cleanFilterVars();
    for(pat::TriggerObjectStandAlone trgObj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames

      // check tau objects
      if( deltaR2(trgObj,tau) < 0.5*0.5 && 
	  ( trgObj.hasTriggerObjectType(trigger::TriggerTau) || trgObj.hasTriggerObjectType(trigger::TriggerL1TauJet) ) ) {
	for(std::string &label : tauFilters_){
	  if( trgObj.hasFilterLabel(label) ) { 
	    treeVars_[label] = trgObj.pt();	    
	    /*
	      std::cout<<"\tFilter \""<<label<<"\" matched to tau, DR="<<deltaR(trgObj,tau)<<std::endl
	      <<"\t\tTrigger object: pt="<<trgObj.pt()<<", eta="<<trgObj.eta()<<", phi="<<trgObj.phi() << std::endl; 
	    */
	  }
	}
      }
    }
    //matching to L1Extra
    //CenJets
    float maxL1Et = 0;
    treeVars_["l1CenJet"] = 0;
    if(l1CenJets.isValid()){
      for(const l1extra::L1JetParticle &l1Tau : *l1CenJets) {
	if(fabs(l1Tau.eta())<2.2 && l1Tau.et()>maxL1Et && deltaR2(l1Tau,tau) < 0.5*0.5){
	  maxL1Et = l1Tau.et();
	  treeVars_["l1CenJet"] = l1Tau.et();
	}
      }
    }
    //TauJets
    maxL1Et = 0;
    treeVars_["l1TauJet"] = 0;
    if(l1TauJets.isValid()){
      for(const l1extra::L1JetParticle &l1Tau : *l1TauJets) {
	if(fabs(l1Tau.eta())<2.2 && l1Tau.et()>maxL1Et && deltaR2(l1Tau,tau) < 0.5*0.5){
	  maxL1Et = l1Tau.et();
	  treeVars_["l1TauJet"] = l1Tau.et();
	}
      }
    }
    //IsoTaus
    maxL1Et = 0;
    treeVars_["l1IsoTau"] = 0;
    if(l1IsoTaus.isValid()){
      for(const l1extra::L1JetParticle &l1Tau : *l1IsoTaus) {
	if(fabs(l1Tau.eta())<2.2 && l1Tau.et()>maxL1Et && deltaR2(l1Tau,tau) < 0.5*0.5){
	  maxL1Et = l1Tau.et();
	  treeVars_["l1IsoTau"] = l1Tau.et();
	}
      }
    }
    //fill tree for each good tau
    tree_->Fill();
  }

  //std::cout << std::endl;
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODSingleTauTriggerAnalyzer);
