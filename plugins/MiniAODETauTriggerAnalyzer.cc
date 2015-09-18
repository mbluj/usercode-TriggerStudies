// Based at example code from https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
// M. Bluj

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

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

class MiniAODETauTriggerAnalyzer : public edm::EDAnalyzer {
public:
  explicit MiniAODETauTriggerAnalyzer(const edm::ParameterSet&);
  ~MiniAODETauTriggerAnalyzer() {}

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  TTree *initTree(edm::Service<TFileService> &fs, std::string name="tree");
  void bookVariable(TTree *t=0, std::string var="foo");
  void cleanFilterVars();
  std::string triggerNameWithoutVersion(const std::string &triggerName);
  float electronIso(const pat::Electron &aEl, float dBetaFactor=0.5, bool allCharged=false);
  bool electronId(const pat::Electron &aEl, 
		  std::string idName="cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight");//FIXME
  bool isGoodElectron(const pat::Electron &aEl, const reco::Vertex & vtx,
		      std::string elId="cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight",
		      float dBetaFactor=0.5, bool allCharged=false);
  bool isGoodTau(const pat::Tau &aTau);

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

  edm::EDGetTokenT<pat::ElectronCollection> els_;
  edm::EDGetTokenT<pat::TauCollection> taus_;
  edm::EDGetTokenT<l1extra::L1JetParticleCollection> l1CenJets_, l1TauJets_, l1IsoTaus_;
  edm::EDGetTokenT<pat::METCollection> mets_;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;

  //std::vector<std::string> electronTriggers_;
  std::set<std::string> electronTriggers_;
  std::vector<std::string> electronFilters_;
  std::vector<std::string> tauFilters_;
  
  std::map<std::string, float> treeVars_;
  
  double minPtEl_, maxEtaEl_, isoEl_;
  double minPtTau_, maxEtaTau_;
  std::vector<std::string> tauIds_;

  std::string elId_;
  bool checkMCMatch_;

  TTree *tree_;
};

MiniAODETauTriggerAnalyzer::MiniAODETauTriggerAnalyzer(const edm::ParameterSet& iConfig):
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  els_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  taus_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  l1CenJets_(consumes<l1extra::L1JetParticleCollection>(iConfig.getParameter<edm::InputTag>("l1CenJets"))),
  l1TauJets_(consumes<l1extra::L1JetParticleCollection>(iConfig.getParameter<edm::InputTag>("l1TauJets"))),
  l1IsoTaus_(consumes<l1extra::L1JetParticleCollection>(iConfig.getParameter<edm::InputTag>("l1IsoTaus"))),
  mets_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("met"))),
  vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  electronFilters_(iConfig.getParameter<std::vector<std::string> >("electronFilters")),
  tauFilters_(iConfig.getParameter<std::vector<std::string> >("tauFilters")),
  minPtEl_(iConfig.getParameter<double>("minPtEl")),
  maxEtaEl_(iConfig.getParameter<double>("maxEtaEl")),
  isoEl_(iConfig.getParameter<double>("isoEl")),
  minPtTau_(iConfig.getParameter<double>("minPtTau")),
  maxEtaTau_(iConfig.getParameter<double>("maxEtaTau")),
  tauIds_(iConfig.getParameter<std::vector<std::string> >("tauIds")),
  elId_(iConfig.getParameter<std::string>("elId")),
  checkMCMatch_(iConfig.getParameter<bool>("checkMCMatch"))
{
  edm::Service<TFileService> fs;
  tree_ = initTree(fs,"eTauTriggerTree");
  bookVariable(tree_,"ePt");
  bookVariable(tree_,"eEta");
  bookVariable(tree_,"ePhi");
  bookVariable(tree_,"eM");
  bookVariable(tree_,"eCharge");
  bookVariable(tree_,"tauPt");
  bookVariable(tree_,"tauEta");
  bookVariable(tree_,"tauPhi");
  bookVariable(tree_,"tauM");
  bookVariable(tree_,"tauCharge");
  bookVariable(tree_,"eGenMatch");
  bookVariable(tree_,"tauGenMatch");
  bookVariable(tree_,"MEt");
  bookVariable(tree_,"MEtPhi");
  bookVariable(tree_,"Mass");
  bookVariable(tree_,"Mt");
  bookVariable(tree_,"Vx");
  bookVariable(tree_,"Vy");
  bookVariable(tree_,"Vz");
  
  std::vector<std::string> eTrgs( iConfig.getParameter<std::vector<std::string> >("electronTriggers") );
  for(unsigned int i=0; i<eTrgs.size(); ++i)
    electronTriggers_.insert( triggerNameWithoutVersion(eTrgs[i]) );
  for(std::set<std::string>::const_iterator it = electronTriggers_.begin(); it != electronTriggers_.end(); ++it)
    bookVariable(tree_,(*it) );
  for(unsigned int i=0; i<electronFilters_.size(); ++i)
    bookVariable(tree_,electronFilters_[i]);
  for(unsigned int i=0; i<tauFilters_.size(); ++i)
    bookVariable(tree_,tauFilters_[i]);
  bookVariable(tree_,"l1CenJet");
  bookVariable(tree_,"l1TauJet");
  bookVariable(tree_,"l1IsoTau");
}

TTree * MiniAODETauTriggerAnalyzer::initTree(edm::Service<TFileService> &fs, std::string name)
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

void MiniAODETauTriggerAnalyzer::bookVariable(TTree *t, std::string var)
{
  if(!t) return;
  treeVars_[var];
  t->Branch(var.c_str(), &treeVars_[var], (var + "/F").c_str());

  return;
}

void MiniAODETauTriggerAnalyzer::cleanFilterVars(){
  for(std::string &label : electronFilters_)
    treeVars_[label] = 0.;
  for(std::string &label : tauFilters_)
    treeVars_[label] = 0.;
  treeVars_["l1CenJet"] = 0.;
  treeVars_["l1TauJet"] = 0.;
  treeVars_["l1IsoTau"] = 0.;
}

std::string MiniAODETauTriggerAnalyzer::triggerNameWithoutVersion(const std::string &triggerName)
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
float MiniAODETauTriggerAnalyzer::electronIso(const pat::Electron &aEl, float dBetaFactor, bool allCharged){

  float iso = std::max((float)0.0,aEl.userIsolation("pfPhotons") + aEl.userIsolation("pfNeutralHadrons") 
		       - dBetaFactor*aEl.userIsolation("pfPUChargedHadrons") );
  if(allCharged)
    iso += aEl.userIsolation("pfChargedAll");
  else
    iso += aEl.userIsolation("pfChargedHadrons");

  return iso;
}

bool MiniAODETauTriggerAnalyzer::electronId(const pat::Electron &aEl, std::string idName){//FIXME

  return aEl.electronID(idName);
}

bool MiniAODETauTriggerAnalyzer::isGoodElectron(const pat::Electron &aEl, const reco::Vertex & vtx, std::string elId, float dBetaFactor, bool allCharged)
{
  //std::cout<<"Checking e: pt="<<aEl.pt()<<", eta="<<aEl.eta()<<", phi="<<aEl.phi()<<std::endl;
  // kinematics
  if(aEl.pt() < minPtEl_) return false;
  if(fabs( aEl.eta() ) > maxEtaEl_) return false;
  // e-Id
  if(elId != "")
    if( !electronId(aEl, elId) ) return false;  //FIXME
  // Vertex
  if( !(fabs( aEl.gsfTrack()->dxy( vtx.position() ) ) < 0.045) ) return false;
  if( !(fabs( aEl.gsfTrack()->dz( vtx.position() ) ) < 0.2) ) return false;
  // Iso
  if( electronIso(aEl,dBetaFactor,allCharged) > isoEl_ ) return false;

  return true;
}

bool MiniAODETauTriggerAnalyzer::isGoodTau(const pat::Tau &aTau)
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

void MiniAODETauTriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  edm::Handle<pat::ElectronCollection> els;
  iEvent.getByToken(els_, els);
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

  bool isTriggered = false;
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  //std::cout << "\n === TRIGGER PATHS === " << std::endl;
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
    if(electronTriggers_.find( triggerNameWithoutVersion( names.triggerName(i) ) ) != electronTriggers_.end() ) {
      /*
      std::cout << "Trigger " << names.triggerName(i) << 
	", prescale " << triggerPrescales->getPrescaleForIndex(i) <<
	": " << (triggerBits->accept(i) ? "PASS" : "fail (or not run)") 
		<< std::endl;
      */
      treeVars_[ triggerNameWithoutVersion( names.triggerName(i) )] = (triggerBits->accept(i) ? 1. : 0.)*triggerPrescales->getPrescaleForIndex(i);
      isTriggered |= triggerBits->accept(i);
    }
  }
  if(!isTriggered) return;

  //At least one vertex
  if(vtxs->size() == 0) return;
  treeVars_["Vx"] = vtxs->at(0).x();
  treeVars_["Vy"] = vtxs->at(0).y();
  treeVars_["Vz"] = vtxs->at(0).z();

  if(mets->size() > 0){
    treeVars_["MEt"] = mets->at(0).pt();
    treeVars_["MEtPhi"] = mets->at(0).phi();
  }


  for(const pat::Electron &el : *els) {
    if( !isGoodElectron(el,vtxs->at(0)) ) continue;
    //std::cout<<"Electron: pt="<<el.pt()<<", eta="<<el.eta()<<", phi="<<el.phi()<<std::endl;
    
    for(const pat::Tau &tau : *taus) {
      if( deltaR2(tau,el) < 0.5*0.5) continue;
      //FIXME OS??
      if( !isGoodTau(tau) ) continue;
      // Vertex
      if( fabs(tau.vertex().z() - vtxs->at(0).z())>0.0001 ) continue;
      // z impact at ECAL surface
      float zImpact = vtxs->at(0).z() + 130./tan(tau.theta());
      if( !(zImpact > 0.5 || zImpact < -1.5) ) continue; //FIXME: what is this? Why assymetric?
      //std::cout<<"Tau: pt="<<tau.pt()<<", eta="<<tau.eta()<<", phi="<<tau.phi()<<std::endl;
      //std::cout<<"DR(el,tau)="<<deltaR(tau,el)<<std::endl;
      treeVars_["ePt"] = el.pt();
      treeVars_["eEta"] = el.eta();
      treeVars_["ePhi"] = el.phi();
      treeVars_["eM"] = el.mass();
      treeVars_["eCharge"] = el.charge();
      treeVars_["tauPt"] = tau.pt();
      treeVars_["tauEta"] = tau.eta();
      treeVars_["tauPhi"] = tau.phi();
      treeVars_["tauM"] = tau.mass();
      treeVars_["tauCharge"] = tau.charge();
      treeVars_["Mass"] = (tau.p4()+el.p4()).mass();
      treeVars_["Mt"] = sqrt( 2.*el.pt()*mets->at(0).pt()*(1.-cos(deltaPhi(el.phi(),mets->at(0).phi()))));

      // gen match
      if(!checkMCMatch_ || el.genLepton() ) 
	treeVars_["eGenMatch"] = 1;
      else
	treeVars_["eGenMatch"] = 0;
      if(!checkMCMatch_ || tau.genJet() ) 
	treeVars_["tauGenMatch"] = 1;
      else
	treeVars_["tauGenMatch"] = 0;

      //std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
      cleanFilterVars();
      for(pat::TriggerObjectStandAlone trgObj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
	// check electron objects
	if( deltaR2(trgObj,el) < 0.5*0.5 && 
	    ( trgObj.hasTriggerObjectType(trigger::TriggerElectron) || trgObj.hasTriggerObjectType(trigger::TriggerPhoton) ) ) {
	  for(std::string &label : electronFilters_){
	    if( trgObj.hasFilterLabel(label) ) { 
	      treeVars_[label] = trgObj.pt();	    
	      /*
	      std::cout<<"\tFilter \""<<label<<"\" matched to electron, DR="<<deltaR(trgObj,el)<<std::endl
		       <<"\t\tTrigger object: pt="<<trgObj.pt()<<", eta="<<trgObj.eta()<<", phi="<<trgObj.phi() << std::endl; 
	      */
	    }
	  }
	}
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
      if(l1IsoTaus.isValid()){
	for(const l1extra::L1JetParticle &l1Tau : *l1IsoTaus) {
	  if(fabs(l1Tau.eta())<2.2 && l1Tau.et()>maxL1Et && deltaR2(l1Tau,tau) < 0.5*0.5){
	    maxL1Et = l1Tau.et();
	    treeVars_["l1IsoTau"] = l1Tau.et();
	  }
	}
      }
      //fill tree for each pair of good electron and good tau
      tree_->Fill();
    }
  }

  //std::cout << std::endl;
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODETauTriggerAnalyzer);
