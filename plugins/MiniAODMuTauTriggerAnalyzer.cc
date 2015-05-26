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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

class MiniAODMuTauTriggerAnalyzer : public edm::EDAnalyzer {
public:
  explicit MiniAODMuTauTriggerAnalyzer(const edm::ParameterSet&);
  ~MiniAODMuTauTriggerAnalyzer() {}

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  TTree *initTree(edm::Service<TFileService> &fs, std::string name="tree");
  void bookVariable(TTree *t=0, std::string var="foo");
  void cleanFilterVars();
  std::string triggerNameWithoutVersion(const std::string &triggerName);
  float muonIso(const pat::Muon &aMu, float dBetaFactor=0.5, bool allCharged=false);
  bool isGoodMuon(const pat::Muon &aMu, const reco::Vertex & vtx, bool useOldId=false, 
		  float dBetaFactor=0.5, bool allCharged=false);
  bool isGoodTau(const pat::Tau &aTau);

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

  edm::EDGetTokenT<pat::MuonCollection> muons_;
  edm::EDGetTokenT<pat::TauCollection> taus_;
  edm::EDGetTokenT<pat::METCollection> mets_;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;

  //std::vector<std::string> muonTriggers_;
  std::set<std::string> muonTriggers_;
  std::vector<std::string> muonFilters_;
  std::vector<std::string> tauFilters_;
  
  std::map<std::string, float> treeVars_;
  
  double minPtMu_, maxEtaMu_, isoMu_;
  double minPtTau_, maxEtaTau_;
  std::vector<std::string> tauIds_;

  bool muTightId_;
  bool checkMCMatch_;

  TTree *tree_;
};

MiniAODMuTauTriggerAnalyzer::MiniAODMuTauTriggerAnalyzer(const edm::ParameterSet& iConfig):
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  muons_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  taus_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  mets_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("met"))),
  vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  muonFilters_(iConfig.getParameter<std::vector<std::string> >("muonFilters")),
  tauFilters_(iConfig.getParameter<std::vector<std::string> >("tauFilters")),
  minPtMu_(iConfig.getParameter<double>("minPtMu")),
  maxEtaMu_(iConfig.getParameter<double>("maxEtaMu")),
  isoMu_(iConfig.getParameter<double>("isoMu")),
  minPtTau_(iConfig.getParameter<double>("minPtTau")),
  maxEtaTau_(iConfig.getParameter<double>("maxEtaTau")),
  tauIds_(iConfig.getParameter<std::vector<std::string> >("tauIds")),
  muTightId_(iConfig.getParameter<bool>("muTightId")),
  checkMCMatch_(iConfig.getParameter<bool>("checkMCMatch"))
{
  edm::Service<TFileService> fs;
  tree_ = initTree(fs,"muTauTriggerTree");
  bookVariable(tree_,"muPt");
  bookVariable(tree_,"muEta");
  bookVariable(tree_,"muPhi");
  bookVariable(tree_,"muM");
  bookVariable(tree_,"muCharge");
  bookVariable(tree_,"tauPt");
  bookVariable(tree_,"tauEta");
  bookVariable(tree_,"tauPhi");
  bookVariable(tree_,"tauM");
  bookVariable(tree_,"tauCharge");
  bookVariable(tree_,"muGenMatch");
  bookVariable(tree_,"tauGenMatch");
  bookVariable(tree_,"MEt");
  bookVariable(tree_,"MEtPhi");
  bookVariable(tree_,"Mass");
  bookVariable(tree_,"Mt");
  bookVariable(tree_,"Vx");
  bookVariable(tree_,"Vy");
  bookVariable(tree_,"Vz");


  std::vector<std::string> muTrgs( iConfig.getParameter<std::vector<std::string> >("muonTriggers") );
  for(unsigned int i=0; i<muTrgs.size(); ++i)
    muonTriggers_.insert( triggerNameWithoutVersion(muTrgs[i]) );
  for(std::set<std::string>::const_iterator it = muonTriggers_.begin(); it != muonTriggers_.end(); ++it)
    bookVariable(tree_,(*it) );
  for(unsigned int i=0; i<muonFilters_.size(); ++i)
    bookVariable(tree_,muonFilters_[i]);
  for(unsigned int i=0; i<tauFilters_.size(); ++i)
    bookVariable(tree_,tauFilters_[i]);
}

TTree * MiniAODMuTauTriggerAnalyzer::initTree(edm::Service<TFileService> &fs, std::string name)
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

void MiniAODMuTauTriggerAnalyzer::bookVariable(TTree *t, std::string var)
{
  if(!t) return;
  treeVars_[var];
  t->Branch(var.c_str(), &treeVars_[var], (var + "/F").c_str());

  return;
}

void MiniAODMuTauTriggerAnalyzer::cleanFilterVars(){
  for(std::string &label : muonFilters_)
    treeVars_[label] = 0.;
  for(std::string &label : tauFilters_)
    treeVars_[label] = 0.;
}

std::string MiniAODMuTauTriggerAnalyzer::triggerNameWithoutVersion(const std::string &triggerName)
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

float MiniAODMuTauTriggerAnalyzer::muonIso(const pat::Muon &aMu, float dBetaFactor, bool allCharged)
{
  float iso = std::max((float)0.0,aMu.userIsolation("pfPhotons") + aMu.userIsolation("pfNeutralHadrons") - dBetaFactor*aMu.userIsolation("pfPUChargedHadrons") );
  if(allCharged)
    iso += aMu.userIsolation("pfChargedAll");
  else
    iso += aMu.userIsolation("pfChargedHadrons");

  return iso;
}
      

bool MiniAODMuTauTriggerAnalyzer::isGoodMuon(const pat::Muon &aMu, const reco::Vertex & vtx, bool useOldId, float dBetaFactor, bool allCharged)
{
  //std::cout<<"Checking mu: pt="<<aMu.pt()<<", eta="<<aMu.eta()<<", phi="<<aMu.phi()<<std::endl;
  // kinematics
  if(aMu.pt() < minPtMu_) return false;
  if(fabs( aMu.eta() ) > maxEtaMu_) return false;
  if(useOldId){
    // Tight mu
    if( !aMu.isGlobalMuon() || !aMu.isPFMuon() ) return false;
    if( !(aMu.normChi2() < 10.) ) return false;
    if( !(aMu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0) ) return false;
    if( !(aMu.numberOfMatchedStations() > 1) ) return false;
    if( !(fabs( aMu.muonBestTrack()->dxy( vtx.position() ) ) < 0.2) ) return false;
    if( !(fabs( aMu.muonBestTrack()->dz( vtx.position() ) ) < 0.5) ) return false;
    if( !(aMu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0) ) return false;
    if( !(aMu.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5) ) return false;
  }
  else {
    // Medium mu
    if( !aMu.isLooseMuon() ) return false;
    bool isGlb = ( aMu.isGlobalMuon() &&
		   aMu.normChi2() < 3 &&
		   aMu.combinedQuality().chi2LocalPosition < 12 && //FIXME: accessible in miniAOD??
		   aMu.combinedQuality().trkKink < 20 ); //FIXME: accessible in miniAOD??
    bool isGood = ( aMu.innerTrack()->validFraction() >= 0.8 &&
		    aMu.segmentCompatibility() >=  (isGlb ? 0.303 : 0.451) );
    if( !isGood) return false;
  }
  // Vertex
  if( !(fabs( aMu.muonBestTrack()->dxy( vtx.position() ) ) < 0.045) ) return false;
  if( !(fabs( aMu.muonBestTrack()->dz( vtx.position() ) ) < 0.2) ) return false;
  // Iso
  if( muonIso(aMu, dBetaFactor=0.5, allCharged) > isoMu_ ) return false;

  return true;
}

bool MiniAODMuTauTriggerAnalyzer::isGoodTau(const pat::Tau &aTau)
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

void MiniAODMuTauTriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muons_, muons);
  edm::Handle<pat::TauCollection> taus;
  iEvent.getByToken(taus_, taus);
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
    if(muonTriggers_.find( triggerNameWithoutVersion( names.triggerName(i) ) ) != muonTriggers_.end() ) {
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

  for(const pat::Muon &mu : *muons) {
    if( !isGoodMuon(mu,vtxs->at(0),muTightId_) ) continue;
    //std::cout<<"Muon: pt="<<mu.pt()<<", eta="<<mu.eta()<<", phi="<<mu.phi()<<std::endl;
    
    for(const pat::Tau &tau : *taus) {
      if( deltaR2(tau,mu) < 0.5*0.5) continue;
      //FIXME OS??
      if( !isGoodTau(tau) ) continue;
      // Vertex
      if( fabs(tau.vertex().z() - vtxs->at(0).z())>0.0001 ) continue;
      //std::cout<<"Tau: pt="<<tau.pt()<<", eta="<<tau.eta()<<", phi="<<tau.phi()<<std::endl;
      //std::cout<<"DR(mu,tau)="<<deltaR(tau,mu)<<std::endl;
      treeVars_["muPt"] = mu.pt();
      treeVars_["muEta"] = mu.eta();
      treeVars_["muPhi"] = mu.phi();
      treeVars_["muM"] = mu.mass();
      treeVars_["muCharge"] = mu.charge();
      treeVars_["tauPt"] = tau.pt();
      treeVars_["tauEta"] = tau.eta();
      treeVars_["tauPhi"] = tau.phi();
      treeVars_["tauM"] = tau.mass();
      treeVars_["tauCharge"] = tau.charge();
      treeVars_["Mass"] = (tau.p4()+mu.p4()).mass();
      treeVars_["Mt"] = fabs( 2.*mu.pt()*mets->at(0).pt()*(1.-cos(deltaPhi(mu.phi(),mets->at(0).phi()))));
      // gen match
      if(!checkMCMatch_ || mu.genLepton() ) 
	treeVars_["muGenMatch"] = 1;
      else
	treeVars_["muGenMatch"] = 0;
      if(!checkMCMatch_ || tau.genJet() ) 
	treeVars_["tauGenMatch"] = 1;
      else
	treeVars_["tauGenMatch"] = 0;

      //std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
      cleanFilterVars();
      for(pat::TriggerObjectStandAlone trgObj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
	// check muon objects
	if( deltaR2(trgObj,mu) < 0.3*0.3 && trgObj.hasTriggerObjectType(trigger::TriggerMuon) ) {
	  for(std::string &label : muonFilters_){
	    if( trgObj.hasFilterLabel(label) ) { 
	      treeVars_[label] = trgObj.pt();	    
	      /*
	      std::cout<<"\tFilter \""<<label<<"\" matched to muon, DR="<<deltaR(trgObj,mu)<<std::endl
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
      //fill tree for each pair of good muon and good tau
      tree_->Fill();
    }
  }

  //std::cout << std::endl;
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODMuTauTriggerAnalyzer);
