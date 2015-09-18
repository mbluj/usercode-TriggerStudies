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
  bool isGoodMuon(const pat::Muon &aMu, const reco::Vertex & vtx, std::string muId="POG_Medium", 
		  float dBetaFactor=0.5, bool allCharged=false);
  bool isGoodTau(const pat::Tau &aTau);

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

  edm::EDGetTokenT<pat::MuonCollection> muons_;
  edm::EDGetTokenT<pat::TauCollection> taus_;
  edm::EDGetTokenT<l1extra::L1JetParticleCollection> l1CenJets_, l1TauJets_, l1IsoTaus_;
  edm::EDGetTokenT<pat::METCollection> mets_;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;
  edm::EDGetTokenT<GenEventInfoProduct> genEvtInfo_;

  //std::vector<std::string> muonTriggers_;
  std::set<std::string> muonTriggers_;
  std::vector<std::string> muonFilters_;
  std::vector<std::string> tauFilters_;
  
  std::map<std::string, float> treeVars_;
  
  double minPtMu_, maxEtaMu_, isoMu_;
  double minPtTau_, maxEtaTau_;
  std::vector<std::string> tauIds_;
  std::vector<std::string> tauIdsForTrees_;

  bool muTightId_;
  bool checkMCMatch_;
  bool isMC_;

  TTree *tree_;
};

MiniAODMuTauTriggerAnalyzer::MiniAODMuTauTriggerAnalyzer(const edm::ParameterSet& iConfig):
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  muons_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  taus_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
  l1CenJets_(consumes<l1extra::L1JetParticleCollection>(iConfig.getParameter<edm::InputTag>("l1CenJets"))),
  l1TauJets_(consumes<l1extra::L1JetParticleCollection>(iConfig.getParameter<edm::InputTag>("l1TauJets"))),
  l1IsoTaus_(consumes<l1extra::L1JetParticleCollection>(iConfig.getParameter<edm::InputTag>("l1IsoTaus"))),
  mets_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("met"))),
  vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  genEvtInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEvtInfo"))),
  muonFilters_(iConfig.getParameter<std::vector<std::string> >("muonFilters")),
  tauFilters_(iConfig.getParameter<std::vector<std::string> >("tauFilters")),
  minPtMu_(iConfig.getParameter<double>("minPtMu")),
  maxEtaMu_(iConfig.getParameter<double>("maxEtaMu")),
  isoMu_(iConfig.getParameter<double>("isoMu")),
  minPtTau_(iConfig.getParameter<double>("minPtTau")),
  maxEtaTau_(iConfig.getParameter<double>("maxEtaTau")),
  tauIds_(iConfig.getParameter<std::vector<std::string> >("tauIds")),
  tauIdsForTrees_(iConfig.getParameter<std::vector<std::string> >("tauIdsForTrees")),
  muTightId_(iConfig.getParameter<bool>("muTightId")),
  checkMCMatch_(iConfig.getParameter<bool>("checkMCMatch")),
  isMC_(iConfig.getParameter<bool>("isMC"))
{
  edm::Service<TFileService> fs;
  tree_ = initTree(fs,"muTauTriggerTree");
  bookVariable(tree_,"muPt");
  bookVariable(tree_,"muEta");
  bookVariable(tree_,"muPhi");
  bookVariable(tree_,"muM");
  bookVariable(tree_,"muCharge");
  bookVariable(tree_,"muIso");
  bookVariable(tree_,"mu2Pt");
  bookVariable(tree_,"mu2Eta");
  bookVariable(tree_,"mu2Phi");
  bookVariable(tree_,"mu2M");
  bookVariable(tree_,"mu2Charge");
  bookVariable(tree_,"mu2Iso");
  bookVariable(tree_,"mu2Id");
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
  bookVariable(tree_,"Mass");
  bookVariable(tree_,"Mt");
  bookVariable(tree_,"nVtx");
  bookVariable(tree_,"weight");
  bookVariable(tree_,"iPair");


  std::vector<std::string> muTrgs( iConfig.getParameter<std::vector<std::string> >("muonTriggers") );
  for(unsigned int i=0; i<muTrgs.size(); ++i)
    muonTriggers_.insert( triggerNameWithoutVersion(muTrgs[i]) );
  for(std::set<std::string>::const_iterator it = muonTriggers_.begin(); it != muonTriggers_.end(); ++it)
    bookVariable(tree_,(*it) );
  for(unsigned int i=0; i<muonFilters_.size(); ++i)
    bookVariable(tree_,muonFilters_[i]);
  for(unsigned int i=0; i<tauFilters_.size(); ++i)
    bookVariable(tree_,tauFilters_[i]);
  for(unsigned int i=0; i<tauIdsForTrees_.size(); ++i)
    bookVariable(tree_,tauIdsForTrees_[i]);
  bookVariable(tree_,"l1CenJet");
  bookVariable(tree_,"l1TauJet");
  bookVariable(tree_,"l1IsoTau");
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

  treeVars_["l1CenJet"] = 0.;
  treeVars_["l1TauJet"] = 0.;
  treeVars_["l1IsoTau"] = 0.;
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
      

bool MiniAODMuTauTriggerAnalyzer::isGoodMuon(const pat::Muon &aMu, const reco::Vertex & vtx, std::string muId, float dBetaFactor, bool allCharged)
{
  //See https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
  if(muId=="POG_Tight"||muId=="Tight"){
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
  else if(muId=="POG_Medium"||muId=="Medium"){
    // Medium mu
    if( !aMu.isLooseMuon() ) return false;
    bool isGlb = ( aMu.isGlobalMuon() &&
		   aMu.normChi2() < 3 &&
		   aMu.combinedQuality().chi2LocalPosition < 12 && //FIXME: accessible in miniAOD??
		   aMu.combinedQuality().trkKink < 20 ); //FIXME: accessible in miniAOD??
    bool isGood = ( aMu.innerTrack()->validFraction() > 0.8 &&
		    aMu.segmentCompatibility() >  (isGlb ? 0.303 : 0.451) );
    if( !isGood) return false;
  }
  else {
    // Basic mu-Id (==loose)
    if( !aMu.isLooseMuon() ) return false; //i.e. isGlobalMuon&&isTrackerMuon&&isPFMuon
  }
  // Vertex
  if( !(fabs( aMu.muonBestTrack()->dxy( vtx.position() ) ) < 0.045) ) return false;
  if( !(fabs( aMu.muonBestTrack()->dz( vtx.position() ) ) < 0.2) ) return false;
  // Iso
  //if( muonIso(aMu, dBetaFactor, allCharged) > isoMu_ ) return false;

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
      //treeVars_[ triggerNameWithoutVersion( names.triggerName(i) )] = (triggerBits->accept(i) ? 1. : 0.)*triggerPrescales->getPrescaleForIndex(i);
      treeVars_[ triggerNameWithoutVersion( names.triggerName(i) )] = (triggerBits->accept(i) ? 1. : 0.);//FIXME: removed prescales for mismatched menu and prescale table (as in samples with rerun HLT) 
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
  int iPair=0;
  //for(const pat::Muon &mu : *muons) {
  for(pat::MuonCollection::const_iterator iMu=muons->begin(); 
      iMu!=muons->end(); ++iMu){
    const pat::Muon &mu = (*iMu);
    // kinematics
    if(mu.pt() < minPtMu_) continue;
    if(fabs( mu.eta() ) > maxEtaMu_) continue;
    std::string muWP = muTightId_ ? "POG_Tight" : "POG_Medium";
    if( !isGoodMuon(mu,vtxs->at(0),muWP) ) continue;
    float muIso = muonIso(mu);
    if(muonIso(mu) > isoMu_) continue;
    //std::cout<<"Muon: pt="<<mu.pt()<<", eta="<<mu.eta()<<", phi="<<mu.phi()<<std::endl;
    
    for(const pat::Tau &tau : *taus) {
      if( deltaR2(tau,mu) < 0.5*0.5) continue;
      if( !isGoodTau(tau) ) continue;
      // Vertex
      if( fabs(tau.vertex().z() - vtxs->at(0).z())>0.0001 ) continue;
      pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(tau.leadChargedHadrCand().get());
      if( fabs( packedLeadTauCand->dz(vtxs->at(0).position()) ) > 0.2) continue; 
      //std::cout<<"Tau: pt="<<tau.pt()<<", eta="<<tau.eta()<<", phi="<<tau.phi()<<std::endl;
      //std::cout<<"DR(mu,tau)="<<deltaR(tau,mu)<<std::endl;
      treeVars_["muPt"] = mu.pt();
      treeVars_["muEta"] = mu.eta();
      treeVars_["muPhi"] = mu.phi();
      treeVars_["muM"] = mu.mass();
      treeVars_["muCharge"] = mu.charge();
      treeVars_["muIso"] = muIso;
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
      treeVars_["Mass"] = (tau.p4()+mu.p4()).mass();
      treeVars_["Mt"] = sqrt( 2.*mu.pt()*mets->at(0).pt()*(1.-cos(deltaPhi(mu.phi(),mets->at(0).phi()))));
      // gen match
      if(!checkMCMatch_ || mu.genLepton() ) 
	treeVars_["muGenMatch"] = 1;
      else
	treeVars_["muGenMatch"] = 0;
      if(!checkMCMatch_ || tau.genJet() ) 
	treeVars_["tauGenMatch"] = 1;
      else
	treeVars_["tauGenMatch"] = 0;
      //check presence of 2nd muon
      treeVars_["mu2Pt"] = 0;
      treeVars_["mu2Eta"] = 0;
      treeVars_["mu2Phi"] = 0;
      treeVars_["mu2M"] = 0;
      treeVars_["mu2Charge"] = 0;
      treeVars_["mu2Iso"] = 999;
      treeVars_["mu2Id"] = -1;
      for(pat::MuonCollection::const_iterator iMu2=iMu; 
	  iMu2!=muons->end(); ++iMu2){
	const pat::Muon &mu2 = (*iMu2);
	if( deltaR2(mu2,mu) < 0.15*0.15) continue;
	if(mu2.pt() < 10) continue;
	if(fabs( mu2.eta() ) > 2.4) continue;
	if( !isGoodMuon(mu2,vtxs->at(0),"") ) continue;
	int mu2Id = 0;
	if( isGoodMuon(mu2,vtxs->at(0),"POG_Tight") )
	  mu2Id += 10;
	if( isGoodMuon(mu2,vtxs->at(0),"POG_Medium") )
	  mu2Id += 1;
	float mu2Iso = muonIso(mu2);
	if(mu2.pt()>treeVars_["mu2Pt"] && mu2Iso<treeVars_["mu2Iso"]){
	  if( (mu2.pt()>15 && mu.charge()*mu2.charge()<0 ) || //as in 2nd mu veto
	      (mu2.pt()>10 && mu2Id%10==1) || //as in 2nd mu veto
	      treeVars_["mu2Id"] < 0 ){ //any additional muon
	    treeVars_["mu2Pt"] = mu2.pt();
	    treeVars_["mu2Eta"] = mu2.eta();
	    treeVars_["mu2Phi"] = mu2.phi();
	    treeVars_["mu2M"] = mu2.mass();
	    treeVars_["mu2Charge"] = mu2.charge();
	    treeVars_["mu2Iso"] = mu2Iso;
	    treeVars_["mu2Id"] = mu2Id;
	  }
	}
      }
      
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
      treeVars_["iPair"] = iPair;
      //fill tree for each pair of good muon and good tau
      tree_->Fill();
      iPair++;
    }
  }

  //std::cout << std::endl;
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODMuTauTriggerAnalyzer);
