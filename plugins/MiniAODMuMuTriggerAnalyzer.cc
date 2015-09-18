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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

class MiniAODMuMuTriggerAnalyzer : public edm::EDAnalyzer {
public:
  explicit MiniAODMuMuTriggerAnalyzer(const edm::ParameterSet&);
  ~MiniAODMuMuTriggerAnalyzer() {}

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  TTree *initTree(edm::Service<TFileService> &fs, std::string name="tree");
  void bookVariable(TTree *t=0, std::string var="foo");
  void cleanFilterVars();
  std::string triggerNameWithoutVersion(const std::string &triggerName);
  float muonIso(const pat::Muon &aMu, float dBetaFactor=0.5, bool allCharged=false);
  bool isGoodMuon(const pat::Muon &aMu, const reco::Vertex & vtx, bool useOldId=false);

  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

  edm::EDGetTokenT<pat::MuonCollection> muons_;
  edm::EDGetTokenT<pat::METCollection> mets_;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;
  edm::EDGetTokenT<GenEventInfoProduct> genEvtInfo_;
  
  std::set<std::string> tagTriggers_;
  std::vector<std::string> tagFilters_;
  std::vector<std::string> probeFilters_;
  
  std::map<std::string, float> treeVars_;
  
  double minPtTag_, maxEtaTag_, isoTag_;
  double minPtProbe_, maxEtaProbe_, isoProbe_;

  bool tagTightId_, probeTightId_;
  bool checkMCMatch_;
  bool isMC_;

  TTree *tree_;
};

MiniAODMuMuTriggerAnalyzer::MiniAODMuMuTriggerAnalyzer(const edm::ParameterSet& iConfig):
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  muons_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  mets_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("met"))),
  vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  genEvtInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEvtInfo"))),
  tagFilters_(iConfig.getParameter<std::vector<std::string> >("tagFilters")),
  probeFilters_(iConfig.getParameter<std::vector<std::string> >("probeFilters")),
  minPtTag_(iConfig.getParameter<double>("minPtTag")),
  maxEtaTag_(iConfig.getParameter<double>("maxEtaTag")),
  isoTag_(iConfig.getParameter<double>("isoTag")),
  minPtProbe_(iConfig.getParameter<double>("minPtProbe")),
  maxEtaProbe_(iConfig.getParameter<double>("maxEtaProbe")),
  isoProbe_(iConfig.getParameter<double>("isoProbe")),
  tagTightId_(iConfig.getParameter<bool>("tagTightId")),
  probeTightId_(iConfig.getParameter<bool>("probeTightId")),
  checkMCMatch_(iConfig.getParameter<bool>("checkMCMatch")),
  isMC_(iConfig.getParameter<bool>("isMC"))
{
  edm::Service<TFileService> fs;
  tree_ = initTree(fs,"muMuTriggerTree");
  bookVariable(tree_,"tagPt");
  bookVariable(tree_,"tagEta");
  bookVariable(tree_,"tagPhi");
  bookVariable(tree_,"tagM");
  bookVariable(tree_,"tagCharge");
  bookVariable(tree_,"probePt");
  bookVariable(tree_,"probeEta");
  bookVariable(tree_,"probePhi");
  bookVariable(tree_,"probeM");
  bookVariable(tree_,"probeCharge");
  bookVariable(tree_,"tagGenMatch");
  bookVariable(tree_,"probeGenMatch");
  bookVariable(tree_,"MEt");
  bookVariable(tree_,"MEtPhi");
  bookVariable(tree_,"Mass");
  bookVariable(tree_,"Vx");
  bookVariable(tree_,"Vy");
  bookVariable(tree_,"Vz");
  bookVariable(tree_,"nVtx");
  bookVariable(tree_,"weight");


  std::vector<std::string> tagTrgs( iConfig.getParameter<std::vector<std::string> >("tagTriggers") );
  for(unsigned int i=0; i<tagTrgs.size(); ++i)
    tagTriggers_.insert( triggerNameWithoutVersion(tagTrgs[i]) );
  for(std::set<std::string>::const_iterator it = tagTriggers_.begin(); it != tagTriggers_.end(); ++it)
    bookVariable(tree_,(*it) );
  for(unsigned int i=0; i<tagFilters_.size(); ++i)
    bookVariable(tree_,"tag_"+tagFilters_[i]);
  for(unsigned int i=0; i<probeFilters_.size(); ++i)
    bookVariable(tree_,"probe_"+probeFilters_[i]);
}

TTree * MiniAODMuMuTriggerAnalyzer::initTree(edm::Service<TFileService> &fs, std::string name)
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

void MiniAODMuMuTriggerAnalyzer::bookVariable(TTree *t, std::string var)
{
  if(!t) return;
  treeVars_[var];
  t->Branch(var.c_str(), &treeVars_[var], (var + "/F").c_str());

  return;
}

void MiniAODMuMuTriggerAnalyzer::cleanFilterVars(){
  for(std::string &label : tagFilters_)
    treeVars_["tag_"+label] = 0.;
  for(std::string &label : probeFilters_)
    treeVars_["probe_"+label] = 0.;
}

std::string MiniAODMuMuTriggerAnalyzer::triggerNameWithoutVersion(const std::string &triggerName)
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
float MiniAODMuMuTriggerAnalyzer::muonIso(const pat::Muon &aMu, float dBetaFactor, bool allCharged)
{
  float iso = std::max((float)0.0,aMu.userIsolation("pfPhotons") + aMu.userIsolation("pfNeutralHadrons") - dBetaFactor*aMu.userIsolation("pfPUChargedHadrons") );
  if(allCharged)
    iso += aMu.userIsolation("pfChargedAll");
  else
    iso += aMu.userIsolation("pfChargedHadrons");

  return iso;
}

bool MiniAODMuMuTriggerAnalyzer::isGoodMuon(const pat::Muon &aMu, const reco::Vertex & vtx, bool useOldId)
{
  //std::cout<<"Checking mu: pt="<<aMu.pt()<<", eta="<<aMu.eta()<<", phi="<<aMu.phi()<<std::endl;
  if(useOldId){
    // Tight mu
    if( !aMu.isGlobalMuon() || !aMu.isPFMuon() ) return false;
    if( !(aMu.normChi2() < 10.) ) return false;
    if( !(aMu.globalTrack()->hitPattern().numberOfValidMuonHits() > 0) ) return false;
    if( !(aMu.numberOfMatchedStations() > 1) ) return false;
    if( !(fabs( aMu.muonBestTrack()->dxy( vtx.position() ) ) < 0.2) ) return false; //FIXME: provide vertex
    if( !(fabs( aMu.muonBestTrack()->dz( vtx.position() ) ) < 0.5) ) return false; //FIXME: provide vertex
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

  return true;
}

void MiniAODMuMuTriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  edm::Handle<pat::MuonCollection> muons;
  iEvent.getByToken(muons_, muons);
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
    if(tagTriggers_.find( triggerNameWithoutVersion( names.triggerName(i) ) ) != tagTriggers_.end() ) {
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

  for(const pat::Muon &tag : *muons) {
    // kinematics
    if(tag.pt() < minPtTag_) continue;
    if(fabs( tag.eta() ) > maxEtaTag_) continue;
    // Iso
    if( muonIso(tag) > isoTag_ ) continue;
    // id
    if( !isGoodMuon(tag,vtxs->at(0),tagTightId_) ) continue;
    //std::cout<<"Tag muon: pt="<<tag.pt()<<", eta="<<tag.eta()<<", phi="<<tag.phi()<<std::endl;
    
    for(const pat::Muon &probe : *muons) {
      // kinematics
      if(probe.pt() < minPtProbe_) continue;
      if(fabs( probe.eta() ) > maxEtaProbe_) continue;
      // Iso
      if( muonIso(probe) > isoProbe_ ) continue;
      //pair
      if( deltaR2(tag,probe) < 0.5*0.5) continue;
      //FIXME OS??
      // id
      if( !isGoodMuon(probe,vtxs->at(0),probeTightId_) ) continue;
      //std::cout<<"Probe muon: pt="<<probe.pt()<<", eta="<<probe.eta()<<", phi="<<probe.phi()<<std::endl;
      //std::cout<<"DR(tag,probe)="<<deltaR(tag,probe)<<std::endl;
      treeVars_["tagPt"] = tag.pt();
      treeVars_["tagEta"] = tag.eta();
      treeVars_["tagPhi"] = tag.phi();
      treeVars_["tagM"] = tag.mass();
      treeVars_["tagCharge"] = tag.charge();
      treeVars_["probePt"] = probe.pt();
      treeVars_["probeEta"] = probe.eta();
      treeVars_["probePhi"] = probe.phi();
      treeVars_["probeM"] = probe.mass();
      treeVars_["probeCharge"] = probe.charge();
      treeVars_["Mass"] = (tag.p4()+probe.p4()).mass();
      // gen match
      if(!checkMCMatch_ || tag.genLepton() ) 
	treeVars_["tagGenMatch"] = 1;
      else
	treeVars_["tagGenMatch"] = 0;
      if(!checkMCMatch_ || probe.genLepton() ) 
	treeVars_["probeGenMatch"] = 1;
      else
	treeVars_["probeGenMatch"] = 0;

      //std::cout << "\n === TRIGGER OBJECTS === " << std::endl;
      cleanFilterVars();
      for(pat::TriggerObjectStandAlone trgObj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
	// check tag muon objects
	if( deltaR2(trgObj,tag) < 0.3*0.3 && trgObj.hasTriggerObjectType(trigger::TriggerMuon) ) {
	  for(std::string &label : tagFilters_){
	    if( trgObj.hasFilterLabel(label) ) { 
	      treeVars_["tag_"+label] = trgObj.pt();	    
	      /*
	      std::cout<<"\tFilter \""<<label<<"\" matched to tag muon, DR="<<deltaR(trgObj,tag)<<std::endl
		       <<"\t\tTrigger object: pt="<<trgObj.pt()<<", eta="<<trgObj.eta()<<", phi="<<trgObj.phi() << std::endl; 
	      */
	    }
	  }
	}
	// check probe muon objects
	if( deltaR2(trgObj,probe) < 0.3*0.3 && 
	    ( trgObj.hasTriggerObjectType(trigger::TriggerMuon) || trgObj.hasTriggerObjectType(trigger::TriggerL1Mu) ) ) {
	  for(std::string &label : probeFilters_){
	    if( trgObj.hasFilterLabel(label) ) { 
	      treeVars_["probe_"+label] = trgObj.pt();	    
	      /*
	      std::cout<<"\tFilter \""<<label<<"\" matched to probe muon, DR="<<deltaR(trgObj,probe)<<std::endl
		       <<"\t\tTrigger object: pt="<<trgObj.pt()<<", eta="<<trgObj.eta()<<", phi="<<trgObj.phi() << std::endl; 
	      */
	    }
	  }
	}
      }
      //fill tree for each pair of good muons (note: two muons can give two good tag-probe pairs)
      tree_->Fill();
    }
  }

  //std::cout << std::endl;
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODMuMuTriggerAnalyzer);
