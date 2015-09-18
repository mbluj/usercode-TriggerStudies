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
#include <utility>

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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

struct CutBasedElId{
  CutBasedElId(){};
  CutBasedElId(std::string name,
	       std::pair<float,float> dEtaInCut,
	       std::pair<float,float> dPhiInCut,
	       std::pair<float,float> sigmaIEtaIEtaCut,
	       std::pair<float,float> HoECut,
	       std::pair<float,float> EoPCut):
    name_(name),
    dEtaInCut_(dEtaInCut),dPhiInCut_(dPhiInCut),sigmaIEtaIEtaCut_(sigmaIEtaIEtaCut),
    HoECut_(HoECut),EoPCut_(EoPCut)
  {};
  CutBasedElId(std::string name,
	       float dEtaInCutEB, float dEtaInCutEE,
	       float dPhiInCutEB, float dPhiInCutEE,
	       float sigmaIEtaIEtaCutEB, float sigmaIEtaIEtaCutEE,
	       float HoECutEB, float HoECutEE,
	       float EoPCutEB, float EoPCutEE):
    name_(name),
    dEtaInCut_(dEtaInCutEB,dEtaInCutEE),dPhiInCut_(dPhiInCutEB,dPhiInCutEE),
    sigmaIEtaIEtaCut_(sigmaIEtaIEtaCutEB,sigmaIEtaIEtaCutEE),
    HoECut_(HoECutEB,HoECutEE),EoPCut_(EoPCutEB,EoPCutEE)
  {};
  ~CutBasedElId(){};
  std::string name_;
  std::pair<float,float> dEtaInCut_,dPhiInCut_,sigmaIEtaIEtaCut_,HoECut_,EoPCut_;
  void setName(std::string name) { name_=name; }
  bool checkIdEB(float &dEtaIn,
	       float &dPhiIn,
	       float &sigmaIEtaIEta,
	       float &HoE,
	       float &EoP){
    if(dEtaIn >= dEtaInCut_.first) return false;
    if(dPhiIn >= dPhiInCut_.first) return false;
    if(sigmaIEtaIEta >= sigmaIEtaIEtaCut_.first) return false;
    if(HoE >= HoECut_.first) return false;
    if(EoP >= EoPCut_.first) return false;
    return true;
  }
  bool checkIdEE(float &dEtaIn,
	       float &dPhiIn,
	       float &sigmaIEtaIEta,
	       float &HoE,
	       float &EoP){
    if(dEtaIn >= dEtaInCut_.second) return false;
    if(dPhiIn >= dPhiInCut_.second) return false;
    if(sigmaIEtaIEta >= sigmaIEtaIEtaCut_.second) return false;
    if(HoE >= HoECut_.second) return false;
    if(EoP >= EoPCut_.second) return false;
    return true;
  }
  bool checkId(float &dEtaIn,
               float &dPhiIn,
               float &sigmaIEtaIEta,
               float &HoE,
               float &EoP,
	       bool isEB){
    if(isEB) return checkIdEB(dEtaIn,dPhiIn,sigmaIEtaIEta,HoE,EoP);
    else return checkIdEE(dEtaIn,dPhiIn,sigmaIEtaIEta,HoE,EoP);
  }
};

class MiniAODEETriggerAnalyzer : public edm::EDAnalyzer {
public:
  explicit MiniAODEETriggerAnalyzer(const edm::ParameterSet&);
  ~MiniAODEETriggerAnalyzer() {}

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  TTree *initTree(edm::Service<TFileService> &fs, std::string name="tree");
  void bookVariable(TTree *t=0, std::string var="foo");
  void cleanFilterVars();
  std::string triggerNameWithoutVersion(const std::string &triggerName);
  float electronIso(const pat::Electron &aEl, float dBetaFactor=0.5, bool allCharged=false);
  bool electronId(const pat::Electron &aEl,
                  std::string idName="cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight");//FIXME
  bool isElIdSupported(std::string idName);
  bool isGoodElectron(const pat::Electron &aEl, const reco::Vertex & vtx, 
		      std::string idName="cutBasedElectronID-CSA14-PU20bx25-V0-standalone-tight");//FIXME


  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescales_;

  edm::EDGetTokenT<pat::ElectronCollection> electrons_;
  edm::EDGetTokenT<l1extra::L1EmParticleCollection> l1NoIsoEGs_, l1IsoEGs_;
  edm::EDGetTokenT<pat::METCollection> mets_;
  edm::EDGetTokenT<reco::VertexCollection> vertices_;

  std::set<std::string> tagTriggers_;
  std::vector<std::string> tagFilters_;
  std::vector<std::string> probeFilters_;
  
  std::map<std::string, float> treeVars_;
  
  double minPtTag_, maxEtaTag_, isoTag_;
  double minPtProbe_, maxEtaProbe_, isoProbe_;

  std::string tagId_, probeId_;
  bool checkMCMatch_;

  std::vector<CutBasedElId> cutBasedIds_;

  TTree *tree_;
};

MiniAODEETriggerAnalyzer::MiniAODEETriggerAnalyzer(const edm::ParameterSet& iConfig):
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  electrons_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
  l1NoIsoEGs_(consumes<l1extra::L1EmParticleCollection>(iConfig.getParameter<edm::InputTag>("l1NoIsoEGs"))),
  l1IsoEGs_(consumes<l1extra::L1EmParticleCollection>(iConfig.getParameter<edm::InputTag>("l1IsoEGs"))),
  mets_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("met"))),  
  vertices_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  tagFilters_(iConfig.getParameter<std::vector<std::string> >("tagFilters")),
  probeFilters_(iConfig.getParameter<std::vector<std::string> >("probeFilters")),
  minPtTag_(iConfig.getParameter<double>("minPtTag")),
  maxEtaTag_(iConfig.getParameter<double>("maxEtaTag")),
  isoTag_(iConfig.getParameter<double>("isoTag")),
  minPtProbe_(iConfig.getParameter<double>("minPtProbe")),
  maxEtaProbe_(iConfig.getParameter<double>("maxEtaProbe")),
  isoProbe_(iConfig.getParameter<double>("isoProbe")),
  tagId_(iConfig.getParameter<std::string>("tagId")),
  probeId_(iConfig.getParameter<std::string>("probeId")),
  checkMCMatch_(iConfig.getParameter<bool>("checkMCMatch"))
{
  edm::Service<TFileService> fs;
  tree_ = initTree(fs,"eETriggerTree");
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


  std::vector<std::string> tagTrgs( iConfig.getParameter<std::vector<std::string> >("tagTriggers") );
  for(unsigned int i=0; i<tagTrgs.size(); ++i)
    tagTriggers_.insert( triggerNameWithoutVersion(tagTrgs[i]) );
  for(std::set<std::string>::const_iterator it = tagTriggers_.begin(); it != tagTriggers_.end(); ++it)
    bookVariable(tree_,(*it) );
  for(unsigned int i=0; i<tagFilters_.size(); ++i)
    bookVariable(tree_,"tag_"+tagFilters_[i]);
  for(unsigned int i=0; i<probeFilters_.size(); ++i)
    bookVariable(tree_,"probe_"+probeFilters_[i]);
  bookVariable(tree_,"probe_l1NoIsoEG");
  bookVariable(tree_,"probe_l1IsoEG");

}

TTree * MiniAODEETriggerAnalyzer::initTree(edm::Service<TFileService> &fs, std::string name)
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

void MiniAODEETriggerAnalyzer::bookVariable(TTree *t, std::string var)
{
  if(!t) return;
  treeVars_[var];
  t->Branch(var.c_str(), &treeVars_[var], (var + "/F").c_str());

  return;
}

void MiniAODEETriggerAnalyzer::cleanFilterVars(){
  for(std::string &label : tagFilters_)
    treeVars_["tag_"+label] = 0.;
  for(std::string &label : probeFilters_)
    treeVars_["probe_"+label] = 0.;
   treeVars_["probe_l1NoIsoEG"] = 0.;
   treeVars_["probe_l1IsoEG"] = 0.;
}

std::string MiniAODEETriggerAnalyzer::triggerNameWithoutVersion(const std::string &triggerName)
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
float MiniAODEETriggerAnalyzer::electronIso(const pat::Electron &aEl, float dBetaFactor, bool allCharged)
{
  float iso = std::max((float)0.0,aEl.userIsolation("pfPhotons") + aEl.userIsolation("pfNeutralHadrons") - dBetaFactor*aEl.userIsolation("pfPUChargedHadrons") );
  if(allCharged)
    iso += aEl.userIsolation("pfChargedAll");
  else
    iso += aEl.userIsolation("pfChargedHadrons");

  return iso;
}
bool MiniAODEETriggerAnalyzer::electronId(const pat::Electron &aEl, std::string idName){//FIXME

  return aEl.electronID(idName);
}

bool MiniAODEETriggerAnalyzer::isElIdSupported(std::string idName){
  std::set<std::string> supportedElId;
  supportedElId.insert("POG_2012_Veto");
  supportedElId.insert("POG_2012_Loose");
  supportedElId.insert("POG_2012_Medium");
  supportedElId.insert("POG_2012_Tight");
  supportedElId.insert("POG_CSA14_25ns_v1_Veto");
  supportedElId.insert("POG_CSA14_25ns_v1_Loose");
  supportedElId.insert("POG_CSA14_25ns_v1_Medium");
  supportedElId.insert("POG_CSA14_25ns_v1_Tight");
  supportedElId.insert("POG_CSA14_50ns_v1_Veto");
  supportedElId.insert("POG_CSA14_50ns_v1_Loose");
  supportedElId.insert("POG_CSA14_50ns_v1_Medium");
  supportedElId.insert("POG_CSA14_50ns_v1_Tight");
  supportedElId.insert("POG_PHYS14_25ns_v1_Veto");
  supportedElId.insert("POG_PHYS14_25ns_v1_Loose");
  supportedElId.insert("POG_PHYS14_25ns_v1_Medium");
  supportedElId.insert("POG_PHYS14_25ns_v1_Tight");
  //firstly check among custom IDs
  if(supportedElId.find(idName) != supportedElId.end()) return true;
  //then inside patObj 
  //how??

  return false;
}

bool MiniAODEETriggerAnalyzer::isGoodElectron(const pat::Electron &aEl, const reco::Vertex & vtx, std::string elId)
{
  //std::cout<<"Checking e: pt="<<aEl.pt()<<", eta="<<aEl.eta()<<", phi="<<aEl.phi()<<std::endl;
  // e-Id
  if(elId != "")
    if( !electronId(aEl, elId) ) return false;  //FIXME
  // Vertex
  if( !(fabs( aEl.gsfTrack()->dxy( vtx.position() ) ) < 0.045) ) return false;
  if( !(fabs( aEl.gsfTrack()->dz( vtx.position() ) ) < 0.2) ) return false;

  return true;
}

void MiniAODEETriggerAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;

  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);

  edm::Handle<pat::ElectronCollection> electrons;
  iEvent.getByToken(electrons_, electrons);
  edm::Handle<l1extra::L1EmParticleCollection> l1NoIsoEGs;
  iEvent.getByToken(l1NoIsoEGs_, l1NoIsoEGs);
  edm::Handle<l1extra::L1EmParticleCollection> l1IsoEGs;
  iEvent.getByToken(l1IsoEGs_, l1IsoEGs);

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

  for(const pat::Electron &tag : *electrons) {
    // kinematics
    if(tag.pt() < minPtTag_) continue;
    if(fabs( tag.eta() ) > maxEtaTag_) continue;
    // Iso
    if( electronIso(tag) > isoTag_ ) continue;
    // id
    if( !isGoodElectron(tag,vtxs->at(0)) ) continue;
    //std::cout<<"Tag electron: pt="<<tag.pt()<<", eta="<<tag.eta()<<", phi="<<tag.phi()<<std::endl;
    
    for(const pat::Electron &probe : *electrons) {
      // kinematics
      if(probe.pt() < minPtProbe_) continue;
      if(fabs( probe.eta() ) > maxEtaProbe_) continue;
      // Iso
      if( electronIso(probe) > isoProbe_ ) continue;
      //pair
      if( deltaR2(tag,probe) < 0.5*0.5) continue;
      //FIXME OS??
      // id
      if( !isGoodElectron(probe,vtxs->at(0)) ) continue;
      //std::cout<<"Probe electron: pt="<<probe.pt()<<", eta="<<probe.eta()<<", phi="<<probe.phi()<<std::endl;
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
	// check tag electron objects
	if( deltaR2(trgObj,tag) < 0.5*0.5 && 
	    ( trgObj.hasTriggerObjectType(trigger::TriggerElectron) || trgObj.hasTriggerObjectType(trigger::TriggerPhoton) ) ) {
	  for(std::string &label : tagFilters_){
	    if( trgObj.hasFilterLabel(label) ) { 
	      treeVars_["tag_"+label] = trgObj.pt();	    
	      /*
	      std::cout<<"\tFilter \""<<label<<"\" matched to tag electron, DR="<<deltaR(trgObj,tag)<<std::endl
		       <<"\t\tTrigger object: pt="<<trgObj.pt()<<", eta="<<trgObj.eta()<<", phi="<<trgObj.phi() << std::endl; 
	      */
	    }
	  }
	}
	// check probe electron objects
	if( deltaR2(trgObj,probe) < 0.5*0.5 && 
	    ( trgObj.hasTriggerObjectType(trigger::TriggerElectron) || trgObj.hasTriggerObjectType(trigger::TriggerPhoton) || 
	      trgObj.hasTriggerObjectType(trigger::TriggerL1IsoEG) || trgObj.hasTriggerObjectType(trigger::TriggerL1NoIsoEG) ) ) {
	  for(std::string &label : probeFilters_){
	    if( trgObj.hasFilterLabel(label) ) { 
	      treeVars_["probe_"+label] = trgObj.pt();	    
	      /*
	      std::cout<<"\tFilter \""<<label<<"\" matched to probe electron, DR="<<deltaR(trgObj,probe)<<std::endl
		       <<"\t\tTrigger object: pt="<<trgObj.pt()<<", eta="<<trgObj.eta()<<", phi="<<trgObj.phi() << std::endl; 
	      */
	    }
	  }
	}
      }
      //matching to L1Extra
      //NoIso 
      float maxL1Et = 0;
      for(const l1extra::L1EmParticle &l1EG : *l1NoIsoEGs) {
	if(l1EG.et()>maxL1Et && deltaR2(l1EG,probe) < 0.5*0.5){
	  maxL1Et = l1EG.et();
	  treeVars_["probe_l1NoIsoEG"] = l1EG.et();
	}
      }
      //NoIso 
      maxL1Et = 0;
      for(const l1extra::L1EmParticle &l1EG : *l1IsoEGs) {
	if(l1EG.et()>maxL1Et && deltaR2(l1EG,probe) < 0.5*0.5){
	  maxL1Et = l1EG.et();
	  treeVars_["probe_l1IsoEG"] = l1EG.et();
	}
      }
      //fill tree for each pair of good electrons (note: two electrons can give two good tag-probe pairs)
      tree_->Fill();
    }
  }

  //std::cout << std::endl;
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODEETriggerAnalyzer);
