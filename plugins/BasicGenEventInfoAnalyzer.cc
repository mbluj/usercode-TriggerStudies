//
// M. Bluj, National Centre for Nuclear Research, Poland
//

// system include files
#include <memory>
#include <cmath>
#include <iostream>
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

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"
#include "TH1F.h"

class BasicGenEventInfoAnalyzer : public edm::EDAnalyzer {
public:
  explicit BasicGenEventInfoAnalyzer(const edm::ParameterSet&);
  ~BasicGenEventInfoAnalyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;

  edm::EDGetTokenT<GenEventInfoProduct> genEvtInfo_;

  bool isMC_;

  //TTree *tree_;
  TH1F *hCounts_;
};

BasicGenEventInfoAnalyzer::BasicGenEventInfoAnalyzer(const edm::ParameterSet& iConfig):
  genEvtInfo_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEvtInfo"))),
  isMC_(iConfig.getParameter<bool>("isMC"))
{
  edm::Service<TFileService> fs;
  hCounts_ = fs->make<TH1F>("hCounts","event counts",10,0,10);
  
}

BasicGenEventInfoAnalyzer::~BasicGenEventInfoAnalyzer() {

  //set errors
  for(int i=0; i<hCounts_->GetNbinsX()+1; ++i){
    double count=hCounts_->GetBinContent(i);
    hCounts_->SetBinError(i,sqrt(count));
    //std::cout<<i<<". "<<count<<" +- "<<hCounts_->GetBinError(i)<<std::endl;
  }
}

void BasicGenEventInfoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  double evtWeight=1;
  if(isMC_){
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    iEvent.getByToken(genEvtInfo_,genEvtInfo);

    //get event weight
    evtWeight = genEvtInfo->weight();
  }
  //counts events
  hCounts_->SetBinContent(1,hCounts_->GetBinContent(1)+1);
  //std::cout<<"No. of events: "<<hCounts_->GetBinContent(1)<<std::endl;

  //sum weights
  hCounts_->SetBinContent(2,hCounts_->GetBinContent(2)+evtWeight);
  //std::cout<<"weight: "<<evtWeight
  //	   <<", sum of weights: "<< hCounts_->GetBinContent(2)
  //	   <<std::endl;

  //std::cout << std::endl;
  
}

//define this as a plug-in
DEFINE_FWK_MODULE(BasicGenEventInfoAnalyzer);
