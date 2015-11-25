#ifndef TriggerStudiesTauPlugins_ShiftedBeamSpotProducer_h
#define TriggerStudiesTauPlugins_ShiftedBeamSpotProducer_h

/**_________________________________________________________________
   class:   ShiftedBeamSpotProducer.h
  
   description: Produces beamspot shifted accordingly to user settings 

   author: Michal Bluj, NCBJ, Poland (michal.bluj@cern.ch)

________________________________________________________________**/

#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"


class ShiftedBeamSpotProducer: public edm::stream::EDProducer<> {

  public:

	/// constructor
	explicit ShiftedBeamSpotProducer(const edm::ParameterSet& iConf);
	/// destructor
	~ShiftedBeamSpotProducer();
	
	/// produce a beam spot class
	virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;

  private:

	const double theMaxZ_;
	double theMaxR2_;
	const edm::EDGetTokenT<reco::BeamSpot> srcToken_;
	const double shiftX_, shiftY_, shiftZ_;

};

#endif
