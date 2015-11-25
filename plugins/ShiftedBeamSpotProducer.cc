
#include "ShiftedBeamSpotProducer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace edm;


ShiftedBeamSpotProducer::ShiftedBeamSpotProducer(const ParameterSet& iconf)
  : theMaxZ_( iconf.getParameter<double>("maxZ") )
  , srcToken_( consumes<reco::BeamSpot> ( iconf.getParameter<InputTag>("src") ) )
  , shiftX_( iconf.getParameter<double>("shiftX") )
  , shiftY_( iconf.getParameter<double>("shiftY") )
  , shiftZ_( iconf.getParameter<double>("shiftZ") )

{

  theMaxR2_ = iconf.getParameter<double>("maxRadius");
  theMaxR2_ *= theMaxR2_;

  produces<reco::BeamSpot>();

} 

ShiftedBeamSpotProducer::~ShiftedBeamSpotProducer() {}

void ShiftedBeamSpotProducer::produce(Event& iEvent, const EventSetup& iSetup)
{

  // get scalar collection
  Handle<reco::BeamSpot> handleSrc;
  iEvent.getByToken(srcToken_, handleSrc);

  reco::BeamSpot::Point aPoint;
  reco::BeamSpot::BeamType type = reco::BeamSpot::Unknown;
  double sigmaZ=0;
  double dxdz=0, dydz=0;
  double beamWidthX=0, beamWidthY=0; 
  reco::BeamSpot::CovarianceMatrix matrix;
  double emittX=0, emittY=0, bStar=0;
  if (handleSrc.isValid()){
    reco::BeamSpot beamSpot = *handleSrc;
    double x = beamSpot.x0()+shiftX_;
    double y = beamSpot.y0()+shiftY_;
    double z = beamSpot.z0()+shiftZ_;
    double r2 = x*x + y*y;
    if(theMaxR2_>0 && r2>theMaxR2_){
      x *= theMaxR2_/r2;
      y *= theMaxR2_/r2;
    }
    if(theMaxZ_>0 && z>theMaxZ_)
      z *= theMaxZ_/z;
    aPoint = reco::BeamSpot::Point(x,y,z);
    sigmaZ = beamSpot.sigmaZ();
    dxdz = beamSpot.dxdz();
    dydz = beamSpot.dydz();
    beamWidthX = beamSpot.BeamWidthX();
    beamWidthY = beamSpot.BeamWidthY();
    matrix =  beamSpot.covariance();
    type = beamSpot.type();
    emittX = beamSpot.emittanceX();
    emittY = beamSpot.emittanceX();
    bStar = beamSpot.betaStar();
  }else{
    edm::LogError("UnusableBeamSpot") << "No beam spot found in Event";
  }

  // product is a reco::BeamSpot object
  std::auto_ptr<reco::BeamSpot> result(new reco::BeamSpot);
  
  reco::BeamSpot aSpot = reco::BeamSpot(aPoint,
					sigmaZ,
					dxdz,
					dydz,
					beamWidthX,
					matrix);    
  aSpot.setBeamWidthY(beamWidthY);
  aSpot.setEmittanceX(emittX);
  aSpot.setEmittanceY(emittY);
  aSpot.setbetaStar(bStar);
  aSpot.setType(type);
    
  *result = aSpot;
  iEvent.put(result);

}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ShiftedBeamSpotProducer);
