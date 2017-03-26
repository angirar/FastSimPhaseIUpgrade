#include <vector>
#include <memory>
#include <algorithm>
#include <TFile.h>
#include <TROOT.h>
#include <TH2.h>
#include <iostream>
#include <fstream>
// framework
#include "FWCore/Framework/interface/Event.h"
#include "MagneticField/UniformEngine/src/UniformMagneticField.h"

// tracking
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/GeomPropagators/interface/HelixArbitraryPlaneCrossing.h"

// fastsim
#include "FastSimulation/TrajectoryManager/interface/InsideBoundsMeasurementEstimator.h" //TODO move this
#include "FastSimulation/SimplifiedGeometryPropagator/interface/Particle.h"
#include "FastSimulation/SimplifiedGeometryPropagator/interface/SimplifiedGeometry.h"
#include "FastSimulation/SimplifiedGeometryPropagator/interface/InteractionModel.h"
#include "FastSimulation/SimplifiedGeometryPropagator/interface/InteractionModelFactory.h"
#include "FastSimulation/SimplifiedGeometryPropagator/interface/Constants.h"

// data formats
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

// other
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "CondFormats/External/interface/DetID.h"
#include "FWCore/Framework/interface/ProducerBase.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/MakerMacros.h"

namespace edm
{
    class ParameterSet;
}

typedef std::pair<const GeomDet*,TrajectoryStateOnSurface> DetWithState;

namespace fastsim
{
    class TrackerSimHitProducer : public InteractionModel
    {
    public:
	TrackerSimHitProducer(const std::string & name,const edm::ParameterSet & cfg);
	~TrackerSimHitProducer(){;}
	void interact(Particle & particle,const SimplifiedGeometry & layer,std::vector<std::unique_ptr<Particle> > & secondaries,const RandomEngineAndDistribution & random) override;
      virtual void registerProducts(edm::ProducerBase & producer) const override;
      virtual void storeProducts(edm::Event & iEvent) override;
	std::pair<double, PSimHit*> createHitOnDetector(const TrajectoryStateOnSurface & particle,int pdgId,int simTrackId,const GeomDet & detector, GlobalPoint & refPos);
    private:
	const float onSurfaceTolerance_;
      std::unique_ptr<edm::PSimHitContainer> simHitContainer_;
      edm::Service<TFileService> FileService;
      TH2F* hitsZPerpfull;
      TH2F* hitsZPerp;
      TH2F* hitsXY;
      TH2F* lochitsZPerp;
    };
}



fastsim::TrackerSimHitProducer::TrackerSimHitProducer(const std::string & name,const edm::ParameterSet & cfg)
    : fastsim::InteractionModel(name)
    , onSurfaceTolerance_(0.01)
    , simHitContainer_(new edm::PSimHitContainer)
{
  edm::Service<TFileService> fs;
  hitsZPerpfull = fs->make<TH2F>("simhitsZPerpfull","",1280,-320,320,1040,-130,130);
  hitsZPerp = fs->make<TH2F>("simhitsZPerp","",600,-60,60,600,-60,60);
  hitsXY = fs->make<TH2F>("simhitsXY","",750,-130,130,750,-130,130);
  lochitsZPerp = fs->make<TH2F>("local_simhitsZPerp","",1280,-320,320,520,0,130);
}

void fastsim::TrackerSimHitProducer::registerProducts(edm::ProducerBase & producer) const
{
    producer.produces<edm::PSimHitContainer>("TrackerHits");
}

void fastsim::TrackerSimHitProducer::storeProducts(edm::Event & iEvent)
{
  
  iEvent.put(std::move(simHitContainer_), "TrackerHits");
  simHitContainer_.reset(new edm::PSimHitContainer);
}

void fastsim::TrackerSimHitProducer::interact(Particle & particle,const SimplifiedGeometry & layer,std::vector<std::unique_ptr<Particle> > & secondaries,const RandomEngineAndDistribution & random)
{
    //
    // check that layer has tracker modules
    //
    if(!layer.getDetLayer())
    {
	return;
    }

    //
    // no material
    //
    if(layer.getThickness(particle.position(), particle.momentum()) < 1E-10)
    {
    return;
    }

    //
    // create the trajectory of the particle
    //
    UniformMagneticField magneticField(layer.getMagneticFieldZ(particle.position())); 
    GlobalPoint  position( particle.position().X(), particle.position().Y(), particle.position().Z());
    GlobalVector momentum( particle.momentum().Px(), particle.momentum().Py(), particle.momentum().Pz());
    auto plane = layer.getDetLayer()->surface().tangentPlane(position);
    TrajectoryStateOnSurface trajectory(GlobalTrajectoryParameters( position, momentum, TrackCharge( particle.charge()), &magneticField), *plane);
    
    //
    // find detectors compatible with the particle's trajectory
    //
    AnalyticalPropagator propagator(&magneticField, anyDirection);
    InsideBoundsMeasurementEstimator est;
    std::vector<DetWithState> compatibleDetectors = layer.getDetLayer()->compatibleDets(trajectory, propagator, est);

    ////////
    // You have to sort the simHits in the order they occur!
    ////////

    // The old algorithm (sorting by distance to IP) doesn't seem to make sense to me (what if particle moves inwards??)

    // Detector layers have to be sorted by proximity to particle.position
    // Doesn't always work! Particle could be already have been propagated in between the layers!
    // Proximity to previous hit also doesn't work since simHits only store the localPosition
    // Propagate particle backwards a bit to make sure it's outside any components (straight line should work well enough) 
    std::map<double, PSimHit*> distAndHits;
    // Position relative to which the hits should be sorted
    GlobalPoint positionOutside(particle.position().x()-particle.momentum().x()/particle.momentum().Rho()*10.,
                                particle.position().y()-particle.momentum().y()/particle.momentum().Rho()*10.,
                                particle.position().z()-particle.momentum().z()/particle.momentum().Rho()*10.);
    //
    // loop over the compatible detectors
    //
    for (const auto & detectorWithState : compatibleDetectors)
    {
    	const GeomDet & detector = *detectorWithState.first;
    	const TrajectoryStateOnSurface & particleState = detectorWithState.second;
    	// if the detector has no components
    	if(detector.isLeaf())
    	{
    	    std::pair<double, PSimHit*> hitPair = createHitOnDetector(particleState,particle.pdgId(),particle.simTrackIndex(),detector,positionOutside);
    	    if(hitPair.second){
    		  	distAndHits.insert(distAndHits.end(), hitPair);
    		  }
    	}
    	else
    	{
    	    // if the detector has components
    	    for(const auto component : detector.components())
    	    {
    		  std::pair<double, PSimHit*> hitPair = createHitOnDetector(particleState,particle.pdgId(),particle.simTrackIndex(),*component,positionOutside);    		  
    		  if(hitPair.second){
    		  	distAndHits.insert(distAndHits.end(), hitPair);
    		  }
    	    }
    	}
    }

    // Fill simHitContainer
    for(std::map<double, PSimHit*>::const_iterator it = distAndHits.begin(); it != distAndHits.end(); it++){
    	simHitContainer_->push_back(*(it->second));
    }
    
}

// Also returns distance to simHit since hits have to be ordered (in time) afterwards. Necessary to do explicit copy of TrajectoryStateOnSurface particle (not call by reference)
std::pair<double, PSimHit*> fastsim::TrackerSimHitProducer::createHitOnDetector(const TrajectoryStateOnSurface & particle, int pdgId, int simTrackId, const GeomDet & detector, GlobalPoint & refPos)
{
    //
    // determine position and momentum of particle in the coordinate system of the detector
    //
    LocalPoint localPosition;
    LocalVector localMomentum;
    // if the particle is close enough, no further propagation is needed
    if ( fabs( detector.toLocal(particle.globalPosition()).z()) < onSurfaceTolerance_) 
    {
	   localPosition = particle.localPosition();
	   localMomentum = particle.localMomentum();
    }
    // else, propagate 
    else 
    {
    	// find crossing of particle with 
    	HelixArbitraryPlaneCrossing crossing(particle.globalPosition().basicVector(),
    					      particle.globalMomentum().basicVector(),
    					      particle.transverseCurvature(),
    					      anyDirection);
    	std::pair<bool,double> path = crossing.pathLength(detector.surface());
    	// case propagation succeeds
    	if (path.first) 	
    	{
    	    localPosition = detector.toLocal( GlobalPoint( crossing.position(path.second)));
    	    localMomentum = detector.toLocal( GlobalVector( crossing.direction(path.second)));
    	    localMomentum = localMomentum.unit() * particle.localMomentum().mag();
    	}
    	// case propagation fails
    	else
    	{
    	    return std::pair<double, PSimHit*>(0, 0);
    	}
    }

    // 
    // find entry and exit point of particle in detector
    //
    const Plane& detectorPlane = detector.surface();
    float halfThick = 0.5*detectorPlane.bounds().thickness();
    float pZ = localMomentum.z();
    LocalPoint entry = localPosition + (-halfThick/pZ) * localMomentum;
    LocalPoint exit = localPosition + halfThick/pZ * localMomentum;
    float tof = particle.globalPosition().mag() / fastsim::Constants::speedOfLight ; // in nanoseconds
    
    //
    // make sure the simhit is physically on the module
    //
    double boundX = detectorPlane.bounds().width()/2.;
    double boundY = detectorPlane.bounds().length()/2.;
    // Special treatment for TID and TEC trapeziodal modules
    unsigned subdet = DetId(detector.geographicalId()).subdetId(); 
    //    std::cout<<"subdet="<<subdet<<std::endl;
    if ( subdet == 4 || subdet == 6 ) 
	boundX *=  1. - localPosition.y()/detectorPlane.position().perp();
    if(fabs(localPosition.x()) > boundX  || fabs(localPosition.y()) > boundY )
    {
	   return std::pair<double, PSimHit*>(0, 0);
    }

    //
    // create the hit
    //
    double energyDeposit = 0.; // do something about the energy deposit

    GlobalPoint hitPos(detector.surface().toGlobal(localPosition));
    //std::cout<<hitPos.x()<<","<<hitPos.y()<<std::endl;
    double r = std::sqrt(hitPos.x()*hitPos.x()+hitPos.y()*hitPos.y());
    if(hitPos.phi()<0){
      hitsZPerpfull->Fill(hitPos.z(),-r);
      if( subdet == 1 || subdet == 2)
	hitsZPerp->Fill(hitPos.z(),-r);
    }
    else{
      hitsZPerpfull->Fill(hitPos.z(),r);
      if( subdet == 1 || subdet == 2)
	hitsZPerp->Fill(hitPos.z(),r);
    }
    hitsXY->Fill(hitPos.x(),hitPos.y());
    lochitsZPerp->Fill(localPosition.z(),std::sqrt(localPosition.x()*localPosition.x()+localPosition.y()*localPosition.y()));

    return std::pair<double, PSimHit*>((hitPos-refPos).mag(),
                                        new PSimHit(entry, exit, localMomentum.mag(), tof, energyDeposit, pdgId,
                                				   detector.geographicalId().rawId(),simTrackId,
                                				   localMomentum.theta(),
                                				   localMomentum.phi()));
}

DEFINE_EDM_PLUGIN(
    fastsim::InteractionModelFactory,
    fastsim::TrackerSimHitProducer,
    "fastsim::TrackerSimHitProducer"
    );
