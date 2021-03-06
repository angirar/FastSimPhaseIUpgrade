// system include files
#include <memory>
#include <TFile.h>
#include <TROOT.h>
#include <TH2.h>
#include <TH1.h>
#include <iostream>
#include <fstream>
// framework stuff
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

// fast tracker rechits
#include "DataFormats/TrackerRecHit2D/interface/FastTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/FastMatchedTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/FastProjectedTrackerRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/FastTrackerRecHitCollection.h"

// geometry stuff
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/GluedGeomDet.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// sim stuff
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

class FastTrackerRecHitMatcher : public edm::stream::EDProducer<>  {

    public:

    explicit FastTrackerRecHitMatcher(const edm::ParameterSet&);
    ~FastTrackerRecHitMatcher(){;}
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
    
    virtual void produce(edm::Event&, const edm::EventSetup&) override;
  TH2F* rechitsRZfull;
  TH2F* rechitsRZ;
  TH2F* rechitsxy;
  TH1F* resol[4];

    // ---------- typedefs -----------------------------
    typedef std::pair<LocalPoint,LocalPoint>                   StripPosition; 

    // ----------internal functions ---------------------------

    // create projected hit
    std::unique_ptr<FastTrackerRecHit> projectOnly( const FastSingleTrackerRecHit *originalRH,
						  const GeomDet * monoDet,
						  const GluedGeomDet* gluedDet,
						  LocalVector& ldir) const;
    
    // create matched hit
    std::unique_ptr<FastTrackerRecHit> match( const FastSingleTrackerRecHit *monoRH,
					    const FastSingleTrackerRecHit *stereoRH,
					    const GluedGeomDet* gluedDet,
					    LocalVector& trackdirection,
					    bool stereLayerFirst) const;
    
    StripPosition project(const GeomDetUnit *det,
			  const GluedGeomDet* glueddet,
			  const StripPosition& strip,
			  const LocalVector& trackdirection)const;
    
    inline const FastSingleTrackerRecHit * _cast2Single(const FastTrackerRecHit * recHit) const{
	if(!recHit->isSingle()){
	    throw cms::Exception("FastTrackerRecHitMatcher") << "all rechits in simHit2RecHitMap must be instances of FastSingleTrackerRecHit. recHit's rtti: " << recHit->rtti() << std::endl;
	}
	return dynamic_cast<const FastSingleTrackerRecHit *>(recHit);
    }

    // ----------member data ---------------------------
    edm::EDGetTokenT<edm::PSimHitContainer> simHitsToken; 
    edm::EDGetTokenT<FastTrackerRecHitRefCollection> simHit2RecHitMapToken;

};

FastTrackerRecHitMatcher::FastTrackerRecHitMatcher(const edm::ParameterSet& iConfig)

{
  edm::Service<TFileService> fs;
  rechitsRZfull = fs->make<TH2F>("matcherrechitsZPerpfull","",1280,-320,320,1040,-130,130);
  rechitsRZ = fs->make<TH2F>("matcherrechitsZPerp","",600,-60,60,600,-60,60);
  rechitsxy = fs->make<TH2F>("matcherrechitsXY","",1500,-250,250,750,-130,130);
  resol[0] = fs->make<TH1F>("Xresol","Resolution plot b/w xpos of Simhit & Rechit",10000,-5,5);
  resol[1] = fs->make<TH1F>("Yresol","Resolution plot b/w ypos of Simhit & Rechit",10000,-5,5);
  resol[2] = fs->make<TH1F>("Zresol","Resolution plot b/w zpos of Simhit & Rechit",10000,-5,5);
  resol[3] = fs->make<TH1F>("Rresol","Resolution plot b/w R of Simhit & Rechit",10000,-5,5);

    simHitsToken = consumes<edm::PSimHitContainer>(iConfig.getParameter<edm::InputTag>("simHits"));
    simHit2RecHitMapToken = consumes<FastTrackerRecHitRefCollection>(iConfig.getParameter<edm::InputTag>("simHit2RecHitMap"));
    
    produces<FastTrackerRecHitCollection>();
    produces<FastTrackerRecHitRefCollection>("simHit2RecHitMap");
}


void FastTrackerRecHitMatcher::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    // services
    edm::ESHandle<TrackerGeometry> geometry;
    iSetup.get<TrackerDigiGeometryRecord> ().get (geometry);

    // input
    edm::Handle<edm::PSimHitContainer> simHits;
    iEvent.getByToken(simHitsToken,simHits);

    edm::Handle<FastTrackerRecHitRefCollection> simHit2RecHitMap;
    iEvent.getByToken(simHit2RecHitMapToken,simHit2RecHitMap);
    
    // output
    auto output_recHits = std::make_unique<FastTrackerRecHitCollection>();
    auto output_simHit2RecHitMap = std::make_unique<FastTrackerRecHitRefCollection>(simHit2RecHitMap->size(),FastTrackerRecHitRef());
    edm::RefProd<FastTrackerRecHitCollection> output_recHits_refProd = iEvent.getRefBeforePut<FastTrackerRecHitCollection>();

    bool skipNext = false;
    for(unsigned simHitCounter = 0;simHitCounter < simHits->size();++simHitCounter){

	// skip hit in case it was matched to previous one
	if(skipNext){
	    skipNext = false;
	    continue;
	}
	skipNext = false;

	// get simHit and associated recHit
	const PSimHit & simHit = (*simHits)[simHitCounter];
	const FastTrackerRecHitRef & recHitRef = (*simHit2RecHitMap)[simHitCounter];

	const LocalPoint& localPoint = simHit.localPosition();
	//std::cout<<"Found localpos"<<std::endl;
	const DetId& theDetId = simHit.detUnitId();
	//std::cout<<"Found det ID"<<std::endl;
	const GeomDet* theGeomDet = geometry->idToDet(theDetId);
	//std::cout<<"Found pointer to Geom det from ID"<<std::endl;
	const GlobalPoint& globalPoint = theGeomDet->toGlobal(localPoint);
	
	//simhit x,y,z
	double simx = globalPoint.x();
        double simy = globalPoint.y();
        double simz = globalPoint.z();
	double simr = globalPoint.mag();
	// skip simHits w/o associated recHit 
	if(recHitRef.isNull())
	    continue;

	// cast
	const FastSingleTrackerRecHit * recHit = _cast2Single(recHitRef.get());
	
	double x = recHit->globalPosition().x();
	double y = recHit->globalPosition().y();
	double z = recHit->globalPosition().z();
	double mag = recHit->globalPosition().mag();
	//resolution plots
	resol[0]->Fill(simx-x);
	resol[1]->Fill(simy-y);
	resol[2]->Fill(simz-z);
	resol[3]->Fill(simr-mag);

	double phi_val = recHit->globalPosition().phi();
	double r = std::sqrt(x*x+y*y);
	
	// get subdetector id
	DetId detid = recHit->geographicalId();
	unsigned int subdet = detid.subdetId();
	
	if(phi_val<0){
	  rechitsRZfull->Fill(z,-r);
	  if( subdet == 1 || subdet == 2)
	    rechitsRZ->Fill(z,-r);
	}
	else{
	  rechitsRZfull->Fill(z,r);
	  if( subdet == 1 || subdet == 2)
	    rechitsRZ->Fill(z,r);
	}
	rechitsxy->Fill(x,y);

	// treat pixel hits
	if(subdet <= 2){ 
	    (*output_simHit2RecHitMap)[simHitCounter] = recHitRef;
	}

	// treat strip hits
	else{

	    StripSubdetector stripSubDetId(detid);

	    // treat regular regular strip hits
	    if(!stripSubDetId.glued()){ 
		(*output_simHit2RecHitMap)[simHitCounter] = recHitRef;
	    }
	    
	    // treat strip hits on glued layers
	    else{ 
		
		// Obtain direction of simtrack at simhit in local coordinates of glued module
		//   - direction of simtrack at simhit, in coordindates of the single module
		LocalVector localSimTrackDir = simHit.localDirection();
		//   - transform to global coordinates
		GlobalVector globalSimTrackDir= recHit->det()->surface().toGlobal(localSimTrackDir);
		//   - transform to local coordinates of glued module
		const GluedGeomDet* gluedDet = (const GluedGeomDet*)geometry->idToDet(DetId(stripSubDetId.glued()));
		LocalVector gluedLocalSimTrackDir = gluedDet->surface().toLocal(globalSimTrackDir);
		
		// check whether next hit is partner
		const FastSingleTrackerRecHit * partnerRecHit = 0;
		//      - there must be a next hit
		if(simHitCounter + 1 < simHits->size()){
		    const FastTrackerRecHitRef & nextRecHitRef = (*simHit2RecHitMap)[simHitCounter + 1];
		    const PSimHit & nextSimHit = (*simHits)[simHitCounter + 1];
		    //  - partner hit must not be null
		    //  - simHit and partner simHit must belong to same simTrack 
		    //  - partner hit must be on the module glued to the module of the hit
		    if( (!nextRecHitRef.isNull())
		        && simHit.trackId() == nextSimHit.trackId()
			&& StripSubdetector( nextRecHitRef->geographicalId() ).partnerDetId() == detid.rawId() ) {
			partnerRecHit = _cast2Single(nextRecHitRef.get());
			skipNext = true;
		    }
		}
		
                std::unique_ptr<FastTrackerRecHit> newRecHit(nullptr);
		
		// if partner found: create a matched hit
		if( partnerRecHit ){
		    newRecHit = match( stripSubDetId.stereo() ? partnerRecHit : recHit,
				       stripSubDetId.stereo() ? recHit : partnerRecHit,
				       gluedDet  , gluedLocalSimTrackDir,
				       stripSubDetId.stereo());
		}
		// else: create projected hit
		else{
		    newRecHit = projectOnly( recHit , geometry->idToDet(detid),gluedDet, gluedLocalSimTrackDir  );
		}
		output_recHits->push_back(std::move(newRecHit));
		(*output_simHit2RecHitMap)[simHitCounter] = FastTrackerRecHitRef(output_recHits_refProd,output_recHits->size()-1);
	    }
	}
    } 
    
    iEvent.put(std::move(output_recHits));
    iEvent.put(std::move(output_simHit2RecHitMap),"simHit2RecHitMap");

}

std::unique_ptr<FastTrackerRecHit> FastTrackerRecHitMatcher::match(const FastSingleTrackerRecHit *monoRH,
								 const FastSingleTrackerRecHit *stereoRH,
								 const GluedGeomDet* gluedDet,
								 LocalVector& trackdirection,
								 bool stereoHitFirst) const
{

    // stripdet = mono
    // partnerstripdet = stereo
    const GeomDetUnit* stripdet = gluedDet->monoDet();
    const GeomDetUnit* partnerstripdet = gluedDet->stereoDet();
    const StripTopology& topol=(const StripTopology&)stripdet->topology();

    LocalPoint position;    

    // position of the initial and final point of the strip (RPHI cluster) in local strip coordinates
    MeasurementPoint RPHIpoint=topol.measurementPosition(monoRH->localPosition());
    MeasurementPoint RPHIpointini=MeasurementPoint(RPHIpoint.x(),-0.5);
    MeasurementPoint RPHIpointend=MeasurementPoint(RPHIpoint.x(),0.5);
  
    // position of the initial and final point of the strip in local coordinates (mono det)
    StripPosition stripmono=StripPosition(topol.localPosition(RPHIpointini),topol.localPosition(RPHIpointend));

    if(trackdirection.mag2()<FLT_MIN){// in case of no track hypothesis assume a track from the origin through the center of the strip
	LocalPoint lcenterofstrip=monoRH->localPosition();
	GlobalPoint gcenterofstrip=(stripdet->surface()).toGlobal(lcenterofstrip);
	GlobalVector gtrackdirection=gcenterofstrip-GlobalPoint(0,0,0);
	trackdirection=(gluedDet->surface()).toLocal(gtrackdirection);
    }
 
    //project mono hit on glued det
    StripPosition projectedstripmono=project(stripdet,gluedDet,stripmono,trackdirection);
    const StripTopology& partnertopol=(const StripTopology&)partnerstripdet->topology();

    //error calculation (the part that depends on mono RH only)
    LocalVector  RPHIpositiononGluedendvector=projectedstripmono.second-projectedstripmono.first;
    double c1=sin(RPHIpositiononGluedendvector.phi()); 
    double s1=-cos(RPHIpositiononGluedendvector.phi());
    MeasurementError errormonoRH=topol.measurementError(monoRH->localPosition(),monoRH->localPositionError());
    double pitch=topol.localPitch(monoRH->localPosition());
    double sigmap12=errormonoRH.uu()*pitch*pitch;

    // position of the initial and final point of the strip (STEREO cluster)
    MeasurementPoint STEREOpoint=partnertopol.measurementPosition(stereoRH->localPosition());

    MeasurementPoint STEREOpointini=MeasurementPoint(STEREOpoint.x(),-0.5);
    MeasurementPoint STEREOpointend=MeasurementPoint(STEREOpoint.x(),0.5);

    // position of the initial and final point of the strip in local coordinates (stereo det)
    StripPosition stripstereo(partnertopol.localPosition(STEREOpointini),partnertopol.localPosition(STEREOpointend));

    //project stereo hit on glued det
    StripPosition projectedstripstereo=project(partnerstripdet,gluedDet,stripstereo,trackdirection);

    //perform the matching
    //(x2-x1)(y-y1)=(y2-y1)(x-x1)
    AlgebraicMatrix22 m; AlgebraicVector2 c, solution;
    m(0,0)=-(projectedstripmono.second.y()-projectedstripmono.first.y()); m(0,1)=(projectedstripmono.second.x()-projectedstripmono.first.x());
    m(1,0)=-(projectedstripstereo.second.y()-projectedstripstereo.first.y()); m(1,1)=(projectedstripstereo.second.x()-projectedstripstereo.first.x());
    c(0)=m(0,1)*projectedstripmono.first.y()+m(0,0)*projectedstripmono.first.x();
    c(1)=m(1,1)*projectedstripstereo.first.y()+m(1,0)*projectedstripstereo.first.x();
    m.Invert(); solution = m * c;
    position=LocalPoint(solution(0),solution(1));


    //
    // temporary fix by tommaso
    //


    LocalError tempError (100,0,100);

    // calculate the error
    LocalVector  stereopositiononGluedendvector=projectedstripstereo.second-projectedstripstereo.first;
    double c2=sin(stereopositiononGluedendvector.phi()); double s2=-cos(stereopositiononGluedendvector.phi());
    MeasurementError errorstereoRH=partnertopol.measurementError(stereoRH->localPosition(), stereoRH->localPositionError());
    pitch=partnertopol.localPitch(stereoRH->localPosition());
    double sigmap22=errorstereoRH.uu()*pitch*pitch;
    double diff=(c1*s2-c2*s1);
    double invdet2=1/(diff*diff);
    float xx=invdet2*(sigmap12*s2*s2+sigmap22*s1*s1);
    float xy=-invdet2*(sigmap12*c2*s2+sigmap22*c1*s1);
    float yy=invdet2*(sigmap12*c2*c2+sigmap22*c1*c1);
    LocalError error=LocalError(xx,xy,yy);

    //Added by DAO to make sure y positions are zero.
    DetId det(monoRH->geographicalId());
    if(det.subdetId() > 2) {
	return std::make_unique<FastMatchedTrackerRecHit>(position, error, *gluedDet, *monoRH, *stereoRH,stereoHitFirst);
	
    }
  
    else {
	throw cms::Exception("FastTrackerRecHitMatcher") << "Matched Pixel!?";
    }
}


FastTrackerRecHitMatcher::StripPosition 
    FastTrackerRecHitMatcher::project(const GeomDetUnit *det,
				      const GluedGeomDet* glueddet,
				      const StripPosition& strip,
				      const LocalVector& trackdirection)const
{

    GlobalPoint globalpointini=(det->surface()).toGlobal(strip.first);
    GlobalPoint globalpointend=(det->surface()).toGlobal(strip.second);

    // position of the initial and final point of the strip in glued local coordinates
    LocalPoint positiononGluedini=(glueddet->surface()).toLocal(globalpointini);
    LocalPoint positiononGluedend=(glueddet->surface()).toLocal(globalpointend);

    //correct the position with the track direction

    float scale=-positiononGluedini.z()/trackdirection.z();

    LocalPoint projpositiononGluedini= positiononGluedini + scale*trackdirection;
    LocalPoint projpositiononGluedend= positiononGluedend + scale*trackdirection;

    return StripPosition(projpositiononGluedini,projpositiononGluedend);
}



std::unique_ptr<FastTrackerRecHit> FastTrackerRecHitMatcher::projectOnly( const FastSingleTrackerRecHit *originalRH,
									const GeomDet * monoDet,
									const GluedGeomDet* gluedDet,
									LocalVector& ldir) const
{
    LocalPoint position(originalRH->localPosition().x(), 0.,0.);
    const BoundPlane& gluedPlane = gluedDet->surface();
    const BoundPlane& hitPlane = monoDet->surface();

    double delta = gluedPlane.localZ( hitPlane.position());

    LocalPoint lhitPos = gluedPlane.toLocal( monoDet->surface().toGlobal(position ) );
    LocalPoint projectedHitPos = lhitPos - ldir * delta/ldir.z();

    LocalVector hitXAxis = gluedPlane.toLocal( hitPlane.toGlobal( LocalVector(1,0,0)));
    LocalError hitErr = originalRH->localPositionError();

    if (gluedPlane.normalVector().dot( hitPlane.normalVector()) < 0) {
	// the two planes are inverted, and the correlation element must change sign
	hitErr = LocalError( hitErr.xx(), -hitErr.xy(), hitErr.yy());
    }
    LocalError rotatedError = hitErr.rotate( hitXAxis.x(), hitXAxis.y());
  
  
    const GeomDetUnit *gluedMonoDet = gluedDet->monoDet();
    const GeomDetUnit *gluedStereoDet = gluedDet->stereoDet();
    int isMono = 0;
    int isStereo = 0;
  
    if(monoDet->geographicalId()==gluedMonoDet->geographicalId()) isMono = 1;
    if(monoDet->geographicalId()==gluedStereoDet->geographicalId()) isStereo = 1;
    //Added by DAO to make sure y positions are zero and correct Mono or stereo Det is filled.
  
    if ((isMono && isStereo)||(!isMono&&!isStereo)) throw cms::Exception("FastTrackerRecHitMatcher") << "Something wrong with DetIds.";
    return std::make_unique<FastProjectedTrackerRecHit>(projectedHitPos, rotatedError, *gluedDet, *originalRH);
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void FastTrackerRecHitMatcher::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FastTrackerRecHitMatcher);
