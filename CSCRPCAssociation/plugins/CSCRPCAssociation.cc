// -*- C++ -*-
//
// Package:    CSCRPCAssociation/CSCRPCAssociation
// Class:      CSCRPCAssociation
//
/**\class CSCRPCAssociation CSCRPCAssociation.cc CSCRPCAssociation/CSCRPCAssociation/plugins/CSCRPCAssociation.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mauricio Thiel
//         Created:  Sun, 06 Jun 2021 21:01:33 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCRangeMapAccessor.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHit.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"

#include "Geometry/RPCGeometry/interface/RPCChamber.h"
#include "Geometry/RPCGeometry/interface/RPCRoll.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include "Geometry/CSCGeometry/interface/CSCChamberSpecs.h"
#include "Geometry/CSCGeometry/interface/CSCLayer.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCChamber.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"
#include "Geometry/CommonTopologies/interface/RectangularStripTopology.h"
#include "Geometry/CommonTopologies/interface/TrapezoidalStripTopology.h"

#include "RecoLocalMuon/RPCRecHit/src/CSCStationIndex.h"
#include "RecoLocalMuon/RPCRecHit/src/CSCObjectMap.h"


//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

//using reco::TrackCollection;

class CSCRPCAssociation : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit CSCRPCAssociation(const edm::ParameterSet&);
  ~CSCRPCAssociation();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;


  // Calling the collections
  edm::EDGetTokenT<RPCRecHitCollection> rpc_token;
  edm::EDGetTokenT<CSCSegmentCollection> cscseg_token;
  edm::Handle<RPCRecHitCollection> pRPCrecHits;
  edm::Handle<CSCSegmentCollection> cscsegments;


  // ----------member data ---------------------------
//  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CSCRPCAssociation::CSCRPCAssociation(const edm::ParameterSet& iConfig) // Taking the collections names from the config file
      : rpc_token(consumes<RPCRecHitCollection>(iConfig.getParameter<edm::InputTag>("rpcRecHitTag")))
      , cscseg_token(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("cscSegTag")))
{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

CSCRPCAssociation::~CSCRPCAssociation() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void CSCRPCAssociation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

        using namespace edm;

        edm::ESHandle<RPCGeometry> rpcGeo;
        edm::ESHandle<CSCGeometry> cscGeo;

        iSetup.get<MuonGeometryRecord>().get(rpcGeo);
        iSetup.get<MuonGeometryRecord>().get(cscGeo);

        std::map<CSCStationIndex,std::set<RPCDetId>> rollstore;
        rollstore.clear();
        for (TrackingGeometry::DetContainer::const_iterator it=rpcGeo->dets().begin();it<rpcGeo->dets().end();it++){
                if(dynamic_cast< const RPCChamber* >( *it ) != nullptr ){
                        auto ch = dynamic_cast< const RPCChamber* >( *it );
                        std::vector< const RPCRoll*> roles = (ch->rolls());
                        for(std::vector<const RPCRoll*>::const_iterator r = roles.begin();r != roles.end(); ++r){
                                RPCDetId rpcId = (*r)->id();
                                int region=rpcId.region();
                                if(region!=0){
                                        int station=rpcId.station();
                                        int ring=rpcId.ring();
                                        int cscring=ring;
                                        int cscstation=station;
                                        RPCGeomServ rpcsrv(rpcId);
                                        int rpcsegment = rpcsrv.segment();
                                        int cscchamber = rpcsegment;
                                        if((station==2||station==3)&&ring==3){
                                                cscring = 2;
                                        }
                                        CSCStationIndex ind(region,cscstation,cscring,cscchamber);
                                        std::set<RPCDetId> myrolls;
                                        if (rollstore.find(ind)!=rollstore.end()) myrolls=rollstore[ind];
                                        myrolls.insert(rpcId);
                                        rollstore[ind]=myrolls;
                                }
                        }
                }
        }


        iEvent.getByToken(cscseg_token, cscsegments);
        iEvent.getByToken(rpc_token, pRPCrecHits);

        double MaxD = 80.;
        float eyr = 1; 
        CSCSegmentCollection::const_iterator segment;
        int n_cscsegments = 0;
	if(cscsegments->begin()!=cscsegments->end()) std::cout << ">>>>>>> LIST OF CSC SEGMENTS AND ITS PROPAGATED POINTS IN RPC <<<<<<< " << std::endl;
        for (segment = cscsegments->begin();segment!=cscsegments->end(); ++segment){
                n_cscsegments++;
                CSCDetId CSCId = segment->cscDetId();
                int cscEndCap = CSCId.endcap();
                int cscStation = CSCId.station();
                int cscRing = CSCId.ring();
                int rpcRegion = 1; if(cscEndCap==2) rpcRegion= -1;
                int rpcRing = cscRing;
                if(cscRing==4)rpcRing =1;
                int rpcStation = cscStation;
                int rpcSegment = CSCId.chamber();
                LocalPoint segmentPosition= segment->localPosition();
                LocalVector segmentDirection=segment->localDirection();
                float dz=segmentDirection.z();
                if((segment->dimension()==4) && (segment->nRecHits()<=10 && segment->nRecHits()>=4)){
                        float Xo=segmentPosition.x();
                        float Yo=segmentPosition.y();
                        float Zo=segmentPosition.z();
                        float dx=segmentDirection.x();
                        float dy=segmentDirection.y();
                        float dz=segmentDirection.z();
                        CSCStationIndex theindex(rpcRegion,rpcStation,rpcRing,rpcSegment);
                        static const std::set<RPCDetId> empty;//
                        std::set<RPCDetId> rollsForThisCSC = (rollstore.find(theindex) == rollstore.end()) ? empty : rollstore.at(theindex);//
                        const CSCChamber* TheChamber=cscGeo->chamber(CSCId);
                        if(rpcRing!=1&&rpcStation!=4){
                                assert(!rollsForThisCSC.empty());
                                for (std::set<RPCDetId>::iterator iteraRoll = rollsForThisCSC.begin();iteraRoll != rollsForThisCSC.end(); iteraRoll++){
                                        const RPCRoll* rollasociated = rpcGeo->roll(*iteraRoll);
                                        RPCDetId rpcId = rollasociated->id();
                                        const BoundPlane & RPCSurface = rollasociated->surface();
                                        GlobalPoint CenterPointRollGlobal = RPCSurface.toGlobal(LocalPoint(0,0,0));
                                        LocalPoint CenterRollinCSCFrame = TheChamber->toLocal(CenterPointRollGlobal);
                                        float D=CenterRollinCSCFrame.z();
                                        float X=Xo+dx*D/dz;
                                        float Y=Yo+dy*D/dz;
                                        float Z=D;
                                        const TrapezoidalStripTopology* top_=dynamic_cast<const TrapezoidalStripTopology*>(&(rollasociated->topology()));
                                        LocalPoint xmin = top_->localPosition(0.);
                                        LocalPoint xmax = top_->localPosition((float)rollasociated->nstrips());
                                        float rsize = fabs( xmax.x()-xmin.x() );
                                        float stripl = top_->stripLength();
                                        float stripw = top_->pitch();
                                        float extrapolatedDistance = sqrt((X-Xo)*(X-Xo)+(Y-Yo)*(Y-Yo)+(Z-Zo)*(Z-Zo));
                                        if(extrapolatedDistance<=MaxD){
                                                GlobalPoint GlobalPointExtrapolated=TheChamber->toGlobal(LocalPoint(X,Y,Z));
                                                LocalPoint PointExtrapolatedRPCFrame = RPCSurface.toLocal(GlobalPointExtrapolated);
                                                if(fabs(PointExtrapolatedRPCFrame.z()) < 1. &&
                                                                fabs(PointExtrapolatedRPCFrame.x()) < rsize*eyr &&
                                                                fabs(PointExtrapolatedRPCFrame.y()) < stripl*eyr){
                                                        RPCRecHit RPCPoint(rpcId,0,PointExtrapolatedRPCFrame);
                                                        RPCDetId rpcId_new = RPCPoint.rpcId();
                                                        RPCGeomServ rpcpointsrv_new(rpcId_new);
                                                        std::cout << "CSCSegment in chamber: " << CSCId.chamberName()  << " prapagated to RPC chamber: " << rpcpointsrv_new.chambername() << " with X position: " << RPCPoint.localPosition().x() << std::endl;
                                                }
                                        }
                                }
                        }
                }
        }

        if(pRPCrecHits->begin()!=pRPCrecHits->end()) std::cout << ">>>>>>> LIST OF RPC RECHIT POINTS <<<<<<< " << std::endl;
        RPCRecHitCollection::const_iterator recHit;
        for(recHit = pRPCrecHits->begin(); recHit != pRPCrecHits->end(); recHit++){
                RPCDetId rpcId = recHit->rpcId();
                LocalPoint rpc_position=recHit->localPosition();
                RPCGeomServ rpcsrv(rpcId);
		std::cout << " RPC rechit in chamber: " << rpcsrv.chambername() << "  with local position: " << rpc_position.x() << std::endl;
	}

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
//  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
//  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void CSCRPCAssociation::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void CSCRPCAssociation::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void CSCRPCAssociation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CSCRPCAssociation);
