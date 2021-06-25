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
// Joao
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH2.h"
#include "TH1.h"
#include <set>
#include <map>
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

///Rumi
#include "DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h"
// Joao
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

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

///Rumi
  edm::EDGetTokenT<CSCCorrelatedLCTDigiCollection> cscCorrDigis_token;
  edm::Handle<CSCCorrelatedLCTDigiCollection> cscCorrDigis;

// Joao
  edm::EDGetTokenT<RPCDigiCollection> rpcDigis_token;
  edm::Handle<RPCDigiCollection> rpcDigis;

  int nomatch;
  int MATCH_csc;
  int MATCH_counter=0;
  int NOMATCH_counter=0;
  TH2F *nomatch_rings;
  TH2F *nomatch_rings1;
  TH2F *nomatch_rings2;
  TH2F *nomatch_rings3;
  TH2F *nomatch_rings4;
  TH2F *nomatch_stations;

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
//Rumi
      , cscCorrDigis_token(consumes<CSCCorrelatedLCTDigiCollection>(iConfig.getParameter<edm::InputTag>("cscCorrDigisTag")))
// Joao
      , rpcDigis_token(consumes<RPCDigiCollection>(iConfig.getParameter<edm::InputTag>("rpcDigisTag")))
{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
   edm::Service<TFileService> fs;
   nomatch_rings1 = fs->make<TH2F>("nomatch_rings1" , "Ring of CSC correlated LCT and segments STATION 1;CSC segment ring;CSC correlated LCT ring" , 5 , 0 , 5,5,0,5 );
   nomatch_rings2 = fs->make<TH2F>("nomatch_rings2" , "Ring of CSC correlated LCT and segments STATION 2;CSC segment ring;CSC correlated LCT ring" , 5 , 0 , 5,5,0,5 );
   nomatch_rings3 = fs->make<TH2F>("nomatch_rings3" , "Ring of CSC correlated LCT and segments STATION 3;CSC segment ring;CSC correlated LCT ring" , 5 , 0 , 5,5,0,5 );
   nomatch_rings4 = fs->make<TH2F>("nomatch_rings4" , "Ring of CSC correlated LCT and segments STATION 4;CSC segment ring;CSC correlated LCT ring" , 5 , 0 , 5,5,0,5 );
   nomatch_rings = fs->make<TH2F>("nomatch_rings" , "Ring of CSC correlated LCT and segments;CSC segment ring;CSC correlated LCT ring" , 5 , 0 , 5,5,0,5 );
   nomatch_stations = fs->make<TH2F>("nomatch_stations" , "Station of CSC correlated LCT and segments;CSC segment station;CSC correlated LCT station" , 4 , 1 , 5,4,1,5 );


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
	
	// Joao
	edm::ESHandle<RPCGeometry> rpcGeo;
	edm::ESHandle<CSCGeometry> cscGeo;

	iSetup.get<MuonGeometryRecord>().get(rpcGeo);
	iSetup.get<MuonGeometryRecord>().get(cscGeo);

	iEvent.getByToken(cscCorrDigis_token, cscCorrDigis);
	iEvent.getByToken(rpc_token, pRPCrecHits);
	iEvent.getByToken(cscseg_token, cscsegments);
	iEvent.getByToken(rpcDigis_token, rpcDigis);

	RPCDigiCollection::DigiRangeIterator detUnitIt;
	for (detUnitIt = rpcDigis->begin(); detUnitIt != rpcDigis->end(); ++detUnitIt) {
	    const RPCDetId id = (*detUnitIt).first;
            if (dynamic_cast<const RPCRoll*>(rpcGeo->roll(id)) != 0){//joao insert this if
               const RPCRoll* roll = dynamic_cast<const RPCRoll*>(rpcGeo->roll(id));
               const RPCDigiCollection::Range range = (*detUnitIt).second;
	       // RPCDetId print-out
	       std::cout << "id: " << id.rawId() << " number of strip " << roll->nstrips() << std::endl;
	       // Loop over the digis of this DetUnit
	       // Take the global coordiantes of the center of the rpc roll (eta partition)
	       const int nstrips = roll->nstrips();
	       float middleStrip = nstrips/2.;
	       const LocalPoint& middleOfRoll = roll->centreOfStrip(middleStrip);
	       const GlobalPoint& globMiddleRol = roll->toGlobal(middleOfRoll);
	       std::cout << "eta " << globMiddleRol.eta() << "\tphi " << globMiddleRol.phi() << std::endl;
	       //take the list of fired strips for a given roll
	       int digisInRoll = 0;
	       for (RPCDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; ++digiIt) {
	  	   //std::cout << " digi " << *digiIt << std::endl;
                   if (digiIt->strip() < 1 || digiIt->strip() > roll->nstrips()) {
		      std::cout << " XXXXXXXXXXXXX Problemt with " << id << std::endl;
                   }  
		   int bx=(*digiIt).bx();
	           digisInRoll++;
		   int strip= (*digiIt).strip();
		   std::cout << "strip " <<  strip << "\tbx " <<  bx << std::endl;     
		   //take the global coordinates of the center of a given strip
	           const LocalPoint& middleOfStrip = roll->centreOfStrip(strip);
	           const GlobalPoint& globMiddleStrip = roll->toGlobal(middleOfStrip);
                   std::cout << "eta strip " << globMiddleStrip.eta() << "\tphi strip " << globMiddleStrip.phi() << std::endl;
  
		}	
	        std::cout << "digis in roll " << digisInRoll << std::endl;
	     }
          }

    //std::cout << "CSC INFO" << std::endl;
    int n_cscsegments = 0;
    int n_csccorrelated = 0;
    MATCH_csc=0;
    CSCSegmentCollection::const_iterator segment;
    for (segment = cscsegments->begin();segment!=cscsegments->end(); ++segment){
        n_cscsegments++;
        CSCDetId CSCId = segment->cscDetId();
        int cscEndCap = CSCId.endcap();
        int cscStation = CSCId.station();
        int cscSector = CSCId.triggerSector();
        int cscLayer = CSCId.layer();
        int cscRing = CSCId.ring();
        if (cscRing==4) cscRing=1;
        int cscChamber =  CSCId.chamber();

        MATCH_csc=0;
        for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc = cscCorrDigis.product()->begin();
            csc != cscCorrDigis.product()->end();
            csc++) {
            CSCCorrelatedLCTDigiCollection::Range myCscRange = cscCorrDigis.product()->get((*csc).first);
            for (CSCCorrelatedLCTDigiCollection::const_iterator lct = myCscRange.first; lct != myCscRange.second; lct++) {
                int endcap = (*csc).first.endcap();
                int station = (*csc).first.station();
                int sector = (*csc).first.triggerSector();
                int ring = (*csc).first.ring();
                int cscId = (*csc).first.triggerCscId();
                int bx = lct->getBX();
		int bxData = lct->getBXData();
                int bx0 = lct->getBX0();
                int chamberId =  (*csc).first.chamberId();
                int layer = (*csc).first.layer();
                int chamber = (*csc).first.chamber();
                //if (n_cscsegments==1) std::cout << "CSCCorrLCT: " << "endcap: " << endcap <<"\tstation: " << station <<  "\tsector: " << sector << "\tchamber: " << chamber << "\tring: " << ring <<   "\tbx: " << bx << "\tcscId: " << cscId <<"\tbxData: " << bxData << "\tbx0: " << bx0 << std::endl;
                n_csccorrelated++;
		//std::cout << "CSC segment" << std::endl;
		//if (cscStation==3&&n_csccorrelated ==1)std::cout << "CSCSEGMENT INFO: "<< "cscEndCap: " << cscEndCap << "\tcscStation: " << cscStation << "\tcscSector: " << cscSector << "\tcscChamber: " << cscChamber << "\tcscLayer: " << cscLayer<< "\tcscRing: " << cscRing << std::endl;  
		if (endcap==cscEndCap && station==cscStation && sector==cscSector && chamber== cscChamber &&  cscRing== ring) MATCH_csc = 1;
            }
        }
	if (MATCH_csc==1) {
           //std::cout << "MATCH" << std::endl; 
           MATCH_counter++;
	   //std::cout << MATCH_counter << std::endl;
	}
	if (MATCH_csc==0) {
            NOMATCH_counter++;
	    //std::cout << "NO MATCH" << std::endl;
	    //std::cout << NOMATCH_counter << std::endl;
	    //std::cout << "CSCCORRELATEDLCT INFO: " << "endcap: " << endcap <<"\tstation: " << station <<  "\tsector: " << sector << "\tchamber: " << chamber << "\tlayer: " << layer << "\tring: " << ring <<  "\tchamberID: " << chamberId  << "\tcscId: " << cscId << "\tbx: " << bx << std::endl; 
            for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc = cscCorrDigis.product()->begin();
            csc != cscCorrDigis.product()->end();
            csc++) {
            CSCCorrelatedLCTDigiCollection::Range myCscRange = cscCorrDigis.product()->get((*csc).first);
                for (CSCCorrelatedLCTDigiCollection::const_iterator lct = myCscRange.first; lct != myCscRange.second; lct++) {
                    int endcap = (*csc).first.endcap();
                    int station = (*csc).first.station();
                    int sector = (*csc).first.triggerSector();
                    int ring = (*csc).first.ring();
                    int cscId = (*csc).first.triggerCscId();
                    int bx = lct->getBX();
                    int chamberId =  (*csc).first.chamberId();
                    int layer = (*csc).first.layer();
                    int chamber = (*csc).first.chamber();
                    if (endcap==cscEndCap && station==cscStation && sector==cscSector && chamber== cscChamber && cscLayer==layer &&  cscRing!= ring) {
                       nomatch_rings->Fill(cscRing,ring);
                       if (station==1) nomatch_rings1->Fill(cscRing,ring);
                       if (station==2) nomatch_rings2->Fill(cscRing,ring);
                       if (station==3) nomatch_rings3->Fill(cscRing,ring);
                       if (station==4) nomatch_rings4->Fill(cscRing,ring);
                    }
                    if (endcap==cscEndCap && station!=cscStation && sector==cscSector && chamber== cscChamber && cscLayer==layer &&  cscRing== ring) {nomatch_stations->Fill(cscStation,station);}
                    //std::cout << "CSCSEGMENT INFO: "<< "cscEndCap: " << cscEndCap << "\tcscStation: " << cscStation << "\tcscSector: " << cscSector << "\tcscChamber: " << cscChamber << "\tcscLayer: " << cscLayer<< "\tcscRing: " << cscRing << std::endl; 
                }       
            }
        }
     }
     //std::cout << "RPC INFO" << std::endl;
     //if(pRPCrecHits->begin()!=pRPCrecHits->end()) std::cout << ">>>>>>> LIST OF RPC RECHIT POINTS <<<<<<< " << std::endl;
     RPCRecHitCollection::const_iterator recHit;
     for(recHit = pRPCrecHits->begin(); recHit != pRPCrecHits->end(); recHit++){
        RPCDetId rpcId = recHit->rpcId();
	//LocalPoint rpc_position=recHit->localPosition();
        RPCGeomServ rpcsrv(rpcId);
        int stationRPC =  rpcId.station();
        int regionRPC  =  rpcId.region();
        if (regionRPC==0) continue;
        int ringRPC    =  rpcId.ring();
        int sectorRPC  =  rpcId.sector();
        int layerRPC   =  rpcId.layer();
        int BXRPC      =  recHit->BunchX();
        int subsectorRPC = rpcId.subsector();
        int rollRPC = rpcId.roll();
        int segmentRPC =  rpcsrv.segment();
        // Skip iRPCs sectors
	if (stationRPC==3 && ringRPC==1) continue;
        if (stationRPC==4 && ringRPC==1) continue;
        for (CSCCorrelatedLCTDigiCollection::DigiRangeIterator csc = cscCorrDigis.product()->begin();
            csc != cscCorrDigis.product()->end();
            csc++) {
            CSCCorrelatedLCTDigiCollection::Range myCscRange = cscCorrDigis.product()->get((*csc).first);
            for (CSCCorrelatedLCTDigiCollection::const_iterator lct = myCscRange.first; lct != myCscRange.second; lct++) {
                int endcapCSC = (*csc).first.endcap();
                int stationCSC = (*csc).first.station();
                int sectorCSC = (*csc).first.triggerSector();
                int ringCSC = (*csc).first.ring();
                int cscIdCSC = (*csc).first.triggerCscId();
                int bxCSC = lct->getBX();
                int bxDataCSC = lct->getBXData();
                int bx0CSC = lct->getBX0();
                int chamberIdCSC =  (*csc).first.chamberId();
                int layerCSC = (*csc).first.layer();
                int chamberCSC = (*csc).first.chamber();
		if (endcapCSC==2) endcapCSC=-1; //for RPC consistency
		if (endcapCSC==regionRPC && stationCSC==stationRPC && chamberCSC==segmentRPC){
		   //std::cout << "CSCCorrLCT: " << "endcap: " << endcapCSC <<"\tstation: " << stationCSC <<  "\tsector: " << sectorCSC << "\tchamber: " << chamberCSC << "\tring: " << ringCSC <<   "\tbx: " << bxCSC << "\tcscId: " << cscIdCSC <<"\tbxData: " << bxDataCSC << "\tbx0: " << bx0CSC << std::endl;
		   //std::cout <<  "RPCrecHit: " << " region: " << regionRPC << "\tstation: " << stationRPC <<  "\tsector: " << sectorRPC << "\tchamber: " << segmentRPC <<"\tring: " << ringRPC << "\tBX: " << BXRPC   << "\tsubsector: " << subsectorRPC << " chambername: " << rpcsrv.chambername() << " roll: " << rollRPC << std::endl;
		   //std::cout << " " << std::endl;		   
		}

	    }
	}

	//std::cout <<  "RPCrecHit: " << " region: " << region << "\tstation: " << station <<  "\tsector: " << sector << "\tchamber: " << segment <<"\tring: " << ring << "\tBX: " << BX   << "\tsubsector: " << subsector << " chambername: " << rpcsrv.chambername() << " roll: " << roll << std::endl;
     }
     std::map<CSCStationIndex,std::set<RPCDetId>> rollstore;
     rollstore.clear();
     //std::cout << ">>>>>>> LIST OF RPC GENERAL INFORMATION <<<<<<< " << std::endl;
     for (TrackingGeometry::DetContainer::const_iterator it=rpcGeo->dets().begin();it<rpcGeo->dets().end();it++){
         if (dynamic_cast< const RPCChamber* >( *it ) != nullptr ){
            auto ch = dynamic_cast< const RPCChamber* >( *it );
            std::vector< const RPCRoll*> roles = (ch->rolls());
            for (std::vector<const RPCRoll*>::const_iterator r = roles.begin();r != roles.end(); ++r){
                RPCDetId rpcId = (*r)->id();
                int rpcregion=rpcId.region();
                if (rpcregion!=0){
                   int rpcstation=rpcId.station();
                   int rpcring=rpcId.ring();
                   int cscring=rpcring;
                   int cscstation=rpcstation;
                   RPCGeomServ rpcsrv(rpcId);
                   int rpcsegment = rpcsrv.segment();
                   int cscchamber = rpcsegment;
                   if ((rpcstation==2||rpcstation==3)&&rpcring==3){
                      cscring = 2;
                   }
                   //if (station==3) std::cout << "region: " << region << "\tstation: " << station << "\tring: " << ring << "\trpcsegment: " << rpcsegment << std::endl;
                   CSCStationIndex ind(rpcregion,cscstation,cscring,cscchamber);
                   std::set<RPCDetId> myrolls;
                   if (rollstore.find(ind)!=rollstore.end()) myrolls=rollstore[ind];
                   myrolls.insert(rpcId);
                   rollstore[ind]=myrolls;
                }
             }
          }
     }

     double MaxD = 80.;
     float eyr = 1;
     //int n_cscsegments = 0;
     //if(cscsegments->begin()!=cscsegments->end()) std::cout << ">>>>>>> LIST OF CSC SEGMENTS AND ITS PROPAGATED POINTS IN RPC <<<<<<< " << std::endl;
     //std::cout << "CSC segment" << std::endl;

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
         // float dz=segmentDirection.z();
         //if((segment->dimension()==4) && (segment->nRecHits()<=10 && segment->nRecHits()>=4)){
         float Xo=segmentPosition.x();
         float Yo=segmentPosition.y();
         float Zo=segmentPosition.z();
         float dx=segmentDirection.x();
         float dy=segmentDirection.y();
         float dz=segmentDirection.z();
         CSCStationIndex theindex(rpcRegion,rpcStation,rpcRing,rpcSegment);
         static const std::set<RPCDetId> empty;
         //std::set<RPCDetId> rollsForThisCSC = (rollstore.find(theindex) == rollstore.end()) ? empty : rollstore.at(theindex);//
         const CSCChamber* TheChamber=cscGeo->chamber(CSCId);
	 //if (cscStation==3) std::cout << "cscEndCap: " << cscEndCap << "\tcscStation: " << cscStation << "\tcscChamber: " << rpcSegment << "\tcscLayer: " << CSCId.layer()<< "\tcscRing: " << cscRing << std::endl;
	 /*if (rpcRing!=1&&rpcStation!=4){
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
                  if (extrapolatedDistance<=MaxD){
                     GlobalPoint GlobalPointExtrapolated=TheChamber->toGlobal(LocalPoint(X,Y,Z));
                     LocalPoint PointExtrapolatedRPCFrame = RPCSurface.toLocal(GlobalPointExtrapolated);
                     if(fabs(PointExtrapolatedRPCFrame.z()) < 1. &&
                        fabs(PointExtrapolatedRPCFrame.x()) < rsize*eyr &&
                        fabs(PointExtrapolatedRPCFrame.y()) < stripl*eyr){
                        RPCRecHit RPCPoint(rpcId,0,PointExtrapolatedRPCFrame);
                        RPCDetId rpcId_new = RPCPoint.rpcId();
                        RPCGeomServ rpcpointsrv_new(rpcId_new);
			//std::cout << "CSCSegment in chamber: " << CSCId.chamberName()  << " propagated to RPC chamber: " << rpcpointsrv_new.chambername() << " with X position: " << RPCPoint.localPosition().x() << std::endl;
		     }
		  }
	      }
	   }
	   //}
	*/
     }
     /*
     if(pRPCrecHits->begin()!=pRPCrecHits->end()) std::cout << ">>>>>>> LIST OF RPC RECHIT POINTS <<<<<<< " << std::endl;
     RPCRecHitCollection::const_iterator recHit;
     for (recHit = pRPCrecHits->begin(); recHit != pRPCrecHits->end(); recHit++){
          RPCDetId rpcId = recHit->rpcId();
          LocalPoint rpc_position=recHit->localPosition();
          RPCGeomServ rpcsrv(rpcId);
          std::cout << " RPC rechit in chamber: " << rpcsrv.chambername() << "  with local position: " << rpc_position.x() << std::endl;
     }
    */


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
