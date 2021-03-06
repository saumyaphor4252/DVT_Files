// -*- C++ -*-
//
// Package:    DatasetValidation/DatasetValidationTool
// Class:      DatasetValidationTool_Tree
// 
/**\class DatasetValidationTool_Tree DatasetValidationTool_Tree.cc DatasetValidation/DatasetValidationTool/plugins/DatasetValidationTool_Tree.cc

 Description: Tool to validate datasets used for CMS Tracker Alignment by producing standard distributions

 Implementation:
     Creates tuple with track variables required to validate the datasets 
*/
//
// Original Author:  Saumya Saumya
//         Created:  Wed, 02 Sep 2020 11:05:57 GMT
//
//


// system include files
#include <memory>
/*
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Framework/interface/EventSetup.h"
*/
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "CondFormats/DataRecord/interface/RunSummaryRcd.h"
#include "CondFormats/DataRecord/interface/SiStripCondDataRecords.h"
#include "CondFormats/RunInfo/interface/RunInfo.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/Phi.h"
#include "DataFormats/GeometryVector/interface/Theta.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Provenance/interface/RunLumiEventNumber.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackResiduals.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include <TTree.h>
#include <TMath.h>
#include <TLorentzVector.h>
//
// class declaration
//
class DatasetValidationTool_Tree: public edm::one::EDAnalyzer<edm::one::WatchRuns, edm::one::WatchLuminosityBlocks, edm::one::SharedResources>
{
   public:
      explicit DatasetValidationTool_Tree(const edm::ParameterSet&);
      ~DatasetValidationTool_Tree();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&) override; 
      virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::TrackCollection> tracksToken_;
      edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
      edm::EDGetTokenT<reco::VertexCollection> vertexToken_;

      const bool isHit2D(const TrackingRecHit &hit);
 
      TTree *treeEvent, *treeRun, *treeLuminosity;
      bool isResonance;
           
      int nTracks,nEvents,nTracksInEvent,nEventsInRun,nTracksInRun,nTracksInLuminosity,nEventsInLuminosity,InvalidHit,nNotDetUnits,nZeroSubDetId;
      int nHits_2D,nHits_PIXEL,nHits_FPIXplus,nHits_FPIXminus,nHits_TIDplus,nHits_TIDminus,nHits_TECplus,nHits_TECminus,nHits_ENDCAP,nHits_ENDCAPplus,nHits_ENDCAPminus;
     
      //-----Kinematic Variables-----// 
      std::vector<double> charge;
      std::vector<double> p;
      std::vector<double> pt;
      std::vector<double> eta;
      std::vector<double> theta;
      std::vector<double> phi;
      std::vector<double> chi2;
      std::vector<double> chi2_ndf;
      std::vector<double> chi2_Prob;
      std::vector<double> d0;
      std::vector<double> dz;
      std::vector<double> dxy;
      std::vector<double> d0PV;
      std::vector<double> dzPV;
      std::vector<double> dxyPV;
      std::vector<double> d0BS;
      std::vector<double> dzBS;
      std::vector<double> dxyBS;
      //-----Hits-----//
      std::vector<int> nh_Valid;
      std::vector<int> nh_2D;
      std::vector<int> nh_PIXEL;
      std::vector<int> nh_BPIX;
      std::vector<int> nh_FPIX;
      std::vector<int> nh_FPIXplus;
      std::vector<int> nh_FPIXminus;
      std::vector<int> nh_TIB;
      std::vector<int> nh_TOB;
      std::vector<int> nh_TID;
      std::vector<int> nh_TIDplus;
      std::vector<int> nh_TIDminus;
      std::vector<int> nh_TEC;
      std::vector<int> nh_TECplus;
      std::vector<int> nh_TECminus;
      std::vector<int> nh_ENDCAP;
      std::vector<int> nh_ENDCAPplus;
      std::vector<int> nh_ENDCAPminus;
      //-----Residuals-----//
      std::vector<double> Res_BPIX_xPrime; 
      std::vector<double> Res_FPIX_xPrime;
//      std::vector<double> Res_FPIXplus_xPrime;
//      std::vector<double> Res_FPIXminus_xPrime;
      std::vector<double> Res_TIB_xPrime;
      std::vector<double> Res_TID_xPrime;
      std::vector<double> Res_TOB_xPrime;
      std::vector<double> Res_TEC_xPrime;
      std::vector<double> Res_BPIX_yPrime;
      std::vector<double> Res_FPIX_yPrime;
//      std::vector<double> Res_FPIXplus_yPrime;
//      std::vector<double> Res_FPIXminus_yPrime;
      std::vector<double> Res_TIB_yPrime;
      std::vector<double> Res_TID_yPrime;
      std::vector<double> Res_TOB_yPrime;
      std::vector<double> Res_TEC_yPrime;
 
     //----- Resonances------//
     std::vector<double> Resonance;
     std::vector<double> Resonance_Eta;
     std::vector<double> Resonance_Pt;

     //-----Multiplicity-----//
     std::vector<int> Tracks_In_Event;
     std::vector<int> Events_In_Run;
     std::vector<int> Tracks_In_Run;
     std::vector<int> Events_In_Luminosity;
     std::vector<int> Tracks_In_Luminosity;
     
     std::vector<double> Temp;
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
DatasetValidationTool_Tree::DatasetValidationTool_Tree(const edm::ParameterSet& iConfig)
: tracksToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"))),
  beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("BS"))),
  vertexToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
{
   //now do what ever initialization is needed
   usesResource("TFileService");
   edm::Service<TFileService> fs;
   treeEvent = fs->make<TTree>("Event", "");
   treeRun = fs->make<TTree>("Run","");
   treeLuminosity = fs->make<TTree>("Luminosity","");
   isResonance = iConfig.getParameter<bool>("IsResonance");

   nTracks=0;nEvents=0,InvalidHit=0,nNotDetUnits=0,nZeroSubDetId=0;
}


DatasetValidationTool_Tree::~DatasetValidationTool_Tree()
{
}

//
// member functions
//

// ------------ method called for each event  ------------
void
DatasetValidationTool_Tree::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   nTracksInEvent=0; 
   charge.clear();
   p.clear();
   pt.clear();
   eta.clear();
   theta.clear();
   phi.clear();
   chi2.clear();
   chi2_ndf.clear();
   chi2_Prob.clear();
   d0.clear();
   dz.clear();
   dxy.clear();
   d0PV.clear();
   dxyPV.clear();
   dzPV.clear();
   dxyBS.clear();
   d0BS.clear();
   dzBS.clear();

   nh_Valid.clear();
   nh_2D.clear();
   nh_PIXEL.clear();
   nh_BPIX.clear();
   nh_FPIX.clear();
   nh_FPIXplus.clear();
   nh_FPIXminus.clear();
   nh_TIB.clear();
   nh_TOB.clear();
   nh_TID.clear();
   nh_TIDplus.clear();
   nh_TIDminus.clear();
   nh_TEC.clear();
   nh_TECplus.clear();
   nh_TECminus.clear();
   nh_ENDCAP.clear();
   nh_ENDCAPplus.clear();
   nh_ENDCAPminus.clear();

   Resonance.clear();
   Resonance_Eta.clear();
   Resonance_Pt.clear();

   Tracks_In_Event.clear();

   using namespace edm;

   edm::Handle<reco::TrackCollection> tracksHandle_;
   iEvent.getByToken(tracksToken_,tracksHandle_);
   const reco::TrackCollection tC = *(tracksHandle_.product());   
 
   // Geometry setup
   edm::ESHandle<TrackerGeometry> geometry;
   iSetup.get<TrackerDigiGeometryRecord>().get(geometry);
   const TrackerGeometry *theGeometry = &(*geometry);

   //Retrieve tracker topology from geometry
   edm::ESHandle<TrackerTopology> tTopoHandle;
   iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
   const TrackerTopology *const tTopo = tTopoHandle.product();

   //Beamspot Handle for dxy(BS)
   reco::BeamSpot beamSpot;
   edm::Handle<reco::BeamSpot> beamSpotHandle_;
   iEvent.getByToken(beamspotToken_, beamSpotHandle_);

   //Vertex Handle for dxy(PV)
   edm::Handle<reco::VertexCollection> vertexHandle_;
   iEvent.getByToken(vertexToken_,vertexHandle_);
    
   for(reco::TrackCollection::const_iterator track=tC.begin(); track!=tC.end(); track++)
   {  
     nTracks++; nTracksInEvent++; nTracksInRun++; nTracksInLuminosity++;

     Res_BPIX_xPrime.clear();
     Res_FPIX_xPrime.clear();
     //Res_FPIXplus_xPrime.clear();
     //Res_FPIXminus_xPrime.clear();
     Res_TIB_xPrime.clear();
     Res_TID_xPrime.clear();
     Res_TOB_xPrime.clear();
     Res_TEC_xPrime.clear();
     Temp.clear();
     Res_BPIX_yPrime.clear();
     Res_FPIX_yPrime.clear();
//   Res_FPIXplus_yPrime.clear();
//   Res_FPIXminus_yPrime.clear();
     Res_TIB_yPrime.clear();
     Res_TID_yPrime.clear();
     Res_TOB_yPrime.clear();
     Res_TEC_yPrime.clear();
           
     charge.push_back(track->charge());
     p.push_back(track->p());
     pt.push_back(track->pt());
     eta.push_back(track->eta());
     theta.push_back(track->theta());
     phi.push_back(track->phi());
     chi2.push_back(track->chi2());
     chi2_ndf.push_back(track->normalizedChi2());
     chi2_Prob.push_back(TMath::Prob(track->chi2(), track->ndof()));
     d0.push_back(track->d0());
     dz.push_back(track->dz());

     if (beamSpotHandle_.isValid())
     {
        beamSpot = *beamSpotHandle_;
        math::XYZPoint BSpoint(beamSpot.x0(), beamSpot.y0(), beamSpot.z0());
        double dxy = track->dxy(BSpoint);
        double dz = track->dz(BSpoint);
        dxyBS.push_back(dxy);
        d0BS.push_back(-dxy);
        dzBS.push_back(dz);
     }

     if(vertexHandle_.isValid())
     {
        reco::Vertex vtx;
        double min_dxy=100., dz=100.;
        for(auto vtx=vertexHandle_->cbegin();vtx!=vertexHandle_->cend();++vtx)
        {
            math::XYZPoint Vpoint(vtx->x(), vtx->y(), vtx->z());
            if(abs(min_dxy) > abs(track->dxy(Vpoint)))
            {
               min_dxy = track->dxy(Vpoint);
               dz = track->dz(Vpoint);
            }
        }
        dxyPV.push_back(min_dxy);
        d0PV.push_back(-min_dxy);
        dzPV.push_back(dz);
     }

     //---------------------------Resonances-----------------------------//
     if(isResonance)
     {
         double InvMass=0;
         for (reco::TrackCollection::const_iterator trk2=track+1; trk2!=tC.end(); trk2++)
         {   
             if(track->charge()*trk2->charge()>0) continue;
             TLorentzVector track1,track2,mother;
             track1.SetPtEtaPhiE(track->pt(),track->eta(),track->phi(),sqrt((track->p()*track->p())+(0.105*0.105)));
             track2.SetPtEtaPhiE(trk2->pt(),trk2->eta(),trk2->phi(),sqrt((trk2->p()*trk2->p())+(0.105*0.105)));
             mother=track1+track2; InvMass=mother.M();
             if(InvMass<0. || InvMass>5.) continue;
             Resonance.push_back(InvMass);
             Resonance_Eta.push_back(mother.Eta());
             Resonance_Pt.push_back(mother.Pt());             
         }
     }

     //-------------------- Hits and Residuals --------------------------//
     nh_Valid.push_back(track->numberOfValidHits());
     nh_BPIX.push_back(track->hitPattern().numberOfValidPixelBarrelHits());
     nh_FPIX.push_back(track->hitPattern().numberOfValidPixelEndcapHits());
     nh_TIB.push_back(track->hitPattern().numberOfValidStripTIBHits());
     nh_TOB.push_back(track->hitPattern().numberOfValidStripTOBHits());
     nh_TID.push_back(track->hitPattern().numberOfValidStripTIDHits());
     nh_TEC.push_back(track->hitPattern().numberOfValidStripTECHits());

     nHits_2D=0;
     nHits_PIXEL=0;
     nHits_FPIXplus=0;
     nHits_FPIXminus=0;
     nHits_TIDplus=0;
     nHits_TIDminus=0;
     nHits_TECplus=0;
     nHits_TECminus=0;
     nHits_ENDCAP=0;
     nHits_ENDCAPplus=0;
     nHits_ENDCAPminus=0;

     auto const &trajParams = track->extra()->trajParams();
     auto const &residuals = track->extra()->residuals();
     assert(trajParams.size() == track->recHitsSize());

     int h_index = 0;
     for(auto iHit = track->recHitsBegin(); iHit!=track->recHitsEnd(); ++iHit,++h_index)
     {   
                  
         const DetId detId=(*iHit)->geographicalId();
         const int SubDetId = detId.subdetId();       
         const GeomDet *geomDet(theGeometry->idToDet(detId));

         if (!(*iHit)->isValid()) { InvalidHit++; continue; } 
         if (SubDetId == 0) { nZeroSubDetId++; continue;} 
         if (!(*iHit)->detUnit())  {  nNotDetUnits++; continue;  }// is it a single physical module?

         if (isHit2D(**iHit)) { ++nHits_2D; }  // 2D Hit

         float uOrientation(-999.F), vOrientation(-999.F);
        
         double resX = residuals.residualX(h_index); 
         double resY = residuals.residualY(h_index);

         // All the transformations here      
         LocalPoint lPModule(0., 0., 0.), lUDirection(1., 0., 0.), lVDirection(0., 1., 0.);
         GlobalPoint gUDirection = geomDet->surface().toGlobal(lUDirection);
         GlobalPoint gVDirection = geomDet->surface().toGlobal(lVDirection);
         GlobalPoint gPModule = geomDet->surface().toGlobal(lPModule);          
        
 
         //-------------------------------------- Filling Hits and Residuals -------------------------------------------//

         if (SubDetId == PixelSubdetector::PixelBarrel || SubDetId == StripSubdetector::TIB || SubDetId == StripSubdetector::TOB) 
         {
            uOrientation = deltaPhi(gUDirection.barePhi(), gPModule.barePhi()) >= 0. ? +1.F : -1.F;
            vOrientation = gVDirection.z() - gPModule.z() >= 0 ? +1.F : -1.F;
            switch(SubDetId)
            {
               case PixelSubdetector::PixelBarrel:
               Res_BPIX_xPrime.push_back(uOrientation*resX*10000);
               Res_BPIX_yPrime.push_back(vOrientation*resY*10000);
               nHits_PIXEL++; 
               break;

               case StripSubdetector::TIB:
               Res_TIB_xPrime.push_back(uOrientation*resX*10000);
               Temp.push_back(uOrientation*resX*10000);
               Res_TIB_yPrime.push_back(vOrientation*resY*10000);
               break;

               case StripSubdetector::TOB:
               Res_TOB_xPrime.push_back(uOrientation*resX*10000);
               Res_TOB_yPrime.push_back(vOrientation*resY*10000);
               break;
            }
        }
        else if ( SubDetId == PixelSubdetector::PixelEndcap || SubDetId == StripSubdetector::TID || SubDetId == StripSubdetector::TEC)     
        {
           uOrientation = deltaPhi(gUDirection.barePhi(), gPModule.barePhi()) >= 0. ? +1.F : -1.F;
           vOrientation = gVDirection.perp() - gPModule.perp() >= 0. ? +1.F : -1.F;
           switch(SubDetId)
           {
              case PixelSubdetector::PixelEndcap:
              Res_FPIX_xPrime.push_back(uOrientation*resX*10000);
              Res_FPIX_yPrime.push_back(vOrientation*resY*10000);
              nHits_PIXEL++; nHits_ENDCAP++;
              if(tTopo->pxfSide(detId)==1) {  nHits_FPIXminus++; nHits_ENDCAPminus++;}
              else {  nHits_FPIXplus++; nHits_ENDCAPplus++; }
              break;

              case StripSubdetector::TID:
              Res_TID_xPrime.push_back(uOrientation*resX*10000);
              Res_TID_yPrime.push_back(vOrientation*resY*10000);
              nHits_ENDCAP++;
              if(tTopo->tidIsZMinusSide(detId)) {  nHits_TIDminus++; nHits_ENDCAPminus++; }
              else {  nHits_TIDplus++; nHits_ENDCAPplus++; }
              break;

              case StripSubdetector::TEC:
              Res_TEC_xPrime.push_back(uOrientation*resX*10000);
              Res_TEC_yPrime.push_back(vOrientation*resY*10000);
              nHits_ENDCAP++;
              if(tTopo->tecIsZMinusSide(detId)) {  nHits_TECminus++; nHits_ENDCAPminus++; }
              else {  nHits_TECplus++; nHits_ENDCAPplus++; }
              break;
           }
        } 

     }  //Hits Loop

     nh_2D.push_back(nHits_2D);
     nh_PIXEL.push_back(nHits_PIXEL);
     nh_FPIXplus.push_back(nHits_FPIXplus);
     nh_FPIXminus.push_back(nHits_FPIXminus);
     nh_TIDplus.push_back(nHits_TIDplus);
     nh_TIDminus.push_back(nHits_TIDminus);
     nh_TECplus.push_back(nHits_TECplus);
     nh_TECminus.push_back(nHits_TECminus);
     nh_ENDCAP.push_back(nHits_ENDCAP);
     nh_ENDCAPplus.push_back(nHits_ENDCAPplus);
     nh_ENDCAPminus.push_back(nHits_ENDCAPminus);

   } //Tracks Loop

   nEvents++;  nEventsInRun++;  nEventsInLuminosity++;
   Tracks_In_Event.push_back(nTracksInEvent);
   treeEvent->Fill();


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
DatasetValidationTool_Tree::beginJob()
{
   treeEvent->Branch("charge", &charge);
   treeEvent->Branch("p", &p);
   treeEvent->Branch("pt", &pt);
   treeEvent->Branch("eta", &eta);
   treeEvent->Branch("theta", &theta);
   treeEvent->Branch("phi", &phi);
   treeEvent->Branch("chi2", &chi2);
   treeEvent->Branch("chi2_ndf", &chi2_ndf);
   treeEvent->Branch("chi2_Prob", &chi2_Prob);
   treeEvent->Branch("d0", &d0);
   treeEvent->Branch("dz", &dz);
   treeEvent->Branch("d0PV", &d0PV);
   treeEvent->Branch("dxyPV", &dxyPV);
   treeEvent->Branch("dzPV", &dzPV);
   treeEvent->Branch("d0BS", &d0BS);
   treeEvent->Branch("dxyBS", &dxyBS);
   treeEvent->Branch("dzBS", &dzBS);

   treeEvent->Branch("nHitsValid", &nh_Valid);
   treeEvent->Branch("nHits2D", &nh_2D);
   treeEvent->Branch("nHitsPIXEL", &nh_PIXEL);
   treeEvent->Branch("nHitsBPIX", &nh_BPIX);
   treeEvent->Branch("nHitsFPIX", &nh_FPIX);
   treeEvent->Branch("nHitsFPIXplus", &nh_FPIXplus);
   treeEvent->Branch("nHitsFPIXminus", &nh_FPIXminus);
   treeEvent->Branch("nHitsTIB", &nh_TIB);
   treeEvent->Branch("nHitsTOB", &nh_TOB);
   treeEvent->Branch("nHitsTID", &nh_TID);
   treeEvent->Branch("nHitsTIDplus", &nh_TIDplus);
   treeEvent->Branch("nHitsTIDminus", &nh_TIDminus);
   treeEvent->Branch("nHitsTEC", &nh_TEC);
   treeEvent->Branch("nHitsTECplus", &nh_TECplus);
   treeEvent->Branch("nHitsTECminus", &nh_TECminus);
   treeEvent->Branch("nHitsENDCAP", &nh_ENDCAP);
   treeEvent->Branch("nHitsENDCAPplus", &nh_ENDCAPplus);
   treeEvent->Branch("nHitsENDCAPminus", &nh_ENDCAPminus);

   treeEvent->Branch("ResBPIXxPrime",&Res_BPIX_xPrime);
   treeEvent->Branch("ResFPIXxPrime", &Res_FPIX_xPrime);
//   treeEvent->Branch("ResFPIXplusxPrime", &Res_FPIXplus_xPrime);
//   treeEvent->Branch("ResFPIXminusxPrime", &Res_FPIXminus_xPrime);
   treeEvent->Branch("ResTIBxPrime", &Res_TIB_xPrime);
   treeEvent->Branch("ResTIDxPrime", &Res_TID_xPrime);
   treeEvent->Branch("ResTOBxPrime", &Res_TOB_xPrime);
   treeEvent->Branch("ResTECxPrime", &Res_TEC_xPrime);
   treeEvent->Branch("ResBPIXyPrime",&Res_BPIX_yPrime);
   treeEvent->Branch("ResFPIXyPrime", &Res_FPIX_yPrime);
//   treeEvent->Branch("ResFPIXplusyPrime", &Res_FPIXplus_yPrime);
//   treeEvent->Branch("ResFPIXminusyPrime", &Res_FPIXminus_yPrime);
   treeEvent->Branch("ResTIByPrime", &Res_TIB_yPrime);
   treeEvent->Branch("ResTIDyPrime", &Res_TID_yPrime);
   treeEvent->Branch("ResTOByPrime", &Res_TOB_yPrime);
   treeEvent->Branch("ResTECyPrime", &Res_TEC_yPrime);

   treeEvent->Branch("ResonanceMass", &Resonance);
   treeEvent->Branch("ResonanceEta", &Resonance_Eta);
   treeEvent->Branch("ResonancePt", &Resonance_Pt);
    
   treeEvent->Branch("TracksPerEvent",&Tracks_In_Event);
   treeRun->Branch("EventsPerRun",&Events_In_Run);
   treeRun->Branch("TracksPerRun",&Tracks_In_Run);
   treeLuminosity->Branch("EventsPerLuminosity",&Events_In_Luminosity);
   treeLuminosity->Branch("TracksPerLuminosity",&Tracks_In_Luminosity);

   treeEvent->Branch("Temp",&Temp);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DatasetValidationTool_Tree::endJob() 
{
   std::cout<<"Events: "<<nEvents<<std::endl;
   std::cout<<"Tracks: "<<nTracks<<std::endl;
   std::cout<<"Invalid Hits: "<<InvalidHit<<std::endl;
   std::cout<<"SubDet0 HIts: "<<nZeroSubDetId<<std::endl;
   std::cout<<"nNotDetUnits: "<<nNotDetUnits<<std::endl;
 }

// ------------ method called when starting to processes a run  ------------
void
DatasetValidationTool_Tree::beginRun(edm::Run const&, edm::EventSetup const&)
{   
   nEventsInRun=0; nTracksInRun=0;
   Tracks_In_Run.clear();
   Events_In_Run.clear();
}

// ------------ method called when ending the processing of a run  ------------
void
DatasetValidationTool_Tree::endRun(edm::Run const&, edm::EventSetup const&)
{  
    Events_In_Run.push_back(nEventsInRun); 
    Tracks_In_Run.push_back(nTracksInRun);
    treeRun->Fill();
} 

// ------------ method called when starting to processes a luminosity block  ------------
void
DatasetValidationTool_Tree::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{ 
    nEventsInLuminosity=0; nTracksInLuminosity=0;
    Tracks_In_Luminosity.clear();
    Events_In_Luminosity.clear();
}

// ------------ method called when ending the processing of a luminosity block  ------------
void
DatasetValidationTool_Tree::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
    Events_In_Luminosity.push_back(nEventsInLuminosity);
    Tracks_In_Luminosity.push_back(nTracksInLuminosity);
    treeLuminosity->Fill();
}

// ------------------------------ 2D Hit Function ---------------------------------------
// WILL BREAK DOWN FOR PHASE 2 !!!!!!!!!!!!!!!!!!   
const bool DatasetValidationTool_Tree::isHit2D(const TrackingRecHit &hit)
{
     bool countStereoHitAs2D_ = true;
     // we count SiStrip stereo modules as 2D if selected via countStereoHitAs2D_
     // (since they provide theta information)
     if (!hit.isValid() || (hit.dimension() < 2 && !countStereoHitAs2D_ && !dynamic_cast<const SiStripRecHit1D *>(&hit)))
     {  return false; // real RecHit1D - but SiStripRecHit1D depends on countStereoHitAs2D_  
     }
     else 
     {
        const DetId detId(hit.geographicalId());
        if (detId.det() == DetId::Tracker)
        {
           if (detId.subdetId() == PixelSubdetector::PixelBarrel || detId.subdetId() == PixelSubdetector::PixelEndcap)
           {  return true;  // pixel is always 2D
           }
           else
           {  // should be SiStrip now
              const SiStripDetId stripId(detId);
              if(stripId.stereo())  return countStereoHitAs2D_;  // stereo modules
              else if (dynamic_cast<const SiStripRecHit1D *>(&hit) || dynamic_cast<const SiStripRecHit2D *>(&hit)) return false;  // rphi modules hit
              //the following two are not used any more since ages...
              else if (dynamic_cast<const SiStripMatchedRecHit2D *>(&hit)) return true;  // matched is 2D       
              else if (dynamic_cast<const ProjectedSiStripRecHit2D *>(&hit))
              {
                 const ProjectedSiStripRecHit2D *pH = static_cast<const ProjectedSiStripRecHit2D *>(&hit);
                 return (countStereoHitAs2D_ && isHit2D(pH->originalHit()));  // depends on original...
              }
              else 
              { // edm::LogError("UnkownType") << "@SUB=DMRChecker::isHit2D"
                //                         << "Tracker hit not in pixel, neither SiStripRecHit[12]D nor "
                //                         << "SiStripMatchedRecHit2D nor ProjectedSiStripRecHit2D.";
                 return false;
              }
           } 
        }
        else 
        {   // not tracker??
            // edm::LogWarning("DetectorMismatch") << "@SUB=DMRChecker::isHit2D"
            //                                  << "Hit not in tracker with 'official' dimension >=2.";
            return true;  // dimension() >= 2 so accept that...
        }
    }// never reached...
}


// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DatasetValidationTool_Tree::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters

 edm::ParameterSetDescription desc;
 desc.setUnknown();
 descriptions.addDefault(desc);

 // edm::ParameterSetDescription desc;
  //desc.setComment("Create tuple with all variables required to workk with");
 // desc.add<edm::InputTag> ("tracksToken_",edm::InputTag("ALCARECOTkAlCosmicsCTF0T"));
 // descriptions.add("datasetValidationTool",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DatasetValidationTool_Tree);
