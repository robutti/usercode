// -*- C++ -*-
//
// Package:    HoughCheckXY
// Class:      HoughCheckXYOnTracks
// 
/**\class HoughCheckXYOnTracks HoughCheckXYOnTracks.cc HoughTest/HoughCheckXY/src/HoughCheckXYOnTracks.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Enrico Robutti
//         Created:  Fri, Jan 25, 2013
// $Id: HoughCheckXYOnTracks.cc
//
//


// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "TH3S.h"
#include "TFile.h"
#include "TTree.h"

//
// class declaration
//

using namespace std;
using namespace edm;

class HoughCheckXYOnTracks : public edm::EDAnalyzer {
public:
  explicit HoughCheckXYOnTracks(const edm::ParameterSet&);
  ~HoughCheckXYOnTracks();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // Parameters
  edm::InputTag trackTags_; //used to select what tracks to read from configuration file
  std::string builderName_;
  const TransientTrackingRecHitBuilder* builder_;
  const float maxDoca_;
  const int nBinsDoca_;
  const float maxKappa_ = 1./45.;  // larger than radius of inner pixel layer
  const int nBinsKappa_;
  const int nBinsPhi_;
  const unsigned int verbosity_;
  const unsigned int printProgressFrequency_;

  // Histogram and tree
  auto_ptr<TH3S> hHoughVotes_;
  auto_ptr<TTree> trackTree_;
  //  vector<float> vDoca_;
  //  vector<float> vKappa_;
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
HoughCheckXYOnTracks::HoughCheckXYOnTracks(const edm::ParameterSet& iConfig)
:
  trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks")),
  builderName_(iConfig.getParameter<std::string>("TTRHBuilder")),
  maxDoca_(iConfig.getUntrackedParameter<double>("maxDoca", 44.)),
  nBinsDoca_(iConfig.getUntrackedParameter<int>("nBinsDoca", 100)),
  maxKappa_(iConfig.getUntrackedParameter<double>("maxKappa", 1./45.)),
  nBinsKappa_(iConfig.getUntrackedParameter<int>("nBinsKappa", 100)),
  nBinsPhi_(iConfig.getUntrackedParameter<int>("nBinsPhi", 100)),
  verbosity_(iConfig.getUntrackedParameter<unsigned int>("verbosity", 0)),
  printProgressFrequency_(iConfig.getUntrackedParameter<unsigned int>("printProgressFrequency", 1000))
{

}


HoughCheckXYOnTracks::~HoughCheckXYOnTracks()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HoughCheckXYOnTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using reco::TrackCollection;
  
  // Get run and event number
  unsigned int runNum = iEvent.id().run();
  unsigned int evtNum = iEvent.id().event();
  
  // Print progress frequency
  if (verbosity_ > 0 && evtNum%printProgressFrequency_ == 0)
    cout << "Run n. " << runNum << ", Event n. " << evtNum << endl << endl;

  Handle<TrackCollection> tracks;
  iEvent.getByLabel(trackTags_,tracks);
  int ntrk = 0 ;
  for(TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack) {
    if (verbosity_ > 1)
      cout << "Track #" << ntrk << " with q=" << itTrack->charge() 
	   << ", pT=" << itTrack->pt() << " GeV, eta=" << itTrack->eta() 
	   << ", Nhits=" << itTrack->recHitsSize() 
	   << ", doca=" << 10.*(itTrack->dxy())
	   << ", kappa=" << 1.139e-3/(itTrack->pt())
	   << ", algo=" << itTrack->algoName(itTrack->algo()).c_str() << endl;
    //    vDoca_.push_back(10.*(itTrack->dxy()));
    //    vKappa_.push_back(1.139e-3/(itTrack->qoverp()));
    int nhit = 0;
    for (trackingRecHit_iterator i=itTrack->recHitsBegin(); i!=itTrack->recHitsEnd(); i++){
      if (verbosity_ > 2)
	cout << "hit #" << nhit;
      TransientTrackingRecHit::RecHitPointer hit = builder_->build(&**i );
//	DetId hitId = hit->geographicalId();
// 	if(hitId.det() == DetId::Tracker) {
// 	  if (hitId.subdetId() == StripSubdetector::TIB )  
// 	    cout << " - TIB " << TIBDetId(hitId).layer();
// 	  else if (hitId.subdetId() == StripSubdetector::TOB ) 
// 	    cout << " - TOB " << TOBDetId(hitId).layer();
// 	  else if (hitId.subdetId() == StripSubdetector::TEC ) 
// 	    cout << " - TEC " << TECDetId(hitId).wheel();
// 	  else if (hitId.subdetId() == StripSubdetector::TID ) 
// 	    cout << " - TID " << TIDDetId(hitId).wheel();
// 	  else if (hitId.subdetId() == StripSubdetector::TID ) 
// 	    cout << " - TID " << TIDDetId(hitId).wheel();
// 	  else if (hitId.subdetId() == (int) PixelSubdetector::PixelBarrel ) 
// 	    cout << " - PixBar " << PXBDetId(hitId).layer();
// 	  else if (hitId.subdetId() == (int) PixelSubdetector::PixelEndcap )
// 	    cout << " - PixFwd " << PXFDetId(hitId).disk();
// 	  else 
// 	    cout << " UNKNOWN TRACKER HIT TYPE ";
// 	}
      if (hit->isValid()) {
	GlobalPoint hitPosition = hit->globalPosition();
	if (verbosity_ > 2)
	  cout << " - globalPos = " << hitPosition << endl;
	double x = 10.*hitPosition.x();
	double y = 10.*hitPosition.y();
	double r = sqrt(x*x + y*y);
	double docaBinWidth = maxDoca_/nBinsDoca_;
	// Loop on allowed histogram bins (first doca, then kappa, then phi) and fill histogram
	for (double doca = docaBinWidth/2; doca < maxDoca_ && doca < r; doca += docaBinWidth) {
	  // Start from first allowed kappa bin (must be k >= -2/(r + doca))
	  double kappaBinWidth = 2*maxKappa_/nBinsKappa_;
	  double firstKappa = -maxKappa_ + kappaBinWidth/2;
	  if (-2/(r + doca) > firstKappa) {
	    int binPosition = (-2/(r + doca))/kappaBinWidth;
	    firstKappa = (binPosition - 0.5)*kappaBinWidth;
	  }
	  for (double kappa = firstKappa; kappa < 2/(r - doca); kappa += kappaBinWidth) {
	    double akappa = (2*doca + kappa*doca*doca + kappa*r*r)/(2*r);
	    double hkappa = sqrt((doca*kappa + 1)*(doca*kappa + 1) - akappa*akappa);
	    double kappax1 = akappa*x/r + hkappa*y/r;
	    double kappax2 = akappa*x/r - hkappa*y/r;
	    double kappay1 = akappa*y/r - hkappa*x/r;
	    double kappay2 = akappa*y/r + hkappa*x/r;
	    hHoughVotes_->Fill(doca, kappa, atan2(kappay1, kappax1));
	    hHoughVotes_->Fill(doca, kappa, atan2(kappay2, kappax2));
	  }
	}
      } else 
	if (verbosity_ > 2)
	  cout << " - invalid hit" << endl;
      nhit++;
    }
    if (verbosity_ > 2)
      cout << endl;
    ntrk++;
  }
  //  trackTree_->Fill();
  //  vDoca_.clear();
  //  vKappa_.clear();
}


// ------------ method called once each job just before starting event loop  ------------
void 
HoughCheckXYOnTracks::beginJob()
{
  // Book 3D histogram
  if (maxDoca_ <= 0. || maxDoca_ > 550) {
    cout << "Invalid maxDoca parameter (" << maxDoca_ << " m). Valid range is 0. < maxDoca <= 550 mm. No histogram booked." << endl;
    return;
  }
  if (nBinsDoca_ <= 0 || nBinsDoca_ > 200) {
    cout << "Invalid nBinsDoca parameter (" << nBinsDoca_ << "). Valid range is 0 < nBinsDoca <= 200. No histogram booked." << endl;
    return;
  }
  if (nBinsKappa_ <= 0 || nBinsKappa_ > 200) {
    cout << "Invalid nBinsKappa parameter (" << nBinsKappa_ << "). Valid range is 0 < nBinsKappa <= 200. No histogram booked." << endl;
    return;
  }
  if (nBinsPhi_ <= 0 || nBinsPhi_ > 200) {
    cout << "Invalid nBinsPhi parameter (" << nBinsPhi_ << "). Valid range is 0 < nBinsPhi <= 200. No histogram booked." << endl;
    return;
  }
  hHoughVotes_.reset(new TH3S("hHoughVotes", "circle Hough transform votes", nBinsDoca_, 0., maxDoca_, nBinsKappa_, -maxKappa_, maxKappa_, nBinsPhi_, -M_PI, M_PI));

  //  // Create tree, with branches
  //  trackTree_.reset(new TTree("trackTree", "Fitted track parameters"));
  //  trackTree_->Branch("doca", &vDoca_);
  //  trackTree_->Branch("kappa", &vKappa_);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HoughCheckXYOnTracks::endJob() 
{
  TFile outRootFile("houghCheckXY_tracks.root", "RECREATE");
  hHoughVotes_->Write();
  //  trackTree_->Write();
  outRootFile.Close();
}

// ------------ method called when starting to processes a run  ------------
void 
HoughCheckXYOnTracks::beginRun(edm::Run const& run, edm::EventSetup const& setup)
{

  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
  setup.get<TransientRecHitRecord>().get(builderName_,theBuilder);
  builder_=theBuilder.product();
}

// ------------ method called when ending the processing of a run  ------------
void 
HoughCheckXYOnTracks::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HoughCheckXYOnTracks::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HoughCheckXYOnTracks::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HoughCheckXYOnTracks::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(HoughCheckXYOnTracks);
