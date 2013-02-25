// -*- C++ -*-
//
// Package:    HoughTransChecks
// Class:      HoughCheckOnTracks
// 
/**\class HoughCheckOnTracks HoughCheckOnTracks.cc HoughTest/HoughTransChecks/src/HoughCheckOnTracks.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Enrico Robutti
//         Created:  Fri, Feb 15, 2013
// $Id: HoughCheckOnTracks.cc
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

#include "ERobutti/HoughTransChecks/interface/TH5.h"
#include "TFile.h"
#include "TTree.h"

//
// class declaration
//

using namespace std;
using namespace edm;

class HoughCheckOnTracks : public edm::EDAnalyzer {
public:
  explicit HoughCheckOnTracks(const edm::ParameterSet&);
  ~HoughCheckOnTracks();
  
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
  const vector<unsigned int> algoSel_;
  const vector<unsigned int> layerSel_;
  const float maxDoca_;
  const int nBinsDoca_;
  const float maxKappa_;
  const int nBinsKappa_;
  const int nBinsPhi_;
  const float maxAbsZ0_;
  const int nBinsZ0_;
  const float maxAbsLambda_;
  const int nBinsLambda_;
  const unsigned int verbosity_;
  const unsigned int printProgressFrequency_;

  // Histogram and tree
  auto_ptr<TH5C> hHoughVotes_;
  auto_ptr<TTree> trackTree_;
  vector<double> vDoca_;
  vector<double> vKappa_;
  vector<double> vPhi_;
  vector<double> vZ0_;
  vector<double> vLambda_;
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
HoughCheckOnTracks::HoughCheckOnTracks(const edm::ParameterSet& iConfig)
:
  trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks")),
  builderName_(iConfig.getParameter<std::string>("TTRHBuilder")),
  algoSel_(iConfig.getParameter<vector<unsigned int> >("algoSel")),
  layerSel_(iConfig.getParameter<vector<unsigned int> >("layerSel")),
  maxDoca_(iConfig.getUntrackedParameter<double>("maxDoca", 44.)),  // radius of inner pixel layer
  nBinsDoca_(iConfig.getUntrackedParameter<int>("nBinsDoca", 100)),
  maxKappa_(iConfig.getUntrackedParameter<double>("maxKappa", 1./45.)),  // larger than radius of inner pixel layer
  nBinsKappa_(iConfig.getUntrackedParameter<int>("nBinsKappa", 100)),
  nBinsPhi_(iConfig.getUntrackedParameter<int>("nBinsPhi", 100)),
  maxAbsZ0_(iConfig.getUntrackedParameter<double>("maxAbsZ0", 265.)),  // pixel barrel half-length
  nBinsZ0_(iConfig.getUntrackedParameter<int>("nBinsZ0", 100)),
  maxAbsLambda_(iConfig.getUntrackedParameter<double>("maxAbsLambda", 1.40)),  // tracker eta acceptance
  nBinsLambda_(iConfig.getUntrackedParameter<int>("nBinsLambda", 100)),
  verbosity_(iConfig.getUntrackedParameter<unsigned int>("verbosity", 0)),
  printProgressFrequency_(iConfig.getUntrackedParameter<unsigned int>("printProgressFrequency", 1000))
{
}


HoughCheckOnTracks::~HoughCheckOnTracks()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HoughCheckOnTracks::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
  int ntrk = 0;
  vDoca_.clear();
  vKappa_.clear();
  vPhi_.clear();
  for(TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack) {
    if (verbosity_ > 1)
      cout << "Track #" << ntrk << " with q=" << itTrack->charge() 
	   << ", pT=" << itTrack->pt() << " GeV, eta=" << itTrack->eta() 
	   << ", Nhits=" << itTrack->recHitsSize()
	   << ", (vx,vy,vz)=(" << itTrack->vx() << "," << itTrack->vy() << "," << itTrack->vz() << ")"
	   << ", doca=" << fabs(10.*(itTrack->dxy()))
	   << ", kappa=" << 1.139e-3/(itTrack->pt())
	   << ", phi=" << atan2(itTrack->vy(), itTrack->vx())
	   << ", algo=" << itTrack->algoName(itTrack->algo()).c_str() << endl;
    if (algoSel_.size() != 0) {  // empty vector means 'keep all tracks'
      bool validAlgo = false;
      for (vector<unsigned int>::const_iterator itAlgo = algoSel_.begin(); itAlgo != algoSel_.end(); itAlgo++) {
	if (itTrack->algo() == *itAlgo) {
	  validAlgo = true;
	  break;
	}
      }
      if (!validAlgo)
	break;
    }
    vDoca_.push_back(10.*(itTrack->d0()));
    vKappa_.push_back(1.139e-3*itTrack->charge()/(itTrack->pt()));
    vPhi_.push_back(itTrack->phi());
    int nhit = 0;
    for (trackingRecHit_iterator i=itTrack->recHitsBegin(); i!=itTrack->recHitsEnd(); i++){
      if (verbosity_ > 2)
	cout << "hit #" << nhit;
      TransientTrackingRecHit::RecHitPointer hit = builder_->build(&**i );
      // Keep only hits from selected layers
      DetId hitId = hit->geographicalId();
      if(hitId.det() != DetId::Tracker)
	break;
      if (layerSel_.size() != 0) {  // empty vector means 'keep all hits'
	bool validLyr = false;
	for (vector<unsigned int>::const_iterator itLyr = layerSel_.begin(); itLyr != layerSel_.end(); itLyr++) {
	  int subdetSel = *itLyr/10;
	  unsigned int lyrSel = *itLyr%10;
 	  if (hitId.subdetId() == subdetSel)
	    switch (subdetSel) {
	    case (PixelSubdetector::PixelBarrel):
	      if (PXBDetId(hitId).layer() == lyrSel) {
		validLyr = true;
		break;
	      }
	    case (PixelSubdetector::PixelEndcap):
	      if (PXFDetId(hitId).disk() == lyrSel) {
		validLyr = true;
		break;
	      }
	    case (StripSubdetector::TIB):
	      if (TIBDetId(hitId).layer() == lyrSel) {
		validLyr = true;
		break;
	      }
	    case (StripSubdetector::TID):
	      if (TIDDetId(hitId).ring() == lyrSel) {
		validLyr = true;
		break;
	      }
	    case (StripSubdetector::TOB):
	      if (TOBDetId(hitId).layer() == lyrSel) {
		validLyr = true;
		break;
	      }
	    case (StripSubdetector::TEC):
	      if (TECDetId(hitId).ring() == lyrSel) {
		validLyr = true;
		break;
	      }
	    default:
	      break;
	    }
 	}
	if (!validLyr)
	  break;
      }
      if (hit->isValid()) {
	GlobalPoint hitPosition = hit->globalPosition();
	if (verbosity_ > 2)
	  cout << " - globalPos = " << hitPosition << endl;
	double x = 10.*hitPosition.x();
	double y = 10.*hitPosition.y();
	double r = sqrt(x*x + y*y);
	double docaBinWidth = maxDoca_/nBinsDoca_;
	double kappaBinWidth = 2*maxKappa_/nBinsKappa_;
	double z0BinWidth = 2*maxAbsZ0_/nBinsZ0_;
	// Loop on allowed histogram bins (first doca, then kappa, then phi) and fill histogram
	for (double docaH = docaBinWidth/2; docaH < maxDoca_ && docaH < r; docaH += docaBinWidth) {
	  // Start from first allowed kappa bin (must be k >= -2/(r + docaH))
	  double firstKappa = -maxKappa_ + kappaBinWidth/2;
	  if (-2/(r + docaH) > firstKappa) {
	    int binPosition = (-2/(r + docaH))/kappaBinWidth;
	    firstKappa = (binPosition - 0.5)*kappaBinWidth;
	  }
	  for (double kappaH = firstKappa; kappaH < 2/(r - docaH); kappaH += kappaBinWidth) {
	    double akappa = (2*docaH + kappaH*docaH*docaH + kappaH*r*r)/(2*r);
	    double hkappa = sqrt((docaH*kappaH + 1)*(docaH*kappaH + 1) - akappa*akappa);
	    for (int sign = -1; sign <= 1; sign += 2) {
	      double kappax = akappa*x/r - sign*hkappa*y/r;
	      double kappay = akappa*y/r + sign*hkappa*x/r;
	      // Convert to perigee parameters
	      double phiHit = atan2(y, x);
	      double phiD = atan2(kappay, kappax);
	      int rot = -2*(int((phiHit - phiD + 2*M_PI)/M_PI)%2) + 1;  // hit position wrt. poca: +1 = anticlockwise; -1 = clockwise
	      double doca = rot*docaH;
	      double kappa = rot*kappaH;
	      double phi = phiD + rot*(M_PI/2.);
	      phi += -2.*M_PI*(int((phi + 3.*M_PI)/(2.*M_PI)) - 1);  // map to range (-PI, PI)
	      double xc = (doca - 1./kappa)*sin(phi);
	      double yc = -(doca - 1./kappa)*cos(phi);
	      double psi = (kappa*doca > 0) ? phi - M_PI/2. : phi + M_PI/2.;
	      for (double z0 = -maxAbsZ0_ + z0BinWidth/2; z0 < maxAbsZ0_; z0 += z0BinWidth) {
		double lambda = 1./kappa*(atan2(kappa*(y - yc), kappa*(x - xc)) - psi);
		hHoughVotes_->Fill(doca, kappa, phi, z0, lambda);
	      }
	    }
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
  trackTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
HoughCheckOnTracks::beginJob()
{
  // Book 5D histogram
  if (maxDoca_ <= 0. || maxDoca_ > 550.) {
    cout << "Invalid maxDoca parameter (" << maxDoca_ << " mm). Valid range is 0. < maxDoca <= 550 mm. No histogram booked." << endl;
    return;
  }
  if (nBinsDoca_ <= 0 || nBinsDoca_ > 200) {
    cout << "Invalid nBinsDoca parameter (" << nBinsDoca_ << "). Valid range is 0 < nBinsDoca <= 200. No histogram booked." << endl;
    return;
  }
  if (maxDoca_ <= 0. || maxKappa_ > 1./45.) {
    cout << "Invalid maxKappa parameter (" << maxKappa_ << " mm^-1). Valid range is 0. < maxKappa <= 1/45 mm. No histogram booked." << endl;
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
  if (maxAbsZ0_ <= 0 || maxAbsZ0_ > 600) {
    cout << "Invalid maxAbsZ0 parameter (" << maxAbsZ0_ << " mm). Valid range is 0 < maxAbsZ0 <= 600. No histogram booked." << endl;
    return;
  }
  if (nBinsZ0_ <= 0 || nBinsZ0_ > 200) {
    cout << "Invalid nBinsZ0 parameter (" << nBinsZ0_ << "). Valid range is 0 < nBinsZ0 <= 200. No histogram booked." << endl;
    return;
  }
  if (maxAbsLambda_ <= 0 || maxAbsLambda_ > 1.4) {
    cout << "Invalid maxAbsLambda parameter (" << maxAbsLambda_ << "). Valid range is 0 < maxAbsLambda <= 1.4. No histogram booked." << endl;
    return;
  }
  if (nBinsLambda_ <= 0 || nBinsLambda_ > 200) {
    cout << "Invalid nBinsLambda parameter (" << nBinsLambda_ << "). Valid range is 0 < nBinsLambda <= 200. No histogram booked." << endl;
    return;
  }
  hHoughVotes_.reset(new TH5C("hHoughVotes", "helix Hough transform votes", nBinsDoca_, -maxDoca_, maxDoca_, nBinsKappa_, -maxKappa_, maxKappa_,
			      nBinsPhi_, -M_PI, M_PI, nBinsZ0_, -maxAbsZ0_, maxAbsZ0_, nBinsLambda_, -maxAbsLambda_, maxAbsLambda_));
  hHoughVotes_->SetDirectory(0);

  // Set content to -128 for all histogram bins, to use full 8-bit range
  for (int iDoca = 1; iDoca <= nBinsDoca_; iDoca++)
    for (int iKappa = 1; iKappa <= nBinsKappa_; iKappa++)
      for (int iPhi = 1; iPhi <= nBinsPhi_; iPhi++)
	for (int iZ0 = 1; iZ0 <= nBinsZ0_; iZ0++)
	  for (int iLambda = 1; iLambda <= nBinsLambda_; iLambda++)
	    hHoughVotes_->SetBinContent(iDoca, iKappa, iPhi, iZ0, iLambda, -128);

  //  // Create tree, with branches
  trackTree_.reset(new TTree("trackTree", "Fitted track parameters"));
  trackTree_->Branch("doca", &vDoca_);
  trackTree_->Branch("kappa", &vKappa_);
  trackTree_->Branch("phi", &vPhi_);
  trackTree_->Branch("z0", &vZ0_);
  trackTree_->Branch("lambda", &vLambda_);
  trackTree_->SetDirectory(0);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HoughCheckOnTracks::endJob() 
{
  TFile outRootFile("houghCheck_tracks.root", "RECREATE");
  hHoughVotes_->Write();
  trackTree_->Write();
  outRootFile.Close();
}

// ------------ method called when starting to processes a run  ------------
void 
HoughCheckOnTracks::beginRun(edm::Run const& run, edm::EventSetup const& setup)
{

  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
  setup.get<TransientRecHitRecord>().get(builderName_,theBuilder);
  builder_=theBuilder.product();
}

// ------------ method called when ending the processing of a run  ------------
void 
HoughCheckOnTracks::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HoughCheckOnTracks::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HoughCheckOnTracks::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HoughCheckOnTracks::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
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
DEFINE_FWK_MODULE(HoughCheckOnTracks);