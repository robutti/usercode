// -*- C++ -*-
//
// Package:    HoughTransChecks
// Class:      HoughCheckAllHits
// 
/**\class HoughCheckAllHits HoughCheckAllHits.cc HoughTest/HoughTransChecks/src/HoughCheckAllHits.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Enrico Robutti
//         Created:  Thu, May 16, 2013
// $Id: HoughCheckAllHits.cc
//
//


// system include files
#include <memory>
#include <vector>
#include <algorithm>

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
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryParametrization/interface/PerigeeTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/PerigeeConversions.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"
#include "RecoTracker/TkSeedingLayers/interface/SeedingLayerSetsBuilder.h"

#include "ERobutti/HoughTransChecks/interface/TH5.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"


//
// class declaration
//

using namespace std;
using namespace edm;

class HoughCheckAllHits : public edm::EDAnalyzer {
public:
  explicit HoughCheckAllHits(const edm::ParameterSet&);
  ~HoughCheckAllHits();
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
  const vector<unsigned int> algoSel_;
  string  layerListName_;
  ctfseeding::SeedingLayerSets layerSets_;
  double minDoca_, maxDoca_;
  const int nBinsDoca_;
  double minDocaScan_, maxDocaScan_;
  double minSqrtK_, maxSqrtK_;
  const int nBinsSqrtK_;
  double minSqrtKScan_, maxSqrtKScan_;
  double minPhi_, maxPhi_;
  const int nBinsPhi_;
  double minZ0_, maxZ0_;
  const int nBinsZ0_;
  double minEta_, maxEta_;
  const int nBinsEta_;
  const unsigned int verbosity_;

  // Histogram and tree
  auto_ptr<TH5C> hHoughVotes_;
  auto_ptr<TTree> trackTree_;
  vector<double> vDoca_;
  vector<double> vKappa_;
  vector<double> vPhi_;
  vector<double> vZ0_;
  vector<double> vTheta_;
  vector<double> vAlgo_;

  // Layer mapping
  static const int lyrMapOffset_[7];
};

//
// constants, enums and typedefs
//
const int HoughCheckAllHits::lyrMapOffset_[7] = {0, 0, 3, 5, 9, 12, 18};

//
// static data member definitions
//

//
// constructors and destructor
//
HoughCheckAllHits::HoughCheckAllHits(const edm::ParameterSet& iConfig)
:
  trackTags_(iConfig.getParameter<edm::InputTag>("tracks")),
  algoSel_(iConfig.getParameter<vector<unsigned int> >("algoSel")),
  layerListName_ (iConfig.getParameter<std::string> ("seedingLayers")),
  nBinsDoca_(iConfig.getUntrackedParameter<int>("nBinsDoca", 100)),
  nBinsSqrtK_(iConfig.getUntrackedParameter<int>("nBinsSqrtK", 100)),
  nBinsPhi_(iConfig.getUntrackedParameter<int>("nBinsPhi", 100)),
  nBinsZ0_(iConfig.getUntrackedParameter<int>("nBinsZ0", 100)),
  nBinsEta_(iConfig.getUntrackedParameter<int>("nBinsEta", 100)),
  verbosity_(iConfig.getUntrackedParameter<unsigned int>("verbosity", 0))
{
  vector<double> rangeDoca = iConfig.getUntrackedParameter<vector<double> >("rangeDoca", {-44., 44.});  // radius of inner pixel layer
  minDoca_ = rangeDoca[0];
  maxDoca_ = rangeDoca[1];
  vector<double> rangeSqrtK = iConfig.getUntrackedParameter<vector<double> >("rangeSqrtK", {-0.075, 0.075});  // corresponding to ~200 MeV (~18 cm radius)
  minSqrtK_ = rangeSqrtK[0];
  maxSqrtK_ = rangeSqrtK[1];
  vector<double> rangePhi = iConfig.getUntrackedParameter<vector<double> >("rangePhi", {-M_PI, M_PI});
  minPhi_ = rangePhi[0];
  maxPhi_ = rangePhi[1];
  vector<double> rangeZ0 = iConfig.getUntrackedParameter<vector<double> >("rangeZ0", {-265., 265.});  // pixel barrel half-length
  minZ0_ = rangeZ0[0];
  maxZ0_ = rangeZ0[1];
  vector<double> rangeEta = iConfig.getUntrackedParameter<vector<double> >("rangeEta", {-2.5, 2.5});  // tracker eta acceptance
  minEta_ = rangeEta[0];
  maxEta_ = rangeEta[1];
}


HoughCheckAllHits::~HoughCheckAllHits()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HoughCheckAllHits::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using reco::TrackCollection;
   
  Handle<TrackCollection> tracks;
  iEvent.getByLabel(trackTags_,tracks);
  int ntrk = 0;
  vDoca_.clear();
  vKappa_.clear();
  vPhi_.clear();
  vZ0_.clear();
  vTheta_.clear();
  vAlgo_.clear();

  double docaBinWidth = (maxDoca_ - minDoca_)/nBinsDoca_;
  double sqrtKBinWidth = (maxSqrtK_ - minSqrtK_)/nBinsSqrtK_;
  double z0BinWidth = (maxZ0_ - minZ0_)/nBinsZ0_;
  
  for(TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack) {

    // Get fitted track parameters in perigee convention
    edm::ESHandle<MagneticField> theMF;
    iSetup.get<IdealMagneticFieldRecord>().get(theMF);
    GlobalPoint refPoint(itTrack->vx(), itTrack->vy(), itTrack->vz());
    GlobalVector pTrack(itTrack->px(), itTrack->py(), itTrack->pz());
    FreeTrajectoryState fts(refPoint, pTrack, TrackCharge(itTrack->charge()), theMF.product());
    TSCPBuilderNoMaterial tscpBuilder;
    TrajectoryStateClosestToPoint tsAtClosestApproach = tscpBuilder(fts, GlobalPoint(0.,0.,0.));
    PerigeeTrajectoryParameters perigeePars = tsAtClosestApproach.perigeeParameters();

    if (verbosity_ > 1)
      cout << "Track #" << ntrk << " with q=" << itTrack->charge() 
	   << ", Nhits=" << itTrack->recHitsSize()
	   << ", (vx,vy,vz)=(" << itTrack->vx() << "," << itTrack->vy() << "," << itTrack->vz() << ")"
           << ", doca=" << perigeePars.transverseImpactParameter()
           << ", kappa=" << perigeePars.transverseCurvature()
	   << ", phi=" << perigeePars.phi()
           << ", z0=" << perigeePars.longitudinalImpactParameter()
	   << ", theta=" << perigeePars.theta()
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
	continue;
    }
    vDoca_.push_back(10.*perigeePars.transverseImpactParameter());
    vKappa_.push_back(perigeePars.transverseCurvature()/10.);
    vPhi_.push_back(perigeePars.phi());
    vZ0_.push_back(10.*perigeePars.longitudinalImpactParameter());
    vTheta_.push_back(perigeePars.theta());
    vAlgo_.push_back(itTrack->algo());
  }

  if (!(hHoughVotes_.get()))
    return;
  // Loop on layers and hits within
  map<int, int> layersInBin;
  int nhit = 0;
  for (ctfseeding::SeedingLayerSets::const_iterator itLS = layerSets_.begin(); itLS != layerSets_.end(); itLS++ ) {
//     ctfseeding::SeedingLayer::Hits lyrSetHits;
    for (ctfseeding::SeedingLayers::const_iterator itLyr = itLS->begin(); itLyr!= itLS->end(); itLyr++ ) {
      ctfseeding::SeedingLayer::Hits hits = itLyr->hits(iEvent, iSetup);
//    lyrSetHits.insert(lyrSetHits.end(), hits.begin(), hits.end());
      for (ctfseeding::SeedingLayer::Hits::const_iterator itHit = hits.begin(); itHit != hits.end(); itHit++) {
	if (verbosity_ > 2)
	  cout << "hit #" << nhit++;
	DetId hitId = (*itHit)->geographicalId();
	if (hitId.det() != DetId::Tracker)
	  continue;
	int hitSubDet = hitId.subdetId();
	int hitLayer = 0;
	if (hitSubDet == PixelSubdetector::PixelBarrel) 
          hitLayer = PXBDetId(hitId).layer();
	else if (hitSubDet == PixelSubdetector::PixelEndcap)
	  hitLayer = PXFDetId(hitId).disk();
	else if (hitSubDet == StripSubdetector::TIB)  
	  hitLayer = TIBDetId(hitId).layer();
	else if (hitSubDet == StripSubdetector::TID) 
	  hitLayer = TIDDetId(hitId).wheel();
	else if (hitSubDet == StripSubdetector::TOB) 
	  hitLayer = TOBDetId(hitId).layer();
	else if (hitSubDet == StripSubdetector::TEC)
         hitLayer = TECDetId(hitId).wheel();
	else 
	  continue;
 	if ((*itHit)->isValid()) {
	  GlobalPoint hitPosition = (*itHit)->globalPosition();
	  if (verbosity_ > 2)
	    cout << " - globalPos = " << hitPosition << endl;
	  double x = 10.*hitPosition.x();
	  double y = 10.*hitPosition.y();
	  double z = 10.*hitPosition.z();
	  double r = sqrt(x*x + y*y);
	  //	  double phiHit = atan2(y, x);
	  // Loop on allowed histogram bins (first doca, then kappa, then phi) and fill histogram
// 	  for (double docaScan = minDocaScan_; docaScan < maxDocaScan_ + 0.5*docaBinWidth; docaScan += docaBinWidth) {
// 	    // Start from first allowed sqrtK bin (must be kappaScan >= -2/(r + docaScan))
// 	    double firstSqrtK = -maxSqrtKScan_;
// 	    if (-2/(r + docaScan) > -firstSqrtK*firstSqrtK) {
// 	      int binPosition = -sqrt(2/(r + docaScan))/sqrtKBinWidth;
// 	      firstSqrtK = (binPosition - 0.5)*sqrtKBinWidth;
// 	    }
// 	    for (double sqrtKScan = firstSqrtK; sqrtKScan < min(maxSqrtKScan_ + 0.5*sqrtKBinWidth, sqrt(2/(r - docaScan))); sqrtKScan += sqrtKBinWidth) {
// 	      if (sqrtKScan > -minSqrtKScan_ + 0.5*sqrtKBinWidth && sqrtKScan < minSqrtKScan_ - 0.5*sqrtKBinWidth)  // skip irrelevant values for histogram
// 		continue;
// 	      double kappaScan = sqrtKScan*fabs(sqrtKScan);
// 	      double akappa = (2*docaScan + kappaScan*docaScan*docaScan + kappaScan*r*r)/(2*r);
// 	      double hkappa = sqrt((docaScan*kappaScan + 1)*(docaScan*kappaScan + 1) - akappa*akappa);
// 	      for (int sign = -1; sign <= 1; sign += 2) {
// 		double kappax = akappa*x/r - sign*hkappa*y/r;
// 		double kappay = akappa*y/r + sign*hkappa*x/r;
// 		// Convert to perigee parameters
// 		double phiD = atan2(kappay, kappax);
// 		int rot = -2*(int((phiHit - phiD + 2*M_PI)/M_PI)%2) + 1;  // hit position wrt. poca: +1 = anticlockwise; -1 = clockwise
// 		double doca = rot*docaScan;
// 		//	      double kappa = -rot*kappaScan;
// 		double kappa = -rot*kappaScan;
// 		if (fabs(kappa) < 1.e-6)
// 		  kappa = (kappa >= 0) ? 1.e-6 : -1.e-6;  // avoid kappa = 0 (maximum curvature radius = 1 km)
// 		double phi = phiD + rot*(M_PI/2.);
// 		phi += -2.*M_PI*(int((phi + 3.*M_PI)/(2.*M_PI)) - 1);  // map to range (-PI, PI)
	  for (int iDoca = 0; iDoca < nBinsDoca_; iDoca++) {
	    double doca = minDoca_ + (iDoca + 0.5)*docaBinWidth;
	    for (int iSqrtK = 0; iSqrtK < nBinsSqrtK_; iSqrtK++) {
	      double sqrtK = minSqrtK_ + (iSqrtK + 0.5)*sqrtKBinWidth;
	      double kappa = sqrtK*fabs(sqrtK);
	      // Check for allowed range -2/(r - doca) <= kappa <= 2/(r + doca)
	      if (kappa < -2/(r - doca) || kappa > 2/(r + doca))
		continue;
	      double akappa = (kappa*doca*doca + kappa*r*r - 2*doca)/(2*r);
	      double hkappa = sqrt((1 - doca*kappa)*(1 - doca*kappa) - akappa*akappa);
	      double kappaxD = -akappa*x/r + hkappa*y/r;
	      double kappayD = -akappa*y/r - hkappa*x/r;
	      double phi = atan2(kappayD, kappaxD) + M_PI/2.;
	      phi += -2.*M_PI*(int((phi + 3.*M_PI)/(2.*M_PI)) - 1);  // map to range (-PI, PI)
	      if (fabs(kappa) < 1.e-6)
		kappa = (kappa >= 0) ? 1.e-6 : -1.e-6;  // avoid kappa = 0 (maximum curvature radius = 1 km)
	      double xc = (doca - 1./kappa)*sin(phi);
	      double yc = -(doca - 1./kappa)*cos(phi);
	      for (int iZ0 = 0; iZ0 < nBinsZ0_; iZ0++) {
		double z0 = minZ0_ + (iZ0 + 0.5)*z0BinWidth;
		double st = 1./kappa*(atan2(kappa*(y - yc), kappa*(x - xc)) - phi + M_PI/2.);
		if (st < 0)  // rotation angle has crossed the +/-PI border
		  st += 2.*M_PI/fabs(kappa);
		else if (st > 2.*M_PI/fabs(kappa))
		  st -= 2.*M_PI/fabs(kappa);
		double lambda = atan((z - z0)/st);
		double eta = -log(tan((M_PI/2. - lambda)/2.));
		//		hHoughVotes_->Fill(doca, kappa, phi, z0, eta);
		int key = hHoughVotes_->FindBin(doca, sqrtK, phi, z0, eta);
		if ((layersInBin[key] & (1 << (lyrMapOffset_[hitSubDet] + hitLayer - 1))) == 0) {
		  layersInBin[key] |= (1 << (lyrMapOffset_[hitSubDet] + hitLayer - 1));
		  hHoughVotes_->Fill(doca, sqrtK, phi, z0, eta);
		}
	      }
	    }
	  }
	} else 
	  if (verbosity_ > 2)
	    cout << " - invalid hit" << endl;
      }
    }
  }

  trackTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
HoughCheckAllHits::beginJob()
{
  // Create tree, with branches
  trackTree_.reset(new TTree("trackTree", "Fitted track parameters"));
  trackTree_->Branch("doca", &vDoca_);
  trackTree_->Branch("kappa", &vKappa_);
  trackTree_->Branch("phi", &vPhi_);
  trackTree_->Branch("z0", &vZ0_);
  trackTree_->Branch("theta", &vTheta_);
  trackTree_->Branch("algo", &vAlgo_);
  trackTree_->SetDirectory(0);

  // Set edges for histogram axes and scan ranges
  if (nBinsDoca_ <= 0 || nBinsDoca_ > 200) {
    cout << "Invalid nBinsDoca parameter (" << nBinsDoca_ << "). Valid range is 0 < nBinsDoca <= 200. No histogram booked." << endl;
    return;
  }
  if (minDoca_ >= maxDoca_) {
    cout << "Invalid doca range: min(doca) >= max(doca). No histogram booked." << endl;
    return;
  }
  if (nBinsSqrtK_ <= 0 || nBinsSqrtK_ > 200) {
    cout << "Invalid nBinsSqrtK parameter (" << nBinsSqrtK_ << "). Valid range is 0 < nBinsSqrtK <= 200. No histogram booked." << endl;
    return;
  }
  if (minSqrtK_ >= maxSqrtK_) {
    cout << "Invalid sqrtK range: min(sqrtK) >= max(sqrtK). No histogram booked." << endl;
    return;
  }
  if (nBinsPhi_ <= 0 || nBinsPhi_ > 200) {
    cout << "Invalid nBinsPhi parameter (" << nBinsPhi_ << "). Valid range is 0 < nBinsPhi <= 200. No histogram booked." << endl;
    return;
  }
  if (minPhi_ >= maxPhi_) {
    cout << "Invalid phi range: min(phi) >= max(phi). No histogram booked." << endl;
    return;
  } else if (minPhi_ < -2*M_PI || maxPhi_ > 2*M_PI) {
    cout << "Invalid phi range. Valid range is -pi <= phi <= pi. No histogram booked." << endl;
    return;
  }
  if (nBinsZ0_ <= 0 || nBinsZ0_ > 200) {
    cout << "Invalid nBinsZ0 parameter (" << nBinsZ0_ << "). Valid range is 0 < nBinsZ0 <= 200. No histogram booked." << endl;
    return;
  }
  if (minZ0_ >= maxZ0_) {
    cout << "Invalid z0 range: min(z0) >= max(z0). No histogram booked." << endl;
    return;
  }
  if (nBinsEta_ <= 0 || nBinsEta_ > 200) {
    cout << "Invalid nBinsEta parameter (" << nBinsEta_ << "). Valid range is 0 < nBinsEta <= 200. No histogram booked." << endl;
    return;
  }
  if (minEta_ >= maxEta_) {
    cout << "Invalid eta range: min(eta) >= max(eta). No histogram booked." << endl;
    return;
  }
  // Book 5D histogram
  hHoughVotes_.reset(new TH5C("hHoughVotes", "helix Hough transform votes", nBinsDoca_, minDoca_, maxDoca_, nBinsSqrtK_, minSqrtK_, maxSqrtK_,
			      nBinsPhi_, minPhi_, maxPhi_, nBinsZ0_, minZ0_, maxZ0_, nBinsEta_, minEta_, maxEta_));
  hHoughVotes_->SetDirectory(0);

  // Set content to -128 for all histogram bins, to use full 8-bit range
  for (int iDoca = 1; iDoca <= nBinsDoca_; iDoca++)
    for (int iSqrtK = 1; iSqrtK <= nBinsSqrtK_; iSqrtK++)
      for (int iPhi = 1; iPhi <= nBinsPhi_; iPhi++)
	for (int iZ0 = 1; iZ0 <= nBinsZ0_; iZ0++)
	  for (int iEta = 1; iEta <= nBinsEta_; iEta++)
	    hHoughVotes_->SetBinContent(iDoca, iSqrtK, iPhi, iZ0, iEta, -128);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HoughCheckAllHits::endJob() 
{
  TFile outRootFile("houghCheck_hits.root", "RECREATE");
  hHoughVotes_->Write();
  trackTree_->Write();
}

// ------------ method called when starting to processes a run  ------------
void 
HoughCheckAllHits::beginRun(edm::Run const& run, edm::EventSetup const& setup)
{
  cout << "SeedingLayerSetsBuilder name: " << layerListName_ << endl;
  ESHandle<SeedingLayerSetsBuilder> layerBuilder;
  setup.get<TrackerDigiGeometryRecord>().get(layerListName_.c_str(), layerBuilder);
  layerSets_ = layerBuilder->layers(setup); 
  int i = 0;
  for (ctfseeding::SeedingLayerSets::const_iterator itLS = layerSets_.begin(); itLS != layerSets_.end(); itLS++) {
    cout << "SeedingLayerSet number " << ++i << endl;
    for (ctfseeding::SeedingLayers::const_iterator itLyr = itLS->begin(); itLyr != itLS->end(); itLyr++) {
      std::cout << "  " << itLyr->name() << std::endl; 
    }
  }
}

// ------------ method called when ending the processing of a run  ------------
void 
HoughCheckAllHits::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HoughCheckAllHits::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HoughCheckAllHits::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HoughCheckAllHits::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<vector<unsigned int> >("algoSel", {});
  desc.add<std::string>("seedingLayers", "HoughTransformSeedLayersAllHitsOneSet");
  desc.addUntracked<int>("nBinsDoca", 100);
  desc.addUntracked<int>("nBinsSqrtK", 100);
  desc.addUntracked<int>("nBinsPhi", 100);
  desc.addUntracked<int>("nBinsZ0", 100);
  desc.addUntracked<int>("nBinsEta", 100);
  desc.addUntracked<unsigned int>("verbosity", 0);
  desc.addUntracked<vector<double> >("rangeDoca", {-44., 44.});  // radius of inner pixel layer
  desc.addUntracked<vector<double> >("rangeSqrtK", {-0.075, 0.075});  // corresponding to ~200 MeV (~18 cm radius)
  desc.addUntracked<vector<double> >("rangePhi", {-M_PI, M_PI});
  desc.addUntracked<vector<double> >("rangeZ0", {-265., 265.});  // pixel barrel half-length
  desc.addUntracked<vector<double> >("rangeEta", {-2.5, 2.5});  // tracker eta acceptance
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HoughCheckAllHits);
