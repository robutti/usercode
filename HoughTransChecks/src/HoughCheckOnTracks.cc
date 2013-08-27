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
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryParametrization/interface/PerigeeTrajectoryParameters.h"
#include "TrackingTools/TrajectoryState/interface/PerigeeConversions.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"

#include "ERobutti/HoughTransChecks/interface/TH5.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "Math/RootFinder.h"
#include "Math/WrappedTF1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TH2S.h"
#include "TEllipse.h"


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

  double fZDoca_(double* var, double* par);
  double fZKappa_(double* var, double* par);
  double fPhi_(double x, double y, int sign);
  double fZ0_(double x, double y, double z);
  double fTheta_(double x, double y, double z);

  // Parameters
  edm::InputTag trackTags_; //used to select what tracks to read from configuration file
  std::string builderName_;
  const TransientTrackingRecHitBuilder* builder_;
  const vector<unsigned int> algoSel_;
  const vector<unsigned int> layerSel_;
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
  const unsigned int printProgressFrequency_;

  // Histogram and tree
  auto_ptr<TH5C> hHoughVotes_;
  auto_ptr<TTree> trackTree_;
  vector<double> vDoca_;
  vector<double> vKappa_;
  vector<double> vPhi_;
  vector<double> vZ0_;
  vector<double> vTheta_;
  vector<double> vAlgo_;
  vector<vector<double> > vDocaExp_;
  vector<vector<double> > vKappaExp_;
  vector<vector<double> > vPhiExp_;
  vector<vector<double> > vZ0Exp_;
  vector<vector<double> > vThetaExp_;
  vector<vector<double> > vXHit_;
  vector<vector<double> > vYHit_;
  vector<vector<double> > vZHit_;
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
  nBinsDoca_(iConfig.getUntrackedParameter<int>("nBinsDoca", 100)),
  nBinsSqrtK_(iConfig.getUntrackedParameter<int>("nBinsSqrtK", 100)),
  nBinsPhi_(iConfig.getUntrackedParameter<int>("nBinsPhi", 100)),
  nBinsZ0_(iConfig.getUntrackedParameter<int>("nBinsZ0", 100)),
  nBinsEta_(iConfig.getUntrackedParameter<int>("nBinsEta", 100)),
  verbosity_(iConfig.getUntrackedParameter<unsigned int>("verbosity", 0)),
  printProgressFrequency_(iConfig.getUntrackedParameter<unsigned int>("printProgressFrequency", 1000))
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
  vZ0_.clear();
  vTheta_.clear();
  vAlgo_.clear();
  vDocaExp_.clear();
  vKappaExp_.clear();
  vPhiExp_.clear();
  vZ0Exp_.clear();
  vThetaExp_.clear();

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
    if (!(hHoughVotes_.get()))
      return;
    int nhit = 0;
    vector<double> vDocaHit;
    vector<double> vKappaHit;
    vector<double> vPhiHit;
    vector<double> vZ0Hit;
    vector<double> vThetaHit;
    vector<double> vX;
    vector<double> vY;
    vector<double> vZ;
    for (trackingRecHit_iterator i = itTrack->recHitsBegin(); i != itTrack->recHitsEnd(); i++){
      if (verbosity_ > 2)
 	cout << "hit #" << nhit++;
      TransientTrackingRecHit::RecHitPointer hit = builder_->build(&**i );
      // Keep only hits from selected layers
      DetId hitId = hit->geographicalId();
      if(hitId.det() != DetId::Tracker)
	continue;
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
	  continue;
      }
      if (hit->isValid()) {
	GlobalPoint hitPosition = hit->globalPosition();
	if (verbosity_ > 2)
	  cout << " - globalPos = " << hitPosition << endl;
	double x = 10.*hitPosition.x();
	vX.push_back(x);
	double y = 10.*hitPosition.y();
	vY.push_back(y);
	double z = 10.*hitPosition.z();
	vZ.push_back(z);
	double r = sqrt(x*x + y*y);
	double phiHit = atan2(y, x);
	// Get expected parameters based on fitted remaining ones
	ROOT::Math::RootFinder rootFinder;
	// Expected DOCA
	TF1 tfZDoca("tfZDoca", this, &HoughCheckOnTracks::fZDoca_, 0., 100., 4, "this", "fZDoca_");
	ROOT::Math::WrappedTF1 wfZDoca(tfZDoca);
	rootFinder.SetFunction(wfZDoca, 0., 100.);
	// Get both solutions and keep closer to fitted value
	double zPar[4];
	zPar[0] = x;
	zPar[1] = y;
	zPar[2] = 1;
	zPar[3] = 1000.;
	wfZDoca.SetParameters(zPar);
	rootFinder.Solve();
	double doca1 = (fabs(wfZDoca.Parameters()[3]) < 0.05) ? rootFinder.Root() : 1.e10;
	doca1 *= vDoca_.back()/fabs(vDoca_.back());  // put "right" sign
	zPar[2] = -1;
	zPar[3] = 1000.;
	wfZDoca.SetParameters(zPar);
	rootFinder.Solve();
	double doca2 = (fabs(wfZDoca.Parameters()[3]) < 0.05) ? rootFinder.Root() : 1.e10;
	doca2 *= vDoca_.back()/fabs(vDoca_.back());  // put "right" sign
	double bestDoca = (fabs(doca1 - vDoca_.back()) < fabs(doca2 - vDoca_.back())) ? doca1 : doca2;
	vDocaHit.push_back(bestDoca);
	// Expected kappa
	double kappaMMin = (vDoca_.back()*vKappa_.back() < 0) ? 1.e-6 : -2./(r + vDoca_.back());
	double kappaMMax = (vDoca_.back()*vKappa_.back() < 0) ? 0.01 : -1.e-6;
	TF1 tfZKappa("tfZKappa", this, &HoughCheckOnTracks::fZKappa_, kappaMMin, kappaMMax, 4, "this", "fZKappa_");
	ROOT::Math::WrappedTF1 wfZKappa(tfZKappa);
	rootFinder.SetFunction(wfZKappa, kappaMMin, kappaMMax);
	// Get both solutions and keep closer to fitted value
	zPar[2] = 1;
	zPar[3] = 1000.;
	wfZKappa.SetParameters(zPar);
	rootFinder.Solve();
	double kappa1 = (fabs(wfZKappa.Parameters()[3]) < 0.05) ? rootFinder.Root() : 1.e10;
	kappa1 = (kappa1*vKappa_.back() > 0) ? kappa1 : -kappa1;  // put "right" sign
	zPar[2] = -1;
	zPar[3] = 1000.;
	wfZKappa.SetParameters(zPar);
	rootFinder.Solve();
	double kappa2 = (fabs(wfZKappa.Parameters()[3]) < 0.05) ? rootFinder.Root() : 1.e10;
	kappa2 = (kappa2*vKappa_.back() > 0) ? kappa2 : -kappa2;  // put "right" sign
	double bestKappa = (fabs(kappa1 - vKappa_.back()) < fabs(kappa2 - vKappa_.back())) ? kappa1 : kappa2;
	vKappaHit.push_back(bestKappa);
	// Expected phi
	double phi1 = fPhi_(x, y, 1);
	double phi2 = fPhi_(x, y, -1);
	double bestPhi = (fabs(phi1 - vPhi_.back()) < fabs(phi2 - vPhi_.back())) ? phi1 : phi2;
	vPhiHit.push_back(bestPhi);
	// Expected z0
	vZ0Hit.push_back(fZ0_(x, y, z));
	// Expected theta
	vThetaHit.push_back(fTheta_(x, y, z));
	// Loop on histogram bins
	double docaBinWidth = (maxDoca_ - minDoca_)/nBinsDoca_;
	double sqrtKBinWidth = (maxSqrtK_ - minSqrtK_)/nBinsSqrtK_;
	double z0BinWidth = (maxZ0_ - minZ0_)/nBinsZ0_;
	// Loop on allowed histogram bins (first doca, then kappa, then phi) and fill histogram
	for (double docaScan = minDocaScan_; docaScan < maxDocaScan_ + 0.5*docaBinWidth; docaScan += docaBinWidth) {
	  // Start from first allowed sqrtK bin (must be kappaScan >= -2/(r + docaScan))
	  double firstSqrtK = -maxSqrtKScan_;
	  if (-2/(r + docaScan) > -firstSqrtK*firstSqrtK) {
	    int binPosition = -sqrt(2/(r + docaScan))/sqrtKBinWidth;
	    firstSqrtK = (binPosition - 0.5)*sqrtKBinWidth;
	  }
	  for (double sqrtKScan = firstSqrtK; sqrtKScan < min(maxSqrtKScan_ + 0.5*sqrtKBinWidth, sqrt(2/(r - docaScan))); sqrtKScan += sqrtKBinWidth) {
	    if (sqrtKScan > -minSqrtKScan_ + 0.5*sqrtKBinWidth && sqrtKScan < minSqrtKScan_ - 0.5*sqrtKBinWidth)  // skip irrelevant values for histogram
	      continue;
	    double kappaScan = sqrtKScan*fabs(sqrtKScan);
	    double akappa = (2*docaScan + kappaScan*docaScan*docaScan + kappaScan*r*r)/(2*r);
	    double hkappa = sqrt((docaScan*kappaScan + 1)*(docaScan*kappaScan + 1) - akappa*akappa);
	    for (int sign = -1; sign <= 1; sign += 2) {
	      double kappax = akappa*x/r - sign*hkappa*y/r;
	      double kappay = akappa*y/r + sign*hkappa*x/r;
	      // Convert to perigee parameters
	      double phiD = atan2(kappay, kappax);
	      int rot = -2*(int((phiHit - phiD + 2*M_PI)/M_PI)%2) + 1;  // hit position wrt. poca: +1 = anticlockwise; -1 = clockwise
	      double doca = rot*docaScan;
	      //	      double kappa = -rot*kappaScan;
	      double kappa = -rot*kappaScan;
	      if (fabs(kappa) < 1.e-6)
		kappa = (kappa >= 0) ? 1.e-6 : -1.e-6;  // avoid kappa = 0 (maximum curvature radius = 1 km)
	      double phi = phiD + rot*(M_PI/2.);
	      phi += -2.*M_PI*(int((phi + 3.*M_PI)/(2.*M_PI)) - 1);  // map to range (-PI, PI)
	      double xc = (doca - 1./kappa)*sin(phi);
	      double yc = -(doca - 1./kappa)*cos(phi);
	      for (double z0 = minZ0_ + z0BinWidth/2; z0 < maxZ0_; z0 += z0BinWidth) {
		double st = 1./kappa*(atan2(kappa*(y - yc), kappa*(x - xc)) - phi + M_PI/2.);
		if (st < 0)  // rotation angle has crossed the +/-PI border
		  st += 2.*M_PI/fabs(kappa);
		else if (st > 2.*M_PI/fabs(kappa))
		  st -= 2.*M_PI/fabs(kappa);
		double lambda = atan((z - z0)/st);
		double eta = -log(tan((M_PI/2. - lambda)/2.));
//		hHoughVotes_->Fill(doca, kappa, phi, z0, eta);
		hHoughVotes_->Fill(doca, -rot*sqrtKScan, phi, z0, eta);
	      }
	    }
	  }
	}
      } else 
	if (verbosity_ > 2)
	  cout << " - invalid hit" << endl;
    }
    vDocaExp_.push_back(vDocaHit);
    vKappaExp_.push_back(vKappaHit);
    vPhiExp_.push_back(vPhiHit);
    vZ0Exp_.push_back(vZ0Hit);
    vThetaExp_.push_back(vThetaHit);
    vXHit_.push_back(vX);
    vYHit_.push_back(vY);
    vZHit_.push_back(vZ);
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
  // Create tree, with branches
  trackTree_.reset(new TTree("trackTree", "Fitted track parameters"));
  trackTree_->Branch("doca", &vDoca_);
  trackTree_->Branch("kappa", &vKappa_);
  trackTree_->Branch("phi", &vPhi_);
  trackTree_->Branch("z0", &vZ0_);
  trackTree_->Branch("theta", &vTheta_);
  trackTree_->Branch("algo", &vAlgo_);
  trackTree_->Branch("docaExp", &vDocaExp_);
  trackTree_->Branch("kappaExp", &vKappaExp_);
  trackTree_->Branch("phiExp", &vPhiExp_);
  trackTree_->Branch("z0Exp", &vZ0Exp_);
  trackTree_->Branch("thetaExp", &vThetaExp_);
  trackTree_->Branch("xHit", &vXHit_);
  trackTree_->Branch("yHit", &vYHit_);
  trackTree_->Branch("zHit", &vZHit_);
  trackTree_->SetDirectory(0);

  // Set edges for histogram axes and scan ranges
  if (nBinsDoca_ <= 0 || nBinsDoca_ > 200) {
    cout << "Invalid nBinsDoca parameter (" << nBinsDoca_ << "). Valid range is 0 < nBinsDoca <= 200. No histogram booked." << endl;
    return;
  }
  float docaBinWidth = (maxDoca_ - minDoca_)/nBinsDoca_; 
  if (docaBinWidth <= 0) {
    cout << "Invalid doca range: min(doca) >= max(doca). No histogram booked." << endl;
    return;
  } else if (minDoca_ < 0 && maxDoca_ > 0) {  // shift the range so that 0 is a bin edge
    float shift = docaBinWidth*int((maxDoca_ + 0.5*docaBinWidth)/docaBinWidth) - maxDoca_;
    minDoca_ += shift;
    maxDoca_ += shift;
    minDocaScan_ = 0.5*docaBinWidth;
    maxDocaScan_ = max(fabs(minDoca_), fabs(maxDoca_)) - 0.5*docaBinWidth;
  } else {
    minDocaScan_ = min(fabs(minDoca_), fabs(maxDoca_)) + 0.5*docaBinWidth;
    maxDocaScan_ = max(fabs(minDoca_), fabs(maxDoca_)) - 0.5*docaBinWidth;
  }
  if (nBinsSqrtK_ <= 0 || nBinsSqrtK_ > 200) {
    cout << "Invalid nBinsSqrtK parameter (" << nBinsSqrtK_ << "). Valid range is 0 < nBinsSqrtK <= 200. No histogram booked." << endl;
    return;
  }
  float sqrtKBinWidth = (maxSqrtK_ - minSqrtK_)/nBinsSqrtK_;
  if (sqrtKBinWidth <= 0) {
    cout << "Invalid sqrtK range: min(sqrtK) >= max(sqrtK). No histogram booked." << endl;
    return;
  } else if (minSqrtK_ < 0 && maxSqrtK_ > 0) {  // shift the range so that 0 is a bin edge
    float shift = sqrtKBinWidth*int((maxSqrtK_ + 0.5*sqrtKBinWidth)/sqrtKBinWidth) - maxSqrtK_;
    minSqrtK_ += shift;
    maxSqrtK_ += shift;
    minSqrtKScan_ = 0.5*sqrtKBinWidth;
    maxSqrtKScan_ = max(fabs(minSqrtK_), fabs(maxSqrtK_)) - 0.5*sqrtKBinWidth;
  } else {
    minSqrtKScan_ = min(fabs(minSqrtK_), fabs(maxSqrtK_)) + 0.5*sqrtKBinWidth;
    maxSqrtKScan_ = max(fabs(minSqrtK_), fabs(maxSqrtK_)) - 0.5*sqrtKBinWidth;
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
HoughCheckOnTracks::endJob() 
{
  TFile outRootFile("houghCheck_tracks.root", "RECREATE");
  hHoughVotes_->Write();
  trackTree_->Write();
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


// ------------ function to zero to get doca given kappa and phi  ------------

double HoughCheckOnTracks::fZDoca_(double* var, double* par)
{
  double x = par[0];
  double y = par[1];
  double r = sqrt(x*x + y*y);
  double sign = par[2];
  int signDoca = vDoca_.back()/fabs(vDoca_.back());
  double kappaM = -signDoca*vKappa_.back();
  double phiM = vPhi_.back() - signDoca*M_PI/2.;
  phiM += -2.*M_PI*(int((phiM + 3.*M_PI)/(2.*M_PI)) - 1);  // map to range (-PI, PI)
  double akappa = (2*var[0] + kappaM*var[0]*var[0] + kappaM*r*r)/(2*r);
  double hkappa = sqrt((var[0]*kappaM + 1)*(var[0]*kappaM + 1) - akappa*akappa);
  double rkappax = akappa*x - sign*hkappa*y;
  double rkappay = akappa*y + sign*hkappa*x;
  par[3] = atan2(rkappay, rkappax) - phiM;  // return parameter to check that angle is not supplementary
  return rkappay/rkappax - tan(phiM);
}


// ------------ function to zero to get kappa given doca and phi  ------------

double HoughCheckOnTracks::fZKappa_(double* var, double* par)
{
  double x = par[0];
  double y = par[1];
  double r = sqrt(x*x + y*y);
  double sign = par[2];
  double docaM = fabs(vDoca_.back());
  int signDoca = vDoca_.back()/fabs(vDoca_.back());
  double phiM = vPhi_.back() - signDoca*M_PI/2.;
  phiM += -2.*M_PI*(int((phiM + 3.*M_PI)/(2.*M_PI)) - 1);  // map to range (-PI, PI)
  double akappa = (2*docaM + var[0]*docaM*docaM + var[0]*r*r)/(2*r);
  double hkappa = sqrt((docaM*var[0] + 1)*(docaM*var[0] + 1) - akappa*akappa);
  double rkappax = akappa*x - sign*hkappa*y;
  double rkappay = akappa*y + sign*hkappa*x;
  par[3] = atan2(rkappay, rkappax) - phiM;  // return parameter to check that angle is not supplementary
  return rkappay/rkappax - tan(phiM);
}


// ------------ function to get phi given doca and kappa  ------------

double HoughCheckOnTracks::fPhi_(double x, double y, int sign)
{
  double r = sqrt(x*x + y*y);
  double docaM = fabs(vDoca_.back());
  int signDoca = vDoca_.back()/fabs(vDoca_.back());
  double kappaM = -signDoca*vKappa_.back();
  double akappa = (2*docaM + kappaM*docaM*docaM + kappaM*r*r)/(2*r);
  double hkappa = sqrt((docaM*kappaM + 1)*(docaM*kappaM + 1) - akappa*akappa);
  double rkappax = akappa*x - sign*hkappa*y;
  double rkappay = akappa*y + sign*hkappa*x;
  double phi = atan2(rkappay, rkappax) + signDoca*M_PI/2.;
  phi += -2.*M_PI*(int((phi + 3.*M_PI)/(2.*M_PI)) - 1);  // map to range (-PI, PI)
  return phi;
}



// ------------ function to get z0 given doca, kappa, phi, eta  ------------

double HoughCheckOnTracks::fZ0_(double x, double y, double z)
{
  double xc = (vDoca_.back() - 1./vKappa_.back())*sin(vPhi_.back());
  double yc = -(vDoca_.back() - 1./vKappa_.back())*cos(vPhi_.back());
  double st = 1./vKappa_.back()*(atan2(vKappa_.back()*(y - yc), vKappa_.back()*(x - xc)) - vPhi_.back() + M_PI/2.);
  if (st < 0)  // rotation angle has crossed the +/-PI border
    st += 2.*M_PI/fabs(vKappa_.back());
  else if (st > 2.*M_PI/fabs(vKappa_.back()))
    st -= 2.*M_PI/fabs(vKappa_.back());
  double lambda = M_PI/2. - vTheta_.back();
  return z - st*tan(lambda);
}


// ------------ function to get theta given doca, kappa, phi, z0  ------------

double HoughCheckOnTracks::fTheta_(double x, double y, double z)
{
  double xc = (vDoca_.back() - 1./vKappa_.back())*sin(vPhi_.back());
  double yc = -(vDoca_.back() - 1./vKappa_.back())*cos(vPhi_.back());
  double st = 1./vKappa_.back()*(atan2(vKappa_.back()*(y - yc), vKappa_.back()*(x - xc)) - vPhi_.back() + M_PI/2.);
  if (st < 0)  // rotation angle has crossed the +/-PI border
    st += 2.*M_PI/fabs(vKappa_.back());
  else if (st > 2.*M_PI/fabs(vKappa_.back()))
    st -= 2.*M_PI/fabs(vKappa_.back());
  return M_PI/2. - atan((z - vZ0_.back())/st);
}


//define this as a plug-in
DEFINE_FWK_MODULE(HoughCheckOnTracks);
