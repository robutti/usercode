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

#include "ERobutti/HoughTransChecks/interface/TH5.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
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

  void plotHoughProjections(int iparX, int iparY);
  
  // Parameters
  edm::InputTag trackTags_; //used to select what tracks to read from configuration file
  std::string builderName_;
  const TransientTrackingRecHitBuilder* builder_;
  const vector<unsigned int> algoSel_;
  const vector<unsigned int> layerSel_;
  double minDoca_, maxDoca_;
  const int nBinsDoca_;
  double minDocaScan_, maxDocaScan_;
  double minKappa_, maxKappa_;
  const int nBinsKappa_;
  double minKappaScan_, maxKappaScan_;
  double minPhi_, maxPhi_;
  const int nBinsPhi_;
  double minZ0_, maxZ0_;
  const int nBinsZ0_;
  double minLambda_, maxLambda_;
  const int nBinsLambda_;
  const int projectionPars_;
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
  vector<double> vAlgo_;
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
  nBinsKappa_(iConfig.getUntrackedParameter<int>("nBinsKappa", 100)),
  nBinsPhi_(iConfig.getUntrackedParameter<int>("nBinsPhi", 100)),
  nBinsZ0_(iConfig.getUntrackedParameter<int>("nBinsZ0", 100)),
  nBinsLambda_(iConfig.getUntrackedParameter<int>("nBinsLambda", 100)),
  projectionPars_(iConfig.getUntrackedParameter<int>("projectionPars", 0)),
  verbosity_(iConfig.getUntrackedParameter<unsigned int>("verbosity", 0)),
  printProgressFrequency_(iConfig.getUntrackedParameter<unsigned int>("printProgressFrequency", 1000))
{
  vector<double> rangeDoca = iConfig.getUntrackedParameter<vector<double> >("rangeDoca", {-44., 44.});  // radius of inner pixel layer
  minDoca_ = rangeDoca[0];
  maxDoca_ = rangeDoca[1];
  vector<double> rangeKappa = iConfig.getUntrackedParameter<vector<double> >("rangeKappa", {-1./45, 1./45});  // larger than radius of inner pixel layer
  minKappa_ = rangeKappa[0];
  maxKappa_ = rangeKappa[1];
  vector<double> rangePhi = iConfig.getUntrackedParameter<vector<double> >("rangePhi", {-M_PI, M_PI});
  minPhi_ = rangePhi[0];
  maxPhi_ = rangePhi[1];
  vector<double> rangeZ0 = iConfig.getUntrackedParameter<vector<double> >("rangeZ0", {-265., 265.});  // pixel barrel half-length
  minZ0_ = rangeZ0[0];
  maxZ0_ = rangeZ0[1];
  vector<double> rangeLambda = iConfig.getUntrackedParameter<vector<double> >("rangeLambda", {-1.40, 1.40});  // tracker eta acceptance
  minLambda_ = rangeLambda[0];
  maxLambda_ = rangeLambda[1];
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
  vLambda_.clear();
  for(TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack) {
    if (verbosity_ > 1)
      cout << "Track #" << ntrk << " with q=" << itTrack->charge() 
	   << ", pT=" << itTrack->pt() << " GeV, eta=" << itTrack->eta() 
	   << ", Nhits=" << itTrack->recHitsSize()
	   << ", (vx,vy,vz)=(" << itTrack->vx() << "," << itTrack->vy() << "," << itTrack->vz() << ")"
	   << ", doca=" << 10.*(itTrack->d0())
	   << ", kappa=" << -1.139e-3*itTrack->charge()/(itTrack->pt())
	   << ", phi=" << itTrack->phi()
	   << ", z0=" << 10.*(itTrack->dz())
	   << ", lambda=" << itTrack->lambda()
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
    vDoca_.push_back(10.*(itTrack->d0()));
    vKappa_.push_back(-1.139e-3*itTrack->charge()/(itTrack->pt()));
    vPhi_.push_back(itTrack->phi());
    vZ0_.push_back(10.*(itTrack->dz()));
    vLambda_.push_back(itTrack->lambda());
    vAlgo_.push_back(itTrack->algo());
    if (!(hHoughVotes_.get()))
      return;
    int nhit = 0;
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
	double y = 10.*hitPosition.y();
	double z = 10.*hitPosition.z();
	double r = sqrt(x*x + y*y);
	double docaBinWidth = (maxDoca_ - minDoca_)/nBinsDoca_;
	double kappaBinWidth = (maxKappa_ - minKappa_)/nBinsKappa_;
	double z0BinWidth = (maxZ0_ - minZ0_)/nBinsZ0_;
	// Loop on allowed histogram bins (first doca, then kappa, then phi) and fill histogram
	for (double docaScan = minDocaScan_; docaScan < maxDocaScan_ + 0.5*docaBinWidth; docaScan += docaBinWidth) {
	  // Start from first allowed kappa bin (must be kappaScan >= -2/(r + docaScan))
	  double firstKappa = -maxKappaScan_;
	  if (-2/(r + docaScan) > firstKappa) {
	    int binPosition = (-2/(r + docaScan))/kappaBinWidth;
	    firstKappa = (binPosition - 0.5)*kappaBinWidth;
	  }
	  for (double kappaScan = firstKappa; kappaScan < min(maxKappaScan_ + 0.5*kappaBinWidth, 2/(r - docaScan)); kappaScan += kappaBinWidth) {
	    if (kappaScan > -minKappaScan_ + 0.5*kappaBinWidth && kappaScan < minKappaScan_ - 0.5*kappaBinWidth)  // skip irrelevant values for histogram
	      continue;
	    double akappa = (2*docaScan + kappaScan*docaScan*docaScan + kappaScan*r*r)/(2*r);
	    double hkappa = sqrt((docaScan*kappaScan + 1)*(docaScan*kappaScan + 1) - akappa*akappa);
	    for (int sign = -1; sign <= 1; sign += 2) {
	      double kappax = akappa*x/r - sign*hkappa*y/r;
	      double kappay = akappa*y/r + sign*hkappa*x/r;
	      // Convert to perigee parameters
	      double phiHit = atan2(y, x);
	      double phiD = atan2(kappay, kappax);
	      int rot = -2*(int((phiHit - phiD + 2*M_PI)/M_PI)%2) + 1;  // hit position wrt. poca: +1 = anticlockwise; -1 = clockwise
	      double doca = rot*docaScan;
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
// 		  cout << "(doca, kappa, phi, z0, psi, st) = "
// 		       << "(" << doca << ", " << kappa << ", " << phi << ", " << z0 << ", " << psi << ", " << st << ")" << endl;
		double lambda = atan((z - z0)/st);
// 		cout << "(doca, kappa, phi, z0, lambda) = "
// 		     << "(" << doca << ", " << kappa << ", " << phi << ", " << z0 << ", " << lambda << ")" << endl;
		hHoughVotes_->Fill(doca, kappa, phi, z0, lambda);
	      }
	    }
	  }
	}
      } else 
	if (verbosity_ > 2)
	  cout << " - invalid hit" << endl;
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
  // Create tree, with branches
  trackTree_.reset(new TTree("trackTree", "Fitted track parameters"));
  trackTree_->Branch("doca", &vDoca_);
  trackTree_->Branch("kappa", &vKappa_);
  trackTree_->Branch("phi", &vPhi_);
  trackTree_->Branch("z0", &vZ0_);
  trackTree_->Branch("lambda", &vLambda_);
  trackTree_->Branch("algo", &vAlgo_);
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
  if (nBinsKappa_ <= 0 || nBinsKappa_ > 200) {
    cout << "Invalid nBinsKappa parameter (" << nBinsKappa_ << "). Valid range is 0 < nBinsKappa <= 200. No histogram booked." << endl;
    return;
  }
  float kappaBinWidth = (maxKappa_ - minKappa_)/nBinsKappa_;
  if (kappaBinWidth <= 0) {
    cout << "Invalid kappa range: min(kappa) >= max(kappa). No histogram booked." << endl;
    return;
  } else if (minKappa_ < 0 && maxKappa_ > 0) {  // shift the range so that 0 is a bin edge
    float shift = kappaBinWidth*int((maxKappa_ + 0.5*kappaBinWidth)/kappaBinWidth) - maxKappa_;
    minKappa_ += shift;
    maxKappa_ += shift;
    minKappaScan_ = 0.5*kappaBinWidth;
    maxKappaScan_ = max(fabs(minKappa_), fabs(maxKappa_)) - 0.5*kappaBinWidth;
  } else {
    minKappaScan_ = min(fabs(minKappa_), fabs(maxKappa_)) + 0.5*kappaBinWidth;
    maxKappaScan_ = max(fabs(minKappa_), fabs(maxKappa_)) - 0.5*kappaBinWidth;
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
  if (nBinsLambda_ <= 0 || nBinsLambda_ > 200) {
    cout << "Invalid nBinsLambda parameter (" << nBinsLambda_ << "). Valid range is 0 < nBinsLambda <= 200. No histogram booked." << endl;
    return;
  }
  if (minLambda_ >= maxLambda_) {
    cout << "Invalid lambda range: min(lambda) >= max(lambda). No histogram booked." << endl;
    return;
  }
  // Book 5D histogram
  hHoughVotes_.reset(new TH5C("hHoughVotes", "helix Hough transform votes", nBinsDoca_, minDoca_, maxDoca_, nBinsKappa_, minKappa_, maxKappa_,
			      nBinsPhi_, minPhi_, maxPhi_, nBinsZ0_, minZ0_, maxZ0_, nBinsLambda_, minLambda_, maxLambda_));
  hHoughVotes_->SetDirectory(0);

  // Set content to -128 for all histogram bins, to use full 8-bit range
  for (int iDoca = 1; iDoca <= nBinsDoca_; iDoca++)
    for (int iKappa = 1; iKappa <= nBinsKappa_; iKappa++)
      for (int iPhi = 1; iPhi <= nBinsPhi_; iPhi++)
	for (int iZ0 = 1; iZ0 <= nBinsZ0_; iZ0++)
	  for (int iLambda = 1; iLambda <= nBinsLambda_; iLambda++)
	    hHoughVotes_->SetBinContent(iDoca, iKappa, iPhi, iZ0, iLambda, -128);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HoughCheckOnTracks::endJob() 
{
  int iparX = projectionPars_/10;
  int iparY = projectionPars_%10;
  if (iparX >= 0 && iparX <= 4 && iparY >= 0 && iparY <= 4 && iparX != iparY)
    plotHoughProjections(iparX, iparY);
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

// ------------ method to plot projection histograms for Hough votes  ------------
void
HoughCheckOnTracks::plotHoughProjections(int iparX, int iparY) {
  // Build all-events vectors from tree (per-event) entries
  vector<double> vPar[5];
  // Loop on tree
  for (int iEv = 0; iEv < trackTree_->GetEntries(); iEv++) {
    long tEntry = trackTree_->LoadTree(iEv);
    trackTree_->GetEntry(tEntry);
    for (unsigned int iTk = 0; iTk < vDoca_.size(); iTk++) {
      vPar[0].push_back(vDoca_.at(iTk));
      vPar[1].push_back(vKappa_.at(iTk));
      vPar[2].push_back(vPhi_.at(iTk));
      vPar[3].push_back(vZ0_.at(iTk));
      vPar[4].push_back(vLambda_.at(iTk));
    }
  }

  // Index names, axes, etc.
  TString axisName[5];
  axisName[0] = "x";
  axisName[1] = "y";
  axisName[2] = "z";
  axisName[3] = "u";
  axisName[4] = "v";
  TString parName[5];
  parName[0] = "doca";
  parName[1] = "kappa";
  parName[2] = "phi";
  parName[3] = "z0";
  parName[4] = "lambda";
  TAxis* axis[5];
  axis[0] = hHoughVotes_->GetXaxis();
  axis[1] = hHoughVotes_->GetYaxis();
  axis[2] = hHoughVotes_->GetZaxis();
  axis[3] = hHoughVotes_->GetUaxis();
  axis[4] = hHoughVotes_->GetVaxis();
  int iProjAxis[3];
  unsigned int nProj = 0;
  for (int iAxis = 0; iAxis < 5; iAxis++)
    if (iAxis != iparX && iAxis != iparY) {
      iProjAxis[nProj++] = iAxis;
      if (nProj == 3)
	break;
    }

  // Set up graphics
  TString psFileName = "houghVoteProjections.ps";
  TString openPsStr = psFileName + "[";
  TString closePsStr = psFileName + "]";
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0001);
  gStyle->SetOptStat(0);
  double maxVote = hHoughVotes_->GetMaximum() + 128;
  TCanvas cHoughProj("cHoughProj", "", 800, 1200);
  TH2S* hProj[6];
  TString hName, hTitle;
  float rX = 0.05*(axis[iparX]->GetXmax() - axis[iparX]->GetXmin());
  float rY = 0.05*(axis[iparY]->GetXmax() - axis[iparY]->GetXmin());
  cHoughProj.Print(openPsStr);

  //
  TString projName = axisName[iparY] + axisName[iparX];
  int iPad = 0;
  cHoughProj.Divide(2, 3);
  for (int iBin0 = 1; iBin0 <= axis[iProjAxis[0]]->GetNbins(); iBin0++) {
    float par0 = axis[iProjAxis[0]]->GetBinCenter(iBin0);
    stringstream par0Str;
    par0Str << setprecision(4) << par0;
    for (int iBin1 = 1; iBin1 <= axis[iProjAxis[1]]->GetNbins(); iBin1++) {
      float par1 = axis[iProjAxis[1]]->GetBinCenter(iBin1);
      stringstream par1Str;
      par1Str << setprecision(4) << par1;
      for (int iBin2 = 1; iBin2 <= axis[iProjAxis[2]]->GetNbins(); iBin2++) {
	cHoughProj.cd(++iPad);
	float par2 = axis[iProjAxis[2]]->GetBinCenter(iBin2);
	stringstream par2Str;
	par2Str << setprecision(4) << par2;
	hName = TString("h") + parName[iparX] + parName[iparY] + 
	  "-" + parName[iProjAxis[0]] + par0Str.str().c_str() + 
	  "-" + parName[iProjAxis[1]] + par1Str.str().c_str() + 
	  "-" + parName[iProjAxis[2]] + par2Str.str().c_str();
	hTitle = parName[iparY] + " vs. " + parName[iparX] + ", " + 
	  parName[iProjAxis[0]] + " = " + par0Str.str().c_str() + ", " +
	  parName[iProjAxis[1]] + " = " + par1Str.str().c_str() + ", " +
	  parName[iProjAxis[2]] + " = " + par2Str.str().c_str();
	axis[iProjAxis[0]]->SetRange(iBin0, iBin0);
	axis[iProjAxis[1]]->SetRange(iBin1, iBin1);
	axis[iProjAxis[2]]->SetRange(iBin2, iBin2);
	hProj[iPad - 1] = (TH2S*)(hHoughVotes_->Project5D(projName));
	hProj[iPad - 1]->SetNameTitle(hName, hTitle);
	// Add 128 to each projection bin
	for (int iBinX = 1; iBinX <= hProj[iPad - 1]->GetNbinsX(); iBinX++)
	  for (int iBinY = 1; iBinY <= hProj[iPad - 1]->GetNbinsY(); iBinY++)
	    hProj[iPad - 1]->SetBinContent(iBinX, iBinY, hProj[iPad - 1]->GetBinContent(iBinX, iBinY) + 128);
	hProj[iPad - 1]->SetMinimum(0.);
	hProj[iPad - 1]->SetMaximum(maxVote);
	hProj[iPad - 1]->Draw("COLZ");
	// Draw circles around fitted track parameters
	for (unsigned int iTkPar = 0; iTkPar < vPar[iProjAxis[0]].size(); iTkPar++) {
	  if (vPar[iProjAxis[0]].at(iTkPar) > axis[iProjAxis[0]]->GetBinLowEdge(iBin0) &&
	      vPar[iProjAxis[0]].at(iTkPar) <= axis[iProjAxis[0]]->GetBinUpEdge(iBin0) &&
	      vPar[iProjAxis[1]].at(iTkPar) > axis[iProjAxis[1]]->GetBinLowEdge(iBin1) &&
	      vPar[iProjAxis[1]].at(iTkPar) <= axis[iProjAxis[1]]->GetBinUpEdge(iBin1) &&
	      vPar[iProjAxis[2]].at(iTkPar) > axis[iProjAxis[2]]->GetBinLowEdge(iBin2) &&
	      vPar[iProjAxis[2]].at(iTkPar) <= axis[iProjAxis[2]]->GetBinUpEdge(iBin2)) {
	    TEllipse* circle = new TEllipse(vPar[iparX].at(iTkPar), vPar[iparY].at(iTkPar), rX, rY);
	    circle->SetLineWidth(2);
	    circle->SetLineColor(kYellow);
	    circle->SetFillStyle(0);
	    circle->Draw();
	  }
	}
	// Print page if complete      
	if (iPad == 6) {
	  cHoughProj.Update();
	  cHoughProj.Print(psFileName);
	  cHoughProj.Clear();
	  cHoughProj.Divide(2, 3);
	  iPad = 0;
	}
      }
    }
  }
  cHoughProj.Update();
  cHoughProj.Print(psFileName);

  cHoughProj.Print(closePsStr);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HoughCheckOnTracks);
