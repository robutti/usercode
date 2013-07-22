// -*- C++ -*-
//
// Package:    HoughTransChecks
// Class:      HoughCheck2Steps
// 
/**\class HoughCheck2Steps HoughCheck2Steps.cc HoughTest/HoughTransChecks/src/HoughCheck2Steps.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Enrico Robutti
//         Created:  Thu, June 19, 2013
// $Id: HoughCheck2Steps.cc
//
//


// system include files
#include <memory>
#include <vector>
#include <algorithm>
#include <map>
#include <utility>
#include <iomanip>

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

#include "TH3S.h"
#include "TH2S.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"


//
// class declaration
//

using namespace std;
using namespace edm;

class HoughCheck2Steps : public edm::EDAnalyzer {
public:
  explicit HoughCheck2Steps(const edm::ParameterSet&);
  ~HoughCheck2Steps();
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
  const edm::InputTag trackTags_; //used to select what tracks to read from configuration file
  string builderName_;
  const TransientTrackingRecHitBuilder* builder_;
  const vector<unsigned int> algoSel_;
  const string layerListName_;
  ctfseeding::SeedingLayerSets layerSets_;
  vector<double> minPar_;
  vector<double> maxPar_;
  vector<int> nBins_;
  const double phiBinOverlap_;
  const double etaBinOverlap_;
  const unsigned int xyVoteThr_;
  const unsigned int xyzVoteThr_;
  const bool cleanupSeeds_;
  const string outRootFile_;
  const unsigned int verbosity_;

  // Histogram and tree 
  auto_ptr<TH3S> hXYHoughVotes_;
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
const int HoughCheck2Steps::lyrMapOffset_[7] = {0, 0, 3, 5, 9, 12, 18};

//
// static data member definitions
//

//
// constructors and destructor
//
HoughCheck2Steps::HoughCheck2Steps(const edm::ParameterSet& iConfig)
:
  trackTags_(iConfig.getParameter<edm::InputTag>("tracks")),
  builderName_(iConfig.getParameter<std::string>("TTRHBuilder")),
  algoSel_(iConfig.getParameter<vector<unsigned int> >("algoSel")),
  layerListName_(iConfig.getParameter<std::string> ("seedingLayers")),
  minPar_(iConfig.getUntrackedParameter<vector<double> >("minPar", {-44., -0.075, -M_PI, -265., -2.5})),
  maxPar_(iConfig.getUntrackedParameter<vector<double> >("maxPar", {44., 0.075, M_PI, 265., 2.5})),
  nBins_(iConfig.getUntrackedParameter<vector<int> >("nBins", {100, 100, 100, 100, 100})),
  phiBinOverlap_(iConfig.getParameter<double>("phiBinOverlap")),
  etaBinOverlap_(iConfig.getParameter<double>("etaBinOverlap")),
  xyVoteThr_(iConfig.getParameter<unsigned int>("xyVoteThr")),
  xyzVoteThr_(iConfig.getParameter<unsigned int>("xyzVoteThr")),
  cleanupSeeds_(iConfig.getUntrackedParameter<bool>("cleanupSeeds", false)),
  outRootFile_(iConfig.getUntrackedParameter<string>("outRootFile", "houghCheck_2steps.root")),
  verbosity_(iConfig.getUntrackedParameter<unsigned int>("verbosity", 0))
{
}


HoughCheck2Steps::~HoughCheck2Steps()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HoughCheck2Steps::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using reco::TrackCollection;
   
  Handle<TrackCollection> tracks;
  iEvent.getByLabel(trackTags_,tracks);
  int ntrk = -1;
  vDoca_.clear();
  vKappa_.clear();
  vPhi_.clear();
  vZ0_.clear();
  vTheta_.clear();
  vAlgo_.clear();

  double docaBinWidth = (maxPar_[0] - minPar_[0])/nBins_[0];
  double sqrtKBinWidth = (maxPar_[1] - minPar_[1])/nBins_[1];
  double phiBinWidth = (maxPar_[2] - minPar_[2])/nBins_[2];
  double phiOverlap = phiBinOverlap_*phiBinWidth;
  double z0BinWidth = (maxPar_[3] - minPar_[3])/nBins_[3];
  double etaBinWidth = (maxPar_[4] - minPar_[4])/nBins_[4];
  double etaOverlap = etaBinOverlap_*etaBinWidth;
  
  vector<vector<TransientTrackingRecHit::RecHitPointer> > vHitsTrk;
  for(TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack) {
    ntrk++;
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

    // Store track hits in auxiliary vector
    vector<TransientTrackingRecHit::RecHitPointer> vTHits;
    for (trackingRecHit_iterator i = itTrack->recHitsBegin(); i != itTrack->recHitsEnd(); i++){
      TransientTrackingRecHit::RecHitPointer hit = builder_->build(&**i );
      if (hit->isValid()) {
	DetId hitId = hit->geographicalId();
	if(hitId.det() != DetId::Tracker)
	  continue;
	vTHits.push_back(hit);
      }
    }
    vHitsTrk.push_back(vTHits);
  }

  if (!(hXYHoughVotes_.get()))
    return;
  // Transverse step.
  // Loop on layers and hits within, store in auxiliary vector
  vector<TransientTrackingRecHit::ConstRecHitPointer> vHits;
  map<unsigned int, map<unsigned int, vector<pair<int, double> > > > xyVotes;
  int nHit = -1;
  for (ctfseeding::SeedingLayerSets::const_iterator itLS = layerSets_.begin(); itLS != layerSets_.end(); itLS++ ) {
    for (ctfseeding::SeedingLayers::const_iterator itLyr = itLS->begin(); itLyr!= itLS->end(); itLyr++ ) {
      ctfseeding::SeedingLayer::Hits hits = itLyr->hits(iEvent, iSetup);
      for (ctfseeding::SeedingLayer::Hits::const_iterator itHit = hits.begin(); itHit != hits.end(); itHit++) {
	nHit++;
	vHits.push_back(*itHit);
	if (verbosity_ > 2)
	  cout << "hit #" << nHit;
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
	  double r = sqrt(x*x + y*y);
	  // Loop on allowed histogram bins (first doca, then kappa, then phi) and fill histogram
	  for (int iDoca = 0; iDoca < nBins_[0]; iDoca++) {
	    double doca = minPar_[0] + (iDoca + 0.5)*docaBinWidth;
	    for (int iSqrtK = 0; iSqrtK < nBins_[1]; iSqrtK++) {
	      double sqrtK = minPar_[1] + (iSqrtK + 0.5)*sqrtKBinWidth;
	      double kappa = sqrtK*fabs(sqrtK);
	      // Check for allowed range -2/(r - doca) <= kappa <= 2/(r + doca)
	      if (kappa < -2/(r - doca) || kappa > 2/(r + doca))
		continue;
	      double akappa = (kappa*doca*doca + kappa*r*r - 2*doca)/(2*r);
	      double hkappa = sqrt((1 - doca*kappa)*(1 - doca*kappa) - akappa*akappa);
	      double kappaxD = -akappa*x/r + hkappa*y/r;
	      double kappayD = -akappa*y/r - hkappa*x/r;
	      double phi = atan2(kappayD, kappaxD) + M_PI/2.;
	      if (phi > M_PI)
		phi -= 2.*M_PI;  // map to range (-PI, PI)
	      // Insert point in vote map, taking overlaps into account
	      int firstPhiBin = max(0, (int)((phi - minPar_[2] - phiOverlap)/phiBinWidth));
	      int lastPhiBin = min(nBins_[2] - 1, (int)((phi - minPar_[2] + phiOverlap)/phiBinWidth));
	      for (int iPhi = firstPhiBin; iPhi <= lastPhiBin; iPhi++) {
		unsigned int keyBin = iDoca | (iSqrtK << 8) | (iPhi << 16);
		map<unsigned int, vector<pair<int, double> > >& lyrMap = xyVotes[keyBin];
		unsigned int keyLyr = (hitSubDet << 4) | hitLayer;
		vector<pair<int, double> >& vLyrHits = lyrMap[keyLyr];
		vLyrHits.push_back(make_pair(nHit, phi));
	      }
	    }
	  }
	} else 
	  if (verbosity_ > 2)
	    cout << " - invalid hit" << endl;
      }
    }
  }

  // Longitudinal step
  map<unsigned long, map<unsigned int, vector<int> > > xyzVotes;
  // Fill 5-par map for 3-par bins over threshold
  for (map<unsigned int, map<unsigned int, vector<pair<int, double> > > >::iterator itBin = xyVotes.begin(); itBin != xyVotes.end(); itBin++) {
    unsigned int keyBin = (*itBin).first;
    int docaBin = (keyBin & 255) + 1;
    int sqrtKBin = ((keyBin >> 8) & 255) + 1;
    int phiBin = ((keyBin >> 16) & 255) + 1;
    unsigned int lyrVotes = (*itBin).second.size();
    // Fill x-y histogram with number of voting layers
    hXYHoughVotes_->SetBinContent(docaBin, sqrtKBin, phiBin, lyrVotes);
    if (lyrVotes > xyVoteThr_) {
      double doca = minPar_[0] + (docaBin + 0.5)*docaBinWidth;
      double sqrtK = minPar_[1] + (sqrtKBin + 0.5)*sqrtKBinWidth;
      double kappa = sqrtK*fabs(sqrtK);
      for (map<unsigned int, vector<pair<int, double> > >::iterator itLyr = (*itBin).second.begin(); itLyr != (*itBin).second.end(); itLyr++) {
	unsigned int keyLyr = (*itLyr).first;
	for (vector<pair<int, double> >::iterator itHit = (*itLyr).second.begin(); itHit != (*itLyr).second.end(); itHit++) {
	  TransientTrackingRecHit::ConstRecHitPointer hit = vHits[(*itHit).first];
	  GlobalPoint hitPosition = hit->globalPosition();
	  double x = 10.*hitPosition.x();
	  double y = 10.*hitPosition.y();
	  double z = 10.*hitPosition.z();
	  double phi = (*itHit).second;
	  double xc = (doca - 1./kappa)*sin(phi);
	  double yc = -(doca - 1./kappa)*cos(phi);
	  for (int iZ0 = 0; iZ0 < nBins_[3]; iZ0++) {
	    double z0 = minPar_[3] + (iZ0 + 0.5)*z0BinWidth;
	    double st = 1./kappa*(atan2(kappa*(y - yc), kappa*(x - xc)) - phi + M_PI/2.);
	    if (st < 0)  // rotation angle has crossed the +/-PI border
	      st += 2.*M_PI/fabs(kappa);
	    else if (st > 2.*M_PI/fabs(kappa))
	      st -= 2.*M_PI/fabs(kappa);
	    double lambda = atan((z - z0)/st);
	    double eta = -log(tan((M_PI/2. - lambda)/2.));
	    // Insert point in 5-parameter vote map, taking overlaps into account
	    int firstEtaBin = max(0, (int)((eta - minPar_[4] - etaOverlap)/etaBinWidth));
	    int lastEtaBin = min(nBins_[4] - 1, (int)((eta - minPar_[4] + etaOverlap)/etaBinWidth));
	    for (int iEta = firstEtaBin; iEta <= lastEtaBin; iEta++) {
	      unsigned long key5Par = keyBin | (iZ0 << 24) | ((unsigned long)iEta << 32);
	      map<unsigned int, vector<int> >& lyrMap = xyzVotes[key5Par];
	      vector<int>& vLyrHits = lyrMap[keyLyr];
	      vLyrHits.push_back((*itHit).first);
	    }
	  }
	}
      }      
    }
  }

  // Clean-up seeds: remove below threshold and duplicates if requested
  map<unsigned long, map<unsigned int, vector<int> > > passVotes;
  if (verbosity_ > 0)
    cout << "Number of seeds before clean-up = " << xyzVotes.size() << endl;
  for (map<unsigned long, map<unsigned int, vector<int> > >::iterator itBin = xyzVotes.begin(); itBin != xyzVotes.end(); itBin++) {
    bool nextBin = false;
    unsigned int lyrVotes = (*itBin).second.size();
    if (lyrVotes > xyzVoteThr_) {
      if (cleanupSeeds_) {
	// Get bin
	unsigned long htBin = (*itBin).first;
	//	cout << "* DEBUG: htBin = 0x" << hex <<htBin << dec << endl;  // debug
	int docaBin = htBin & 255;
	int sqrtKBin = (htBin >> 8) & 255;
	int phiBin = (htBin >> 16) & 255;
	int z0Bin = (htBin >> 24) & 255;
	int etaBin = (htBin >> 32) & 255;
	int deltaBin = 1;  // width of bin range where to look for duplicates
	for (int iDoca = max(0, docaBin - deltaBin); !nextBin && iDoca < min(nBins_[0], docaBin + deltaBin + 1); iDoca++)
	  for (int iSqrtK = max(0, sqrtKBin - deltaBin); !nextBin && iSqrtK < min(nBins_[1], sqrtKBin + deltaBin + 1); iSqrtK++)
	    for (int iPhi = max(0, phiBin - deltaBin); !nextBin && iPhi < min(nBins_[2], phiBin + deltaBin + 1); iPhi++)
	      for (int iZ0 = max(0, z0Bin - deltaBin); !nextBin && iZ0 < min(nBins_[3], z0Bin + deltaBin + 1); iZ0++)
		for (int iEta = max(0, etaBin - deltaBin); !nextBin && iEta < min(nBins_[4], etaBin + deltaBin + 1); iEta++) {
		  unsigned long key = iDoca | (iSqrtK << 8) | (iPhi << 16) | ((unsigned long)iZ0 << 24) | ((unsigned long)iEta << 32);
		  //		  cout << "  key = 0x" << hex << key << dec << endl;
		  if (key == htBin || xyzVotes.find(key) == xyzVotes.end())
		    continue;
		  bool hitFound = true;
		  for (map<unsigned int, vector<int> >::iterator itLyr1 = (*itBin).second.begin(); hitFound && itLyr1 != (*itBin).second.end(); itLyr1++) {
		    for (vector<int>::iterator itHit1 = (*itLyr1).second.begin(); hitFound && itHit1 != (*itLyr1).second.end(); itHit1++) {
		      hitFound = false;
		      for (map<unsigned int, vector<int> >::iterator itLyr2 = xyzVotes[key].begin(); !hitFound && itLyr2 != xyzVotes[key].end(); itLyr2++) {
			if ((*itLyr2).first != (*itLyr1).first)
			  continue;
			for (vector<int>::iterator itHit2 = (*itLyr2).second.begin(); itHit2 != (*itLyr2).second.end(); itHit2++) {
			  if ((*itHit1) == (*itHit2)) {
			    hitFound = true;
			    break;
			  }
			}
		      }
		    }
		  }
		  if (hitFound) {  // all hits found in a "nearby" bin
		    (*itBin).second.clear();  // remove all hits from the bin
		    nextBin = true;
		  }
		}
      }
      if (!nextBin)  // no nearby bin containing all hits (or no clean-up)
	passVotes[(*itBin).first] = (*itBin).second;
    }
  }
  unsigned int nSeeds = passVotes.size();
  if (verbosity_ > 0)
    cout << "Number of seeds after clean-up = " << nSeeds << endl;

  // Compare content of 5-par bin over thresholds with track hits
  // Vector with tracks found in HT histogram
  vector<bool> trackFound(vHitsTrk.size(), false);
  // Number of HT "seeds" above threshold and "seeds" corresponding to tracks
  unsigned int nGoodSeeds = 0;
  for (map<unsigned long, map<unsigned int, vector<int> > >::iterator itBin = passVotes.begin(); itBin != passVotes.end(); itBin++) {
    vector<map<unsigned int, bool> > lyrFound(vHitsTrk.size());  // vector with layers with association found
    for (map<unsigned int, vector<int> >::iterator itLyr = (*itBin).second.begin(); itLyr != (*itBin).second.end(); itLyr++) {
      unsigned int lyrId = (*itLyr).first;
      for (vector<int>::iterator itHit = (*itLyr).second.begin(); itHit != (*itLyr).second.end(); itHit++) {
	DetId hitId = vHits[(*itHit)]->geographicalId();
	GlobalPoint hitPosition = vHits[(*itHit)]->globalPosition();
	for (unsigned int iTrk = 0; iTrk < vHitsTrk.size(); iTrk++){
	  for (vector<TransientTrackingRecHit::RecHitPointer>::iterator itTrkHit = vHitsTrk[iTrk].begin(); itTrkHit != vHitsTrk[iTrk].end(); itTrkHit++) {
	    GlobalPoint trkHitPosition = (*itTrkHit)->globalPosition();
	    DetId trkHitId = (*itTrkHit)->geographicalId();
	    if ((hitId.rawId() & (~3)) == (trkHitId.rawId() & ~3)) {
	      double hitPhi = atan2(hitPosition.y(), hitPosition.x());
	      double trkHitPhi = atan2(trkHitPosition.y(), trkHitPosition.x());
	      double hitR = sqrt(hitPosition.x()*hitPosition.x() + hitPosition.y()*hitPosition.y());
	      switch (trkHitId.subdetId()) {
	      case PixelSubdetector::PixelBarrel:
		// Insert appropriate condition
		continue;
		break;
	      case PixelSubdetector::PixelEndcap:
		// Insert appropriate condition
		continue;
		break;
	      case StripSubdetector::TIB:
		if (((TIBDetId)hitId).isDoubleSide() && ((TIBDetId)trkHitId).isRPhi()) {
		  if (fabs(hitR*(hitPhi - trkHitPhi) > 0.02))
		    continue;
		} else
		  // Add other cases
		  continue;
		break;
	      case StripSubdetector::TID:
		if (((TIDDetId)hitId).isDoubleSide() && ((TIDDetId)trkHitId).isRPhi()) {
		  if (fabs(hitR*(hitPhi - trkHitPhi) > 0.02))
		    continue;
		} else
		  // Add other cases
		  break;
	      case StripSubdetector::TOB:
		if (((TOBDetId)hitId).isDoubleSide() && ((TOBDetId)trkHitId).isRPhi()) {
		  if (fabs(hitR*(hitPhi - trkHitPhi) > 0.02))
		    continue;
		} else
		  // Add other cases
		  continue;
		break;
	      case StripSubdetector::TEC:
		if (((TECDetId)hitId).isDoubleSide() && ((TECDetId)trkHitId).isRPhi()) {
		  if (fabs(hitR*(hitPhi - trkHitPhi) > 0.02))
		    continue;
		} else
		  // Add other cases
		  continue;
		break;
	      default:
		break;
	      }
	      if (lyrFound[iTrk].find(lyrId) == lyrFound[iTrk].end())  // new layer associated
		(lyrFound[iTrk])[lyrId] = true;
	      break;
	    }
	  }
	}
      }
    }
    // Loop again on tracks to check whether some has at least three hits associated to HT votes in bin
    for (unsigned int iTrk = 0; iTrk < vHitsTrk.size(); iTrk++) {
      if ((lyrFound[iTrk]).size() >= 3) {
	trackFound[iTrk] = true;
	nGoodSeeds++;
	break;
      }
    }
  }

  // Final loop on tracks to count how many were "found" by HT
  unsigned int nFoundTracks = 0;
  for (unsigned int iTrk = 0; iTrk < vHitsTrk.size(); iTrk++)
    if (trackFound[iTrk])
      nFoundTracks++;
  // Print results
  cout << "Total number of seeds: " << nSeeds << endl;
  cout << "Seeds associated to tracks: " << nGoodSeeds << endl;
  cout << "Total number of tracks: " << trackFound.size() << endl;
  cout << "Tracks associated to seeds: " << nFoundTracks << endl;

  // Fill track tree
  trackTree_->Fill();
}


// ------------ method called once each job just before starting event loop  ------------
void 
HoughCheck2Steps::beginJob()
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

  // Check parameter ranges and binnings
  for (int iPar = 0; iPar < 5; iPar++) {
    if (nBins_[iPar] <= 0 || nBins_[iPar] > 256) {
      cout << "Invalid nBins[" << iPar << "] parameter (" << nBins_[iPar] << "). Valid range is 0 < nBins <= 256. No histogram booked." << endl;
      return;
    }
    if (minPar_[iPar] >= maxPar_[iPar]) {
      cout << "Invalid parameter " << iPar << " range: min(par" << iPar << ") >= max(par" << iPar << "). No histogram booked." << endl;
    return;
    }
  }
  if (minPar_[2] < -M_PI || maxPar_[2] > M_PI) {
    cout << "Invalid phi range. Valid range is -pi <= phi <= pi. No histogram booked." << endl;
    return;
  }

  // Book 3D histogram
  hXYHoughVotes_.reset(new TH3S("hXYHoughVotes", "helix Hough transform votes", nBins_[0], minPar_[0], maxPar_[0], nBins_[1], minPar_[1], maxPar_[1], nBins_[2], minPar_[2], maxPar_[2]));
  hXYHoughVotes_->SetDirectory(0);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HoughCheck2Steps::endJob() 
{
  TFile outRootFile(outRootFile_.c_str(), "RECREATE");
  hXYHoughVotes_->Write();
  trackTree_->Write();
}

// ------------ method called when starting to processes a run  ------------
void 
HoughCheck2Steps::beginRun(edm::Run const& run, edm::EventSetup const& setup)
{
  // TTRH builder
  edm::ESHandle<TransientTrackingRecHitBuilder> theBuilder;
  setup.get<TransientRecHitRecord>().get(builderName_,theBuilder);
  builder_=theBuilder.product();

  // Seeding layers
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
HoughCheck2Steps::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HoughCheck2Steps::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HoughCheck2Steps::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HoughCheck2Steps::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<std::string>("TTRHBuilder");
  desc.add<vector<unsigned int> >("algoSel", {});
  desc.add<std::string>("seedingLayers", "HoughTransformSeedLayersAllHitsOneSet");
  desc.addUntracked<vector<double> >("minPar", {-44., -0.075, -M_PI, -265., -2.5});
  desc.addUntracked<vector<double> >("maxPar", {44., 0.075, M_PI, 265., 2.5});
  desc.addUntracked<vector<int> >("nBins", {100, 100, 100, 100, 100});
  desc.add<double>("phiBinOverlap", 0.);
  desc.add<double>("etaBinOverlap", 0.);
  desc.add<unsigned int>("xyVoteThr", 0);
  desc.add<unsigned int>("xyzVoteThr", 0);
  desc.addUntracked<bool>("cleanupSeeds", false);
  desc.addUntracked<string>("outRootFile", "houghCheck_2steps.root");
  desc.addUntracked<unsigned int>("verbosity", 0);
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HoughCheck2Steps);
