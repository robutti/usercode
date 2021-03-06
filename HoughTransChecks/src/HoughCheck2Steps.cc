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
#include "TH1I.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include "ERobutti/HoughTransChecks/interface/LayerSequence.h"

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
  auto_ptr<TH1I> hNVotes_;
  auto_ptr<TH1I> hNVotesMatched_;
  auto_ptr<TTree> trackTree_;
  vector<double> vDoca_;
  vector<double> vKappa_;
  vector<double> vPhi_;
  vector<double> vZ0_;
  vector<double> vTheta_;
  vector<double> vAlgo_;
  unsigned int nSeeds_;
  unsigned int nGoodSeeds_;
  unsigned int nTriedSeeds_;
  unsigned int nAssSeeds_;
  unsigned int nTracks_;
  unsigned int nTracksAlgo_;
  unsigned int nAccTracks_;
  unsigned int nAssTracks_;
  unsigned int nFoundTracks_;

  // Allowed Layer sequences
  LayerSequence layerSeq_;
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
  verbosity_(iConfig.getUntrackedParameter<unsigned int>("verbosity", 0)),
  layerSeq_(iConfig.getParameter<unsigned int>("xyzVoteThr"))
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
   
  double docaBinWidth = (maxPar_[0] - minPar_[0])/nBins_[0];
  double sqrtKBinWidth = (maxPar_[1] - minPar_[1])/nBins_[1];
  double phiBinWidth = (maxPar_[2] - minPar_[2])/nBins_[2];
  double phiOverlap = phiBinOverlap_*phiBinWidth;
  double z0BinWidth = (maxPar_[3] - minPar_[3])/nBins_[3];
  double etaBinWidth = (maxPar_[4] - minPar_[4])/nBins_[4];
  double etaOverlap = etaBinOverlap_*etaBinWidth;
  
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
	  cout << "hit #" << nHit << ":";
	DetId hitId = (*itHit)->geographicalId();
	if (hitId.det() != DetId::Tracker)
	  continue;
	int hitSubDet = hitId.subdetId();
	int hitLayer = 0;
	if (hitSubDet == PixelSubdetector::PixelBarrel) {
	  if (verbosity_ > 2)
	    cout << "hitId = " << PXBDetId(hitId);
          hitLayer = PXBDetId(hitId).layer();
	}
	else if (hitSubDet == PixelSubdetector::PixelEndcap) {
	  if (verbosity_ > 2)
	    cout << "hitId = " << PXFDetId(hitId);
	  hitLayer = PXFDetId(hitId).disk();
	}
	else if (hitSubDet == StripSubdetector::TIB) {
	  if (verbosity_ > 2)
	    cout << "hitId = " << TIBDetId(hitId);  
	  hitLayer = TIBDetId(hitId).layer();
	}
	else if (hitSubDet == StripSubdetector::TID) {
	  if (verbosity_ > 2)
	    cout << "hitId = " << TIDDetId(hitId);
	  hitLayer = TIDDetId(hitId).wheel();
	}
	else if (hitSubDet == StripSubdetector::TOB) {
	  if (verbosity_ > 2)
	    cout << "hitId = " << TOBDetId(hitId); 
 	  hitLayer = TOBDetId(hitId).layer();
	}
	else if (hitSubDet == StripSubdetector::TEC) {
	  if (verbosity_ > 2)
	    cout << "hitId = " << TECDetId(hitId);
         hitLayer = TECDetId(hitId).wheel();
	}
	else 
	  continue;
 	if ((*itHit)->isValid()) {
	  GlobalPoint hitPosition = (*itHit)->globalPosition();
	  if (verbosity_ > 2)
	    cout << ", position = " << hitPosition << endl;
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
    if (layerSeq_.testPattern((*itBin).second)) {
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
  nSeeds_ = xyzVotes.size();

  // Clean-up seeds: remove below threshold and duplicates if requested
  map<unsigned long, map<unsigned int, vector<int> > > passVotes;
  for (map<unsigned long, map<unsigned int, vector<int> > >::iterator itBin = xyzVotes.begin(); itBin != xyzVotes.end(); itBin++) {
    bool nextBin = false;
    if (layerSeq_.testPattern((*itBin).second)) {
      if (cleanupSeeds_) {
	// Get bin
	unsigned long htBin = (*itBin).first;
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
  nGoodSeeds_ = passVotes.size();

  // Associate track hits and get parameters
  Handle<TrackCollection> tracks;
  iEvent.getByLabel(trackTags_,tracks);
  nTracks_ = 0;
  nTracksAlgo_ = 0;
  nAccTracks_ = 0;
  vDoca_.clear();
  vKappa_.clear();
  vPhi_.clear();
  vZ0_.clear();
  vTheta_.clear();
  vAlgo_.clear();
  vector<vector<int> > vHitsTrk;
  for(TrackCollection::const_iterator itTrack = tracks->begin(); itTrack != tracks->end(); ++itTrack) {
    nTracks_++;
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
      cout << "Track #" << nTracks_ - 1 << " with q=" << itTrack->charge() 
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
    nTracksAlgo_++;
    vDoca_.push_back(10.*perigeePars.transverseImpactParameter());
    vKappa_.push_back(perigeePars.transverseCurvature()/10.);
    vPhi_.push_back(perigeePars.phi());
    vZ0_.push_back(10.*perigeePars.longitudinalImpactParameter());
    vTheta_.push_back(perigeePars.theta());
    vAlgo_.push_back(itTrack->algo());
    // Check if track is in histogram acceptance
    if (vDoca_.back() < minPar_[0] || vDoca_.back() > maxPar_[0] ||
        vKappa_.back()/sqrt(fabs(vKappa_.back())) < minPar_[1] || vKappa_.back()/sqrt(fabs(vKappa_.back())) > maxPar_[1] ||
        vPhi_.back() < minPar_[2] || vPhi_.back() > maxPar_[2] ||
        vZ0_.back() < minPar_[3] || vZ0_.back() > maxPar_[3] ||
        -log(tan(0.5*vTheta_.back())) < minPar_[4] || -log(tan(0.5*vTheta_.back())) > maxPar_[4])
      continue;
    nAccTracks_++;      
    vector<int> vTHits;
    unsigned int nHitsOrig = 0;
    for (trackingRecHit_iterator i = itTrack->recHitsBegin(); i != itTrack->recHitsEnd(); i++){
      nHitsOrig++;
      int assHit = -1;
      double maxPixelDxy = 0.01;
      double maxPixelDz = 0.01;
      double maxStripDrphi = 0.02;
      double deltaZ = 1000;
      TransientTrackingRecHit::RecHitPointer trkHit = builder_->build(&**i );
      if (trkHit->isValid()) {
	DetId trkHitId = trkHit->geographicalId();
	if(trkHitId.det() != DetId::Tracker)
	  continue;
	GlobalPoint trkHitPosition = trkHit->globalPosition();
        if (verbosity_ > 2)
          switch (trkHitId.subdetId()) {
          case PixelSubdetector::PixelBarrel:
            cout << "Hit detId = " << (PXBDetId)trkHitId << ", position = " << trkHitPosition;
            break;
          case PixelSubdetector::PixelEndcap:
            cout << "Hit detId = " << (PXFDetId)trkHitId << ", position = " << trkHitPosition;
            break;
          case StripSubdetector::TIB:
            cout << "Hit detId = " << (TIBDetId)trkHitId << ", position = " << trkHitPosition;
            break;
          case StripSubdetector::TID:
            cout << "Hit detId = " << (TIDDetId)trkHitId << ", position = " << trkHitPosition;
            break;
          case StripSubdetector::TOB:
            cout << "Hit detId = " << (TOBDetId)trkHitId << ", position = " << trkHitPosition;
            break;
          case StripSubdetector::TEC:
            cout << "Hit detId = " << (TECDetId)trkHitId << ", position = " << trkHitPosition;
            break;
          default:
            break;
          }
	double trkHitPhi = atan2(trkHitPosition.y(), trkHitPosition.x());
	for (unsigned int iHit = 0; iHit < vHits.size(); iHit++) {
	  DetId hitId = (vHits[iHit])->geographicalId();
	  if ((hitId.rawId() & (~3)) == (trkHitId.rawId() & ~3)) {
	    GlobalPoint hitPosition = (vHits[iHit])->globalPosition();
	    double hitPhi = atan2(hitPosition.y(), hitPosition.x());
	    double hitR = sqrt(hitPosition.x()*hitPosition.x() + hitPosition.y()*hitPosition.y());
	    switch (trkHitId.subdetId()) {
	    case PixelSubdetector::PixelBarrel:
              if (fabs(hitPosition.x() - trkHitPosition.x()) > maxPixelDxy ||
                  fabs(hitPosition.y() - trkHitPosition.y()) > maxPixelDxy ||
                  fabs(hitPosition.z() - trkHitPosition.z()) > maxPixelDz)
                continue;
	      break;
	    case PixelSubdetector::PixelEndcap:
              if (fabs(hitPosition.x() - trkHitPosition.x()) > maxPixelDxy ||
                  fabs(hitPosition.y() - trkHitPosition.y()) > maxPixelDxy ||
                  fabs(hitPosition.z() - trkHitPosition.z()) > maxPixelDz)
                continue;
	      break;
	    case StripSubdetector::TIB:
	      if (((TIBDetId)hitId).isDoubleSide()) {
		if (!(((TIBDetId)trkHitId).isRPhi()) || hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi || fabs(hitPosition.z() - trkHitPosition.z()) > deltaZ)
		  continue;
	      } else if (hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi)
                continue;
	      break;
	    case StripSubdetector::TID:
	      if (((TIDDetId)hitId).isDoubleSide()) {
		if (!(((TIDDetId)trkHitId).isRPhi()) || hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi || fabs(hitPosition.z() - trkHitPosition.z()) > deltaZ)
		  continue;
	      } else if (hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi)
		continue;
	      break;
	    case StripSubdetector::TOB:
	      if (((TOBDetId)hitId).isDoubleSide()) {
		if (!(((TOBDetId)trkHitId).isRPhi()) || hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi || fabs(hitPosition.z() - trkHitPosition.z()) > deltaZ)
		  continue;
	      } else if (hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi)
		continue;
	      break;
	    case StripSubdetector::TEC:
	      if (((TECDetId)hitId).isDoubleSide()) {
		if (!(((TECDetId)trkHitId).isRPhi()) || hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi || fabs(hitPosition.z() - trkHitPosition.z()) > deltaZ)
		  continue;
	      } else if (hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi)
		continue;
	      break;
	    default:
	      continue;
	      break;
	    }
	    deltaZ = fabs(hitPosition.z() - trkHitPosition.z());
	    assHit = iHit;
            if (deltaZ < 4.)
              break;
	  }
	}
        if (assHit >= 0) {
          vTHits.push_back(assHit);
          if (verbosity_ > 2)
            cout << ": associated to hit " << assHit;
        }
        if (verbosity_ > 2)
          cout << endl;
      }
    }
    if (verbosity_ > 2)
      cout << "Number of hits in track: " << nHitsOrig << "; number of associated hits: " << vTHits.size() << endl;
    if (vTHits.size() >= 3)
      vHitsTrk.push_back(vTHits);
  }
  nAssTracks_ = vHitsTrk.size();

  // Compare content of 5-par bin over thresholds with track hits
  map<int, bool> assHits;  // map with hits in found tracks
  vector<bool> trackFound(nAssTracks_, false);  // vector with tracks found in HT histogram
  // Number of HT "seeds" surviving after using hits and "seeds" corresponding to tracks
  nTriedSeeds_ = 0;
  nAssSeeds_ = 0;
  for (map<unsigned long, map<unsigned int, vector<int> > >::iterator itBin = passVotes.begin(); itBin != passVotes.end(); itBin++) {
    map<unsigned int, bool> lyrTried;  // map with layers tried
    vector<map<unsigned int, bool> > lyrFound(nAssTracks_);  // vector with layers with association found
    unsigned int nVotes = 0;
    for (map<unsigned int, vector<int> >::iterator itLyr = (*itBin).second.begin(); itLyr != (*itBin).second.end(); itLyr++) {
      nVotes += (*itLyr).second.size();
      unsigned int lyrId = (*itLyr).first;
      for (vector<int>::iterator itHit = (*itLyr).second.begin(); itHit != (*itLyr).second.end(); itHit++) {
	if (assHits.find(*itHit) == assHits.end()) {
	  if (lyrTried.find(lyrId) == lyrTried.end())  // new layer tried
	    lyrTried[lyrId] = true;
	  for (unsigned int iTrk = 0; iTrk < nAssTracks_; iTrk++){
	    for (vector<int>::iterator itTrkHit = vHitsTrk[iTrk].begin(); itTrkHit != vHitsTrk[iTrk].end(); itTrkHit++) {
	      if ((*itHit) == (*itTrkHit)) {
		if (lyrFound[iTrk].find(lyrId) == lyrFound[iTrk].end())  // new layer associated
		  (lyrFound[iTrk])[lyrId] = true;
		break;
	      }
	    }
	  }
	}
      }
    }
    hNVotes_->Fill(nVotes);
    if (lyrTried.size() >= xyzVoteThr_) {
      nTriedSeeds_++;
      // Loop again on tracks to check whether some has at least three hits associated to HT votes in bin
      bool seedAss = false;
      for (unsigned int iTrk = 0; iTrk < nAssTracks_; iTrk++) {
	if ((lyrFound[iTrk]).size() >= 3) {
	  trackFound[iTrk] = true;
          if (!seedAss) {
            nAssSeeds_++;
            hNVotesMatched_->Fill(nVotes);
            seedAss = true;
          }
	  // Flag as associated all track hits
	  for (vector<int>::iterator itTrkHit = vHitsTrk[iTrk].begin(); itTrkHit != vHitsTrk[iTrk].end(); itTrkHit++)
	    assHits[*itTrkHit] = true;
	}
      }
    }
  }

  // Final loop on tracks to count how many were "found" by HT
  nFoundTracks_ = 0;
  for (unsigned int iTrk = 0; iTrk < nAssTracks_; iTrk++)
    if (trackFound[iTrk])
      nFoundTracks_++;
  // Print results
  cout << "Total number of seeds: " << nSeeds_ << endl;
  cout << "Seeds after clean-up: " << nGoodSeeds_ << endl;
  cout << "Seeds surviving after removing used hits: " << nTriedSeeds_ << endl;
  cout << "Seeds associated to tracks: " << nAssSeeds_ << endl;
  cout << "Total number of tracks in event: " << nTracks_ << endl;
  cout << "Tracks found by selected iterations: " << nTracksAlgo_ << endl;
  cout << "Tracks in parameter acceptance: " << nAccTracks_ << endl;
  cout << "Tracks with at least 3 associated hits: " << nAssTracks_ << endl;
  cout << "Tracks associated to seeds: " << nFoundTracks_ << endl;

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
  trackTree_->Branch("nSeeds", &nSeeds_);
  trackTree_->Branch("nGoodSeeds", &nGoodSeeds_);
  trackTree_->Branch("nTriedSeeds", &nTriedSeeds_);
  trackTree_->Branch("nAssSeeds", &nAssSeeds_);
  trackTree_->Branch("nTracks", &nTracks_);
  trackTree_->Branch("nTracksAlgo", &nTracksAlgo_);
  trackTree_->Branch("nAccTracks", &nAccTracks_);
  trackTree_->Branch("nAssTracks", &nAssTracks_);
  trackTree_->Branch("nFoundTracks", &nFoundTracks_);

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

  // Book histograms
  hXYHoughVotes_.reset(new TH3S("hXYHoughVotes", "helix Hough transform votes", nBins_[0], minPar_[0], maxPar_[0], nBins_[1], minPar_[1], maxPar_[1], nBins_[2], minPar_[2], maxPar_[2]));
  hXYHoughVotes_->SetDirectory(0);
  hNVotes_.reset(new TH1I("hNVotes", "number of votes per bin", 200, 0., 200.));
  hNVotes_->SetDirectory(0);
  hNVotesMatched_.reset(new TH1I("hNVotesMatched", "number of votes per track-matched bin", 200, 0., 200.));
  hNVotesMatched_->SetDirectory(0);
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HoughCheck2Steps::endJob() 
{
  TFile outRootFile(outRootFile_.c_str(), "RECREATE");
  hXYHoughVotes_->Write();
  hNVotes_->Write();
  hNVotesMatched_->Write();
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
  if (!layerListName_.empty()) {
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
