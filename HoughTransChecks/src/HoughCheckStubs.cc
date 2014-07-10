// -*- C++ -*-
//
// Package:    HoughTransChecks
// Class:      HoughCheckStubs
// 
/**\class HoughCheckStubs HoughCheckStubs.cc ERobutti/HoughTransChecks/src/HoughCheckStubs.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Enrico Robutti
//         Created:  Tue, October 15, 2013
// $Id: HoughCheckStubs.cc
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

#include "ERobutti/HoughTransChecks/interface/TH5.h"
#include "TH1I.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"


//
// class declaration
//

using namespace std;
using namespace edm;

class HoughCheckStubs : public edm::EDAnalyzer {
public:
  explicit HoughCheckStubs(const edm::ParameterSet&);
  ~HoughCheckStubs();
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
  double binOverlap_;
  double binWidth_[5];
  double parOverlap_[5];
  double minKappaScan_, maxKappaScan_;
  const unsigned int voteThr_;
  const bool cleanupSeeds_;
  const string outRootFile_;
  const unsigned int verbosity_;

  // Histogram and tree 
  auto_ptr<TH5C> hHoughVotes_;
  auto_ptr<TH1I> hNVotes_;
  auto_ptr<TH1I> hNVotesMatched_;
  auto_ptr<TTree> trackTree_;
  vector<int> vSubdet_;
  vector<int> vLayer_;
  vector<double> vXHit_;
  vector<double> vYHit_;
  vector<double> vZHit_;
  vector<pair<int, int> > vPair_;
  vector<double> vDoca_;
  vector<double> vKappa_;
  vector<double> vPhi_;
  vector<double> vZ0_;
  vector<double> vTheta_;
  vector<double> vAlgo_;
  vector<vector<int> > vTHit_;
  
  unsigned int nSeeds_;
  unsigned int nGoodSeeds_;
  unsigned int nTriedSeeds_;
  unsigned int nAssSeeds_;
  unsigned int nTracks_;
  unsigned int nTracksAlgo_;
  unsigned int nAccTracks_;
  unsigned int nAssTracks_;
  unsigned int nFoundTracks_;

  // Layer mapping and maximum hit distance in pair
  static const int lyrMapOffset_[7];
  map<pair<int, int>, float> maxPairDist_;
};

//
// constants, enums and typedefs
//
const int HoughCheckStubs::lyrMapOffset_[7] = {0, 0, 3, 5, 9, 12, 18};

//
// static data member definitions
//

//
// constructors and destructor
//
HoughCheckStubs::HoughCheckStubs(const edm::ParameterSet& iConfig)
:
  trackTags_(iConfig.getParameter<edm::InputTag>("tracks")),
  builderName_(iConfig.getParameter<std::string>("TTRHBuilder")),
  algoSel_(iConfig.getParameter<vector<unsigned int> >("algoSel")),
  layerListName_(iConfig.getParameter<std::string> ("seedingLayers")),
  minPar_(iConfig.getUntrackedParameter<vector<double> >("minPar", {-20, -0.1, -M_PI, -50., -2.5})),
  maxPar_(iConfig.getUntrackedParameter<vector<double> >("maxPar", {20, 0.1, M_PI, 50., 2.5})),
  nBins_(iConfig.getUntrackedParameter<vector<int> >("nBins", {100, 100, 100, 100, 100})),
  binOverlap_(iConfig.getParameter<double>("binOverlap")),
  voteThr_(iConfig.getParameter<unsigned int>("voteThr")),
  cleanupSeeds_(iConfig.getUntrackedParameter<bool>("cleanupSeeds", false)),
  outRootFile_(iConfig.getUntrackedParameter<string>("outRootFile", "houghCheck_stubs.root")),
  verbosity_(iConfig.getUntrackedParameter<unsigned int>("verbosity", 0))
{
}


HoughCheckStubs::~HoughCheckStubs()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HoughCheckStubs::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using reco::TrackCollection;

  // Reset vectors
  vSubdet_.clear();
  vLayer_.clear();
  vXHit_.clear();
  vYHit_.clear();
  vZHit_.clear();
  vPair_.clear();

  if (!(hHoughVotes_.get()))
    return;
  // Loop on layers and hits within, store in auxiliary vectors and create allowed pairs
  map<unsigned int, vector<pair<int, int> > > votes;
  for (ctfseeding::SeedingLayerSets::const_iterator itLS = layerSets_.begin(); itLS != layerSets_.end(); itLS++ ) {
    for (ctfseeding::SeedingLayers::const_iterator itLyr = itLS->begin(); itLyr!= itLS->end(); itLyr++ ) {
      ctfseeding::SeedingLayer::Hits hits = itLyr->hits(iEvent, iSetup);
      for (ctfseeding::SeedingLayer::Hits::const_iterator itHit = hits.begin(); itHit != hits.end(); itHit++) {
        if (!((*itHit)->isValid()))
          continue;
        DetId hitId = (*itHit)->geographicalId();
        if (hitId.det() != DetId::Tracker)
          continue;
        int hitSubDet = hitId.subdetId();
        int hitLayer = 0;
        if (hitSubDet == PixelSubdetector::PixelBarrel) {
          if (verbosity_ > 2)
            cout << "hitId = " << PXBDetId(hitId) << endl;
          hitLayer = PXBDetId(hitId).layer();
        }
        else if (hitSubDet == PixelSubdetector::PixelEndcap) {
          if (verbosity_ > 2)
            cout << "hitId = " << PXFDetId(hitId) << endl;
          hitLayer = PXFDetId(hitId).disk();
        }
        else if (hitSubDet == StripSubdetector::TIB) {
          if (verbosity_ > 2)
            cout << "hitId = " << TIBDetId(hitId) << endl;
          hitLayer = TIBDetId(hitId).layer();
        }
        else if (hitSubDet == StripSubdetector::TID) {
          if (verbosity_ > 2)
            cout << "hitId = " << TIDDetId(hitId) << endl;
          hitLayer = TIDDetId(hitId).wheel();
        }
        else if (hitSubDet == StripSubdetector::TOB) {
          if (verbosity_ > 2)
            cout << "hitId = " << TOBDetId(hitId) << endl; 
          hitLayer = TOBDetId(hitId).layer();
        }
        else if (hitSubDet == StripSubdetector::TEC) {
          if (verbosity_ > 2)
            cout << "hitId = " << TECDetId(hitId) << endl;
          hitLayer = TECDetId(hitId).wheel();
        }
        else 
          continue;
        vSubdet_.push_back(hitSubDet);
        vLayer_.push_back(hitLayer);
        GlobalPoint hitPosition = (*itHit)->globalPosition();
        vXHit_.push_back(hitPosition.x());
        vYHit_.push_back(hitPosition.y());
        vZHit_.push_back(hitPosition.z());
        // int lyrId = lyrMapOffset_[hitSubDet] + hitLayer - 1;
        // vHits.push_back(make_pair(lyrId, *itHit));
      }
    }
    // Loop on hit pairs and transform
    if (verbosity_ > 2)
      cout << "Total number of good hits: " << vSubdet_.size() << endl;
    for (unsigned int iHit1 = 0; iHit1 < vSubdet_.size(); iHit1++) {
      if (verbosity_ > 2)
        cout << "hit #" << iHit1 << ": ";
      int lyr1 = lyrMapOffset_[vSubdet_[iHit1]] + vLayer_[iHit1] - 1;
      // TransientTrackingRecHit::ConstRecHitPointer hit1 = vHits[iHit1].second;
      // GlobalPoint hitPosition1 = hit1->globalPosition();
      double x1 = vXHit_[iHit1];
      double y1 = vYHit_[iHit1];
      double z1 = vZHit_[iHit1];
      if (verbosity_ > 2)
        cout << "position = (" << x1 << ", " << y1 << ", " << z1 << ")" << endl;
      for (unsigned int iHit2 = iHit1 + 1; iHit2 < vSubdet_.size(); iHit2++) {
        int lyr2 = lyrMapOffset_[vSubdet_[iHit2]] + vLayer_[iHit2] - 1;
        // TransientTrackingRecHit::ConstRecHitPointer hit2 = vHits[iHit2].second;
        // GlobalPoint hitPosition2 = hit2->globalPosition();
        double x2 = vXHit_[iHit2];
        double y2 = vYHit_[iHit2];
        double z2 = vZHit_[iHit2];
        double pairDist = sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1) + (z2 - z1)*(z2 - z1));
        if (maxPairDist_.count(make_pair(lyr1, lyr2)) == 0 || pairDist > maxPairDist_.at(make_pair(lyr1, lyr2)))
          continue;
        vPair_.push_back(make_pair(iHit1, iHit2));
        if (verbosity_ > 2)
          cout << "Pair (" << iHit1 << ", " << iHit2 << ") accepted (distance = " << pairDist << ")." << endl;
        double xM = 0.5*(x1 + x2);
        double yM = 0.5*(y1 + y2);
        double a = (x1 - x2)/(y2 - y1);
        double c2 = (x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1);
	// Loop on kappa and fill histogram
        for (double kappa = minKappaScan_; kappa < maxKappaScan_; kappa += binWidth_[1]) {
          double h2 = 1./(kappa*kappa) - c2/4.;
          double dxM = sqrt(h2/(1 + a*a));
          for (int iSign = -1; iSign <= 1; iSign += 2) {
            int signh = ((y1 - y2) > 0) ? iSign : -iSign;
            double xC = xM + signh*dxM;
            double yC = yM + a*signh*dxM;
            //            double lnr2 = log(xC*xC + yC*yC);
            double doca = iSign*(1./kappa - sqrt(xC*xC + yC*yC));
            //            double phiC = atan2(yC, xC);
            double phi = atan2(-iSign*xC, iSign*yC);
            double Ds = 2./kappa*asin(kappa*sqrt(c2)/2.);
            double theta = atan2(Ds, z2 - z1);
            double eta = -log(tan(theta/2.));
            double s1 = iSign*(atan2(y1 - yC, x1 - xC) - atan2(-yC, -xC))/kappa;
            if (s1 < 0)
              s1 += 2*M_PI/kappa;
            double z0 = z1 - (z2 - z1)*s1/Ds;
            hHoughVotes_->Fill(doca, iSign*kappa, phi, z0, eta);
          }
        }
      }
    }
  }

  //             // Insert point in vote map, taking overlaps into account
  //             int firstPhiBin = max(0, (int)((phi - minPar_[2] - phiOverlap)/phiBinWidth));
  //             int lastPhiBin = min(nBins_[2] - 1, (int)((phi - minPar_[2] + phiOverlap)/phiBinWidth));
  //             for (int iPhi = firstPhiBin; iPhi <= lastPhiBin; iPhi++) {
  //       	unsigned int keyBin = iDoca | (iSqrtK << 8) | (iPhi << 16);
  //       	map<unsigned int, vector<pair<int, double> > >& lyrMap = xyVotes[keyBin];
  //       	unsigned int keyLyr = (hitSubDet << 4) | hitLayer;
  //       	vector<pair<int, double> >& vLyrHits = lyrMap[keyLyr];
  //       	vLyrHits.push_back(make_pair(nHit, phi));
  //             }
  //           }
  //         }
  //       } else
  //         if (verbosity_ > 2)
  //           cout << " - invalid hit" << endl;
  //     }

  // // Longitudinal step
  // map<unsigned long, map<unsigned int, vector<int> > > xyzVotes;
  // // Fill 5-par map for 3-par bins over threshold
  // for (map<unsigned int, map<unsigned int, vector<pair<int, double> > > >::iterator itBin = xyVotes.begin(); itBin != xyVotes.end(); itBin++) {
  //   unsigned int keyBin = (*itBin).first;
  //   int docaBin = (keyBin & 255) + 1;
  //   int sqrtKBin = ((keyBin >> 8) & 255) + 1;
  //   int phiBin = ((keyBin >> 16) & 255) + 1;
  //   unsigned int lyrVotes = (*itBin).second.size();
  //   // Fill x-y histogram with number of voting layers
  //   hXYHoughVotes_->SetBinContent(docaBin, sqrtKBin, phiBin, lyrVotes);
  //   if (lyrVotes >= xyVoteThr_) {
  //     double doca = minPar_[0] + (docaBin + 0.5)*docaBinWidth;
  //     double sqrtK = minPar_[1] + (sqrtKBin + 0.5)*sqrtKBinWidth;
  //     double kappa = sqrtK*fabs(sqrtK);
  //     for (map<unsigned int, vector<pair<int, double> > >::iterator itLyr = (*itBin).second.begin(); itLyr != (*itBin).second.end(); itLyr++) {
  //       unsigned int keyLyr = (*itLyr).first;
  //       for (vector<pair<int, double> >::iterator itHit = (*itLyr).second.begin(); itHit != (*itLyr).second.end(); itHit++) {
  //         TransientTrackingRecHit::ConstRecHitPointer hit = vHits[(*itHit).first];
  //         GlobalPoint hitPosition = hit->globalPosition();
  //         double x = 10.*hitPosition.x();
  //         double y = 10.*hitPosition.y();
  //         double z = 10.*hitPosition.z();
  //         double phi = (*itHit).second;
  //         double xc = (doca - 1./kappa)*sin(phi);
  //         double yc = -(doca - 1./kappa)*cos(phi);
  //         for (int iZ0 = 0; iZ0 < nBins_[3]; iZ0++) {
  //           double z0 = minPar_[3] + (iZ0 + 0.5)*z0BinWidth;
  //           double st = 1./kappa*(atan2(kappa*(y - yc), kappa*(x - xc)) - phi + M_PI/2.);
  //           if (st < 0)  // rotation angle has crossed the +/-PI border
  //             st += 2.*M_PI/fabs(kappa);
  //           else if (st > 2.*M_PI/fabs(kappa))
  //             st -= 2.*M_PI/fabs(kappa);
  //           double lambda = atan((z - z0)/st);
  //           double eta = -log(tan((M_PI/2. - lambda)/2.));
  //           // Insert point in 5-parameter vote map, taking overlaps into account
  //           int firstEtaBin = max(0, (int)((eta - minPar_[4] - etaOverlap)/etaBinWidth));
  //           int lastEtaBin = min(nBins_[4] - 1, (int)((eta - minPar_[4] + etaOverlap)/etaBinWidth));
  //           for (int iEta = firstEtaBin; iEta <= lastEtaBin; iEta++) {
  //             unsigned long key5Par = keyBin | (iZ0 << 24) | ((unsigned long)iEta << 32);
  //             map<unsigned int, vector<int> >& lyrMap = xyzVotes[key5Par];
  //             vector<int>& vLyrHits = lyrMap[keyLyr];
  //             vLyrHits.push_back((*itHit).first);
  //           }
  //         }
  //       }
  //     }      
  //   }
  // }
  // nSeeds_ = xyzVotes.size();

  // // Clean-up seeds: remove below threshold and duplicates if requested
  map<unsigned long, map<unsigned int, vector<int> > > passVotes;
  // for (map<unsigned long, map<unsigned int, vector<int> > >::iterator itBin = xyzVotes.begin(); itBin != xyzVotes.end(); itBin++) {
  //   bool nextBin = false;
    // unsigned int lyrVotes = (*itBin).second.size();
    // if (lyrVotes >= xyzVoteThr_) {
    //   if (cleanupSeeds_) {
    //     // Get bin
    //     unsigned long htBin = (*itBin).first;
    //     int docaBin = htBin & 255;
    //     int sqrtKBin = (htBin >> 8) & 255;
    //     int phiBin = (htBin >> 16) & 255;
    //     int z0Bin = (htBin >> 24) & 255;
    //     int etaBin = (htBin >> 32) & 255;
    //     int deltaBin = 1;  // width of bin range where to look for duplicates
    //     for (int iDoca = max(0, docaBin - deltaBin); !nextBin && iDoca < min(nBins_[0], docaBin + deltaBin + 1); iDoca++)
    //       for (int iSqrtK = max(0, sqrtKBin - deltaBin); !nextBin && iSqrtK < min(nBins_[1], sqrtKBin + deltaBin + 1); iSqrtK++)
    //         for (int iPhi = max(0, phiBin - deltaBin); !nextBin && iPhi < min(nBins_[2], phiBin + deltaBin + 1); iPhi++)
    //           for (int iZ0 = max(0, z0Bin - deltaBin); !nextBin && iZ0 < min(nBins_[3], z0Bin + deltaBin + 1); iZ0++)
    //     	for (int iEta = max(0, etaBin - deltaBin); !nextBin && iEta < min(nBins_[4], etaBin + deltaBin + 1); iEta++) {
    //     	  unsigned long key = iDoca | (iSqrtK << 8) | (iPhi << 16) | ((unsigned long)iZ0 << 24) | ((unsigned long)iEta << 32);
    //     	  if (key == htBin || xyzVotes.find(key) == xyzVotes.end())
    //     	    continue;
    //     	  bool hitFound = true;
    //     	  for (map<unsigned int, vector<int> >::iterator itLyr1 = (*itBin).second.begin(); hitFound && itLyr1 != (*itBin).second.end(); itLyr1++) {
    //     	    for (vector<int>::iterator itHit1 = (*itLyr1).second.begin(); hitFound && itHit1 != (*itLyr1).second.end(); itHit1++) {
    //     	      hitFound = false;
    //     	      for (map<unsigned int, vector<int> >::iterator itLyr2 = xyzVotes[key].begin(); !hitFound && itLyr2 != xyzVotes[key].end(); itLyr2++) {
    //     		if ((*itLyr2).first != (*itLyr1).first)
    //     		  continue;
    //     		for (vector<int>::iterator itHit2 = (*itLyr2).second.begin(); itHit2 != (*itLyr2).second.end(); itHit2++) {
    //     		  if ((*itHit1) == (*itHit2)) {
    //     		    hitFound = true;
    //     		    break;
    //     		  }
    //     		}
    //     	      }
    //     	    }
    //     	  }
    //     	  if (hitFound) {  // all hits found in a "nearby" bin
    //     	    (*itBin).second.clear();  // remove all hits from the bin
  //       	    nextBin = true;
  //       	  }
  //       	}
  //     }
  //     if (!nextBin)  // no nearby bin containing all hits (or no clean-up)
  //       passVotes[(*itBin).first] = (*itBin).second;
  //   }
  // }
  // nGoodSeeds_ = passVotes.size();

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
    vDoca_.push_back(perigeePars.transverseImpactParameter());
    vKappa_.push_back(perigeePars.transverseCurvature());
    vPhi_.push_back(perigeePars.phi());
    vZ0_.push_back(perigeePars.longitudinalImpactParameter());
    vTheta_.push_back(perigeePars.theta());
    vAlgo_.push_back(itTrack->algo());
    // Check if track is in histogram acceptance
    // if (vDoca_.back() < minPar_[0] || vDoca_.back() > maxPar_[0] ||
    //     vKappa_.back()/sqrt(fabs(vKappa_.back())) < minPar_[1] || vKappa_.back()/sqrt(fabs(vKappa_.back())) > maxPar_[1] ||
    //     vPhi_.back() < minPar_[2] || vPhi_.back() > maxPar_[2] ||
    //     vZ0_.back() < minPar_[3] || vZ0_.back() > maxPar_[3] ||
    //     -log(tan(0.5*vTheta_.back())) < minPar_[4] || -log(tan(0.5*vTheta_.back())) > maxPar_[4])
    //   continue;
    nAccTracks_++;      
    vector<int> vTHits;
    // vector<double> vTSubdet;
    // vector<double> vTLayer;
    // vector<double> vTXHit;
    // vector<double> vTYHit;
    // vector<double> vTZHit;
    unsigned int nHitsOrig = 0;
    for (trackingRecHit_iterator i = itTrack->recHitsBegin(); i != itTrack->recHitsEnd(); i++) {
      nHitsOrig++;
      int assHit = -1;
      // double maxPixelDxy = 0.01;
      // double maxPixelDz = 0.01;
      // double maxStripDrphi = 0.02;
      // double deltaZ = 1000;
      TransientTrackingRecHit::RecHitPointer trkHit = builder_->build(&**i );
      if (trkHit->isValid()) {
	DetId trkHitId = trkHit->geographicalId();
	if(trkHitId.det() != DetId::Tracker)
	  continue;
	GlobalPoint trkHitPosition = trkHit->globalPosition();
        //	vTSubdet.push_back(trkHitId.subdetId());
	switch (trkHitId.subdetId()) {
	case PixelSubdetector::PixelBarrel:
          //	  vTLayer.push_back(PXBDetId(trkHitId).layer());
	  if (verbosity_ > 2)
	    cout << "Hit detId = " << PXBDetId(trkHitId) << ", position = " << trkHitPosition;
	  break;
        case PixelSubdetector::PixelEndcap:
          //	  vTLayer.push_back(PXFDetId(trkHitId).disk());
	  if (verbosity_ > 2)
            cout << "Hit detId = " << PXFDetId(trkHitId) << ", position = " << trkHitPosition;
	  break;
        case StripSubdetector::TIB:
          //	  vTLayer.push_back(TIBDetId(trkHitId).layer());
	  if (verbosity_ > 2)
            cout << "Hit detId = " << TIBDetId(trkHitId) << ", position = " << trkHitPosition;
	  break;
        case StripSubdetector::TID:
          //	  vTLayer.push_back(TIDDetId(trkHitId).wheel());
	  if (verbosity_ > 2)
            cout << "Hit detId = " << TIDDetId(trkHitId) << ", position = " << trkHitPosition;
	  break;
        case StripSubdetector::TOB:
          //	  vTLayer.push_back(TOBDetId(trkHitId).layer());
	  if (verbosity_ > 2)
            cout << "Hit detId = " << TOBDetId(trkHitId) << ", position = " << trkHitPosition;
	  break;
        case StripSubdetector::TEC:
          //	  vTLayer.push_back(TECDetId(trkHitId).wheel());
	  if (verbosity_ > 2)
	    cout << "Hit detId = " << TECDetId(trkHitId) << ", position = " << trkHitPosition;
	  break;
        default:
	  break;
	}
        //	vTXHit.push_back(trkHitPosition.x());
        //	vTYHit.push_back(trkHitPosition.y());
        //	vTZHit.push_back(trkHitPosition.z());
        //	double trkHitPhi = atan2(trkHitPosition.y(), trkHitPosition.x());
	//	for (unsigned int iHit = 0; iHit < vHits.size(); iHit++) {
	//   DetId hitId = (vHits[iHit])->geographicalId();
	//   if ((hitId.rawId() & (~3)) == (trkHitId.rawId() & ~3)) {
	//     GlobalPoint hitPosition = (vHits[iHit])->globalPosition();
	//     double hitPhi = atan2(hitPosition.y(), hitPosition.x());
	//     double hitR = sqrt(hitPosition.x()*hitPosition.x() + hitPosition.y()*hitPosition.y());
	//     switch (trkHitId.subdetId()) {
	//     case PixelSubdetector::PixelBarrel:
        //       if (fabs(hitPosition.x() - trkHitPosition.x()) > maxPixelDxy ||
        //           fabs(hitPosition.y() - trkHitPosition.y()) > maxPixelDxy ||
        //           fabs(hitPosition.z() - trkHitPosition.z()) > maxPixelDz)
        //         continue;
	//       break;
	//     case PixelSubdetector::PixelEndcap:
        //       if (fabs(hitPosition.x() - trkHitPosition.x()) > maxPixelDxy ||
        //           fabs(hitPosition.y() - trkHitPosition.y()) > maxPixelDxy ||
        //           fabs(hitPosition.z() - trkHitPosition.z()) > maxPixelDz)
        //         continue;
	//       break;
	//     case StripSubdetector::TIB:
	//       if (((TIBDetId)hitId).isDoubleSide()) {
	// 	if (!(((TIBDetId)trkHitId).isRPhi()) ||  hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi || fabs(hitPosition.z() - trkHitPosition.z()) > deltaZ)
	// 	  continue;
	//       } else if (hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi)
        //         continue;
	//       break;
	//     case StripSubdetector::TID:
	//       if (((TIDDetId)hitId).isDoubleSide()) {
	// 	if (!(((TIDDetId)trkHitId).isRPhi()) || hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi || fabs(hitPosition.z() - trkHitPosition.z()) > deltaZ)
	// 	  continue;
	//       } else if (hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi)
	// 	continue;
	//       break;
	//     case StripSubdetector::TOB:
	//       if (((TOBDetId)hitId).isDoubleSide()) {
	// 	if (!(((TOBDetId)trkHitId).isRPhi()) || hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi || fabs(hitPosition.z() - trkHitPosition.z()) > deltaZ)
	// 	  continue;
	//       } else if (hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi)
	// 	continue;
	//       break;
	//     case StripSubdetector::TEC:
	//       if (((TECDetId)hitId).isDoubleSide()) {
	// 	if (!(((TECDetId)trkHitId).isRPhi()) || hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi || fabs(hitPosition.z() - trkHitPosition.z()) > deltaZ)
	// 	  continue;
	//       } else if (hitR*fabs(hitPhi - trkHitPhi) > maxStripDrphi)
	// 	continue;
	//       break;
	//     default:
	//       continue;
	//       break;
	//     }
	//     deltaZ = fabs(hitPosition.z() - trkHitPosition.z());
	//     assHit = iHit;
        //     if (deltaZ < 4.)
        //       break;
	//   }
	// }
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
    // vSubdet_.push_back(vTSubdet);
    // vLayer_.push_back(vTLayer);
    // vXHit_.push_back(vTXHit);
    // vYHit_.push_back(vTYHit);
    // vZHit_.push_back(vTZHit);
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
	  for (unsigned int iTrk = -1; iTrk < nAssTracks_; iTrk++){
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
    if (lyrTried.size() >= voteThr_) {
      nTriedSeeds_++;
      // Loop again on tracks to check whether some has at least three hits associated to HT votes in bin
      bool seedAss = false;
      for (unsigned int iTrk = 0; iTrk < nAssTracks_; iTrk++) {
	if ((lyrFound[iTrk]).size() >= 3) {
	  trackFound[iTrk] = true;
          if (!seedAss) {
            nAssSeeds_++;
            hNVotesMatched_->Fill(nVotes);
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
HoughCheckStubs::beginJob()
{

  // Create tree, with branches
  trackTree_.reset(new TTree("trackTree", "Fitted track parameters"));
  trackTree_->Branch("subdet", &vSubdet_);
  trackTree_->Branch("layer", &vLayer_);
  trackTree_->Branch("xHit", &vXHit_);
  trackTree_->Branch("yHit", &vYHit_);
  trackTree_->Branch("zHit", &vZHit_);
  trackTree_->Branch("pair", &vPair_);
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
    if (iPar == 2 && (minPar_[iPar] < -M_PI || maxPar_[iPar] > M_PI)) {
      cout << "Invalid phi range. Valid range is -pi <= phi <= pi. No histogram booked." << endl;
      return;
    }
    binWidth_[iPar] = (maxPar_[iPar] - minPar_[iPar])/nBins_[iPar];
    if (iPar == 1) {
      if (minPar_[iPar] < 0 && maxPar_[iPar] > 0) {  // shift the range so that 0 is a bin edge
        float shift = binWidth_[iPar]*int((maxPar_[iPar] + 0.5*binWidth_[iPar])/binWidth_[iPar]) - maxPar_[iPar];
        minPar_[iPar] += shift;
        maxPar_[iPar] += shift;
        minKappaScan_ = 0.5*binWidth_[iPar];
        maxKappaScan_ = max(fabs(minPar_[iPar]), fabs(maxPar_[iPar]));
      } else {
        minKappaScan_ = min(fabs(minPar_[iPar]), fabs(maxPar_[iPar])) + 0.5*binWidth_[iPar];
        maxKappaScan_ = max(fabs(minPar_[iPar]), fabs(maxPar_[iPar]));
      }
    } else {
      parOverlap_[iPar] = binOverlap_*binWidth_[iPar];
    }
  }
  if (verbosity_ > 2)
    cout << "minKappaScan_, maxKappaScan_, binWidth_[1] = " << minKappaScan_ << ", " << maxKappaScan_ << ", " << binWidth_[1] << endl;

  // Book histograms
  hHoughVotes_.reset(new TH5C("hHoughVotes", "helix Hough transform votes", nBins_[0], minPar_[0], maxPar_[0], nBins_[1], minPar_[1], maxPar_[1], nBins_[2], minPar_[2], maxPar_[2], nBins_[3], minPar_[3], maxPar_[3], nBins_[4], minPar_[4], maxPar_[4]));
  hHoughVotes_->SetDirectory(0);
  hNVotes_.reset(new TH1I("hNVotes", "number of votes per bin", 200, 0., 200.));
  hNVotes_->SetDirectory(0);
  hNVotesMatched_.reset(new TH1I("hNVotes", "number of votes per track-matched bin", 200, 0., 200.));
  hNVotesMatched_->SetDirectory(0);

  // Set content to -128 for all histogram bins, to use full 8-bit range
  int nBinsTot = 1;
  for (int iPar = 0; iPar < 5; iPar++)
    nBinsTot *= nBins_[iPar] + 2;
  for (int iBin = 0; iBin < nBinsTot; iBin++)
    hHoughVotes_->SetBinContent(iBin, -128);

  // Initialize map with maximum distance for consecutive hits
  maxPairDist_[(make_pair(0, 1))] = 14.;
  maxPairDist_[(make_pair(0, 2))] = 10.;
  maxPairDist_[(make_pair(0, 3))] = 24.;
  maxPairDist_[(make_pair(1, 2))] = 10.;
  maxPairDist_[(make_pair(1, 3))] = 20.;
  maxPairDist_[(make_pair(1, 5))] = 36.;
  maxPairDist_[(make_pair(2, 3))] = 12.;
  maxPairDist_[(make_pair(2, 5))] = 44.;
  maxPairDist_[(make_pair(2, 6))] = 40.;
  maxPairDist_[(make_pair(2, 9))] = 64.;
  maxPairDist_[(make_pair(3, 4))] = 16.;
  maxPairDist_[(make_pair(3, 5))] = 30.;
  maxPairDist_[(make_pair(3, 9))] = 52.;
  maxPairDist_[(make_pair(3, 10))] = 64.;
  maxPairDist_[(make_pair(3, 11))] = 76.;
  maxPairDist_[(make_pair(3, 18))] = 94.;
  maxPairDist_[(make_pair(3, 19))] = 108.;
  maxPairDist_[(make_pair(4, 9))] = 36.;
  maxPairDist_[(make_pair(4, 10))] = 50.;
  maxPairDist_[(make_pair(4, 11))] = 64.;
  maxPairDist_[(make_pair(4, 18))] = 90.;
  maxPairDist_[(make_pair(4, 19))] = 104.;
  maxPairDist_[(make_pair(4, 20))] = 110.;
  maxPairDist_[(make_pair(5, 6))] = 28.;
  maxPairDist_[(make_pair(5, 7))] = 28.;
  maxPairDist_[(make_pair(5, 8))] = 30.;
  maxPairDist_[(make_pair(5, 9))] = 42.;
  maxPairDist_[(make_pair(5, 10))] = 34.;
  maxPairDist_[(make_pair(6, 7))] = 26.;
  maxPairDist_[(make_pair(6, 8))] = 26.;
  maxPairDist_[(make_pair(6, 9))] = 26.;
  maxPairDist_[(make_pair(6, 10))] = 32.;
  maxPairDist_[(make_pair(7, 8))] = 22.;
  maxPairDist_[(make_pair(7, 9))] = 22.;
  maxPairDist_[(make_pair(7, 12))] = 48.;
  maxPairDist_[(make_pair(8, 12))] = 30.;
  maxPairDist_[(make_pair(8, 13))] = 20.;
  maxPairDist_[(make_pair(9, 10))] = 20.;
  maxPairDist_[(make_pair(9, 11))] = 24.;
  maxPairDist_[(make_pair(9, 12))] = 30.;
  maxPairDist_[(make_pair(9, 18))] = 72.;
  maxPairDist_[(make_pair(10, 11))] = 22.;
  maxPairDist_[(make_pair(10, 18))] = 56.;
  maxPairDist_[(make_pair(11, 18))] = 42.;
  maxPairDist_[(make_pair(11, 19))] = 40.;
  maxPairDist_[(make_pair(12, 13))] = 24.;
  maxPairDist_[(make_pair(12, 18))] = 46.;
  maxPairDist_[(make_pair(13, 14))] = 26.;
  maxPairDist_[(make_pair(13, 15))] = 22.;
  maxPairDist_[(make_pair(13, 18))] = 56.;
  maxPairDist_[(make_pair(14, 15))] = 24.;
  maxPairDist_[(make_pair(14, 16))] = 28.;
  maxPairDist_[(make_pair(14, 17))] = 38.;
  maxPairDist_[(make_pair(14, 18))] = 70.;
  maxPairDist_[(make_pair(14, 19))] = 56.;
  maxPairDist_[(make_pair(15, 16))] = 28.;
  maxPairDist_[(make_pair(15, 17))] = 32.;
  maxPairDist_[(make_pair(15, 18))] = 38.;
  maxPairDist_[(make_pair(16, 17))] = 28.;
  maxPairDist_[(make_pair(18, 19))] = 26.;
  maxPairDist_[(make_pair(18, 20))] = 22.;
  maxPairDist_[(make_pair(19, 20))] = 24.;
  maxPairDist_[(make_pair(19, 21))] = 26.;
  maxPairDist_[(make_pair(20, 21))] = 26.;
  maxPairDist_[(make_pair(20, 22))] = 34.;
  maxPairDist_[(make_pair(20, 23))] = 50.;
  maxPairDist_[(make_pair(21, 22))] = 28.;
  maxPairDist_[(make_pair(22, 23))] = 30.;
  maxPairDist_[(make_pair(22, 24))] = 36.;
  maxPairDist_[(make_pair(23, 24))] = 32.;
  maxPairDist_[(make_pair(23, 25))] = 44.;
  maxPairDist_[(make_pair(24, 25))] = 34.;
  maxPairDist_[(make_pair(24, 26))] = 42.;
  maxPairDist_[(make_pair(25, 26))] = 36.;
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HoughCheckStubs::endJob() 
{
  TFile outRootFile(outRootFile_.c_str(), "RECREATE");
  hHoughVotes_->Write();
  hNVotes_->Write();
  hNVotesMatched_->Write();
  trackTree_->Write();
}

// ------------ method called when starting to processes a run  ------------
void 
HoughCheckStubs::beginRun(edm::Run const& run, edm::EventSetup const& setup)
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
HoughCheckStubs::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
HoughCheckStubs::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
HoughCheckStubs::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HoughCheckStubs::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("tracks", edm::InputTag("generalTracks"));
  desc.add<std::string>("TTRHBuilder");
  desc.add<vector<unsigned int> >("algoSel", {});
  desc.add<std::string>("seedingLayers", "HoughTransformSeedLayersAllHitsOneSet");
  desc.addUntracked<vector<double> >("minPar", {-20., -0.1, -M_PI, -50., -2.5});
  desc.addUntracked<vector<double> >("maxPar", {20., 0.1, M_PI, 50., 2.5});
  desc.addUntracked<vector<int> >("nBins", {100, 100, 100, 100, 100});
  desc.add<double>("binOverlap", 0.);
  desc.add<unsigned int>("voteThr", 0);
  desc.addUntracked<bool>("cleanupSeeds", false);
  desc.addUntracked<string>("outRootFile", "houghCheck_stubs.root");
  desc.addUntracked<unsigned int>("verbosity", 0);
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HoughCheckStubs);
