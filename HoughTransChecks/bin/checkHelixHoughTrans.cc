#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <math.h>

#include <TFile.h>
#include <TString.h>
#include <TH3S.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraph2D.h>
#include "../interface/TH5.h"

#define M_PI 3.14159265358979323846

using namespace std;

void getPhiTheta(double hitPos[], double htPar[])
{
  if (fabs(htPar[1]) < 1.e-6)
    htPar[1] = (htPar[1] >= 0) ? 1.e-6 : -1.e-6;  // avoid kappa = 0 (maximum curvature radius = 1 km)
  double r = sqrt(hitPos[0]*hitPos[0] + hitPos[1]*hitPos[1]);
  double akappa = (htPar[1]*htPar[0]*htPar[0] + htPar[1]*r*r - 2*htPar[0])/(2*r);
  double hkappa = sqrt((1 - htPar[0]*htPar[1])*(1 - htPar[0]*htPar[1]) - akappa*akappa);
  double kappaxD = -akappa*hitPos[0]/r + hkappa*hitPos[1]/r;
  double kappayD = -akappa*hitPos[1]/r - hkappa*hitPos[0]/r;
  htPar[2] = atan2(kappayD, kappaxD) + M_PI/2.;
  htPar[2] += -2.*M_PI*(int((htPar[2] + 3.*M_PI)/(2.*M_PI)) - 1);  // map to range (-PI, PI)
  double xc = (htPar[0] - 1./htPar[1])*sin(htPar[2]);
  double yc = -(htPar[0] - 1./htPar[1])*cos(htPar[2]);
  double st = 1./htPar[1]*(atan2(htPar[1]*(hitPos[1] - yc), htPar[1]*(hitPos[0] - xc)) - htPar[2] + M_PI/2.);
  if (st < 0)  // rotation angle has crossed the +/-PI border
    st += 2.*M_PI/fabs(htPar[1]);
  else if (st > 2.*M_PI/fabs(htPar[1]))
    st -= 2.*M_PI/fabs(htPar[1]);
  htPar[4] = atan2(st, hitPos[2] - htPar[3]);
}

int main(int argc, char* argv[])
{
  // Parameters
  const int nBins = 24;
  const double minAbsDoca = 0.;
  const double maxAbsDoca = 45.;
  const double minPt = 0.05;
  const double maxPt = 10.;
  const double minZ0 = -100.;
  const double maxZ0 = 100.;
  const double minTheta = 0.16;
  const double maxTheta = M_PI - 0.16;
  const double maxSt = 1000.;

  // Check arguments
  if (argc != 3) {
    cout << "Usage: checkHelixHoughTrans <n. of tracks> <n. of hits>" << endl;
    return 1;
  }

  // Get number of tracks
  int nTracks;
  TString nTracksStr(argv[1]);
  if (nTracksStr.IsDec())
    nTracks = nTracksStr.Atoi();
  else {
    cout << "First argument must be integer number. Abort." << endl;
    return 1;
  }
  
  // Get number of hits per track
  int nHits;
  TString nHitsStr(argv[2]);
  if (nHitsStr.IsDec())
    nHits = nHitsStr.Atoi();
  else {
    cout << "Second argument must be integer number. Abort." << endl;
    return 1;
  }
  
  // Create tree, with branches
  TTree tTrackPar("trackTree", "Fitted track parameters");
  vector<double> vDoca;
  vector<double> vKappa;
  vector<double> vPhi;
  vector<double> vZ0;
  vector<double> vTheta;
  vector<double> vAlgo;
  tTrackPar.Branch("doca", &vDoca);
  tTrackPar.Branch("kappa", &vKappa);
  tTrackPar.Branch("phi", &vPhi);
  tTrackPar.Branch("z0", &vZ0);
  tTrackPar.Branch("theta", &vTheta);
  tTrackPar.Branch("algo", &vAlgo);
  tTrackPar.SetDirectory(0);

  // Book 5D histogram
  double maxAbsKappa = 1.139e-3/minPt;
  TH5C hHoughVotes("hHoughVotes", "helix Hough transform votes", nBins, -maxAbsDoca, maxAbsDoca,
		   nBins, -maxAbsKappa, maxAbsKappa, nBins, -M_PI, M_PI, nBins, minZ0, maxZ0, nBins, minTheta, maxTheta);
  hHoughVotes.SetDirectory(0);
  // Set content to -128 for all histogram bins, to use full 8-bit range
  for (int iDoca = 1; iDoca <= nBins; iDoca++)
    for (int iKappa = 1; iKappa <= nBins; iKappa++)
      for (int iPhi = 1; iPhi <= nBins; iPhi++)
	for (int iZ0 = 1; iZ0 <= nBins; iZ0++)
	  for (int iTheta = 1; iTheta <= nBins; iTheta++)
	    hHoughVotes.SetBinContent(iDoca, iKappa, iPhi, iZ0, iTheta, -128);

  // Graph with helix points
  TGraph2D* grHelix;

  // Auxiliary arrays
  double hitPos[3];
  double htPar[5];

  // Loop on tracks
  srand((unsigned int)time(0));  // Initialise random generator seed
  for (int iTrack = 0; iTrack < nTracks; iTrack++) {

    // Generate random parameters
    double absDocaGen = minAbsDoca + (maxAbsDoca - minAbsDoca)*(double)rand()/RAND_MAX;
    double phiPoca = 2*M_PI*(double)rand()/RAND_MAX;
    double pt = exp(log(minPt) + log(maxPt/minPt)*(double)rand()/RAND_MAX);
    int charge = 2*(rand()%2) - 1;
    int sameSide = 2*(rand()%2) - 1;  // helix axis on same (+1) or opposite (-1) side wrt. POCA
    double z0Gen = minZ0 + (maxZ0 - minZ0)*(double)rand()/RAND_MAX;
    double thetaGen = minTheta + (maxTheta - minTheta)*(double)rand()/RAND_MAX;

    // Get helix axis
    double absKappaGen = 1.139e-3/pt;
    double xC = (absDocaGen + sameSide/absKappaGen)*cos(phiPoca);
    double yC = (absDocaGen + sameSide/absKappaGen)*sin(phiPoca);

    // Fill parameter vectors
    htPar[0]  = charge*sameSide*absDocaGen;
    vDoca.push_back(htPar[0]);
    htPar[1] = -charge*absKappaGen;
    vKappa.push_back(htPar[1]);
    double phiGen = phiPoca + charge*sameSide*M_PI/2.;
    if (phiGen >= M_PI)
      phiGen -= 2*M_PI;  // map to range (-PI, PI)
    vPhi.push_back(phiGen);
    htPar[3] = z0Gen;
    vZ0.push_back(htPar[3]);
    vTheta.push_back(thetaGen);
    vAlgo.push_back(-1);
    cout << endl << "Track n. " << iTrack << ": generated parameters = ( "
	 << vDoca.at(iTrack) << ", " << vKappa.at(iTrack) << ", " << vPhi.at(iTrack) << ", "
	 << vZ0.at(iTrack) << ", " << vTheta.at(iTrack) << ")" << endl;

    // Generate hits on track
    double x[nHits];
    double y[nHits];
    double z[nHits];
    for (int iHit = 0; iHit < nHits; iHit++) {
      double st = min(maxSt, M_PI/absKappaGen)*(double)rand()/RAND_MAX;
      double phiC = phiPoca + M_PI*(sameSide + 1)/2 - charge*absKappaGen*st;
      x[iHit] = xC + cos(phiC)/absKappaGen;
      y[iHit] = yC + sin(phiC)/absKappaGen;
      z[iHit] = z0Gen + st/tan(thetaGen);
      hitPos[0] = x[iHit];
      hitPos[1] = y[iHit];
      hitPos[2] = z[iHit];
      // Do the Hough transform to get back original parameters
      getPhiTheta(hitPos, htPar);
      cout << "Hit n." << iHit << ": phi = " << htPar[2] << "; theta = " << htPar[4] << endl;
    }
    
    // Populate helix graph
    if (iTrack == 0) {
      grHelix = new TGraph2D("grHelix", "helix trajectory", nHits, x, y, z);
      grHelix->SetDirectory(0);
    }

    // Loop on generated hits and fill vote histogram
    for (int iHit = 0; iHit < nHits; iHit++) {
      double r = sqrt(x[iHit]*x[iHit] + y[iHit]*y[iHit]);
      double docaBinWidth = 2*maxAbsDoca/nBins;
      double kappaBinWidth = 2*maxAbsKappa/nBins;
      double z0BinWidth = (maxZ0 - minZ0)/nBins;
      // Loop on allowed histogram bins (first doca, then kappa, then z0) and fill histogram
      for (htPar[0] = -maxAbsDoca + 0.5*docaBinWidth; htPar[0] < maxAbsDoca; htPar[0] += docaBinWidth) {
      // Start from first allowed kappa bin (must be kappaScan >= -2/(r + docaScan))
	double firstKappa = -maxAbsKappa;
	if (-2/(r + htPar[0]) > firstKappa) {
	  int binPosition = (-2/(r + htPar[0]))/kappaBinWidth;
	  firstKappa = (binPosition - 0.5)*kappaBinWidth;
	}
	for (htPar[1] = firstKappa; htPar[1] < min(maxAbsKappa, 2/(r - htPar[0])); htPar[1] += kappaBinWidth) {
	  for (htPar[3] = minZ0 + z0BinWidth/2; htPar[3] < maxZ0; htPar[3] += z0BinWidth) {
	    hitPos[0] = x[iHit];
	    hitPos[1] = y[iHit];
	    hitPos[2] = z[iHit];
	    getPhiTheta(hitPos, htPar);
	    hHoughVotes.Fill(htPar[0], htPar[1], htPar[2], htPar[3], htPar[4]);
	  }
	}
      }
    }
  }
  tTrackPar.Fill();

  // Write ROOT objects to file
  TFile outRootFile("houghCheck_genHits.root", "RECREATE");
  hHoughVotes.Write();
  tTrackPar.Write();
  grHelix->Write();

  return 0;
}
