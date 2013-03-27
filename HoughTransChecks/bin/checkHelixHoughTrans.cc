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

void getPhiLambda(double hitPos[], double docaScan, double kappaScan, double z0, double htPar[][2])
{
  double r = sqrt(hitPos[0]*hitPos[0] + hitPos[1]*hitPos[1]);
  double akappa = (2*docaScan + kappaScan*docaScan*docaScan + kappaScan*r*r)/(2*r);
  double hkappa = sqrt((docaScan*kappaScan + 1)*(docaScan*kappaScan + 1) - akappa*akappa);
  for (int iSol = 0; iSol < 2; iSol ++) {
    double kappax = (akappa*hitPos[0] - (2*iSol - 1)*hkappa*hitPos[1])/r;
    double kappay = (akappa*hitPos[1] + (2*iSol - 1)*hkappa*hitPos[0])/r;
    // Convert to perigee parameters
    double phiHit = atan2(hitPos[1], hitPos[0]);
    double phiD = atan2(kappay, kappax);
    int rot = -2*(int((phiHit - phiD + 2*M_PI)/M_PI)%2) + 1;  // hit position wrt. poca: +1 = anticlockwise; -1 = clockwise
    htPar[0][iSol] = rot*docaScan;
    htPar[1][iSol] = -rot*kappaScan;
    if (fabs(htPar[1][iSol]) < 1.e-6)
      htPar[1][iSol] = (htPar[1][iSol] >= 0) ? 1.e-6 : -1.e-6;  // avoid kappa = 0 (maximum curvature radius = 1 km)
    htPar[2][iSol] = phiD + rot*(M_PI/2.);
    htPar[2][iSol] += -2.*M_PI*(int((htPar[2][iSol] + 3.*M_PI)/(2.*M_PI)) - 1);  // map to range (-PI, PI)
    double xc = (htPar[0][iSol] - 1./htPar[1][iSol])*sin(htPar[2][iSol]);
    double yc = -(htPar[0][iSol] - 1./htPar[1][iSol])*cos(htPar[2][iSol]);
    double st = 1./htPar[1][iSol]*(atan2(htPar[1][iSol]*(hitPos[1] - yc), htPar[1][iSol]*(hitPos[0] - xc)) - htPar[2][iSol] + M_PI/2.);
    if (st < 0)  // rotation angle has crossed the +/-PI border
      st += 2.*M_PI/fabs(htPar[1][iSol]);
    else if (st > 2.*M_PI/fabs(htPar[1][iSol]))
      st -= 2.*M_PI/fabs(htPar[1][iSol]);
    htPar[3][iSol] = z0;
    htPar[4][iSol] = atan((hitPos[2] - z0)/st);
  }
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
  const double minLambda = -1.4;
  const double maxLambda = 1.4;
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
  vector<double> vLambda;
  vector<double> vAlgo;
  tTrackPar.Branch("doca", &vDoca);
  tTrackPar.Branch("kappa", &vKappa);
  tTrackPar.Branch("phi", &vPhi);
  tTrackPar.Branch("z0", &vZ0);
  tTrackPar.Branch("lambda", &vLambda);
  tTrackPar.Branch("algo", &vAlgo);
  tTrackPar.SetDirectory(0);

  // Book 5D histogram
  double maxAbsKappa = 1.139e-3/minPt;
  TH5C hHoughVotes("hHoughVotes", "helix Hough transform votes", nBins, -maxAbsDoca, maxAbsDoca,
		   nBins, -maxAbsKappa, maxAbsKappa, nBins, -M_PI, M_PI, nBins, minZ0, maxZ0, nBins, minLambda, maxLambda);
  hHoughVotes.SetDirectory(0);
  // Set content to -128 for all histogram bins, to use full 8-bit range
  for (int iDoca = 1; iDoca <= nBins; iDoca++)
    for (int iKappa = 1; iKappa <= nBins; iKappa++)
      for (int iPhi = 1; iPhi <= nBins; iPhi++)
	for (int iZ0 = 1; iZ0 <= nBins; iZ0++)
	  for (int iLambda = 1; iLambda <= nBins; iLambda++)
	    hHoughVotes.SetBinContent(iDoca, iKappa, iPhi, iZ0, iLambda, -128);

  // Graph with helix points
  TGraph2D* grHelix;

  // Auxiliary arrays
  double hitPos[3];
  double htPar[5][2];

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
    double lambdaGen = minLambda + (maxLambda - minLambda)*(double)rand()/RAND_MAX;

    // Get helix axis
    double absKappaGen = 1.139e-3/pt;
    double xC = (absDocaGen + sameSide/absKappaGen)*cos(phiPoca);
    double yC = (absDocaGen + sameSide/absKappaGen)*sin(phiPoca);

    // Fill parameter vectors
    vDoca.push_back(charge*sameSide*absDocaGen);
    vKappa.push_back(-charge*absKappaGen);
    double phiGen = phiPoca + charge*sameSide*M_PI/2.;
    if (phiGen >= M_PI)
      phiGen -= 2*M_PI;  // map to range (-PI, PI)
    vPhi.push_back(phiGen);
    vZ0.push_back(z0Gen);
    vLambda.push_back(lambdaGen);
    vAlgo.push_back(-1);
    cout << endl << "Track n. " << iTrack << ": generated parameters = ( "
	 << vDoca.at(iTrack) << ", " << vKappa.at(iTrack) << ", " << vPhi.at(iTrack) << ", "
	 << vZ0.at(iTrack) << ", " << vLambda.at(iTrack) << ")" << endl;

    // Generate hits on track
    double x[nHits];
    double y[nHits];
    double z[nHits];
    for (int iHit = 0; iHit < nHits; iHit++) {
      double st = min(maxSt, M_PI/absKappaGen)*(double)rand()/RAND_MAX;
      double phiC = phiPoca + M_PI*(sameSide + 1)/2 - charge*absKappaGen*st;
      x[iHit] = xC + cos(phiC)/absKappaGen;
      y[iHit] = yC + sin(phiC)/absKappaGen;
      z[iHit] = z0Gen + st*tan(lambdaGen);
      hitPos[0] = x[iHit];
      hitPos[1] = y[iHit];
      hitPos[2] = z[iHit];
      // Do the Hough transform to get back original parameters
      getPhiLambda(hitPos, absDocaGen, sameSide*absKappaGen, z0Gen, htPar);
      cout << "Hit n." << iHit << ": phi = [" << htPar[2][0] << ", " << htPar[2][1] << "]"
	   << "; lambda = [" << htPar[4][0] << ", " << htPar[4][1] << "]" << endl;
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
      for (double docaScan = -maxAbsDoca + 0.5*docaBinWidth; docaScan < maxAbsDoca; docaScan += docaBinWidth) {
      // Start from first allowed kappa bin (must be kappaScan >= -2/(r + docaScan))
	double firstKappa = -maxAbsKappa;
	if (-2/(r + docaScan) > firstKappa) {
	  int binPosition = (-2/(r + docaScan))/kappaBinWidth;
	  firstKappa = (binPosition - 0.5)*kappaBinWidth;
	}
	for (double kappaScan = firstKappa; kappaScan < min(maxAbsKappa, 2/(r - docaScan)); kappaScan += kappaBinWidth) {
	  for (double z0 = minZ0 + z0BinWidth/2; z0 < maxZ0; z0 += z0BinWidth) {
	    hitPos[0] = x[iHit];
	    hitPos[1] = y[iHit];
	    hitPos[2] = z[iHit];
	    getPhiLambda(hitPos, docaScan, kappaScan, z0, htPar);
	    for (int iSol = 0; iSol < 2; iSol++)
	      hHoughVotes.Fill(htPar[0][iSol], htPar[1][iSol], htPar[2][iSol], htPar[3][iSol], htPar[4][iSol]);
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
