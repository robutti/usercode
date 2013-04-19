#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <math.h>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH2F.h>
#include <TProfile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>

#define M_PI 3.14159265358979323846

using namespace std;

int main(int argc, char* argv[])
{
  // Check arguments
  if (argc != 2) {
    cout << "Usage: plotParamRes <input file name>" << endl;
    return 1;
  }
  char* inFileName = argv[1];

 
  // Open input file and get objects
  TFile* inFile = new TFile(inFileName);
  if (inFile == 0) {
    cout << "Error opening input file " << inFileName << ". Abort." << endl;
    return 1;
  }

  // Get tree and set branches
  TTree* tHelixPar = (TTree*)(inFile->Get("trackTree"));
  vector<double>* vDoca = 0;
  vector<vector <double> >* vDocaExp = 0;
  vector<double>* vKappa = 0;
  vector<vector <double> >* vKappaExp = 0;
  vector<double>* vPhi = 0;
  vector<vector <double> >* vPhiExp = 0;
  vector<double>* vZ0 = 0;
  vector<vector <double> >* vZ0Exp = 0;
  vector<double>* vTheta = 0;
  vector<vector <double> >* vThetaExp = 0;
  TBranch* bDoca;
  TBranch* bDocaExp;
  TBranch* bKappa;
  TBranch* bKappaExp;
  TBranch* bPhi;
  TBranch* bPhiExp;
  TBranch* bZ0;
  TBranch* bZ0Exp;
  TBranch* bTheta;
  TBranch* bThetaExp;
  tHelixPar->SetBranchAddress("doca", &vDoca, &bDoca);
  tHelixPar->SetBranchAddress("docaExp", &vDocaExp, &bDocaExp);
  tHelixPar->SetBranchAddress("kappa", &vKappa, &bKappa);
  tHelixPar->SetBranchAddress("kappaExp", &vKappaExp, &bKappaExp);
  tHelixPar->SetBranchAddress("phi", &vPhi, &bPhi);
  tHelixPar->SetBranchAddress("phiExp", &vPhiExp, &bPhiExp);
  tHelixPar->SetBranchAddress("z0", &vZ0, &bZ0);
  tHelixPar->SetBranchAddress("z0Exp", &vZ0Exp, &bZ0Exp);
  tHelixPar->SetBranchAddress("theta", &vTheta, &bTheta);
  tHelixPar->SetBranchAddress("thetaExp", &vThetaExp, &bThetaExp);

  // Book histograms
  TH2F hDoca("hDoca", "HT vs. fitted doca", 100, -100., 100., 100, -100., 100.);
  //  TH2F hKappa("hKappa", "HT vs. fitted kappa", 100, -0.006, 0.006, 100, -0.006, 0.006);
  TH2F hSqrtKappa("hSqrtKappa", "HT vs. fitted sqrt(kappa)", 100, -0.075, 0.075, 100, -0.075, 0.075);
  TH2F hPhi("hPhi", "HT vs. fitted phi", 100, -M_PI, M_PI, 100, -M_PI, M_PI);
  TH2F hZ0("hZ0", "HT vs. fitted z0", 100, -1000., 1000., 100, -1000., 1000);
  TH2F hEta("hEta", "HT vs. fitted eta", 100, -2.5, 2.5, 100, -2.5, 2.5);
  TProfile profDeltaDoca("profDeltaDoca", "HT - fitted doca", 40, -100., 100., "s");
  //  TProfile profDeltaKappa("profDeltaKappa", "HT - fitted kappa", 40, -0.006, 0.006, "s");
  TProfile profDeltaSqrtKappa("profDeltaSqrtKappa", "HT - fitted sqrt(kappa)", 40, -0.075, 0.075, "s");
  TProfile profDeltaPhi("profDeltaPhi", "HT - fitted phi", 40, -M_PI, M_PI, "s");
  TProfile profDeltaZ0("profDeltaZ0", "HT - fitted z0", 40, -1000., 1000., "s");
  TProfile profDeltaEta("profDeltaEta", "HT - fitted eta", 40, -2.5, 2.5, "s");
  
  // Loop on tree and fill histograms
  for (int iEv = 0; iEv < tHelixPar->GetEntries(); iEv++) {
    long tEntry = tHelixPar->LoadTree(iEv);
    tHelixPar->GetEntry(tEntry);
    for (int iTrack = 0; iTrack < vDoca->size(); iTrack++) {
      double doca = vDoca->at(iTrack);
      double kappa = vKappa->at(iTrack);
      double phi = vPhi->at(iTrack);
      double z0 = vZ0->at(iTrack);
      double eta = -log(tan((vTheta->at(iTrack))/2.));
      for (int iHit = 0; iHit < vDocaExp->at(iTrack).size(); iHit++) {
	hDoca.Fill(doca, vDocaExp->at(iTrack).at(iHit));
	if (fabs(vDocaExp->at(iTrack).at(iHit) - doca) < 100.)
	  profDeltaDoca.Fill(doca, vDocaExp->at(iTrack).at(iHit) - doca);
	//	hKappa.Fill(kappa, vKappaExp->at(iTrack).at(iHit));
	if (fabs(vKappaExp->at(iTrack).at(iHit)/sqrt(fabs(vKappaExp->at(iTrack).at(iHit))) - kappa/sqrt(fabs(kappa))) < 0.02) {
	  hSqrtKappa.Fill(kappa/sqrt(fabs(kappa)), vKappaExp->at(iTrack).at(iHit)/sqrt(fabs(vKappaExp->at(iTrack).at(iHit))));
	  //	  profDeltaKappa.Fill(kappa, vKappaExp->at(iTrack).at(iHit) - kappa);
	  profDeltaSqrtKappa.Fill(kappa/fabs(kappa)*sqrt(fabs(kappa)), vKappaExp->at(iTrack).at(iHit)/fabs(vKappaExp->at(iTrack).at(iHit))*sqrt(fabs(vKappaExp->at(iTrack).at(iHit))) - kappa/fabs(kappa)*sqrt(fabs(kappa)));
	}
	hPhi.Fill(phi, vPhiExp->at(iTrack).at(iHit));
	if (fabs(vPhiExp->at(iTrack).at(iHit) - phi) < 0.5)
	  profDeltaPhi.Fill(phi, vPhiExp->at(iTrack).at(iHit) - phi);
	hZ0.Fill(z0, vZ0Exp->at(iTrack).at(iHit));
	if (fabs(vZ0Exp->at(iTrack).at(iHit) - z0) < 500)
	  profDeltaZ0.Fill(z0, vZ0Exp->at(iTrack).at(iHit) - z0);
	hEta.Fill(eta, -log(tan((vThetaExp->at(iTrack).at(iHit))/2.)));
	if (fabs(-log(tan(vThetaExp->at(iTrack).at(iHit)/2.)) - eta) < 0.5)
	  profDeltaEta.Fill(eta, -log(tan((vThetaExp->at(iTrack).at(iHit))/2.)) - eta);
      }
    }
  }

  // Set up graphics
  TString psFileName = "houghParamRes.ps";
  TString openPsStr = psFileName + "[";
  TString closePsStr = psFileName + "]";
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0001);
  gStyle->SetOptStat(0);
  TCanvas cParamRes("cParamRes", "", 900, 1500);

  // Plot all histograms (fill Delta histograms)
  cParamRes.Print(openPsStr);
  cParamRes.Divide(3, 5);
  cParamRes.cd(1);
  hDoca.Draw("COLZ"); 
  cParamRes.cd(2);
  TProfile* profDoca = hDoca.ProfileX("profDoca", 1, -1, "s");
  profDoca->Draw();
  cParamRes.cd(3);
  profDeltaDoca.SetAxisRange(-40., 40., "Y");
  profDeltaDoca.Draw();
  cParamRes.cd(4);
  //  hKappa.Draw("COLZ");  
  hSqrtKappa.Draw("COLZ");  
  cParamRes.cd(5);
  //  TProfile* profKappa = hKappa.ProfileX("profKappa", 1, -1, "s");
  TProfile* profSqrtKappa = hSqrtKappa.ProfileX("profSqrtKappa", 1, -1, "s");
  //  profKappa->Draw();  
  profSqrtKappa->Draw();  
  cParamRes.cd(6);
  //  profDeltaKappa.SetAxisRange(-1.e-3, 1.e-3, "Y");
  profDeltaSqrtKappa.SetAxisRange(-5.e-3, 5.e-3, "Y");
  //  profDeltaKappa.Draw();
  profDeltaSqrtKappa.Draw();
  cParamRes.cd(7);
  hPhi.Draw("COLZ");  
  cParamRes.cd(8);
  TProfile* profPhi = hPhi.ProfileX("profPhi", 1, -1, "s");
  profPhi->Draw();
  cParamRes.cd(9);
  profDeltaPhi.SetAxisRange(-0.2, 0.2, "Y");
  profDeltaPhi.Draw();
  cParamRes.cd(10);
  hZ0.Draw("COLZ");  
  cParamRes.cd(11);
  TProfile* profZ0 = hZ0.ProfileX("profZ0", 1, -1, "s");
  profZ0->Draw();  
  cParamRes.cd(12);
  profDeltaZ0.SetAxisRange(-500., 500., "Y");
  profDeltaZ0.Draw();
  cParamRes.cd(13);
  hEta.Draw("COLZ");  
  cParamRes.cd(14);
  TProfile* profEta = hEta.ProfileX("profEta", 1, -1, "s");
  profEta->Draw();
  cParamRes.cd(15);
  profDeltaEta.SetAxisRange(-0.6, 0.6, "Y");
  profDeltaEta.Draw();

  cParamRes.Update();
  cParamRes.Print(psFileName);
  cParamRes.Print(closePsStr);

  return 0;
}
