#include <iostream>
#include <sstream>
#include <iomanip>
#include <vector>

#include <TFile.h>
#include <TString.h>
#include <TH3S.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH2S.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include "../interface/TH5.h"

using namespace std;

int main(int argc, char* argv[])
{
  // Check arguments
  if (argc != 4) {
    cout << "Usage: plotHoughProjections <input file name> <parX> <parY>" << endl;
    return 1;
  }
  char* inFileName = argv[1];
  
  // Open input file and get objects
  TFile* inFile = new TFile(inFileName);
  if (inFile == 0) {
    cout << "Error opening input file " << inFileName << ". Abort." << endl;
    return 1;
  }
  TH5C* hHoughVotes = (TH5C*)(inFile->Get("hHoughVotes"));
  if (hHoughVotes == 0) {
    cout << "Error getting votes histogram. Abort." << endl;
    return 1;
  }
  TTree* tTrackPar = (TTree*)(inFile->Get("trackTree"));
  if (tTrackPar == 0) {
    cout << "Error getting track parameter tree. Abort." << endl;
    return 1;
  }

  // Get projection parameters
  TString strParX(argv[2]);
  TString strParY(argv[3]);
  int iParX = -1, iParY = -1;
  TString parName[5] = {"doca", "kappa", "phi", "z0", "lambda"};
  TString axisName[5] = {"x", "y", "z", "u", "v"};
  for (int iPar = 0; iPar < 5; iPar++) {
    if (strParX.Atoi() == iPar || strParX == parName[iPar] || strParX == axisName[iPar])
      iParX = iPar;
    if (strParY.Atoi() == iPar || strParY == parName[iPar] || strParY == axisName[iPar])
      iParY = iPar;
  }
  if (iParX < 0 || iParX > 4 || iParY < 0 || iParY > 4 || iParX == iParY) {
    cout << "Invalid paramater specification. Abort." << endl;
    return 2;
  }

  // Build all-events vectors from tree (per-event) entries
  vector<double> vPar[5];
  vector<double>* vDocaEv = 0;
  vector<double>* vKappaEv = 0;
  vector<double>* vPhiEv = 0;
  vector<double>* vZ0Ev = 0;
  vector<double>* vLambdaEv = 0;
  TBranch* bDocaEv = 0;
  TBranch* bKappaEv = 0;
  TBranch* bPhiEv = 0;
  TBranch* bZ0Ev = 0;
  TBranch* bLambdaEv = 0;
  tTrackPar->SetBranchAddress("doca", &vDocaEv, &bDocaEv);
  tTrackPar->SetBranchAddress("kappa", &vKappaEv, &bKappaEv);
  tTrackPar->SetBranchAddress("phi", &vPhiEv, &bPhiEv);
  tTrackPar->SetBranchAddress("z0", &vZ0Ev, &bZ0Ev);
  tTrackPar->SetBranchAddress("lambda", &vLambdaEv, &bLambdaEv);
  // Loop on tree
  for (int iEv = 0; iEv < tTrackPar->GetEntries(); iEv++) {
    long tEntry = tTrackPar->LoadTree(iEv);
    tTrackPar->GetEntry(tEntry);
    for (unsigned int iTk = 0; iTk < vDocaEv->size(); iTk++) {
      vPar[0].push_back(vDocaEv->at(iTk));
      vPar[1].push_back(vKappaEv->at(iTk));
      vPar[2].push_back(vPhiEv->at(iTk));
      vPar[3].push_back(vZ0Ev->at(iTk));
      vPar[4].push_back(vLambdaEv->at(iTk));
    }
  }

  // Index axes
  TAxis* axis[5];
  axis[0] = hHoughVotes->GetXaxis();
  axis[1] = hHoughVotes->GetYaxis();
  axis[2] = hHoughVotes->GetZaxis();
  axis[3] = hHoughVotes->GetUaxis();
  axis[4] = hHoughVotes->GetVaxis();
  int iProjAxis[3];
  unsigned int nProj = 0;
  for (int iAxis = 0; iAxis < 5; iAxis++)
    if (iAxis != iParX && iAxis != iParY) {
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
  double maxVote = hHoughVotes->GetMaximum() + 128;
  TCanvas cHoughProj("cHoughProj", "", 800, 1200);
  TH2S* hProj[6];
  TString hName, hTitle;
  float rX = 0.05*(axis[iParX]->GetXmax() - axis[iParX]->GetXmin());
  float rY = 0.05*(axis[iParY]->GetXmax() - axis[iParY]->GetXmin());
  cHoughProj.Print(openPsStr);

  //
  TString projName = axisName[iParY] + axisName[iParX];
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
	hName = TString("h") + parName[iParX] + parName[iParY] + 
	  "-" + parName[iProjAxis[0]] + par0Str.str().c_str() + 
	  "-" + parName[iProjAxis[1]] + par1Str.str().c_str() + 
	  "-" + parName[iProjAxis[2]] + par2Str.str().c_str();
	hTitle = parName[iParY] + " vs. " + parName[iParX] + ", " + 
	  parName[iProjAxis[0]] + " = " + par0Str.str().c_str() + ", " +
	  parName[iProjAxis[1]] + " = " + par1Str.str().c_str() + ", " +
	  parName[iProjAxis[2]] + " = " + par2Str.str().c_str();
	axis[iProjAxis[0]]->SetRange(iBin0, iBin0);
	axis[iProjAxis[1]]->SetRange(iBin1, iBin1);
	axis[iProjAxis[2]]->SetRange(iBin2, iBin2);
	hProj[iPad - 1] = (TH2S*)(hHoughVotes->Project5D(projName));
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
	    TEllipse* circle = new TEllipse(vPar[iParX].at(iTkPar), vPar[iParY].at(iTkPar), rX, rY);
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

  return 0;
}
