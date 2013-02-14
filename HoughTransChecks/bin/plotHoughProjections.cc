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

using namespace std;

int main(int argc, char* argv[])
{
  // Check arguments
  if (argc != 2) {
    cout << "Usage: plotHoughProjections <input file name>" << endl;
    return 1;
  }
  char* inFileName = argv[1];
  
  // Open input file and get objects
  TFile* inFile = new TFile(inFileName);
  if (inFile == 0) {
    cout << "Error opening input file " << inFileName << ". Abort." << endl;
    return 1;
  }
  TH3S* hHoughVotes = (TH3S*)(inFile->Get("hHoughVotes"));
  if (hHoughVotes == 0) {
    cout << "Error getting votes histogram. Abort." << endl;
    return 1;
  }
  TTree* tTrackPar = (TTree*)(inFile->Get("trackTree"));
  if (tTrackPar == 0) {
    cout << "Error getting track parameter tree. Abort." << endl;
    return 1;
  }

  // Build all-events vectors from tree (per-event) entries
  vector<double> vDoca, vKappa, vPhi;
  vector<double>* vDocaEv = 0;
  vector<double>* vKappaEv = 0;
  vector<double>* vPhiEv = 0;
  TBranch* bDocaEv = 0;
  TBranch* bKappaEv = 0;
  TBranch* bPhiEv = 0;
  tTrackPar->SetBranchAddress("doca", &vDocaEv, &bDocaEv);
  tTrackPar->SetBranchAddress("kappa", &vKappaEv, &bKappaEv);
  tTrackPar->SetBranchAddress("phi", &vPhiEv, &bPhiEv);
  // Loop on tree
  for (int iEv = 0; iEv < tTrackPar->GetEntries(); iEv++) {
    long tEntry = tTrackPar->LoadTree(iEv);
    tTrackPar->GetEntry(tEntry);
    for (int iTk = 0; iTk < vDocaEv->size(); iTk++) {
      vDoca.push_back(vDocaEv->at(iTk));
      vKappa.push_back(vKappaEv->at(iTk));
      vPhi.push_back(vPhiEv->at(iTk));
    }
  }

  // Set up graphics
  TString psFileName = "houghVoteProjections.ps";
  TString openPsStr = psFileName + "[";
  TString closePsStr = psFileName + "]";
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0001);
  gStyle->SetOptStat(0);
  double maxVote = hHoughVotes->GetMaximum();
  TCanvas cHoughProj("cHoughProj", "", 800, 1200);
  TH2S* hProj[6];
  TString hName, hTitle;
  float rDoca =0.05*(hHoughVotes->GetXaxis()->GetXmax() - hHoughVotes->GetXaxis()->GetXmin());
  float rKappa =0.05*(hHoughVotes->GetYaxis()->GetXmax() - hHoughVotes->GetYaxis()->GetXmin());
  float rPhi =0.05*(hHoughVotes->GetZaxis()->GetXmax() - hHoughVotes->GetZaxis()->GetXmin());
  cHoughProj.Print(openPsStr);

  // phi slices
  int iPad = 0;
  cHoughProj.Divide(2, 3);
  for (int iBin = 1; iBin <= hHoughVotes->GetNbinsZ(); iBin++) {
    cHoughProj.cd(++iPad);
    float phi = hHoughVotes->GetZaxis()->GetBinCenter(iBin);
    stringstream phiStr;
    phiStr << setprecision(4) << phi;
    hName = TString("hDocaKappa_phi") + phiStr.str().c_str();
    hTitle = TString("kappa vs. doca, phi = ") + phiStr.str().c_str();
    hHoughVotes->GetZaxis()->SetRange(iBin, iBin);
    hProj[iPad - 1] = (TH2S*)(hHoughVotes->Project3D("yx"));
    hProj[iPad - 1]->SetNameTitle(hName, hTitle);
    hProj[iPad - 1]->SetMinimum(0.);
    hProj[iPad - 1]->SetMaximum(maxVote);
    hProj[iPad - 1]->Draw("COLZ");
    // Draw circles around fitted track parameters
    for (int iTkPar = 0; iTkPar < vPhi.size(); iTkPar++) {
      if (vPhi.at(iTkPar) > hHoughVotes->GetZaxis()->GetBinLowEdge(iBin) &&
	  vPhi.at(iTkPar) <= hHoughVotes->GetZaxis()->GetBinUpEdge(iBin)) {
	TEllipse* circle = new TEllipse(vDoca.at(iTkPar), vKappa.at(iTkPar), rDoca, rKappa);
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
  cHoughProj.Update();
  cHoughProj.Print(psFileName);
  hHoughVotes->GetZaxis()->SetRange(1, hHoughVotes->GetNbinsZ());

  // kappa slices
  iPad = 0;
  cHoughProj.Clear();
  cHoughProj.Divide(2, 3);
  for (int iBin = 1; iBin <= hHoughVotes->GetNbinsY(); iBin++) {
    cHoughProj.cd(++iPad);
    float kappa = hHoughVotes->GetYaxis()->GetBinCenter(iBin);
    stringstream kappaStr;
    kappaStr << setprecision(4) << kappa;
    hName = TString("hDocaPhi_kappa") + kappaStr.str().c_str();
    hTitle = TString("phi vs. doca, kappa = ") + kappaStr.str().c_str();
    hHoughVotes->GetYaxis()->SetRange(iBin, iBin);
    hProj[iPad - 1] = (TH2S*)(hHoughVotes->Project3D("zx"));
    hProj[iPad - 1]->SetNameTitle(hName, hTitle);
    hProj[iPad - 1]->SetMinimum(0.);
    hProj[iPad - 1]->SetMaximum(maxVote);
    hProj[iPad - 1]->Draw("COLZ");
    // Draw circles around fitted track parameters
    for (int iTkPar = 0; iTkPar < vKappa.size(); iTkPar++) {
      if (vKappa.at(iTkPar) > hHoughVotes->GetYaxis()->GetBinLowEdge(iBin) &&
	  vKappa.at(iTkPar) <= hHoughVotes->GetYaxis()->GetBinUpEdge(iBin)) {
	TEllipse* circle = new TEllipse(vDoca.at(iTkPar), vPhi.at(iTkPar), rDoca, rPhi);
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
  cHoughProj.Update();
  cHoughProj.Print(psFileName);
  hHoughVotes->GetYaxis()->SetRange(1, hHoughVotes->GetNbinsY());

  // phi slices
  iPad = 0;
  cHoughProj.Clear();
  cHoughProj.Divide(2, 3);
  for (int iBin = 1; iBin <= hHoughVotes->GetNbinsX(); iBin++) {
    cHoughProj.cd(++iPad);
    float doca = hHoughVotes->GetXaxis()->GetBinCenter(iBin);
    stringstream docaStr;
    docaStr << setprecision(4) << doca;
    hName = TString("hKappaPhi_doca") + docaStr.str().c_str();
    hTitle = TString("phi vs. kappa, doca = ") + docaStr.str().c_str();
    hHoughVotes->GetXaxis()->SetRange(iBin, iBin);
    hProj[iPad - 1] = (TH2S*)(hHoughVotes->Project3D("zy"));
    hProj[iPad - 1]->SetNameTitle(hName, hTitle);
    hProj[iPad - 1]->SetMinimum(0.);
    hProj[iPad - 1]->SetMaximum(maxVote);
    hProj[iPad - 1]->Draw("COLZ");
    // Draw circles around fitted track parameters
    for (int iTkPar = 0; iTkPar < vDoca.size(); iTkPar++) {
      if (vDoca.at(iTkPar) > hHoughVotes->GetXaxis()->GetBinLowEdge(iBin) &&
	  vDoca.at(iTkPar) <= hHoughVotes->GetXaxis()->GetBinUpEdge(iBin)) {
	TEllipse* circle = new TEllipse(vKappa.at(iTkPar), vPhi.at(iTkPar), rKappa, rPhi);
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
  cHoughProj.Update();
  cHoughProj.Print(psFileName);
  hHoughVotes->GetXaxis()->SetRange(1, hHoughVotes->GetNbinsX());

  cHoughProj.Print(closePsStr);
}
