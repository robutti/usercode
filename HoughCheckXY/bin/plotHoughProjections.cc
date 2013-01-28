#include <iostream>
#include <sstream>
#include <iomanip>

#include <TFile.h>
#include <TString.h>
#include <TH3S.h>
#include <TH2S.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

int main(int argc, char* argv[])
{
  // Check arguments
  if (argc != 2) {
    cout << "Usage: plotHoughProjections <input file name>" << endl;
    return 1;
  }
  char* inFileName = argv[1];
  
  // Open input file and get 3D histogram
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
  
  // Set up graphics
  TString psFileName = "houghVoteProjections.ps";
  TString openPsStr = psFileName + "[";
  TString closePsStr = psFileName + "]";
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0001);
  gStyle->SetOptStat(0);
  TCanvas cHoughProj("cHoughProj", "", 800, 1200);
  TH2S* hProj[6];
  TString hName, hTitle;
  
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
    hProj[iPad - 1]->Draw("COLZ");
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
    hProj[iPad - 1]->Draw("COLZ");
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
    hProj[iPad - 1]->Draw("COLZ");
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
