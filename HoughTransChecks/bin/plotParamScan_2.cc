#include <iostream>
#include <dirent.h>
#include <stdlib.h>
#include <math.h>
#include <map>

#include <TGraph.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TPad.h>

using namespace std;

int main(int argc, char* argv[])
{
  // Check arguments
  if (argc != 2) {
    cout << "Usage: plotParamScan <input directory name>" << endl;
    return 1;
  }
  char* dirName = argv[1];

  // Tree variables
  unsigned int nSeeds, nGoodSeeds, nTriedSeeds, nAssSeeds;
  unsigned int nTracks, nTracksAlgo, nAccTracks, nAssTracks, nFoundTracks;
  TBranch *bNSeeds, *bNGoodSeeds, *bNTriedSeeds, *bNAssSeeds;
  TBranch *bNTracks, *bNTracksAlgo, *bNAccTracks, *bNAssTracks, *bNFoundTracks;;

  // Graphs
  TGraph grSeedEff[5];
  TGraph grTrackEff[5];
  
  // Open the directory and look for ROOT files
  DIR* dir = opendir(dirName);
  if (dir == NULL) {
    cout << "Error opening directory " << dirName << ". Abort." << endl;
    return 1;
  }
  // First find default values for nBins
  unsigned int defaultNBins[5]= {0, 0, 0, 0, 0};
  map<unsigned int, bool> foundNBins[5];
  unsigned int nFoundDef = 0;
  while (struct dirent*  dirEntry = readdir(dir)) {
    string fileName = dirEntry->d_name;
    if (fileName.length() > 5 && fileName.compare(fileName.length() - 5, 5, string(".root")) == 0) {
      // Get run parameters
      int iParStr = fileName.find_last_of(string("_"));
      if (iParStr == string::npos) {
        cout << "Incorrect file name " << fileName << ". Abort." << endl;
        return 2;
      }
      for (int iPar = 0; iPar < 5; iPar++) {
        iParStr++;
        int iValStr = fileName.find_first_not_of("0123456789", iParStr);
        if (iValStr == string::npos) {
          cout << "Incorrect file name " << fileName << ". Abort." << endl;
          return 2;
        }
        unsigned int nBins = atoi((fileName.substr(iParStr, iValStr - iParStr)).c_str());
        if (foundNBins[iPar].find(nBins) == foundNBins[iPar].end())
          foundNBins[iPar][nBins] = true;
        else if (defaultNBins[iPar] == 0) {
          defaultNBins[iPar] = nBins;
          nFoundDef++;
          if (nFoundDef == 5)
            break;
        }
        iParStr = iValStr;
      }
      if (nFoundDef == 5)
        break;
    }
  }   
  // Now re-scan the directory and get data from trees
  unsigned int nBins[5];
  while (struct dirent*  dirEntry = readdir(dir)) {
    string fileName = dirEntry->d_name;
    if (fileName.length() > 5 && fileName.compare(fileName.length() - 5, 5, string(".root")) == 0) {
      string filePath = string(dirName) + "/" + fileName;
      // Get run parameters
      int iParStr = fileName.find_last_of(string("_"));
      if (iParStr == string::npos) {
        cout << "Incorrect file name " << fileName << ". Abort." << endl;
        return 2;
      }
      for (int iPar = 0; iPar < 5; iPar++) {
        iParStr++;
        int iValStr = fileName.find_first_not_of("0123456789", iParStr);
        if (iValStr == string::npos) {
          cout << "Incorrect file name " << fileName << ". Abort." << endl;
          return 2;
        }
        nBins[iPar] = atoi((fileName.substr(iParStr, iValStr - iParStr)).c_str());
        iParStr = iValStr;
      }
      // Get run statistics from tree
      TFile fHT(filePath.c_str());
      if (!fHT.IsOpen()) {
        cout << "Error opening file " << filePath << ". Skipping." << endl;
        continue;
      }
      TTree* tHT = (TTree*)(fHT.Get("trackTree"));
      if (tHT == NULL) {
        cout << "Error getting HT tree from file" << filePath << ". Skipping." << endl;
        continue;
      }
      tHT->SetBranchAddress("nSeeds", &nSeeds, &bNSeeds);
      tHT->SetBranchAddress("nGoodSeeds", &nGoodSeeds, &bNGoodSeeds);
      tHT->SetBranchAddress("nTriedSeeds", &nTriedSeeds, &bNTriedSeeds);
      tHT->SetBranchAddress("nAssSeeds", &nAssSeeds, &bNAssSeeds);
      tHT->SetBranchAddress("nTracks", &nTracks, &bNTracks);
      tHT->SetBranchAddress("nTracksAlgo", &nTracksAlgo, &bNTracksAlgo);
      tHT->SetBranchAddress("nAccTracks", &nAccTracks, &bNAccTracks);
      tHT->SetBranchAddress("nAssTracks", &nAssTracks, &bNAssTracks);
      tHT->SetBranchAddress("nFoundTracks", &nFoundTracks, &bNFoundTracks);

      // Loop on tree
      unsigned int nTriedSeedsTot = 0;
      unsigned int nAssSeedsTot = 0;
      unsigned int nAssTracksTot = 0;
      unsigned int nFoundTracksTot = 0;
      for (int iEv = 0; iEv < tHT->GetEntries(); iEv++) {
        long tEntry = tHT->LoadTree(iEv);
        tHT->GetEntry(tEntry);
        nTriedSeedsTot += nTriedSeeds;
        nAssSeedsTot += nAssSeeds;
        nAssTracksTot += nAssTracks;
        nFoundTracksTot += nFoundTracks;
      }
      //      cout << "Seed efficiency = " << (float)nAssSeeds/nTriedSeeds << "; Track efficiency = " << (float)nFoundTracks/nAssTracks << endl;
      for (int iPar = 0; iPar < 5; iPar++) {
        bool validPoint = true;
        for (int iOtherPar = 1; iOtherPar <= 4; iOtherPar++) {
          int otherPar = (iPar + iOtherPar)%5;
          if (nBins[otherPar] != defaultNBins[otherPar]) {
            validPoint = false;
            break;
          }
        }
        if (validPoint) {
          grSeedEff[iPar].SetPoint(grSeedEff[iPar].GetN(), nBins[iPar], (float)nAssSeeds/nTriedSeeds);
          grTrackEff[iPar].SetPoint(grSeedEff[iPar].GetN(), nBins[iPar], (float)nFoundTracks/nAssTracks);
        }
      }
    }
  }

  // Set up graphics
  TString psFileName = "paramBinScan_2.ps";
  TString openPsStr = psFileName + "[";
  TString closePsStr = psFileName + "]";
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0001);
  gStyle->SetOptStat(0);
  TCanvas cParamScan("cParamScan", "", 1000, 1414);

  // Plot all graphs
  TH2F* hDummy[5];  // Dummy 2D histograms to set axes range
  cParamScan.Print(openPsStr);
  cParamScan.Clear();
  cParamScan.Divide(2, 3);
  for (int iPar = 0; iPar < 5; iPar++) {
    float xMin = grSeedEff[iPar].GetXaxis()->GetXmin();
    float xMax = grSeedEff[iPar].GetXaxis()->GetXmax();
    cParamScan.cd(iPar + 1);
    hDummy[iPar] = new TH2F("hDummy", "", 1, xMin, xMax, 1, 1.e-6, 1.e-3);
    hDummy[iPar]->SetXTitle("#font[52]{N}_{bins}");
    hDummy[iPar]->SetYTitle("Seed eff.");
    gPad->SetLogx();
    gPad->SetLogy();
    hDummy[iPar]->Draw();
    grSeedEff[iPar].SetMarkerStyle(20);
    grSeedEff[iPar].Draw("P");
  }
  cParamScan.Update();
  cParamScan.Print(psFileName);
  cParamScan.Clear();
  cParamScan.Divide(2, 3);
  for (int iPar = 0; iPar < 5; iPar++) {
    cParamScan.cd(iPar + 1);
    hDummy[iPar]->GetYaxis()->SetLimits(-0.1, 1.1);
    hDummy[iPar]->SetYTitle("Track eff.");
    gPad->SetLogx();
    hDummy[iPar]->Draw();
    grTrackEff[iPar].SetMarkerStyle(20);
    grTrackEff[iPar].Draw("P");
  }
  cParamScan.Update();
  cParamScan.Print(psFileName);

  cParamScan.Print(closePsStr);

  return 0;

}
