#include <iostream>
#include <dirent.h>
#include <stdlib.h>
#include <math.h>

#include <TGraph.h>
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TPad.h>
#include <TGaxis.h>

using namespace std;

int main(int argc, char* argv[])
{
  // Check arguments
  if (argc != 2) {
    cout << "Usage: plotParamScan <input directory name>" << endl;
    return 1;
  }
  char* dirName = argv[1];

  // Scan parameters (customized by case)
  const unsigned int nBinsDoca[3] = {4, 8, 16};
  const unsigned int nBinsSqrtK[3] = {20, 40, 80};
  const unsigned int nBinsPhi[3] = {64, 128, 256};
  const unsigned int nBinsZ0[3] = {3, 6, 12};
  const unsigned int nBinsEta[3] = {25, 50, 100};
  const unsigned int* nBinsAll[5] = {nBinsDoca, nBinsSqrtK, nBinsPhi, nBinsZ0, nBinsEta};
  const float phiBinOverlap[3] = {0., 0.25, 0.5};
  const float etaBinOverlap[3] = {0., 0.25, 0.5};
  const float* binOverlapAll[5] = {phiBinOverlap, etaBinOverlap};

  // Graphical parameters
  const int markerType[3] = {20, 21, 22};
  const int markerColor[3] = {kRed, kGreen, kBlue};
  const int markerBrightness[3] = {2, 0, -9};
  const float markerShift[3] = {-0.02, 0., 0.02};

  // Tree variables
  unsigned int nSeeds, nGoodSeeds, nTriedSeeds, nAssSeeds;
  unsigned int nTracks, nTracksAlgo, nAccTracks, nAssTracks, nFoundTracks;
  TBranch *bNSeeds, *bNGoodSeeds, *bNTriedSeeds, *bNAssSeeds;
  TBranch *bNTracks, *bNTracksAlgo, *bNAccTracks, *bNAssTracks, *bNFoundTracks;;

  // Graphs
  TGraph grSeedEff[5][9];
  TGraph grTrackEff[5][9];
  TGraph2D grEffVsFake[9];
  
  // Open the directory and look for ROOT files
  DIR* dir = opendir(dirName);
  if (dir == NULL) {
    cout << "Error opening directory " << dirName << ". Abort." << endl;
    return 1;
  }
  unsigned int nBins[5];
  float overlap[2];
  int step;
  while (struct dirent*  dirEntry = readdir(dir)) {
    string fileName = dirEntry->d_name;
    if (fileName.length() > 5 && fileName.compare(fileName.length() - 5, 5, string(".root")) == 0) {
      string filePath = string(dirName) + "/" + fileName;
      //      cout << filePath << endl;
      // Get run parameters
      for (int iPar = 0; iPar < 5; iPar++) {
        step = atoi((fileName.substr(fileName.length() - 12 + iPar, 1)).c_str());
        if (step == 0) {
          cout << "Incorrect file name " << fileName << ". Abort." << endl;
          return 2;
        }
        nBins[iPar] = nBinsAll[iPar][step - 1];
      }
      unsigned int iOvlp = 0;
      for (int iOLPar = 0; iOLPar < 2; iOLPar++) {
        step = atoi((fileName.substr(fileName.length() - 7 + iOLPar, 1)).c_str());
        if (step == 0) {
          cout << "Incorrect file name " << fileName << ". Abort." << endl;
          return 2;
        }
        iOvlp += pow(3, iOLPar)*(step - 1);
        overlap[iOLPar] = binOverlapAll[iOLPar][step - 1];
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
      //      tHT->Scan("nTriedSeeds:nAssSeeds:nAssTracks:nFoundTracks");

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
        grSeedEff[iPar][iOvlp].SetPoint(grSeedEff[iPar][iOvlp].GetN(), nBins[iPar], (float)nAssSeeds/nTriedSeeds);
        grTrackEff[iPar][iOvlp].SetPoint(grSeedEff[iPar][iOvlp].GetN(), nBins[iPar], (float)nFoundTracks/nAssTracks);
        grEffVsFake[iOvlp].SetPoint(grEffVsFake[iOvlp].GetN(), (float)nAssSeeds/nTriedSeeds, (float)nFoundTracks/nAssTracks, nBins[0]*nBins[1]*nBins[2]*nBins[3]*nBins[4]);
      }
    }
  }

  // Set up graphics
  TString psFileName = "paramBinScan.ps";
  TString openPsStr = psFileName + "[";
  TString closePsStr = psFileName + "]";
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0001);
  gStyle->SetOptStat(0);
  TCanvas cParamScan("cParamScan", "", 1000, 1414);

  // Plot all graphs
  TH2F* hDummy[9];  // Dummy 2D histograms to set axes range
  cParamScan.Print(openPsStr);
  for (int iPar = 0; iPar < 5; iPar++) {
    cParamScan.Clear();
    cParamScan.Divide(3, 3);
    float xMin = nBinsAll[iPar][0] - 0.1*(nBinsAll[iPar][2] - nBinsAll[iPar][0]);
    float xMax = nBinsAll[iPar][2] + 0.1*(nBinsAll[iPar][2] - nBinsAll[iPar][0]);
    for (int iOvlp = 0; iOvlp < 9; iOvlp++) {
      cParamScan.cd(iOvlp + 1);
      hDummy[iOvlp] = new TH2F("hDummy", "", 1, xMin, xMax, 1, 1.e-6, 1.e-3);
      hDummy[iOvlp]->SetXTitle("#font[52]{N}_{bins}");
      hDummy[iOvlp]->SetYTitle("Seed eff.");
      gPad->SetLogy();
      hDummy[iOvlp]->Draw();
      grTrackEff[iPar][iOvlp].SetMarkerStyle(20);
      grSeedEff[iPar][iOvlp].SetMarkerStyle(20);
      grSeedEff[iPar][iOvlp].Draw("P");
    }
    cParamScan.Update();
    cParamScan.Print(psFileName);
    for (int iOvlp = 0; iOvlp < 9; iOvlp++)
      delete hDummy[iOvlp];
    cParamScan.Clear();
    cParamScan.Divide(3, 3);
    for (int iOvlp = 0; iOvlp < 9; iOvlp++) {
      cParamScan.cd(iOvlp + 1);
      hDummy[iOvlp] = new TH2F("hDummy", "", 1, xMin, xMax, 1, -0.1, 1.1);
      hDummy[iOvlp]->SetXTitle("#font[52]{N}_{bins}");
      hDummy[iOvlp]->SetYTitle("Track eff.");
      //      gPad->SetLogy();
      hDummy[iOvlp]->Draw();
      grTrackEff[iPar][iOvlp].SetMarkerStyle(20);
      grTrackEff[iPar][iOvlp].Draw("P");
    }
    cParamScan.Update();
    cParamScan.Print(psFileName);
    for (int iOvlp = 0; iOvlp < 9; iOvlp++)
      delete hDummy[iOvlp];
  }
  TFile fTest("fTest.root", "RECREATE");
  cParamScan.Clear();
  cParamScan.Divide(3, 3);
  for (int iOvlp = 0; iOvlp < 9; iOvlp++) {
    cParamScan.cd(iOvlp + 1);
    gPad->SetMargin(0.1, 0.15, 0.235, 0.235);
    hDummy[iOvlp] = new TH2F("hDummy", "", 1, 1.e-6, 1.e-3, 1, 0., 1.);
    hDummy[iOvlp]->SetXTitle("Seed eff.");
    hDummy[iOvlp]->SetYTitle("TrackEff eff.");
    gPad->SetLogx();
    //    hDummy[iOvlp]->Draw();
    gPad->SetLogz();
    gPad->SetPhi(0.);
    gPad->SetTheta(90.);
    grEffVsFake[iOvlp].SetMarkerStyle(20);
    // grEffVsFake[iOvlp].GetYaxis()->SetLabelOffset(999);
    // grEffVsFake[iOvlp].GetYaxis()->SetTickLength(0);
    grEffVsFake[iOvlp].Draw("PCOLZ");
    // TGaxis newYAxis(gPad->GetUxmin(), gPad->GetUymin(), gPad->GetUxmin(), gPad->GetUymax(), grEffVsFake[iOvlp].GetYaxis()->GetXmin(), grEffVsFake[iOvlp].GetYaxis()->GetXmax(), 510, "G");
    // newYAxis.SetLabelFont(grEffVsFake[iOvlp].GetXaxis()->GetLabelFont());
    // newYAxis.SetLabelSize(grEffVsFake[iOvlp].GetXaxis()->GetLabelSize());
    // newYAxis.Draw();
    grEffVsFake[iOvlp].Write();
  }
  cParamScan.Update();
  cParamScan.Print(psFileName);
  for (int iOvlp = 0; iOvlp < 9; iOvlp++)
    delete hDummy[iOvlp];

  fTest.Close();

  cParamScan.Print(closePsStr);

  return 0;

}
