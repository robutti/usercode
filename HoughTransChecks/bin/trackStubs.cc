#include <vector>
#include <iostream>
#include <stdlib.h>
#include <math.h>

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH2I.h>
#include <TH2F.h>
#include <TPad.h>
#include <TString.h>
#include <TPDF.h>
#include <TPaletteAxis.h>

using namespace std;

int main(int argc, char* argv[])
{
  // Parameters
  int lyrMapOffset[7] = {0, 0, 3, 5, 9, 12, 18};

  // Check arguments
  if (argc != 2) {
    cout << "Usage: trackStubs <input file name>" << endl;
    return 1;
  }
  char* inFileName = argv[1];

  // Open input file and get objects
  TFile* inFile = new TFile(inFileName);
  if (inFile == 0) {
    cout << "Error opening input file " << inFileName << ". Abort." << endl;
    return 1;
  }

  // Histograms
  TH2I* hNextLayer = new TH2I("hNextLayer", "Next layer in track", 27, 0., 27., 27, 0., 27.);
  TH2F* hDistNext[27];
  for (int iLyr = 0; iLyr < 27; iLyr++) {
    TString hName = TString("hDistNext_") + (unsigned long)iLyr;
    TString hTitle = TString("Distance of next hit for layer ") + (unsigned long)iLyr;
    hDistNext[iLyr] = new TH2F(hName, hTitle, 75, 0., 150., 27, 0., 27.);
  }

  // Tree variables
  vector<vector<double> >* subdet = 0;
  vector<vector<double> >* layer = 0;
  vector<vector<double> >* xHit = 0;
  vector<vector<double> >* yHit = 0;
  vector<vector<double> >* zHit = 0;
  unsigned int nTracks;

   // Branches
  TBranch* b_subdet;
  TBranch* b_layer;
  TBranch* b_xHit;
  TBranch* b_yHit;
  TBranch* b_zHit;
  TBranch* b_nTracks;

  // Get tree and set branches
  TTree* tTracks = (TTree*)(inFile->Get("trackTree"));
  tTracks->SetBranchAddress("subdet", &subdet, &b_subdet);
  tTracks->SetBranchAddress("layer", &layer, &b_layer);
  tTracks->SetBranchAddress("xHit", &xHit, &b_xHit);
  tTracks->SetBranchAddress("yHit", &yHit, &b_yHit);
  tTracks->SetBranchAddress("zHit", &zHit, &b_zHit);
  tTracks->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
  
  // Loop on tree
  for (int iEv = 0; iEv < tTracks->GetEntries(); iEv++) {
    long tEntry = tTracks->LoadTree(iEv);
    tTracks->GetEntry(tEntry);
    cout << "Event n. " << iEv << ": " << nTracks << " tracks found." << endl;
    for (int iTrack = 0; iTrack < nTracks; iTrack++) {
      vector<double> trkSubdet = subdet->at(iTrack); 
      vector<double> trkLayer = layer->at(iTrack); 
      vector<double> trkXHit = xHit->at(iTrack); 
      vector<double> trkYHit = yHit->at(iTrack); 
      vector<double> trkZHit = zHit->at(iTrack); 
      int nHits = trkSubdet.size();
      cout << "  Track n. " << iTrack << ": " << nHits << " hits." << endl;
      for (int iHit = 0; iHit < nHits; iHit++) {
        float minDist = 1.e6;
        int nextLayer = -1;
        for (int iNextHit = 0; iNextHit < nHits; iNextHit++) {
          if (iNextHit != iHit &&
              (trkSubdet[iNextHit] > trkSubdet[iHit] || 
               ((trkSubdet[iNextHit] == trkSubdet[iHit] && trkLayer[iNextHit] > trkLayer[iHit])))) {
            float dist = sqrt((trkXHit[iNextHit] - trkXHit[iHit])*(trkXHit[iNextHit] - trkXHit[iHit]) +
                              (trkYHit[iNextHit] - trkYHit[iHit])*(trkYHit[iNextHit] - trkYHit[iHit]) +
                              (trkZHit[iNextHit] - trkZHit[iHit])*(trkZHit[iNextHit] - trkZHit[iHit]));
            if (dist < minDist) {
              minDist = dist;
              nextLayer = lyrMapOffset[int(trkSubdet[iNextHit])] + trkLayer[iNextHit] - 1;
            }
          }
        }
        if (nextLayer >= 0) {
          int layer = lyrMapOffset[int(trkSubdet[iHit])] + trkLayer[iHit] - 1;
          hNextLayer->Fill(layer, nextLayer);
          hDistNext[layer]->Fill(minDist, nextLayer);
        }
      }
    }
  }

  // Plot histograms on PDF file
  // Set up graphics
  TString pdfFileName = "trackStubs.pdf";
  TString openPdfStr = pdfFileName + "[";
  TString closePdfStr = pdfFileName + "]";
  gROOT->SetStyle("plain");
  gStyle->SetPalette(1);
  gStyle->SetOptFit(0001);
  gStyle->SetOptStat(0);
  TCanvas cNextHit("cNextHit", "Next hit in track", 1000, 1414);
  cNextHit.Print(openPdfStr);
  // Next layer histogram
  cNextHit.Divide(1, 2, 0., 0.);
  cNextHit.cd(1);
  gPad->SetMargin(0.2172, 0.2172, 0.1, 0.1);
  hNextLayer->Draw("COLZ");
  gPad->Update();
  TPaletteAxis* palette = (TPaletteAxis*)(hNextLayer->GetListOfFunctions()->FindObject("palette"));
  palette->SetX1NDC(0.8);
  palette->SetX2NDC(0.82);
  palette->SetY1NDC(0.1);
  palette->SetY2NDC(0.9);
  hNextLayer->Draw("COLZ");
  cNextHit.Update();
  cNextHit.cd(2);
  gPad->SetMargin(0.2172, 0.2172, 0.1, 0.1);
  gPad->SetLogz();
  hNextLayer->Draw("COLZ");
  gPad->Update();
  palette = (TPaletteAxis*)(hNextLayer->GetListOfFunctions()->FindObject("palette"));
  palette->SetX1NDC(0.8);
  palette->SetX2NDC(0.82);
  palette->SetY1NDC(0.1);
  palette->SetY2NDC(0.9);
  hNextLayer->Draw("COLZ");
  cNextHit.Update();
  cNextHit.Print(pdfFileName);
  // Next hit distance histograms
  cNextHit.Clear();
  cNextHit.Divide(1, 4);
  int iPad = 0;
  for (int iLyr = 0; iLyr < 27; iLyr++) {
    cNextHit.cd(++iPad);
    double norm = 1./(hDistNext[iLyr]->Integral());
    hDistNext[iLyr]->Scale(norm);
    gPad->SetLogz();
    hDistNext[iLyr]->Draw("COLZ");
    // Print page if complete      
    if (iPad == 4) {
      cNextHit.Update();
      cNextHit.Print(pdfFileName);
      cNextHit.Clear();
      cNextHit.Divide(1, 4);
      iPad = 0;
    }
  }
  cNextHit.Update();
  cNextHit.Print(pdfFileName);
  // Close PDF file
  cNextHit.Print(closePdfStr);

  // Write histograms to file
  TFile outFile("trackStubs.root", "RECREATE");
  hNextLayer->Write();
  for (int iLyr = 0; iLyr < 27; iLyr++)
    hDistNext[iLyr]->Write();

  return 0;

}
