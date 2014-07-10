#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <TString.h>
#include <TFile.h>
#include <TH2F.h>

using namespace std;

int main(int argc, char* argv[])
{
  // Check arguments
  if (argc != 3) {
    cout << "Usage: writePairDistCode <input ROOT file> <hit fraction threshold>" << endl;
    return 1;
  }

  // Input file
  char* inFileName = argv[1];
  TFile* inFile = new TFile(inFileName);
  if (inFile == 0) {
    cout << "Error opening input file " << inFileName << ". Abort." << endl;
    return 1;
  }
  
  // Hit fraction threshold
  double fracThr;
  TString fracThrStr(argv[2]);
  if (!fracThrStr.IsFloat() || fracThrStr.Atof() <= 0. || fracThrStr.Atof() > 1.) {
    cout << "Second argument must be a number between 0 and 1. Abort." << endl;
    return 1;
  } else {
   fracThr = fracThrStr.Atof();
  }

  // Loop on layer histograms and extract max distance parameters
  ofstream codeFile("pairDistCode.cc");
  TH2F* hDistNext[27];
  for (int iLyr = 0; iLyr < 27; iLyr++) {
    TString hName = TString("hDistNext_") + (unsigned long)iLyr;
    hDistNext[iLyr] = (TH2F*)(inFile->Get(hName));
    if (hDistNext[iLyr] == 0) {
      cout << "Histogram " << hName << " not found!" << endl;
      continue;
    }
    for (int iNext = iLyr + 1; iNext < 27; iNext++) {
      double integral = 0.;
      for (int iBin = hDistNext[iLyr]->GetNbinsX(); iBin > 0; iBin--) {
        integral += hDistNext[iLyr]->GetBinContent(iBin, iNext + 1);
        if (integral > fracThr) {
          double maxDist = hDistNext[iLyr]->GetXaxis()->GetBinUpEdge(iBin);
          codeFile << "  maxPairDist_[(make_pair(" << iLyr << ", " << iNext << "))] = " << maxDist << ";" << endl;
          break;
        }
      }
    }
  }

}
