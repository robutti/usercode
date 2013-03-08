// -*- C++ -*-
//
// Package:    HoughTransChecks
// Class:      TH5
// 
/**\class TH5 TH5.cc HoughTest/HoughTransChecks/src/TH5.cc

 Description: 5-D histogram classes

 Implementation:
     Only char specialization TH5C (one byte per bin) implemented so far.
*/
//
// Original Author:  Enrico Robutti (adapted from ROOT TH3)
//         Created:  Fri, Feb 20, 2013
// $Id: TH5
//
//

#include "TROOT.h"
#include "TClass.h"
#include "THashList.h"
#include "../interface/TH5.h"
//#include "ERobutti/HoughTransChecks/interface/TH5.h"
#include "TProfile2D.h"
#include "TH2.h"
#include "TF1.h"
#include "TVirtualPad.h"
#include "TVirtualHistPainter.h"
#include "THLimitsFinder.h"
#include "TRandom.h"
#include "TError.h"
#include "TMath.h"
#include "TObjString.h"
#include "TStyle.h"

ClassImp(TH5)
//______________________________________________________________________________
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
//*-*
//*-*  The 5-D histogram classes derived from the 1-D histogram classes.
//*-*
//
//  TH5C a 5-D histogram with one byte per cell (char)
//
//*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

//______________________________________________________________________________
TH5::TH5()
{
   // Default constructor.
  fDimension = 5;
  fUaxis.SetName("uaxis");
  fVaxis.SetName("vaxis");
  fUaxis.SetParent(this);
  fVaxis.SetParent(this);
  fTsumwy = fTsumwy2 = fTsumwxy = 0;
  fTsumwz = fTsumwz2 = fTsumwxz = fTsumwyz = 0;
  fTsumwu = fTsumwu2 = fTsumwxu = fTsumwyu = fTsumwzu = 0;
  fTsumwv = fTsumwv2 = fTsumwxv = fTsumwyv = fTsumwzv = fTsumwuv = 0;
}

//______________________________________________________________________________
TH5::TH5(const char *name, const char *title,
	 Int_t nbinsx, Double_t xlow, Double_t xup,
	 Int_t nbinsy, Double_t ylow, Double_t yup,
	 Int_t nbinsz, Double_t zlow, Double_t zup,
	 Int_t nbinsu, Double_t ulow, Double_t uup,
	 Int_t nbinsv, Double_t vlow, Double_t vup)
:TH1(name,title,nbinsx,xlow,xup)
{
//*-*-*-*-*-*-*-*-*Normal constructor for fixed bin size 5-D histograms*-*-*-*-*
//*-*              ==================================================

  fDimension = 5;
  fUaxis.SetName("uaxis");
  fVaxis.SetName("vaxis");
  fUaxis.SetParent(this);
  fVaxis.SetParent(this);
  if (nbinsy <= 0) nbinsy = 1;
  if (nbinsz <= 0) nbinsz = 1;
  if (nbinsu <= 0) nbinsu = 1;
  if (nbinsv <= 0) nbinsv = 1;
  fYaxis.Set(nbinsy, ylow, yup);
  fZaxis.Set(nbinsz, zlow, zup);
  fUaxis.Set(nbinsu, ulow, uup);
  fVaxis.Set(nbinsv, vlow, vup);
  fNcells = (nbinsx + 2)*(nbinsy + 2)*(nbinsz + 2)*(nbinsu + 2)*(nbinsv + 2);
  fTsumwy = fTsumwy2 = fTsumwxy = 0;
  fTsumwz = fTsumwz2 = fTsumwxz = fTsumwyz = 0;
  fTsumwu = fTsumwu2 = fTsumwxu = fTsumwyu = fTsumwzu = 0;
  fTsumwv = fTsumwv2 = fTsumwxv = fTsumwyv = fTsumwzv = fTsumwuv = 0;
  UseCurrentStyle();
}

//______________________________________________________________________________
TH5::TH5(const char *name, const char *title,
	 Int_t nbinsx, const Float_t *xbins,
	 Int_t nbinsy, const Float_t *ybins,
	 Int_t nbinsz, const Float_t *zbins,
	 Int_t nbinsu, const Float_t *ubins,
	 Int_t nbinsv, const Float_t *vbins)
:TH1(name,title,nbinsx,xbins)
{
//*-*-*-*-*-*-*-*Normal constructor for variable bin size 5-D histograms*-*-*-*
//*-*            =======================================================
  fDimension = 5;
  fUaxis.SetName("uaxis");
  fVaxis.SetName("vaxis");
  fUaxis.SetParent(this);
  fVaxis.SetParent(this);
  if (nbinsy <= 0) nbinsy = 1;
  if (nbinsz <= 0) nbinsz = 1;
  if (nbinsu <= 0) nbinsu = 1;
  if (nbinsv <= 0) nbinsv = 1;
  if (ybins)
    fYaxis.Set(nbinsy, ybins);
  else
    fYaxis.Set(nbinsy, 0, 1);
  if (zbins)
    fZaxis.Set(nbinsz, zbins);
  else
    fZaxis.Set(nbinsz, 0, 1);
  if (ubins)
    fUaxis.Set(nbinsu, ubins);
  else
    fUaxis.Set(nbinsu, 0, 1);
  if (vbins)
    fVaxis.Set(nbinsv, vbins);
  else
    fVaxis.Set(nbinsv, 0, 1);
  fNcells = (nbinsx + 2)*(nbinsy + 2)*(nbinsz + 2)*(nbinsu + 2)*(nbinsv + 2);
  fTsumwy = fTsumwy2 = fTsumwxy = 0;
  fTsumwz = fTsumwz2 = fTsumwxz = fTsumwyz = 0;
  fTsumwu = fTsumwu2 = fTsumwxu = fTsumwyu = fTsumwzu = 0;
  fTsumwv = fTsumwv2 = fTsumwxv = fTsumwyv = fTsumwzv = fTsumwuv = 0;
  UseCurrentStyle();
}

//______________________________________________________________________________
TH5::TH5(const char *name, const char *title,
	 Int_t nbinsx, const Double_t *xbins,
	 Int_t nbinsy, const Double_t *ybins,
	 Int_t nbinsz, const Double_t *zbins,
	 Int_t nbinsu, const Double_t *ubins,
	 Int_t nbinsv, const Double_t *vbins)
:TH1(name, title, nbinsx, xbins)
{
//*-*-*-*-*-*-*-*Normal constructor for variable bin size 5-D histograms*-*-*-*
//*-*            =======================================================
  fDimension = 5;
  fUaxis.SetName("uaxis");
  fVaxis.SetName("vaxis");
  fUaxis.SetParent(this);
  fVaxis.SetParent(this);
  if (nbinsy <= 0) nbinsy = 1;
  if (nbinsz <= 0) nbinsz = 1;
  if (nbinsu <= 0) nbinsu = 1;
  if (nbinsv <= 0) nbinsv = 1;
  if (ybins)
    fYaxis.Set(nbinsy, ybins);
  else
    fYaxis.Set(nbinsy, 0, 1);
  if (zbins)
    fZaxis.Set(nbinsz, zbins);
  else
    fZaxis.Set(nbinsz, 0, 1);
  if (ubins)
    fUaxis.Set(nbinsu, ubins);
  else
    fUaxis.Set(nbinsu, 0, 1);
  if (vbins)
    fVaxis.Set(nbinsv, vbins);
  else
    fVaxis.Set(nbinsv, 0, 1);
  fNcells = (nbinsx + 2)*(nbinsy + 2)*(nbinsz + 2)*(nbinsu + 2)*(nbinsv + 2);
  fTsumwy = fTsumwy2 = fTsumwxy = 0;
  fTsumwz = fTsumwz2 = fTsumwxz = fTsumwyz = 0;
  fTsumwu = fTsumwu2 = fTsumwxu = fTsumwyu = fTsumwzu = 0;
  fTsumwv = fTsumwv2 = fTsumwxv = fTsumwyv = fTsumwzv = fTsumwuv = 0;
  UseCurrentStyle();
}

//______________________________________________________________________________
TH5::TH5(const TH5& h) : TH1()
{
  // Copy constructor.
  // The list of functions is not copied. (Use Clone if needed)
  
  ((TH5&)h).Copy(*this);
}

//______________________________________________________________________________
TH5::~TH5()
{
  // Destructor.
}

//______________________________________________________________________________
void TH5::Copy(TObject &obj) const
{
  // Copy.
  
  TH1::Copy(obj);
  ((TH5&)obj).fTsumwy = fTsumwy;
  ((TH5&)obj).fTsumwy2 = fTsumwy2;
  ((TH5&)obj).fTsumwxy = fTsumwxy;
  ((TH5&)obj).fTsumwz = fTsumwz;
  ((TH5&)obj).fTsumwz2 = fTsumwz2;
  ((TH5&)obj).fTsumwxz = fTsumwxz;
  ((TH5&)obj).fTsumwyz = fTsumwyz;
  ((TH5&)obj).fTsumwu = fTsumwu;
  ((TH5&)obj).fTsumwu2 = fTsumwu2;
  ((TH5&)obj).fTsumwxu = fTsumwxu;
  ((TH5&)obj).fTsumwyu = fTsumwyu;
  ((TH5&)obj).fTsumwzu = fTsumwzu;
  ((TH5&)obj).fTsumwv = fTsumwv;
  ((TH5&)obj).fTsumwv2 = fTsumwv2;
  ((TH5&)obj).fTsumwxv = fTsumwxv;
  ((TH5&)obj).fTsumwyv = fTsumwyv;
  ((TH5&)obj).fTsumwzv = fTsumwzv;
  ((TH5&)obj).fTsumwuv = fTsumwuv;
}

//______________________________________________________________________________
Int_t TH5::BufferEmpty(Int_t action)
{
  // Fill histogram with all entries in the buffer.
  // action = -1 histogram is reset and refilled from the buffer (called by THistPainter::Paint)
  // action =  0 histogram is filled from the buffer
  // action =  1 histogram is filled and buffer is deleted
  //             The buffer is automatically deleted when the number of entries
  //             in the buffer is greater than the number of entries in the histogram
  
  // do we need to compute the bin size?
  if (!fBuffer) return 0;
  Int_t nbentries = (Int_t)fBuffer[0];
  if (!nbentries) return 0;
  Double_t* buffer = fBuffer;
  if (nbentries < 0) {
    if (action == 0) return 0;
    nbentries = -nbentries;
    fBuffer = 0;
    Reset("ICES");  
    fBuffer = buffer;
  }
  if (TestBit(kCanRebin) || fXaxis.GetXmax() <= fXaxis.GetXmin() ||
      fYaxis.GetXmax() <= fYaxis.GetXmin() ||
      fZaxis.GetXmax() <= fZaxis.GetXmin() ||
      fUaxis.GetXmax() <= fUaxis.GetXmin() ||
      fVaxis.GetXmax() <= fVaxis.GetXmin()) {
    //find min, max of entries in buffer
    Double_t xmin = fBuffer[2];
    Double_t xmax = xmin;
    Double_t ymin = fBuffer[3];
    Double_t ymax = ymin;
    Double_t zmin = fBuffer[4];
    Double_t zmax = zmin;
    Double_t umin = fBuffer[5];
    Double_t umax = umin;
    Double_t vmin = fBuffer[6];
    Double_t vmax = vmin;
    for (Int_t i = 1; i < nbentries; i++) {
      Double_t x = fBuffer[6*i + 2];
      if (x < xmin) xmin = x;
      if (x > xmax) xmax = x;
      Double_t y = fBuffer[6*i + 3];
      if (y < ymin) ymin = y;
      if (y > ymax) ymax = y;
      Double_t z = fBuffer[6*i + 4];
      if (z < zmin) zmin = z;
      if (z > zmax) zmax = z;
      Double_t u = fBuffer[6*i + 5];
      if (u < umin) umin = u;
      if (u > umax) umax = u;
      Double_t v = fBuffer[6*i + 6];
      if (v < vmin) vmin = v;
      if (v > vmax) vmax = v;
    }
    if (fXaxis.GetXmax() <= fXaxis.GetXmin() ||
	fYaxis.GetXmax() <= fYaxis.GetXmin() ||
	fZaxis.GetXmax() <= fZaxis.GetXmin() ||
	fUaxis.GetXmax() <= fUaxis.GetXmin() ||
	fVaxis.GetXmax() <= fVaxis.GetXmin()) {
      //    THLimitsFinder::GetLimitsFinder()->FindGoodLimits(this,xmin,xmax,ymin,ymax,zmin,zmax);
    } else {
      fBuffer = 0;
      Int_t keep = fBufferSize; fBufferSize = 0;
      if (xmin < fXaxis.GetXmin()) RebinAxis(xmin, &fXaxis);
      if (xmax >= fXaxis.GetXmax()) RebinAxis(xmax, &fXaxis);
      if (ymin < fYaxis.GetXmin()) RebinAxis(ymin, &fYaxis);
      if (ymax >= fYaxis.GetXmax()) RebinAxis(ymax, &fYaxis);
      if (zmin < fZaxis.GetXmin()) RebinAxis(zmin, &fZaxis);
      if (zmax >= fZaxis.GetXmax()) RebinAxis(zmax, &fZaxis);
      if (umin < fUaxis.GetXmin()) RebinAxis(umin, &fUaxis);
      if (umax >= fUaxis.GetXmax()) RebinAxis(umax, &fUaxis);
      if (vmin < fVaxis.GetXmin()) RebinAxis(vmin, &fVaxis);
      if (vmax >= fVaxis.GetXmax()) RebinAxis(vmax, &fVaxis);
      fBuffer = buffer;
      fBufferSize = keep;
    }
  }
  fBuffer = 0;
  
  for (Int_t i=0;i < nbentries; i++) {
    Fill(buffer[6*i+2], buffer[6*i + 3], buffer[6*i + 4], buffer[6*i + 5], buffer[6*i + 6], buffer[6*i + 1]);
  }
  fBuffer = buffer;

  if (action > 0) { delete [] fBuffer; fBuffer = 0; fBufferSize = 0;}
  else {
    if (nbentries == (Int_t)fEntries)
      fBuffer[0] = -nbentries;
    else
      fBuffer[0] = 0;
  }
  return nbentries;
}

//______________________________________________________________________________
Int_t TH5::BufferFill(Double_t x, Double_t y, Double_t z, Double_t u, Double_t v, Double_t w)
{
  // accumulate arguments in buffer. When buffer is full, empty the buffer
  // fBuffer[0] = number of entries in buffer
  // fBuffer[1] = w of first entry
  // fBuffer[2] = x of first entry
  // fBuffer[3] = y of first entry
  // fBuffer[4] = z of first entry
  // fBuffer[5] = u of first entry
  // fBuffer[6] = v of first entry

  if (!fBuffer) return -3;
  Int_t nbentries = (Int_t)fBuffer[0];
  if (nbentries < 0) {
    nbentries  = -nbentries;
    fBuffer[0] =  nbentries;
    if (fEntries > 0) {
      Double_t *buffer = fBuffer; fBuffer=0;
      Reset("ICES");  
      fBuffer = buffer;
    }
  }
  if (6*nbentries + 6 >= fBufferSize) {
    BufferEmpty(1);
    return Fill(x, y ,z, u, v, w);
  }
  fBuffer[6*nbentries + 1] = w;
  fBuffer[6*nbentries + 2] = x;
  fBuffer[6*nbentries + 3] = y;
  fBuffer[6*nbentries + 4] = z;
  fBuffer[6*nbentries + 5] = u;
  fBuffer[6*nbentries + 6] = v;
  fBuffer[0] += 1;
  return -3;
}

//______________________________________________________________________________
Int_t TH5::Fill(Double_t )
{
   // Invalid Fill method
   Error("Fill", "Invalid signature - do nothing"); 
   return -1;
}
//______________________________________________________________________________
Int_t TH5::Fill(Double_t x, Double_t y, Double_t z, Double_t u, Double_t v)
{
   //*-*-*-*-*-*-*-*-*-*-*Increment cell defined by x, y, z, u, v by 1 *-*-*-*-*
   //*-*                  ====================================
   //*-*
   //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  if (fBuffer) return BufferFill(x, y, z, u, v, 1);
  
  Int_t binx, biny, binz, binu, binv, bin;
  fEntries++;
  binx = fXaxis.FindBin(x);
  biny = fYaxis.FindBin(y);
  binz = fZaxis.FindBin(z);
  binu = fUaxis.FindBin(u);
  binv = fVaxis.FindBin(v);
  if (binx < 0 || biny < 0 || binz < 0 || binu < 0 || binv < 0) return -1;
  bin = GetBin(binx, biny, binz, binu, binv);
  AddBinContent(bin);
  if (fSumw2.fN) ++fSumw2.fArray[bin];
  if (binx == 0 || binx > fXaxis.GetNbins())
    if (!fgStatOverflows) return -1;
  if (biny == 0 || biny > fYaxis.GetNbins())
    if (!fgStatOverflows) return -1;
  if (binz == 0 || binz > fZaxis.GetNbins())
    if (!fgStatOverflows) return -1;
  if (binu == 0 || binu > fUaxis.GetNbins())
    if (!fgStatOverflows) return -1;
  if (binv == 0 || binv > fVaxis.GetNbins())
    if (!fgStatOverflows) return -1;

  ++fTsumw;
  ++fTsumw2;
  fTsumwx  += x;
  fTsumwx2 += x*x;
  fTsumwy  += y;
  fTsumwy2 += y*y;
  fTsumwxy += x*y;
  fTsumwz  += z;
  fTsumwz2 += z*z;
  fTsumwxz += x*z;
  fTsumwyz += y*z;
  fTsumwu  += u;
  fTsumwu2 += u*u;
  fTsumwxu += x*u;
  fTsumwyu += y*u;
  fTsumwzu += z*u;
  fTsumwv  += v;
  fTsumwv2 += v*v;
  fTsumwxv += x*v;
  fTsumwyv += y*v;
  fTsumwzv += z*v;
  fTsumwuv += u*v;
  return bin;
}

//______________________________________________________________________________
Int_t TH5::Fill(Double_t x, Double_t y, Double_t z, Double_t u, Double_t v, Double_t w)
{
  //*-*-*-*-*-*-*-*-*-*-*Increment cell defined by x, y, z, u, v by a weight w*-*-*-*-*
  //*-*                     //*-*
  //*-* If the storage of the sum of squares of weights has been triggered,
  //*-* via the function Sumw2, then the sum of the squares of weights is incremented
  //*-* by w^2 in the cell corresponding to x, y, z, u, v.
  //*-*
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  if (fBuffer) return BufferFill(x, y, z, u, v, w);

  Int_t binx, biny, binz, binu, binv, bin;
  fEntries++;
  binx = fXaxis.FindBin(x);
  biny = fYaxis.FindBin(y);
  binz = fZaxis.FindBin(z);
  binu = fUaxis.FindBin(u);
  binv = fVaxis.FindBin(v);
  if (binx < 0 || biny < 0 || binz < 0 || binu < 0 || binv < 0) return -1;
  bin = GetBin(binx, biny, binz, binu, binv);
  AddBinContent(bin,w);
  if (fSumw2.fN) fSumw2.fArray[bin] += w*w;
  if (binx == 0 || binx > fXaxis.GetNbins()) {
    if (!fgStatOverflows) return -1;
  }
  if (biny == 0 || biny > fYaxis.GetNbins()) {
    if (!fgStatOverflows) return -1;
  }
  if (binz == 0 || binz > fZaxis.GetNbins()) {
    if (!fgStatOverflows) return -1;
  }
  if (binu == 0 || binu > fUaxis.GetNbins()) {
    if (!fgStatOverflows) return -1;
  }
  if (binv == 0 || binv > fVaxis.GetNbins()) {
    if (!fgStatOverflows) return -1;
  }
   fTsumw   += w;
   fTsumw2  += w*w;
   fTsumwx  += w*x;
   fTsumwx2 += w*x*x;
   fTsumwy  += w*y;
   fTsumwy2 += w*y*y;
   fTsumwxy += w*x*y;
   fTsumwz  += w*z;
   fTsumwz2 += w*z*z;
   fTsumwxz += w*x*z;
   fTsumwyz += w*y*z;
   fTsumwu  += w*u;
   fTsumwu2 += w*u*u;
   fTsumwxu += w*x*u;
   fTsumwyu += w*y*u;
   fTsumwzu += w*z*u;
   fTsumwu  += w*v;
   fTsumwv2 += w*v*v;
   fTsumwxv += w*x*v;
   fTsumwyv += w*y*v;
   fTsumwzv += w*z*v;
   fTsumwuv += w*u*v;
   return bin;
}

//______________________________________________________________________________
// void TH5::FillRandom(const char *fname, Int_t ntimes)
// {
//    //*-*-*-*-*-*-*Fill histogram following distribution in function fname*-*-*-*
//    //*-*          =======================================================
//    //*-*
//    //*-*   The distribution contained in the function fname (TF1) is integrated
//    //*-*   over the channel contents.
//    //*-*   It is normalized to 1.
//    //*-*   Getting one random number implies:
//    //*-*     - Generating a random number between 0 and 1 (say r1)
//    //*-*     - Look in which bin in the normalized integral r1 corresponds to
//    //*-*     - Fill histogram channel
//    //*-*   ntimes random numbers are generated
//    //*-*
//    //*-*  One can also call TF1::GetRandom to get a random variate from a function.
//    //*-*
//    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*

//    Int_t bin, binx, biny, binz, ibin, loop;
//    Double_t r1, x, y,z, xv[3];
//    //*-*- Search for fname in the list of ROOT defined functions
//    TF1 *f1 = (TF1*)gROOT->GetFunction(fname);
//    if (!f1) { Error("FillRandom", "Unknown function: %s",fname); return; }

//    //*-*- Allocate temporary space to store the integral and compute integral
//    Int_t nbinsx = GetNbinsX();
//    Int_t nbinsy = GetNbinsY();
//    Int_t nbinsz = GetNbinsZ();
//    Int_t nxy    = nbinsx*nbinsy;
//    Int_t nbins  = nxy*nbinsz;

//    Double_t *integral = new Double_t[nbins+1];
//    ibin = 0;
//    integral[ibin] = 0;
//    for (binz=1;binz<=nbinsz;binz++) {
//       xv[2] = fZaxis.GetBinCenter(binz);
//       for (biny=1;biny<=nbinsy;biny++) {
//          xv[1] = fYaxis.GetBinCenter(biny);
//          for (binx=1;binx<=nbinsx;binx++) {
//             xv[0] = fXaxis.GetBinCenter(binx);
//             ibin++;
//             integral[ibin] = integral[ibin-1] + f1->Eval(xv[0],xv[1],xv[2]);
//          }
//       }
//    }

//    //*-*- Normalize integral to 1
//    if (integral[nbins] == 0 ) {
//       delete [] integral;
//       Error("FillRandom", "Integral = zero"); return;
//    }
//    for (bin=1;bin<=nbins;bin++)  integral[bin] /= integral[nbins];

//    //*-*--------------Start main loop ntimes
//    if (fDimension < 2) nbinsy = -1;
//    if (fDimension < 3) nbinsz = -1;
//    for (loop=0;loop<ntimes;loop++) {
//       r1 = gRandom->Rndm(loop);
//       ibin = TMath::BinarySearch(nbins,&integral[0],r1);
//       binz = ibin/nxy;
//       biny = (ibin - nxy*binz)/nbinsx;
//       binx = 1 + ibin - nbinsx*(biny + nbinsy*binz);
//       if (nbinsz) binz++;
//       if (nbinsy) biny++;
//       x    = fXaxis.GetBinCenter(binx);
//       y    = fYaxis.GetBinCenter(biny);
//       z    = fZaxis.GetBinCenter(binz);
//       Fill(x,y,z, 1.);
//    }
//    delete [] integral;
// }

//______________________________________________________________________________
// void TH5::FillRandom(TH1 *h, Int_t ntimes)
// {
//    //*-*-*-*-*-*-*Fill histogram following distribution in histogram h*-*-*-*
//    //*-*          ====================================================
//    //*-*
//    //*-*   The distribution contained in the histogram h (TH5) is integrated
//    //*-*   over the channel contents.
//    //*-*   It is normalized to 1.
//    //*-*   Getting one random number implies:
//    //*-*     - Generating a random number between 0 and 1 (say r1)
//    //*-*     - Look in which bin in the normalized integral r1 corresponds to
//    //*-*     - Fill histogram channel
//    //*-*   ntimes random numbers are generated
//    //*-*
//    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*

//    if (!h) { Error("FillRandom", "Null histogram"); return; }
//    if (fDimension != h->GetDimension()) {
//       Error("FillRandom", "Histograms with different dimensions"); return;
//    }

//    if (h->ComputeIntegral() == 0) return;

//    TH5 *h3 = (TH5*)h;
//    Int_t loop;
//    Double_t x,y,z;
//    for (loop=0;loop<ntimes;loop++) {
//       h3->GetRandom3(x,y,z);
//       Fill(x,y,z,1.);
//    }
// }

//______________________________________________________________________________
Double_t TH5::GetCorrelationFactor(Int_t axis1, Int_t axis2) const
{
  //*-*-*-*-*-*-*-*Return correlation factor between axis1 and axis2*-*-*-*-*
  //*-*            ====================================================
  if (axis1 < 1 || axis2 < 1 || axis1 > 5 || axis2 > 5) {
    Error("GetCorrelationFactor","Wrong parameters");
    return 0;
  }
  if (axis1 == axis2) return 1;
  Double_t rms1 = GetRMS(axis1);
  if (rms1 == 0) return 0;
  Double_t rms2 = GetRMS(axis2);
  if (rms2 == 0) return 0;
  return GetCovariance(axis1,axis2)/rms1/rms2;
}

//______________________________________________________________________________
Double_t TH5::GetCovariance(Int_t axis1, Int_t axis2) const
{
  //*-*-*-*-*-*-*-*Return covariance between axis1 and axis2*-*-*-*-*
  //*-*            ====================================================
  
  if (axis1 < 1 || axis2 < 1 || axis1 > 5 || axis2 > 5) {
    Error("GetCovariance","Wrong parameters");
    return 0;
  }
  Double_t stats[kNstat5];
  GetStats(stats);
  Double_t sumw   = stats[0];
  Double_t sumw2  = stats[1];
  Double_t sumwx  = stats[2];
  Double_t sumwx2 = stats[3];
  Double_t sumwy  = stats[4];
  Double_t sumwy2 = stats[5];
  Double_t sumwxy = stats[6];
  Double_t sumwz  = stats[7];
  Double_t sumwz2 = stats[8];
  Double_t sumwxz = stats[9];
  Double_t sumwyz = stats[10];
  Double_t sumwu  = stats[11];
  Double_t sumwu2 = stats[12];
  Double_t sumwxu = stats[13];
  Double_t sumwyu = stats[14];
  Double_t sumwzu = stats[15];
  Double_t sumwv  = stats[16];
  Double_t sumwv2 = stats[17];
  Double_t sumwxv = stats[18];
  Double_t sumwyv = stats[19];
  Double_t sumwzv = stats[20];
  Double_t sumwuv = stats[21];

   if (sumw == 0) return 0;
   if (axis1 == 1 && axis2 == 1) {
      return TMath::Abs(sumwx2/sumw - sumwx*sumwx/sumw2);
   }
   if (axis1 == 2 && axis2 == 2) {
      return TMath::Abs(sumwy2/sumw - sumwy*sumwy/sumw2);
   }
   if (axis1 == 3 && axis2 == 3) {
      return TMath::Abs(sumwz2/sumw - sumwz*sumwz/sumw2);
   }
   if (axis1 == 4 && axis2 == 4) {
      return TMath::Abs(sumwu2/sumw - sumwu*sumwu/sumw2);
   }
   if (axis1 == 5 && axis2 == 5) {
      return TMath::Abs(sumwv2/sumw - sumwv*sumwv/sumw2);
   }
   if ((axis1 == 1 && axis2 == 2) || (axis1 == 2 && axis1 == 1)) {
      return sumwxy/sumw - sumwx/sumw*sumwy/sumw;
   }
   if ((axis1 == 1 && axis2 == 3) || (axis1 == 3 && axis2 == 1)) {
      return sumwxz/sumw - sumwx/sumw*sumwz/sumw;
   }
   if ((axis1 == 1 && axis2 == 4) || (axis1 == 4 && axis2 == 1)) {
      return sumwxu/sumw - sumwx/sumw*sumwu/sumw;
   }
   if ((axis1 == 1 && axis2 == 5) || (axis1 == 5 && axis2 == 1)) {
      return sumwxv/sumw - sumwx/sumw*sumwv/sumw;
   }
   if ((axis1 == 2 && axis2 == 3) || (axis1 == 3 && axis2 == 2)) {
      return sumwyz/sumw - sumwy/sumw*sumwz/sumw;
   }
   if ((axis1 == 2 && axis2 == 4) || (axis1 == 4 && axis2 == 2)) {
      return sumwyu/sumw - sumwy/sumw*sumwu/sumw;
   }
   if ((axis1 == 2 && axis2 == 5) || (axis1 == 5 && axis2 == 2)) {
      return sumwyv/sumw - sumwy/sumw*sumwv/sumw;
   }
   if ((axis1 == 3 && axis2 == 4) || (axis1 == 4 && axis2 == 3)) {
      return sumwzu/sumw - sumwz/sumw*sumwu/sumw;
   }
   if ((axis1 == 3 && axis2 == 5) || (axis1 == 5 && axis2 == 3)) {
      return sumwzv/sumw - sumwz/sumw*sumwv/sumw;
   }
   if ((axis1 == 4 && axis2 == 5) || (axis1 == 5 && axis2 == 4)) {
      return sumwuv/sumw - sumwu/sumw*sumwv/sumw;
   }

   return 0;
}


//______________________________________________________________________________
Double_t TH5::GetMaximum(Double_t maxval) const
{
  //  Return maximum value smaller than maxval of bins in the range,
  //  unless the value has been overridden by TH1::SetMaximum,
  //  in which case it returns that value. (This happens, for example,
  //  when the histogram is drawn and the y or z axis limits are changed
  //
  //  To get the maximum value of bins in the histogram regardless of
  //  whether the value has been overridden, use
  //      h->GetBinContent(h->GetMaximumBin())
  
  if (fMaximum != -1111) return fMaximum;
  Int_t bin, binx, biny, binz, binu, binv;
  Int_t xfirst  = fXaxis.GetFirst();
  Int_t xlast   = fXaxis.GetLast();
  Int_t yfirst  = fYaxis.GetFirst();
  Int_t ylast   = fYaxis.GetLast();
  Int_t zfirst  = fZaxis.GetFirst();
  Int_t zlast   = fZaxis.GetLast();
  Int_t ufirst  = fUaxis.GetFirst();
  Int_t ulast   = fUaxis.GetLast();
  Int_t vfirst  = fVaxis.GetFirst();
  Int_t vlast   = fVaxis.GetLast();
  Double_t maximum = -FLT_MAX, value;
  for (binv = vfirst; binv <= vlast; binv++) {
    for (binu = ufirst; binu <= ulast; binu++) {
      for (binz = zfirst; binz <= zlast; binz++) {
	for (biny = yfirst; biny <= ylast; biny++) {
	  for (binx = xfirst; binx <= xlast; binx++) {
            bin = GetBin(binx, biny, binz, binu, binv);
            value = GetBinContent(bin);
            if (value > maximum && value < maxval) maximum = value;
	  }
	}
      }
    }
  }
  return maximum;
}

//______________________________________________________________________________
// void TH5::GetRandom5(Double_t &x, Double_t &y, Double_t &z, Double_t &u, Double_t &v)
// {
//    // return 3 random numbers along axis x , y and z distributed according
//    // the cellcontents of a 3-dim histogram

//    Int_t nbinsx = GetNbinsX();
//    Int_t nbinsy = GetNbinsY();
//    Int_t nbinsz = GetNbinsZ();
//    Int_t nxy    = nbinsx*nbinsy;
//    Int_t nbins  = nxy*nbinsz;
//    Double_t integral;
//    if (fIntegral) {
//       if (fIntegral[nbins+1] != fEntries) integral = ComputeIntegral();
//    } else {
//       integral = ComputeIntegral();
//       if (integral == 0 || fIntegral == 0) return;
//    }
//    Double_t r1 = gRandom->Rndm();
//    Int_t ibin = TMath::BinarySearch(nbins,fIntegral,(Double_t) r1);
//    Int_t binz = ibin/nxy;
//    Int_t biny = (ibin - nxy*binz)/nbinsx;
//    Int_t binx = ibin - nbinsx*(biny + nbinsy*binz);
//    x = fXaxis.GetBinLowEdge(binx+1);
//    if (r1 > fIntegral[ibin]) x +=
//       fXaxis.GetBinWidth(binx+1)*(r1-fIntegral[ibin])/(fIntegral[ibin+1] - fIntegral[ibin]);
//    y = fYaxis.GetBinLowEdge(biny+1) + fYaxis.GetBinWidth(biny+1)*gRandom->Rndm();
//    z = fZaxis.GetBinLowEdge(binz+1) + fZaxis.GetBinWidth(binz+1)*gRandom->Rndm();
// }

//______________________________________________________________________________
void TH5::GetStats(Double_t *stats) const
{
  // fill the array stats from the contents of this histogram
  // The array stats must be correctly dimensionned in the calling program.
  // stats[0] = sumw
  // stats[1] = sumw2
  // stats[2] = sumwx
  // stats[3] = sumwx2
  // stats[4] = sumwy
  // stats[5] = sumwy2
  // stats[6] = sumwxy
  // stats[7] = sumwz
  // stats[8] = sumwz2
  // stats[9] = sumwxz
  // stats[10]= sumwyz
  // stats[11] = sumwu
  // stats[12] = sumwu2
  // stats[13] = sumwxu
  // stats[14] = sumwyu
  // stats[15] = sumwzu
  // stats[16] = sumwv
  // stats[17] = sumwv2
  // stats[18] = sumwxv
  // stats[19] = sumwyv
  // stats[20] = sumwzv
  // stats[21] = sumwuv

  if (fBuffer) ((TH5*)this)->BufferEmpty();
  
  Int_t bin, binx, biny, binz, binu, binv;
  Double_t w,err;
  Double_t x,y,z, u, v;
  if ((fTsumw == 0 && fEntries > 0) || fXaxis.TestBit(TAxis::kAxisRange) || fYaxis.TestBit(TAxis::kAxisRange) ||
      fZaxis.TestBit(TAxis::kAxisRange) || fUaxis.TestBit(TAxis::kAxisRange) || fVaxis.TestBit(TAxis::kAxisRange)) {
    for (bin = 0; bin < 9; bin++) stats[bin] = 0;
    
    Int_t firstBinX = fXaxis.GetFirst();
    Int_t lastBinX  = fXaxis.GetLast();
    Int_t firstBinY = fYaxis.GetFirst();
    Int_t lastBinY  = fYaxis.GetLast();
    Int_t firstBinZ = fZaxis.GetFirst();
    Int_t lastBinZ  = fZaxis.GetLast();
    Int_t firstBinU = fUaxis.GetFirst();
    Int_t lastBinU  = fUaxis.GetLast();
    Int_t firstBinV = fVaxis.GetFirst();
    Int_t lastBinV  = fVaxis.GetLast();
    // include underflow/overflow if TH1::StatOverflows(kTRUE) in case no range is set on the axis
    if (fgStatOverflows) {
      if (!fXaxis.TestBit(TAxis::kAxisRange)) {
	if (firstBinX == 1) firstBinX = 0;
	if (lastBinX ==  fXaxis.GetNbins() ) lastBinX += 1;
      }
      if ( !fYaxis.TestBit(TAxis::kAxisRange) ) {
	if (firstBinY == 1) firstBinY = 0;
	if (lastBinY ==  fYaxis.GetNbins() ) lastBinY += 1;
      }
      if ( !fZaxis.TestBit(TAxis::kAxisRange) ) {
	if (firstBinZ == 1) firstBinZ = 0;
	if (lastBinZ ==  fZaxis.GetNbins() ) lastBinZ += 1;
      }
      if ( !fUaxis.TestBit(TAxis::kAxisRange) ) {
	if (firstBinU == 1) firstBinU = 0;
	if (lastBinU ==  fUaxis.GetNbins() ) lastBinU += 1;
      }
      if ( !fVaxis.TestBit(TAxis::kAxisRange) ) {
	if (firstBinV == 1) firstBinV = 0;
	if (lastBinV ==  fVaxis.GetNbins() ) lastBinV += 1;
      }
    }
    for (binv = firstBinV; binv <= lastBinV; binv++) {
      v = fVaxis.GetBinCenter(binv);
      for (binu = firstBinU; binu <= lastBinU; binu++) {
	u = fUaxis.GetBinCenter(binu);
	for (binz = firstBinZ; binz <= lastBinZ; binz++) {
	  z = fZaxis.GetBinCenter(binz);
	  for (biny = firstBinY; biny <= lastBinY; biny++) {
	    y = fYaxis.GetBinCenter(biny);
	    for (binx = firstBinX; binx <= lastBinX; binx++) {
	      bin = GetBin(binx, biny, binz, binu, binv);
	      x   = fXaxis.GetBinCenter(binx);
	      //w   = TMath::Abs(GetBinContent(bin));
	      w   = GetBinContent(bin);
	      err = TMath::Abs(GetBinError(bin));
	      stats[0] += w;
	      stats[1] += err*err;
	      stats[2] += w*x;
	      stats[3] += w*x*x;
	      stats[4] += w*y;
	      stats[5] += w*y*y;
	      stats[6] += w*x*y;
	      stats[7] += w*z;
	      stats[8] += w*z*z;
	      stats[9] += w*x*z;
	      stats[10] += w*y*z;
	      stats[11] += w*u;
	      stats[12] += w*u*u;
	      stats[13] += w*x*u;
	      stats[14] += w*y*u;
	      stats[15] += w*z*u;
	      stats[16] += w*v;
	      stats[17] += w*v*v;
	      stats[18] += w*x*v;
	      stats[19] += w*y*v;
	      stats[20] += w*z*v;
	      stats[21] += w*u*v;
	    }
	  }
	}
      }
    }
  } else {
    stats[0] = fTsumw;
    stats[1] = fTsumw2;
    stats[2] = fTsumwx;
    stats[3] = fTsumwx2;
    stats[4] = fTsumwy;
    stats[5] = fTsumwy2;
    stats[6] = fTsumwxy;
    stats[7] = fTsumwz;
    stats[8] = fTsumwz2;
    stats[9] = fTsumwxz;
    stats[10] = fTsumwyz;
    stats[11] = fTsumwu;
    stats[12] = fTsumwu2;
    stats[13] = fTsumwxu;
    stats[14] = fTsumwyu;
    stats[15] = fTsumwzu;
    stats[16] = fTsumwv;
    stats[17] = fTsumwv2;
    stats[18] = fTsumwxv;
    stats[19] = fTsumwyv;
    stats[20] = fTsumwzv;
    stats[21] = fTsumwuv;
  }
}

//______________________________________________________________________________
TAxis* TH5::GetUaxis() const
{
  // return a pointer to the U axis object

  return &((TH5*)this)->fUaxis;
}

//______________________________________________________________________________
TAxis* TH5::GetVaxis() const
{
  // return a pointer to the V axis object

  return &((TH5*)this)->fVaxis;
}

//______________________________________________________________________________
// Double_t TH5::Integral(Option_t *option) const
// {
//    //Return integral of bin contents. Only bins in the bins range are considered.
//    // By default the integral is computed as the sum of bin contents in the range.
//    // if option "width" is specified, the integral is the sum of
//    // the bin contents multiplied by the bin width in x, y and in z.

//    return Integral(fXaxis.GetFirst(),fXaxis.GetLast(),
//       fYaxis.GetFirst(),fYaxis.GetLast(),
//       fZaxis.GetFirst(),fZaxis.GetLast(),option);
// }

//______________________________________________________________________________
//Double_t TH5::Integral(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Option_t *option) const
// {
//    //Return integral of bin contents in range [binx1,binx2],[biny1,biny2],[binz1,binz2]
//    // for a 3-D histogram
//    // By default the integral is computed as the sum of bin contents in the range.
//    // if option "width" is specified, the integral is the sum of
//    // the bin contents multiplied by the bin width in x, y and in z.

//    Double_t err = 0; 
//    return DoIntegral(binx1,binx2,biny1,biny2,binz1,binz2,err,option);
// }
//______________________________________________________________________________
//Double_t TH5::IntegralAndError(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Double_t & error, Option_t *option) const
// {
//    //Return integral of bin contents in range [binx1,binx2],[biny1,biny2],[binz1,binz2]
//    // for a 3-D histogram. Calculates also the integral error using error propagation 
//    // from the bin errors assumming that all the bins are uncorrelated. 
//    // By default the integral is computed as the sum of bin contents in the range.
//    // if option "width" is specified, the integral is the sum of
//    // the bin contents multiplied by the bin width in x, y and in z.
//    return DoIntegral(binx1,binx2,biny1,biny2,binz1,binz2,error,option,kTRUE);
// }

//______________________________________________________________________________
// Double_t TH5::Interpolate(Double_t)
// {
   
//    //Not yet implemented
//    Error("Interpolate","This function must be called with 3 arguments for a TH5");
//    return 0;
//}

//______________________________________________________________________________
// Double_t TH5::Interpolate(Double_t, Double_t)
// {
   
//    //Not yet implemented
//    Error("Interpolate","This function must be called with 3 arguments for a TH5");
//    return 0;
//}

//______________________________________________________________________________
// Double_t TH5::Interpolate(Double_t x, Double_t y, Double_t z)
// {
//    // Given a point P(x,y,z), Interpolate approximates the value via trilinear interpolation
//    // based on the 8 nearest bin center points ( corner of the cube surronding the points) 
//    // The Algorithm is described in http://en.wikipedia.org/wiki/Trilinear_interpolation
//    // The given values (x,y,z) must be between first bin center and  last bin center for each coordinate: 
//    //
//    //   fXAxis.GetBinCenter(1) < x  < fXaxis.GetBinCenter(nbinX)     AND
//    //   fYAxis.GetBinCenter(1) < y  < fYaxis.GetBinCenter(nbinY)     AND
//    //   fZAxis.GetBinCenter(1) < z  < fZaxis.GetBinCenter(nbinZ) 

//    Int_t ubx = fXaxis.FindBin(x);
//    if ( x < fXaxis.GetBinCenter(ubx) ) ubx -= 1;
//    Int_t obx = ubx + 1;

//    Int_t uby = fYaxis.FindBin(y);
//    if ( y < fYaxis.GetBinCenter(uby) ) uby -= 1;
//    Int_t oby = uby + 1;

//    Int_t ubz = fZaxis.FindBin(z);
//    if ( z < fZaxis.GetBinCenter(ubz) ) ubz -= 1;
//    Int_t obz = ubz + 1;


////    if ( IsBinUnderflow(GetBin(ubx, uby, ubz)) ||
////         IsBinOverflow (GetBin(obx, oby, obz)) ) {
//    if (ubx <=0 || uby <=0 || ubz <= 0 ||
//        obx > fXaxis.GetNbins() || oby > fYaxis.GetNbins() || obz > fZaxis.GetNbins() ) {
//       Error("Interpolate","Cannot interpolate outside histogram domain.");
//       return 0;
//    }

//    Double_t xw = fXaxis.GetBinCenter(obx) - fXaxis.GetBinCenter(ubx);
//    Double_t yw = fYaxis.GetBinCenter(oby) - fYaxis.GetBinCenter(uby);
//    Double_t zw = fZaxis.GetBinCenter(obz) - fZaxis.GetBinCenter(ubz);

//    Double_t xd = (x - fXaxis.GetBinCenter(ubx)) / xw;
//    Double_t yd = (y - fYaxis.GetBinCenter(uby)) / yw;
//    Double_t zd = (z - fZaxis.GetBinCenter(ubz)) / zw;


//    Double_t v[] = { GetBinContent( ubx, uby, ubz ), GetBinContent( ubx, uby, obz ),
//                     GetBinContent( ubx, oby, ubz ), GetBinContent( ubx, oby, obz ),
//                     GetBinContent( obx, uby, ubz ), GetBinContent( obx, uby, obz ),
//                     GetBinContent( obx, oby, ubz ), GetBinContent( obx, oby, obz ) };


//    Double_t i1 = v[0] * (1 - zd) + v[1] * zd;
//    Double_t i2 = v[2] * (1 - zd) + v[3] * zd;
//    Double_t j1 = v[4] * (1 - zd) + v[5] * zd;
//    Double_t j2 = v[6] * (1 - zd) + v[7] * zd;


//    Double_t w1 = i1 * (1 - yd) + i2 * yd;
//    Double_t w2 = j1 * (1 - yd) + j2 * yd;


//    Double_t result = w1 * (1 - xd) + w2 * xd;

//    return result;
// }

//______________________________________________________________________________
//Double_t TH5::KolmogorovTest(const TH1 *h2, Option_t *option) const
//{
//    //  Statistical test of compatibility in shape between
//    //  THIS histogram and h2, using Kolmogorov test.
//    //     Default: Ignore under- and overflow bins in comparison
//    //
//    //     option is a character string to specify options
//    //         "U" include Underflows in test
//    //         "O" include Overflows
//    //         "N" include comparison of normalizations
//    //         "D" Put out a line of "Debug" printout
//    //         "M" Return the Maximum Kolmogorov distance instead of prob
//    //
//    //   The returned function value is the probability of test
//    //       (much less than one means NOT compatible)
//    //
//    //   The KS test uses the distance between the pseudo-CDF's obtained 
//    //   from the histogram. Since in more than 1D the order for generating the pseudo-CDF is 
//    //   arbitrary, we use the pseudo-CDF's obtained from all the possible 6 combinatons of the 3 axis. 
//    //   The average of all the maximum  distances obtained is used in the tests.  

//    TString opt = option;
//    opt.ToUpper();

//    Double_t prb = 0;
//    TH1 *h1 = (TH1*)this;
//    if (h2 == 0) return 0;
//    TAxis *xaxis1 = h1->GetXaxis();
//    TAxis *xaxis2 = h2->GetXaxis();
//    TAxis *yaxis1 = h1->GetYaxis();
//    TAxis *yaxis2 = h2->GetYaxis();
//    TAxis *zaxis1 = h1->GetZaxis();
//    TAxis *zaxis2 = h2->GetZaxis();
//    Int_t ncx1   = xaxis1->GetNbins();
//    Int_t ncx2   = xaxis2->GetNbins();
//    Int_t ncy1   = yaxis1->GetNbins();
//    Int_t ncy2   = yaxis2->GetNbins();
//    Int_t ncz1   = zaxis1->GetNbins();
//    Int_t ncz2   = zaxis2->GetNbins();

//    // Check consistency of dimensions
//    if (h1->GetDimension() != 3 || h2->GetDimension() != 3) {
//       Error("KolmogorovTest","Histograms must be 3-D\n");
//       return 0;
//    }

//    // Check consistency in number of channels
//    if (ncx1 != ncx2) {
//       Error("KolmogorovTest","Number of channels in X is different, %d and %d\n",ncx1,ncx2);
//       return 0;
//    }
//    if (ncy1 != ncy2) {
//       Error("KolmogorovTest","Number of channels in Y is different, %d and %d\n",ncy1,ncy2);
//       return 0;
//    }
//    if (ncz1 != ncz2) {
//       Error("KolmogorovTest","Number of channels in Z is different, %d and %d\n",ncz1,ncz2);
//       return 0;
//    }

//    // Check consistency in channel edges
//    Bool_t afunc1 = kFALSE;
//    Bool_t afunc2 = kFALSE;
//    Double_t difprec = 1e-5;
//    Double_t diff1 = TMath::Abs(xaxis1->GetXmin() - xaxis2->GetXmin());
//    Double_t diff2 = TMath::Abs(xaxis1->GetXmax() - xaxis2->GetXmax());
//    if (diff1 > difprec || diff2 > difprec) {
//       Error("KolmogorovTest","histograms with different binning along X");
//       return 0;
//    }
//    diff1 = TMath::Abs(yaxis1->GetXmin() - yaxis2->GetXmin());
//    diff2 = TMath::Abs(yaxis1->GetXmax() - yaxis2->GetXmax());
//    if (diff1 > difprec || diff2 > difprec) {
//       Error("KolmogorovTest","histograms with different binning along Y");
//       return 0;
//    }
//    diff1 = TMath::Abs(zaxis1->GetXmin() - zaxis2->GetXmin());
//    diff2 = TMath::Abs(zaxis1->GetXmax() - zaxis2->GetXmax());
//    if (diff1 > difprec || diff2 > difprec) {
//       Error("KolmogorovTest","histograms with different binning along Z");
//       return 0;
//    }

//    //   Should we include Uflows, Oflows?
//    Int_t ibeg = 1, jbeg = 1, kbeg = 1;
//    Int_t iend = ncx1, jend = ncy1, kend = ncz1;
//    if (opt.Contains("U")) {ibeg = 0; jbeg = 0; kbeg = 0;}
//    if (opt.Contains("O")) {iend = ncx1+1; jend = ncy1+1; kend = ncz1+1;}

//    Int_t i,j,k,bin;
//    Double_t sum1  = 0;
//    Double_t sum2  = 0;
//    Double_t w1    = 0;
//    Double_t w2    = 0;
//    for (i = ibeg; i <= iend; i++) {
//       for (j = jbeg; j <= jend; j++) {
//          for (k = kbeg; k <= kend; k++) {
//             bin = h1->GetBin(i,j,k);
//             sum1 += h1->GetBinContent(bin);
//             sum2 += h2->GetBinContent(bin);
//             Double_t ew1   = h1->GetBinError(bin);
//             Double_t ew2   = h2->GetBinError(bin);
//             w1   += ew1*ew1;
//             w2   += ew2*ew2;
//          }
//       }
//    }


//    //    Check that both scatterplots contain events
//    if (sum1 == 0) {
//       Error("KolmogorovTest","Integral is zero for h1=%s\n",h1->GetName());
//       return 0;
//    }
//    if (sum2 == 0) {
//       Error("KolmogorovTest","Integral is zero for h2=%s\n",h2->GetName());
//       return 0;
//    }
//    // calculate the effective entries.  
//    // the case when errors are zero (w1 == 0 or w2 ==0) are equivalent to 
//    // compare to a function. In that case the rescaling is done only on sqrt(esum2) or sqrt(esum1) 
//    Double_t esum1 = 0, esum2 = 0; 
//    if (w1 > 0) 
//       esum1 = sum1 * sum1 / w1; 
//    else 
//       afunc1 = kTRUE;    // use later for calculating z
   
//    if (w2 > 0) 
//       esum2 = sum2 * sum2 / w2; 
//    else 
//       afunc2 = kTRUE;    // use later for calculating z
   
//    if (afunc2 && afunc1) { 
//       Error("KolmogorovTest","Errors are zero for both histograms\n");
//       return 0;
//    }

//    //   Find Kolmogorov distance
//    //   order is arbitrary take average of all possible 6 starting orders x,y,z 
//    int order[3] = {0,1,2};
//    int binbeg[3]; 
//    int binend[3]; 
//    int ibin[3];
//    binbeg[0] = ibeg; binbeg[1] = jbeg; binbeg[2] = kbeg; 
//    binend[0] = iend; binend[1] = jend; binend[2] = kend; 
//    Double_t vdfmax[6]; // there are in total 6 combinations 
//    int icomb = 0; 
//    Double_t s1 = 1./(6.*sum1);
//    Double_t s2 = 1./(6.*sum2);
//    Double_t rsum1=0, rsum2=0;
//    do { 
//       // loop on bins
//       Double_t dmax = 0;
//       for (i = binbeg[order[0] ]; i <= binend[order[0] ]; i++) {
//          for ( j = binbeg[order[1] ]; j <= binend[order[1] ]; j++) {
//             for ( k = binbeg[order[2] ]; k <= binend[order[2] ]; k++) {
//                   ibin[ order[0] ] = i;
//                   ibin[ order[1] ] = j;
//                   ibin[ order[2] ] = k;
//                   bin = h1->GetBin(ibin[0],ibin[1],ibin[2]);
//                   rsum1 += s1*h1->GetBinContent(bin);
//                   rsum2 += s2*h2->GetBinContent(bin);
//                   dmax   = TMath::Max(dmax, TMath::Abs(rsum1-rsum2));
//             }
//          }
//       }
//       vdfmax[icomb] = dmax; 
//       icomb++;
//    } while (TMath::Permute(3,order)  );


//    // get average of distances 
//    Double_t dfmax = TMath::Mean(6,vdfmax);
   
//    //    Get Kolmogorov probability
//    Double_t factnm;
//    if (afunc1)      factnm = TMath::Sqrt(sum2);
//    else if (afunc2) factnm = TMath::Sqrt(sum1);
//    else             factnm = TMath::Sqrt(sum1*sum2/(sum1+sum2));
//    Double_t z  = dfmax*factnm;

//    prb = TMath::KolmogorovProb(z); 

//    Double_t prb1 = 0, prb2 = 0; 
//    // option N to combine normalization makes sense if both afunc1 and afunc2 are false
//    if (opt.Contains("N")  && !(afunc1 || afunc2 ) ) { 
//       // Combine probabilities for shape and normalization
//       prb1   = prb;
//       Double_t d12    = esum1-esum2;
//       Double_t chi2   = d12*d12/(esum1+esum2);
//       prb2   = TMath::Prob(chi2,1);
//       //     see Eadie et al., section 11.6.2
//       if (prb > 0 && prb2 > 0) prb = prb*prb2*(1-TMath::Log(prb*prb2));
//       else                     prb = 0;
//    }

//    //    debug printout
//    if (opt.Contains("D")) {
//       printf(" Kolmo Prob  h1 = %s, sum1=%g\n",h1->GetName(),sum1);
//       printf(" Kolmo Prob  h2 = %s, sum2=%g\n",h2->GetName(),sum2);
//       printf(" Kolmo Probabil = %f, Max Dist = %g\n",prb,dfmax);
//       if (opt.Contains("N"))
//          printf(" Kolmo Probabil = %f for shape alone, =%f for normalisation alone\n",prb1,prb2);
//    }
//    // This numerical error condition should never occur:
//    if (TMath::Abs(rsum1-1) > 0.002) Warning("KolmogorovTest","Numerical problems with h1=%s\n",h1->GetName());
//    if (TMath::Abs(rsum2-1) > 0.002) Warning("KolmogorovTest","Numerical problems with h2=%s\n",h2->GetName());

//    if(opt.Contains("M"))      return dfmax;  // return avergae of max distance

//    return prb;
// }


//______________________________________________________________________________
// Long64_t TH5::Merge(TCollection *list)
// {
//    //Add all histograms in the collection to this histogram.
//    //This function computes the min/max for the axes,
//    //compute a new number of bins, if necessary,
//    //add bin contents, errors and statistics.
//    //If overflows are present and limits are different the function will fail.
//    //The function returns the total number of entries in the result histogram
//    //if the merge is successfull, -1 otherwise.
//    //
//    //IMPORTANT remark. The 2 axis x and y may have different number
//    //of bins and different limits, BUT the largest bin width must be
//    //a multiple of the smallest bin width and the upper limit must also
//    //be a multiple of the bin width.

//    if (!list) return 0;
//    if (list->IsEmpty()) return (Long64_t) GetEntries();

//    TList inlist;
//    inlist.AddAll(list);

//    TAxis newXAxis;
//    TAxis newYAxis;
//    TAxis newZAxis;
//    Bool_t initialLimitsFound = kFALSE;
//    Bool_t allSameLimits = kTRUE;
//    Bool_t allHaveLimits = kTRUE;
//    Bool_t firstNonEmptyHist = kTRUE;

//    TIter next(&inlist);
//    TH5* h = this;
//    do { 
//       // skip empty histograms 
//       if (h->fTsumw == 0 && h->GetEntries() == 0) continue;

//       Bool_t hasLimits = h->GetXaxis()->GetXmin() < h->GetXaxis()->GetXmax();
//       allHaveLimits = allHaveLimits && hasLimits;

//       if (hasLimits) {
//          h->BufferEmpty();

//          // this is done in case the first histograms are empty and 
//          // the histogram have different limits
//          if (firstNonEmptyHist ) { 
//             // set axis limits in the case the first histogram was empty
//             if (h != this ) { 
//                if (!SameLimitsAndNBins(fXaxis, *(h->GetXaxis())) ) 
//                   fXaxis.Set(h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax());
//                if (!SameLimitsAndNBins(fYaxis, *(h->GetYaxis())) ) 
//                   fYaxis.Set(h->GetYaxis()->GetNbins(), h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
//                if (!SameLimitsAndNBins(fZaxis, *(h->GetZaxis())) ) 
//                   fZaxis.Set(h->GetZaxis()->GetNbins(), h->GetZaxis()->GetXmin(),h->GetZaxis()->GetXmax());
//             }
//             firstNonEmptyHist = kFALSE;     
//          }

//          if (!initialLimitsFound) {
//             // this is executed the first time an histogram with limits is found
//             // to set some initial values on the new axes
//             initialLimitsFound = kTRUE;
//             newXAxis.Set(h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
//                h->GetXaxis()->GetXmax());
//             newYAxis.Set(h->GetYaxis()->GetNbins(), h->GetYaxis()->GetXmin(),
//                h->GetYaxis()->GetXmax());
//             newZAxis.Set(h->GetZaxis()->GetNbins(), h->GetZaxis()->GetXmin(),
//                h->GetZaxis()->GetXmax());
//          }
//          else {
//             // check first if histograms have same bins 
//             if (!SameLimitsAndNBins(newXAxis, *(h->GetXaxis())) || 
//                 !SameLimitsAndNBins(newYAxis, *(h->GetYaxis())) ||
//                 !SameLimitsAndNBins(newZAxis, *(h->GetZaxis())) ) { 
               
//                allSameLimits = kFALSE;
 
//                if (!RecomputeAxisLimits(newXAxis, *(h->GetXaxis()))) {
//                   Error("Merge", "Cannot merge histograms - limits are inconsistent:\n "
//                   "first: (%d, %f, %f), second: (%d, %f, %f)",
//                         newXAxis.GetNbins(), newXAxis.GetXmin(), newXAxis.GetXmax(),
//                         h->GetXaxis()->GetNbins(), h->GetXaxis()->GetXmin(),
//                         h->GetXaxis()->GetXmax());
//                   return -1; 
//                }
//                if (!RecomputeAxisLimits(newYAxis, *(h->GetYaxis()))) {
//                   Error("Merge", "Cannot merge histograms - limits are inconsistent:\n "
//                         "first: (%d, %f, %f), second: (%d, %f, %f)",
//                         newYAxis.GetNbins(), newYAxis.GetXmin(), newYAxis.GetXmax(),
//                         h->GetYaxis()->GetNbins(), h->GetYaxis()->GetXmin(),
//                         h->GetYaxis()->GetXmax());
//                   return -1; 
//                }
//                if (!RecomputeAxisLimits(newZAxis, *(h->GetZaxis()))) {
//                   Error("Merge", "Cannot merge histograms - limits are inconsistent:\n "
//                         "first: (%d, %f, %f), second: (%d, %f, %f)",
//                         newZAxis.GetNbins(), newZAxis.GetXmin(), newZAxis.GetXmax(),
//                         h->GetZaxis()->GetNbins(), h->GetZaxis()->GetXmin(),
//                         h->GetZaxis()->GetXmax());
//                   return -1; 
//                }
//             }
//          }
//       }
//    } while ( ( h = dynamic_cast<TH5*> ( next() ) ) != NULL );
//    if (!h && (*next) ) {
//       Error("Merge","Attempt to merge object of class: %s to a %s",
//             (*next)->ClassName(),this->ClassName());
//       return -1;
//    }
//    next.Reset();

//    // In the case of histogram with different limits
//    // newX(Y)Axis will now have the new found limits
//    // but one needs first to clone this histogram to perform the merge
//    // The clone is not needed when all histograms have the same limits
//    TH5 * hclone = 0;
//    if (!allSameLimits) { 
//       // We don't want to add the clone to gDirectory,
//       // so remove our kMustCleanup bit temporarily
//       Bool_t mustCleanup = TestBit(kMustCleanup);
//       if (mustCleanup) ResetBit(kMustCleanup);
//       hclone = (TH5*)IsA()->New();
//       hclone->SetDirectory(0);
//       Copy(*hclone);
//       if (mustCleanup) SetBit(kMustCleanup);
//       BufferEmpty(1);         // To remove buffer.
//       Reset();                // BufferEmpty sets limits so we can't use it later.
//       SetEntries(0);
//       inlist.AddFirst(hclone);
//    }

//    if (!allSameLimits && initialLimitsFound) {
//       SetBins(newXAxis.GetNbins(), newXAxis.GetXmin(), newXAxis.GetXmax(),
//               newYAxis.GetNbins(), newYAxis.GetXmin(), newYAxis.GetXmax(),
//               newZAxis.GetNbins(), newZAxis.GetXmin(), newZAxis.GetXmax());
//    }

//    if (!allHaveLimits) {
//       // fill this histogram with all the data from buffers of histograms without limits
//       while ( (h = dynamic_cast<TH5*> (next())) ) {
//          if (h->GetXaxis()->GetXmin() >= h->GetXaxis()->GetXmax() && h->fBuffer) {
//             // no limits
//             Int_t nbentries = (Int_t)h->fBuffer[0];
//             for (Int_t i = 0; i < nbentries; i++)
//                Fill(h->fBuffer[4*i + 2], h->fBuffer[4*i + 3],
//                h->fBuffer[4*i + 4], h->fBuffer[4*i + 1]);
//             // Entries from buffers have to be filled one by one
//             // because FillN doesn't resize histograms.
//          }
//       }
//       if (!initialLimitsFound) {
//          if (hclone) { 
//             inlist.Remove(hclone);
//             delete hclone; 
//          }
//          return (Long64_t) GetEntries();  // all histograms have been processed
//       }
//       next.Reset();
//    }

//    //merge bin contents and errors
//    Double_t stats[kNstat5], totstats[kNstat5];
//    for (Int_t i=0;i<kNstat5;i++) {totstats[i] = stats[i] = 0;}
//    GetStats(totstats);
//    Double_t nentries = GetEntries();
//    Int_t binx, biny, binz, ix, iy, iz, nx, ny, nz, bin, ibin;
//    Double_t cu;
//    Int_t nbix = fXaxis.GetNbins();
//    Int_t nbiy = fYaxis.GetNbins();
//    Bool_t canRebin=TestBit(kCanRebin);
//    ResetBit(kCanRebin); // reset, otherwise setting the under/overflow will rebin

//    while ( (h=(TH5*)next()) ) {
//       // process only if the histogram has limits; otherwise it was processed before
//       if (h->GetXaxis()->GetXmin() < h->GetXaxis()->GetXmax()) {
//          // import statistics
//          h->GetStats(stats);
//          for (Int_t i = 0; i < kNstat5; i++)
//             totstats[i] += stats[i];
//          nentries += h->GetEntries();

//          nx = h->GetXaxis()->GetNbins();
//          ny = h->GetYaxis()->GetNbins();
//          nz = h->GetZaxis()->GetNbins();
         
//          // mantain loop in separate binz, biny and binz to avoid 
//          // callinig FindBin(x,y,z) for every bin
//          for (binz = 0; binz <= nz + 1; binz++) {
//             if (!allSameLimits)
//                iz = fZaxis.FindBin(h->GetZaxis()->GetBinCenter(binz));
//             else
//                iz = binz;
//             for (biny = 0; biny <= ny + 1; biny++) {
//                if (!allSameLimits)
//                   iy = fYaxis.FindBin(h->GetYaxis()->GetBinCenter(biny));
//                else 
//                   iy = biny; 

//                for (binx = 0; binx <= nx + 1; binx++) {
//                   bin = binx +(nx+2)*(biny + (ny+2)*binz);
//                   cu  = h->GetBinContent(bin);
//                   if (!allSameLimits) { 
//                      // look at non-empty unerflow/overflows
//                      if (cu != 0 && ( h->IsBinUnderflow(bin) || h->IsBinOverflow(bin) )) {
//                         Error("Merge", "Cannot merge histograms - the histograms have"
//                               " different limits and undeflows/overflows are present."
//                               " The initial histogram is now broken!");
//                         return -1;
//                      }
//                      ix = fXaxis.FindBin(h->GetXaxis()->GetBinCenter(binx));
//                   }
//                   else { 
//                      // case histograms have same limits 
//                      ix = binx; 
//                   }
                  
//                   ibin = ix +(nbix+2)*(iy + (nbiy+2)*iz);
//                   if (ibin <0) continue;
//                   AddBinContent(ibin,cu);
//                   if (fSumw2.fN) {
//                      Double_t error1 = h->GetBinError(bin);
//                      fSumw2.fArray[ibin] += error1*error1;
//                   }
//                }
//             }
//          }
//       }
//    }
//    if (canRebin) SetBit(kCanRebin);

//    //copy merged stats
//    PutStats(totstats);
//    SetEntries(nentries);
//    if (hclone) { 
//       inlist.Remove(hclone);
//       delete hclone;
//    }
//    return (Long64_t)nentries;
// }

//______________________________________________________________________________
TH1D *TH5::ProjectionX(const char *name, Int_t iymin, Int_t iymax, Int_t izmin, Int_t izmax,
		       Int_t iumin, Int_t iumax, Int_t ivmin, Int_t ivmax, Option_t *option) const
{
  //*-*-*-*-*Project a 5-D histogram into a 1-D histogram along X*-*-*-*-*-*-*
  //*-*      ====================================================
  //
  //   The projection is always of the type TH1D.
  //   The projection is made from the cells along the X axis
  //   ranging from iymin to iymax, izmin to izmax,... included.
  //   By default, underflow and overflows are included
  //   By Setting iymin=1 and iymax=NbinsY the underflow and/or overflow will be excluded
  //
  //   if option "e" is specified, the errors are computed.
  //   if option "d" is specified, the projection is drawn in the current pad.
  //   if option "o" original axis range of the target axes will be
  //   kept, but only bins inside the selected range will be filled.
  //
  //   NOTE that if a TH1D named "name" exists in the current directory or pad 
  //   the histogram is reset and filled again with the projected contents of the TH5.
  //
  //  implemented using Project5D
  
  
  TString opt = option;
  opt.ToLower();
  
  Int_t piymin = GetYaxis()->GetFirst();
  Int_t piymax = GetYaxis()->GetLast();
  Int_t pizmin = GetZaxis()->GetFirst();
  Int_t pizmax = GetZaxis()->GetLast();   
  Int_t piumin = GetUaxis()->GetFirst();
  Int_t piumax = GetUaxis()->GetLast();   
  Int_t pivmin = GetVaxis()->GetFirst();
  Int_t pivmax = GetVaxis()->GetLast();   
  
  GetYaxis()->SetRange(iymin,iymax);
  GetZaxis()->SetRange(izmin,izmax);
  GetUaxis()->SetRange(iumin,iumax);
  GetVaxis()->SetRange(ivmin,ivmax);
  
  // exclude underflow/overflow by forcing the axis range bit
  // due to limitation of TAxis::SetRange cannot select only underflow or overflow cannot have underflow or overflow only   
  if (iymin == 1 && iymax ==  GetNbinsY() ) GetYaxis()->SetBit(TAxis::kAxisRange); 
  if (izmin == 1 && izmax ==  GetNbinsZ() ) GetZaxis()->SetBit(TAxis::kAxisRange); 
  if (iumin == 1 && iumax ==  GetNbinsU() ) GetUaxis()->SetBit(TAxis::kAxisRange); 
  if (ivmin == 1 && ivmax ==  GetNbinsV() ) GetVaxis()->SetBit(TAxis::kAxisRange); 
  // I can use the useUF or useOF flag for exclude/include the underflow/overflow separately
  // (-1 in the imax is considered as an overflow)
  Bool_t useUF = (iymin == 0 || izmin == 0 || iumin == 0 || ivmin == 0); 
  Bool_t useOF = (iymax < 0) || (iymax > GetNbinsY()) || (izmax < 0) || (izmax > GetNbinsZ()) ||
                 (iumax < 0) || (iumax > GetNbinsU()) || (ivmax < 0) || (ivmax > GetNbinsV()); 
		  
  Bool_t computeErrors = GetSumw2N();
  if (opt.Contains("e") ) { 
    computeErrors = kTRUE;
    opt.Remove(opt.First("e"),1);
  }
  Bool_t originalRange = kFALSE;
  if (opt.Contains('o') ) { 
    originalRange = kTRUE; 
    opt.Remove(opt.First("o"),1);
  }
  
  TH1D * h1 = DoProject1D(name, GetTitle(), this->GetXaxis(), computeErrors, originalRange, useUF, useOF);
  
  // restore original range
  GetYaxis()->SetRange(piymin,piymax);
  GetZaxis()->SetRange(pizmin,pizmax);
  GetUaxis()->SetRange(piumin,piumax);
  GetVaxis()->SetRange(pivmin,pivmax);
  
  // draw in current pad 
  if (h1 && opt.Contains("d")) {
    opt.Remove(opt.First("d"),1);
    TVirtualPad *padsav = gPad;
    TVirtualPad *pad = gROOT->GetSelectedPad();
    if (pad) pad->cd();
    if (!gPad->FindObject(h1)) {
      h1->Draw(opt);
    } else {
      h1->Paint(opt);
    }
    if (padsav) padsav->cd();
  }
  
  return h1;
}

//______________________________________________________________________________
TH1D *TH5::ProjectionY(const char *name, Int_t ixmin, Int_t ixmax, Int_t izmin, Int_t izmax,
		       Int_t iumin, Int_t iumax, Int_t ivmin, Int_t ivmax, Option_t *option) const
{
  //*-*-*-*-*Project a 5-D histogram into a 1-D histogram along Y*-*-*-*-*-*-*
  //*-*      ====================================================
  //
  //   The projection is always of the type TH1D.
  //   The projection is made from the cells along the Y axis
  //   ranging from ixmin to ixmax, izmin to izmax,... included.
  //   By default, underflow and overflows are included
  //   By Setting ixmin=1 and ixmax=NbinsX the underflow and/or overflow will be excluded
  //
  //   if option "e" is specified, the errors are computed.
  //   if option "d" is specified, the projection is drawn in the current pad.
  //   if option "o" original axis range of the target axes will be
  //   kept, but only bins inside the selected range will be filled.
  //
  //   NOTE that if a TH1D named "name" exists in the current directory or pad 
  //   the histogram is reset and filled again with the projected contents of the TH5.
  //
  //  implemented using Project5D
  
  
  TString opt = option;
  opt.ToLower();
  
  Int_t pixmin = GetXaxis()->GetFirst();
  Int_t pixmax = GetXaxis()->GetLast();
  Int_t pizmin = GetZaxis()->GetFirst();
  Int_t pizmax = GetZaxis()->GetLast();   
  Int_t piumin = GetUaxis()->GetFirst();
  Int_t piumax = GetUaxis()->GetLast();   
  Int_t pivmin = GetVaxis()->GetFirst();
  Int_t pivmax = GetVaxis()->GetLast();   
  
  GetXaxis()->SetRange(ixmin,ixmax);
  GetZaxis()->SetRange(izmin,izmax);
  GetUaxis()->SetRange(iumin,iumax);
  GetVaxis()->SetRange(ivmin,ivmax);
  
  // exclude underflow/overflow by forcing the axis range bit
  // due to limitation of TAxis::SetRange cannot select only underflow or overflow cannot have underflow or overflow only   
  if (ixmin == 1 && ixmax ==  GetNbinsX() ) GetXaxis()->SetBit(TAxis::kAxisRange); 
  if (izmin == 1 && izmax ==  GetNbinsZ() ) GetZaxis()->SetBit(TAxis::kAxisRange); 
  if (iumin == 1 && iumax ==  GetNbinsU() ) GetUaxis()->SetBit(TAxis::kAxisRange); 
  if (ivmin == 1 && ivmax ==  GetNbinsV() ) GetVaxis()->SetBit(TAxis::kAxisRange); 
  // I can use the useUF or useOF flag for exclude/include the underflow/overflow separately
  // (-1 in the imax is considered as an overflow)
  Bool_t useUF = (ixmin == 0 || izmin == 0 || iumin == 0 || ivmin == 0); 
  Bool_t useOF = (ixmax < 0) || (ixmax > GetNbinsX()) || (izmax < 0) || (izmax > GetNbinsZ()) ||
                 (iumax < 0) || (iumax > GetNbinsU()) || (ivmax < 0) || (ivmax > GetNbinsV()); 
		  
  Bool_t computeErrors = GetSumw2N();
  if (opt.Contains("e") ) { 
    computeErrors = kTRUE;
    opt.Remove(opt.First("e"),1);
  }
  Bool_t originalRange = kFALSE;
  if (opt.Contains('o') ) { 
    originalRange = kTRUE; 
    opt.Remove(opt.First("o"),1);
  }
  
  TH1D * h1 = DoProject1D(name, GetTitle(), this->GetYaxis(), computeErrors, originalRange, useUF, useOF);
  
  // restore original range
  GetXaxis()->SetRange(pixmin,pixmax);
  GetZaxis()->SetRange(pizmin,pizmax);
  GetUaxis()->SetRange(piumin,piumax);
  GetVaxis()->SetRange(pivmin,pivmax);
  
  // draw in current pad 
  if (h1 && opt.Contains("d")) {
    opt.Remove(opt.First("d"),1);
    TVirtualPad *padsav = gPad;
    TVirtualPad *pad = gROOT->GetSelectedPad();
    if (pad) pad->cd();
    if (!gPad->FindObject(h1)) {
      h1->Draw(opt);
    } else {
      h1->Paint(opt);
    }
    if (padsav) padsav->cd();
  }
  
  return h1;
}

//______________________________________________________________________________
TH1D *TH5::ProjectionZ(const char *name, Int_t ixmin, Int_t ixmax, Int_t iymin, Int_t iymax,
		       Int_t iumin, Int_t iumax, Int_t ivmin, Int_t ivmax, Option_t *option) const
{
  //*-*-*-*-*Project a 5-D histogram into a 1-D histogram along Z*-*-*-*-*-*-*
  //*-*      ====================================================
  //
  //   The projection is always of the type TH1D.
  //   The projection is made from the cells along the Z axis
  //   ranging from ixmin to ixmax, iymin to iymax,... included.
  //   By default, underflow and overflows are included
  //   By Setting ixmin=1 and ixmax=NbinsX the underflow and/or overflow will be excluded
  //
  //   if option "e" is specified, the errors are computed.
  //   if option "d" is specified, the projection is drawn in the current pad.
  //   if option "o" original axis range of the target axes will be
  //   kept, but only bins inside the selected range will be filled.
  //
  //   NOTE that if a TH1D named "name" exists in the current directory or pad 
  //   the histogram is reset and filled again with the projected contents of the TH5.
  //
  //  implemented using Project5D
  
  
  TString opt = option;
  opt.ToLower();
  
  Int_t pixmin = GetXaxis()->GetFirst();
  Int_t pixmax = GetXaxis()->GetLast();
  Int_t piymin = GetYaxis()->GetFirst();
  Int_t piymax = GetYaxis()->GetLast();   
  Int_t piumin = GetUaxis()->GetFirst();
  Int_t piumax = GetUaxis()->GetLast();   
  Int_t pivmin = GetVaxis()->GetFirst();
  Int_t pivmax = GetVaxis()->GetLast();   
  
  GetXaxis()->SetRange(ixmin,ixmax);
  GetYaxis()->SetRange(iymin,iymax);
  GetUaxis()->SetRange(iumin,iumax);
  GetVaxis()->SetRange(ivmin,ivmax);
  
  // exclude underflow/overflow by forcing the axis range bit
  // due to limitation of TAxis::SetRange cannot select only underflow or overflow cannot have underflow or overflow only   
  if (ixmin == 1 && ixmax ==  GetNbinsX() ) GetXaxis()->SetBit(TAxis::kAxisRange); 
  if (iymin == 1 && iymax ==  GetNbinsY() ) GetYaxis()->SetBit(TAxis::kAxisRange); 
  if (iumin == 1 && iumax ==  GetNbinsU() ) GetUaxis()->SetBit(TAxis::kAxisRange); 
  if (ivmin == 1 && ivmax ==  GetNbinsV() ) GetVaxis()->SetBit(TAxis::kAxisRange); 
  // I can use the useUF or useOF flag for exclude/include the underflow/overflow separately
  // (-1 in the imax is considered as an overflow)
  Bool_t useUF = (ixmin == 0 || iymin == 0 || iumin == 0 || ivmin == 0); 
  Bool_t useOF = (ixmax < 0) || (ixmax > GetNbinsX()) || (iymax < 0) || (iymax > GetNbinsY()) ||
                 (iumax < 0) || (iumax > GetNbinsU()) || (ivmax < 0) || (ivmax > GetNbinsV()); 
		  
  Bool_t computeErrors = GetSumw2N();
  if (opt.Contains("e") ) { 
    computeErrors = kTRUE;
    opt.Remove(opt.First("e"),1);
  }
  Bool_t originalRange = kFALSE;
  if (opt.Contains('o') ) { 
    originalRange = kTRUE; 
    opt.Remove(opt.First("o"),1);
  }
  
  TH1D * h1 = DoProject1D(name, GetTitle(), this->GetZaxis(), computeErrors, originalRange, useUF, useOF);
  
  // restore original range
  GetXaxis()->SetRange(pixmin,pixmax);
  GetYaxis()->SetRange(piymin,piymax);
  GetUaxis()->SetRange(piumin,piumax);
  GetVaxis()->SetRange(pivmin,pivmax);
  
  // draw in current pad 
  if (h1 && opt.Contains("d")) {
    opt.Remove(opt.First("d"),1);
    TVirtualPad *padsav = gPad;
    TVirtualPad *pad = gROOT->GetSelectedPad();
    if (pad) pad->cd();
    if (!gPad->FindObject(h1)) {
      h1->Draw(opt);
    } else {
      h1->Paint(opt);
    }
    if (padsav) padsav->cd();
  }
  
  return h1;
}

//______________________________________________________________________________
TH1D *TH5::ProjectionU(const char *name, Int_t ixmin, Int_t ixmax, Int_t iymin, Int_t iymax,
		       Int_t izmin, Int_t izmax, Int_t ivmin, Int_t ivmax, Option_t *option) const
{
  //*-*-*-*-*Project a 5-D histogram into a 1-D histogram along U*-*-*-*-*-*-*
  //*-*      ====================================================
  //
  //   The projection is always of the type TH1D.
  //   The projection is made from the cells along the U axis
  //   ranging from ixmin to ixmax, iymin to iymax,... included.
  //   By default, underflow and overflows are included
  //   By Setting ixmin=1 and ixmax=NbinsX the underflow and/or overflow will be excluded
  //
  //   if option "e" is specified, the errors are computed.
  //   if option "d" is specified, the projection is drawn in the current pad.
  //   if option "o" original axis range of the target axes will be
  //   kept, but only bins inside the selected range will be filled.
  //
  //   NOTE that if a TH1D named "name" exists in the current directory or pad 
  //   the histogram is reset and filled again with the projected contents of the TH5.
  //
  //  implemented using Project5D
  
  
  TString opt = option;
  opt.ToLower();
  
  Int_t pixmin = GetXaxis()->GetFirst();
  Int_t pixmax = GetXaxis()->GetLast();
  Int_t piymin = GetYaxis()->GetFirst();
  Int_t piymax = GetYaxis()->GetLast();   
  Int_t pizmin = GetZaxis()->GetFirst();
  Int_t pizmax = GetZaxis()->GetLast();   
  Int_t pivmin = GetVaxis()->GetFirst();
  Int_t pivmax = GetVaxis()->GetLast();   
  
  GetXaxis()->SetRange(ixmin,ixmax);
  GetYaxis()->SetRange(iymin,iymax);
  GetZaxis()->SetRange(izmin,izmax);
  GetVaxis()->SetRange(ivmin,ivmax);
  
  // exclude underflow/overflow by forcing the axis range bit
  // due to limitation of TAxis::SetRange cannot select only underflow or overflow cannot have underflow or overflow only   
  if (ixmin == 1 && ixmax ==  GetNbinsX() ) GetXaxis()->SetBit(TAxis::kAxisRange); 
  if (iymin == 1 && iymax ==  GetNbinsY() ) GetYaxis()->SetBit(TAxis::kAxisRange); 
  if (izmin == 1 && izmax ==  GetNbinsZ() ) GetZaxis()->SetBit(TAxis::kAxisRange); 
  if (ivmin == 1 && ivmax ==  GetNbinsV() ) GetVaxis()->SetBit(TAxis::kAxisRange); 
  // I can use the useUF or useOF flag for exclude/include the underflow/overflow separately
  // (-1 in the imax is considered as an overflow)
  Bool_t useUF = (ixmin == 0 || iymin == 0 || izmin == 0 || ivmin == 0); 
  Bool_t useOF = (ixmax < 0) || (ixmax > GetNbinsX()) || (iymax < 0) || (iymax > GetNbinsY()) ||
                 (izmax < 0) || (izmax > GetNbinsZ()) || (ivmax < 0) || (ivmax > GetNbinsV()); 
		  
  Bool_t computeErrors = GetSumw2N();
  if (opt.Contains("e") ) { 
    computeErrors = kTRUE;
    opt.Remove(opt.First("e"),1);
  }
  Bool_t originalRange = kFALSE;
  if (opt.Contains('o') ) { 
    originalRange = kTRUE; 
    opt.Remove(opt.First("o"),1);
  }
  
  TH1D * h1 = DoProject1D(name, GetTitle(), this->GetUaxis(), computeErrors, originalRange, useUF, useOF);
  
  // restore original range
  GetXaxis()->SetRange(pixmin,pixmax);
  GetYaxis()->SetRange(piymin,piymax);
  GetZaxis()->SetRange(pizmin,pizmax);
  GetVaxis()->SetRange(pivmin,pivmax);
  
  // draw in current pad 
  if (h1 && opt.Contains("d")) {
    opt.Remove(opt.First("d"),1);
    TVirtualPad *padsav = gPad;
    TVirtualPad *pad = gROOT->GetSelectedPad();
    if (pad) pad->cd();
    if (!gPad->FindObject(h1)) {
      h1->Draw(opt);
    } else {
      h1->Paint(opt);
    }
    if (padsav) padsav->cd();
  }
  
  return h1;
}

//______________________________________________________________________________
TH1D *TH5::ProjectionV(const char *name, Int_t ixmin, Int_t ixmax, Int_t iymin, Int_t iymax,
		       Int_t izmin, Int_t izmax, Int_t iumin, Int_t iumax, Option_t *option) const
{
  //*-*-*-*-*Project a 5-D histogram into a 1-D histogram along V*-*-*-*-*-*-*
  //*-*      ====================================================
  //
  //   The projection is always of the type TH1D.
  //   The projection is made from the cells along the V axis
  //   ranging from ixmin to ixmax, iymin to iymax,... included.
  //   By default, underflow and overflows are included
  //   By Setting ixmin=1 and ixmax=NbinsX the underflow and/or overflow will be excluded
  //
  //   if option "e" is specified, the errors are computed.
  //   if option "d" is specified, the projection is drawn in the current pad.
  //   if option "o" original axis range of the target axes will be
  //   kept, but only bins inside the selected range will be filled.
  //
  //   NOTE that if a TH1D named "name" exists in the current directory or pad 
  //   the histogram is reset and filled again with the projected contents of the TH5.
  //
  //  implemented using Project5D
  
  
  TString opt = option;
  opt.ToLower();
  
  Int_t pixmin = GetXaxis()->GetFirst();
  Int_t pixmax = GetXaxis()->GetLast();
  Int_t piymin = GetYaxis()->GetFirst();
  Int_t piymax = GetYaxis()->GetLast();   
  Int_t pizmin = GetZaxis()->GetFirst();
  Int_t pizmax = GetZaxis()->GetLast();   
  Int_t piumin = GetUaxis()->GetFirst();
  Int_t piumax = GetUaxis()->GetLast();   
  
  GetXaxis()->SetRange(ixmin,ixmax);
  GetYaxis()->SetRange(iymin,iymax);
  GetZaxis()->SetRange(izmin,izmax);
  GetUaxis()->SetRange(iumin,iumax);
  
  // exclude underflow/overflow by forcing the axis range bit
  // due to limitation of TAxis::SetRange cannot select only underflow or overflow cannot have underflow or overflow only   
  if (ixmin == 1 && ixmax ==  GetNbinsX() ) GetXaxis()->SetBit(TAxis::kAxisRange); 
  if (iymin == 1 && iymax ==  GetNbinsY() ) GetYaxis()->SetBit(TAxis::kAxisRange); 
  if (izmin == 1 && izmax ==  GetNbinsZ() ) GetZaxis()->SetBit(TAxis::kAxisRange); 
  if (iumin == 1 && iumax ==  GetNbinsU() ) GetUaxis()->SetBit(TAxis::kAxisRange); 
  // I can use the useUF or useOF flag for exclude/include the underflow/overflow separately
  // (-1 in the imax is considered as an overflow)
  Bool_t useUF = (ixmin == 0 || iymin == 0 || izmin == 0 || iumin == 0); 
  Bool_t useOF = (ixmax < 0) || (ixmax > GetNbinsX()) || (iymax < 0) || (iymax > GetNbinsY()) ||
                 (izmax < 0) || (izmax > GetNbinsZ()) || (iumax < 0) || (iumax > GetNbinsU()); 
		  
  Bool_t computeErrors = GetSumw2N();
  if (opt.Contains("e") ) { 
    computeErrors = kTRUE;
    opt.Remove(opt.First("e"),1);
  }
  Bool_t originalRange = kFALSE;
  if (opt.Contains('o') ) { 
    originalRange = kTRUE; 
    opt.Remove(opt.First("o"),1);
  }
  
  TH1D * h1 = DoProject1D(name, GetTitle(), this->GetVaxis(), computeErrors, originalRange, useUF, useOF);
  
  // restore original range
  GetXaxis()->SetRange(pixmin,pixmax);
  GetYaxis()->SetRange(piymin,piymax);
  GetZaxis()->SetRange(pizmin,pizmax);
  GetUaxis()->SetRange(piumin,piumax);
  
  // draw in current pad 
  if (h1 && opt.Contains("d")) {
    opt.Remove(opt.First("d"),1);
    TVirtualPad *padsav = gPad;
    TVirtualPad *pad = gROOT->GetSelectedPad();
    if (pad) pad->cd();
    if (!gPad->FindObject(h1)) {
      h1->Draw(opt);
    } else {
      h1->Paint(opt);
    }
    if (padsav) padsav->cd();
  }
  
  return h1;
}

//______________________________________________________________________________
TH1D *TH5::DoProject1D(const char* name, const char* title, TAxis* projX, 
		       bool computeErrors, bool originalRange, 
		       bool useUF, bool useOF) const
{
  // internal methdod performing the projection to 1D histogram
  // called from TH5::Project5D
  
  // Create the projection histogram
  TH1D *h1 = 0;
  
  // Get range to use as well as bin limits
  Int_t ixmin = projX->GetFirst();
  Int_t ixmax = projX->GetLast();
  if (ixmin == 0 && ixmax == 0) { ixmin = 1; ixmax = projX->GetNbins(); }
  Int_t nx = ixmax-ixmin+1;
  
  // Create the histogram, either reseting a preexisting one 
  TObject *h1obj = gROOT->FindObject(name);
  if (h1obj && h1obj->InheritsFrom(TH1::Class())) {
    if (h1obj->IsA() != TH1D::Class() ) { 
      Error("DoProject1D","Histogram with name %s must be a TH1D and is a %s",name,h1obj->ClassName());
      return 0; 
    }
    h1 = (TH1D*)h1obj;
    // reset histogram and re-set the axis in any case
    h1->Reset();
    const TArrayD *bins = projX->GetXbins();
    if ( originalRange )
      {
	if (bins->fN == 0) {
	  h1->SetBins(projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
	} else {
	  h1->SetBins(projX->GetNbins(),bins->fArray);
	}
      } else {
      if (bins->fN == 0) {
	h1->SetBins(nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
      } else {
	h1->SetBins(nx,&bins->fArray[ixmin-1]);
      }
    }
  }
  
  if (!h1) { 
    const TArrayD *bins = projX->GetXbins();
    if ( originalRange )
      {
	if (bins->fN == 0) {
	  h1 = new TH1D(name,title,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
	} else {
	  h1 = new TH1D(name,title,projX->GetNbins(),bins->fArray);
	}
      } else {
      if (bins->fN == 0) {
	h1 = new TH1D(name,title,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
      } else {
	h1 = new TH1D(name,title,nx,&bins->fArray[ixmin-1]);
      }
    }
  }
  
  // Copy the axis attributes and the axis labels if needed.
  h1->GetXaxis()->ImportAttributes(projX);
  THashList* labels = projX->GetLabels();
  if (labels) {
    TIter iL(labels);
    TObjString* lb;
    Int_t i = 1;
    while ((lb=(TObjString*)iL())) {
      h1->GetXaxis()->SetBinLabel(i,lb->String().Data());
      i++;
    }
  }
  h1->SetLineColor(this->GetLineColor());
  h1->SetFillColor(this->GetFillColor());
  h1->SetMarkerColor(this->GetMarkerColor());
  h1->SetMarkerStyle(this->GetMarkerStyle());
  
  // Activate errors
  if ( computeErrors ) h1->Sumw2();
  
  // Set references to the axis, so that the bucle has no branches.
  TAxis* out1 = 0;
  TAxis* out2 = 0;
  TAxis* out3 = 0;
  TAxis* out4 = 0;
  if (projX == GetXaxis()) {
    out1 = GetYaxis();
    out2 = GetZaxis();
    out3 = GetUaxis();
    out4 = GetVaxis();
  } else if (projX == GetYaxis()) {
    out1 = GetXaxis();
    out2 = GetZaxis();
    out3 = GetUaxis();
    out4 = GetVaxis();
  } else if (projX == GetZaxis()) {
    out1 = GetXaxis();
    out2 = GetYaxis();
    out3 = GetUaxis();
    out4 = GetVaxis();
  } else if (projX == GetUaxis()) {
    out1 = GetXaxis();
    out2 = GetYaxis();
    out3 = GetZaxis();
    out4 = GetVaxis();
  } else {
    out1 = GetXaxis();
    out2 = GetYaxis();
    out3 = GetZaxis();
    out4 = GetUaxis();
  }
  
  Int_t *refX = 0, *refY = 0, *refZ = 0, *refU = 0, *refV = 0;
  Int_t ixbin, out1bin, out2bin, out3bin, out4bin;
  if (projX == GetXaxis()) {refX = &ixbin;   refY = &out1bin; refZ = &out2bin; refU = &out3bin; refV = &out4bin;}
  if (projX == GetYaxis()) {refX = &out1bin; refY = &ixbin;   refZ = &out2bin; refU = &out3bin; refV = &out4bin;}
  if (projX == GetZaxis()) {refX = &out1bin; refY = &out2bin; refZ = &ixbin;   refU = &out3bin; refV = &out4bin;}
  if (projX == GetUaxis()) {refX = &out1bin; refY = &out2bin; refZ = &out3bin; refU = &ixbin;   refV = &out4bin;}
  if (projX == GetVaxis()) {refX = &out1bin; refY = &out2bin; refZ = &out3bin; refU = &out4bin; refV = &ixbin;}
  R__ASSERT (refX != 0 && refY != 0 && refZ != 0 && refU != 0 && refV != 0); 
  
  // Fill the projected histogram excluding underflow/overflows if considered in the option
  // if specified in the option (by default they considered)
  Double_t totcont  = 0;
  
  Int_t out1min = out1->GetFirst(); 
  Int_t out1max = out1->GetLast(); 
  // GetFirst(), GetLast() can return (0,0) when the range bit is set artifically (see TAxis::SetRange)
  if (out1min == 0 && out1max == 0) {out1min = 1; out1max = out1->GetNbins();}
  // correct for underflow/overflows
  if (useUF && !out1->TestBit(TAxis::kAxisRange) )  out1min -= 1; 
  if (useOF && !out1->TestBit(TAxis::kAxisRange) )  out1max += 1; 
  Int_t out2min = out2->GetFirst(); 
  Int_t out2max = out2->GetLast(); 
  if (out2min == 0 && out2max == 0) {out2min = 1; out2max = out2->GetNbins();}
  if (useUF && !out2->TestBit(TAxis::kAxisRange) )  out2min -= 1; 
  if (useOF && !out2->TestBit(TAxis::kAxisRange) )  out2max += 1; 
  Int_t out3min = out3->GetFirst(); 
  Int_t out3max = out3->GetLast(); 
  if (out3min == 0 && out3max == 0) {out3min = 1; out3max = out3->GetNbins();}
  if (useUF && !out3->TestBit(TAxis::kAxisRange) )  out3min -= 1; 
  if (useOF && !out3->TestBit(TAxis::kAxisRange) )  out3max += 1; 
  Int_t out4min = out4->GetFirst(); 
  Int_t out4max = out4->GetLast(); 
  if (out4min == 0 && out4max == 0) {out4min = 1; out4max = out4->GetNbins();}
  if (useUF && !out4->TestBit(TAxis::kAxisRange) )  out4min -= 1; 
  if (useOF && !out4->TestBit(TAxis::kAxisRange) )  out4max += 1; 
  
  for (ixbin=0; ixbin <= 1 + projX->GetNbins(); ixbin++){
    if (projX->TestBit(TAxis::kAxisRange) && (ixbin < ixmin || ixbin > ixmax)) continue;
    
    Double_t cont = 0; 
    Double_t err2 = 0; 
    
    // loop on the bins to be integrated (outbin should be called inbin)
    for (out1bin = out1min; out1bin <= out1max; out1bin++){
      for (out2bin = out2min; out2bin <= out2max; out2bin++){
	for (out3bin = out3min; out3bin <= out3max; out3bin++){
	  for (out4bin = out4min; out4bin <= out4max; out4bin++){
	    
            Int_t bin = GetBin(*refX, *refY, *refZ, *refU, *refV);

            // sum the bin contents and errors if needed
            cont += GetBinContent(bin);
            if (computeErrors) { 
	      Double_t exyzuv = GetBinError(bin);
	      err2 += exyzuv*exyzuv;
            }
	  }
	}
      }
    }
    Int_t ix = h1->FindBin(projX->GetBinCenter(ixbin));
    h1->SetBinContent(ix ,cont);
    if (computeErrors) h1->SetBinError(ix, TMath::Sqrt(err2) ); 
    // sum all content
    totcont += cont;
    
  }
  
  // since we use a combination of fill and SetBinError we need to reset and recalculate the statistics
  // for weighted histograms otherwise sumw2 will be wrong. 
  // We  can keep the original statistics from the TH5 if the projected sumw is consistent with original one
  // i.e. when no events are thrown away  
  bool resetStats = true; 
  double eps = 1.E-12;
  //  if (IsA() == TH5F::Class() ) eps = 1.E-6;
  if (fTsumw != 0 && TMath::Abs( fTsumw - totcont) <  TMath::Abs(fTsumw) * eps) resetStats = false; 
  
  bool resetEntries = resetStats; 
  // entries are calculated using underflow/overflow. If excluded entries must be reset
  resetEntries |= !useUF || !useOF; 
  
  
  if (!resetStats) {
    Double_t stats[kNstat5];
    GetStats(stats); 
    if (projX == GetYaxis()) {
      stats[2] = stats[4];
      stats[3] = stats[5]; 
    } else if (projX == GetZaxis()) {
      stats[2] = stats[7];
      stats[3] = stats[8]; 
    } else if (projX == GetUaxis()) {
      stats[2] = stats[11];
      stats[3] = stats[12]; 
    } else if (projX == GetVaxis()) {
      stats[2] = stats[16];
      stats[3] = stats[17]; 
    }
    h1->PutStats(stats);
  } else {
    // reset statistics 
    h1->ResetStats();
  }
  if (resetEntries) { 
    // in case of error calculation (i.e. when Sumw2() is set) 
    // use the effective entries for the entries
    // since this  is the only way to estimate them
    Double_t entries =  TMath::Floor( totcont + 0.5); // to avoid numerical rounding         
    if (computeErrors) entries = h1->GetEffectiveEntries();
    h1->SetEntries( entries );  
  }
  else { 
    h1->SetEntries( fEntries );  
  }
  
  return h1;
 }

//______________________________________________________________________________
TH2D *TH5::DoProject2D(const char* name, const char * title, TAxis* projX, TAxis* projY,  
                    bool computeErrors, bool originalRange,
                    bool useUF, bool useOF) const
{
  // internal method performing the projection to a 2D histogram
  // called from TH5::Project5D 
  
  TH2D *h2 = 0;
  
  // Get range to use as well as bin limits
  Int_t ixmin = projX->GetFirst();
  Int_t ixmax = projX->GetLast();
  Int_t iymin = projY->GetFirst();
  Int_t iymax = projY->GetLast();
  if (ixmin == 0 && ixmax == 0) { ixmin = 1; ixmax = projX->GetNbins(); }
  if (iymin == 0 && iymax == 0) { iymin = 1; iymax = projY->GetNbins(); }
  Int_t nx = ixmax-ixmin+1;
  Int_t ny = iymax-iymin+1;
  
  // Create the histogram, either reseting a preexisting one 
  //  or creating one from scratch.
  // Does an object with the same name exists?
  TObject *h2obj = gROOT->FindObject(name);
  if (h2obj && h2obj->InheritsFrom(TH1::Class())) {
    if ( h2obj->IsA() != TH2D::Class() ) { 
      Error("DoProject2D","Histogram with name %s must be a TH2D and is a %s",name,h2obj->ClassName());
      return 0; 
    }
    h2 = (TH2D*)h2obj;
    // reset histogram and its axes
    h2->Reset();
    const TArrayD *xbins = projX->GetXbins();
    const TArrayD *ybins = projY->GetXbins();
    if ( originalRange ) {
      h2->SetBins(projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
		  ,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
      // set bins for mixed axis do not exists - need to set afterwards the variable bins
      if (ybins->fN != 0) 
	h2->GetXaxis()->Set(projY->GetNbins(),&ybins->fArray[iymin-1]);
      if (xbins->fN != 0)
	h2->GetYaxis()->Set(projX->GetNbins(),&xbins->fArray[ixmin-1]);
    } else {
      h2->SetBins(ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
		  ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
      if (ybins->fN != 0) 
	h2->GetXaxis()->Set(ny,&ybins->fArray[iymin-1]);
      if (xbins->fN != 0) 
	h2->GetYaxis()->Set(nx,&xbins->fArray[ixmin-1]);
    }
   }
  
  
  if (!h2) { 
    const TArrayD *xbins = projX->GetXbins();
    const TArrayD *ybins = projY->GetXbins();
    if ( originalRange )
      {
	if (xbins->fN == 0 && ybins->fN == 0) {
	  h2 = new TH2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
			,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
	} else if (ybins->fN == 0) {
	  h2 = new TH2D(name,title,projY->GetNbins(),projY->GetXmin(),projY->GetXmax()
			,projX->GetNbins(),&xbins->fArray[ixmin-1]);
	} else if (xbins->fN == 0) {
	  h2 = new TH2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1]
			,projX->GetNbins(),projX->GetXmin(),projX->GetXmax());
	} else {
	  h2 = new TH2D(name,title,projY->GetNbins(),&ybins->fArray[iymin-1],projX->GetNbins(),&xbins->fArray[ixmin-1]);
	}
      } else {
      if (xbins->fN == 0 && ybins->fN == 0) {
	h2 = new TH2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
		      ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
      } else if (ybins->fN == 0) {
	h2 = new TH2D(name,title,ny,projY->GetBinLowEdge(iymin),projY->GetBinUpEdge(iymax)
		      ,nx,&xbins->fArray[ixmin-1]);
      } else if (xbins->fN == 0) {
	h2 = new TH2D(name,title,ny,&ybins->fArray[iymin-1]
		      ,nx,projX->GetBinLowEdge(ixmin),projX->GetBinUpEdge(ixmax));
      } else {
	h2 = new TH2D(name,title,ny,&ybins->fArray[iymin-1],nx,&xbins->fArray[ixmin-1]);
      }
    }
  }

  // Copy the axis attributes and the axis labels if needed.
  THashList* labels1 = 0;
  THashList* labels2 = 0;
  // "xy"
  h2->GetXaxis()->ImportAttributes(projY);
  h2->GetYaxis()->ImportAttributes(projX);
  labels1 = projY->GetLabels();
  labels2 = projX->GetLabels();
  if (labels1) {
    TIter iL(labels1);
    TObjString* lb;
    Int_t i = 1;
    while ((lb=(TObjString*)iL())) {
      h2->GetXaxis()->SetBinLabel(i,lb->String().Data());
      i++;
    }
  }
  if (labels2) {
    TIter iL(labels2);
    TObjString* lb;
    Int_t i = 1;
    while ((lb=(TObjString*)iL())) {
      h2->GetYaxis()->SetBinLabel(i,lb->String().Data());
      i++;
    }
  }
  h2->SetLineColor(this->GetLineColor());
  h2->SetFillColor(this->GetFillColor());
  h2->SetMarkerColor(this->GetMarkerColor());
  h2->SetMarkerStyle(this->GetMarkerStyle());
  
  // Activate errors
  if ( computeErrors) h2->Sumw2();
  
  // Set references to the axis, so that the bucle has no branches.
  TAxis *out1 = 0, *out2 = 0, *out3 = 0;
  Int_t *refX = 0, *refY = 0, *refZ = 0, *refU = 0, *refV = 0;
  Int_t ixbin, iybin, out1bin, out2bin, out3bin;
  if (projX == GetXaxis() && projY == GetYaxis()) {
    out1 = GetZaxis(); out2 = GetUaxis(); out3 = GetVaxis();
    refX = &ixbin; refY = &iybin; refZ = &out1bin; refU = &out2bin; refV = &out3bin;
  }
  if (projX == GetYaxis() && projY == GetXaxis()) {
    out1 = GetZaxis(); out2 = GetUaxis(); out3 = GetVaxis();
    refX = &iybin; refY = &ixbin; refZ = &out1bin; refU = &out2bin; refV = &out3bin;
  }
  if (projX == GetXaxis() && projY == GetZaxis()) {
    out1 = GetYaxis(); out2 = GetUaxis(); out3 = GetVaxis();
    refX = &ixbin; refY = &out1bin; refZ = &iybin; refU = &out2bin; refV = &out3bin;
  }
  if (projX == GetZaxis() && projY == GetXaxis()) {
    out1 = GetYaxis(); out2 = GetUaxis(); out3 = GetVaxis();
    refX = &iybin; refY = &out1bin; refZ = &ixbin; refU = &out2bin; refV = &out3bin;
  }
  if (projX == GetXaxis() && projY == GetUaxis()) {
    out1 = GetYaxis(); out2 = GetZaxis(); out3 = GetVaxis();
    refX = &ixbin; refY = &out1bin; refZ = &out2bin; refU = &iybin; refV = &out3bin;
  }
  if (projX == GetUaxis() && projY == GetXaxis()) {
    out1 = GetYaxis(); out2 = GetZaxis(); out3 = GetVaxis();
    refX = &iybin; refY = &out1bin; refZ = &out2bin; refU = &ixbin; refV = &out3bin;
  }
  if (projX == GetXaxis() && projY == GetVaxis()) {
    out1 = GetYaxis(); out2 = GetZaxis(); out3 = GetUaxis();
    refX = &ixbin; refY = &out1bin; refZ = &out2bin; refU = &out3bin; refV = &iybin;
  }
  if (projX == GetVaxis() && projY == GetXaxis()) {
    out1 = GetYaxis(); out2 = GetZaxis(); out3 = GetUaxis();
    refX = &iybin; refY = &out1bin; refZ = &out2bin; refU = &out3bin; refV = &ixbin;
  }
  if (projX == GetYaxis() && projY == GetZaxis()) {
    out1 = GetXaxis(); out2 = GetUaxis(); out3 = GetVaxis();
    refX = &out1bin; refY = &ixbin; refZ = &iybin; refU = &out2bin; refV = &out3bin;
  }
  if (projX == GetZaxis() && projY == GetYaxis()) {
    out1 = GetXaxis(); out2 = GetUaxis(); out3 = GetVaxis();
    refX = &out1bin; refY = &iybin; refZ = &ixbin; refU = &out2bin; refV = &out3bin;
  }
  if (projX == GetYaxis() && projY == GetUaxis()) {
    out1 = GetXaxis(); out2 = GetZaxis(); out3 = GetVaxis();
    refX = &out1bin; refY = &ixbin; refZ = &out2bin; refU = &iybin; refV = &out3bin;
  }
  if (projX == GetUaxis() && projY == GetYaxis()) {
    out1 = GetXaxis(); out2 = GetZaxis(); out3 = GetVaxis();
    refX = &out1bin; refY = &iybin; refZ = &out2bin; refU = &ixbin; refV = &out3bin;
  }
  if (projX == GetYaxis() && projY == GetVaxis()) {
    out1 = GetXaxis(); out2 = GetZaxis(); out3 = GetUaxis();
    refX = &out1bin; refY = &ixbin; refZ = &out2bin; refU = &out3bin; refV = &iybin;
  }
  if (projX == GetVaxis() && projY == GetYaxis()) {
    out1 = GetXaxis(); out2 = GetZaxis(); out3 = GetUaxis();
    refX = &out1bin; refY = &iybin; refZ = &out2bin; refU = &out3bin; refV = &ixbin;
  }
  if (projX == GetZaxis() && projY == GetUaxis()) {
    out1 = GetXaxis(); out2 = GetYaxis(); out3 = GetVaxis();
    refX = &out1bin; refY = &out2bin; refZ = &ixbin; refU = &iybin; refV = &out3bin;
  }
  if (projX == GetUaxis() && projY == GetZaxis()) {
    out1 = GetXaxis(); out2 = GetYaxis(); out3 = GetVaxis();
    refX = &out1bin; refY = &out2bin; refZ = &iybin; refU = &ixbin; refV = &out3bin;
  }
  if (projX == GetZaxis() && projY == GetVaxis()) {
    out1 = GetXaxis(); out2 = GetYaxis(); out3 = GetUaxis();
    refX = &out1bin; refY = &out2bin; refZ = &ixbin; refU = &out3bin; refV = &iybin;
  }
  if (projX == GetVaxis() && projY == GetZaxis()) {
    out1 = GetXaxis(); out2 = GetYaxis(); out3 = GetUaxis();
    refX = &out1bin; refY = &out2bin; refZ = &iybin; refU = &out3bin; refV = &ixbin;
  }
  if (projX == GetUaxis() && projY == GetVaxis()) {
    out1 = GetXaxis(); out2 = GetYaxis(); out3 = GetZaxis();
    refX = &out1bin; refY = &out2bin; refZ = &out3bin; refU = &ixbin; refV = &iybin;
  }
  if (projX == GetVaxis() && projY == GetUaxis()) {
    out1 = GetXaxis(); out2 = GetYaxis(); out3 = GetZaxis();
    refX = &out1bin; refY = &out2bin; refZ = &out3bin; refU = &iybin; refV = &ixbin;
  }
  R__ASSERT (refX != 0 && refY != 0 && refZ != 0 && refU != 0 && refV != 0); 

  // Fill the projected histogram excluding underflow/overflows if considered in the option
  // if specified in the option (by default they considered)
  Double_t totcont  = 0;
  
  Int_t out1min = out1->GetFirst(); 
  Int_t out1max = out1->GetLast(); 
  // GetFirst(), GetLast() can return (0,0) when the range bit is set artifically (see TAxis::SetRange)
  if (out1min == 0 && out1max == 0) {out1min = 1; out1max = out1->GetNbins();}
  // correct for underflow/overflows
  if (useUF && !out1->TestBit(TAxis::kAxisRange) )  out1min -= 1; 
  if (useOF && !out1->TestBit(TAxis::kAxisRange) )  out1max += 1; 
  Int_t out2min = out2->GetFirst(); 
  Int_t out2max = out2->GetLast(); 
  if (out2min == 0 && out2max == 0) {out2min = 1; out2max = out2->GetNbins();}
  if (useUF && !out2->TestBit(TAxis::kAxisRange) )  out2min -= 1; 
  if (useOF && !out2->TestBit(TAxis::kAxisRange) )  out2max += 1; 
  Int_t out3min = out3->GetFirst(); 
  Int_t out3max = out3->GetLast(); 
  if (out3min == 0 && out3max == 0) {out3min = 1; out3max = out3->GetNbins();}
  if (useUF && !out3->TestBit(TAxis::kAxisRange) )  out3min -= 1; 
  if (useOF && !out3->TestBit(TAxis::kAxisRange) )  out3max += 1; 

  for (ixbin=0; ixbin <= 1 + projX->GetNbins(); ixbin++){
    if (projX->TestBit(TAxis::kAxisRange) && (ixbin < ixmin || ixbin > ixmax)) continue;
    Int_t ix = h2->GetYaxis()->FindBin(projX->GetBinCenter(ixbin));
    
    for (iybin=0;iybin<=1+projY->GetNbins();iybin++){
      if ( projY->TestBit(TAxis::kAxisRange) && ( iybin < iymin || iybin > iymax )) continue;
      Int_t iy = h2->GetXaxis()->FindBin(projY->GetBinCenter(iybin));

      Double_t cont = 0; 
      Double_t err2 = 0;
      
      // loop on the bins to be integrated (outbin should be called inbin)
      for (out1bin = out1min; out1bin <= out1max; out1bin++){
	for (out2bin = out2min; out2bin <= out2max; out2bin++){
	  for (out3bin = out3min; out3bin <= out3max; out3bin++){
	
	    Int_t bin = GetBin(*refX, *refY, *refZ, *refU, *refV);
	    
	    // sum the bin contents and errors if needed
	    cont += GetBinContent(bin);
	    if (computeErrors) { 
	      Double_t exyzuv = GetBinError(bin);
	      err2 += exyzuv*exyzuv;
	    }
	  }
	}
      }
      
      // remember axis are inverted 
      h2->SetBinContent(iy , ix, cont);
      if (computeErrors) h2->SetBinError(iy, ix, TMath::Sqrt(err2) ); 
      // sum all content
      totcont += cont;
      
    }
  }
  
  // since we use fill we need to reset and recalculate the statistics (see comment in DoProject1D )
  // or keep original statistics if consistent sumw2  
  bool resetStats = true; 
  double eps = 1.E-12;
  if (fTsumw != 0 && TMath::Abs(fTsumw - totcont) <  TMath::Abs(fTsumw) * eps) resetStats = false; 

  bool resetEntries = resetStats;
  // entries are calculated using underflow/overflow. If excluded entries must be reset
  resetEntries |= !useUF || !useOF; 
  
  if (!resetStats) {
    Double_t stats[kNstat5];
    Double_t oldst[kNstat5]; // old statistics
    for (Int_t i = 0; i < kNstat5; ++i) {oldst[i] = 0;}
    GetStats(oldst); 
    std::copy(oldst,oldst+kNstat5,stats);
    // note that projX refer to Y axis and projX refer to the X axis of projected histogram
    // nothing to do for projection in Y vs X
    if (projY == GetXaxis()) {
      if (projX == GetZaxis()) {  // case XZ
	stats[4] = oldst[7];
	stats[5] = oldst[8];
	stats[6] = oldst[9];
      } else if (projX == GetUaxis()) {  // case XU
	stats[4] = oldst[11];
	stats[5] = oldst[12];
	stats[6] = oldst[13];
      } else if (projX == GetVaxis()) {  // case XV
	stats[4] = oldst[16];
	stats[5] = oldst[17];
	stats[6] = oldst[18];
      }
    } else if (projY == GetYaxis()) {
      stats[2] = oldst[4];
      stats[3] = oldst[5];
      if (projX == GetXaxis()) {  // case YX
	stats[4] = oldst[2]; 
	stats[5] = oldst[3];
      } else if (projX == GetZaxis()) {  // case YZ
	stats[4] = oldst[7]; 
	stats[5] = oldst[8];
	stats[6] = oldst[10];
      } else if (projX == GetUaxis()) {  // case YU
	stats[4] = oldst[11];
	stats[5] = oldst[12];
	stats[6] = oldst[14];
      } else if (projX == GetVaxis()) {  // case YV
	stats[4] = oldst[16];
	stats[5] = oldst[17];
	stats[6] = oldst[19];
      }
    } else if (projY == GetZaxis()) {
      stats[2] = oldst[7];
      stats[3] = oldst[8]; 
      if (projX == GetXaxis()) {  // case ZX
	stats[4] = oldst[2]; 
	stats[5] = oldst[3];
	stats[6] = oldst[9];
      } else if (projX == GetYaxis())  { // case ZY
	stats[4] = oldst[4]; 
	stats[5] = oldst[5];
	stats[6] = oldst[10];
      } else if (projX == GetUaxis())  { // case ZU
	stats[4] = oldst[11]; 
	stats[5] = oldst[12];
	stats[6] = oldst[15];
      } else if (projX == GetVaxis())  { // case ZV
	stats[4] = oldst[16]; 
	stats[5] = oldst[17];
	stats[6] = oldst[20];
      }
    } else if (projY == GetUaxis()) {
      stats[2] = oldst[11];
      stats[3] = oldst[12]; 
      if (projX == GetXaxis()) {  // case UX
	stats[4] = oldst[2]; 
	stats[5] = oldst[3];
	stats[6] = oldst[13];
      } else if (projX == GetYaxis())  { // case UY
	stats[4] = oldst[4]; 
	stats[5] = oldst[5];
	stats[6] = oldst[14];
      } else if (projX == GetZaxis())  { // case UZ
	stats[4] = oldst[7]; 
	stats[5] = oldst[8];
	stats[6] = oldst[15];
      } else if (projX == GetVaxis())  { // case UV
	stats[4] = oldst[16]; 
	stats[5] = oldst[17];
	stats[6] = oldst[21];
      }
    } else if (projY == GetVaxis()) {
      stats[2] = oldst[16];
      stats[3] = oldst[17]; 
      if (projX == GetXaxis()) {  // case VX
	stats[4] = oldst[2]; 
	stats[5] = oldst[3];
	stats[6] = oldst[18];
      } else if (projX == GetYaxis())  { // case VY
	stats[4] = oldst[4]; 
	stats[5] = oldst[5];
	stats[6] = oldst[19];
      } else if (projX == GetZaxis())  { // case VZ
	stats[4] = oldst[7]; 
	stats[5] = oldst[8];
	stats[6] = oldst[20];
      } else if (projX == GetUaxis())  { // case VU
	stats[4] = oldst[11]; 
	stats[5] = oldst[12];
	stats[6] = oldst[21];
      }
    }
    // set the new statistics 
    h2->PutStats(stats);
  } else { 
    // recalculate the statistics
    h2->ResetStats(); 
  }

  if (resetEntries) { 
    // use the effective entries for the entries
    // since this  is the only way to estimate them
    Double_t entries =  h2->GetEffectiveEntries();
    if (!computeErrors) entries = TMath::Floor( entries + 0.5); // to avoid numerical rounding
    h2->SetEntries( entries );  
  }
  else { 
    h2->SetEntries( fEntries );  
  }
  
  return h2;
 }


//______________________________________________________________________________
TH1 *TH5::Project5D(Option_t *option) const
{
  // Project a 5-d histogram into 1 or 2-d histograms depending on the
  // option parameter
  // option may contain a combination of the characters x, y, z, u, v, e
  // option = "x" return the x projection into a TH1D histogram
  // option = "y" return the y projection into a TH1D histogram
  // ...
  // option = "xy" return the x versus y projection into a TH2D histogram
  // option = "yx" return the y versus x projection into a TH2D histogram
  // ...
  // NB: the notation "a vs b" means "a" vertical and "b" horizontal
  //
  // option = "o" original axis range of the target axes will be
  //   kept, but only bins inside the selected range will be filled.
  //
  // If option contains the string "e", errors are computed
  //
  // The projection is made for the selected bins only.
  // To select a bin range along an axis, use TAxis::SetRange, eg
  //    h3.GetYaxis()->SetRange(23,56);
  //
  // NOTE 1: The generated histogram is named th3name + option
  // eg if the TH5* h histogram is named "myhist", then
  // h->Project5D("xy"); produces a TH2D histogram named "myhist_xy"
  // if a histogram of the same type already exists, it is overwritten.
  // The following sequence
  //    h->Project5D("xy");
  //    h->Project5D("xy2");
  //  will generate two TH2D histograms named "myhist_xy" and "myhist_xy2"
  //  A different name can be generated by attaching a string to the option
  //  For example h->Project5D("name_xy") will generate an histogram with the name:  h3dname_name_xy. 
  //
  //  NOTE 2: If an histogram of the same type already exists, 
  //  the histogram is reset and filled again with the projected contents of the TH5.
  //
  //  NOTE 3: The number of entries in the projected histogram is estimated from the number of 
  //  effective entries for all the cells included in the projection. 
  //
  //  NOTE 4: underflow/overflow are included by default in the projection 
  //  To exclude underflow and/or overflow (for both axis in case of a projection to a 1D histogram) use option "NUF" and/or "NOF"
  //  With SetRange() you can have all bins except underflow/overflow only if you set the axis bit range as 
  //  following after having called SetRange: 
  //    axis->SetRange(1, axis->GetNbins());
  //    axis->SetBit(TAxis::kAxisRange);  
  //          
  
  TString opt = option; opt.ToLower();
  Int_t pcase = 0;
  TString ptype; 
  if (opt.Contains("x"))  { pcase = 1; ptype = "x"; }
  if (opt.Contains("y"))  { pcase = 2; ptype = "y"; }
  if (opt.Contains("z"))  { pcase = 3; ptype = "z"; }
  if (opt.Contains("u"))  { pcase = 4; ptype = "u"; }
  if (opt.Contains("v"))  { pcase = 5; ptype = "v"; }
  if (opt.Contains("xy")) { pcase = 6; ptype = "xy"; }
  if (opt.Contains("yx")) { pcase = 7; ptype = "yx"; }
  if (opt.Contains("xz")) { pcase = 8; ptype = "xz"; }
  if (opt.Contains("zx")) { pcase = 9; ptype = "zx"; }
  if (opt.Contains("xu")) { pcase = 10; ptype = "xu"; }
  if (opt.Contains("ux")) { pcase = 11; ptype = "ux"; }
  if (opt.Contains("xv")) { pcase = 12; ptype = "xv"; }
  if (opt.Contains("vx")) { pcase = 13; ptype = "vx"; }
  if (opt.Contains("yz")) { pcase = 14; ptype = "yz"; }
  if (opt.Contains("zy")) { pcase = 15; ptype = "zy"; }
  if (opt.Contains("yu")) { pcase = 16; ptype = "yu"; }
  if (opt.Contains("uy")) { pcase = 17; ptype = "uy"; }
  if (opt.Contains("yv")) { pcase = 18; ptype = "yv"; }
  if (opt.Contains("vy")) { pcase = 19; ptype = "vy"; }
  if (opt.Contains("zu")) { pcase = 20; ptype = "zu"; }
  if (opt.Contains("uz")) { pcase = 21; ptype = "uz"; }
  if (opt.Contains("zv")) { pcase = 22; ptype = "zv"; }
  if (opt.Contains("vz")) { pcase = 23; ptype = "vz"; }
  if (opt.Contains("uv")) { pcase = 24; ptype = "uv"; }
  if (opt.Contains("vu")) { pcase = 25; ptype = "vu"; }
   
  if (pcase == 0) { 
    Error("Project5D","No projection axis specified - return a NULL pointer"); 
    return 0; 
  }
  // do not remove ptype from opt to use later in the projected histo name
  
  Bool_t computeErrors = GetSumw2N();
  if (opt.Contains("e") ) { 
    computeErrors = kTRUE;
    opt.Remove(opt.First("e"),1);
  }
  
  Bool_t useUF = kTRUE;
  Bool_t useOF = kTRUE;
  if (opt.Contains("nuf") ) { 
    useUF = kFALSE;
    opt.Remove(opt.Index("nuf"),3);
  }
  if (opt.Contains("nof") ) { 
    useOF = kFALSE;
    opt.Remove(opt.Index("nof"),3);
  }
  
  Bool_t originalRange = kFALSE;
  if (opt.Contains('o') ) { 
    originalRange = kTRUE; 
    opt.Remove(opt.First("o"),1);
  }

  
  // Create the projection histogram
  TH1 *h = 0;
  
  TString name = GetName();
  TString title = GetTitle();
  name  += "_";   name  += opt;  // opt may include a user defined name
  title += " ";   title += ptype; title += " projection";   
  
  switch (pcase) {
  case 1: 
    // "x"
    h = DoProject1D(name, title, this->GetXaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 2:
    // "y"
    h = DoProject1D(name, title, this->GetYaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 3:
    // "z"
    h = DoProject1D(name, title, this->GetZaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 4:
    // "u"
    h = DoProject1D(name, title, this->GetUaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 5:
    // "v"
    h = DoProject1D(name, title, this->GetVaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 6:
    // "xy"
    h = DoProject2D(name, title, this->GetXaxis(),this->GetYaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 7:
    // "yx"
    h = DoProject2D(name, title, this->GetYaxis(),this->GetXaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 8:
    // "xz"
    h = DoProject2D(name, title, this->GetXaxis(),this->GetZaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 9:
    // "zx"
    h = DoProject2D(name, title, this->GetZaxis(),this->GetXaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 10:
    // "xu"
    h = DoProject2D(name, title, this->GetXaxis(),this->GetUaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 11:
    // "ux"
    h = DoProject2D(name, title, this->GetUaxis(),this->GetXaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 12:
    // "xv"
    h = DoProject2D(name, title, this->GetXaxis(),this->GetVaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 13:
    // "vx"
    h = DoProject2D(name, title, this->GetVaxis(),this->GetXaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 14:
    // "yz"
    h = DoProject2D(name, title, this->GetYaxis(),this->GetZaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 15:
    // "zy"
    h = DoProject2D(name, title, this->GetZaxis(),this->GetYaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 16:
    // "yu"
    h = DoProject2D(name, title, this->GetYaxis(),this->GetUaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 17:
    // "uy"
    h = DoProject2D(name, title, this->GetUaxis(),this->GetYaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 18:
    // "yv"
    h = DoProject2D(name, title, this->GetYaxis(),this->GetVaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 19:
    // "vy"
    h = DoProject2D(name, title, this->GetVaxis(),this->GetYaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 20:
    // "zu"
    h = DoProject2D(name, title, this->GetZaxis(),this->GetUaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 21:
    // "uz"
    h = DoProject2D(name, title, this->GetUaxis(),this->GetZaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 22:
    // "zv"
    h = DoProject2D(name, title, this->GetZaxis(),this->GetVaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 23:
    // "vz"
    h = DoProject2D(name, title, this->GetVaxis(),this->GetZaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 24:
    // "uv"
    h = DoProject2D(name, title, this->GetUaxis(),this->GetVaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  case 25:
    // "vu"
    h = DoProject2D(name, title, this->GetVaxis(),this->GetUaxis(), computeErrors, originalRange, useUF, useOF);
    break;
  }

  // draw in current pad 
  if (h && opt.Contains("d")) {
    opt.Remove(opt.First("d"),1);
    TVirtualPad *padsav = gPad;
    TVirtualPad *pad = gROOT->GetSelectedPad();
    if (pad) pad->cd();
    if (!gPad->FindObject(h)) {
      h->Draw(opt);
    } else {
      h->Paint(opt);
    }
    if (padsav) padsav->cd();
  }
  
  return h;
}

//______________________________________________________________________________
void TH5::PutStats(Double_t *stats)
{
  // Replace current statistics with the values in array stats
  
  TH1::PutStats(stats);
  fTsumwy = stats[4];
  fTsumwy2 = stats[5];
  fTsumwxy = stats[6];
  fTsumwz = stats[7];
  fTsumwz2 = stats[8];
  fTsumwxz = stats[9];
  fTsumwyz = stats[10];
  fTsumwu = stats[11];
  fTsumwu2 = stats[12];
  fTsumwxu = stats[13];
  fTsumwyu = stats[14];
  fTsumwzu = stats[15];
  fTsumwv = stats[16];
  fTsumwv2 = stats[17];
  fTsumwxv = stats[18];
  fTsumwyv = stats[19];
  fTsumwzv = stats[20];
  fTsumwuv = stats[21];
}

//______________________________________________________________________________
// TH5 *TH5::RebinX(Int_t ngroup, const char *newname)
// {
//   // Rebin only the X axis
//   // see Rebin3D
//   return Rebin3D(ngroup, 1, 1, newname);
// }

// //______________________________________________________________________________
// TH5 *TH5::RebinY(Int_t ngroup, const char *newname)
// {
//   // Rebin only the Y axis
//   // see Rebin3D
//   return Rebin3D(1, ngroup, 1, newname);
// }

// //______________________________________________________________________________
// TH5 *TH5::RebinZ(Int_t ngroup, const char *newname)
// {
//   // Rebin only the Z axis
//   // see Rebin3D
//   return Rebin3D(1, 1, ngroup, newname);
  
// }

// //______________________________________________________________________________
// TH5 *TH5::Rebin3D(Int_t nxgroup, Int_t nygroup, Int_t nzgroup, const char *newname)
// {
//    //   -*-*-*Rebin this histogram grouping nxgroup/nygroup/nzgroup bins along the xaxis/yaxis/zaxis together*-*-*-*-
//    //         =================================================================================
//    //   if newname is not blank a new temporary histogram hnew is created.
//    //   else the current histogram is modified (default)
//    //   The parameter nxgroup/nygroup indicate how many bins along the xaxis/yaxis of this
//    //   have to me merged into one bin of hnew
//    //   If the original histogram has errors stored (via Sumw2), the resulting
//    //   histograms has new errors correctly calculated.
//    //
//    //   examples: if hpxpy is an existing TH5 histogram with 40 x 40 x 40 bins
//    //     hpxpypz->Rebin3D();  // merges two bins along the xaxis and yaxis in one in hpxpypz
//    //                          // Carefull: previous contents of hpxpy are lost
//    //     hpxpypz->RebinX(5);  //merges five bins along the xaxis in one in hpxpypz
//    //     TH5 *hnew = hpxpypz->RebinY(5,"hnew"); // creates a new histogram hnew
//    //                                          // merging 5 bins of h1 along the yaxis in one bin
//    //
//    //   NOTE : If nxgroup/nygroup is not an exact divider of the number of bins,
//    //          along the xaxis/yaxis the top limit(s) of the rebinned histogram
//    //          is changed to the upper edge of the xbin=newxbins*nxgroup resp.
//    //          ybin=newybins*nygroup and the corresponding bins are added to
//    //          the overflow bin.
//    //          Statistics will be recomputed from the new bin contents.

//    Int_t i,j,k,xbin,ybin,zbin;
//    Int_t nxbins  = fXaxis.GetNbins();
//    Int_t nybins  = fYaxis.GetNbins();
//    Int_t nzbins  = fZaxis.GetNbins();
//    Double_t xmin  = fXaxis.GetXmin();
//    Double_t xmax  = fXaxis.GetXmax();
//    Double_t ymin  = fYaxis.GetXmin();
//    Double_t ymax  = fYaxis.GetXmax();
//    Double_t zmin  = fZaxis.GetXmin();
//    Double_t zmax  = fZaxis.GetXmax();
//    if ((nxgroup <= 0) || (nxgroup > nxbins)) {
//       Error("Rebin", "Illegal value of nxgroup=%d",nxgroup);
//       return 0;
//    }
//    if ((nygroup <= 0) || (nygroup > nybins)) {
//       Error("Rebin", "Illegal value of nygroup=%d",nygroup);
//       return 0;
//    }
//    if ((nzgroup <= 0) || (nzgroup > nzbins)) {
//       Error("Rebin", "Illegal value of nzgroup=%d",nzgroup);
//       return 0;
//    }

//    Int_t newxbins = nxbins/nxgroup;
//    Int_t newybins = nybins/nygroup;
//    Int_t newzbins = nzbins/nzgroup;

//    // Save old bin contents into a new array
//    Double_t entries = fEntries;
//    Double_t *oldBins = new Double_t[fNcells];
//    for (Int_t ibin = 0; ibin < fNcells; ibin++) {
//       oldBins[ibin] = GetBinContent(ibin);
//    }
//    Double_t *oldSumw2 = 0;
//    if (fSumw2.fN != 0) { 
//       oldSumw2 = new Double_t[fNcells];
//       for (Int_t ibin = 0; ibin < fNcells; ibin++) {
//          oldSumw2[ibin] = fSumw2.fArray[ibin];
//       }   
//    }

//    // create a clone of the old histogram if newname is specified
//    TH5 *hnew = this;
//    if (newname && strlen(newname)) {
//       hnew = (TH5*)Clone();
//       hnew->SetName(newname);
//    }

//    // save original statistics
//    Double_t stat[kNstat5];
//    GetStats(stat);
//    bool resetStat = false;


//    // change axis specs and rebuild bin contents array
//    if(newxbins*nxgroup != nxbins) {
//       xmax = fXaxis.GetBinUpEdge(newxbins*nxgroup);
//       resetStat = true; //stats must be reset because top bins will be moved to overflow bin
//    }
//    if(newybins*nygroup != nybins) {
//       ymax = fYaxis.GetBinUpEdge(newybins*nygroup);
//       resetStat = true; //stats must be reset because top bins will be moved to overflow bin
//    }
//    if(newzbins*nzgroup != nzbins) {
//       zmax = fZaxis.GetBinUpEdge(newzbins*nzgroup);
//       resetStat = true; //stats must be reset because top bins will be moved to overflow bin
//    }
//    // save the TAttAxis members (reset by SetBins) for x axis
//    Int_t    nXdivisions  = fXaxis.GetNdivisions();
//    Color_t  xAxisColor   = fXaxis.GetAxisColor();
//    Color_t  xLabelColor  = fXaxis.GetLabelColor();
//    Style_t  xLabelFont   = fXaxis.GetLabelFont();
//    Float_t  xLabelOffset = fXaxis.GetLabelOffset();
//    Float_t  xLabelSize   = fXaxis.GetLabelSize();
//    Float_t  xTickLength  = fXaxis.GetTickLength();
//    Float_t  xTitleOffset = fXaxis.GetTitleOffset();
//    Float_t  xTitleSize   = fXaxis.GetTitleSize();
//    Color_t  xTitleColor  = fXaxis.GetTitleColor();
//    Style_t  xTitleFont   = fXaxis.GetTitleFont();
//    // save the TAttAxis members (reset by SetBins) for y axis
//    Int_t    nYdivisions  = fYaxis.GetNdivisions();
//    Color_t  yAxisColor   = fYaxis.GetAxisColor();
//    Color_t  yLabelColor  = fYaxis.GetLabelColor();
//    Style_t  yLabelFont   = fYaxis.GetLabelFont();
//    Float_t  yLabelOffset = fYaxis.GetLabelOffset();
//    Float_t  yLabelSize   = fYaxis.GetLabelSize();
//    Float_t  yTickLength  = fYaxis.GetTickLength();
//    Float_t  yTitleOffset = fYaxis.GetTitleOffset();
//    Float_t  yTitleSize   = fYaxis.GetTitleSize();
//    Color_t  yTitleColor  = fYaxis.GetTitleColor();
//    Style_t  yTitleFont   = fYaxis.GetTitleFont();
//    // save the TAttAxis members (reset by SetBins) for z axis
//    Int_t    nZdivisions  = fZaxis.GetNdivisions();
//    Color_t  zAxisColor   = fZaxis.GetAxisColor();
//    Color_t  zLabelColor  = fZaxis.GetLabelColor();
//    Style_t  zLabelFont   = fZaxis.GetLabelFont();
//    Float_t  zLabelOffset = fZaxis.GetLabelOffset();
//    Float_t  zLabelSize   = fZaxis.GetLabelSize();
//    Float_t  zTickLength  = fZaxis.GetTickLength();
//    Float_t  zTitleOffset = fZaxis.GetTitleOffset();
//    Float_t  zTitleSize   = fZaxis.GetTitleSize();
//    Color_t  zTitleColor  = fZaxis.GetTitleColor();
//    Style_t  zTitleFont   = fZaxis.GetTitleFont();

//    // copy merged bin contents (ignore under/overflows)
//    if (nxgroup != 1 || nygroup != 1 || nzgroup != 1) {
//       if(fXaxis.GetXbins()->GetSize() > 0 || fYaxis.GetXbins()->GetSize() > 0 || fZaxis.GetXbins()->GetSize() > 0){
//     	 // variable bin sizes in x or y, don't treat both cases separately
//          Double_t *xbins = new Double_t[newxbins+1];
//          for(i = 0; i <= newxbins; ++i) xbins[i] = fXaxis.GetBinLowEdge(1+i*nxgroup);
//          Double_t *ybins = new Double_t[newybins+1];
//          for(i = 0; i <= newybins; ++i) ybins[i] = fYaxis.GetBinLowEdge(1+i*nygroup);
//          Double_t *zbins = new Double_t[newzbins+1];
//          for(i = 0; i <= newzbins; ++i) zbins[i] = fZaxis.GetBinLowEdge(1+i*nzgroup);
//          hnew->SetBins(newxbins,xbins, newybins, ybins, newzbins, zbins);//changes also errors array (if any)
//          delete [] xbins;
//          delete [] ybins;
//          delete [] zbins;
//       } else {
// 	hnew->SetBins(newxbins, xmin, xmax, newybins, ymin, ymax, newzbins, zmin, zmax);//changes also errors array
//       }

//       Double_t binContent, binSumw2;
//       Int_t oldxbin = 1;
//       Int_t oldybin = 1;
//       Int_t oldzbin = 1;
//       Int_t bin;
//       for (xbin = 1; xbin <= newxbins; xbin++) {
// 	oldybin=1;
// 	oldzbin=1;
// 	for (ybin = 1; ybin <= newybins; ybin++) {
// 	  oldzbin=1;
// 	  for (zbin = 1; zbin <= newzbins; zbin++) {
// 	    binContent = 0;
// 	    binSumw2   = 0;
// 	    for (i = 0; i < nxgroup; i++) {
// 	      if (oldxbin+i > nxbins) break;
// 	      for (j =0; j < nygroup; j++) {
// 		if (oldybin+j > nybins) break;
// 		for (k =0; k < nzgroup; k++) {
// 		  if (oldzbin+k > nzbins) break;
// 		  //get global bin (same conventions as in TH1::GetBin(xbin,ybin)
// 		  bin = oldxbin + i + (oldybin + j)*(nxbins + 2) + (oldzbin + k)*(nxbins + 2)*(nybins + 2);
// 		  binContent += oldBins[bin];
// 		  if (oldSumw2) binSumw2 += oldSumw2[bin];
// 		}
// 	      }
// 	    }
//             Int_t ibin = hnew->GetBin(xbin,ybin,zbin);  // new bin number 
// 	    hnew->SetBinContent(ibin, binContent);
// 	    if (oldSumw2) hnew->fSumw2.fArray[ibin] = binSumw2;
// 	    oldzbin += nzgroup;
// 	  }
// 	  oldybin += nygroup;
// 	}
// 	oldxbin += nxgroup;
//       }

//       // compute new underflow/overflows for the 8 vertices 
//       for (Int_t xover = 0; xover <= 1; xover++) { 
//          for (Int_t yover = 0; yover <= 1; yover++) { 
//             for (Int_t zover = 0; zover <= 1; zover++) { 
//                binContent = 0;
//                binSumw2 = 0;
//                // make loop in case of only underflow/overflow
//                for(xbin = xover*oldxbin; xbin <= xover*(nxbins+1); xbin++) {
//                   for(ybin = yover*oldybin; ybin <= yover*(nybins+1); ybin++){
//                      for(zbin = zover*oldzbin; zbin <= zover*(nzbins+1); zbin++){
//                         bin = GetBin(xbin,ybin,zbin);
//                         binContent += oldBins[bin];
//                         if(oldSumw2) binSumw2 += oldSumw2[bin];
//                      }
//                   }
//                }
//                Int_t binNew = hnew->GetBin( xover *(newxbins+1),
//                                             yover*(newybins+1), zover*(newzbins+1) );
//                hnew->SetBinContent(binNew,binContent);                     
//                if (oldSumw2) hnew->fSumw2.fArray[binNew] = binSumw2;
//             }
//          }
//       }         

//       Double_t binContent0, binContent2, binContent3, binContent4;
//       Double_t binError0, binError2, binError3, binError4;
//       Int_t oldxbin2, oldybin2, oldzbin2;
//       Int_t ufbin, ofbin, ofbin2, ofbin3, ofbin4;

//       //  recompute under/overflow contents in y for the new  x and z bins
//       oldxbin2 = 1;
//       oldybin2 = 1;
//       oldzbin2 = 1;
//       for (xbin = 1; xbin<=newxbins; xbin++) {
// 	oldzbin2 = 1;
// 	for (zbin = 1; zbin<=newzbins; zbin++) {
// 	  binContent0 = binContent2 = 0;
// 	  binError0 = binError2 = 0;
// 	  for (i=0; i<nxgroup; i++) {
//             if (oldxbin2+i > nxbins) break;
// 	    for (k=0; k<nzgroup; k++) {
// 	      if (oldzbin2+k > nzbins) break;
	      
// 	      //old underflow bin (in y)
// 	      ufbin = oldxbin2 + i + (nxbins+2)*(nybins+2)*(oldzbin2+k);
// 	      binContent0 += oldBins[ufbin];
// 	      if(oldSumw2) binError0 += oldSumw2[ufbin];
// 	      for(ybin = oldybin; ybin <= nybins + 1; ybin++){
// 		//old overflow bin (in y)
// 		ofbin = ufbin + ybin*(nxbins+2);
// 		binContent2 += oldBins[ofbin];
// 		if(oldSumw2) binError2 += oldSumw2[ofbin];
// 	      }
// 	    }
// 	  }
// 	  hnew->SetBinContent(xbin,0,zbin,binContent0);
// 	  hnew->SetBinContent(xbin,newybins+1,zbin,binContent2);
// 	  if (oldSumw2) {
//             hnew->SetBinError(xbin,0,zbin,TMath::Sqrt(binError0));
//             hnew->SetBinError(xbin,newybins+1,zbin,TMath::Sqrt(binError2) );
// 	  }
// 	  oldzbin2 += nzgroup;
// 	}
// 	oldxbin2 += nxgroup;
//       }

//       //  recompute under/overflow contents in x for the new  y and z bins
//       oldxbin2 = 1;
//       oldybin2 = 1;
//       oldzbin2 = 1;
//       for (ybin = 1; ybin<=newybins; ybin++) {
// 	oldzbin2 = 1;
// 	for (zbin = 1; zbin<=newzbins; zbin++) {
// 	  binContent0 = binContent2 = 0;
// 	  binError0 = binError2 = 0;
// 	  for (j=0; j<nygroup; j++) {
//             if (oldybin2+j > nybins) break;
// 	    for (k=0; k<nzgroup; k++) {
// 	      if (oldzbin2+k > nzbins) break;
	      
// 	      //old underflow bin (in y)
// 	      ufbin = (oldybin2 + j)*(nxbins+2) + (nxbins+2)*(nybins+2)*(oldzbin2+k);
// 	      binContent0 += oldBins[ufbin];
// 	      if(oldSumw2) binError0 += oldSumw2[ufbin];
// 	      for(xbin = oldxbin; xbin <= nxbins + 1; xbin++){
// 		//old overflow bin (in x)
// 		ofbin = ufbin + xbin;
// 		binContent2 += oldBins[ofbin];
// 		if(oldSumw2) binError2 += oldSumw2[ofbin];
// 	      }
// 	    }
// 	  }
// 	  hnew->SetBinContent(0,ybin,zbin,binContent0);
// 	  hnew->SetBinContent(newxbins+1,ybin,zbin,binContent2);
// 	  if (oldSumw2) {
//             hnew->SetBinError(0,ybin,zbin,TMath::Sqrt(binError0));
//             hnew->SetBinError(newxbins+1,ybin,zbin,TMath::Sqrt(binError2) );
// 	  }
// 	  oldzbin2 += nzgroup;
// 	}
// 	oldybin2 += nygroup;
//       }

//       //  recompute under/overflow contents in z for the new  x and y bins
//       oldxbin2 = 1;
//       oldybin2 = 1;
//       oldzbin2 = 1;
//       for (xbin = 1; xbin<=newxbins; xbin++) {
// 	oldybin2 = 1;
// 	for (ybin = 1; ybin<=newybins; ybin++) {
// 	  binContent0 = binContent2 = 0;
// 	  binError0 = binError2 = 0;
// 	  for (i=0; i<nxgroup; i++) {
//             if (oldxbin2+i > nxbins) break;
// 	    for (j=0; j<nygroup; j++) {
// 	      if (oldybin2+j > nybins) break;
	      
// 	      //old underflow bin (in z)
// 	      ufbin = oldxbin2 + i + (nxbins+2)*(oldybin2+j);
// 	      binContent0 += oldBins[ufbin];
// 	      if(oldSumw2) binError0 += oldSumw2[ufbin];
// 	      for(zbin = oldzbin; zbin <= nzbins + 1; zbin++){
// 		//old overflow bin (in z)
// 		ofbin = ufbin + (nxbins+2)*(nybins+2)*zbin;
// 		binContent2 += oldBins[ofbin];
// 		if(oldSumw2) binError2 += oldSumw2[ofbin];
// 	      }
// 	    }
// 	  }
// 	  hnew->SetBinContent(xbin,ybin,0,binContent0);
// 	  hnew->SetBinContent(xbin,ybin,newzbins+1,binContent2);
// 	  if (oldSumw2) {
//             hnew->SetBinError(xbin,ybin,0,TMath::Sqrt(binError0));
//             hnew->SetBinError(xbin,ybin,newzbins+1,TMath::Sqrt(binError2) );
// 	  }
// 	  oldybin2 += nygroup;
// 	}
// 	oldxbin2 += nxgroup;
//       }

//       //  recompute under/overflow contents in y, z for the new  x
//       oldxbin2 = 1;
//       oldybin2 = 1;
//       oldzbin2 = 1;
//       for (xbin = 1; xbin<=newxbins; xbin++) {
// 	  binContent0 = 0;
// 	  binContent2 = 0;
// 	  binContent3 = 0;
// 	  binContent4 = 0;
// 	  binError0 = 0;
// 	  binError2 = 0;
// 	  binError3 = 0;
// 	  binError4 = 0;
// 	  for (i=0; i<nxgroup; i++) {
//               if (oldxbin2+i > nxbins) break;	      
// 	      ufbin = oldxbin2 + i; //
// 	      binContent0 += oldBins[ufbin];
// 	      if(oldSumw2) binError0 += oldSumw2[ufbin];
	      
// 	      for(ybin = oldybin; ybin <= nybins + 1; ybin++){
// 		ofbin3 =  ufbin+ybin*(nxbins+2);
// 		binContent3 += oldBins[ ofbin3 ];
// 		if (oldSumw2)  binError3 += oldSumw2[ofbin3];
// 		for(zbin = oldzbin; zbin <= nzbins + 1; zbin++){
// 		  //old overflow bin (in z)
// 		  ofbin4 =   oldxbin2 + i + ybin*(nxbins+2) + (nxbins+2)*(nybins+2)*zbin;
// 		  binContent4 += oldBins[ofbin4];
// 		  if(oldSumw2) binError4 += oldSumw2[ofbin4];
// 		}
// 	      }
// 	      for(zbin = oldzbin; zbin <= nzbins + 1; zbin++){
// 		ofbin2 =  ufbin+zbin*(nxbins+2)*(nybins+2);
// 		binContent2 += oldBins[ ofbin2 ];
// 		if (oldSumw2)  binError2 += oldSumw2[ofbin2];
// 	      }
// 	  }
// 	  hnew->SetBinContent(xbin,0,0,binContent0);
// 	  hnew->SetBinContent(xbin,0,newzbins+1,binContent2);
// 	  hnew->SetBinContent(xbin,newybins+1,0,binContent3);
// 	  hnew->SetBinContent(xbin,newybins+1,newzbins+1,binContent4);
// 	  if (oldSumw2) {
//             hnew->SetBinError(xbin,0,0,TMath::Sqrt(binError0));	    
//             hnew->SetBinError(xbin,0,newzbins+1,TMath::Sqrt(binError2) );
//             hnew->SetBinError(xbin,newybins+1,0,TMath::Sqrt(binError3) );
//             hnew->SetBinError(xbin,newybins+1,newzbins+1,TMath::Sqrt(binError4) );
// 	  }
// 	  oldxbin2 += nxgroup;
//       }
   

//       //  recompute under/overflow contents in x, y for the new z
//       oldxbin2 = 1;
//       oldybin2 = 1;
//       oldzbin2 = 1;
//       for (zbin = 1; zbin<=newzbins; zbin++) {
// 	  binContent0 = 0;
// 	  binContent2 = 0;
// 	  binContent3 = 0;
// 	  binContent4 = 0;
// 	  binError0 = 0;
// 	  binError2 = 0;
// 	  binError3 = 0;
// 	  binError4 = 0;
// 	  for (i=0; i<nzgroup; i++) {
//               if (oldzbin2+i > nzbins) break;	      
// 	      ufbin = (oldzbin2 + i)*(nxbins+2)*(nybins+2); //
// 	      binContent0 += oldBins[ufbin];
// 	      if(oldSumw2) binError0 += oldSumw2[ufbin];
// 	      for(ybin = oldybin; ybin <= nybins + 1; ybin++){
// 		ofbin3 =  ufbin+ybin*(nxbins+2);
// 		binContent3 += oldBins[ ofbin3 ];
// 		if (oldSumw2)  binError3 += oldSumw2[ofbin3];
// 		for(xbin = oldxbin; xbin <= nxbins + 1; xbin++){
// 		  //old overflow bin (in z)
// 		  ofbin4 = ufbin + xbin + ybin*(nxbins+2);
// 		  binContent4 += oldBins[ofbin4];
// 		  if(oldSumw2) binError4 += oldSumw2[ofbin4];
// 		}
// 	      }
// 	      for(xbin = oldxbin; xbin <= nxbins + 1; xbin++){
// 		ofbin2 =  xbin +(oldzbin2+i)*(nxbins+2)*(nybins+2);
// 		binContent2 += oldBins[ ofbin2 ];
// 		if (oldSumw2)  binError2 += oldSumw2[ofbin2];
// 	      }
// 	  }
// 	  hnew->SetBinContent(0,0,zbin,binContent0);
// 	  hnew->SetBinContent(0,newybins+1,zbin,binContent3);
// 	  hnew->SetBinContent(newxbins+1,0,zbin,binContent2);
// 	  hnew->SetBinContent(newxbins+1,newybins+1,zbin,binContent4);
// 	  if (oldSumw2) {
//             hnew->SetBinError(0,0,zbin,TMath::Sqrt(binError0));	    
//             hnew->SetBinError(0,newybins+1,zbin,TMath::Sqrt(binError3) );
//             hnew->SetBinError(newxbins+1,0,zbin,TMath::Sqrt(binError2) );
//             hnew->SetBinError(newxbins+1,newybins+1,zbin,TMath::Sqrt(binError4) );
// 	  }
// 	  oldzbin2 += nzgroup;
//       }
   

//       //  recompute under/overflow contents in x, z for the new  y
//       oldxbin2 = 1;
//       oldybin2 = 1;
//       oldzbin2 = 1;
//       for (ybin = 1; ybin<=newybins; ybin++) {
// 	  binContent0 = 0;
// 	  binContent2 = 0;
// 	  binContent3 = 0;
// 	  binContent4 = 0;
// 	  binError0 = 0;
// 	  binError2 = 0;
// 	  binError3 = 0;
// 	  binError4 = 0;
// 	  for (i=0; i<nygroup; i++) {
//               if (oldybin2+i > nybins) break;	      
// 	      ufbin = (oldybin2 + i)*(nxbins+2); //
// 	      binContent0 += oldBins[ufbin];
// 	      if(oldSumw2) binError0 += oldSumw2[ufbin];
// 	      for(xbin = oldxbin; xbin <= nxbins + 1; xbin++){
// 		ofbin3 =  ufbin+xbin;
// 		binContent3 += oldBins[ ofbin3 ];
// 		if (oldSumw2)  binError3 += oldSumw2[ofbin3];
// 		for(zbin = oldzbin; zbin <= nzbins + 1; zbin++){
// 		  //old overflow bin (in z)
// 		  ofbin4 = xbin + (nxbins+2)*(nybins+2)*zbin+(oldybin2+i)*(nxbins+2);
// 		  binContent4 += oldBins[ofbin4];
// 		  if(oldSumw2) binError4 += oldSumw2[ofbin4];
// 		}
// 	      }
// 	      for(zbin = oldzbin; zbin <= nzbins + 1; zbin++){
// 		ofbin2 =  (oldybin2+i)*(nxbins+2)+zbin*(nxbins+2)*(nybins+2);
// 		binContent2 += oldBins[ ofbin2 ];
// 		if (oldSumw2)  binError2 += oldSumw2[ofbin2];
// 	      }
// 	  }
// 	  hnew->SetBinContent(0,ybin,0,binContent0);
// 	  hnew->SetBinContent(0,ybin,newzbins+1,binContent2);
// 	  hnew->SetBinContent(newxbins+1,ybin,0,binContent3);
// 	  hnew->SetBinContent(newxbins+1,ybin,newzbins+1,binContent4);
// 	  if (oldSumw2) {
//             hnew->SetBinError(0,ybin,0,TMath::Sqrt(binError0));	    
//             hnew->SetBinError(0,ybin,newzbins+1,TMath::Sqrt(binError2) );
//             hnew->SetBinError(newxbins+1,ybin,0,TMath::Sqrt(binError3) );
//             hnew->SetBinError(newxbins+1,ybin,newzbins+1,TMath::Sqrt(binError4) );
// 	  }
// 	  oldybin2 += nygroup;
//       }
//    }

//    // Restore x axis attributes
//    fXaxis.SetNdivisions(nXdivisions);
//    fXaxis.SetAxisColor(xAxisColor);
//    fXaxis.SetLabelColor(xLabelColor);
//    fXaxis.SetLabelFont(xLabelFont);
//    fXaxis.SetLabelOffset(xLabelOffset);
//    fXaxis.SetLabelSize(xLabelSize);
//    fXaxis.SetTickLength(xTickLength);
//    fXaxis.SetTitleOffset(xTitleOffset);
//    fXaxis.SetTitleSize(xTitleSize);
//    fXaxis.SetTitleColor(xTitleColor);
//    fXaxis.SetTitleFont(xTitleFont);
//    // Restore y axis attributes
//    fYaxis.SetNdivisions(nYdivisions);
//    fYaxis.SetAxisColor(yAxisColor);
//    fYaxis.SetLabelColor(yLabelColor);
//    fYaxis.SetLabelFont(yLabelFont);
//    fYaxis.SetLabelOffset(yLabelOffset);
//    fYaxis.SetLabelSize(yLabelSize);
//    fYaxis.SetTickLength(yTickLength);
//    fYaxis.SetTitleOffset(yTitleOffset);
//    fYaxis.SetTitleSize(yTitleSize);
//    fYaxis.SetTitleColor(yTitleColor);
//    fYaxis.SetTitleFont(yTitleFont);
//    // Restore z axis attributes
//    fZaxis.SetNdivisions(nZdivisions);
//    fZaxis.SetAxisColor(zAxisColor);
//    fZaxis.SetLabelColor(zLabelColor);
//    fZaxis.SetLabelFont(zLabelFont);
//    fZaxis.SetLabelOffset(zLabelOffset);
//    fZaxis.SetLabelSize(zLabelSize);
//    fZaxis.SetTickLength(zTickLength);
//    fZaxis.SetTitleOffset(zTitleOffset);
//    fZaxis.SetTitleSize(zTitleSize);
//    fZaxis.SetTitleColor(zTitleColor);
//    fZaxis.SetTitleFont(zTitleFont);

//    //restore statistics and entries  modified by SetBinContent
//    hnew->SetEntries(entries);
//    if (!resetStat) hnew->PutStats(stat);

//    delete [] oldBins;
//    if (oldSumw2) delete [] oldSumw2;
//    return hnew;
// }


//______________________________________________________________________________
void TH5::Reset(Option_t *option)
{
  //*-*-*-*-*-*-*-*Reset this histogram: contents, errors, etc*-*-*-*-*-*-*-*
  //*-*            ===========================================
  
  TH1::Reset(option);
  TString opt = option;
  opt.ToUpper();
  if (opt.Contains("ICE") && !opt.Contains("S")) return;
  fTsumwy  = 0;
  fTsumwy2 = 0;
  fTsumwxy = 0;
  fTsumwz  = 0;
  fTsumwz2 = 0;
  fTsumwxz = 0;
  fTsumwyz = 0;
  fTsumwyz = 0;
  fTsumwu = 0;
  fTsumwu2 = 0;
  fTsumwxu = 0;
  fTsumwyu = 0;
  fTsumwzu = 0;
  fTsumwv = 0;
  fTsumwv2 = 0;
  fTsumwxv = 0;
  fTsumwyv = 0;
  fTsumwzv = 0;
  fTsumwuv = 0;
}

//______________________________________________________________________________
// void TH5::Streamer(TBuffer &R__b)
// {
//   // Stream an object of class TH5.
  
//   if (R__b.IsReading()) {
//     UInt_t R__s, R__c;
//     Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
//     if (R__v > 2) {
//       R__b.ReadClassBuffer(TH5::Class(), this, R__v, R__s, R__c);
//       return;
//     }
//     //====process old versions before automatic schema evolution
//     TH1::Streamer(R__b);
//     TAtt3D::Streamer(R__b);
//     R__b.CheckByteCount(R__s, R__c, TH5::IsA());
//     //====end of old versions
    
//   } else {
//     R__b.WriteClassBuffer(TH5::Class(),this);
//   }
// }
 
//______________________________________________________________________________
void TH5::UseCurrentStyle()
{
  //   Copy current attributes from/to current style

  TH1::UseCurrentStyle();  
  if (!gStyle) return;
  if (gStyle->IsReading()) {
    fUaxis.ResetAttAxis("X");  // TStyle only knows x, y, z, axes
    fVaxis.ResetAttAxis("X");  // TStyle only knows x, y, z, axes
  }
}

//______________________________________________________________________________
//                     TH5C methods
//  TH5C a 3-D histogram with one byte per cell (char)
//______________________________________________________________________________
 
ClassImp(TH5C)
   
//______________________________________________________________________________
TH5C::TH5C(): TH5(), TArrayC()
{
  // Constructor.
  SetBinsLength(243);
  if (fgDefaultSumw2) Sumw2();
}
  
//______________________________________________________________________________
TH5C::~TH5C()
{
  // Destructor.
}
 
//______________________________________________________________________________
TH5C::TH5C(const char *name, const char *title,
	   Int_t nbinsx, Double_t xlow, Double_t xup,
	   Int_t nbinsy, Double_t ylow, Double_t yup,
	   Int_t nbinsz, Double_t zlow, Double_t zup,
	   Int_t nbinsu, Double_t ulow, Double_t uup,
	   Int_t nbinsv, Double_t vlow, Double_t vup):
  TH5(name, title, nbinsx, xlow, xup, nbinsy, ylow, yup, nbinsz, zlow, zup, nbinsu, ulow, uup, nbinsv, vlow, vup)
{
  //*-*-*-*-*-*-*-*-*Normal constructor for fix bin size 5-D histograms*-*-*-*-*
  //*-*              ==================================================

  TArrayC::Set(fNcells);
  if (fgDefaultSumw2) Sumw2();
  
  if (xlow >= xup || ylow >= yup || zlow >= zup || ulow >= uup || vlow >= vup) SetBuffer(fgBufferSize);
 }

//______________________________________________________________________________
TH5C::TH5C(const char *name, const char *title, 
	   Int_t nbinsx, const Float_t *xbins,
	   Int_t nbinsy, const Float_t *ybins, 
	   Int_t nbinsz, const Float_t *zbins,
	   Int_t nbinsu, const Float_t *ubins,
	   Int_t nbinsv, const Float_t *vbins):
  TH5(name, title, nbinsx, xbins, nbinsy, ybins, nbinsz, zbins, nbinsu, ubins, nbinsv, vbins)
{
  //*-*-*-*-*-*-*-*Normal constructor for variable bin size 5-D histograms*-*-*-*
  //*-*            =======================================================
  TArrayC::Set(fNcells);
  if (fgDefaultSumw2) Sumw2();
}

//______________________________________________________________________________
TH5C::TH5C(const char *name, const char *title,
	   Int_t nbinsx, const Double_t *xbins,
	   Int_t nbinsy, const Double_t *ybins,
	   Int_t nbinsz, const Double_t *zbins,
	   Int_t nbinsu, const Double_t *ubins,
	   Int_t nbinsv, const Double_t *vbins):
  TH5(name, title, nbinsx, xbins, nbinsy, ybins, nbinsz, zbins, nbinsu, ubins, nbinsv, vbins)
{
  //*-*-*-*-*-*-*-*Normal constructor for variable bin size 3-D histograms*-*-*-*
  //*-*            =======================================================
  TArrayC::Set(fNcells);
  if (fgDefaultSumw2) Sumw2();
}

//______________________________________________________________________________
TH5C::TH5C(const TH5C &h3c): TH5(), TArrayC()
{
  // Copy constructor.
  
  ((TH5C&)h3c).Copy(*this);
}

//______________________________________________________________________________
void TH5C::AddBinContent(Int_t bin)
{
  //*-*-*-*-*-*-*-*-*-*Increment bin content by 1*-*-*-*-*-*-*-*-*-*-*-*-*-*
  //*-*                ==========================

   if (fArray[bin] < 127) fArray[bin]++;
}

//______________________________________________________________________________
void TH5C::AddBinContent(Int_t bin, Double_t w)
{
  //*-*-*-*-*-*-*-*-*-*Increment bin content by w*-*-*-*-*-*-*-*-*-*-*-*-*-*
  //*-*                ==========================
  
  Int_t newval = fArray[bin] + Int_t(w);
  if (newval > -128 && newval < 128) {fArray[bin] = Char_t(newval); return;}
  if (newval < -127) fArray[bin] = -127;
  if (newval >  127) fArray[bin] =  127;
}

//______________________________________________________________________________
 void TH5C::Copy(TObject &newth3) const
{
  //*-*-*-*-*-*-*Copy this 5-D histogram structure to newth5*-*-*-*-*-*-*-*-*-*
  //*-*          ===========================================
  
  TH5::Copy((TH5C&)newth3);
}

//______________________________________________________________________________
TH1 *TH5C::DrawCopy(Option_t *option) const
{
  // Draw copy.
  
  TString opt = option;
  opt.ToLower();
  if (gPad && !opt.Contains("same")) gPad->Clear();
  TH5C *newth5 = (TH5C*)Clone();
  newth5->SetDirectory(0);
  newth5->SetBit(kCanDelete);
  newth5->AppendPad(option);
  return newth5;
}

//______________________________________________________________________________
Double_t TH5C::GetBinContent(Int_t bin) const
{
  // Get bin content.
  
  if (fBuffer) ((TH5C*)this)->BufferEmpty();
  if (bin < 0) bin = 0;
  if (bin >= fNcells) bin = fNcells-1;
  if (!fArray) return 0;
  return Double_t (fArray[bin]);
}

//______________________________________________________________________________
Int_t TH5::GetBin(Int_t binx, Int_t biny, Int_t binz, Int_t binu, Int_t binv) const
{
  //   -*-*-*-*Return Global bin number corresponding to binx, y, z, u, v*-*-*-*-*-*-*
  //           ==================================================
  //
  //      5-D histograms are represented with a one dimensional
  //      structure.
  //      This has the advantage that all existing functions, such as
  //        GetBinContent, GetBinError, GetBinFunction work for all dimensions.
  //
  //      Convention for numbering bins
  //      =============================
  //      For all histogram types: nbins, xlow, xup
  //        bin = 0;       underflow bin
  //        bin = 1;       first bin with low-edge xlow INCLUDED
  //        bin = nbins;   last bin with upper-edge xup EXCLUDED
  //        bin = nbins+1; overflow bin
  //   -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  Int_t nx, ny, nz, nu, nv;
  nx = fXaxis.GetNbins()+2;
  if (binx < 0)
    binx = 0;
  else if (binx >= nx) 
    binx = nx - 1;
  ny = fYaxis.GetNbins()+2;
  if (biny < 0)
    biny = 0;
  else  if (biny >= ny)
    biny = ny - 1;
  nz = fZaxis.GetNbins()+2;
  if (binz < 0)
    binz = 0;
  else if (binz >= nz)
    binz = nz - 1;
  nu = fUaxis.GetNbins()+2;
  if (binu < 0)
    binu = 0;
  else if (binu >= nu)
    binu = nu - 1;
  nv = fVaxis.GetNbins()+2;
  if (binv < 0)
    binv = 0;
  else if (binv >= nv)
    binv = nv - 1;

  return binx + nx*(biny +ny*(binz + nz*(binu + nu*binv)));
}

//______________________________________________________________________________
void TH5C::Reset(Option_t *option)
{
  //*-*-*-*-*-*-*-*Reset this histogram: contents, errors, etc*-*-*-*-*-*-*-*
  //*-*            ===========================================
  
  TH5::Reset(option);
  TArrayC::Reset();
  // should also reset statistics once statistics are implemented for TH5
}

//______________________________________________________________________________
void TH5C::SetBinsLength(Int_t n)
{
  // Set total number of bins including under/overflow
  // Reallocate bin contents array
  
  if (n < 0) n = (fXaxis.GetNbins()+2)*(fYaxis.GetNbins()+2)*(fZaxis.GetNbins()+2);
  fNcells = n;
  TArrayC::Set(n);
}

//______________________________________________________________________________
void TH5C::SetBinContent(Int_t bin, Double_t content)
{
  // Set bin content
  fEntries++;
  fTsumw = 0;
  if (bin < 0) return;
  if (bin >= fNcells) return;
  fArray[bin] = Char_t (content);
}


//______________________________________________________________________________
// void TH5::SetShowProjection(const char *option,Int_t nbins)
// {
//    // When the mouse is moved in a pad containing a 3-d view of this histogram
//    // a second canvas shows a projection type given as option.
//    // To stop the generation of the projections, delete the canvas
//    // containing the projection.
//    // option may contain a combination of the characters x,y,z,e
//    // option = "x" return the x projection into a TH1D histogram
//    // option = "y" return the y projection into a TH1D histogram
//    // option = "z" return the z projection into a TH1D histogram
//    // option = "xy" return the x versus y projection into a TH2D histogram
//    // option = "yx" return the y versus x projection into a TH2D histogram
//    // option = "xz" return the x versus z projection into a TH2D histogram
//    // option = "zx" return the z versus x projection into a TH2D histogram
//    // option = "yz" return the y versus z projection into a TH2D histogram
//    // option = "zy" return the z versus y projection into a TH2D histogram
//    // option can also include the drawing option for the projection, eg to draw
//    // the xy projection using the draw option "box" do
//    //   myhist.SetShowProjection("xy box");
//    // This function is typically called from the context menu.
//    // NB: the notation "a vs b" means "a" vertical and "b" horizontal

//    GetPainter();

//    if (fPainter) fPainter->SetShowProjection(option,nbins);
// }

//______________________________________________________________________________
//void TH5C::Streamer(TBuffer &R__b)
// {
//    // Stream an object of class TH5C.

//    if (R__b.IsReading()) {
//       UInt_t R__s, R__c;
//       if (R__b.GetParent() && R__b.GetVersionOwner() < 22300) return;
//       Version_t R__v = R__b.ReadVersion(&R__s, &R__c);
//       if (R__v > 2) {
//          R__b.ReadClassBuffer(TH5C::Class(), this, R__v, R__s, R__c);
//          return;
//       }
//       //====process old versions before automatic schema evolution
//       if (R__v < 2) {
//          R__b.ReadVersion();
//          TH1::Streamer(R__b);
//          TArrayC::Streamer(R__b);
//          R__b.ReadVersion(&R__s, &R__c);
//          TAtt3D::Streamer(R__b);
//       } else {
//          TH5::Streamer(R__b);
//          TArrayC::Streamer(R__b);
//          R__b.CheckByteCount(R__s, R__c, TH5C::IsA());
//       }
//       //====end of old versions

//    } else {
//       R__b.WriteClassBuffer(TH5C::Class(),this);
//    }
// }

//______________________________________________________________________________
TH5C& TH5C::operator=(const TH5C &h1)
{
  // Operator =
  
  if (this != &h1)  ((TH5C&)h1).Copy(*this);
  return *this;
}

//______________________________________________________________________________
TH5C operator*(Float_t c1, TH5C &h1)
{
  // Operator *
  
  TH5C hnew = h1;
  hnew.Scale(c1);
  hnew.SetDirectory(0);
  return hnew;
}

//______________________________________________________________________________
TH5C operator+(TH5C &h1, TH5C &h2)
{
  // Operator +
  
  TH5C hnew = h1;
  hnew.Add(&h2, 1);
  hnew.SetDirectory(0);
  return hnew;
}

//______________________________________________________________________________
TH5C operator-(TH5C &h1, TH5C &h2)
{
  // Operator -
  
  TH5C hnew = h1;
  hnew.Add(&h2, -1);
  hnew.SetDirectory(0);
  return hnew;
}

//______________________________________________________________________________
TH5C operator*(TH5C &h1, TH5C &h2)
{
  // Operator *
  
  TH5C hnew = h1;
  hnew.Multiply(&h2);
  hnew.SetDirectory(0);
  return hnew;
}

//______________________________________________________________________________
TH5C operator/(TH5C &h1, TH5C &h2)
{
  // Operator /
  
  TH5C hnew = h1;
  hnew.Divide(&h2);
  hnew.SetDirectory(0);
  return hnew;
}


