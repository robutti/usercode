// -*- C++ -*-
//
// Package:    HoughTransChecks
// Class:      TH5
// 
/**\class TH5 TH5.h HoughTest/HoughTransChecks/interface/TH5.h

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

#ifndef ROOT_TH5
#define ROOT_TH5

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// TH5                                                                  //
//                                                                      //
// 5-Dim histogram base class.                                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#ifndef ROOT_TH1
#include "TH1.h"
#endif

//#ifndef ROOT_TAtt3D
//#include "TAtt3D.h"
//#endif

class TH2D; 
class TProfile2D;

class TH5 : public TH1 {

protected:
  Double_t fTsumwy;  // Total Sum of weight*Y
  Double_t fTsumwy2;  // Total Sum of weight*Y*Y
  Double_t fTsumwxy;  // Total Sum of weight*X*Y
  Double_t fTsumwz;  // Total Sum of weight*Z
  Double_t fTsumwz2;  // Total Sum of weight*Z*Z
  Double_t fTsumwxz;  // Total Sum of weight*X*Z
  Double_t fTsumwyz;  // Total Sum of weight*Y*Z
  Double_t fTsumwu;  // Total Sum of weight*U
  Double_t fTsumwu2;  // Total Sum of weight*U*U
  Double_t fTsumwxu;  // Total Sum of weight*X*U
  Double_t fTsumwyu;  // Total Sum of weight*Y*U
  Double_t fTsumwzu;  // Total Sum of weight*Z*U
  Double_t fTsumwv;  // Total Sum of weight*V
  Double_t fTsumwv2;  // Total Sum of weight*V*V
  Double_t fTsumwxv;  // Total Sum of weight*X*V
  Double_t fTsumwyv;  // Total Sum of weight*Y*V
  Double_t fTsumwzv;  // Total Sum of weight*Z*V
  Double_t fTsumwuv;  // Total Sum of weight*U*V
  TAxis fUaxis;  // U axis descriptor
  TAxis fVaxis;  // V axis descriptor

  TH5();
  TH5(const char *name, const char *title,
      Int_t nbinsx, Double_t xlow, Double_t xup,
      Int_t nbinsy, Double_t ylow, Double_t yup,
      Int_t nbinsz, Double_t zlow, Double_t zup,
      Int_t nbinsu, Double_t ulow, Double_t uup,
      Int_t nbinsv, Double_t vlow, Double_t vup);
  TH5(const char *name, const char *title,
      Int_t nbinsx, const Float_t *xbins,
      Int_t nbinsy, const Float_t *ybins,
      Int_t nbinsz, const Float_t *zbins,
      Int_t nbinsu, const Float_t *ubins,
      Int_t nbinsv, const Float_t *vbins);
  TH5(const char *name, const char *title,
      Int_t nbinsx, const Double_t *xbins,
      Int_t nbinsy, const Double_t *ybins,
      Int_t nbinsz, const Double_t *zbins,
      Int_t nbinsu, const Double_t *ubins,
      Int_t nbinsv, const Double_t *vbins);
  virtual Int_t BufferFill(Double_t x, Double_t y, Double_t z, Double_t u, Double_t v, Double_t w);
  
  Int_t Fill(Double_t);  // MayNotUse
  Int_t Fill(Double_t,Double_t) {return Fill(0.);}  // MayNotUse
  Int_t Fill(const char*, Double_t) {return Fill(0);}  // MayNotUse
  
public:
  enum {kNstat5 = 24};
  TH5(const TH5&);
  virtual ~TH5();
  virtual Int_t BufferEmpty(Int_t action=0);
  virtual void Copy(TObject &hnew) const;
  virtual Int_t Fill(Double_t x, Double_t y, Double_t z, Double_t u, Double_t v);
  virtual Int_t Fill(Double_t x, Double_t y, Double_t z, Double_t u, Double_t v, Double_t w);
  
//  virtual void FillRandom(const char *fname, Int_t ntimes=5000);
//  virtual void FillRandom(TH1 *h, Int_t ntimes=5000);
  
  Int_t FindBin(Double_t x, Double_t y, Double_t z, Double_t u, Double_t v);
  Int_t GetBin(Int_t binx, Int_t biny, Int_t binz, Int_t binu, Int_t binv) const;
  using TH1::GetBinErrorLow;
  using TH1::GetBinErrorUp;
//   virtual Double_t GetBinErrorLow(Int_t binx, Int_t biny, Int_t binz) { return TH1::GetBinErrorLow( GetBin(binx, biny, binz) ); }
//   virtual Double_t GetBinErrorUp(Int_t binx, Int_t biny, Int_t binz)  { return TH1::GetBinErrorUp( GetBin(binx, biny, binz) ); }
  virtual Double_t GetCorrelationFactor(Int_t axis1=1,Int_t axis2=2) const;
  virtual Double_t GetCovariance(Int_t axis1=1,Int_t axis2=2) const;
  virtual Double_t GetMaximum(Double_t maxval=FLT_MAX) const;
  Int_t GetNbinsU() const {return fUaxis.GetNbins();}
  Int_t GetNbinsV() const {return fVaxis.GetNbins();}
//  virtual void     GetRandom5(Double_t &x, Double_t &y, Double_t &z, Double_t &u, Double_t &v);
  virtual void     GetStats(Double_t *stats) const;
  TAxis* GetUaxis() const;
  TAxis* GetVaxis() const;
//   virtual Double_t Integral(Option_t *option="") const;
  using TH1::Integral;
//  virtual Double_t Integral(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Int_t binu1, Int_t binu2, Int_t binv1, Int_t binv2, Option_t *option="") const;
  using TH1::IntegralAndError;
//  virtual Double_t IntegralAndError(Int_t binx1, Int_t binx2, Int_t biny1, Int_t biny2, Int_t binz1, Int_t binz2, Int_t binu1, Int_t binu2, Int_t binv1, Int_t binv2, Double_t & err, Option_t *option="") const;
//  virtual Double_t Interpolate(Double_t x);
//  virtual Double_t Interpolate(Double_t x, Double_t y);
//  virtual Double_t Interpolate(Double_t x, Double_t y, Double_t z);
//  virtual Double_t Interpolate(Double_t x, Double_t y, Double_t z, Double_t u, Double_t v);
//  virtual Double_t KolmogorovTest(const TH1 *h2, Option_t *option="") const;
//  virtual Long64_t Merge(TCollection *list);
  TH1D* ProjectionX(const char *name="_px", Int_t firstybin=0, Int_t lastybin=-1, Int_t firstzbin=0, Int_t lastzbin=-1,
		    Int_t firstubin=0, Int_t lastubin=-1, Int_t firstvbin=0, Int_t lastvbin=-1, Option_t *option="") const; // *MENU*
  TH1D* ProjectionY(const char *name="_py", Int_t firstxbin=0, Int_t lastxbin=-1, Int_t firstzbin=0, Int_t lastzbin=-1,
		    Int_t firstubin=0, Int_t lastubin=-1, Int_t firstvbin=0, Int_t lastvbin=-1, Option_t *option="") const; // *MENU*
  TH1D* ProjectionZ(const char *name="_pz", Int_t firstxbin=0, Int_t lastxbin=-1, Int_t firstybin=0, Int_t lastybin=-1,
		    Int_t firstubin=0, Int_t lastubin=-1, Int_t firstvbin=0, Int_t lastvbin=-1, Option_t *option="") const; // *MENU*
  TH1D* ProjectionU(const char *name="_pu", Int_t firstxbin=0, Int_t lastxbin=-1, Int_t firstybin=0, Int_t lastybin=-1,
		    Int_t firstzbin=0, Int_t lastzbin=-1, Int_t firstvbin=0, Int_t lastvbin=-1, Option_t *option="") const; // *MENU*
  TH1D* ProjectionV(const char *name="_pv", Int_t firstxbin=0, Int_t lastxbin=-1, Int_t firstybin=0, Int_t lastybin=-1,
		    Int_t firstzbin=0, Int_t lastzbin=-1, Int_t firstubin=0, Int_t lastubin=-1, Option_t *option="") const; // *MENU*
  TH1* Project5D(Option_t *option="x") const; // *MENU*
  virtual void PutStats(Double_t *stats);
//   virtual TH5* RebinX(Int_t ngroup = 2, const char *newname = "");
//   virtual TH5* RebinY(Int_t ngroup = 2, const char *newname = "");
//   virtual TH5* RebinZ(Int_t ngroup = 2, const char *newname = "");
//   virtual TH5* RebinU(Int_t ngroup = 2, const char *newname = "");
//   virtual TH5* RebinV(Int_t ngroup = 2, const char *newname = "");
//  virtual TH5* Rebin3D(Int_t nxgroup = 2, Int_t nygroup = 2, Int_t nzgroup = 2, Int_t nugroup = 2, Int_t nvgroup = 2, const char *newname = "");
  virtual void Reset(Option_t *option="");
  //  virtual void SetShowProjection(const char *option="xy",Int_t nbins=1);   // *MENU*

  virtual void UseCurrentStyle();

protected:
  TH1D* DoProject1D(const char* name, const char * title, TAxis* projX, 
                    bool computeErrors, bool originalRange,
		    bool useUF, bool useOF) const;
  TH2D* DoProject2D(const char* name, const char * title, TAxis* projX, TAxis* projY, 
		    bool computeErrors, bool originalRange,
		    bool useUF, bool useOF) const;
   
  ClassDef(TH5, 1)  // 5-Dim histogram base class
};

//________________________________________________________________________

class TH5C : public TH5, public TArrayC {
public:
  TH5C();
  TH5C(const char *name, const char *title,
       Int_t nbinsx, Double_t xlow, Double_t xup,
       Int_t nbinsy, Double_t ylow, Double_t yup,
       Int_t nbinsz, Double_t zlow, Double_t zup,
       Int_t nbinsu, Double_t ulow, Double_t uup,
       Int_t nbinsv, Double_t vlow, Double_t vup);
  TH5C(const char *name, const char *title, 
       Int_t nbinsx, const Float_t *xbins,
       Int_t nbinsy, const Float_t *ybins, 
       Int_t nbinsz, const Float_t *zbins,
       Int_t nbinsu, const Float_t *ubins,
       Int_t nbinsv, const Float_t *vbins);
  TH5C(const char *name, const char *title,
       Int_t nbinsx, const Double_t *xbins,
       Int_t nbinsy, const Double_t *ybins,
       Int_t nbinsz, const Double_t *zbins,
       Int_t nbinsu, const Double_t *ubins,
       Int_t nbinsv, const Double_t *vbins);
  TH5C(const TH5C &h5c);
  virtual ~TH5C();
  virtual void AddBinContent(Int_t bin);
  virtual void AddBinContent(Int_t bin, Double_t w);
  virtual void Copy(TObject &hnew) const;
  virtual TH1* DrawCopy(Option_t *option="") const;
  virtual Double_t GetBinContent(Int_t bin) const;
  virtual Double_t GetBinContent(Int_t binx, Int_t biny) const {return GetBinContent(binx);}
  virtual Double_t GetBinContent(Int_t binx, Int_t biny, Int_t binz) const {return GetBinContent(binx);}
  Double_t GetBinContent(Int_t binx, Int_t biny, Int_t binz, Int_t binu, Int_t binv) const {return GetBinContent(GetBin(binx, biny, binz, binu, binv));}
  virtual void Reset(Option_t *option="");
  virtual void SetBinContent(Int_t bin, Double_t content);
  virtual void SetBinContent(Int_t binx, Int_t biny, Double_t content) {SetBinContent(binx, content);}
  virtual void SetBinContent(Int_t binx, Int_t biny, Int_t binz, Double_t content) {SetBinContent(binx, content);}
  void SetBinContent(Int_t binx, Int_t biny, Int_t binz, Int_t binu, Int_t binv, Double_t content) {SetBinContent(GetBin(binx, biny, binz, binu, binv), content);}
  virtual void SetBinsLength(Int_t n = -1);
  TH5C& operator=(const TH5C &h1);
  friend TH5C operator*(Float_t c1, TH5C &h1);
  friend TH5C operator*(TH5C &h1, Float_t c1) {return operator*(c1,h1);}
  friend TH5C operator+(TH5C &h1, TH5C &h2);
  friend TH5C operator-(TH5C &h1, TH5C &h2);
  friend TH5C operator*(TH5C &h1, TH5C &h2);
  friend TH5C operator/(TH5C &h1, TH5C &h2);

  ClassDef(TH5C, 1)  // 5-Dim histograms (one char per channel)

};

//________________________________________________________________________

//class TH5S : public TH5, public TArrayS {
//public:
//   TH5S();
//   TH5S(const char *name,const char *title,Int_t nbinsx,Double_t xlow,Double_t xup
//                                  ,Int_t nbinsy,Double_t ylow,Double_t yup
//                                  ,Int_t nbinsz,Double_t zlow,Double_t zup);
//   TH5S(const char *name,const char *title,Int_t nbinsx,const Float_t *xbins
//                                          ,Int_t nbinsy,const Float_t *ybins
//                                          ,Int_t nbinsz,const Float_t *zbins);
//   TH5S(const char *name,const char *title,Int_t nbinsx,const Double_t *xbins
//                                          ,Int_t nbinsy,const Double_t *ybins
//                                          ,Int_t nbinsz,const Double_t *zbins);
//   TH5S(const TH5S &h3s);
//   virtual ~TH5S();
//   virtual void      AddBinContent(Int_t bin);
//   virtual void      AddBinContent(Int_t bin, Double_t w);
//   virtual void      Copy(TObject &hnew) const;
//   virtual TH1      *DrawCopy(Option_t *option="") const;
//   virtual Double_t  GetBinContent(Int_t bin) const;
//   virtual Double_t  GetBinContent(Int_t bin, Int_t) const {return GetBinContent(bin);}
//   virtual Double_t  GetBinContent(Int_t binx, Int_t biny, Int_t binz) const {return GetBinContent(GetBin(binx,biny,binz));}
//   virtual void      Reset(Option_t *option="");
//   virtual void      SetBinContent(Int_t bin, Double_t content);
//   virtual void      SetBinContent(Int_t bin, Int_t, Double_t content) {SetBinContent(bin,content);}
//   virtual void      SetBinContent(Int_t binx, Int_t biny, Int_t binz, Double_t content) {SetBinContent(GetBin(binx,biny,binz),content);}
//   virtual void      SetBinsLength(Int_t n=-1);
//           TH5S&     operator=(const TH5S &h1);
//   friend  TH5S      operator*(Float_t c1, TH5S &h1);
//   friend  TH5S      operator*(TH5S &h1, Float_t c1) {return operator*(c1,h1);}
//   friend  TH5S      operator+(TH5S &h1, TH5S &h2);
//   friend  TH5S      operator-(TH5S &h1, TH5S &h2);
//   friend  TH5S      operator*(TH5S &h1, TH5S &h2);
//   friend  TH5S      operator/(TH5S &h1, TH5S &h2);
//
//   ClassDef(TH5S,3)  //3-Dim histograms (one short per channel)
//};

#endif
