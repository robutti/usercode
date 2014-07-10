//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb 11 17:04:25 2014 by ROOT version 5.34/07
// from TTree trackTree/Fitted track parameters
// found on file: houghCheck_stubs.root
//////////////////////////////////////////////////////////

#ifndef trackStubs_h
#define trackStubs_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

class trackStubs {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   vector<double>  *doca;
   vector<double>  *kappa;
   vector<double>  *phi;
   vector<double>  *z0;
   vector<double>  *theta;
   vector<double>  *algo;
   /* vector<vector<int> > *subdet; */
   /* vector<vector<int> > *layer; */
   /* vector<vector<double> > *xHit; */
   /* vector<vector<double> > *yHit; */
   /* vector<vector<double> > *zHit; */
   UInt_t          nSeeds;
   UInt_t          nGoodSeeds;
   UInt_t          nTriedSeeds;
   UInt_t          nAssSeeds;
   UInt_t          nTracks;
   UInt_t          nTracksAlgo;
   UInt_t          nAccTracks;
   UInt_t          nAssTracks;
   UInt_t          nFoundTracks;

   // List of branches
   TBranch        *b_doca;   //!
   TBranch        *b_kappa;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_z0;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_algo;   //!
   /* TBranch        *b_subdet;   //! */
   /* TBranch        *b_layer;   //! */
   /* TBranch        *b_xHit;   //! */
   /* TBranch        *b_yHit;   //! */
   /* TBranch        *b_zHit;   //! */
   TBranch        *b_nSeeds;   //!
   TBranch        *b_nGoodSeeds;   //!
   TBranch        *b_nTriedSeeds;   //!
   TBranch        *b_nAssSeeds;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_nTracksAlgo;   //!
   TBranch        *b_nAccTracks;   //!
   TBranch        *b_nAssTracks;   //!
   TBranch        *b_nFoundTracks;   //!

   trackStubs(TTree *tree=0);
   virtual ~trackStubs();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef trackStubs_cxx
trackStubs::trackStubs(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("houghCheck_stubs.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("houghCheck_stubs.root");
      }
      f->GetObject("trackTree",tree);

   }
   Init(tree);
}

trackStubs::~trackStubs()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t trackStubs::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t trackStubs::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void trackStubs::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   doca = 0;
   kappa = 0;
   phi = 0;
   z0 = 0;
   theta = 0;
   algo = 0;
   /* subdet = 0; */
   /* layer = 0; */
   /* xHit = 0; */
   /* yHit = 0; */
   /* zHit = 0; */
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("doca", &doca, &b_doca);
   fChain->SetBranchAddress("kappa", &kappa, &b_kappa);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("z0", &z0, &b_z0);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("algo", &algo, &b_algo);
   /* fChain->SetBranchAddress("subdet", &subdet, &b_subdet); */
   /* fChain->SetBranchAddress("layer", &layer, &b_layer); */
   /* fChain->SetBranchAddress("xHit", &xHit, &b_xHit); */
   /* fChain->SetBranchAddress("yHit", &yHit, &b_yHit); */
   /* fChain->SetBranchAddress("zHit", &zHit, &b_zHit); */
   fChain->SetBranchAddress("nSeeds", &nSeeds, &b_nSeeds);
   fChain->SetBranchAddress("nGoodSeeds", &nGoodSeeds, &b_nGoodSeeds);
   fChain->SetBranchAddress("nTriedSeeds", &nTriedSeeds, &b_nTriedSeeds);
   fChain->SetBranchAddress("nAssSeeds", &nAssSeeds, &b_nAssSeeds);
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   fChain->SetBranchAddress("nTracksAlgo", &nTracksAlgo, &b_nTracksAlgo);
   fChain->SetBranchAddress("nAccTracks", &nAccTracks, &b_nAccTracks);
   fChain->SetBranchAddress("nAssTracks", &nAssTracks, &b_nAssTracks);
   fChain->SetBranchAddress("nFoundTracks", &nFoundTracks, &b_nFoundTracks);
   Notify();
}

Bool_t trackStubs::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void trackStubs::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t trackStubs::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef trackStubs_cxx
