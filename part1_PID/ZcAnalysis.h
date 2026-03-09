//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jan 11 13:04:48 2023 by ROOT version 6.24/02
// from TTree Zc3900/Zc3900 Tuple
// found on file: Zc3900Alg_11259576_0.root
//////////////////////////////////////////////////////////

#ifndef ZcAnalysis_h
#define ZcAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

// Header file for the classes stored in the TTree if any.
class ZcAnalysis {
public:
    TTree* fChain;  //!pointer to the analyzed TTree or TChain
    Int_t fCurrent; //!current Tree number in a TChain

    // Fixed size dimensions of array or collections stored in the TTree if any.

    // Declaration of leaf types
    Int_t           intRunNumber;
    Int_t           intEventNumber;
    Double_t        dblCmsEnergy;
    Int_t           intNumberChargedTracks;
    Int_t           intNumberShowers;
    Int_t           intNumberTotalTracks;
    Int_t           intNumberGoodChargedTracks;
    Int_t           intTracksIndex;
    Double_t        dblTracksECal[4];   //[intTracksIndex]
    Double_t        dblTracksPx[4];   //[intTracksIndex]
    Double_t        dblTracksPy[4];   //[intTracksIndex]
    Double_t        dblTracksPz[4];   //[intTracksIndex]
    Double_t        dblTracksEoP[4];   //[intTracksIndex]
    Double_t        dblTracksCharge[4];   //[intTracksIndex]
    Double_t        dblTracksMucDepth[4];   //[intTracksIndex]
    Double_t        dblTracksProbElectron[4];   //[intTracksIndex]
    Double_t        dblTracksProbMuon[4];   //[intTracksIndex]
    Double_t        dblTracksProbPion[4];   //[intTracksIndex]
    Double_t        dblTracksProbKaon[4];   //[intTracksIndex]
    Double_t        dblTracksProbProton[4];   //[intTracksIndex]
    Int_t           intNumberGoodPhotons;
    Int_t           intPhotonsIndex;
    Double_t        dblPhotonsE[20];   //[intPhotonsIndex]
    Double_t        dblPhotonsPx[20];   //[intPhotonsIndex]
    Double_t        dblPhotonsPy[20];   //[intPhotonsIndex]
    Double_t        dblPhotonsPz[20];   //[intPhotonsIndex]
    Double_t        dblPhotonsEmcTime[20];   //[intPhotonsIndex]
    Double_t        dblPhotonsSecondMoment[20];   //[intPhotonsIndex]
    Double_t        dblPhotonsLateralMoment[20];   //[intPhotonsIndex]
    Double_t        dblPhotonsA20Moment[20];   //[intPhotonsIndex]
    Double_t        dblPhotonsA42Moment[20];   //[intPhotonsIndex]
    Double_t        dblPhotonsE9oE25[20];   //[intPhotonsIndex]
    Double_t        dblVertexPositionX;
    Double_t        dblVertexPositionY;
    Double_t        dblVertexPositionZ;
    Int_t           intNumberPions;
    Int_t           intNumberElectrons;
    Int_t           intNumberMuons;
    Double_t        dblKinFit4CChiSq;
    Int_t           intKinFit4CMomentaIndex;
    Double_t        dblKinFit4CMomentaE[4];   //[intKinFit4CMomentaIndex]
    Double_t        dblKinFit4CMomentaPx[4];   //[intKinFit4CMomentaIndex]
    Double_t        dblKinFit4CMomentaPy[4];   //[intKinFit4CMomentaIndex]
    Double_t        dblKinFit4CMomentaPz[4];   //[intKinFit4CMomentaIndex]
    Int_t           intNumberTrueParticles;
    Int_t           intTrueIndex;
    Int_t           intTruePdgId[50];   //[intTrueIndex]
    Double_t        dblTrueE[50];   //[intTrueIndex]
    Double_t        dblTruePx[50];   //[intTrueIndex]
    Double_t        dblTruePy[50];   //[intTrueIndex]
    Double_t        dblTruePz[50];   //[intTrueIndex]
    Int_t           intTruePrimary[50];   //[intTrueIndex]
    Int_t           intTrueLeaf[50];   //[intTrueIndex]
    Int_t           intTrueDecayFromGenerator[50];   //[intTrueIndex]
    Int_t           intTrueDecayInFlight[50];   //[intTrueIndex]
    Int_t           intTrueTrackIndex[50];   //[intTrueIndex]
    Int_t           intTrueVertexIndex0[50];   //[intTrueIndex]
    Int_t           intTrueVertexIndex1[50];   //[intTrueIndex]
    Int_t           intTrueMotherTrackIndex[50];   //[intTrueIndex]
    Int_t           intTrueMotherPdgId[50];   //[intTrueIndex]

    // List of branches
    TBranch* b_intRunNumber;
    TBranch* b_intEventNumber;
    TBranch* b_dblCmsEnergy;
    TBranch* b_intNumberChargedTracks;
    TBranch* b_intNumberShowers;
    TBranch* b_intNumberTotalTracks;
    TBranch* b_intNumberGoodChargedTracks;
    TBranch* b_intTracksIndex;
    TBranch* b_dblTracksECal;
    TBranch* b_dblTracksPx;
    TBranch* b_dblTracksPy;
    TBranch* b_dblTracksPz;
    TBranch* b_dblTracksEoP;
    TBranch* b_dblTracksCharge;
    TBranch* b_dblTracksMucDepth;
    TBranch* b_dblTracksProbElectron;
    TBranch* b_dblTracksProbMuon;
    TBranch* b_dblTracksProbPion;
    TBranch* b_dblTracksProbKaon;
    TBranch* b_dblTracksProbProton;
    TBranch* b_intNumberGoodPhotons;
    TBranch* b_intPhotonsIndex;
    TBranch* b_dblPhotonsE;
    TBranch* b_dblPhotonsPx;
    TBranch* b_dblPhotonsPy;
    TBranch* b_dblPhotonsPz;
    TBranch* b_dblPhotonsEmcTime;
    TBranch* b_dblPhotonsSecondMoment;
    TBranch* b_dblPhotonsLateralMoment;
    TBranch* b_dblPhotonsA20Moment;
    TBranch* b_dblPhotonsA42Moment;
    TBranch* b_dblPhotonsE9oE25;
    TBranch* b_dblVertexPositionX;
    TBranch* b_dblVertexPositionY;
    TBranch* b_dblVertexPositionZ;
    TBranch* b_intNumberPions;
    TBranch* b_intNumberElectrons;
    TBranch* b_intNumberMuons;
    TBranch* b_dblKinFit4CChiSq;
    TBranch* b_intKinFit4CMomentaIndex;
    TBranch* b_dblKinFit4CMomentaE;
    TBranch* b_dblKinFit4CMomentaPx;
    TBranch* b_dblKinFit4CMomentaPy;
    TBranch* b_dblKinFit4CMomentaPz;
    TBranch* b_intNumberTrueParticles;
    TBranch* b_intTrueIndex;
    TBranch* b_intTruePdgId;
    TBranch* b_dblTrueE;
    TBranch* b_dblTruePx;
    TBranch* b_dblTruePy;
    TBranch* b_dblTruePz;
    TBranch* b_intTruePrimary;
    TBranch* b_intTrueLeaf;
    TBranch* b_intTrueDecayFromGenerator;
    TBranch* b_intTrueDecayInFlight;
    TBranch* b_intTrueTrackIndex;
    TBranch* b_intTrueVertexIndex0;
    TBranch* b_intTrueVertexIndex1;
    TBranch* b_intTrueMotherTrackIndex;
    TBranch* b_intTrueMotherPdgId;

    // Functions
    ZcAnalysis(TString dataPath, TString savePath);
    virtual ~ZcAnalysis();
    virtual Int_t    Cut(Long64_t entry);
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init(TTree* tree);
    virtual void     Loop(TString savePath);
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef ZcAnalysis_cxx
ZcAnalysis::ZcAnalysis(TString dataPath, TString savePath) : fChain(0)
{
    TFile* f = new TFile(dataPath); // Load file from dataPath
    TTree* tree = (TTree*) f->Get("Zc3900"); // Get tree called "Zc3900" from file
    Init(tree); // Set branch addresses and branch pointers of tree
    Loop(savePath); // Execute Loop function with current savePath
}

ZcAnalysis::~ZcAnalysis()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t ZcAnalysis::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t ZcAnalysis::LoadTree(Long64_t entry)
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

void ZcAnalysis::Init(TTree* tree)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    // Set branch addresses and branch pointers
    if (!tree) return;
    fChain = tree;
    fCurrent = -1;
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("intRunNumber", &intRunNumber, &b_intRunNumber);
    fChain->SetBranchAddress("intEventNumber", &intEventNumber, &b_intEventNumber);
    fChain->SetBranchAddress("dblCmsEnergy", &dblCmsEnergy, &b_dblCmsEnergy);
    fChain->SetBranchAddress("intNumberChargedTracks", &intNumberChargedTracks, &b_intNumberChargedTracks);
    fChain->SetBranchAddress("intNumberShowers", &intNumberShowers, &b_intNumberShowers);
    fChain->SetBranchAddress("intNumberTotalTracks", &intNumberTotalTracks, &b_intNumberTotalTracks);
    fChain->SetBranchAddress("intNumberGoodChargedTracks", &intNumberGoodChargedTracks, &b_intNumberGoodChargedTracks);
    fChain->SetBranchAddress("intTracksIndex", &intTracksIndex, &b_intTracksIndex);
    fChain->SetBranchAddress("dblTracksECal", dblTracksECal, &b_dblTracksECal);
    fChain->SetBranchAddress("dblTracksPx", dblTracksPx, &b_dblTracksPx);
    fChain->SetBranchAddress("dblTracksPy", dblTracksPy, &b_dblTracksPy);
    fChain->SetBranchAddress("dblTracksPz", dblTracksPz, &b_dblTracksPz);
    fChain->SetBranchAddress("dblTracksEoP", dblTracksEoP, &b_dblTracksEoP);
    fChain->SetBranchAddress("dblTracksCharge", dblTracksCharge, &b_dblTracksCharge);
    fChain->SetBranchAddress("dblTracksMucDepth", dblTracksMucDepth, &b_dblTracksMucDepth);
    fChain->SetBranchAddress("dblTracksProbElectron", dblTracksProbElectron, &b_dblTracksProbElectron);
    fChain->SetBranchAddress("dblTracksProbMuon", dblTracksProbMuon, &b_dblTracksProbMuon);
    fChain->SetBranchAddress("dblTracksProbPion", dblTracksProbPion, &b_dblTracksProbPion);
    fChain->SetBranchAddress("dblTracksProbKaon", dblTracksProbKaon, &b_dblTracksProbKaon);
    fChain->SetBranchAddress("dblTracksProbProton", dblTracksProbProton, &b_dblTracksProbProton);
    fChain->SetBranchAddress("intNumberGoodPhotons", &intNumberGoodPhotons, &b_intNumberGoodPhotons);
    fChain->SetBranchAddress("intPhotonsIndex", &intPhotonsIndex, &b_intPhotonsIndex);
    fChain->SetBranchAddress("dblPhotonsE", dblPhotonsE, &b_dblPhotonsE);
    fChain->SetBranchAddress("dblPhotonsPx", dblPhotonsPx, &b_dblPhotonsPx);
    fChain->SetBranchAddress("dblPhotonsPy", dblPhotonsPy, &b_dblPhotonsPy);
    fChain->SetBranchAddress("dblPhotonsPz", dblPhotonsPz, &b_dblPhotonsPz);
    fChain->SetBranchAddress("dblPhotonsEmcTime", dblPhotonsEmcTime, &b_dblPhotonsEmcTime);
    fChain->SetBranchAddress("dblPhotonsSecondMoment", dblPhotonsSecondMoment, &b_dblPhotonsSecondMoment);
    fChain->SetBranchAddress("dblPhotonsLateralMoment", dblPhotonsLateralMoment, &b_dblPhotonsLateralMoment);
    fChain->SetBranchAddress("dblPhotonsA20Moment", dblPhotonsA20Moment, &b_dblPhotonsA20Moment);
    fChain->SetBranchAddress("dblPhotonsA42Moment", dblPhotonsA42Moment, &b_dblPhotonsA42Moment);
    fChain->SetBranchAddress("dblPhotonsE9oE25", dblPhotonsE9oE25, &b_dblPhotonsE9oE25);
    fChain->SetBranchAddress("dblVertexPositionX", &dblVertexPositionX, &b_dblVertexPositionX);
    fChain->SetBranchAddress("dblVertexPositionY", &dblVertexPositionY, &b_dblVertexPositionY);
    fChain->SetBranchAddress("dblVertexPositionZ", &dblVertexPositionZ, &b_dblVertexPositionZ);
    fChain->SetBranchAddress("intNumberPions", &intNumberPions, &b_intNumberPions);
    fChain->SetBranchAddress("intNumberElectrons", &intNumberElectrons, &b_intNumberElectrons);
    fChain->SetBranchAddress("intNumberMuons", &intNumberMuons, &b_intNumberMuons);
    fChain->SetBranchAddress("dblKinFit4CChiSq", &dblKinFit4CChiSq, &b_dblKinFit4CChiSq);
    fChain->SetBranchAddress("intKinFit4CMomentaIndex", &intKinFit4CMomentaIndex, &b_intKinFit4CMomentaIndex);
    fChain->SetBranchAddress("dblKinFit4CMomentaE", dblKinFit4CMomentaE, &b_dblKinFit4CMomentaE);
    fChain->SetBranchAddress("dblKinFit4CMomentaPx", dblKinFit4CMomentaPx, &b_dblKinFit4CMomentaPx);
    fChain->SetBranchAddress("dblKinFit4CMomentaPy", dblKinFit4CMomentaPy, &b_dblKinFit4CMomentaPy);
    fChain->SetBranchAddress("dblKinFit4CMomentaPz", dblKinFit4CMomentaPz, &b_dblKinFit4CMomentaPz);
    fChain->SetBranchAddress("intNumberTrueParticles", &intNumberTrueParticles, &b_intNumberTrueParticles);
    fChain->SetBranchAddress("intTrueIndex", &intTrueIndex, &b_intTrueIndex);
    fChain->SetBranchAddress("intTruePdgId", intTruePdgId, &b_intTruePdgId);
    fChain->SetBranchAddress("dblTrueE", dblTrueE, &b_dblTrueE);
    fChain->SetBranchAddress("dblTruePx", dblTruePx, &b_dblTruePx);
    fChain->SetBranchAddress("dblTruePy", dblTruePy, &b_dblTruePy);
    fChain->SetBranchAddress("dblTruePz", dblTruePz, &b_dblTruePz);
    fChain->SetBranchAddress("intTruePrimary", intTruePrimary, &b_intTruePrimary);
    fChain->SetBranchAddress("intTrueLeaf", intTrueLeaf, &b_intTrueLeaf);
    fChain->SetBranchAddress("intTrueDecayFromGenerator", intTrueDecayFromGenerator, &b_intTrueDecayFromGenerator);
    fChain->SetBranchAddress("intTrueDecayInFlight", intTrueDecayInFlight, &b_intTrueDecayInFlight);
    fChain->SetBranchAddress("intTrueTrackIndex", intTrueTrackIndex, &b_intTrueTrackIndex);
    fChain->SetBranchAddress("intTrueVertexIndex0", intTrueVertexIndex0, &b_intTrueVertexIndex0);
    fChain->SetBranchAddress("intTrueVertexIndex1", intTrueVertexIndex1, &b_intTrueVertexIndex1);
    fChain->SetBranchAddress("intTrueMotherTrackIndex", intTrueMotherTrackIndex, &b_intTrueMotherTrackIndex);
    fChain->SetBranchAddress("intTrueMotherPdgId", intTrueMotherPdgId, &b_intTrueMotherPdgId);
    Notify();
}

Bool_t ZcAnalysis::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void ZcAnalysis::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}

Int_t ZcAnalysis::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif // #ifdef ZcAnalysis_cxx
