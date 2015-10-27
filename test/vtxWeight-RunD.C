#include "TMath.h" 
#include <math.h> 
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"

void vtxWeight(){

  double scale[5] = {
    1,     //data
    2008.4*3, //DY x-sec [pb]
    20508.9*3, //W x-sec [pb]
    831.76,  //TTbar x-sec [pb]
    720648000*0.00042 //QCD 
  };
  TH1D *hNVtx[2];

  TFile *f[4];
  f[0]=TFile::Open("effTrees/DYToLL_M50_Eff_MuMu_74X_25ns_v2/DYToLL_M50Count.root");
  f[1]=TFile::Open("effTrees/WToLNu_Eff_MuMu_74X_25ns_v2/WToLNuCount.root");
  f[2]=TFile::Open("effTrees/TTbar_Eff_MuMu_74X_25ns_v2/TTbarCount.root");
  f[3]=TFile::Open("effTrees/QCD-Pt20-MuEnrPt15_Eff_MuMu_74X_25ns_v2/QCD-Pt20-MuEnrPt15Count.root");
  for(int i=0; i<4; ++i){
    TH1D *h = (TH1D*)f[i]->Get("evtCounter/hCounts");
    float w = h->GetBinContent(2);
    scale[i+1] = scale[i+1]/w;
    if(i==0){
      hNVtx[1] = (TH1D*)f[i]->Get("evtCounter/hNVtx")->Clone("hNVtx_1");
      hNVtx[1]->Reset();
    }
    hNVtx[1]->Add( (TH1D*)f[i]->Get("evtCounter/hNVtx")->Clone(Form("hNVtx_1%i",i)), scale[i+1]/w);
  }

  TFile *fData[2];
  fData[0]=TFile::Open("effTrees/Run2015D_05Oct2015-v1_Cert_246908-258750_Eff_MuMu_74X_25ns_v1/DataCount.root");
  fData[1]=TFile::Open("effTrees/Run2015D_PromptReco-v4_Cert_246908-258750_Eff_MuMu_74X_25ns_v1/DataCount.root");
  for(int i=0; i<2; ++i){
    if(i==0){
      hNVtx[0] = (TH1D*)fData[i]->Get("evtCounter/hNVtx")->Clone("hNVtx_0");
      hNVtx[0]->Reset();
    }
    hNVtx[0]->Add( (TH1D*)fData[i]->Get("evtCounter/hNVtx")->Clone(Form("hNVtx_0%i",i)) );
  }

  hNVtx[0]->Scale(1./hNVtx[0]->Integral(0,hNVtx[0]->GetNbinsX()+1));
  hNVtx[1]->Scale(1./hNVtx[1]->Integral(0,hNVtx[1]->GetNbinsX()+1));
  TH1D *hNVtxNorm=(TH1D*)hNVtx[0]->Clone("hNVtxNorm");
  hNVtxNorm->Divide(hNVtx[1]);

  TFile *fOut=TFile::Open("vtxNorm-RunD.root","RECREATE");
  hNVtxNorm->Write();
  hNVtx[0]->Write();
  hNVtx[1]->Write();
  fOut->Close();

  for(int i=0; i<4; ++i)
    f[i]->Close();
  for(int i=0; i<2; ++i)
    fData[i]->Close();

}
