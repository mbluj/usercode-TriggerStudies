#include "TMath.h" 
#include <math.h> 
#include "TFile.h"
#include "TChain.h"
#include "TH1D.h"


void vtxWeight(){

  TChain *t[5];
  //data
  t[0] = new TChain("muMu/muMuTriggerTree");
  t[0]->Add("effTrees/Run2015D_05Oct2015-v1_Cert_246908-258750_Eff_MuMu_74X_25ns_v1/muMuTrgAna*.root"); //552.67pb-1 (after brilcalc tool)
  t[0]->Add("effTrees/Run2015D_PromptReco-v4_Cert_246908-258750_Eff_MuMu_74X_25ns_v1/muMuTrgAna*.root"); //711.21pb-1 (after brilcalc tool)
  //DY
  t[1] = new TChain("muMu/muMuTriggerTree");
  t[1]->Add("effTrees/DYToLL_M50_Eff_MuMu_74X_25ns_v2/muMuTrgAna*.root");
  //W
  t[2] = new TChain("muMu/muMuTriggerTree");
  t[2]->Add("effTrees/WToLNu_Eff_MuMu_74X_25ns_v2/muMuTrgAna*.root");
  //TTbar
  t[3] = new TChain("muMu/muMuTriggerTree");
  t[3]->Add("effTrees/TTbar_Eff_MuMu_74X_25ns_v2/muMuTrgAna*.root");
  //QCD
  t[4] = new TChain("muMu/muMuTriggerTree");
  t[4]->Add("effTrees/QCD-Pt20-MuEnrPt15_Eff_MuMu_74X_25ns_v2/muMuTrgAna*.root");

  double scale[5] = {
    1,     //data
    2008.4*3, //DY x-sec [pb]
    20508.9*3, //W x-sec [pb]
    831.76,  //TTbar x-sec [pb]
    720648000*0.00042 //QCD 
  };

  TFile *f[4];
  f[0]=TFile::Open("effTrees/DYToLL_M50_Eff_MuMu_74X_25ns_v2/DYToLL_M50Count.root");
  f[1]=TFile::Open("effTrees/WToLNu_Eff_MuMu_74X_25ns_v2/WToLNuCount.root");
  f[2]=TFile::Open("effTrees/TTbar_Eff_MuMu_74X_25ns_v2/TTbarCount.root");
  f[3]=TFile::Open("effTrees/QCD-Pt20-MuEnrPt15_Eff_MuMu_74X_25ns_v2/QCD-Pt20-MuEnrPt15Count.root");
  for(int i=0; i<4; ++i){
    TH1D *h = (TH1D*)f[i]->Get("evtCounter/hCounts");
    float w = h->GetBinContent(2);
    scale[i+1] = scale[i+1]/w;
    f[i]->Close();
  }

  TH1D *hNVtx[2];
  for(int i=0; i<2; ++i){
    //hNVtx[i] = new TH1D(Form("hNVtx_%i",i)," ;no. of vertices; enties",25,0,50);
    hNVtx[i] = new TH1D(Form("hNVtx_%i",i)," ;no. of vertices; enties",50,0,50);
    hNVtx[i]->Sumw2();
    hNVtx[i]->SetMarkerStyle(20);
    //hNVtx[i]->SetMarkerSize(0.7);
    hNVtx[i]->SetStats(0);
    //hNVtx[i]->SetMinimum(0.1);
  }  

  t[0]->Project(Form("hNVtx_%i",0),"nVtx","","");
  for(int i=1; i<5; ++i){
    t[i]->Project(Form("hNVtx_%i",1),"nVtx","weight","");
    hNVtx[1]->Scale(scale[i]);
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

}
