#include <math.h> 
#include <string.h>
#include <iostream>

#include "TROOT.h"
#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TAxis.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TFrame.h"

#include "tdrstyle.C" //sic! *.C included
#include "CMS_lumi.C" //sic! *.C included

//Implementation of CBErf from HTT TWiki by J.Swanson, R.Lane
double myCBErf(double m, double m0, double sigma, double alpha,
	     double n, double norm){
  //Useful constants
  const double sqrtPiOver2 = sqrt(TMath::PiOver2()); //1.2533141373;
  const double sqrt2 = sqrt(2.); //1.4142135624;

  
  double sig = fabs((double)sigma);
  double t = (m - m0)/sig;
  if(alpha < 0)
    t = -t;
  double absAlpha = fabs(alpha/sig);
  double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
  double b = absAlpha - n/absAlpha;
  double ApproxErf;
  double arg = absAlpha / sqrt2;
  if(arg > 5.) ApproxErf = 1.;
  else if(arg < -5.) ApproxErf = -1.;
  else ApproxErf = TMath::Erf(arg);
  double leftArea = (1. + ApproxErf) * sqrtPiOver2;
  double rightArea = ( a * 1./TMath::Power(absAlpha-b, n-1) ) / (n - 1);
  double area = leftArea + rightArea;
  if( t <= absAlpha ){
    arg = t / sqrt2;
    if(arg > 5.) ApproxErf = 1.;
    else if(arg < -5.) ApproxErf = -1.;
    else ApproxErf = TMath::Erf(arg);
    return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
  }
  else{
    return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) -
				   1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area;
  }
}
/////////////////
void CMSPrelim(const char* dataset, const char* channel, const char* preliminary, double lowX, double lowY)
{
  //TPaveText* cmsprel  = new TPaveText(lowX, lowY+0.06, lowX+0.30, lowY+0.16, "NDC");
  //TPaveText* cmsprel  = new TPaveText(lowX+0.08, lowY+0.061, lowX+0.45, lowY+0.161, "NDC");
  //v1TPaveText* cmsprel  = new TPaveText(lowX+0.07, lowY-0.043, lowX+0.45, lowY-0.023, "NDC");
  TPaveText* cmsprel  = new TPaveText(lowX+0.07, lowY-0.053, lowX+0.45, lowY-0.033, "NDC");
  cmsprel->SetBorderSize(   0 );
  cmsprel->SetFillStyle(    0 );
  cmsprel->SetTextAlign(   12 );
  //cmsprel->SetTextSize ( 0.03 );
  cmsprel->SetTextSize ( 0.04 );
  cmsprel->SetTextColor(    1 );
  cmsprel->SetTextFont (   62 );//52
  cmsprel->AddText("CMS");
  cmsprel->Draw();

  TPaveText* cmsprel1 = new TPaveText(lowX+0.16, lowY-0.057, lowX+0.54, lowY-0.037, "NDC");
  cmsprel1->SetBorderSize(   0 );
  cmsprel1->SetFillStyle(    0 );
  cmsprel1->SetTextAlign(   12 );
  //cmsprel1->SetTextSize ( 0.03 );
  cmsprel1->SetTextSize ( cmsprel->GetTextSize() );
  cmsprel1->SetTextColor(    1 );
  cmsprel1->SetTextFont (   52 );//62
  cmsprel1->AddText("Preliminary");
  cmsprel1->Draw();

  TPaveText* cmsprel2  = new TPaveText(cmsprel->GetX1(), cmsprel->GetY1()-1.1*cmsprel->GetTextSize(), 
				       cmsprel->GetX2(), cmsprel->GetY2()-1.1*cmsprel->GetTextSize(), "NDC");
  cmsprel2->SetBorderSize(   0 );
  cmsprel2->SetFillStyle(    0 );
  cmsprel2->SetTextAlign(   12 );
  //cmsprel2->SetTextSize ( 0.03 );
  cmsprel2->SetTextSize ( cmsprel->GetTextSize() );
  cmsprel2->SetTextColor(    1 );
  cmsprel2->SetTextFont (   52 );//62
  cmsprel2->AddText(preliminary); //
  cmsprel2->Draw();

  //TPaveText* lumi  = new TPaveText(lowX+0.08, lowY+0.061, lowX+0.45, lowY+0.161, "NDC");
  TPaveText* lumi  = new TPaveText(cmsprel->GetX1(), lowY, cmsprel->GetX2(), lowY+0.05, "NDC");
  lumi->SetBorderSize(   0 );
  lumi->SetFillStyle(    0 );
  lumi->SetTextAlign(   12 );
  //lumi->SetTextSize ( 0.03 );
  lumi->SetTextSize ( 0.04 );
  lumi->SetTextColor(    1 );
  lumi->SetTextFont (   62 );
  lumi->AddText(dataset);
  lumi->Draw();

  //TPaveText* chan     = new TPaveText(lowX+0.68, lowY+0.061, lowX+0.80, lowY+0.161, "NDC");
  TPaveText* chan     = new TPaveText(cmsprel->GetX1()+/*0.55*//*0.45*/0.35, cmsprel->GetY1(), 
				      cmsprel->GetX1()+0.90, cmsprel->GetY2(), "NDC");
  chan->SetBorderSize(   0 );
  chan->SetFillStyle(    0 );
  chan->SetTextAlign(   12 );
  //chan->SetTextSize ( 0.05 );
  //chan->SetTextSize ( 0.03 );
  chan->SetTextSize ( 0.04 );
  chan->SetTextColor(    1 );
  chan->SetTextFont (   62 );
  chan->AddText(channel);
  chan->Draw();
}
////////////////////////
TFile *fVtxW=TFile::Open("vtxNorm-RunC.root");
TH1D *hNVtxNorm=(TH1D*)fVtxW->Get("hNVtxNorm")->Clone();
//fVtxW->Close();
double vtxWeight(int nvtx, bool doIt=true){
  if(!doIt) return 1;
  int iBin=hNVtxNorm->GetXaxis()->FindBin(double(nvtx+0.5));
  return hNVtxNorm->GetBinContent(iBin);
}

////////////////////////

void eff(double lumi=1, /*pb-1*/
	 bool fit=false,
	 bool doVtxReWeight=true){

  //set style
  setTDRStyle();
  tdrGrid(false, gStyle);//switch off the grid
  gStyle->SetPadRightMargin(0.05);
  //gStyle->SetErrorX(0); //To not display x-error bars of TGraph
  //std::cout<<"ErrX: "<< gStyle->GetErrorX()<<std::endl;
  bool zeroErrX=false;
  gStyle->SetOptFit(0); //To not display fit results
  bool trueDYtt=false;  //To draw efficiency for true taus from DY

  TChain *t[5];
  //data
  t[0] = new TChain("muLooseTau/muTauTriggerTree");
  //t[0]->Add("effTrees/Run2015C_Cert_254231-254879_Eff_74X_v3/muTauTrgAna*.root"); //8.4pb-1 (after brilcalc tool)
  t[0]->Add("effTrees/Run2015C_Cert_246908-255031_Eff_74X_v3.1/muTauTrgAna*.root"); //16.1pb-1 (after brilcalc tool)
  t[0]->Add("effTrees/Run2015C_Cert_254833_Eff_74X_v3.1/muTauTrgAna*.root"); //23.2pb-1 (after brilcalc tool)
  //DY
  t[1] = new TChain("muLooseTau/muTauTriggerTree");
  t[1]->Add("effTrees/DYToLL_M50_Eff_74X_25ns_v3/muTauTrgAna*.root");
  //W
  t[2] = new TChain("muLooseTau/muTauTriggerTree");
  t[2]->Add("effTrees/WToLNu_Eff_74X_25ns_v3/muTauTrgAna*.root");
  //TTbar
  t[3] = new TChain("muLooseTau/muTauTriggerTree");
  t[3]->Add("effTrees/TTbar_Eff_74X_25ns_v3/muTauTrgAna*.root");
  //QCD
  t[4] = new TChain("muLooseTau/muTauTriggerTree");
  t[4]->Add("effTrees/QCD-Pt20-MuEnrPt15_Eff_74X_25ns_v3/muTauTrgAna*.root");

  double scale[5] = {
        1,     //data
     2008.4*3, //DY x-sec [pb]
    20508.9*3, //W x-sec [pb]
      831.76,  //TTbar x-sec [pb]
720648000*0.00042 //QCD 
  };
  TFile *f[4];
  f[0]=TFile::Open("effTrees/DYToLL_M50_Eff_74X_25ns_v3/DYToLL_M50Count.root");
  f[1]=TFile::Open("effTrees/WToLNu_Eff_74X_25ns_v3/WToLNuCount.root");
  f[2]=TFile::Open("effTrees/TTbar_Eff_74X_25ns_v3/TTbarCount.root");
  f[3]=TFile::Open("effTrees/QCD-Pt20-MuEnrPt15_Eff_74X_25ns_v3/QCD-Pt20-MuEnrPt15Count.root");
  for(unsigned int i=0; i<4; ++i){
    TH1D *h = (TH1D*)f[i]->Get("evtCounter/hCounts");
    double w = h->GetBinContent(2);
    scale[i+1] = lumi*scale[i+1]/w;
    f[i]->Close();
  }


  //double ptBins[]={0,2,4,6,8,10,12,14,16,18,20,22,24,28,30,35,40,45,50,55,60,70,80,100,120,140,170,200,250,300,350,400};
  //hPt = new TH1D("hPt"," ;p_{T}^{#tau offline} [GeV]; Events",24,ptBins);//upto 120
  double ptBins[]={0,2,4,6,8,10,12,14,16,18,20,22,26,30,35,40,50,70,100,150,200,250,300,350,400};
  TH1D *hPt[7];
  TH1D *hEta[7];
  TH1D *hMVis[7];
  TH1D *hMTau[7];
  TH1D *hMt[7];
  TH1D *hNVtx[7];
  for(unsigned int i=0; i<7; ++i){
    hPt[i] = new TH1D(Form("hPt_%i",i)," ;p_{T}^{#tau offline} [GeV]; Events/width [1/GeV]",18,ptBins);//upto 100
    hPt[i]->Sumw2();
    hPt[i]->SetMarkerStyle(20);
    //hPt[i]->SetMarkerSize(0.7);
    hPt[i]->SetStats(0);
    //hPt[i]->SetMinimum(0.1);

    hEta[i] = new TH1D(Form("hEta_%i",i)," ;#eta^{#tau offline}; Events",25,-2.5,2.5);
    hEta[i]->Sumw2();
    hEta[i]->SetMarkerStyle(20);
    //hEta[i]->SetMarkerSize(0.7);
    hEta[i]->SetStats(0);
    //hMVis[i]->SetMinimum(0.1);

    hMVis[i] = new TH1D(Form("hMVis_%i",i)," ;M_{vis}(#mu,#tau) [GeV]; Events",20,0,200);
    hMVis[i]->Sumw2();
    hMVis[i]->SetMarkerStyle(20);
    //hMVis[i]->SetMarkerSize(0.7);
    hMVis[i]->SetStats(0);
    //hMVis[i]->SetMinimum(0.1);

    hMTau[i] = new TH1D(Form("hMTau_%i",i)," ;m_{#tau} [GeV]; Events",18,0,1.8);
    hMTau[i]->Sumw2();
    hMTau[i]->SetMarkerStyle(20);
    //hMVis[i]->SetMarkerSize(0.7);
    hMTau[i]->SetStats(0);
    //hMTau[i]->SetMinimum(0.1);

    hMt[i] = new TH1D(Form("hMt_%i",i)," ;M_{T}(#mu,E_{T}^{miss}) [GeV]; Events",16,0,160);
    hMt[i]->Sumw2();
    hMt[i]->SetMarkerStyle(20);
    //hMt[i]->SetMarkerSize(0.7);
    hMt[i]->SetStats(0);
    //hMt[i]->SetMinimum(0.1);

    hNVtx[i] = new TH1D(Form("hNVtx_%i",i)," ;No. of vertices; Events",50,0,50);
    hNVtx[i]->Sumw2();
    hNVtx[i]->SetMarkerStyle(20);
    //hNVtx[i]->SetMarkerSize(0.7);
    hNVtx[i]->SetStats(0);
    //hNVtx[i]->SetMinimum(0.1);
  }  

  std::string muTag = "(HLT_IsoMu20_eta2p1_v>0 && hltL3crIsoL1sMu16Eta2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p09>0)";
  std::string probe1 = "hltOverlapFilterSingleIsoMu17LooseIsoPFTau20>0";
  std::string probe2 = "hltL1sMu16erTauJet20er>0 && hltOverlapFilterIsoMu17LooseIsoPFTau20>0";
  std::string probe20 = "hltL1sMu16erTauJet20er>0";
  std::string mu2veto = "!(mu2Pt>15&&mu2Iso<0.3&&muCharge*mu2Charge<0)";//||abs(mu2Eta)>2.4)";
  std::string mu3veto = /*"!(mu2Pt>10&&mu2Iso<0.3&&mu2Id>0)";*/"!(mu2Pt>10&&mu2Iso<0.3&&(mu2Id>10.9||(mu2Id>0.9&&mu2Id<1.1)))";
  std::string off   = "muPt>21 && abs(muEta)<2.1 && muIso<0.1 && tauPt>18 && abs(tauEta)<2.3 && muCharge*tauCharge<0 && decayModeFindingNewDMs>0.5 && abs(tauCharge)<1.001 && abs(tauCharge)>0.999 && byCombinedIsolationDeltaBetaCorrRaw3Hits<1.5 && againstElectronVLooseMVA5>0.5 && againstMuonTight3>0.5";
  std::string offSS = "muPt>21 && abs(muEta)<2.1 && muIso<0.1 && tauPt>18 && abs(tauEta)<2.3 && muCharge*tauCharge>0 && decayModeFindingNewDMs>0.5 && abs(tauCharge)<1.001 && abs(tauCharge)>0.999 && byCombinedIsolationDeltaBetaCorrRaw3Hits<1.5 && againstElectronVLooseMVA5>0.5 && againstMuonTight3>0.5";
  off   = off   + " && " + mu2veto + " && " + mu3veto;
  offSS = offSS + " && " + mu2veto + " && " + mu3veto;
  std::string window = /*"1";*/"Mass>40";//"(Mass>40 && Mass<80)";
  std::string sel = muTag + " && " + off + " && " + window;
  std::string mTcut = "Mt<40";
  std::string mTHcut = "Mt>80";//"Mt>70";
  std::string selMt = sel + " && " + mTcut;
  std::string selNoM = muTag + " && " + off  + " && " + mTcut;
  std::string vtxReWeight="*vtxWeight(nVtx)";
  if(!doVtxReWeight)
    vtxReWeight="*vtxWeight(nVtx,false)";

  std::string cut, cutMt, cutNoM;
  unsigned int ii=0;
  for(unsigned int i=0; i<6; ++i){
    if(i>5){
      std::cout<<"What?!"<<std::endl;
      continue;
    }
    cut = "weight*("+sel+")";
    cutMt = "weight*("+selMt+")";
    cutNoM = "weight*("+selNoM+")";
    ii=i;
    if(i==1){ //DY->tt
      cut = "weight*("+sel+" && tauGenMatch>0.001"+")";
      cutMt = "weight*("+selMt+" && tauGenMatch>0.001"+")";
      cutNoM = "weight*("+selNoM+" && tauGenMatch>0.001"+")";
    }
    else if(i==4){ //DY fake-tau
      ii=1;
      cut = "weight*("+sel+" && !(tauGenMatch>0.001)"+")";
      cutMt = "weight*("+selMt+" && !(tauGenMatch>0.001)"+")";
      cutNoM = "weight*("+selNoM+" && !(tauGenMatch>0.001)"+")";
    }
    else if(i>4){ //DY fake-tau
      ii=i-1;
    }
    if(i!=0){
      cut = cut+vtxReWeight;
      cutMt = cutMt+vtxReWeight;
      cutNoM = cutNoM+vtxReWeight;
    }
    std::cout<<i<<", "<<ii<<std::endl;
    std::cout<<cut<<std::endl;

    t[ii]->Project(Form("hPt_%i",i),"tauPt",cutMt.c_str(),"");
    std::cout<<"Expected: "<<hPt[i]->Integral(0,hPt[i]->GetNbinsX()+1)*scale[ii]
	     <<" for L="<<lumi<<"pb-1"<<std::endl;
    hPt[i]->Scale(scale[ii],"width");
    std::cout<<hPt[i]->Integral(0,hPt[i]->GetNbinsX()+1,"width")<<std::endl;
    //hPt[i]->SetMinimum(0.1);
    //t[ii]->Project(Form("hMVis_%i",i),"Mass",cutMt.c_str(),"");
    t[ii]->Project(Form("hMVis_%i",i),"Mass",cutNoM.c_str(),"");
    //hMVis[i]->Scale(scale[ii],"width");
    hMVis[i]->Scale(scale[ii],"");
    //hMVis[i]->SetMinimum(0.1);
    t[ii]->Project(Form("hEta_%i",i),"tauEta",cutMt.c_str(),"");
    //hEta[i]->Scale(scale[ii],"width");
    hEta[i]->Scale(scale[ii],"");
    std::cout<<hEta[i]->Integral(0,hEta[i]->GetNbinsX()+1)<<std::endl;
    //hEta[i]->SetMinimum(0.1);
    t[ii]->Project(Form("hMTau_%i",i),"tauM",cutMt.c_str(),"");
    //hMTau[i]->Scale(scale[ii],"width");
    hMTau[i]->Scale(scale[ii],"");
    //hMTau[i]->SetMinimum(0.1);
    t[ii]->Project(Form("hMt_%i",i),"Mt",cut.c_str(),"");
    //hMt[i]->Scale(scale[ii],"width");
    hMt[i]->Scale(scale[ii],"");
    //hMt[i]->SetMinimum(0.1);
    t[ii]->Project(Form("hNVtx_%i",i),"nVtx",cut.c_str(),"");
    //hNVtx[i]->Scale(scale[ii],"width");
    hNVtx[i]->Scale(scale[ii]);
    //hNVtx[i]->SetMinimum(0.1);
    std::cout<<"end of loop: "<<i<<std::endl;
  }
  //W normalization Mt>70
  std::cout<<"High Mt for W\n";
  std::string cutHMt, cutHMtSS;//kuku
  cutHMt     = "weight*(" + muTag + " && " + off   + " && " + mTHcut + " && " + window +")";
  cutHMtSS   = "weight*(" + muTag + " && " + offSS + " && " + mTHcut + " && " + window +")";
  TH1D *hEtaMtHighTmp=(TH1D*)hEta[5]->Clone("hEtaMtHighTmp");
  hEtaMtHighTmp->Reset();
  TH1D *hEtaMtHigh[2];
  for(int i=0; i<2; ++i){
    hEtaMtHigh[i]=(TH1D*)hEta[5]->Clone(Form("hEtaMtHigh_%i",i));
    hEtaMtHigh[i]->Reset();  
  }
  //high-Mt OS
  //data
  std::cout<<cutHMt<<std::endl;
  hEtaMtHighTmp->Reset();
  hEtaMtHigh[0]->Reset();
  t[0]->Project("hEtaMtHighTmp","tauEta",cutHMt.c_str(),"");
  std::cout<<"High Mt :"<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1);
  std::cout<<"+"<<hEtaMtHighTmp->Integral(0,hEtaMtHighTmp->GetNbinsX()+1);
  hEtaMtHigh[0]->Add(hEtaMtHighTmp,1);
  std::cout<<"="<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1)<<std::endl;
  //DY
  cutHMt = cutHMt+vtxReWeight;
  hEtaMtHighTmp->Reset();
  t[1]->Project("hEtaMtHighTmp","tauEta",cutHMt.c_str(),"");
  std::cout<<"High Mt :"<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1);
  hEtaMtHighTmp->Scale(scale[1]);
  std::cout<<"-"<<hEtaMtHighTmp->Integral(0,hEtaMtHighTmp->GetNbinsX()+1);
  hEtaMtHigh[0]->Add(hEtaMtHighTmp,-1);
  std::cout<<"="<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1)<<std::endl;
  //TT
  hEtaMtHighTmp->Reset();
  t[3]->Project("hEtaMtHighTmp","tauEta",cutHMt.c_str(),"");
  std::cout<<"High Mt :"<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1);
  hEtaMtHighTmp->Scale(scale[3]);
  std::cout<<"-"<<hEtaMtHighTmp->Integral(0,hEtaMtHighTmp->GetNbinsX()+1);
  hEtaMtHigh[0]->Add(hEtaMtHighTmp,-1);
  std::cout<<"="<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1)<<std::endl;
  //W
  hEtaMtHigh[1]->Reset();
  t[2]->Project("hEtaMtHigh_1","tauEta",cutHMt.c_str(),"");
  hEtaMtHigh[1]->Scale(scale[2]);
  double wNorm=hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1)/hEtaMtHigh[1]->Integral(0,hEtaMtHigh[1]->GetNbinsX()+1);
  std::cout<<"wNorm: "<<wNorm
	   <<"="<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1)
	   <<"/"<<hEtaMtHigh[1]->Integral(0,hEtaMtHigh[1]->GetNbinsX()+1)
	   <<std::endl;
  hPt[2]->Scale(wNorm);
  hEta[2]->Scale(wNorm);
  hMVis[2]->Scale(wNorm);
  hMTau[2]->Scale(wNorm);
  hMt[2]->Scale(wNorm);
  hNVtx[2]->Scale(wNorm);
  //high-Mt SS
  //data
  std::cout<<cutHMtSS<<std::endl;
  hEtaMtHighTmp->Reset();
  hEtaMtHigh[0]->Reset();
  t[0]->Project("hEtaMtHighTmp","tauEta",cutHMtSS.c_str(),"");
  std::cout<<"High Mt :"<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1);
  std::cout<<"+"<<hEtaMtHighTmp->Integral(0,hEtaMtHighTmp->GetNbinsX()+1);
  hEtaMtHigh[0]->Add(hEtaMtHighTmp,1);
  std::cout<<"="<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1)<<std::endl;
  //DY
  cutHMtSS = cutHMtSS+vtxReWeight;
  hEtaMtHighTmp->Reset();
  t[1]->Project("hEtaMtHighTmp","tauEta",cutHMtSS.c_str(),"");
  std::cout<<"High Mt :"<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1);
  hEtaMtHighTmp->Scale(scale[1]);
  std::cout<<"-"<<hEtaMtHighTmp->Integral(0,hEtaMtHighTmp->GetNbinsX()+1);
  hEtaMtHigh[0]->Add(hEtaMtHighTmp,-1);
  std::cout<<"="<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1)<<std::endl;
  //TT
  hEtaMtHighTmp->Reset();
  t[3]->Project("hEtaMtHighTmp","tauEta",cutHMtSS.c_str(),"");
  std::cout<<"High Mt :"<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1);
  hEtaMtHighTmp->Scale(scale[3]);
  std::cout<<"-"<<hEtaMtHighTmp->Integral(0,hEtaMtHighTmp->GetNbinsX()+1);
  hEtaMtHigh[0]->Add(hEtaMtHighTmp,-1);
  std::cout<<"="<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1)<<std::endl;
  //W
  hEtaMtHigh[1]->Reset();
  t[2]->Project("hEtaMtHigh_1","tauEta",cutHMtSS.c_str(),"");
  hEtaMtHigh[1]->Scale(scale[2]);
  double wNormSS=hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1)/hEtaMtHigh[1]->Integral(0,hEtaMtHigh[1]->GetNbinsX()+1);
  std::cout<<"wNormSS: "<<wNormSS
	   <<"="<<hEtaMtHigh[0]->Integral(0,hEtaMtHigh[0]->GetNbinsX()+1)
	   <<"/"<<hEtaMtHigh[1]->Integral(0,hEtaMtHigh[1]->GetNbinsX()+1)
	   <<std::endl;

  //QCD
  std::cout<<"QCD from SS\n";
  std::string cutSS, cutMtSS, cutNoMSS;
  TH1D *hPtTmp=(TH1D*)hPt[6]->Clone("hPtTmp");
  TH1F *hPtSS[4];
  for(int i=0; i<4; ++i){
    hPtSS[i] = (TH1F*)hPt[6]->Clone(Form("hPtSS_%i",i));
    hPtSS[i]->Reset();
  }
  TH1D *hEtaTmp=(TH1D*)hEta[6]->Clone("hEtaTmp");
  TH1D *hMVisTmp=(TH1D*)hMVis[6]->Clone("hMVisTmp");
  TH1D *hMTauTmp=(TH1D*)hMTau[6]->Clone("hMTauTmp");
  TH1D *hMtTmp=(TH1D*)hMt[6]->Clone("hMtTmp");
  TH1D *hNVtxTmp=(TH1D*)hNVtx[6]->Clone("hNVtxTmp");
  for(unsigned int i=0; i<4; ++i){
    cutSS    = "weight*(" + muTag + " && " + offSS + " && " + window +")";
    cutMtSS  = "weight*(" + muTag + " && " + offSS + " && " + window + " && " + mTcut +")";
    cutNoMSS = "weight*(" + muTag + " && " + offSS +                   " && " + mTcut +")";
    std::cout<<cutSS<<std::endl;
    if(i!=0){
      cutSS = cutSS+vtxReWeight;
      cutMtSS = cutMtSS+vtxReWeight;
      cutNoMSS = cutNoMSS+vtxReWeight;
    }
    hPtTmp->Reset();
    t[i]->Project("hPtTmp","tauPt",cutMtSS.c_str(),"");
    hPtTmp->Scale(scale[i],"width");
    if(i==2)//W
      hPtTmp->Scale(wNormSS,"");
    std::cout<<i<<". QCD Pt :"<<hPt[6]->Integral(0,hPt[6]->GetNbinsX()+1,"width");
    if(i==0){
      hPt[6]->Add(hPtTmp,1);
      hPtSS[0]->Add(hPtTmp,1);
      std::cout<<" + ";
    }
    else{
      hPt[6]->Add(hPtTmp,-1);
      std::cout<<" - ";
    }
    std::cout<<hPtTmp->Integral(0,hPtTmp->GetNbinsX()+1,"width");
    std::cout<<" = "<<hPt[6]->Integral(0,hPt[6]->GetNbinsX()+1,"width")<<std::endl;

    hEtaTmp->Reset();
    t[i]->Project("hEtaTmp","tauEta",cutMtSS.c_str(),"");
    //hEtaTmp->Scale(scale[i],"width");
    hEtaTmp->Scale(scale[i],"");
    if(i==2)//W
      hEtaTmp->Scale(wNormSS,"");
    std::cout<<i<<". QCD :"<<hEta[6]->Integral(0,hEta[6]->GetNbinsX()+1);
    if(i==0){
      hEta[6]->Add(hEtaTmp,1);
      std::cout<<" + ";
    }
    else{
      hEta[6]->Add(hEtaTmp,-1);
      std::cout<<" - ";
    }
    std::cout<<hEtaTmp->Integral(0,hEtaTmp->GetNbinsX()+1);
    std::cout<<" = "<<hEta[6]->Integral(0,hEta[6]->GetNbinsX()+1)<<std::endl;

    hMVisTmp->Reset();
    t[i]->Project("hMVisTmp","Mass",cutNoMSS.c_str(),"");
    //hMVisTmp->Scale(scale[i],"width");
    hMVisTmp->Scale(scale[i],"");
    if(i==2)//W
      hMVisTmp->Scale(wNormSS,"");
    if(i==0){
      hMVis[6]->Add(hMVisTmp,1);
    }
    else{
      hMVis[6]->Add(hMVisTmp,-1);
    }

    hMTauTmp->Reset();
    t[i]->Project("hMTauTmp","tauM",cutMtSS.c_str(),"");
    //hMTauTmp->Scale(scale[i],"width");
    hMTauTmp->Scale(scale[i],"");
    if(i==2)//W
      hMTauTmp->Scale(wNormSS,"");
    if(i==0)
      hMTau[6]->Add(hMTauTmp,1);
    else
      hMTau[6]->Add(hMTauTmp,-1);

    hMtTmp->Reset();
    t[i]->Project("hMtTmp","Mt",cutSS.c_str(),"");
    //hMtTmp->Scale(scale[i],"width");
    hMtTmp->Scale(scale[i],"");
    if(i==2)//W
      hMtTmp->Scale(wNormSS,"");
    if(i==0)
      hMt[6]->Add(hMtTmp,1);
    else
      hMt[6]->Add(hMtTmp,-1);

    hNVtxTmp->Reset();
    t[i]->Project("hNVtxTmp","nVtx",cutSS.c_str(),"");
    //hNVtxTmp->Scale(scale[i],"width");
    hNVtxTmp->Scale(scale[i]);
    if(i==2)//W
      hNVtxTmp->Scale(wNormSS,"");
    if(i==0)
      hNVtx[6]->Add(hNVtxTmp,1);
    else
      hNVtx[6]->Add(hNVtxTmp,-1);
  }
  hPt[6]->Scale(1.06);
  hEta[6]->Scale(1.06);
  hMVis[6]->Scale(1.06);
  hMTau[6]->Scale(1.06);
  hMt[6]->Scale(1.06);
  hNVtx[6]->Scale(1.06);

  double QCDnorm=hEta[6]->Integral(0,hEta[6]->GetNbinsX()+1)/hEta[5]->Integral(0,hEta[5]->GetNbinsX()+1);
  double QCDnormNoM=hMVis[6]->Integral(0,hMVis[6]->GetNbinsX()+1)/hMVis[5]->Integral(0,hMVis[5]->GetNbinsX()+1);
  std::cout<<"QCDnorm "<<QCDnorm
	   <<"="<<hEta[6]->Integral(0,hEta[6]->GetNbinsX()+1,"")
	   <<"/"<<hEta[5]->Integral(0,hEta[5]->GetNbinsX()+1,"")
	   <<std::endl;
  std::cout<<"QCDnormNoM "<<QCDnormNoM<<std::endl;
  hPt[5]->Scale(QCDnorm);
  hEta[5]->Scale(QCDnorm);
  hMVis[5]->Scale(QCDnormNoM);
  hMTau[5]->Scale(QCDnorm);
  hMt[5]->Scale(QCDnorm);
  hNVtx[5]->Scale(QCDnorm);
  ////

  //stacked
  int colMap[] = {kOrange-2,kRed+2,kBlue-1,kGreen+3,kMagenta-10};
  TH1D *hSPt[5];
  TH1D *hSEta[5];
  TH1D *hSMVis[5];
  TH1D *hSMTau[5];
  TH1D *hSMt[5];
  TH1D *hSNVtx[5];
  for(unsigned int i=0; i<5; ++i){
    hSPt[i]=(TH1D*)hPt[1]->Clone(Form("hSPt_%i",i));
    hSPt[i]->Reset();
    hSPt[i]->SetFillColor(colMap[i]);
    hSEta[i]=(TH1D*)hEta[1]->Clone(Form("hSEta_%i",i));
    hSEta[i]->Reset();
    hSEta[i]->SetFillColor(colMap[i]);
    hSMVis[i]=(TH1D*)hMVis[1]->Clone(Form("hSMVis_%i",i));
    hSMVis[i]->Reset();
    hSMVis[i]->SetFillColor(colMap[i]);
    hSMTau[i]=(TH1D*)hMTau[1]->Clone(Form("hSMTau_%i",i));
    hSMTau[i]->Reset();
    hSMTau[i]->SetFillColor(colMap[i]);
    hSMt[i]=(TH1D*)hMt[1]->Clone(Form("hSMt_%i",i));
    hSMt[i]->Reset();
    hSMt[i]->SetFillColor(colMap[i]);
    hSNVtx[i]=(TH1D*)hNVtx[1]->Clone(Form("hSNVtx_%i",i));
    hSNVtx[i]->Reset();
    hSNVtx[i]->SetFillColor(colMap[i]);
    for(unsigned int j=i+1; j<6; ++j){
      hSPt[i]->Add(hPt[j]);
      hSEta[i]->Add(hEta[j]);
      hSMVis[i]->Add(hMVis[j]);
      hSMTau[i]->Add(hMTau[j]);
      hSMt[i]->Add(hMt[j]);
      hSNVtx[i]->Add(hNVtx[j]);
    }
  }
  std::cout<<"Data: "<<hEta[0]->Integral(0,hEta[0]->GetNbinsX()+1)
	   <<", MC: "<<hSEta[0]->Integral(0,hSEta[0]->GetNbinsX()+1)<<std::endl;
  std::cout<<"Data (NoM): "<<hMVis[0]->Integral(0,hMVis[0]->GetNbinsX()+1)
	   <<", MC (NoM): "<<hSMVis[0]->Integral(0,hSMVis[0]->GetNbinsX()+1)<<std::endl;

  ///Effs
  TH1D *hhPt[3][3];
  TH1D *hhPtEff[3][3];
  TGraphAsymmErrors *grPtEff[3][3];
  TH1D *hPtSSEff[3];
  TGraphAsymmErrors *grPtSSEff[3];
  //TH2D *hFrame = new TH2D("hFrame",";p_{T}^{#tau offline} [GeV]; Efficiency",2,0,100,2,0.01,1.10);
  //hFrame->SetStats(0);

  for(unsigned int i=0; i<3; ++i){
    for(unsigned int j=0; j<3; ++j){
      hhPt[i][j]=(TH1D*)hPt[1]->Clone(Form("hhPt_%i%i",i,j)); 
      //hhPt[i][j]->Sumw2();
      hhPt[i][j]->Reset();     

      hhPtEff[i][j]=(TH1D*)hPt[1]->Clone(Form("hhPtEff_%i%i",i,j)); 
      //hhPtEff[i][j]->Sumw2();
      hhPtEff[i][j]->SetTitle(";p_{T}^{#tau offline} [GeV]; Efficiency");
      hhPtEff[i][j]->Reset();       
    }
    hPtSSEff[i]=(TH1D*)hPtSS[0]->Clone(Form("hPtSSEff_%i",i)); 
    //hhPtEff[i]->Sumw2();
    hPtSSEff[i]->SetTitle(";p_{T}^{#tau offline} [GeV]; Efficiency");
    hPtSSEff[i]->Reset();       
  }
  //Data and all MC's
  cut = "weight*("+probe1+"&&"+selMt+")";
  std::cout<<cut<<std::endl;
  t[0]->Project("hhPt_00","tauPt",cut.c_str(),"");
  cut = cut+vtxReWeight;
  for(unsigned int i=1; i<5; ++i){
    hPtTmp->Reset();
    t[i]->Project("hPtTmp","tauPt",cut.c_str(),"");
    hPtTmp->Scale(scale[i],"width");
    if(i==2)//W
      hPtTmp->Scale(wNorm,"");
    else if(i==4)//QCD
      hPtTmp->Scale(QCDnorm,"");
    hhPt[0][2]->Add(hPtTmp);
    std::cout<<hhPt[0][2]->GetName()<<": "
	     <<hhPt[0][2]->Integral(0,hhPt[0][2]->GetNbinsX()+1)<<std::endl;
  }
  cut = "weight*("+probe2+"&&"+selMt+")";
  t[0]->Project("hhPt_10","tauPt",cut.c_str(),"");
  cut = cut+vtxReWeight;
  for(unsigned int i=1; i<5; ++i){
    hPtTmp->Reset();
    t[i]->Project("hPtTmp","tauPt",cut.c_str(),"");
    hPtTmp->Scale(scale[i],"width");
    if(i==2)//W
      hPtTmp->Scale(wNorm,"");
    else if(i==4)//QCD
      hPtTmp->Scale(QCDnorm,"");
    hhPt[1][2]->Add(hPtTmp);
    std::cout<<hhPt[1][2]->GetName()<<": "
	     <<hhPt[1][2]->Integral(0,hhPt[1][2]->GetNbinsX()+1)<<std::endl;
  }
  cut = "weight*("+probe20+"&&"+selMt+")";
  t[0]->Project("hhPt_20","tauPt",cut.c_str(),"");
  cut = cut+vtxReWeight;
  for(unsigned int i=1; i<5; ++i){
    hPtTmp->Reset();
    t[i]->Project("hPtTmp","tauPt",cut.c_str(),"");
    hPtTmp->Scale(scale[i],"width");
    if(i==2)//W
      hPtTmp->Scale(wNorm,"");
    else if(i==4)//QCD
      hPtTmp->Scale(QCDnorm,"");
    hhPt[2][2]->Add(hPtTmp);
    std::cout<<hhPt[2][2]->GetName()<<": "
	     <<hhPt[2][2]->Integral(0,hhPt[2][2]->GetNbinsX()+1)<<std::endl;
  }
  //matched DY
  cut = "weight*("+probe1+"&&"+selMt+" && tauGenMatch>0"+")";
  cut = cut+vtxReWeight;
  std::cout<<cut<<std::endl;
  t[1]->Project("hhPt_01","tauPt",cut.c_str(),"");
  std::cout<<hhPt[0][1]->GetName()<<": "
	   <<hhPt[0][1]->Integral(0,hhPt[0][1]->GetNbinsX()+1)<<std::endl;
  cut = "weight*("+probe2+"&&"+selMt+" && tauGenMatch>0"+")";
  cut = cut+vtxReWeight;
  t[1]->Project("hhPt_11","tauPt",cut.c_str(),"");
  cut = "weight*("+probe20+"&&"+selMt+" && tauGenMatch>0"+")";
  cut = cut+vtxReWeight;
  t[1]->Project("hhPt_21","tauPt",cut.c_str(),"");
  //SS data
  cut = "weight*("+probe1+"&&"+muTag+" && "+offSS+" && "+window+" && "+mTcut +")";
  std::cout<<cut<<std::endl;
  t[0]->Project("hPtSS_1","tauPt",cut.c_str(),"");
  std::cout<<hPtSS[1]->GetName()<<": "
	   <<hPtSS[1]->Integral(0,hPtSS[1]->GetNbinsX()+1)<<std::endl;
  cut = "weight*("+probe2+"&&"+muTag+" && "+offSS+" && "+window+" && "+mTcut +")";
  t[0]->Project("hPtSS_2","tauPt",cut.c_str(),"");
  cut = "weight*("+probe20+"&&"+muTag+" && "+offSS+" && "+window+" && "+mTcut +")";
  t[0]->Project("hPtSS_3","tauPt",cut.c_str(),"");


  TH1D *hPtBkg = (TH1D*)hPt[1]->Clone("hPtBkg");
  //hPtBkg->Sumw2();
  hPtBkg->Reset();
  hPtBkg->Add(hPt[1]);
  hPtBkg->Add(hPt[2]);
  hPtBkg->Add(hPt[3]);
  hPtBkg->Add(hPt[4]);
  hPtBkg->Add(hPt[5]);//QCD from MC

  /*
  //emulate QCD with Wjets
  //simple scaling does not work well as Pt for Wjets harder than for QCD 
  //so trigger matched distribution is weighted bin-by-bin
  TH1D* hQCDNorm=(TH1D*)hPt[2]->Clone("hQCDNorm");  
  for(int ib=0;ib<hPt[5]->GetNbinsX()+2;++ib){
    if(hPt[2]->GetBinContent(ib)>0){
      hQCDNorm->SetBinContent(ib,hPt[5]->GetBinContent(ib)/hPt[2]->GetBinContent(ib));
    }
    else{
      hQCDNorm->SetBinContent(ib,0);
    }
  }
  hPtBkg->Add(hPt[5]);
  cut = "weight*("+probe1+"&&"+selMt+")";
  cut = cut+vtxReWeight;
  hPtTmp->Reset();
  t[2]->Project("hPtTmp","tauPt",cut.c_str(),"");
  hPtTmp->Scale(scale[2],"width");
  for(int ib=0;ib<hPtTmp->GetNbinsX()+2;++ib){
    hPtTmp->SetBinContent(ib,hPtTmp->GetBinContent(ib)*hQCDNorm->GetBinContent(ib));
    hPtTmp->SetBinError(ib,hPtTmp->GetBinError(ib)*hQCDNorm->GetBinContent(ib));
  }
  hhPt[0][2]->Add(hPtTmp);
  cut = "weight*("+probe2+"&&"+selMt+")";
  cut = cut+vtxReWeight;
  hPtTmp->Reset();
  t[2]->Project("hPtTmp","tauPt",cut.c_str(),"");
  hPtTmp->Scale(scale[2],"width");
  for(int ib=0;ib<hPtTmp->GetNbinsX()+2;++ib){
    hPtTmp->SetBinContent(ib,hPtTmp->GetBinContent(ib)*hQCDNorm->GetBinContent(ib));
    hPtTmp->SetBinError(ib,hPtTmp->GetBinError(ib)*hQCDNorm->GetBinContent(ib));
  }
  hhPt[1][2]->Add(hPtTmp);
  cut = "weight*("+probe20+"&&"+selMt+")";
  cut = cut+vtxReWeight;
  hPtTmp->Reset();
  t[2]->Project("hPtTmp","tauPt",cut.c_str(),"");
  hPtTmp->Scale(scale[2],"width");
  for(int ib=0;ib<hPtTmp->GetNbinsX()+2;++ib){
    hPtTmp->SetBinContent(ib,hPtTmp->GetBinContent(ib)*hQCDNorm->GetBinContent(ib));
    hPtTmp->SetBinError(ib,hPtTmp->GetBinError(ib)*hQCDNorm->GetBinContent(ib));
  }
  hhPt[2][2]->Add(hPtTmp);
  */


  for(unsigned int i=0; i<3; ++i){
    //Data
    hhPt[i][0]->Scale(scale[0],"width");
    //reset underflows
    hhPt[i][0]->SetBinContent(hhPt[i][0]->GetNbinsX()+1,0);
    hhPt[i][0]->SetBinError(hhPt[i][0]->GetNbinsX()+1,0);
    hPt[0]->SetBinContent(hPt[0]->GetNbinsX()+1,0);
    hPt[0]->SetBinError(hPt[0]->GetNbinsX()+1,0);
    hhPtEff[i][0]->Divide(hhPt[i][0],hPt[0],1,1,"b"); 
    grPtEff[i][0] = new TGraphAsymmErrors(hhPt[i][0],hPt[0],"cl=0.683 b(1,1) mode");
    hhPtEff[i][0]->SetMinimum(0.01); 
    hhPtEff[i][0]->SetMaximum(1.10);
    hhPtEff[i][0]->SetLineColor(kBlack);
    hhPtEff[i][0]->SetMarkerColor(kBlack);
    grPtEff[i][0]->SetLineColor(kBlack);
    grPtEff[i][0]->SetMarkerColor(kBlack);
    grPtEff[i][0]->SetMarkerStyle(grPtEff[i][0]->GetMarkerStyle());

    //SS data
    hPtSS[i+1]->Scale(scale[0],"width");
    //reset underflows
    hPtSS[0]->SetBinContent(hPtSS[0]->GetNbinsX()+1,0);
    hPtSS[0]->SetBinError(hPtSS[0]->GetNbinsX()+1,0);    
    hPtSS[i+1]->SetBinContent(hPtSS[i+1]->GetNbinsX()+1,0);
    hPtSS[i+1]->SetBinError(hPtSS[i+1]->GetNbinsX()+1,0);
    hPtSSEff[i]->Divide(hPtSS[i+1],hPtSS[0],1,1,"b"); 
    hPtSSEff[i]->SetMinimum(0.01); 
    hPtSSEff[i]->SetMaximum(1.10);
    grPtSSEff[i] = new TGraphAsymmErrors(hPtSS[i+1],hPtSS[0],"cl=0.683 b(1,1) mode");
    hPtSSEff[i]->SetLineColor(kRed);
    hPtSSEff[i]->SetMarkerColor(kRed);
    grPtSSEff[i]->SetLineColor(kRed);
    grPtSSEff[i]->SetMarkerColor(kRed);
    grPtSSEff[i]->SetMarkerStyle(grPtEff[i][0]->GetMarkerStyle());

    //matched DY
    hhPt[i][1]->Scale(scale[1],"width"); 
    //reset underflows
    hhPt[i][1]->SetBinContent(hhPt[i][1]->GetNbinsX()+1,0);
    hhPt[i][1]->SetBinError(hhPt[i][1]->GetNbinsX()+1,0);
    hPt[1]->SetBinContent(hPt[1]->GetNbinsX()+1,0);
    hPt[1]->SetBinError(hPt[1]->GetNbinsX()+1,0);
    std::cout<<hhPt[i][1]->GetName()<<": "
	     <<hhPt[i][1]->Integral(0,hhPt[i][1]->GetNbinsX()+1,"weight")<<std::endl;
    //
    hhPtEff[i][1]->Divide(hhPt[i][1],hPt[1],1,1,"b"); 
    grPtEff[i][1] = new TGraphAsymmErrors(hhPt[i][1],hPt[1],"cl=0.683 b(1,1) mode");
    hhPtEff[i][1]->SetMinimum(0.01); 
    hhPtEff[i][1]->SetMaximum(1.10);
    hhPtEff[i][1]->SetLineColor(kGreen+3);
    hhPtEff[i][1]->SetMarkerColor(kGreen+3);
    hhPtEff[i][1]->SetMarkerStyle(21);
    grPtEff[i][1]->SetLineColor(kGreen+3);
    grPtEff[i][1]->SetMarkerColor(kGreen+3);
    grPtEff[i][1]->SetMarkerStyle(grPtEff[i][1]->GetMarkerStyle());
    //all MC (QCD ommited)
    //hhPt[i][2]->Scale(xxx,"width");// already scaled 
    //reset underflows
    hhPt[i][2]->SetBinContent(hhPt[i][2]->GetNbinsX()+1,0);
    hhPt[i][2]->SetBinError(hhPt[i][2]->GetNbinsX()+1,0);
    hPtBkg->SetBinContent(hPtBkg->GetNbinsX()+1,0);
    hPtBkg->SetBinError(hPtBkg->GetNbinsX()+1,0);
    for(int ij=0; ij<hPtBkg->GetNbinsX()+2; ++ij){
      double ratio = hPtBkg->GetBinContent(ij)!=0. ? hhPt[i][2]->GetBinContent(ij)/hPtBkg->GetBinContent(ij) : -100.*hhPt[i][2]->GetBinContent(ij);
      //std::cout<<"hhPt["<<i<<"][2]/hPtBkg("<<ij<<") "<<hhPt[i][2]->GetBinContent(ij)<<"/"
      //	       <<hPtBkg->GetBinContent(ij)<<" = "<<ratio<<std::endl;
      if(ratio>1) {	
	std::cout<<" >1 !"<<std::endl;
	std::cout<<"hhPt["<<i<<"][2]/hPtBkg("<<ij<<") "<<hhPt[i][2]->GetBinContent(ij)<<"/"
		 <<hPtBkg->GetBinContent(ij)<<" = "<<ratio<<std::endl;
	if(ratio-1<0.005)//add hoc 'renormalization' for safety
	  hhPt[i][2]->SetBinContent(ij,hPtBkg->GetBinContent(ij));
      }
    }
    hhPtEff[i][2]->Divide(hhPt[i][2],hPtBkg,1,1,"b"); 
    grPtEff[i][2] = new TGraphAsymmErrors(hhPt[i][2],hPtBkg,"cl=0.683 b(1,1) mode");
    hhPtEff[i][2]->SetMinimum(0.01); 
    hhPtEff[i][2]->SetMaximum(1.10);
    hhPtEff[i][2]->SetLineColor(kRed);
    hhPtEff[i][2]->SetMarkerColor(kRed);
    hhPtEff[i][2]->SetMarkerStyle(21);
    grPtEff[i][2]->SetLineColor(kRed);
    grPtEff[i][2]->SetMarkerColor(kRed);
    grPtEff[i][2]->SetMarkerStyle(grPtEff[i][2]->GetMarkerStyle());
  }

  TCanvas *ca = new TCanvas("ca","myEfficiency",600,600);

  TLegend *leg = new TLegend(0.65,0.53,0.9,0.81);
  TLegend *leg2 = new TLegend(0.65,0.27,0.9,0.40);
  //leg2->SetTextSize(0.03);
  leg2->SetTextSize(0.04);
  leg2->SetEntrySeparation(0.0);
  leg2->SetColumnSeparation(0.05);
  //leg->SetFillColor(kWhite);
  //leg->SetLineColor(kWhite);
  ////leg->SetFillStyle(4000);
  leg->SetFillStyle (0);
  leg->SetFillColor (0);
  leg->SetBorderSize(0);
  leg2->SetFillStyle (0);
  leg2->SetFillColor (0);
  leg2->SetBorderSize(0);

  float sizeScale  = 1.3;
  lumiTextSize     = 0.6*sizeScale;
  lumiTextOffset   = 0.2*sizeScale;
  cmsTextSize      = 0.75*sizeScale;
  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"
  lumi_8TeV  = "19.7 fb^{-1}"; // default is "19.7 fb^{-1}"
  lumi_7TeV  = "4.9 fb^{-1}";  // default is "5.1 fb^{-1}"
  lumi_13TeV = Form("%.1f pb^{-1}",lumi);
  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 4=13TeV, 7=7+8+13TeV 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
  int iPos=11;// top-left, left-aligned
  //int iPos=33;// top-right, right-aligned
  //int iPos=22; center, centered

  hPt[0]->Draw();
  hSPt[0]->Draw("hist same");
  for(unsigned int i=1; i<5; ++i)
    hSPt[i]->Draw("same hist");
  hPt[0]->Draw("same");
  leg->AddEntry(hPt[0], " Data","PLE");
  leg->AddEntry(hSPt[0]," Z#rightarrow#tau#tau","F");
  leg->AddEntry(hSPt[1]," W+jets","F");
  leg->AddEntry(hSPt[2]," t#bar{t}","F");
  leg->AddEntry(hSPt[3]," Z#rightarrow#mu#mu","F");
  leg->AddEntry(hSPt[4]," QCD","F");
  leg->Draw();

  CMS_lumi(ca, iPeriod, iPos);
  ca->Update();
  ca->RedrawAxis();
  ca->GetFrame()->Draw();

  ca->SaveAs((std::string("tauEff32-RunC_Pt")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("tauEff32-RunC_Pt")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("tauEff32-RunC_Pt")+std::string(".png")).c_str());

  hEta[0]->Draw();
  hSEta[0]->Draw("hist same");
  for(unsigned int i=1; i<5; ++i)
    hSEta[i]->Draw("same hist");
  hEta[0]->Draw("same");
  leg->Draw();

  CMS_lumi(ca, iPeriod, iPos);
  ca->Update();
  ca->RedrawAxis();
  ca->GetFrame()->Draw();

  ca->SaveAs((std::string("tauEff32-RunC_Eta")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("tauEff32-RunC_Eta")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("tauEff32-RunC_Eta")+std::string(".png")).c_str());

  //hMVis[0]->Draw();
  hSMVis[0]->Draw("hist");//hMVis[0]->Draw();
  hSMVis[0]->Draw("hist same");
  for(unsigned int i=1; i<5; ++i)
    hSMVis[i]->Draw("same hist");
  hMVis[0]->Draw("same");
  leg->Draw();

  CMS_lumi(ca, iPeriod, iPos);
  ca->Update();
  ca->RedrawAxis();
  ca->GetFrame()->Draw();

  ca->SaveAs((std::string("tauEff32-RunC_MVis")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("tauEff32-RunC_MVis")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("tauEff32-RunC_MVis")+std::string(".png")).c_str());

  hMTau[0]->Draw();
  hSMTau[0]->Draw("hist same");
  for(unsigned int i=1; i<5; ++i)
    hSMTau[i]->Draw("same hist");
  hMTau[0]->Draw("same");
  leg->Draw();

  CMS_lumi(ca, iPeriod, iPos);
  ca->Update();
  ca->RedrawAxis();
  ca->GetFrame()->Draw();

  ca->SaveAs((std::string("tauEff32-RunC_MTau")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("tauEff32-RunC_MTau")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("tauEff32-RunC_MTau")+std::string(".png")).c_str());

  hMt[0]->Draw();
  hSMt[0]->Draw("hist same");
  for(unsigned int i=1; i<5; ++i)
    hSMt[i]->Draw("same hist");
  hMt[0]->Draw("same");
  leg->Draw();

  CMS_lumi(ca, iPeriod, iPos);
  ca->Update();
  ca->RedrawAxis();
  ca->GetFrame()->Draw();

  ca->SaveAs((std::string("tauEff32-RunC_Mt")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("tauEff32-RunC_Mt")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("tauEff32-RunC_Mt")+std::string(".png")).c_str());

  hNVtx[0]->Draw();
  hSNVtx[0]->Draw("hist same");
  for(unsigned int i=1; i<5; ++i)
    hSNVtx[i]->Draw("same hist");
  hNVtx[0]->Draw("same");
  leg->Draw();

  CMS_lumi(ca, iPeriod, iPos);
  ca->Update();
  ca->RedrawAxis();
  ca->GetFrame()->Draw();

  ca->SaveAs((std::string("tauEff32-RunC_NVtx")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("tauEff32-RunC_NVtx")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("tauEff32-RunC_NVtx")+std::string(".png")).c_str());

  std::string names[] = {"HLT","L1+HLT","L1"};

  TF1 *myErf[3];
  /*
  myErf[0] = new TF1("myErf0","[0]*(0.5+0.5*TMath::Erf((x-[2])/(sqrt(2)*[1])))",0,100);
  myErf[0]->SetParNames("eff","res","thr");
  myErf[0]->SetParameters(0.8,1.0,20);
  */
  /**/
  myErf[0] = new TF1("myErf0","myCBErf(x,[0],[1],[2],[3],[4])",0,100);
  myErf[0]->SetParNames("m0","sigma","alpha","n","norm");
  myErf[0]->SetParameters(20,1.0,0.5,7,0.8);
  /**/
  myErf[0]->SetLineColor(kBlack);  
  myErf[1]=(TF1*)myErf[0]->Clone("myErf1");
  myErf[1]->SetLineColor(kGreen+3);
  myErf[2]=(TF1*)myErf[0]->Clone("myErf2");
  myErf[2]->SetLineColor(kRed);  

  TH1D *hFrame=(TH1D*)hhPtEff[0][0]->Clone("hFrame");
  hFrame->Reset();
  //zero x-errors by hand
  if(zeroErrX){
    for(unsigned int i=0; i<3; ++i){
      for(unsigned int j=0; j<3; ++j){
	std::cout<<i<<","<<j<<","<<grPtEff[i][j]->GetN()<<std::endl;
	for(int ib=0; ib<grPtEff[i][j]->GetN()-1; ++ib){
	  //std::cout<<i<<","<<j<<","<<ib<<std::endl;
	  grPtEff[i][j]->SetPointEXlow(ib+1,0);
	  grPtEff[i][j]->SetPointEXhigh(ib+1,0);
	}
      }
    }
  }
  for(unsigned int i=0; i<3; ++i){
    /*
    hhPtEff[i][1]->Draw("e");
    hhPtEff[i][0]->Draw("e same");
    hhPtEff[i][2]->Draw("e same");
    */
    hFrame->Draw();
    if(trueDYtt)
      grPtEff[i][1]->Draw("pz same");
    grPtEff[i][0]->Draw("pz same");
    grPtEff[i][2]->Draw("pz same");
    if(fit){
      //hhPtEff[i]->Fit("myErf","IN","",18,120);
      /*
      myErf[0]->SetParameters(0.8,1.0,20);
      hhPtEff[i][0]->Fit("myErf0","","",18,100);
      if(trueDYtt){
        myErf[1]->SetParameters(0.8,1.0,20);
        hhPtEff[i][1]->Fit("myErf1","","",18,100);
      }
      myErf[2]->SetParameters(0.8,1.0,20);
      hhPtEff[i][2]->Fit("myErf2","","",18,100);
      hhPtEff[i][0]->Draw("e, same");
      if(trueDYtt)
        hhPtEff[i][1]->Draw("e, same");
      hhPtEff[i][2]->Draw("e, same");
      */
      myErf[0]->SetParameters(20,1.0,0.5,7,0.8);
      grPtEff[i][0]->Fit("myErf0","M EX0","",18,100);
      //grPtEff[i][0]->Fit("myErf0","M","",18,100);
      if(trueDYtt){
	myErf[1]->SetParameters(20,1.0,0.5,7,0.8);
	grPtEff[i][1]->Fit("myErf1","M EX0","",18,100);
        //grPtEff[i][1]->Fit("myErf1","M","",18,100);
      }
      myErf[2]->SetParameters(20,1.0,0.5,7,0.8);
      grPtEff[i][2]->Fit("myErf2","M EX0","",18,100);
      //grPtEff[i][2]->Fit("myErf2","M","",18,100);
      grPtEff[i][0]->Draw("pz same");
      if(trueDYtt)
	grPtEff[i][1]->Draw("pz same");
      grPtEff[i][2]->Draw("pz same");
    }
    leg2->Clear();
    leg2->AddEntry(grPtEff[i][0],"Data","PLE");
    leg2->AddEntry(grPtEff[i][2],"Simulation","PLE");
    if(trueDYtt)
      leg2->AddEntry(grPtEff[i][2],"Z#rightarrow#tau#tau simulation","PLE");
    leg2->Draw();
    //CMSPrelim(Form(".1%f pb^{-1} at 13TeV",lumi),"", "", 0.10, 0.95);
    //
    /*
    TPaveText* cmsprel  = new TPaveText(0.10+0.07, 0.95-0.053, 0.10+0.45, 0.95-0.033, "NDC");
    cmsprel->SetBorderSize(   0 );
    cmsprel->SetFillStyle(    0 );
    cmsprel->SetTextAlign(   12 );
    //cmsprel->SetTextSize ( 0.03 );
    cmsprel->SetTextSize ( 0.04 );
    cmsprel->SetTextColor(    1 );
    cmsprel->SetTextFont (   62 );//52
    cmsprel->AddText("CMS");
    cmsprel->Draw();
    */
    CMS_lumi(ca, iPeriod, iPos);
    ca->Update();
    ca->RedrawAxis();
    ca->GetFrame()->Draw();

    ca->SaveAs((std::string("tauEff32-RunC_")+names[i]+std::string(".eps")).c_str());  
    ca->SaveAs((std::string("tauEff32-RunC_")+names[i]+std::string(".pdf")).c_str());  
    ca->SaveAs((std::string("tauEff32-RunC_")+names[i]+std::string(".png")).c_str());  
    //////
    //Data vs SS data
    /*
    hhPtEff[i][1]->Draw("e");
    hhPtSSEff[i]->Draw("e same");
    */
    hFrame->Draw();
    grPtEff[i][0]->Draw("pz same");
    grPtSSEff[i]->Draw("pz same");
    if(fit){
      /*
      myErf[0]->SetParameters(0.8,1.0,20);
      hhPtEff[i][0]->Fit("myErf0","","",18,100);
      myErf[2]->SetParameters(0.8,1.0,20);
      hPtSSEff[i]->Fit("myErf2","","",18,100);
      hhPtEff[i][0]->Draw("e, same");
      hPtSSEff[i]->Draw("e, same");
      */
      myErf[0]->SetParameters(20,1.0,0.5,7,0.8);
      grPtEff[i][0]->Fit("myErf0","M EX0","",18,100);
      //grPtEff[i][0]->Fit("myErf0","M","",18,100);
      myErf[2]->SetParameters(20,1.0,0.5,7,0.8);
      grPtSSEff[i]->Fit("myErf2","M EX0","",18,100);
      //grPtSSEff[i]->Fit("myErf1","M","",18,100);

      grPtEff[i][0]->Draw("pz same");
      grPtSSEff[i]->Draw("pz same");
    }
    leg2->Clear();
    leg2->AddEntry(grPtEff[i][0],"OS data","PLE");
    leg2->AddEntry(grPtSSEff[i], "SS data","PLE");
    leg2->Draw();

    CMS_lumi(ca, iPeriod, iPos);
    ca->Update();
    ca->RedrawAxis();
    ca->GetFrame()->Draw();

    ca->SaveAs((std::string("tauEff32-RunC_SS_")+names[i]+std::string(".eps")).c_str());  
    ca->SaveAs((std::string("tauEff32-RunC_SS_")+names[i]+std::string(".pdf")).c_str());  
    ca->SaveAs((std::string("tauEff32-RunC_SS_")+names[i]+std::string(".png")).c_str());  

  }

}


