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

  TChain *t[5];
  //data
  t[0] = new TChain("muMu/muMuTriggerTree");
  t[0]->Add("effTrees/Run2015C_Cert_246908-255031_Eff_MuMu_74X_v1/muMuTrgAna*.root"); //16.1pb-1 (after brilcalc tool)
  t[0]->Add("effTrees/Run2015C_Cert_254833_Eff_MuMu_74X_v1/muMuTrgAna*.root"); //23.2pb-1 (after brilcalc tool)
  //DY
  t[1] = new TChain("muMu/muMuTriggerTree");
  t[1]->Add("effTrees/DYToLL_M50_Eff_MuMu_74X_25ns_v1/muMuTrgAna*.root");
  //W
  t[2] = new TChain("muMu/muMuTriggerTree");
  t[2]->Add("effTrees/WToLNu_Eff_MuMu_74X_25ns_v1/muMuTrgAna*.root");
  //TTbar
  t[3] = new TChain("muMu/muMuTriggerTree");
  t[3]->Add("effTrees/TTbar_Eff_MuMu_74X_25ns_v1/muMuTrgAna*.root");
  //QCD
  t[4] = new TChain("muMu/muMuTriggerTree");
  t[4]->Add("effTrees/QCD-Pt20-MuEnrPt15_Eff_MuMu_74X_25ns_v1/muMuTrgAna*.root");

  double scale[5] = {
        1,     //data
     2008.4*3, //DY x-sec [pb]
    20508.9*3, //W x-sec [pb]
      831.76,  //TTbar x-sec [pb]
720648000*0.00042 //QCD 
  };
  TFile *f[4];
  f[0]=TFile::Open("effTrees/DYToLL_M50_Eff_MuMu_74X_25ns_v1/DYToLL_M50Count.root");
  f[1]=TFile::Open("effTrees/WToLNu_Eff_MuMu_74X_25ns_v1/WToLNuCount.root");
  f[2]=TFile::Open("effTrees/TTbar_Eff_MuMu_74X_25ns_v1/TTbarCount.root");
  f[3]=TFile::Open("effTrees/QCD-Pt20-MuEnrPt15_Eff_MuMu_74X_25ns_v1/QCD-Pt20-MuEnrPt15Count.root");
  for(unsigned int i=0; i<4; ++i){
    TH1D *h = (TH1D*)f[i]->Get("evtCounter/hCounts");
    double w = h->GetBinContent(2);
    scale[i+1] = lumi*scale[i+1]/w;
    f[i]->Close();
  }


  double ptBins[]={0,1,2,3,4,5,6,7,8,9,
		   10,11,12,13,14,15,16,17,18,19,
		   20,21,22,23,24,25,26,27,28,29,
		   30,31,32,33,34,35,36,37,38,39,
		   40,42,44,46,48,
		   50,55,
		   60,65,
		   70,80,90,100
  };
  TH1D *hPt[4][5];
  TH1D *hEta[5];
  TH1D *hM[5];
  TH1D *hNVtx[5];
  for(unsigned int i=0; i<5; ++i){
    for(unsigned int j=0; j<4; ++j){
      hPt[j][i] = new TH1D(Form("hPt_%i%i",j,i)," ;p_{T}^{#mu offline} [GeV]; Events/width [1/GeV]",52,ptBins);//upto 100
      //hPt[j][i] = new TH1D(Form("hPt_%i%i",j,i)," ;p_{T}^{#mu offline} [GeV]; Events",100,0,100);//upto 100        
      hPt[j][i]->Sumw2();
      hPt[j][i]->SetMarkerStyle(20);
      //hPt[j][i]->SetMarkerSize(0.7);
      hPt[j][i]->SetStats(0);
      //hPt[j][i]->SetMinimum(0.1);
    }
    hEta[i] = new TH1D(Form("hEta_%i",i)," ;#eta^{#mu offline}; Events",50,-2.5,2.5);
    hEta[i]->Sumw2();
    hEta[i]->SetMarkerStyle(20);
    //hEta[i]->SetMarkerSize(0.7);
    hEta[i]->SetStats(0);
    //hEta[i]->SetMinimum(0.1);

    hM[i] = new TH1D(Form("hM_%i",i)," ;M(#mu,#mu) [GeV]; Events",100,50,150);
    hM[i]->Sumw2();
    hM[i]->SetMarkerStyle(20);
    //hM[i]->SetMarkerSize(0.7);
    hM[i]->SetStats(0);
    hM[i]->SetMinimum(0.01);

    hNVtx[i] = new TH1D(Form("hNVtx_%i",i)," ;No. of vertices; Events",50,0,50);
    hNVtx[i]->Sumw2();
    hNVtx[i]->SetMarkerStyle(20);
    //hNVtx[i]->SetMarkerSize(0.7);
    hNVtx[i]->SetStats(0);
    //hNVtx[i]->SetMinimum(0.1);
  }  

  //std::string muTag = "(HLT_IsoMu20_eta2p1_v>0 && tag_hltL3crIsoL1sMu16Eta2p1L1f0L2f10QL3f20QL3trkIsoFiltered0p09>0)";
  std::string muTag = "(HLT_IsoMu24_eta2p1_v>0 && tag_hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09>0)";
  //std::string off   = "tagPt>21 && abs(tagEta)<2.1 && tagIso<0.1 && probePt>15 && abs(probeEta)<2.1 && probeIso<0.1 && tagCharge*probeCharge<0";
  std::string off   = "tagPt>25 && abs(tagEta)<2.1 && tagIso<0.1 && probePt>15 && abs(probeEta)<2.1 && probeIso<0.1 && tagCharge*probeCharge<0";
  std::string window = "abs(Mass-91)<5";
  std::string sel = muTag + " && " + off + " && " + window;
  std::string selNoM = muTag + " && " + off;
  std::string vtxReWeight="*vtxWeight(nVtx)";

  //filters to check
  std::vector<std::string> probe;
  std::vector<std::string> names;
  probe.push_back("probe_hltL3crIsoL1sDoubleMu125L1f16erL2f10QL3f17QL3Dz0p2L3crIsoRhoFiltered0p15IterTrk02>0");//DoubleIsoMu17_eta2p1
  names.push_back("DoubleIsoMu17_eta2p1");
  probe.push_back("probe_hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09>0");//IsoMu17_eta2p1 - tau leg of signle seeded mu+tau
  names.push_back("IsoMu17_eta2p1");
  // regions
  std::string region_name[] = {"inclusive","barrel","overlap","endcap"};

  if(!doVtxReWeight)
    vtxReWeight="*vtxWeight(nVtx,false)";

  std::string cut, cutNoM;
  for(unsigned int i=0; i<5; ++i){
    if(i>4){
      std::cout<<"What?!"<<std::endl;
      continue;
    }
    cut = "weight*("+sel+")";
    cutNoM = "weight*("+selNoM+")";
    if(i!=0){
      cut = cut+vtxReWeight;
      cutNoM = cutNoM+vtxReWeight;
    }
    std::cout<<i<<std::endl;
    std::cout<<cut<<std::endl;

    for(unsigned int j=0; j<4; ++j){
      std::string cutEta=sel;
      if(j==1)
	cutEta = cutEta + " && abs(probeEta)<0.9";//barrel
      else if(j==2)
	cutEta = cutEta + " && abs(probeEta)<1.2 && abs(probeEta)>0.9";//overlap
      else if(j==3)
	cutEta = cutEta + " && abs(probeEta)<2.1 && abs(probeEta)>1.2";//endcap
      cutEta = "weight*("+cutEta+")";
      if(i!=0)
	cutEta = cutEta+vtxReWeight;
      std::cout<<"\t"<<j<<". cut with eta: "<<cutEta<<std::endl;
      t[i]->Project(Form("hPt_%i%i",j,i),"probePt",cutEta.c_str(),"");
      hPt[j][i]->Scale(scale[i],"width");
    }
    std::cout<<"Expected: "<<hPt[0][i]->Integral(0,hPt[0][i]->GetNbinsX()+1)
	     <<" for L="<<lumi<<"pb-1"<<std::endl;
    std::cout<<hPt[0][i]->Integral(0,hPt[0][i]->GetNbinsX()+1,"")<<std::endl;

    t[i]->Project(Form("hM_%i",i),"Mass",cutNoM.c_str(),"");
    hM[i]->Scale(scale[i],"");

    t[i]->Project(Form("hEta_%i",i),"probeEta",cut.c_str(),"");
    hEta[i]->Scale(scale[i],"");
    std::cout<<hEta[i]->Integral(0,hEta[i]->GetNbinsX()+1)<<std::endl;

    t[i]->Project(Form("hNVtx_%i",i),"nVtx",cut.c_str(),"");
    hNVtx[i]->Scale(scale[i]);

    std::cout<<"end of loop: "<<i<<std::endl;
  }

  //stacked
  //              DY       ,W+jets,ttbar  ,QCD
  int colMap[] = {kOrange-2,kRed+2,kBlue-1,kMagenta-10};
  TH1D *hSPt[4][4];
  TH1D *hSEta[4];
  TH1D *hSM[4];
  TH1D *hSNVtx[4];
  for(unsigned int i=0; i<4; ++i){
    for(unsigned int j=0; j<4; ++j){
      hSPt[j][i]=(TH1D*)hPt[0][1]->Clone(Form("hSPt_%i%i",j,i));
      hSPt[j][i]->Reset();
      hSPt[j][i]->SetFillColor(colMap[i]);
    }
    hSEta[i]=(TH1D*)hEta[1]->Clone(Form("hSEta_%i",i));
    hSEta[i]->Reset();
    hSEta[i]->SetFillColor(colMap[i]);
    hSM[i]=(TH1D*)hM[1]->Clone(Form("hSM_%i",i));
    hSM[i]->Reset();
    hSM[i]->SetFillColor(colMap[i]);
    hSNVtx[i]=(TH1D*)hNVtx[1]->Clone(Form("hSNVtx_%i",i));
    hSNVtx[i]->Reset();
    hSNVtx[i]->SetFillColor(colMap[i]);
    for(unsigned int j=i+1; j<5; ++j){
      for(unsigned int k=0; k<4; ++k){
	hSPt[k][i]->Add(hPt[k][j]);
      }
      hSEta[i]->Add(hEta[j]);
      hSM[i]->Add(hM[j]);
      hSNVtx[i]->Add(hNVtx[j]);
    }
  }
  std::cout<<"Data: "<<hEta[0]->Integral(0,hEta[0]->GetNbinsX()+1)
	   <<", MC: "<<hSEta[0]->Integral(0,hSEta[0]->GetNbinsX()+1)<<std::endl;
  std::cout<<"Data (NoM): "<<hM[0]->Integral(0,hM[0]->GetNbinsX()+1)
	   <<", MC (NoM): "<<hSM[0]->Integral(0,hSM[0]->GetNbinsX()+1)<<std::endl;
  
  ///Effs
  // 2(data&MC) x 4 regions in eta x 2 probes  
  TH1D *hhPt[2][4][probe.size()];
  TH1D *hhPtEff[2][4][probe.size()];
  TGraphAsymmErrors *grPtEff[2][4][probe.size()];
  //TH2D *hFrame = new TH2D("hFrame",";p_{T}^{#tau offline} [GeV]; Efficiency",2,0,100,2,0.01,1.10);
  //hFrame->SetStats(0);

  for(unsigned int iProb=0; iProb<probe.size(); ++iProb){
    for(unsigned int i=0; i<2; ++i){
      for(unsigned int j=0; j<4; ++j){
	hhPt[i][j][iProb]=(TH1D*)hPt[1][0]->Clone(Form("hhPt_%i%i%i",i,j,iProb)); 
	//hhPt[i][j][iProb]->Sumw2();
	hhPt[i][j][iProb]->Reset();     

	hhPtEff[i][j][iProb]=(TH1D*)hPt[1][0]->Clone(Form("hhPtEff_%i%i%i",i,j,iProb)); 
	//hhPtEff[i][j][iProb]->Sumw2();
	hhPtEff[i][j][iProb]->SetTitle(";p_{T}^{#mu offline} [GeV]; Efficiency");
	hhPtEff[i][j][iProb]->Reset();       
      }
    }
  }
  //Data and all MC's
  TH1D *hhPtTmp=(TH1D*)hPt[0][1]->Clone("hhPtTmp"); 
  //hhPtTmp->Sumw2();
  hhPtTmp->Reset();
  for(unsigned int iProb=0; iProb<probe.size(); ++iProb){
    for(unsigned int i=0; i<5; ++i){
      for(unsigned int j=0; j<4; ++j){
	std::string cutEta=sel + " && " + probe[iProb];
	if(j==1)
	  cutEta = cutEta + " && abs(probeEta)<0.9" + " && " + probe[iProb];//barrel
	else if(j==2)
	  cutEta = cutEta + " && abs(probeEta)<1.2 && abs(probeEta)>0.9" + "&& " + probe[iProb];//overlap
	else if(j==3)
	  cutEta = cutEta + " && abs(probeEta)<2.1 && abs(probeEta)>1.2" + " && " + probe[iProb];//endcap
	cutEta = "weight*("+cutEta+")";
	if(i!=0){
	  cutEta = cutEta+vtxReWeight;
	  std::cout<<"\t"<<j<<". cut with eta and probed trigger : "<<cutEta<<std::endl;
	  hhPtTmp->Reset();
	  t[i]->Project("hhPtTmp","probePt",cutEta.c_str(),"");
	  hhPtTmp->Scale(scale[i],"width");
	  hhPt[1][j][iProb]->Add(hhPtTmp);
	}
	else{
	  std::cout<<"\t"<<j<<". cut with eta and probed trigger : "<<cutEta<<std::endl;
	  t[i]->Project(Form("hhPt_%i%i%i",0,j,iProb),"probePt",cutEta.c_str(),"");
	  hhPt[0][j][iProb]->Scale(scale[i],"width");
	}
      }
    }
  }

  for(unsigned int j=0; j<4; ++j){
    //Data
    //reset overflows
    hPt[j][0]->SetBinContent(hPt[j][0]->GetNbinsX()+1,0);
    hPt[j][0]->SetBinError(hPt[j][0]->GetNbinsX()+1,0);
    for(unsigned int iProb=0; iProb<probe.size(); ++iProb){
      //reset overflows
      hhPt[0][j][iProb]->SetBinContent(hhPt[0][j][iProb]->GetNbinsX()+1,0);
      hhPt[0][j][iProb]->SetBinError(hhPt[0][j][iProb]->GetNbinsX()+1,0);

      //efficiency
      hhPtEff[0][j][iProb]->Divide(hhPt[0][j][iProb],hPt[j][0],1,1,"b"); 
      grPtEff[0][j][iProb] = new TGraphAsymmErrors(hhPt[0][j][iProb],hPt[j][0],"cl=0.683 b(1,1) mode");
      
      //style
      hhPtEff[0][j][iProb]->SetMinimum(0.01); 
      hhPtEff[0][j][iProb]->SetMaximum(1.10);
      hhPtEff[0][j][iProb]->SetLineColor(kBlack);
      hhPtEff[0][j][iProb]->SetMarkerColor(kBlack);
      hhPtEff[0][j][iProb]->SetMarkerColor(20);
      grPtEff[0][j][iProb]->SetLineColor(kBlack);
      grPtEff[0][j][iProb]->SetMarkerColor(kBlack);
      grPtEff[0][j][iProb]->SetMarkerStyle(hhPtEff[0][j][iProb]->GetMarkerStyle());
    }
    //all MC
    TH1D *hPtMC = (TH1D*)hSPt[j][0]->Clone(Form("hPtMC_%i",j));
    //reset overflows
    hPtMC->SetBinContent(hPtMC->GetNbinsX()+1,0);
    hPtMC->SetBinError(hPtMC->GetNbinsX()+1,0);
    
    for(unsigned int iProb=0; iProb<probe.size(); ++iProb){
      //reset overflows
      hhPt[1][j][iProb]->SetBinContent(hhPt[1][j][iProb]->GetNbinsX()+1,0);
      hhPt[1][j][iProb]->SetBinError(hhPt[1][j][iProb]->GetNbinsX()+1,0);
      
      //add hoc safety checks
      for(int ij=0; ij<hPtMC->GetNbinsX()+2; ++ij){
	double ratio = hPtMC->GetBinContent(ij)!=0. ? hhPt[1][j][iProb]->GetBinContent(ij)/hPtMC->GetBinContent(ij) : -100.*hhPt[1][j][iProb]->GetBinContent(ij);
	//std::cout<<"hhPt[1]["<<j<<"]["<<iProb<<"]/hPtMC("<<ij<<") "<<hhPt[1][j][iProb]->GetBinContent(ij)<<"/"
	//	       <<hPtMC->GetBinContent(ij)<<" = "<<ratio<<std::endl;
	if(ratio>1) {	
	  std::cout<<" >1 !"<<std::endl;
	  std::cout<<"hhPt[1]["<<j<<"]["<<iProb<<"]/hPtMC("<<ij<<") "<<hhPt[1][j][iProb]->GetBinContent(ij)<<"/"
		   <<hPtMC->GetBinContent(ij)<<" = "<<ratio<<std::endl;
	  if(ratio-1<0.005)//add hoc 'renormalization' for safety
	    hhPt[1][j][iProb]->SetBinContent(ij,hPtMC->GetBinContent(ij));
	}
      }

      //efficiency
      hhPtEff[1][j][iProb]->Divide(hhPt[1][j][iProb],hPtMC,1,1,"b"); 
      grPtEff[1][j][iProb] = new TGraphAsymmErrors(hhPt[1][j][iProb],hPtMC,"cl=0.683 b(1,1) mode");
      
      //style
      hhPtEff[1][j][iProb]->SetMinimum(0.01); 
      hhPtEff[1][j][iProb]->SetMaximum(1.10);
      hhPtEff[1][j][iProb]->SetLineColor(kRed);
      hhPtEff[1][j][iProb]->SetMarkerColor(kRed);
      hhPtEff[1][j][iProb]->SetMarkerStyle(21);
      grPtEff[1][j][iProb]->SetLineColor(kRed);
      grPtEff[1][j][iProb]->SetMarkerColor(kRed);
      grPtEff[1][j][iProb]->SetMarkerStyle(hhPtEff[1][j][iProb]->GetMarkerStyle());
      //grPtEff[1][j][iProb]->SetMarkerStyle(21);
    }
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

  // Pt 
  for(unsigned int j=0; j<4; ++j){
    if( hPt[j][0]->GetBinContent(hPt[j][0]->GetMaximumBin()) > hSPt[j][0]->GetBinContent(hSPt[j][0]->GetMaximumBin()) )
      hPt[j][0]->Draw();
    else
      hSPt[j][0]->Draw("hist");
    hSPt[j][0]->Draw("hist same");
    for(unsigned int i=1; i<4; ++i)
      hSPt[j][i]->Draw("same hist");
    hPt[j][0]->Draw("same");
    if(j==0){//do it once
      leg->AddEntry(hPt[j][0], " Data","PLE");
      leg->AddEntry(hSPt[j][0]," Z#rightarrow#mu#mu","F");
      leg->AddEntry(hSPt[j][1]," W+jets","F");
      leg->AddEntry(hSPt[j][2]," t#bar{t}","F");
      leg->AddEntry(hSPt[j][3]," QCD","F");
    }
    leg->Draw();

    CMS_lumi(ca, iPeriod, iPos);
    ca->Update();
    ca->RedrawAxis();
    ca->GetFrame()->Draw();

    ca->SaveAs((std::string("muEff-RunC_Pt_")+region_name[j]+std::string(".eps")).c_str());  
    ca->SaveAs((std::string("muEff-RunC_Pt_")+region_name[j]+std::string(".pdf")).c_str());
    ca->SaveAs((std::string("muEff-RunC_Pt_")+region_name[j]+std::string(".png")).c_str());
  }

  if( hEta[0]->GetBinContent(hEta[0]->GetMaximumBin()) > hSEta[0]->GetBinContent(hSEta[0]->GetMaximumBin()) )
    hEta[0]->Draw();
  else
    hSEta[0]->Draw("hist");
  hSEta[0]->Draw("hist same");
  for(unsigned int i=1; i<4; ++i)
    hSEta[i]->Draw("same hist");
  hEta[0]->Draw("same");
  leg->Draw();

  CMS_lumi(ca, iPeriod, iPos);
  ca->Update();
  ca->RedrawAxis();
  ca->GetFrame()->Draw();

  ca->SaveAs((std::string("muEff-RunC_Eta")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("muEff-RunC_Eta")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("muEff-RunC_Eta")+std::string(".png")).c_str());

  if( hM[0]->GetBinContent(hM[0]->GetMaximumBin()) > hSM[0]->GetBinContent(hSM[0]->GetMaximumBin()) )
    hM[0]->Draw();
  else
    hSM[0]->Draw("hist");
  hSM[0]->Draw("hist same");
  for(unsigned int i=1; i<4; ++i)
    hSM[i]->Draw("same hist");
  hM[0]->Draw("same");
  leg->Draw();

  CMS_lumi(ca, iPeriod, iPos);
  ca->Update();
  ca->RedrawAxis();
  ca->GetFrame()->Draw();

  ca->SaveAs((std::string("muEff-RunC_Mass")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("muEff-RunC_Mass")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("muEff-RunC_Mass")+std::string(".png")).c_str());
  
  //mass in log
  ca->SetLogy(1);
  if( hM[0]->GetBinContent(hM[0]->GetMaximumBin()) > hSM[0]->GetBinContent(hSM[0]->GetMaximumBin()) )
    hM[0]->Draw();
  else
    hSM[0]->Draw("hist");
  hSM[0]->Draw("hist same");
  for(unsigned int i=1; i<4; ++i)
    hSM[i]->Draw("same hist");
  hM[0]->Draw("same");
  leg->Draw();

  CMS_lumi(ca, iPeriod, iPos);
  ca->Update();
  ca->RedrawAxis();
  ca->GetFrame()->Draw();

  ca->SaveAs((std::string("muEff-RunC_Mass_Log")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("muEff-RunC_Mass_Log")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("muEff-RunC_Mass_Log")+std::string(".png")).c_str());
  ca->SetLogy(0);
  

  if( hNVtx[0]->GetBinContent(hNVtx[0]->GetMaximumBin()) > hSNVtx[0]->GetBinContent(hSNVtx[0]->GetMaximumBin()) )
    hNVtx[0]->Draw();
  else
    hSNVtx[0]->Draw("hist");
  hSNVtx[0]->Draw("hist same");
  for(unsigned int i=1; i<4; ++i)
    hSNVtx[i]->Draw("same hist");
  hNVtx[0]->Draw("same");
  leg->Draw();

  CMS_lumi(ca, iPeriod, iPos);
  ca->Update();
  ca->RedrawAxis();
  ca->GetFrame()->Draw();

  ca->SaveAs((std::string("muEff-RunC_NVtx")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("muEff-RunC_NVtx")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("muEff-RunC_NVtx")+std::string(".png")).c_str());


  //Effs
  TF1 *myErf[2];
  /**
  myErf[0] = new TF1("myErf0","[0]*(0.5+0.5*TMath::Erf((x-[2])/(sqrt(2)*[1])))",0,100);
  myErf[0]->SetParNames("eff","res","thr");
  myErf[0]->SetParameters(0.8,1.0,20);
  **/
  /**/
  myErf[0] = new TF1("myErf0","myCBErf(x,[0],[1],[2],[3],[4])",0,100);
  myErf[0]->SetParNames("m0","sigma","alpha","n","norm");
  myErf[0]->SetParameters(18,5,5,20,0.95);
  //myErf[0]->SetParameters(17,0.55,0.81,1.55,0.9);
  /**/
  myErf[0]->SetLineColor(kBlack);  
  myErf[1]=(TF1*)myErf[0]->Clone("myErf1");
  myErf[1]->SetLineColor(kRed);  

  TH1D *hFrame=(TH1D*)hhPtEff[0][0][0]->Clone("hFrame");
  hFrame->Reset();
  //zero x-errors by hand
  if(zeroErrX){
    for(unsigned int i=0; i<2; ++i){
      for(unsigned int j=0; j<4; ++j){
	for(unsigned int iProb=0; iProb<probe.size(); ++iProb){
	  std::cout<<i<<","<<j<<","<<iProb<<","<<grPtEff[i][j][iProb]->GetN()<<std::endl;
	  for(int ib=0; ib<grPtEff[i][j][iProb]->GetN()-1; ++ib){
	    //std::cout<<i<<","<<j<<","<<ib<<std::endl;
	    grPtEff[i][j][iProb]->SetPointEXlow(ib+1,0);
	    grPtEff[i][j][iProb]->SetPointEXhigh(ib+1,0);
	  }
	}
      }
    }
  }
  for(unsigned int iProb=0; iProb<probe.size(); ++iProb){
    for(unsigned int j=0; j<4; ++j){
      hFrame->Draw();
      grPtEff[0][j][iProb]->Draw("pz same");
      grPtEff[1][j][iProb]->Draw("pz same");
      if(fit){
	/*
	  myErf[0]->SetParameters(0.8,1.0,20);
	  hhPtEff[0][j][iProb]->Fit("myErf0","","",18,100);
	  myErf[1]->SetParameters(0.8,1.0,20);
	  hhPtEff[1][j][iProb]->Fit("myErf1","","",18,100);
	  hhPtEff[0][j][iProb]->Draw("e, same");
	  hhPtEff[1][j][iProb]->Draw("e, same");
	*/
	//myErf[0]->SetParameters(17,0.55,0.81,1.55,0.9);
	myErf[0]->SetParameters(18,5,5,20,0.95);
	grPtEff[0][j][iProb]->Fit("myErf0","M EX0","",18,100);
	//grPtEff[0][j][iProb]->Fit("myErf0","M","",18,100);
	//myErf[1]->SetParameters(17,0.55,0.81,1.55,0.9);
	myErf[1]->SetParameters(18,5,5,20,0.95);
	grPtEff[1][j][iProb]->Fit("myErf1","M EX0","",18,100);
	//grPtEff[1][j][iProb]->Fit("myErf1","M","",18,100);
	grPtEff[0][j][iProb]->Draw("pz same");
	grPtEff[1][j][iProb]->Draw("pz same");
      }
      leg2->Clear();
      leg2->AddEntry(grPtEff[0][j][iProb],"Data","PLE");
      leg2->AddEntry(grPtEff[1][j][iProb],"Simulation","PLE");
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

      ca->SaveAs((std::string("muEff-RunC_")+names[iProb]+std::string("_")+region_name[j]+std::string(".eps")).c_str());  
      ca->SaveAs((std::string("muEff-RunC_")+names[iProb]+std::string("_")+region_name[j]+std::string(".pdf")).c_str());  
      ca->SaveAs((std::string("muEff-RunC_")+names[iProb]+std::string("_")+region_name[j]+std::string(".png")).c_str());  
    }
  }

  return;
}


