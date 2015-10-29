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
#include "TFitResult.h"
#include "TFitResultPtr.h"

#include "tdrstyle.C" //sic! *.C included
#include "CMS_lumi.C" //sic! *.C included

//Implementation of CBErf from HTT TWiki by J.Swanson, R.Lane
double myCBErf(double m, double m0, double sigma, double alpha,
		  double n, double norm){
  //Useful constants
  const double sqrtPiOver2 = sqrt(TMath::PiOver2());
  const double sqrt2 = sqrt(2.);

  
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
TFile *fVtxW=TFile::Open("vtxNorm-RunD.root");
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
  for(unsigned int i=0; i<4; ++i){
    TH1D *h = (TH1D*)f[i]->Get("evtCounter/hCounts");
    double w = h->GetBinContent(2);
    scale[i+1] = lumi*scale[i+1]/w;
    f[i]->Close();
  }

  // regions
  std::vector<std::string> region_name;
  std::vector<std::string> region_def;
  region_name.push_back("inclusive");
  region_def.push_back("1");
  region_name.push_back("barrel");
  region_def.push_back("abs(probeEta)<0.8 && abs(probeEta)>=0.0");
  region_name.push_back("overlap");
  region_def.push_back("abs(probeEta)<1.2 && abs(probeEta)>=0.8");
  region_name.push_back("endcap");
  region_def.push_back("abs(probeEta)<2.1 && abs(probeEta)>=1.2");

  double ptBins[]={0,1,2,3,4,5,6,7,8,9,
		   10,11,12,13,14,15,16,17,18,19,
		   20,21,22,23,24,25,26,27,28,29,
		   30,31,32,33,34,35,36,37,38,39,
		   40,42,44,46,48,
		   50,55,
		   60,65,
		   70,80,90,100
  };
  double etaBins[]={-2.4,-2.1,-1.6,-1.2,-0.8,-0.3,
		    -0.2, 0.2,
		    0.3, 0.8, 1.2, 1.6, 2.1, 2.4};
  TH1D *hPt[region_def.size()][5];
  TH1D *hEta[5];
  TH1D *hEta20[5];
  TH1D *hM[5];
  TH1D *hNVtx[5];
  TH1D *hNVtx20[5];
  for(unsigned int i=0; i<5; ++i){
    for(unsigned int j=0; j<region_def.size(); ++j){
      hPt[j][i] = new TH1D(Form("hPt_%i%i",j,i)," ;p_{T}^{#mu offline} [GeV]; Events/width [1/GeV]",52,ptBins);//upto 100
      //hPt[j][i] = new TH1D(Form("hPt_%i%i",j,i)," ;p_{T}^{#mu offline} [GeV]; Events",100,0,100);//upto 100        
      hPt[j][i]->Sumw2();
      hPt[j][i]->SetMarkerStyle(20);
      //hPt[j][i]->SetMarkerSize(0.7);
      hPt[j][i]->SetStats(0);
      //hPt[j][i]->SetMinimum(0.1);
    }
    hEta[i] = new TH1D(Form("hEta_%i",i)," ;#eta^{#mu offline}; Events",48,-2.4,2.4);
    hEta[i]->Sumw2();
    hEta[i]->SetMarkerStyle(20);
    //hEta[i]->SetMarkerSize(0.7);
    hEta[i]->SetStats(0);
    //hEta[i]->SetMinimum(0.1);

    hEta20[i] = new TH1D(Form("hEta20_%i",i)," ;#eta^{#mu offline}; Events",13,etaBins);
    hEta20[i]->Sumw2();
    hEta20[i]->SetMarkerStyle(20);
    //hEta20[i]->SetMarkerSize(0.7);
    hEta20[i]->SetStats(0);
    //hEta20[i]->SetMinimum(0.1);

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

    hNVtx20[i] = new TH1D(Form("hNVtx20_%i",i)," ;No. of vertices; Events",50,0,50);
    hNVtx20[i]->Sumw2();
    hNVtx20[i]->SetMarkerStyle(20);
    //hNVtx20[i]->SetMarkerSize(0.7);
    hNVtx20[i]->SetStats(0);
    //hNVtx20[i]->SetMinimum(0.1);
  }  

  std::string muTag_MC   = "(HLT_IsoMu17_eta2p1_v>0 && tag_hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09>18)";
  std::string muTag_data = "(HLT_IsoMu18_v>0 && tag_hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09>18)";
  std::string off   = "tagPt>20 && abs(tagEta)<2.1 && tagIso<0.1 && probePt>10 && abs(probeEta)<2.1 && probeIso<0.1 && tagCharge*probeCharge<0";
  std::string window = "abs(Mass-91)<5";
  std::string vtxReWeight="*vtxWeight(nVtx)";

  //filters to check
  std::vector<std::string> probeMC, probeData;
  std::vector<std::string> names;
  probeMC.push_back("probe_hltL3crIsoL1sDoubleMu125L1f16erL2f10QL3f17QL3Dz0p2L3crIsoRhoFiltered0p15IterTrk02>0");//DoubleIsoMu17_eta2p1
  probeData.push_back("probe_hltL3crIsoL1sDoubleMu125L1f16erL2f10QL3f17QL3Dz0p2L3crIsoRhoFiltered0p15IterTrk02>0");//DoubleIsoMu17_eta2p1
  names.push_back("DoubleIsoMu17_eta2p1");
  probeMC.push_back("probe_hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09>0");//IsoMu17_eta2p1 - tau leg of signle seeded mu+tau
  probeData.push_back("probe_hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09>0");//IsoMu17_eta2p1 - tau leg of signle seeded mu+tau
  names.push_back("IsoMu17_eta2p1");
  probeMC.push_back("probe_hltL3crIsoL1sSingleMu16erL1f0L2f10QL3f17QL3trkIsoFiltered0p09>18");//IsoMu17_eta2p1 - tau leg of signle seeded mu+tau titghtened to 18GeV
  probeData.push_back("probe_hltL3crIsoL1sMu16L1f0L2f10QL3f18QL3trkIsoFiltered0p09>0");//IsoMu18
  names.push_back("IsoMu18");

  if(!doVtxReWeight)
    vtxReWeight="*vtxWeight(nVtx,false)";

  std::string cut, cutNoM, cutPt;
  for(unsigned int i=0; i<5; ++i){
    if(i>4){
      std::cout<<"What?!"<<std::endl;
      continue;
    }
    if(i!=0){
      cut    = "weight*("+muTag_MC+" && "+off+" && "+window+")"+vtxReWeight;
      cutPt  = "weight*("+muTag_MC+" && "+off+" && "+window+" && probePt>20)"+vtxReWeight;
      cutNoM = "weight*("+muTag_MC+" && "+off+              ")"+vtxReWeight;
    }
    else{
      cut    = "weight*("+muTag_data+" && "+off+" && "+window+")";
      cutPt  = "weight*("+muTag_data+" && "+off+" && "+window+" && probePt>20)";
      cutNoM = "weight*("+muTag_data+" && "+off+              ")";
    }
    std::cout<<i<<std::endl;
    std::cout<<cut<<std::endl;

    for(unsigned int j=0; j<region_def.size(); ++j){
      std::string cutEta = muTag_MC+" && "+off+" && "+window;
      if(i==0){
	cutEta = muTag_data+" && "+off+" && "+window;
      }
      cutEta = "weight*("+cutEta + " && " + region_def[j]+")";
      if(i!=0)
	cutEta = cutEta+vtxReWeight;
      std::cout<<"\t"<<j<<". cut with eta: "<<cutEta<<std::endl;
      t[i]->Project(Form("hPt_%i%i",j,i),"probePt",cutEta.c_str(),"");
      hPt[j][i]->Scale(scale[i],"width");
    }
    std::cout<<"Expected: "<<hPt[0][i]->Integral(0,hPt[0][i]->GetNbinsX()+1,"width")
	     <<" for L="<<lumi<<"pb-1"<<std::endl;
    std::cout<<hPt[0][i]->Integral(0,hPt[0][i]->GetNbinsX()+1,"width")<<std::endl;

    t[i]->Project(Form("hM_%i",i),"Mass",cutNoM.c_str(),"");
    hM[i]->Scale(scale[i],"");

    t[i]->Project(Form("hEta_%i",i),"probeEta",cut.c_str(),"");
    hEta[i]->Scale(scale[i],"");
    std::cout<<hEta[i]->Integral(0,hEta[i]->GetNbinsX()+1)<<std::endl;

    t[i]->Project(Form("hEta20_%i",i),"probeEta",cutPt.c_str(),"");
    hEta20[i]->Scale(scale[i],"width");

    t[i]->Project(Form("hNVtx_%i",i),"nVtx",cut.c_str(),"");
    hNVtx[i]->Scale(scale[i]);

    t[i]->Project(Form("hNVtx20_%i",i),"nVtx",cutPt.c_str(),"");
    hNVtx20[i]->Scale(scale[i]);

    std::cout<<"end of loop: "<<i<<std::endl;
  }

  //stacked
  //              DY       ,W+jets,ttbar  ,QCD
  int colMap[] = {kOrange-2,kRed+2,kBlue-1,kMagenta-10};
  TH1D *hSPt[region_def.size()][4];
  TH1D *hSEta[4];
  TH1D *hSEta20[4];
  TH1D *hSM[4];
  TH1D *hSNVtx[4];
  TH1D *hSNVtx20[4];
  for(unsigned int i=0; i<4; ++i){
    for(unsigned int j=0; j<region_def.size(); ++j){
      hSPt[j][i]=(TH1D*)hPt[0][1]->Clone(Form("hSPt_%i%i",j,i));
      hSPt[j][i]->Reset();
      hSPt[j][i]->SetFillColor(colMap[i]);
    }
    hSEta[i]=(TH1D*)hEta[1]->Clone(Form("hSEta_%i",i));
    hSEta[i]->Reset();
    hSEta[i]->SetFillColor(colMap[i]);
    hSEta20[i]=(TH1D*)hEta20[1]->Clone(Form("hSEta20_%i",i));
    hSEta20[i]->Reset();
    hSEta20[i]->SetFillColor(colMap[i]);
    hSM[i]=(TH1D*)hM[1]->Clone(Form("hSM_%i",i));
    hSM[i]->Reset();
    hSM[i]->SetFillColor(colMap[i]);
    hSNVtx[i]=(TH1D*)hNVtx[1]->Clone(Form("hSNVtx_%i",i));
    hSNVtx[i]->Reset();
    hSNVtx[i]->SetFillColor(colMap[i]);
    hSNVtx20[i]=(TH1D*)hNVtx20[1]->Clone(Form("hSNVtx20_%i",i));
    hSNVtx20[i]->Reset();
    hSNVtx20[i]->SetFillColor(colMap[i]);
    for(unsigned int j=i+1; j<5; ++j){
      for(unsigned int k=0; k<region_def.size(); ++k){
	hSPt[k][i]->Add(hPt[k][j]);
      }
      hSEta[i]->Add(hEta[j]);
      hSEta20[i]->Add(hEta20[j]);
      hSM[i]->Add(hM[j]);
      hSNVtx[i]->Add(hNVtx[j]);
      hSNVtx20[i]->Add(hNVtx20[j]);
    }
  }
  std::cout<<"Data: "<<hEta[0]->Integral(0,hEta[0]->GetNbinsX()+1)
	   <<", MC: "<<hSEta[0]->Integral(0,hSEta[0]->GetNbinsX()+1)<<std::endl;
  std::cout<<"Data (NoM): "<<hM[0]->Integral(0,hM[0]->GetNbinsX()+1)
	   <<", MC (NoM): "<<hSM[0]->Integral(0,hSM[0]->GetNbinsX()+1)<<std::endl;
  
  ///Effs
  // 2(data&MC) x 4 regions in eta x 2 probes  
  TH1D *hhPt[2][4][probeMC.size()];
  TH1D *hhEta[2][probeMC.size()];
  TH1D *hhNVtx[2][probeMC.size()];
  TH1D *hhPtEff[2][4][probeMC.size()];
  TGraphAsymmErrors *grPtEff[2][4][probeMC.size()];
  TGraphAsymmErrors *grEtaEff[2][probeMC.size()];
  TGraphAsymmErrors *grNVtxEff[2][probeMC.size()];
  //TH2D *hFrame = new TH2D("hFrame",";p_{T}^{#tau offline} [GeV]; Efficiency",2,0,100,2,0.01,1.10);
  //hFrame->SetStats(0);

  for(unsigned int iProb=0; iProb<probeMC.size(); ++iProb){
    for(unsigned int i=0; i<2; ++i){
      for(unsigned int j=0; j<region_def.size(); ++j){
	hhPt[i][j][iProb]=(TH1D*)hPt[1][0]->Clone(Form("hhPt_%i%i%i",i,j,iProb)); 
	//hhPt[i][j][iProb]->Sumw2();
	hhPt[i][j][iProb]->Reset();     

	hhPtEff[i][j][iProb]=(TH1D*)hPt[1][0]->Clone(Form("hhPtEff_%i%i%i",i,j,iProb)); 
	//hhPtEff[i][j][iProb]->Sumw2();
	hhPtEff[i][j][iProb]->SetTitle(";p_{T}^{#mu offline} [GeV]; Efficiency");
	hhPtEff[i][j][iProb]->Reset();       
      }
      hhEta[i][iProb]=(TH1D*)hEta20[0]->Clone(Form("hhEta_%i%i",i,iProb)); 
      //hhEta[i][iProb]->Sumw2();
      hhEta[i][iProb]->Reset();     
      hhNVtx[i][iProb]=(TH1D*)hNVtx20[0]->Clone(Form("hhNVtx_%i%i",i,iProb)); 
      //hhNVtx[i][iProb]->Sumw2();
      hhNVtx[i][iProb]->Reset();     
    }
  }
  //Data and all MC's
  TH1D *hhPtTmp=(TH1D*)hPt[0][1]->Clone("hhPtTmp"); 
  //hhPtTmp->Sumw2();
  hhPtTmp->Reset();
  TH1D *hhEtaTmp=(TH1D*)hEta20[0]->Clone("hhEtaTmp"); 
  //hhEtaTmp->Sumw2();
  hhEtaTmp->Reset();
  TH1D *hhNVtxTmp=(TH1D*)hNVtx20[0]->Clone("hhNVtxTmp"); 
  //hhNVtxTmp->Sumw2();
  hhNVtxTmp->Reset();
  for(unsigned int iProb=0; iProb<probeMC.size(); ++iProb){
    for(unsigned int i=0; i<5; ++i){
      std::string cutProbe=off+" && "+window+" && probePt>20";
      if(i!=0){
	cutProbe="weight*("+muTag_MC+" && "+cutProbe+" && "+probeMC[iProb]+")"+vtxReWeight;
	hhEtaTmp->Reset();
	t[i]->Project("hhEtaTmp","probeEta",cutProbe.c_str(),"");
	hhEtaTmp->Scale(scale[i],"width");
	hhEta[1][iProb]->Add(hhEtaTmp);
	hhNVtxTmp->Reset();
	t[i]->Project("hhNVtxTmp","nVtx",cutProbe.c_str(),"");
	hhNVtxTmp->Scale(scale[i]);
	hhNVtx[1][iProb]->Add(hhNVtxTmp);
      }
      else{
	cutProbe="weight*("+muTag_data+" && "+cutProbe+" && "+probeData[iProb]+")";
	t[i]->Project(Form("hhEta_%i%i",0,iProb),"probeEta",cutProbe.c_str(),"");
	hhEta[0][iProb]->Scale(scale[i],"width");
	t[i]->Project(Form("hhNVtx_%i%i",0,iProb),"nVtx",cutProbe.c_str(),"");
	hhNVtx[0][iProb]->Scale(scale[i]);
      }
      for(unsigned int j=0; j<region_def.size(); ++j){
	std::string cutEta=muTag_MC+" && "+off+" && "+window;
	if(i==0)
	  cutEta=muTag_data+" && "+off+" && "+window;
	cutEta = cutEta + " && " + region_def[j];
	if(i!=0){
	  cutEta = "weight*("+cutEta+" && "+probeMC[iProb]+")"+vtxReWeight;
	  std::cout<<"\t"<<j<<". cut with eta and probed trigger : "<<cutEta<<std::endl;
	  hhPtTmp->Reset();
	  t[i]->Project("hhPtTmp","probePt",cutEta.c_str(),"");
	  hhPtTmp->Scale(scale[i],"width");
	  hhPt[1][j][iProb]->Add(hhPtTmp);
	}
	else{
	  cutEta = "weight*("+cutEta+" && "+probeData[iProb]+")";
	  std::cout<<"\t"<<j<<". cut with eta and probed trigger : "<<cutEta<<std::endl;
	  t[i]->Project(Form("hhPt_%i%i%i",0,j,iProb),"probePt",cutEta.c_str(),"");
	  hhPt[0][j][iProb]->Scale(scale[i],"width");
	}
      }
    }
  }

  //Data
  //reset overflows
  hNVtx20[0]->SetBinError(hNVtx20[0]->GetNbinsX()+1,0);
  for(unsigned int iProb=0; iProb<probeMC.size(); ++iProb){
    //reset overflows
    hhNVtx[0][iProb]->SetBinContent(hhNVtx[0][iProb]->GetNbinsX()+1,0);
    hhNVtx[0][iProb]->SetBinError(hhNVtx[0][iProb]->GetNbinsX()+1,0);

    //efficiency
    grEtaEff[0][iProb] = new TGraphAsymmErrors(hhEta[0][iProb],hEta20[0],"cl=0.683 b(1,1) mode");
    grNVtxEff[0][iProb] = new TGraphAsymmErrors(hhNVtx[0][iProb],hNVtx20[0],"cl=0.683 b(1,1) mode");
      
    //style
    grEtaEff[0][iProb]->SetLineColor(kBlack);
    grNVtxEff[0][iProb]->SetLineColor(kBlack);
    grEtaEff[0][iProb]->SetMarkerColor(kBlack);
    grNVtxEff[0][iProb]->SetMarkerColor(kBlack);
    grEtaEff[0][iProb]->SetMarkerStyle(20);
    grNVtxEff[0][iProb]->SetMarkerStyle(20);
  }
  //all MC
  TH1D *hEtaMC = (TH1D*)hSEta20[0]->Clone("hEtaMC");
  TH1D *hNVtxMC = (TH1D*)hSNVtx20[0]->Clone("hNVtxMC");
  //reset overflows
  hNVtxMC->SetBinContent(hNVtxMC->GetNbinsX()+1,0);
  hNVtxMC->SetBinError(hNVtxMC->GetNbinsX()+1,0);
    
  for(unsigned int iProb=0; iProb<probeMC.size(); ++iProb){
    //reset overflows
    hhNVtx[1][iProb]->SetBinContent(hhNVtx[1][iProb]->GetNbinsX()+1,0);
    hhNVtx[1][iProb]->SetBinError(hhNVtx[1][iProb]->GetNbinsX()+1,0);
      
    //add hoc safety checks
    for(int ij=0; ij<hEtaMC->GetNbinsX()+2; ++ij){
      double ratio = hEtaMC->GetBinContent(ij)!=0. ? hhEta[1][iProb]->GetBinContent(ij)/hEtaMC->GetBinContent(ij) : -100.*hhEta[1][iProb]->GetBinContent(ij);
      if(ratio>1) {	
	std::cout<<" >1 !"<<std::endl;
	std::cout<<"hhEta[1]["<<iProb<<"]/hEtaMC("<<ij<<") "<<hhEta[1][iProb]->GetBinContent(ij)<<"/"
		 <<hEtaMC->GetBinContent(ij)<<" = "<<ratio<<std::endl;
	if(ratio-1<0.005)//add hoc 'renormalization' for safety
	  hhEta[1][iProb]->SetBinContent(ij,hEtaMC->GetBinContent(ij));
      }
    }
    for(int ij=0; ij<hNVtxMC->GetNbinsX()+2; ++ij){
      double ratio = hNVtxMC->GetBinContent(ij)!=0. ? hhNVtx[1][iProb]->GetBinContent(ij)/hNVtxMC->GetBinContent(ij) : -100.*hhNVtx[1][iProb]->GetBinContent(ij);
      if(ratio>1) {	
	std::cout<<" >1 !"<<std::endl;
	std::cout<<"hhNVtx[1]["<<iProb<<"]/hNVtxMC("<<ij<<") "<<hhNVtx[1][iProb]->GetBinContent(ij)<<"/"
		 <<hNVtxMC->GetBinContent(ij)<<" = "<<ratio<<std::endl;
	if(ratio-1<0.005)//add hoc 'renormalization' for safety
	  hhNVtx[1][iProb]->SetBinContent(ij,hNVtxMC->GetBinContent(ij));
      }
    }

    //efficiency
    grEtaEff[1][iProb] = new TGraphAsymmErrors(hhEta[1][iProb],hEtaMC,"cl=0.683 b(1,1) mode");
    grNVtxEff[1][iProb] = new TGraphAsymmErrors(hhNVtx[1][iProb],hNVtxMC,"cl=0.683 b(1,1) mode");
      
    //style
    grEtaEff[1][iProb]->SetLineColor(kRed);
    grNVtxEff[1][iProb]->SetLineColor(kRed);
    grEtaEff[1][iProb]->SetMarkerColor(kRed);
    grNVtxEff[1][iProb]->SetMarkerColor(kRed);
    grEtaEff[1][iProb]->SetMarkerStyle(21);
    grNVtxEff[1][iProb]->SetMarkerStyle(21);
  }
  /////
  for(unsigned int j=0; j<region_def.size(); ++j){
    //Data
    //reset overflows
    hPt[j][0]->SetBinContent(hPt[j][0]->GetNbinsX()+1,0);
    hPt[j][0]->SetBinError(hPt[j][0]->GetNbinsX()+1,0);
    for(unsigned int iProb=0; iProb<probeMC.size(); ++iProb){
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
      hhPtEff[0][j][iProb]->SetMarkerStyle(20);
      grPtEff[0][j][iProb]->SetLineColor(kBlack);
      grPtEff[0][j][iProb]->SetMarkerColor(kBlack);
      grPtEff[0][j][iProb]->SetMarkerStyle(hhPtEff[0][j][iProb]->GetMarkerStyle());
    }
    //all MC
    TH1D *hPtMC = (TH1D*)hSPt[j][0]->Clone(Form("hPtMC_%i",j));
    //reset overflows
    hPtMC->SetBinContent(hPtMC->GetNbinsX()+1,0);
    hPtMC->SetBinError(hPtMC->GetNbinsX()+1,0);
    
    for(unsigned int iProb=0; iProb<probeMC.size(); ++iProb){
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
  if(lumi<1000.)
    lumi_13TeV = Form("%.1f pb^{-1}",lumi);
  else
    lumi_13TeV = Form("%.2f fb^{-1}",lumi/1000.);
  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 4=13TeV, 7=7+8+13TeV 
  //   iPos = 10*(alignement 1/2/3) + position (1/2/3 = left/center/right)
  int iPos=11;// top-left, left-aligned
  //int iPos=33;// top-right, right-aligned
  //int iPos=22; center, centered

  // Pt 
  for(unsigned int j=0; j<region_def.size(); ++j){
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

  ca->SaveAs((std::string("muEff-RunD_Eta")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("muEff-RunD_Eta")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("muEff-RunD_Eta")+std::string(".png")).c_str());

  if( hEta20[0]->GetBinContent(hEta20[0]->GetMaximumBin()) > hSEta20[0]->GetBinContent(hSEta20[0]->GetMaximumBin()) )
    hEta20[0]->Draw();
  else
    hSEta20[0]->Draw("hist");
  hSEta20[0]->Draw("hist same");
  for(unsigned int i=1; i<4; ++i)
    hSEta20[i]->Draw("same hist");
  hEta20[0]->Draw("same");
  leg->Draw();

  CMS_lumi(ca, iPeriod, iPos);
  ca->Update();
  ca->RedrawAxis();
  ca->GetFrame()->Draw();

  ca->SaveAs((std::string("muEff-RunD_EtaPt20")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("muEff-RunD_EtaPt20")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("muEff-RunD_EtaPt20")+std::string(".png")).c_str());

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

  ca->SaveAs((std::string("muEff-RunD_Mass")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("muEff-RunD_Mass")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("muEff-RunD_Mass")+std::string(".png")).c_str());
  
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

  ca->SaveAs((std::string("muEff-RunD_Mass_Log")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("muEff-RunD_Mass_Log")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("muEff-RunD_Mass_Log")+std::string(".png")).c_str());
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

  ca->SaveAs((std::string("muEff-RunD_NVtx")+std::string(".eps")).c_str());  
  ca->SaveAs((std::string("muEff-RunD_NVtx")+std::string(".pdf")).c_str());
  ca->SaveAs((std::string("muEff-RunD_NVtx")+std::string(".png")).c_str());


  //Effs
  TF1 *myErf[2];
  TFitResultPtr res[2];
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
      for(unsigned int j=0; j<region_def.size(); ++j){
	for(unsigned int iProb=0; iProb<probeMC.size(); ++iProb){
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
  for(unsigned int iProb=0; iProb<probeMC.size(); ++iProb){
    for(unsigned int j=0; j<region_def.size(); ++j){
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
	myErf[0]->SetParameters(17.5,5,5,20,0.95);
	//myErf[0]->SetParameters(17.5,0.2,0.01,1.5,0.95);
	//myErf[0]->SetParameters(17.5,0.2,1.0,5.0,0.95);
	//myErf[0]->SetParLimits(0, 12,23);//turn-on point/threshold (m0)
	myErf[0]->SetParLimits(4, 0.8,1.0);//eff at plateau (norm)
	res[0]= grPtEff[0][j][iProb]->Fit("myErf0","M EX0 S","",18,100);
	//res[0] = grPtEff[0][j][iProb]->Fit("myErf0","M S","",16,100);
	std::cout<<"*********************************"<<std::endl
		 <<"Data: "<<names[iProb]<<", "<<region_name[j]<<std::endl
		 <<"Chi2/Ndf+"<<res[0]->Chi2()<<"/"<<res[0]->Ndf()<<"="<<res[0]->Chi2()/res[0]->Ndf()
		 <<", Prob="<<res[0]->Prob()
		 <<std::endl;
	for(unsigned int iPar=0; iPar<res[0]->NPar(); ++iPar){
	  std::cout<<res[0]->ParName(iPar)<<": "
		   <<res[0]->Parameter(iPar)<<"+-"<<res[0]->ParError(iPar)
		   <<std::endl;
	}
	std::cout<<"*********************************"<<std::endl<<std::endl;

	//myErf[1]->SetParameters(17,0.55,0.81,1.55,0.9);	
	myErf[1]->SetParameters(17.5,5,5,20,0.95);
	//myErf[1]->SetParameters(17.5,0.2,0.01,1.5,0.95);
	//myErf[1]->SetParameters(17.5,0.2,1.0,5.0,0.95);
	//myErf[1]->SetParLimits(0, 12,23);//turn-on point/threshold (m0)
	myErf[1]->SetParLimits(4, 0.8,1.0);//eff at plateau (norm)
	res[1] = grPtEff[1][j][iProb]->Fit("myErf1","M EX0 S","",18,100);
	//res[1] = grPtEff[1][j][iProb]->Fit("myErf1","M S","",16,100);
	std::cout<<"*********************************"<<std::endl
		 <<"MC: "<<names[iProb]<<", "<<region_name[j]<<std::endl
		 <<"Chi2/Ndf="<<res[1]->Chi2()<<"/"<<res[1]->Ndf()<<"="<<res[1]->Chi2()/res[1]->Ndf()
		 <<", Prob="<<res[1]->Prob()
		 <<std::endl;
	for(unsigned int iPar=0; iPar<res[1]->NPar(); ++iPar){
	  std::cout<<res[1]->ParName(iPar)<<": "
		   <<res[1]->Parameter(iPar)<<"+-"<<res[1]->ParError(iPar)
		   <<std::endl;
	}
	std::cout<<"*********************************"<<std::endl<<std::endl;

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

      ca->SaveAs((std::string("muEff-RunD_")+names[iProb]+std::string("_")+region_name[j]+std::string(".eps")).c_str());  
      ca->SaveAs((std::string("muEff-RunD_")+names[iProb]+std::string("_")+region_name[j]+std::string(".pdf")).c_str());  
      ca->SaveAs((std::string("muEff-RunD_")+names[iProb]+std::string("_")+region_name[j]+std::string(".png")).c_str());  
    }
  }
  //Eff in eta
  TH1D *hFrameEta=(TH1D*)hhEta[0][0]->Clone("hFrameEta");
  hFrameEta->Reset();
  hFrameEta->SetMinimum(0.01); 
  hFrameEta->SetMaximum(1.10);
  hFrameEta->SetTitle(";#eta^{#mu offline}; Efficiency");

  for(unsigned int iProb=0; iProb<probeMC.size(); ++iProb){
    hFrameEta->Draw();
    grEtaEff[0][iProb]->Draw("pz same");
    grEtaEff[1][iProb]->Draw("pz same");
    leg2->Clear();
    leg2->AddEntry(grEtaEff[0][iProb],"Data","PLE");
    leg2->AddEntry(grEtaEff[1][iProb],"Simulation","PLE");
    leg2->Draw();

    CMS_lumi(ca, iPeriod, iPos);
    ca->Update();
    ca->RedrawAxis();
    ca->GetFrame()->Draw();

    ca->SaveAs((std::string("muEff-RunD_EtaEff_")+names[iProb]+std::string(".eps")).c_str());
    ca->SaveAs((std::string("muEff-RunD_EtaEff_")+names[iProb]+std::string(".pdf")).c_str());
    ca->SaveAs((std::string("muEff-RunD_EtaEff_")+names[iProb]+std::string(".png")).c_str());
  }

  //Eff in NVtx
  TH1D *hFrameNVtx=(TH1D*)hhNVtx[0][0]->Clone("hFrameNVtx");
  hFrameNVtx->Reset();
  hFrameNVtx->SetMinimum(0.01); 
  hFrameNVtx->SetMaximum(1.10);
  hFrameNVtx->SetTitle(";No. of vertices; Efficiency");

  for(unsigned int iProb=0; iProb<probeMC.size(); ++iProb){
    hFrameNVtx->Draw();
    grNVtxEff[0][iProb]->Draw("pz same");
    grNVtxEff[1][iProb]->Draw("pz same");
    leg2->Clear();
    leg2->AddEntry(grNVtxEff[0][iProb],"Data","PLE");
    leg2->AddEntry(grNVtxEff[1][iProb],"Simulation","PLE");
    leg2->Draw();

    CMS_lumi(ca, iPeriod, iPos);
    ca->Update();
    ca->RedrawAxis();
    ca->GetFrame()->Draw();

    ca->SaveAs((std::string("muEff-RunD_NVtxEff_")+names[iProb]+std::string(".eps")).c_str());
    ca->SaveAs((std::string("muEff-RunD_NVtxEff_")+names[iProb]+std::string(".pdf")).c_str());
    ca->SaveAs((std::string("muEff-RunD_NVtxEff_")+names[iProb]+std::string(".png")).c_str());
  }

  return;
}


