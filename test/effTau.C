#include "TROOT.h"

#include "TMath.h"
#include "TChain.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TPaveText.h"

#include "tdrstyle.C"

#include <string>
#include <vector>
#include <utility>

void effTau(TTree *t=0, bool genMatch=false, float ptMin=25, std::string type="loose");

int main() {
  //gROOT->SetBatch(kTRUE);
  //gROOT->SetStyle("Plain");
  setTDRStyle();
  gStyle->SetPadRightMargin(0.05);

  bool genMatch=true;

  TChain *t = new TChain("muLooseTau/muTauTriggerTree");
  //t->Add("muTauTrgAna*.root");
  t->Add("./prod/Eff_GRunV47_v1/GGH125ToTauTau_Eff_v1/res/muTauTrgAna*.root");
  t->Add("./prod/Eff_GRunV47_v1/VBFH125ToTauTau_Eff_v1/res/muTauTrgAna*.root");
  t->Add("./prod/Eff_GRunV47_v1/ZprimeToTauTau_Eff_v1/res/muTauTrgAna*.root");
  //t->ls();
  effTau(t,genMatch,25,"looseHighPt");

  TChain *t2 = new TChain("muLooseTau/muTauTriggerTree");
  t2->Add("./prod/Eff_GRunV47_v1/GGH125ToTauTau_Eff_Std_v1/res/muTauTrgAna*.root");
  t2->Add("./prod/Eff_GRunV47_v1/VBFH125ToTauTau_Eff_Std_v1/res/muTauTrgAna*.root");
  t2->Add("./prod/Eff_GRunV47_v1/ZprimeToTauTau_Eff_Std_v1/res/muTauTrgAna*.root");
  //t2->ls();
  effTau(t2,genMatch,25,"looseOldHighPt");

  TChain *t3 = new TChain("muMediumTau/muTauTriggerTree");
  t3->Add("./prod/Eff_GRunV47_v1/GGH125ToTauTau_Eff_v1/res/muTauTrgAna*.root");
  t3->Add("./prod/Eff_GRunV47_v1/VBFH125ToTauTau_Eff_v1/res/muTauTrgAna*.root");
  t3->Add("./prod/Eff_GRunV47_v1/ZprimeToTauTau_Eff_v1/res/muTauTrgAna*.root");
  //t3->ls();
  effTau(t3,genMatch,45,"mediumHighPt");

  TChain *t4 = new TChain("muLooseTau/muTauTriggerTree");
  t4->Add("./prod/Eff_GRunV47_v1/GGH125ToTauTau_Eff_v1/res/muTauTrgAna*.root");
  t4->Add("./prod/Eff_GRunV47_v1/VBFH125ToTauTau_Eff_v1/res/muTauTrgAna*.root");
  //t4->ls();
  effTau(t4,genMatch,25,"looseLowPt");

  TChain *t5 = new TChain("muLooseTau/muTauTriggerTree");
  t5->Add("./prod/Eff_GRunV47_v1/GGH125ToTauTau_Eff_Std_v1/res/muTauTrgAna*.root");
  t5->Add("./prod/Eff_GRunV47_v1/VBFH125ToTauTau_Eff_Std_v1/res/muTauTrgAna*.root");
  //t5->ls();
  effTau(t5,genMatch,25,"looseOldLowPt");

  TChain *t6 = new TChain("muMediumTau/muTauTriggerTree");
  t6->Add("./prod/Eff_GRunV47_v1/GGH125ToTauTau_Eff_v1/res/muTauTrgAna*.root");
  t6->Add("./prod/Eff_GRunV47_v1/VBFH125ToTauTau_Eff_v1/res/muTauTrgAna*.root");
  //t6->ls();
  effTau(t6,genMatch,45,"mediumLowPt");

  return 0;
}

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

void effTau(TTree *t, bool genMatch, float ptMin, std::string type){

  if(!t) return;

  //define histos
  //float ptBins[]={0,10,20,30,35,40,45,50,60,70,85,100,120};
  //TH1F *hPt = new TH1F("hPt"," ;p_{T}^{#tau offline} [GeV]; enties/width",12,ptBins);
  //float ptBins[]={0,2,4,6,8,10,12,14,16,18,20,22,24,28,30,35,40,45,50,55,60,70,80,100,120};
  //TH1F *hPt = new TH1F("hPt"," ;p_{T}^{#tau offline} [GeV]; enties/width",24,ptBins);
  float ptBins[]={0,2,4,6,8,10,12,14,16,18,20,22,24,28,30,35,40,45,50,55,60,70,80,100,120,140,170,200,250,300,350,400};
  TH1F *hPt;
  if(type.find("low")!=std::string::npos || type.find("Low")!=std::string::npos)
    hPt = new TH1F("hPt"," ;p_{T}^{#tau offline} [GeV]; enties/width",25,ptBins); //upto 140GeV
  else
    hPt = new TH1F("hPt"," ;p_{T}^{#tau offline} [GeV]; enties/width",31,ptBins);
  hPt->Sumw2();
  hPt->SetMarkerStyle(20);
  //hPt->SetMarkerSize(0.7);
  hPt->SetStats(0);
  hPt->SetMinimum(0.1);
  //TH1F *hLepPt = new TH1F("hPt"," ;p_{T}^{l offline} [GeV]; enties/width",12,ptBins);
  //hLepPt->Sumw2();
  //TH1F *hEta = new TH1F("hEta"," ;#eta^{#tau offline}; enties/width",16,-2.4,2.4);
  TH1F *hEta = new TH1F("hEta"," ;#eta^{#tau offline}; enties/width",23,-2.3,2.3);
  hEta->Sumw2();
  hEta->SetMarkerStyle(20);
  //hEta->SetMarkerSize(0.7);
  hEta->SetStats(0);
  hEta->SetMinimum(0.1);
  //TH1F *hLepEta = new TH1F("hEta"," ;#eta^{l offline}; enties/width",16,-2.4,2.4);
  //hLepEta->Sumw2();
  //TH1F *hPhi = new TH1F("hPhi"," ;#phi^{#tau offline}; enties/width",16,-3.2,3.2);
  TH1F *hPhi = new TH1F("hPhi"," ;#phi^{#tau offline}; enties/width",12,-3.1416,3.1416);
  hPhi->Sumw2();
  hPhi->SetMarkerStyle(20);
  //hPhi->SetMarkerSize(0.7);
  hPhi->SetStats(0);
  hPhi->SetMinimum(0.1);
  //TH1F *hLepPhi = new TH1F("hPhi"," ;#phi^{l offline}; enties/width",16,-3.2,3.2);
  //hLepPhi->Sumw2();

  //
  std::string muTag =
    "(HLT_IsoMu24_eta2p1_IterTrk02_v>0 && hltL3crIsoL1sMu20Eta2p1L1f0L2f20QL3f24QL3crIsoRhoFiltered0p15IterTrk02>0)";
  muTag += "|| (HLT_IsoMu24_eta2p1_v>0 && hltL3crIsoL1sMu20Eta2p1L1f0L2f10QL3f24QL3trkIsoFiltered0p09>0)";
  muTag = "("+muTag+")";
  std::string match = (genMatch ? "muGenMatch && tauGenMatch" : "1");

  std::string sel = muTag + "&&" + match;
  
  std::vector<std::pair<std::string,std::string> > probs;
  if(type.find("loose")!=std::string::npos){
    probs.push_back( std::make_pair("hltOverlapFilterSingleIsoMu17LooseIsoPFTau20>0","LooseIsoPFTau20") );
    probs.push_back( std::make_pair("hltOverlapFilterIsoMu24LooseIsoPFTau20>0","LooseIsoPFTau20_2") );
    probs.push_back( std::make_pair("hltL1sMu16erTauJet20er>0","L1Tau20er") );
    probs.push_back( std::make_pair("hltOverlapFilterIsoMu17LooseIsoPFTau20>0 && hltL1sMu16erTauJet20er>0","L1Tau20erANDLooseIsoPFTau20") );
  }
  unsigned int iMed = probs.size();
  probs.push_back( std::make_pair("hltL1sMu16erIsoTau36er>0","L1IsoTau36er") );
  probs.push_back( std::make_pair("hltL2Tau35eta2p2>0 && hltL1sMu16erIsoTau36er>0","L1IsoTau36erANDL2Tau35er") );
  probs.push_back( std::make_pair("hltL2IsoTau35eta2p2>0 && hltL2Tau35eta2p2>0 && hltL1sMu16erIsoTau36er>0","L1IsoTau36erANDL2Tau35erANDL2IsoTau35er") );
  probs.push_back( std::make_pair("hltOverlapFilterIsoMu17L2IsoTau35>0 && hltL2IsoTau35eta2p2>0 && hltL2Tau35eta2p2>0 && hltL1sMu16erIsoTau36er>0","L1IsoTau36erANDL2Tau35erANDL2IsoTau35er_2") );
  probs.push_back( std::make_pair("hltOverlapFilterIsoMu17MediumIsoPFTau40Reg>0 && hltOverlapFilterIsoMu17L2IsoTau35>0 && hltL2IsoTau35eta2p2>0 && hltL2Tau35eta2p2>0 && hltL1sMu16erIsoTau36er>0","L1IsoTau36erANDL2Tau35erANDL2IsoTau35erANDMediumIsoPFTau40er") );

  //
  t->Project("hPt","tauPt",sel.c_str(),"");
  hPt->Scale(1,"width");
  hPt->SetMinimum(0.1);
  t->Project("hEta","tauEta",(sel+"&& tauPt>"+std::string(Form("%f ",ptMin))).c_str(),"");
  hEta->Scale(1,"width");
  hEta->SetMinimum(0.1);
  t->Project("hPhi","tauPhi",(sel+"&& tauPt>"+std::string(Form("%f ",ptMin))).c_str(),"");
  hPhi->Scale(1,"width");
  hPhi->SetMinimum(0.1);
  std::vector<TH1F*> hhPt;
  std::vector<TH1F*> hhPtEff;
  std::vector<TGraphAsymmErrors*> ggPtEff;
  std::vector<TH1F*> hhEta;
  std::vector<TH1F*> hhEtaEff;
  std::vector<TGraphAsymmErrors*> ggEtaEff;
  std::vector<TH1F*> hhPhi;
  std::vector<TH1F*> hhPhiEff;
  std::vector<TGraphAsymmErrors*> ggPhiEff;

  for(unsigned int i=0; i<probs.size(); ++i){
    //TH1F *hPt_t = (TH1F*)hPt->Clone(Form("hPt_%i",i));
    TH1F *hPt_t = (TH1F*)hPt->Clone((std::string("hPt_")+(probs[i]).second).c_str());
    TH1F *hEta_t = (TH1F*)hEta->Clone((std::string("hEta_")+(probs[i]).second).c_str());
    TH1F *hPhi_t = (TH1F*)hPhi->Clone((std::string("hPhi_")+(probs[i]).second).c_str());
    hPt_t->SetMarkerColor(kRed);
    hPt_t->Reset();
    hEta_t->SetMarkerColor(kRed);
    hEta_t->Reset();
    hPhi_t->SetMarkerColor(kRed);
    hPhi_t->Reset();
    //TH1F *hPtEff_t = (TH1F*)hPt->Clone(Form("hPtEff_%i",i));
    TH1F *hPtEff_t = (TH1F*)hPt->Clone((std::string("hPtEff_")+(probs[i]).second).c_str());
    hPtEff_t->Reset();
    hPtEff_t->SetTitle(";p_{T}^{#tau offline} [GeV]; Efficiency");
    TH1F *hEtaEff_t = (TH1F*)hEta->Clone((std::string("hEtaEff_")+(probs[i]).second).c_str());
    hEtaEff_t->Reset();
    hEtaEff_t->SetTitle(";#eta^{#tau offline}; Efficiency");
    TH1F *hPhiEff_t = (TH1F*)hPhi->Clone((std::string("hPhiEff_")+(probs[i]).second).c_str());
    hPhiEff_t->Reset();
    hPhiEff_t->SetTitle(";#phi^{#tau offline} [GeV]; Efficiency");
    //t->Project(Form("hPt_%i",i),"tauPt",(sel+" && "+probs[i]).c_str(),"");
    t->Project((std::string("hPt_")+(probs[i]).second).c_str(),"tauPt",(sel+" && "+(probs[i]).first).c_str(),"");
    hPt_t->Scale(1,"width");
    hPt_t->SetMinimum(0.1);
    hPtEff_t->Divide(hPt_t,hPt,1,1,"b");
    hPtEff_t->SetMinimum(0.01);
    hPtEff_t->SetMaximum(/*1.10*/1.10);
    TGraphAsymmErrors *gPtEff_t = new TGraphAsymmErrors(hPt_t,hPt,"cl=0.683 b(1,1) median");
    gPtEff_t->SetMarkerStyle(hPt_t->GetMarkerStyle());
    hhPt.push_back(hPt_t);
    hhPtEff.push_back(hPtEff_t);
    ggPtEff.push_back(gPtEff_t);
    t->Project((std::string("hEta_")+(probs[i]).second).c_str(),"tauEta",(sel+"&& tauPt>"+std::string(Form("%f ",ptMin))+" && "+(probs[i]).first).c_str(),"");
    hEta_t->Scale(1,"width");
    hEta_t->SetMinimum(0.1);
    hEtaEff_t->Divide(hEta_t,hEta,1,1,"b");
    hEtaEff_t->SetMinimum(0.01);
    hEtaEff_t->SetMaximum(/*1.10*/1.10);
    TGraphAsymmErrors *gEtaEff_t = new TGraphAsymmErrors(hEta_t,hEta,"cl=0.683 b(1,1) median");
    gEtaEff_t->SetMarkerStyle(hEta_t->GetMarkerStyle());
    hhEta.push_back(hEta_t);
    hhEtaEff.push_back(hEtaEff_t);    
    ggEtaEff.push_back(gEtaEff_t);    
    t->Project((std::string("hPhi_")+(probs[i]).second).c_str(),"tauPhi",(sel+"&& tauPt>"+std::string(Form("%f ",ptMin))+" && "+(probs[i]).first).c_str(),"");
    hPhi_t->Scale(1,"width");
    hPhi_t->SetMinimum(0.1);
    hPhiEff_t->Divide(hPhi_t,hPhi,1,1,"b");
    TGraphAsymmErrors *gPhiEff_t = new TGraphAsymmErrors(hPhi_t,hPhi,"cl=0.683 b(1,1) median");
    hPhiEff_t->SetMinimum(0.01);
    hPhiEff_t->SetMaximum(/*1.10*/1.10);
    gPhiEff_t->SetMarkerStyle(hPhi_t->GetMarkerStyle());
    hhPhi.push_back(hPhi_t);
    hhPhiEff.push_back(hPhiEff_t);    
    ggPhiEff.push_back(gPhiEff_t);    
  }
  //
  TF1 *myErf = new TF1("myErf","[0]*(0.5+0.5*TMath::Erf((x-[2])/(sqrt(2)*[1])))",0,120);
  myErf->SetParNames("eff","res","thr");
  myErf->SetParameters(0.5,10,20);
  myErf->SetLineColor(kRed);

  // Draw all frames on a canvas
  TCanvas *ca = new TCanvas("ca","myEfficiency",1200,600);
  ca->Divide(2);
  //ca->cd(1); gPad->SetLeftMargin(0.15); 
  //ca->cd(2); gPad->SetLeftMargin(0.15); 

  TFile *out = TFile::Open((std::string("tauEff_")+type+std::string("_out.root")).c_str(),"RECREATE");
  out->cd();
  hPt->Write();
  hEta->Write();
  hPhi->Write();

  for(unsigned int i=0; i<probs.size(); ++i){
    //
    hhPt[i]->Write();
    hhPtEff[i]->Write();
    hhEta[i]->Write();
    hhEtaEff[i]->Write();
    hhPhi[i]->Write();
    hhPhiEff[i]->Write();

    //Pt
    ca->cd(1);
    /*hPt->GetYaxis()->SetTitleOffset(1.6);*/ hPt->Draw("e"); 
    /*hhPt[i]->GetYaxis()->SetTitleOffset(1.6);*/ hhPt[i]->Draw("e same");
    
    //
    ca->cd(2);
    //hhPtEff[i]->GetYaxis()->SetTitleOffset(1.6);
    hhPtEff[i]->Draw("e");

    hhPtEff[i]->Fit("myErf","IN","",18,120);
    hhPtEff[i]->Draw("e, same");

    //
    ca->SaveAs((std::string("tauEff_")+(probs[i]).second+type+std::string(".eps")).c_str());  
    ca->SaveAs((std::string("tauEff_")+(probs[i]).second+type+std::string(".png")).c_str());  

    //Eta
    ca->cd(1);
    /*hEta->GetYaxis()->SetTitleOffset(1.6);*/ hEta->Draw("e"); 
    /*hhEta[i]->GetYaxis()->SetTitleOffset(1.6);*/ hhEta[i]->Draw("e same");
    
    //
    ca->cd(2);
    //hhEtaEff[i]->GetYaxis()->SetTitleOffset(1.6);
    hhEtaEff[i]->Draw("e");

    //
    ca->SaveAs((std::string("tauEtaEff_")+(probs[i]).second+type+std::string(".eps")).c_str());  
    ca->SaveAs((std::string("tauEtaEff_")+(probs[i]).second+type+std::string(".png")).c_str());  

    //Phi
    ca->cd(1)->SetGrid();
    /*hPhi->GetYaxis()->SetTitleOffset(1.6);*/ hPhi->Draw("e"); 
    /*hhPhi[i]->GetYaxis()->SetTitleOffset(1.6);*/ hhPhi[i]->Draw("e same");
    
    //
    ca->cd(2)->SetGrid();
    //hhPhiEff[i]->GetYaxis()->SetTitleOffset(1.6);
    hhPhiEff[i]->Draw("e");

    //
    ca->SaveAs((std::string("tauPhiEff_")+(probs[i]).second+type+std::string(".eps")).c_str());  
    ca->SaveAs((std::string("tauPhiEff_")+(probs[i]).second+type+std::string(".png")).c_str());  
  }

  out->Close();

  //additional plots
  TCanvas *cb = new TCanvas("cb","myEfficiency",600,600);
  cb->cd()->SetGrid();
  std::cout<<"Canvas left and right margin: "<<cb->GetLeftMargin()<<", "<<cb->GetRightMargin()<<std::endl;
  //gPad->SetLeftMargin(0.15); 
  TLegend *leg = new TLegend(0.33,0.12,0.88,0.3);
  //TLegend *leg2 = new TLegend(0.75,0.30,0.9,0.40);
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
  TH1F *hFrame=(TH1F*)hhPtEff[0]->Clone("hFrame");
  hFrame->Reset();
  if(type.find("loose")!=std::string::npos){
    //compare LoosePFTau without and with L1
    std::cout<<hhPtEff[0]->GetName()<<" vs "
	     <<hhPtEff[3]->GetName()<<std::endl;

    hhPtEff[0]->Draw("e");
    hhPtEff[3]->SetMarkerColor(kRed);
    hhPtEff[3]->Draw("e same");
    //
    /* Plot with graph
    hFrame->Draw();
    ggPtEff[0]->SetMarkerColor(kBlack);
    ggPtEff[0]->Draw("pz same");
    ggPtEff[3]->SetMarkerColor(kRed);
    ggPtEff[3]->Draw("pz same");
    */
    leg->Clear();
    leg->AddEntry(hhPtEff[3],"L1Tau20er & LooseIsoPFTau20","p");
    leg->AddEntry(hhPtEff[0],"LooseIsoPFTau20","p");
    //leg->Draw();
    leg2->Clear();
    //leg2->AddEntry(hhPtEff[3],"L1 Tau p_{T}>20GeV & HLT Tau p_{T}>20GeV","pe");
    //leg2->AddEntry(hhPtEff[0],"HLT Tau p_{T}>20GeV","pe");
    leg2->AddEntry(hhPtEff[3],"L1+HLT","pe");
    //leg2->AddEntry(hhPtEff[0],"HLT","pe");
    //leg2->AddEntry(hhPtEff[0],"#splitline{HLT}{w/o underlying L1}","pe");
    leg2->AddEntry(hhPtEff[0],"HLT (w/o L1)","pe");
    leg2->Draw();    
    if(type.find("low")!=std::string::npos || type.find("Low")!=std::string::npos)
      //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
      //CMSPrelim("2015, 13TeV, PU20, 25ns                H(125)#rightarrow#tau#tau","#tau_{h} leg of #mu+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      //CMSPrelim("H(125)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #mu+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  
    else
      //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)+Z'(1000)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
      //CMSPrelim("2015, 13TeV, PU20, 25ns   H(125)+Z'(1000)#rightarrow#tau#tau","#tau_{h} leg of #mu+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      //CMSPrelim("H(125)+Z'(1000)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #mu+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

    cb->SaveAs((std::string("tauEff_LoosePFTau20withAndWithoutL1")+type+std::string(".eps")).c_str());
    cb->SaveAs((std::string("tauEff_LoosePFTau20withAndWithoutL1")+type+std::string(".pdf")).c_str());
    cb->SaveAs((std::string("tauEff_LoosePFTau20withAndWithoutL1")+type+std::string(".png")).c_str());
    hhEtaEff[0]->Draw("e");
    hhEtaEff[3]->SetMarkerColor(kRed);
    hhEtaEff[3]->Draw("e same");
    //leg->Draw();
    leg2->Draw();
    if(type.find("low")!=std::string::npos || type.find("Low")!=std::string::npos)
      //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
      //CMSPrelim("2015, 13TeV, PU20, 25ns                H(125)#rightarrow#tau#tau","#tau_{h} leg of #mu+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      //CMSPrelim("H(125)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #mu+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  
    else
      //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)+Z'(1000)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
      //CMSPrelim("2015, 13TeV, PU20, 25ns   H(125)+Z'(1000)#rightarrow#tau#tau","#tau_{h} leg of #mu+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      //CMSPrelim("H(125)+Z'(1000)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #mu+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

    cb->SaveAs((std::string("tauEtaEff_LoosePFTau20withAndWithoutL1")+type+std::string(".eps")).c_str());
    cb->SaveAs((std::string("tauEtaEff_LoosePFTau20withAndWithoutL1")+type+std::string(".pdf")).c_str());
    cb->SaveAs((std::string("tauEtaEff_LoosePFTau20withAndWithoutL1")+type+std::string(".png")).c_str());
    hhPhiEff[0]->Draw("e");
    hhPhiEff[3]->SetMarkerColor(kRed);
    hhPhiEff[3]->Draw("e same");
    //leg->Draw();
    leg2->Draw();
    if(type.find("low")!=std::string::npos || type.find("Low")!=std::string::npos)
      //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
      //CMSPrelim("2015, 13TeV, PU20, 25ns                H(125)#rightarrow#tau#tau","#tau_{h} leg of #mu+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      //CMSPrelim("H(125)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #mu+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  
    else
      //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)+Z'(1000)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
      //CMSPrelim("2015, 13TeV, PU20, 25ns   H(125)+Z'(1000)#rightarrow#tau#tau","#tau_{h} leg of #mu+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      //CMSPrelim("H(125)+Z'(1000)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #mu+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

    cb->SaveAs((std::string("tauPhiEff_LoosePFTau20withAndWithoutL1")+type+std::string(".eps")).c_str());
    cb->SaveAs((std::string("tauPhiEff_LoosePFTau20withAndWithoutL1")+type+std::string(".pdf")).c_str());
    cb->SaveAs((std::string("tauPhiEff_LoosePFTau20withAndWithoutL1")+type+std::string(".png")).c_str());
    
    // "stacked L1 and L1+HLT
    std::cout<<hhPtEff[2]->GetName()<<" vs "
	     <<hhPtEff[3]->GetName()<<std::endl;
    hhPtEff[2]->Draw("e");
    hhPtEff[3]->SetMarkerColor(kRed);
    hhPtEff[3]->Draw("e same");
    leg->Clear();
    leg->AddEntry(hhPtEff[2],"L1Tau20er","p");
    leg->AddEntry(hhPtEff[3],"L1Tau20er & LooseIsoPFTau20","p");
    //leg->Draw();
    leg2->Clear();
    //leg2->AddEntry(hhPtEff[0],"L1 Tau p_{T}>20GeV","pe");
    //leg2->AddEntry(hhPtEff[3],"L1 Tau p_{T}>20GeV & HLT Tau p_{T}>20GeV","pe");
    leg2->AddEntry(hhPtEff[0],"L1","pe");
    leg2->AddEntry(hhPtEff[3],"L1+HLT","pe");
    leg2->Draw();
    if(type.find("low")!=std::string::npos || type.find("Low")!=std::string::npos)
      //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
      //CMSPrelim("2015, 13TeV, PU20, 25ns                H(125)#rightarrow#tau#tau","#tau_{h} leg of e(#mu)+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      //CMSPrelim("H(125)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of e(#mu)+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  
    else
      //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)+Z'(1000)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
      //CMSPrelim("2015, 13TeV, PU20, 25ns   H(125)+Z'(1000)#rightarrow#tau#tau","#tau_{h} leg of e(#mu)+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      //CMSPrelim("H(125)+Z'(1000)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of e(#mu)+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

    cb->SaveAs((std::string("tauEff_LoosePFTau20withL1")+type+std::string(".eps")).c_str());
    cb->SaveAs((std::string("tauEff_LoosePFTau20withL1")+type+std::string(".pdf")).c_str());
    cb->SaveAs((std::string("tauEff_LoosePFTau20withL1")+type+std::string(".png")).c_str());
    hhEtaEff[2]->Draw("e");
    hhEtaEff[3]->SetMarkerColor(kRed);
    hhEtaEff[3]->Draw("e same");
    //leg->Draw();
    leg2->Draw();
    if(type.find("low")!=std::string::npos || type.find("Low")!=std::string::npos)
      //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
      //CMSPrelim("2015, 13TeV, PU20, 25ns                H(125)#rightarrow#tau#tau","#tau_{h} leg of e(#mu)+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      //CMSPrelim("H(125)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of e(#mu)+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

    else
      //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)+Z'(1000)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
      //CMSPrelim("2015, 13TeV, PU20, 25ns   H(125)+Z'(1000)#rightarrow#tau#tau","#tau_{h} leg of e(#mu)+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      //CMSPrelim("H(125)+Z'(1000)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of e(#mu)+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

    cb->SaveAs((std::string("tauEtaEff_LoosePFTau20withL1")+type+std::string(".eps")).c_str());
    cb->SaveAs((std::string("tauEtaEff_LoosePFTau20withL1")+type+std::string(".pdf")).c_str());
    cb->SaveAs((std::string("tauEtaEff_LoosePFTau20withL1")+type+std::string(".png")).c_str());
    hhPhiEff[2]->Draw("e");
    hhPhiEff[3]->SetMarkerColor(kRed);
    hhPhiEff[3]->Draw("e same");
    //leg->Draw();
    leg2->Draw();
    if(type.find("low")!=std::string::npos || type.find("Low")!=std::string::npos)
      //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
      //CMSPrelim("2015, 13TeV, PU20, 25ns                H(125)#rightarrow#tau#tau","#tau_{h} leg of e(#mu)+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      //CMSPrelim("H(125)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of e(#mu)+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

    else
      //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)+Z'(1000)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
      //CMSPrelim("2015, 13TeV, PU20, 25ns   H(125)+Z'(1000)#rightarrow#tau#tau","#tau_{h} leg of e(#mu)+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      //CMSPrelim("H(125)+Z'(1000)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of e(#mu)+#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

    cb->SaveAs((std::string("tauPhiEff_LoosePFTau20withL1")+type+std::string(".eps")).c_str());
    cb->SaveAs((std::string("tauPhiEff_LoosePFTau20withL1")+type+std::string(".pdf")).c_str());
    cb->SaveAs((std::string("tauPhiEff_LoosePFTau20withL1")+type+std::string(".png")).c_str());
  }

  // "stacked L1, and L1+HLT
  hhPtEff[iMed]->Draw("e");
  hhPtEff[iMed+1]->SetMarkerColor(kGreen+3);
  //hhPtEff[iMed+1]->Draw("e same");
  hhPtEff[iMed+3]->SetMarkerColor(kBlue);
  //hhPtEff[iMed+3]->Draw("e same");
  hhPtEff[iMed+4]->SetMarkerColor(kRed);
  hhPtEff[iMed+4]->Draw("e same");
  leg->Clear();
  leg->AddEntry(hhPtEff[iMed],"L1IsoTau36er","p");
  leg->AddEntry(hhPtEff[iMed+1],"L1IsoTau36er & L2Tau35er","p");
  leg->AddEntry(hhPtEff[iMed+3],"L1IsoTau36er & L2Tau35er & L2IsoTau35er","p");
  leg->AddEntry(hhPtEff[iMed+4],"L1IsoTau36er & L2Tau35er & L2IsoTau35er & MediumIsoPFTau40er","p");
  //leg->Draw();
  leg2->Clear();
  //leg2->AddEntry(hhPtEff[iMed],"L1 IsoTau p_{T}>36GeV","pe");
  ////leg2->AddEntry(hhPtEff[iMed+1],"L1 IsoTau p_{T}>36GeV & L2 Tau p_{T}>35GeV","pe");
  //leg2->AddEntry(hhPtEff[iMed+3],"L1 IsoTau p_{T}>36GeV & L2 IsoTau p_{T}>35GeV","pe");
  //leg2->AddEntry(hhPtEff[iMed+4],"L1 IsoTau p_{T}>36GeV & L2 IsoTau p_{T}>35GeV & HLT Tau p_{T}>40GeV","pe");
  leg2->AddEntry(hhPtEff[iMed],"L1","pe");
  //leg2->AddEntry(hhPtEff[iMed+1],"L1+L2","pe");
  //leg2->AddEntry(hhPtEff[iMed+3],"L1+L2+L2.5","pe");
  //leg2->AddEntry(hhPtEff[iMed+3],"L1+L2.5","pe");
  //leg2->AddEntry(hhPtEff[iMed+4],"L1+L2+L2.5+HLT","pe");
  leg2->AddEntry(hhPtEff[iMed+4],"L1+HLT","pe");
  leg2->Draw();
  if(type.find("low")!=std::string::npos || type.find("Low")!=std::string::npos)
    //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
    //CMSPrelim("2015, 13TeV, PU20, 25ns                H(125)#rightarrow#tau#tau","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
    //CMSPrelim("H(125)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

  else
    //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)+Z'(1000)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
    //CMSPrelim("2015, 13TeV, PU20, 25ns   H(125)+Z'(1000)#rightarrow#tau#tau","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
    //CMSPrelim("H(125)+Z'(1000)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

  cb->SaveAs((std::string("tauEff_MediumPFTau40withL1HLT")+type+std::string(".eps")).c_str());
  cb->SaveAs((std::string("tauEff_MediumPFTau40withL1HLT")+type+std::string(".pdf")).c_str());
  cb->SaveAs((std::string("tauEff_MediumPFTau40withL1HLT")+type+std::string(".png")).c_str());
  hhEtaEff[iMed]->Draw("e");
  hhEtaEff[iMed+1]->SetMarkerColor(kGreen+3);
  //hhEtaEff[iMed+1]->Draw("e same");
  hhEtaEff[iMed+3]->SetMarkerColor(kBlue);
  hhEtaEff[iMed+3]->Draw("e same");
  hhEtaEff[iMed+4]->SetMarkerColor(kRed);
  hhEtaEff[iMed+4]->Draw("e same");
  //leg->Draw();
  leg2->Draw();
  if(type.find("low")!=std::string::npos || type.find("Low")!=std::string::npos)
    //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
    //CMSPrelim("2015, 13TeV, PU20, 25ns                H(125)#rightarrow#tau#tau","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
    //CMSPrelim("H(125)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  
  else
    //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)+Z'(1000)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
    //CMSPrelim("2015, 13TeV, PU20, 25ns   H(125)+Z'(1000)#rightarrow#tau#tau","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
    //CMSPrelim("H(125)+Z'(1000)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

  cb->SaveAs((std::string("tauEtaEff_MediumPFTau40withL1HLT")+type+std::string(".eps")).c_str());
  cb->SaveAs((std::string("tauEtaEff_MediumPFTau40withL1HLT")+type+std::string(".pdf")).c_str());
  cb->SaveAs((std::string("tauEtaEff_MediumPFTau40withL1HLT")+type+std::string(".png")).c_str());
  hhPhiEff[iMed]->Draw("e");
  hhPhiEff[iMed+1]->SetMarkerColor(kGreen+3);
  //hhPhiEff[iMed+1]->Draw("e same");
  hhPhiEff[iMed+3]->SetMarkerColor(kBlue);
  hhPhiEff[iMed+3]->Draw("e same");
  hhPhiEff[iMed+4]->SetMarkerColor(kRed);
  hhPhiEff[iMed+4]->Draw("e same");
  //leg->Draw();
  leg2->Draw();
  if(type.find("low")!=std::string::npos || type.find("Low")!=std::string::npos)
    //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
    //CMSPrelim("2015, 13TeV, PU20, 25ns                H(125)#rightarrow#tau#tau","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
    //CMSPrelim("H(125)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  
  else
    //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)+Z'(1000)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
    //CMSPrelim("2015, 13TeV, PU20, 25ns   H(125)+Z'(1000)#rightarrow#tau#tau","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
    //CMSPrelim("H(125)+Z'(1000)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

  cb->SaveAs((std::string("tauPhiEff_MediumPFTau40withL1HLT")+type+std::string(".eps")).c_str());
  cb->SaveAs((std::string("tauPhiEff_MediumPFTau40withL1HLT")+type+std::string(".pdf")).c_str());
  cb->SaveAs((std::string("tauPhiEff_MediumPFTau40withL1HLT")+type+std::string(".png")).c_str());

  // "stacked L1, L1+L2, L1+L2+L2.5 and L1+HLT
  // std::cout<<hhPtEff[iMed]->GetName()<<" vs "
  //	   <<hhPtEff[iMed+1]->GetName()<<" vs "
  //	   <<hhPtEff[iMed+3]->GetName()<<" vs "
  //	   <<hhPtEff[iMed+4]->GetName()<<std::endl;
  hhPtEff[iMed]->Draw("e");
  hhPtEff[iMed+1]->SetMarkerColor(kGreen+3);
  //hhPtEff[iMed+1]->Draw("e same");
  hhPtEff[iMed+3]->SetMarkerColor(kBlue);
  hhPtEff[iMed+3]->Draw("e same");
  hhPtEff[iMed+4]->SetMarkerColor(kRed);
  hhPtEff[iMed+4]->Draw("e same");
  leg->Clear();
  leg->AddEntry(hhPtEff[iMed],"L1IsoTau36er","p");
  leg->AddEntry(hhPtEff[iMed+1],"L1IsoTau36er & L2Tau35er","p");
  leg->AddEntry(hhPtEff[iMed+3],"L1IsoTau36er & L2Tau35er & L2IsoTau35er","p");
  leg->AddEntry(hhPtEff[iMed+4],"L1IsoTau36er & L2Tau35er & L2IsoTau35er & MediumIsoPFTau40er","p");
  //leg->Draw();
  leg2->Clear();
  //leg2->AddEntry(hhPtEff[iMed],"L1 IsoTau p_{T}>36GeV","pe");
  ////leg2->AddEntry(hhPtEff[iMed+1],"L1 IsoTau p_{T}>36GeV & L2 Tau p_{T}>35GeV","pe");
  //leg2->AddEntry(hhPtEff[iMed+3],"L1 IsoTau p_{T}>36GeV & L2 IsoTau p_{T}>35GeV","pe");
  //leg2->AddEntry(hhPtEff[iMed+4],"L1 IsoTau p_{T}>36GeV & L2 IsoTau p_{T}>35GeV & HLT Tau p_{T}>40GeV","pe");
  leg2->AddEntry(hhPtEff[iMed],"L1","pe");
  //leg2->AddEntry(hhPtEff[iMed+1],"L1+L2","pe");
  //leg2->AddEntry(hhPtEff[iMed+3],"L1+L2+L2.5","pe");
  leg2->AddEntry(hhPtEff[iMed+3],"L1+L2.5","pe");
  //leg2->AddEntry(hhPtEff[iMed+4],"L1+L2+L2.5+HLT","pe");
  leg2->AddEntry(hhPtEff[iMed+4],"L1+L2.5+L3","pe");
  leg2->Draw();
  if(type.find("low")!=std::string::npos || type.find("Low")!=std::string::npos)
    //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
    //CMSPrelim("2015, 13TeV, PU20, 25ns                H(125)#rightarrow#tau#tau","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
    //CMSPrelim("H(125)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  
  else
    //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)+Z'(1000)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
    //CMSPrelim("2015, 13TeV, PU20, 25ns   H(125)+Z'(1000)#rightarrow#tau#tau","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
    //CMSPrelim("H(125)+Z'(1000)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

  cb->SaveAs((std::string("tauEff_MediumPFTau40withL1L2L2p5")+type+std::string(".eps")).c_str());
  cb->SaveAs((std::string("tauEff_MediumPFTau40withL1L2L2p5")+type+std::string(".pdf")).c_str());
  cb->SaveAs((std::string("tauEff_MediumPFTau40withL1L2L2p5")+type+std::string(".png")).c_str());
  hhEtaEff[iMed]->Draw("e");
  hhEtaEff[iMed+1]->SetMarkerColor(kGreen+3);
  //hhEtaEff[iMed+1]->Draw("e same");
  hhEtaEff[iMed+3]->SetMarkerColor(kBlue);
  hhEtaEff[iMed+3]->Draw("e same");
  hhEtaEff[iMed+4]->SetMarkerColor(kRed);
  hhEtaEff[iMed+4]->Draw("e same");
  //leg->Draw();
  leg2->Draw();
  if(type.find("low")!=std::string::npos || type.find("Low")!=std::string::npos)
    //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
    //CMSPrelim("2015, 13TeV, PU20, 25ns                H(125)#rightarrow#tau#tau","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
    //CMSPrelim("H(125)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  
  else
    //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)+Z'(1000)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
    //CMSPrelim("2015, 13TeV, PU20, 25ns   H(125)+Z'(1000)#rightarrow#tau#tau","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
    //CMSPrelim("H(125)+Z'(1000)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

  cb->SaveAs((std::string("tauEtaEff_MediumPFTau40withL1L2L2p5")+type+std::string(".eps")).c_str());
  cb->SaveAs((std::string("tauEtaEff_MediumPFTau40withL1L2L2p5")+type+std::string(".pdf")).c_str());
  cb->SaveAs((std::string("tauEtaEff_MediumPFTau40withL1L2L2p5")+type+std::string(".png")).c_str());
  hhPhiEff[iMed]->Draw("e");
  hhPhiEff[iMed+1]->SetMarkerColor(kGreen+3);
  //hhPhiEff[iMed+1]->Draw("e same");
  hhPhiEff[iMed+3]->SetMarkerColor(kBlue);
  hhPhiEff[iMed+3]->Draw("e same");
  hhPhiEff[iMed+4]->SetMarkerColor(kRed);
  hhPhiEff[iMed+4]->Draw("e same");
  //leg->Draw();
  leg2->Draw();
  if(type.find("low")!=std::string::npos || type.find("Low")!=std::string::npos)
    //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
    //CMSPrelim("2015, 13TeV, PU20, 25ns                H(125)#rightarrow#tau#tau","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
    //CMSPrelim("H(125)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  
  else
    //CMSPrelim("2015, 13TeV, PU20, 25ns","H(125)+Z'(1000)#rightarrow#tau#tau","Simulation", 0.10, 0.95);  
    //CMSPrelim("2015, 13TeV, PU20, 25ns   H(125)+Z'(1000)#rightarrow#tau#tau","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
    //CMSPrelim("H(125)+Z'(1000)#rightarrow#tau#tau, 13TeV, PU20, 25ns","#tau_{h} leg of #tau_{h}#tau_{h} trigger", "Simulation", 0.10, 0.95);  
      /*CMSPrelim("2015, 13TeV","", "Simulation", 0.10, 0.95);*/CMSPrelim("2015, 13TeV","", "Simulation ", 0.10, 0.95);  

  cb->SaveAs((std::string("tauPhiEff_MediumPFTau40withL1L2L2p5")+type+std::string(".eps")).c_str());
  cb->SaveAs((std::string("tauPhiEff_MediumPFTau40withL1L2L2p5")+type+std::string(".pdf")).c_str());
  cb->SaveAs((std::string("tauPhiEff_MediumPFTau40withL1L2L2p5")+type+std::string(".png")).c_str());

  /***********/
  // 2012 trigger from Bastian's files
  TFile *f2012 = TFile::Open("turnOnEtaBins_Run2012BCD_forTriggerPaper.root");
  TH1F *hPtEff2012[4];
  for(int ii=0;ii<4;++ii){
    hPtEff2012[ii]=(TH1F*)f2012->Get("denominator_Total")->Clone(Form("hPtEff2012_%i",ii));
    hPtEff2012[ii]->Reset();
    hPtEff2012[ii]->SetTitle(";p_{T}^{#tau offline} [GeV]; Efficiency");
    //hPtEff2012[ii]->Sumw2();
    hPtEff2012[ii]->SetMarkerStyle(24);
  }
  hPtEff2012[0]->Divide((TH1F*)f2012->Get("numeratorL1_Total"),(TH1F*)f2012->Get("denominator_Total"),1,1,"b");
  hPtEff2012[1]->Divide((TH1F*)f2012->Get("numeratorL1L2_Total"),(TH1F*)f2012->Get("denominator_Total"),1,1,"b");
  hPtEff2012[2]->Divide((TH1F*)f2012->Get("numeratorL1L2L2p5_Total"),(TH1F*)f2012->Get("denominator_Total"),1,1,"b");
  hPtEff2012[3]->Divide((TH1F*)f2012->Get("numeratorL1_HLT_Total"),(TH1F*)f2012->Get("denominator_Total"),1,1,"b");
  for(int ii=0;ii<4;++ii){
    hPtEff2012[ii]->SetMinimum(0.01);
    hPtEff2012[ii]->SetMaximum(/*1.10*/1.10);
  }
  
  hPtEff2012[0]->Draw("e");
  hPtEff2012[1]->SetMarkerColor(kGreen+3);
  //hPtEff2012[1]->Draw("e same");
  hPtEff2012[2]->SetMarkerColor(kBlue);
  hPtEff2012[2]->Draw("e same");
  hPtEff2012[3]->SetMarkerColor(kRed);
  hPtEff2012[3]->Draw("e same");
  leg2->Clear();
  leg2->AddEntry(hPtEff2012[0],"L1","pe");
  //leg2->AddEntry(hPtEff2012[1],"L1+L2","pe");
  //leg2->AddEntry(hPtEff2012[2],"L1+L2+L2.5","pe");
  leg2->AddEntry(hPtEff2012[2],"L1+L2.5","pe");
  //leg2->AddEntry(hPtEff2012[3],"L1+L2+L2.5+HLT","pe");
  leg2->AddEntry(hPtEff2012[3],"L1+L2.5+L3","pe");
  leg2->Draw();
  CMSPrelim("2012, 8TeV, L=20/fb","","", 0.10, 0.95);  
  cb->SaveAs((std::string("tauEff2012_MediumPFTau35withL1L2L2p5")+type+std::string(".eps")).c_str());
  cb->SaveAs((std::string("tauEff2012_MediumPFTau35withL1L2L2p5")+type+std::string(".pdf")).c_str());
  cb->SaveAs((std::string("tauEff2012_MediumPFTau35withL1L2L2p5")+type+std::string(".png")).c_str());

  hhPtEff[iMed]->Draw("e");
  hhPtEff[iMed+1]->SetMarkerColor(kGreen+3);
  //hhPtEff[iMed+1]->Draw("e same");
  hhPtEff[iMed+3]->SetMarkerColor(kBlue);
  hhPtEff[iMed+3]->Draw("e same");
  hhPtEff[iMed+4]->SetMarkerColor(kRed);
  hhPtEff[iMed+4]->Draw("e same");
  hPtEff2012[0]->Draw("e same");
  hPtEff2012[1]->SetMarkerColor(kGreen+3);
  //hPtEff2012[1]->Draw("e same");
  hPtEff2012[2]->SetMarkerColor(kBlue);
  hPtEff2012[2]->Draw("e same");
  hPtEff2012[3]->SetMarkerColor(kRed);
  hPtEff2012[3]->Draw("e same");
  leg2->Clear();
  leg2->AddEntry((TObject*)0,"L1","");
  leg2->AddEntry(hhPtEff[iMed],"2015","pe");
  leg2->AddEntry(hPtEff2012[0],"2012","pe");
  leg2->AddEntry((TObject*)0,"L1+L2.5","");
  leg2->AddEntry(hhPtEff[iMed+3],"2015","pe");
  leg2->AddEntry(hPtEff2012[2],"2012","pe");
  leg2->AddEntry((TObject*)0,"L1+L2.5+L3","");
  leg2->AddEntry(hhPtEff[iMed+4],"2015","pe");
  leg2->AddEntry(hPtEff2012[3],"2012","pe");
  leg2->Draw();
  CMSPrelim("2012 data vs 2015 Simulation","","", 0.10, 0.95);  
  cb->SaveAs((std::string("tauEff2015vs2012_MediumPFTauWithL1L2L2p5")+type+std::string(".eps")).c_str());
  cb->SaveAs((std::string("tauEff2015vs2012_MediumPFTauWithL1L2L2p5")+type+std::string(".pdf")).c_str());
  cb->SaveAs((std::string("tauEff2015vs2012_MediumPFTauWithL1L2L2p5")+type+std::string(".png")).c_str());

  hhPtEff[iMed+4]->Draw("e");
  hPtEff2012[3]->Draw("e same");
  leg2->Clear();
  leg2->AddEntry((TObject*)0,"L1+L2.5+L3","");
  leg2->AddEntry(hhPtEff[iMed+4],"2015","pe");
  leg2->AddEntry(hPtEff2012[3],"2012","pe");
  leg2->Draw();
  CMSPrelim("2012 data vs 2015 Simulation","","", 0.10, 0.95);  
  cb->SaveAs((std::string("tauEff2015vs2012_MediumPFTau")+type+std::string(".eps")).c_str());
  cb->SaveAs((std::string("tauEff2015vs2012_MediumPFTau")+type+std::string(".pdf")).c_str());
  cb->SaveAs((std::string("tauEff2015vs2012_MediumPFTau")+type+std::string(".png")).c_str());

  f2012->Close();
  /*********/

  delete hPt;
  delete hEta;
  delete hPhi;
  delete myErf;
  for(unsigned int i=0; i<hhPt.size(); ++i)
    delete hhPt[i];
  for(unsigned int i=0; i<hhPtEff.size(); ++i)  
    delete hhPtEff[i];
  for(unsigned int i=0; i<hhEta.size(); ++i)
    delete hhEta[i];
  for(unsigned int i=0; i<hhEtaEff.size(); ++i)
    delete hhEtaEff[i];
  for(unsigned int i=0; i<hhPhi.size(); ++i)
    delete hhPhi[i];
  for(unsigned int i=0; i<hhPhiEff.size(); ++i)
    delete hhPhiEff[i];
  for(unsigned int i=0; i<ggPtEff.size(); ++i)
    delete ggPtEff[i];
  for(unsigned int i=0; i<ggEtaEff.size(); ++i)
    delete ggEtaEff[i];
  for(unsigned int i=0; i<ggPhiEff.size(); ++i)
    delete ggPhiEff[i];
  delete ca;
  delete cb;
  delete leg;
  delete leg2;
}


/*
*/ 

