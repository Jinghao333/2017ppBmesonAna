#include "TROOT.h"
#include "TStyle.h"
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom.h"
#include <iostream>
#include <iomanip>
#include <fstream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"

#include "TH1.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TLine.h"
#include "../../parameter.h"



using namespace std;

void RAAratio(){

  gSystem->mkdir("RAAPlots" ,true );
    gStyle->SetOptStat(0);

	TFile * FileBP = new TFile("../../BP/EffAna/FinalFiles/BPPPCorrYieldPT.root");
	TFile * FileBs = new TFile("../../Bs/EffAna/FinalFiles/BsPPCorrYieldPT.root");
	TH1D * BPCross = (TH1D *) FileBP->Get("hPtSigma");
	TH1D * BsCross = (TH1D *) FileBs->Get("hPtSigma");

	const int NBins = ptbinsvec.size() - 1;
	const int NBins_low = 1;
	const int NBins_high = 3;

	float BsXsecPPY[NBins];
    float BPXsecPPY[NBins];
	    float BsXSecPPYErrUp[NBins];
	    float BsXSecPPYErrDown[NBins];
        float BPXSecPPYErrUp[NBins];
	    float BPXSecPPYErrDown[NBins];
        float BPXSecPPYErrUpRatio[NBins];
		float BPXSecPPYErrDownRatio[NBins];
	    float BsXSecPPYErrUpPercent[NBins];
	    float BsXSecPPYErrDownPercent[NBins];

	float BXsecX[NBins];
	    float BXSecXErrUp[NBins];
	    float BXSecXErrDown[NBins];

	for(int i = 0; i < NBins; i++){
        //Bin X position (CENTERED for now)
        BXsecX[i] = (ptbinsvec[i] + ptbinsvec[i + 1]) / 2;
        BXSecXErrUp[i] = (ptbinsvec[i+1] - ptbinsvec[i]) / 2;
        BXSecXErrDown[i] = BXSecXErrUp[i];
        //Bin X position (CENTERED for now)
        //BS
        BsXsecPPY[i] = BsCross->GetBinContent(i+1);
		BsXSecPPYErrUp[i] = BsCross->GetBinError(i+1);
		BsXSecPPYErrDown[i] = BsCross->GetBinError(i+1);
		BsXSecPPYErrUpPercent[i] = BsXSecPPYErrUp[i]/BsXsecPPY[i];
		BsXSecPPYErrDownPercent[i] = BsXSecPPYErrDown[i]/BsXsecPPY[i];
        //Bs
        //BP
		BPXsecPPY[i] = BPCross->GetBinContent(i+1);
		BPXSecPPYErrUp[i] = BPCross->GetBinError(i+1);
		BPXSecPPYErrDown[i] = BPCross->GetBinError(i+1);
        BPXSecPPYErrUpRatio[i] = BPXSecPPYErrUp[i]/BPXsecPPY[i];
		BPXSecPPYErrDownRatio[i] = BPXSecPPYErrDown[i]/BPXsecPPY[i];
        //BP
	}

//UNCERTANTIES

  TFile fError_bs("../../2DMapSyst/OutFiles/BsError2D.root");
  TFile fError("../../2DMapSyst/OutFiles/BPError2D.root");

	TH1D * TnPSyst_bs = (TH1D *) fError_bs.Get("TnPSyst");
	TH1D * BptSyst_bs = (TH1D *) fError_bs.Get("BptSyst");
	TH1D * MCDataSyst_bs = (TH1D *) fError_bs.Get("MCDataSyst");
  	if (!MCDataSyst_bs) MCDataSyst_bs = (TH1D *) fError_bs.Get("BDTSyst");
  TFile fPdfError_bs("../../bs_pdf.root");
  TGraph* pdfSyst_bs = (TGraph *) fPdfError_bs.Get("bs_error");
  TFile fTrackSelError_bs("../../syst_track_sel.root");
  TGraph* trackSelSyst_bs = (TGraph *) fTrackSelError_bs.Get("bs_track_sel_error");

	TH1D * TnPSyst = (TH1D *) fError.Get("TnPSyst");
	TH1D * BptSyst = (TH1D *) fError.Get("BptSyst");
	TH1D * MCDataSyst = (TH1D *) fError.Get("MCDataSyst");
  if (!MCDataSyst) MCDataSyst = (TH1D *) fError.Get("BDTSyst");
  TFile fPdfError("../../bp_pdf.root");
  TGraph* pdfSyst = (TGraph *) fPdfError.Get("bp_error");
  TFile fTrackSelError("../../syst_track_sel.root");
  TGraph* trackSelSyst = (TGraph *) fTrackSelError.Get("bp_track_sel_error");

  // percent error
	float BPTrackingSyst[NBins] = {[0 ... NBins - 1] = 2.4};
	float BsTrackingSyst[NBins] = {[0 ... NBins - 1] = 4.8};

	float BPMCDataSyst[NBins];
	float BPPDFSyst[NBins];
	float BPTrackSelSyst[NBins];
	float BPPtShapeSyst[NBins];
	float BPTnPSystDown[NBins];
	float BPTnPSystUp[NBins];
	float BsMCDataSyst[NBins];
	float BsPtShapeSyst[NBins];
	float BsPDFSyst[NBins];
	float BsTrackSelSyst[NBins];
	float BsTnPSystDown[NBins];
	float BsTnPSystUp[NBins];

  // Get systematics from input files
    for (auto ibin = 0; ibin < NBins; ++ibin) {
    BsMCDataSyst[ibin] = MCDataSyst_bs->GetBinContent(ibin + 1);
    BsPtShapeSyst[ibin] = BptSyst_bs->GetBinContent(ibin + 1);
    BsTnPSystDown[ibin] = TnPSyst_bs->GetBinContent(ibin + 1);
    BPMCDataSyst[ibin] = MCDataSyst->GetBinContent(ibin + 1);
    BPPtShapeSyst[ibin] = BptSyst->GetBinContent(ibin + 1);
    BPTnPSystDown[ibin] = TnPSyst->GetBinContent(ibin + 1);
    // TnP systematics are symmetric in the binned pT case
    BsTnPSystUp[ibin] = BsTnPSystDown[ibin];
    BsPDFSyst[ibin] = pdfSyst_bs->GetY()[ibin];
    BsTrackSelSyst[ibin] = trackSelSyst_bs->GetY()[ibin];  
    BPTnPSystUp[ibin] = BPTnPSystDown[ibin];
    BPPDFSyst[ibin] = pdfSyst->GetY()[ibin];
    BPTrackSelSyst[ibin] = trackSelSyst->GetY()[ibin];
    
  }

  // RMS of all the errors
	float BsTotalSystDownRatio[NBins];
	float BsTotalSystUpRatio[NBins];
	float BPTotalSystDownRatio[NBins];
	float BPTotalSystUpRatio[NBins];

    float BP_TOTAL_SYST_DoubleRAA[NBins];
    float Bs_TOTAL_SYST_DoubleRAA[NBins];

	for(int i = 0; i < NBins; i++){
		BsTotalSystDownRatio[i] = TMath::Sqrt(TMath::Power(BsTrackingSyst[i], 2) + TMath::Power(BsMCDataSyst[i], 2) +
                                          TMath::Power(BsPDFSyst[i], 2) + TMath::Power(BsTrackSelSyst[i], 2) +
                                          TMath::Power(BsPtShapeSyst[i], 2) + TMath::Power(BsTnPSystDown[i], 2)) / 100;
		BsTotalSystUpRatio[i] = BsTotalSystDownRatio[i];

		BPTotalSystDownRatio[i] = TMath::Sqrt(TMath::Power(BPTrackingSyst[i], 2) + TMath::Power(BPMCDataSyst[i], 2) +
                                          TMath::Power(BPPDFSyst[i], 2) + TMath::Power(BPTrackSelSyst[i], 2) +
                                          TMath::Power(BPPtShapeSyst[i], 2) + TMath::Power(BPTnPSystDown[i], 2)) / 100;
		BPTotalSystUpRatio[i] = BPTotalSystDownRatio[i];


//for Bs one of the track unc survives -> use bp unc which only has 1 track
        Bs_TOTAL_SYST_DoubleRAA[i] = TMath::Sqrt(/*TMath::Power(BPTrackingSyst[i], 2)+*/ TMath::Power(BsMCDataSyst[i], 2) +
                                          TMath::Power(BsPDFSyst[i], 2) + TMath::Power(BsTrackSelSyst[i], 2) +
                                          TMath::Power(BsPtShapeSyst[i], 2) + TMath::Power(BsTnPSystDown[i], 2)) / 100;
        BP_TOTAL_SYST_DoubleRAA[i] = TMath::Sqrt( TMath::Power(BPMCDataSyst[i], 2) +
                                          TMath::Power(BPPDFSyst[i], 2) + TMath::Power(BPTrackSelSyst[i], 2) +
                                          TMath::Power(BPPtShapeSyst[i], 2) + TMath::Power(BPTnPSystDown[i], 2)) / 100;
	}


//UNCERTANTIES


// RAA VALUES    
//x position of the RAA BIN
	float RAAX_high[NBins_high] = {BXsecX[1],BXsecX[2],BXsecX[3]};
	float RAAX_low[NBins_low] = {BXsecX[0]};
	float RAAXErrUp_high[NBins_high] = {BXSecXErrUp[1],BXSecXErrUp[2],BXSecXErrUp[3]};
	float RAAXErrDown_high[NBins_high] = {BXSecXErrUp[1],BXSecXErrUp[2],BXSecXErrUp[3]};
	float RAAXErrUp_low[NBins_low] = {BXSecXErrUp[0]};
	float RAAXErrDown_low[NBins_low] = {BXSecXErrUp[0]};
//x position of the RAA BIN


    float BsRAAY[NBins];
	float BsRAAYErrUp[NBins];
	float BsRAAYErrDown[NBins];
	float BsRAAYSystDown[NBins];
	float BsRAAYSystUp[NBins];
    float BsRAAYErrUpRatio[NBins];
    float BsRAAYErrDownRatio[NBins];
    float BsRAAYSystUpRatio[NBins];
    float BsRAAYSystDownRatio[NBins];

    float BPRAAY[NBins];
	float BPRAAYErrUp[NBins];
	float BPRAAYErrDown[NBins];
    float BPRAAYSystUp[NBins] ;
	float BPRAAYSystDown[NBins];
	float BPRAAYSystUpRatio[NBins] ;
	float BPRAAYSystDownRatio[NBins];
	float BPRAAYErrUpRatio[NBins];
	float BPRAAYErrDownRatio[NBins];

    float BPRAAYSystUp_DR;
    float BPRAAYSystDown_DR;
    float BsRAAYSystUp_DR;
    float BsRAAYSystDown_DR;

    float BPRAAYSyst_upDR[NBins];
    float BsRAAYSyst_upDR[NBins];
    float BPRAAYSyst_downDR[NBins];
    float BsRAAYSyst_downDR[NBins];

	for(int i = 0; i < NBins; i++){

		BPRAAY[i] = BPXsecPbPbY[i] / BPXsecPPY[i];
		BsRAAY[i] = BsXsecPbPbY[i] / BsXsecPPY[i];

		BPRAAYErrUp[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYErrUpRatio[i] * BPXSecPbPbYErrUpRatio[i] + BPXSecPPYErrUpRatio[i] * BPXSecPPYErrUpRatio[i]);
		BPRAAYErrDown[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYErrDownRatio[i] * BPXSecPbPbYErrDownRatio[i] + BPXSecPPYErrDownRatio[i] * BPXSecPPYErrDownRatio[i]);
		BPRAAYSystDown[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYSystDownRatio[i] * BPXSecPbPbYSystDownRatio[i] + BPTotalSystDownRatio[i] * BPTotalSystDownRatio[i]);
		BPRAAYSystUp[i] = BPRAAY[i] * TMath::Sqrt(BPXSecPbPbYSystUpRatio[i] * BPXSecPbPbYSystUpRatio[i] + BPTotalSystUpRatio[i] * BPTotalSystUpRatio[i]);
        BPRAAYErrUpRatio[i] = BPRAAYErrUp[i]/BPRAAY[i];
        BPRAAYErrDownRatio[i] = BPRAAYErrDown[i]/BPRAAY[i];
        BPRAAYSystUpRatio[i] = BPRAAYSystUp[i]/BPRAAY[i];
        BPRAAYSystDownRatio[i] = BPRAAYSystDown[i]/BPRAAY[i];

		BsRAAYErrUp[i] = BsRAAY[i] * TMath::Sqrt(BsXSecPbPbYErrUpPercent[i] * BsXSecPbPbYErrUpPercent[i] + BsXSecPPYErrUpPercent[i] * BsXSecPPYErrUpPercent[i]);
		BsRAAYErrDown[i]  = BsRAAY[i] * TMath::Sqrt(BsXSecPbPbYErrDownPercent[i] * BsXSecPbPbYErrDownPercent[i] + BsXSecPPYErrDownPercent[i] * BsXSecPPYErrDownPercent[i]);
		BsRAAYSystDown[i] = BsRAAY[i] * TMath::Sqrt(BsXSecPbPbYSystDownPercent[i] * BsXSecPbPbYSystDownPercent[i] + BsTotalSystDownRatio[i] * BsTotalSystDownRatio[i]);
		BsRAAYSystUp[i]   = BsRAAY[i] * TMath::Sqrt(BsXSecPbPbYSystUpPercent[i] * BsXSecPbPbYSystUpPercent[i] + BsTotalSystUpRatio[i] * BsTotalSystUpRatio[i]);
        BsRAAYErrUpRatio[i] = BsRAAYErrUp[i]/BsRAAY[i];
        BsRAAYErrDownRatio[i] = BsRAAYErrDown[i]/BsRAAY[i];
        BsRAAYSystUpRatio[i] = BsRAAYSystUp[i]/BsRAAY[i];
        BsRAAYSystDownRatio[i] = BsRAAYSystDown[i]/BsRAAY[i];

//for double ratios
		BPRAAYSyst_upDR[i] = BPRAAY[i] * TMath::Sqrt( BPXSecPbPbYSystUp_DR_Percent[i] * BPXSecPbPbYSystUp_DR_Percent[i] + BP_TOTAL_SYST_DoubleRAA[i] * BP_TOTAL_SYST_DoubleRAA[i]);
		BsRAAYSyst_upDR[i] = BsRAAY[i] * TMath::Sqrt( BsXSecPbPbYSystUp_DR_Percent[i] * BsXSecPbPbYSystUp_DR_Percent[i] + Bs_TOTAL_SYST_DoubleRAA[i] * Bs_TOTAL_SYST_DoubleRAA[i]);
        BPRAAYSyst_downDR[i] = BPRAAY[i] * TMath::Sqrt( BPXSecPbPbYSystDown_DR_Percent[i] * BPXSecPbPbYSystDown_DR_Percent[i] + BP_TOTAL_SYST_DoubleRAA[i] * BP_TOTAL_SYST_DoubleRAA[i]);
		BsRAAYSyst_downDR[i] = BsRAAY[i] * TMath::Sqrt( BsXSecPbPbYSystDown_DR_Percent[i] * BsXSecPbPbYSystDown_DR_Percent[i] + Bs_TOTAL_SYST_DoubleRAA[i] * Bs_TOTAL_SYST_DoubleRAA[i]);
	}
// RAA VALUES  

//DIVE INTO BINS ACCORDING TO FID REGION
float BsRAAY_low[NBins_low] = {BsRAAY[0]};
float BsRAAY_high[NBins_high] = {BsRAAY[1],BsRAAY[2],BsRAAY[3]};
float BsRAAYErrDown_high[NBins_high] = {BsRAAYErrDown[1],BsRAAYErrDown[2],BsRAAYErrDown[3]};
float BsRAAYErrDown_low[NBins_low] = {BsRAAYErrDown[0]};
float BsRAAYSystDown_high[NBins_high] = {BsRAAYSystDown[1],BsRAAYSystDown[2],BsRAAYSystDown[3]};
float BsRAAYSystDown_low[NBins_low] = {BsRAAYSystDown[0]};
float BsRAAYErrUp_high[NBins_high] = {BsRAAYErrUp[1],BsRAAYErrUp[2],BsRAAYErrUp[3]};
float BsRAAYErrUp_low[NBins_low] = {BsRAAYErrUp[0]};
float BsRAAYSystUp_high[NBins_high] = {BsRAAYSystUp[1],BsRAAYSystUp[2],BsRAAYSystUp[3]};
float BsRAAYSystUp_low[NBins_low] = {BsRAAYSystUp[0]};
float BPRAAY_low[NBins_low] = {BPRAAY[0]};
float BPRAAY_high[NBins_high] = {BPRAAY[1],BPRAAY[2],BPRAAY[3]};
float BPRAAYErrDown_high[NBins_high] = {BPRAAYErrDown[1],BPRAAYErrDown[2],BPRAAYErrDown[3]};
float BPRAAYErrDown_low[NBins_low] = {BPRAAYErrDown[0]};
float BPRAAYSystDown_high[NBins_high] = {BPRAAYSystDown[1],BPRAAYSystDown[2],BPRAAYSystDown[3]};
float BPRAAYSystDown_low[NBins_low] = {BPRAAYSystDown[0]};
float BPRAAYErrUp_high[NBins_high] = {BPRAAYErrUp[1],BPRAAYErrUp[2],BPRAAYErrUp[3]};
float BPRAAYErrUp_low[NBins_low] = {BPRAAYErrUp[0]};
float BPRAAYSystUp_high[NBins_high] = {BPRAAYSystUp[1],BPRAAYSystUp[2],BPRAAYSystUp[3]};
float BPRAAYSystUp_low[NBins_low] = {BPRAAYSystUp[0]};

float zero1[1] = {0};


//PLOTS PLOTS PLOTS
	TCanvas * c2 = new TCanvas("c","c",600,600);
	c2->cd();
	c2->SetLeftMargin(0.16);

	TH2D * HisEmptyRAA = new TH2D("HisEmptyRAA","",100,7,50,100,0,1.2);
	HisEmptyRAA->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	HisEmptyRAA->GetYaxis()->SetTitle("R_{AA} = #frac{1}{T_{AA}} #frac{dN_{PbPb}/dp_{T}}{d#sigma_{pp}/dp_{T}}");
	HisEmptyRAA->GetXaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->CenterTitle();
	HisEmptyRAA->GetYaxis()->SetTitleOffset(1.8);

    TH2D * HisEmptyRAA_bs = new TH2D("HisEmptyRAA","",100,7,50,100,0,2.2);
	HisEmptyRAA_bs->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	HisEmptyRAA_bs->GetYaxis()->SetTitle("R_{AA} = #frac{1}{T_{AA}} #frac{dN_{PbPb}/dp_{T}}{d#sigma_{pp}/dp_{T}}");
	HisEmptyRAA_bs->GetXaxis()->CenterTitle();
	HisEmptyRAA_bs->GetYaxis()->CenterTitle();
	HisEmptyRAA_bs->GetYaxis()->SetTitleOffset(1.8);

    TH2D * HisEmptyRAA_bs_bp = new TH2D("HisEmptyRAA","",100,7,50,100,0,1.8);
	HisEmptyRAA_bs_bp->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	HisEmptyRAA_bs_bp->GetYaxis()->SetTitle("R_{AA} = #frac{1}{T_{AA}} #frac{dN_{PbPb}/dp_{T}}{d#sigma_{pp}/dp_{T}}");
	HisEmptyRAA_bs_bp->GetXaxis()->CenterTitle();
	HisEmptyRAA_bs_bp->GetYaxis()->CenterTitle();
	HisEmptyRAA_bs_bp->GetYaxis()->SetTitleOffset(1.8);

    TH2D * HisEmptyRAAr = new TH2D("HisEmptyRAA","",100,7,50,100, 0.5 ,3.8);
	HisEmptyRAAr->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	HisEmptyRAAr->GetYaxis()->SetTitle("R_{AA}^{B^{0}_{s}} / R_{AA}^{B^{+}}");
	HisEmptyRAAr->GetXaxis()->CenterTitle();
	HisEmptyRAAr->GetYaxis()->CenterTitle();
	HisEmptyRAAr->GetYaxis()->SetTitleOffset(1.8);

    TH2D * HisEmptyRAA_Dr = new TH2D("HisEmptyRAA","",100,7,50,100, 0.5 ,6);
	HisEmptyRAA_Dr->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	HisEmptyRAA_Dr->GetYaxis()->SetTitle("R_{AA}^{B^{0}_{s}} / R_{AA}^{B^{+}}");
	HisEmptyRAA_Dr->GetXaxis()->CenterTitle();
	HisEmptyRAA_Dr->GetYaxis()->CenterTitle();
	HisEmptyRAA_Dr->GetYaxis()->SetTitleOffset(1.8);

    TLine * Unity = new TLine(7,1,50,1);
	Unity->SetLineWidth(2);
	Unity->SetLineStyle(2);
	Unity->SetLineColor(1);

	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.025); 
	lat->SetTextFont(42);

//PLOTS PLOTS PLOTS

// DRAW B+ RAA vs 2015 RAA
//  B+ RAA  
TGraphAsymmErrors *BPRAAGraph = new TGraphAsymmErrors(NBins_high, RAAX_high, BPRAAY_high, RAAXErrDown_high, RAAXErrUp_high, BPRAAYErrDown_high,BPRAAYErrUp_high);
TGraphAsymmErrors *BPRAAGraphSyst = new TGraphAsymmErrors(NBins_high, RAAX_high, BPRAAY_high, RAAXErrDown_high, RAAXErrUp_high, BPRAAYSystDown_high,BPRAAYSystUp_high);
TGraphAsymmErrors *BPRAAGraph_low = new TGraphAsymmErrors(NBins_low, RAAX_low, BPRAAY_low, RAAXErrDown_low, RAAXErrUp_low,BPRAAYErrDown_low, BPRAAYErrUp_low);
TGraphAsymmErrors *BPRAAGraph_low_just_marker = new TGraphAsymmErrors(NBins_low, RAAX_low, BPRAAY_low, zero1, zero1,zero1,zero1);
TGraphAsymmErrors *BPRAAGraphSyst_low = new TGraphAsymmErrors(NBins_low, RAAX_low, BPRAAY_low, RAAXErrDown_low, RAAXErrUp_low,BPRAAYSystDown_low, BPRAAYSystUp_low);
//  B+ RAA
//  2015
TGraphAsymmErrors *BPRAAGraph2015 = new TGraphAsymmErrors(NBins2015, BPRAAX2015, BPRAAY2015,BPRAAXErrDown2015, BPRAAXErrUp2015,BPRAAYErrDown2015,BPRAAYErrUp2015);
TGraphAsymmErrors *BPRAAGraphSyst2015 = new TGraphAsymmErrors(NBins2015, BPRAAX2015, BPRAAY2015,BPRAAXErrDown2015, BPRAAXErrUp2015,BPRAAYSystDown2015,BPRAAYSystUp2015);
//  2015
BPRAAGraph_low_just_marker ->SetMarkerStyle(20);
BPRAAGraph_low_just_marker ->SetMarkerSize(0.9);
BPRAAGraph_low_just_marker ->SetMarkerColor(kWhite);
BPRAAGraph_low ->SetMarkerColor(kGreen+2);
BPRAAGraph_low ->SetLineColor(kGreen+2);
BPRAAGraph_low ->SetMarkerStyle(24);
BPRAAGraph_low ->SetMarkerSize(1);
BPRAAGraphSyst_low ->SetFillColorAlpha(kGreen-7,0.5);
BPRAAGraphSyst_low ->SetLineColor(kGreen-7);

BPRAAGraph ->SetMarkerColor(kGreen+2);
BPRAAGraph ->SetLineColor(kGreen+2);
BPRAAGraph ->SetMarkerStyle(20);
BPRAAGraph ->SetMarkerSize(1);
BPRAAGraphSyst ->SetFillColorAlpha(kGreen-7,0.5);
BPRAAGraphSyst ->SetLineColor(kGreen-7);

BPRAAGraph2015->SetLineColor(kOrange+1);
BPRAAGraph2015->SetMarkerColor(kOrange+1);	
BPRAAGraph2015->SetMarkerStyle(21);
BPRAAGraph2015->SetMarkerSize(1);
BPRAAGraphSyst2015->SetFillColorAlpha(kOrange+1,0.5);
BPRAAGraphSyst2015->SetLineColor(kOrange+1);


	HisEmptyRAA->Draw();
	BPRAAGraphSyst_low->Draw("5same");
	BPRAAGraphSyst->Draw("5same");
	BPRAAGraph->Draw("ep");
	BPRAAGraph_low->Draw("epSAME");
	BPRAAGraph_low_just_marker->Draw("epSAME");
	Unity->Draw("SAME");

	c2->SaveAs("RAAPlots/BPRAA.pdf");

	HisEmptyRAA->Draw();
	BPRAAGraphSyst2015->Draw("5same");
	BPRAAGraph2015->Draw("epSAME");
	BPRAAGraphSyst_low->Draw("5same");
	BPRAAGraphSyst->Draw("5same");
	BPRAAGraph->Draw("epSAME");
	BPRAAGraph_low->Draw("epSAME");
	BPRAAGraph_low_just_marker->Draw("epSAME");
	Unity->Draw("SAME");
	
	lat->DrawLatex(0.65,0.58 ,Form("global Unc. #pm %.1f%%",2.4));
	lat->DrawLatex(0.15,0.92 , "CMS work in progress");

        TLegend* leg = new TLegend(0.6,0.62,0.9,0.85,NULL,"brNDC");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry((TObject*)0, "B^{+}", "C");
        leg->AddEntry((TObject*)0, "2018 PbPb + 2017 pp", "");
        leg->AddEntry((TObject*)0, "Cent. 0-90%", "");
        leg->AddEntry(BPRAAGraph,"|y|<2.4","P");
        leg->AddEntry(BPRAAGraph_low,"|y|>1.5","P");
        leg->AddEntry((TObject*)0, "", "");
        leg->AddEntry((TObject*)0, "2015 PbPb + 2015 pp", "");
        leg->AddEntry((TObject*)0, "Cent. 0-100%", "");
	    leg->AddEntry(BPRAAGraph2015,"|y|<2.4","P");
        leg->SetTextSize(0.025);
        leg->Draw("same");      

	c2->SaveAs("RAAPlots/BPRAACompairson.pdf");
// DRAW B+ RAA vs 2015 RAA

// DRAW Bs RAA vs 2015 RAA
//  Bs RAA
    TGraphAsymmErrors *BsRAAGraph_bs = new TGraphAsymmErrors(NBins_high, RAAX_high, BsRAAY_high,RAAXErrDown_high, RAAXErrUp_high,BsRAAYErrDown_high,BsRAAYErrUp_high);
	TGraphAsymmErrors *BsRAAGraphSyst_bs = new TGraphAsymmErrors(NBins_high, RAAX_high, BsRAAY_high,RAAXErrDown_high, RAAXErrUp_high,BsRAAYSystDown_high,BsRAAYSystUp_high);

	TGraphAsymmErrors *BsRAAGraph_low_bs = new TGraphAsymmErrors(NBins_low, RAAX_low, BsRAAY_low, RAAXErrDown_low, RAAXErrUp_low,BsRAAYErrDown_low,BsRAAYErrUp_low);
	TGraphAsymmErrors *BsRAAGraph_low_just_marker_bs = new TGraphAsymmErrors(NBins_low, RAAX_low, BsRAAY_low, zero1, zero1, zero1, zero1);
	TGraphAsymmErrors *BsRAAGraphSyst_low_bs = new TGraphAsymmErrors(NBins_low, RAAX_low, BsRAAY_low, RAAXErrDown_low, RAAXErrUp_low,BsRAAYSystDown_low,BsRAAYSystUp_low);
//  Bs RAA
//  2015	
    TGraphAsymmErrors *BsRAAGraph2015_bs = new TGraphAsymmErrors(NBins2015, BsRAAX2015, BsRAAY2015,BsRAAXErrDown2015, BsRAAXErrUp2015,BsRAAYErrDown2015,BsRAAYErrUp2015);
	TGraphAsymmErrors *BsRAAGraphSyst2015_bs = new TGraphAsymmErrors(NBins2015, BsRAAX2015, BsRAAY2015,BsRAAXErrDown2015, BsRAAXErrUp2015,BsRAAYSystDown2015,BsRAAYSystUp2015);
//  2015	

	BsRAAGraph_low_just_marker_bs ->SetMarkerStyle(20);
	BsRAAGraph_low_just_marker_bs ->SetMarkerSize(0.9);
	BsRAAGraph_low_just_marker_bs ->SetMarkerColor(kWhite);
	BsRAAGraph_low_bs ->SetMarkerColor(kBlue+2);
	BsRAAGraph_low_bs ->SetLineColor(kBlue+2);
	BsRAAGraph_low_bs ->SetMarkerStyle(24);
	BsRAAGraph_low_bs ->SetMarkerSize(1);
	BsRAAGraphSyst_low_bs ->SetFillColorAlpha(kBlue-3,0.5);
	BsRAAGraphSyst_low_bs ->SetLineColor(kBlue-3);
	
    BsRAAGraph_bs->SetLineColor(kBlue+2);
	BsRAAGraph_bs->SetMarkerStyle(20);
	BsRAAGraph_bs->SetMarkerSize(1);
	BsRAAGraph_bs->SetMarkerColor(kBlue+2);
	BsRAAGraphSyst_bs->SetFillColorAlpha(kBlue-3,0.5);
	BsRAAGraphSyst_bs->SetLineColor(kBlue-3);
	
    BsRAAGraph2015_bs->SetLineColor(kOrange+1);
	BsRAAGraph2015_bs->SetMarkerColor(kOrange+1);
	BsRAAGraph2015_bs->SetMarkerStyle(21);
	BsRAAGraph2015_bs->SetMarkerSize(1);
	BsRAAGraphSyst2015_bs->SetFillColorAlpha(kOrange+1,0.5);
	BsRAAGraphSyst2015_bs->SetLineColor(kOrange+1);

    	HisEmptyRAA_bs_bp->Draw(); //vertical range is ok
        BsRAAGraphSyst_low_bs->Draw("5same");
        BsRAAGraphSyst_bs->Draw("5same");
        BsRAAGraph_bs->Draw("epSAME");
        BsRAAGraph_low_bs->Draw("epSAME");
        BsRAAGraph_low_just_marker_bs->Draw("epSAME");
	    Unity->Draw("SAME");

	    c2->SaveAs("RAAPlots/BsRAA.pdf");

        HisEmptyRAA_bs->Draw();
        BsRAAGraphSyst2015_bs->Draw("5same");
        BsRAAGraph2015_bs->Draw("epSAME");
        BsRAAGraphSyst_low_bs->Draw("5same");
        BsRAAGraphSyst_bs->Draw("5same");
        BsRAAGraph_bs->Draw("epSAME");
        BsRAAGraph_low_bs->Draw("epSAME");
        BsRAAGraph_low_just_marker_bs->Draw("epSAME");
        Unity->Draw("SAME");

        lat->DrawLatex(0.65,0.57 ,Form("global Unc. #pm %.1f%%",2.4));
	lat->DrawLatex(0.15,0.92 , "CMS work in progress");

        TLegend* leg2 = new TLegend(0.6,0.62,0.9,0.85,NULL,"brNDC");
        leg2->SetBorderSize(0);
        leg2->SetFillStyle(0);
        leg2->AddEntry((TObject*)0, "B_{s}^{0}", "C");
        leg2->AddEntry((TObject*)0, "2018 PbPb + 2017 pp", "");
        leg2->AddEntry((TObject*)0, "Cent. 0-90%", "");
        leg2->AddEntry(BsRAAGraph_bs,"|y|<2.4","P");
        leg2->AddEntry(BsRAAGraph_low_bs,"|y|>1.5","P");
        leg2->AddEntry((TObject*)0, "", "");
        leg2->AddEntry((TObject*)0, "2015 PbPb + 2015 pp", "");
        leg2->AddEntry((TObject*)0, "Cent. 0-100%", "");
	    leg2->AddEntry(BsRAAGraph2015_bs,"|y|<2.4","P");
        leg2->SetTextSize(0.025);
        leg2->Draw("same");      


        c2->SaveAs("RAAPlots/BsRAAComparison.pdf");
// DRAW Bs RAA vs 2015 RAA

// Bs0 vs B+ 2017 //

    //distinguish plots if black and white
    BsRAAGraph_low_bs ->SetMarkerStyle(25);
    BsRAAGraph_low_just_marker_bs ->SetMarkerStyle(21);
	BsRAAGraph_bs->SetMarkerStyle(21);
    //distinguish plots if black and white

        HisEmptyRAA_bs_bp->Draw();
	    BPRAAGraphSyst_low->Draw("5same");
	    BPRAAGraphSyst->Draw("5same");
	    BPRAAGraph->Draw("ep");
	    BPRAAGraph_low->Draw("epSAME");
	    BPRAAGraph_low_just_marker->Draw("epSAME");
        BsRAAGraphSyst_low_bs->Draw("5same");
        BsRAAGraphSyst_bs->Draw("5same");
        BsRAAGraph_bs->Draw("epSAME");
        BsRAAGraph_low_bs->Draw("epSAME");
        BsRAAGraph_low_just_marker_bs->Draw("epSAME");
        Unity->Draw("SAME");


        TLegend* leg3 = new TLegend(0.55,0.59,0.9,0.85,NULL,"brNDC");
        leg3->SetBorderSize(0);
        leg3->SetFillStyle(0);
        leg3->AddEntry((TObject*)0, "B^{0}_{s} 2018 PbPb + 2017 pp", "");
        leg3->AddEntry((TObject*)0, "Cent. 0-90%", "");
        leg3->AddEntry(BsRAAGraph_bs,"|y|<2.4","P");
        leg3->AddEntry(BsRAAGraph_low_bs,"|y|>1.5","P");
        leg3->AddEntry((TObject*)0, "", "");
        leg3->AddEntry((TObject*)0, "B^{+} 2018 PbPb + 2017 pp", "");
        leg3->AddEntry((TObject*)0, "Cent. 0-90%", "");
        leg3->AddEntry(BPRAAGraph,"|y|<2.4","P");
        leg3->AddEntry(BPRAAGraph_low,"|y|>1.5","P");        
        leg3->SetTextSize(0.025);
        leg3->Draw("same");      

        lat->DrawLatex(0.6,0.5 ,Form("global Unc. #pm %.1f%%",2.4));
	lat->DrawLatex(0.15,0.92 , "CMS work in progress");
    
        
        c2->SaveAs("RAAPlots/Bs_BP_RAA.pdf");

// Bs0 vs B+ 2017 //


// BDOUBLE RATIOS //  (Bs / B+)
//Compute the Ratios
float DR_RAAY_high[NBins_high];
float DR_RAAY_low[NBins_low] = { BsRAAY_low[0]/BPRAAY_low[0] };
float DR_RAAYErrDown_high[NBins_high];
float DR_RAAYErrUp_high[NBins_high];
float DR_RAAYSystDown_high[NBins_high];
float DR_RAAYSystUp_high[NBins_high];
float DR_RAAYErrDown_low[NBins_low] ={sqrtf( pow(DR_RAAY_low[0],2) * ( pow(BsRAAYErrDown_low[0]/BsRAAY_low[0] ,2) + pow(BPRAAYErrDown_low[0]/BPRAAY_low[0] ,2) ))};
float DR_RAAYErrUp_low[NBins_low]   ={sqrtf( pow(DR_RAAY_low[0],2) * ( pow(BsRAAYErrUp_low[0]/BsRAAY_low[0]   ,2) + pow(BPRAAYErrUp_low[0]/BPRAAY_low[0]   ,2) ))};


float DR_RAAYSystDown_low[NBins_low]={sqrtf( pow(DR_RAAY_low[0],2) * ( pow(BsRAAYSyst_downDR[0]/BsRAAY_low[0],2) + pow(BPRAAYSyst_downDR[0]/BPRAAY_low[0],2) ))};
float DR_RAAYSystUp_low[NBins_low]  ={sqrtf( pow(DR_RAAY_low[0],2) * ( pow(BsRAAYSyst_upDR[0]/BsRAAY_low[0]  ,2) + pow(BPRAAYSyst_upDR[0]/BPRAAY_low[0]  ,2) ))};

for (int i=0; i<NBins_high; i++){
    DR_RAAY_high[i]=BsRAAY_high[i]/BPRAAY_high[i];
    DR_RAAYErrDown_high[i] = sqrtf( pow(DR_RAAY_high[i],2) * ( pow(BsRAAYErrDown_high[i]/BsRAAY_high[i] ,2) + pow(BPRAAYErrDown_high[i]/BPRAAY_high[i] ,2) ));
    DR_RAAYErrUp_high[i]   = sqrtf( pow(DR_RAAY_high[i],2) * ( pow(BsRAAYErrUp_high[i]/BsRAAY_high[i]   ,2) + pow(BPRAAYErrUp_high[i]/BPRAAY_high[i]   ,2) ));
    
    
    DR_RAAYSystDown_high[i]= sqrtf( pow(DR_RAAY_high[i],2) * ( pow(BsRAAYSyst_downDR[i+1]/BsRAAY_high[i],2) + pow(BPRAAYSyst_downDR[i+1]/BPRAAY_high[i],2) ));
    DR_RAAYSystUp_high[i]  = sqrtf( pow(DR_RAAY_high[i],2) * ( pow(BsRAAYSyst_upDR[1+i]/BsRAAY_high[i]  ,2) + pow(BPRAAYSyst_upDR[i+1]/BPRAAY_high[i]  ,2) ));
}
//Compute the Ratios


TGraphAsymmErrors * DOUBLERAAGraph     = new TGraphAsymmErrors(NBins_high, RAAX_high, DR_RAAY_high, RAAXErrDown_high, RAAXErrUp_high, DR_RAAYErrDown_high, DR_RAAYErrUp_high);
TGraphAsymmErrors * DOUBLERAAGraph_low     = new TGraphAsymmErrors(NBins_low, RAAX_low, DR_RAAY_low, RAAXErrDown_low, RAAXErrUp_low,DR_RAAYErrDown_low, DR_RAAYErrUp_low);
TGraphAsymmErrors * DOUBLERAAGraph_low_just_marker = new TGraphAsymmErrors(NBins_low, RAAX_low, DR_RAAY_low, zero1, zero1, zero1, zero1);

TGraphAsymmErrors * DOUBLERAAGraphSyst = new TGraphAsymmErrors(NBins_high, RAAX_high, DR_RAAY_high, RAAXErrDown_high, RAAXErrUp_high, DR_RAAYSystDown_high, DR_RAAYSystUp_high);
TGraphAsymmErrors * DOUBLERAAGraphSyst_low = new TGraphAsymmErrors(NBins_low, RAAX_low, DR_RAAY_low, RAAXErrDown_low, RAAXErrUp_low,DR_RAAYSystDown_low, DR_RAAYSystUp_low);


DOUBLERAAGraph_low_just_marker ->SetMarkerStyle(20);
DOUBLERAAGraph_low_just_marker ->SetMarkerSize(0.9);
DOUBLERAAGraph_low_just_marker ->SetMarkerColor(kWhite);
DOUBLERAAGraph_low ->SetMarkerColor(kRed+2);
DOUBLERAAGraph_low ->SetLineColor(kRed+2);
DOUBLERAAGraph_low ->SetMarkerStyle(24);
DOUBLERAAGraph_low ->SetMarkerSize(1);
DOUBLERAAGraphSyst_low ->SetFillColorAlpha(kRed+1,0.5);
DOUBLERAAGraphSyst_low ->SetLineColor(kRed+1);

DOUBLERAAGraph ->SetMarkerColor(kRed+2);
DOUBLERAAGraph ->SetLineColor(kRed+2);
DOUBLERAAGraph ->SetMarkerStyle(20);
DOUBLERAAGraph ->SetMarkerSize(1);
DOUBLERAAGraphSyst ->SetFillColorAlpha(kRed+1,0.5);
DOUBLERAAGraphSyst ->SetLineColor(kRed+1);

	lat->DrawLatex(0.15,0.92 , "CMS work in progress");

        HisEmptyRAAr->Draw();
        DOUBLERAAGraphSyst_low->Draw("5same");
        DOUBLERAAGraphSyst->Draw("5same");
        DOUBLERAAGraph->Draw("epSAME");
        DOUBLERAAGraph_low->Draw("epSAME");
        DOUBLERAAGraph_low_just_marker->Draw("epSAME");
        Unity->Draw("SAME");

        TLegend* leg4 = new TLegend(0.65,0.77,0.9,0.85,NULL,"brNDC");
        leg4->AddEntry((TObject*)0, "y region", "");
        leg4->SetBorderSize(0);
        leg4->SetFillStyle(0);
        leg4->AddEntry(DOUBLERAAGraph,"|y|<2.4","P");
        leg4->AddEntry(DOUBLERAAGraph_low,"|y|>1.5","P");
	    leg4->SetTextSize(0.025);
        leg4->Draw("same");

        c2->SaveAs("RAAPlots/DOUBLE_RATIOS.pdf");

// BDOUBLE RATIOS //



// BDOUBLE RATIOS vs 2015 //

const int nb15_DR = 2;
float DR15_Y[nb15_DR] ={4.0,1.8};
float DR15_Y_Err[nb15_DR] ={1.8,0.7};
float DR15_Y_SYS[nb15_DR] ={1.3,0.3};
float DR15_X_value[nb15_DR]={11,32.5};
float DR15_X_UP_value[nb15_DR]={4,17.5};
float DR15_X_DOWN_value[nb15_DR]={4,17.5};

TGraphAsymmErrors * DOUBLERAAGraph_comp     = new TGraphAsymmErrors(nb15_DR, DR15_X_value, DR15_Y, DR15_X_DOWN_value, DR15_X_UP_value, DR15_Y_Err, DR15_Y_Err);
TGraphAsymmErrors * DOUBLERAAGraphSyst_comp = new TGraphAsymmErrors(nb15_DR, DR15_X_value, DR15_Y, DR15_X_DOWN_value, DR15_X_UP_value, DR15_Y_SYS, DR15_Y_SYS);

DOUBLERAAGraph_comp ->SetMarkerColor(kMagenta-1);
DOUBLERAAGraph_comp ->SetLineColor(kMagenta-1);
DOUBLERAAGraph_comp ->SetMarkerStyle(33);
DOUBLERAAGraph_comp ->SetMarkerSize(1.5);
DOUBLERAAGraphSyst_comp ->SetFillColorAlpha(kMagenta-8,0.5);
DOUBLERAAGraphSyst_comp ->SetLineColor(kMagenta-8);

        HisEmptyRAA_Dr->Draw();
        DOUBLERAAGraphSyst_comp->Draw("5same");
        DOUBLERAAGraph_comp->Draw("epSAME");
        DOUBLERAAGraphSyst_low->Draw("5same");
        DOUBLERAAGraphSyst->Draw("5same");
        DOUBLERAAGraph->Draw("epSAME");
        DOUBLERAAGraph_low->Draw("epSAME");
        DOUBLERAAGraph_low_just_marker->Draw("epSAME");
        Unity->Draw("SAME");

        	lat->DrawLatex(0.15,0.92 , "CMS work in progress");

        TLegend* leg5 = new TLegend(0.6,0.65,0.9,0.85,NULL,"brNDC");
        leg5->SetBorderSize(0);
        leg5->SetFillStyle(0);
        leg5->AddEntry((TObject*)0, "2018 PbPb + 2017 pp", "");
        leg5->AddEntry((TObject*)0, "Cent. 0-90%", "");
        leg5->AddEntry(DOUBLERAAGraph,"|y|<2.4","P");
        leg5->AddEntry(DOUBLERAAGraph_low,"|y|>1.5","P");
        leg5->AddEntry((TObject*)0, "", "");
        leg5->AddEntry((TObject*)0, "2015 PbPb + 2015 pp", "");
        leg5->AddEntry((TObject*)0, "Cent. 0-100%", "");
	    leg5->AddEntry(DOUBLERAAGraph_comp,"|y|<2.4","P");
        leg5->SetTextSize(0.025);
        leg5->Draw("same");        
        
        c2->SaveAs("RAAPlots/DOUBLE_R_COMP_17_15.pdf");



// BDOUBLE RATIOS vs 2015 //











	TFile * fout = new TFile("OutFile/BPRAA.root","RECREATE");
	fout->cd();
	BPRAAGraph->Write();
	fout->Close();
    TFile * fout_bs = new TFile("OutFile/BsRAA.root","RECREATE");
	fout_bs->cd();
	BsRAAGraph_bs->Write();
	fout_bs->Close();

//B+ save things
{
  double lumiUncertainty = TMath::Sqrt(TMath::Power(0.019, 2) + TMath::Power(0.015, 2));
  std::vector<float> globUncert(NBins, lumiUncertainty);
  // summary of errors (in ratio, not percent)

  gSystem->mkdir("../../MakeFinalPlots/NominalPlots/RAA/dataSource/" ,true );
  string outFile = "../../MakeFinalPlots/NominalPlots/RAA/dataSource/RAA_pt_Bp_New.txt";
  ofstream out;
  out.open(outFile);
  unsigned columnWidth = 14;
  out << std::left << std::setw(columnWidth) <<
    "ptmin" << std::setw(columnWidth) << "ptmax" << std::setw(columnWidth) <<
    "central_val" << std::setw(columnWidth) <<
    "statUp" << std::setw(columnWidth) << "statDown" << std::setw(columnWidth) <<
    "systUp" << std::setw(columnWidth) << "systDown" << std::setw(columnWidth) <<
    "glbUp" << std::setw(columnWidth) << "glbDown" << std::setw(columnWidth) <<
    "abscissae" << endl;
  for (auto i = 0; i < NBins; ++i ) {
    out << std::setw(columnWidth) <<
      ptbinsvec[i] << std::setw(columnWidth) <<
      ptbinsvec[i + 1] << std::setw(columnWidth) <<
      setprecision(3) << BPRAAY[i] << std::setw(columnWidth) <<
      setprecision(2) << std::defaultfloat <<
      BPRAAYErrUpRatio[i] << std::setw(columnWidth) <<
      BPRAAYErrDownRatio[i] << std::setw(columnWidth) <<
      BPRAAYSystUpRatio[i] << std::setw(columnWidth) <<
      BPRAAYSystDownRatio[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      setprecision(3) << abscissae[i] << "\n";
  }
  out.close();
}
//B+ save things

//Bs save things
{  std::vector<float> globUncert(NBins, 0.077);
  // summary of errors (in ratio, not percent)

  gSystem->mkdir("../../MakeFinalPlots/NominalPlots/RAA/dataSource/" ,true );
  string outFile = "../../MakeFinalPlots/NominalPlots/RAA/dataSource/RAA_pt_Bs_New.txt";
  ofstream out;
  out.open(outFile);
  unsigned columnWidth = 14;
  out << std::left << std::setw(columnWidth) <<
    "ptmin" << std::setw(columnWidth) << "ptmax" << std::setw(columnWidth) <<
    "central_val" << std::setw(columnWidth) <<
    "statUp" << std::setw(columnWidth) << "statDown" << std::setw(columnWidth) <<
    "systUp" << std::setw(columnWidth) << "systDown" << std::setw(columnWidth) <<
    "glbUp" << std::setw(columnWidth) << "glbDown" << std::setw(columnWidth) <<
    "abscissae" << endl;
  for (auto i = 0; i < NBins; ++i ) {
    out << std::setw(columnWidth) <<
      ptbinsvec[i] << std::setw(columnWidth) << 
      ptbinsvec[i + 1] << std::setw(columnWidth) <<
      setprecision(3) << BsRAAY[i] << std::setw(columnWidth) <<
      setprecision(2) << std::defaultfloat <<
      BsRAAYErrUpRatio[i] << std::setw(columnWidth) <<
      BsRAAYErrDownRatio[i] << std::setw(columnWidth) <<
      BsRAAYSystUpRatio[i] << std::setw(columnWidth) <<
      BsRAAYSystDownRatio[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      globUncert[i] << std::setw(columnWidth) <<
      setprecision(3) << abscissae[i] << "\n";
  }
  out.close();
}
//Bs save things







}