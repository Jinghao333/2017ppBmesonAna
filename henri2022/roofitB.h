#include "TAxis.h"
#include "uti.h"

#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooChi2Var.h"
#include "RooHist.h"
#include "RooProdPdf.h"
#include "RooAddition.h"
#include "RooProduct.h"
#include <RooBifurGauss.h>
#include <RooCmdArg.h>
#include <fstream>
#include <string>
#include <iomanip>
#include "RooMCStudy.h"
#include <RooMinuit.h>



// draw legend and suppress parameters
const bool drawLegend = false;
using namespace RooFit;
using namespace std;

#define BS_MASS 5.36682
#define BP_MASS 5.27915

bool drawSup = 0;
float bkgd;
void clean0 (TH1D* h);
double _ErrCor=1;
double setparam0=100.;
double setparam2=0.02;
double setparam3=0.06;
double fixparam1=BS_MASS;
double Significance;

double real_significance;
double minhisto=5.;
double maxhisto=6.;
int nbinsmasshisto=50;

TString seldata;
TString selmc;
TString selmcgen;
Int_t _count=0;
RooWorkspace* inputw = new RooWorkspace();
RooWorkspace* outputw = new RooWorkspace("w");
RooWorkspace* w_val= new RooWorkspace("w_vl");


//This first function is only for obtaining parameters related to Monte Carlo, for the actual fit, please use 
//the one below this one.
RooFitResult *fitMC(TString variation, TString pdf,TString tree, TCanvas* c, TCanvas* cMC, RooDataSet* ds, RooDataSet* dsMC, RooDataHist* dh, RooDataHist* dhMC, RooRealVar* mass, RooPlot* &outframe, Double_t ptmin, Double_t ptmax, int isMC, TString npfit)
{
	//extern TString pdf_ext = pdf;	

	cout<<"total data: "<<ds->numEntries()<<endl;
	TH1* h = dh->createHistogram("Bmass");
	h->Sumw2(kFALSE);
	h->SetBinErrorOption(TH1::kPoisson);
	h->SetMarkerSize(1);
	h->SetMarkerStyle(20);
	h->SetLineColor(1);
	h->SetLineWidth(2);
	RooPlot* frameMC = mass->frame();
	frameMC->SetTitle("");
	if(tree=="ntKp")frameMC->SetXTitle("m_{J/#psiK^{#pm}} (GeV/c^{2})");
	if(tree=="ntphi")frameMC->SetXTitle("m_{J/#psiK^{+}K^{-}} (GeV/c^{2})");

	frameMC->SetYTitle("Events / (20 MeV/c^{2})");
	frameMC->GetXaxis()->CenterTitle();
	frameMC->GetXaxis()->SetTitleOffset(1.0);
	frameMC->GetYaxis()->SetTitleOffset(2.);
	frameMC->GetXaxis()->SetTitleSize(0.055);
	frameMC->GetYaxis()->SetTitleSize(0.035);
	frameMC->GetXaxis()->SetTitleFont(42);
	frameMC->GetYaxis()->SetTitleFont(42);
	frameMC->GetXaxis()->SetLabelFont(42);
	frameMC->GetYaxis()->SetLabelFont(42);
	frameMC->GetXaxis()->SetLabelSize(0.055);
	frameMC->GetYaxis()->SetLabelSize(0.035);
	frameMC->SetStats(0);
	frameMC->GetXaxis()->SetNdivisions(-50205);

	cMC->cd();
	double init_mean;
	if(tree=="ntphi") init_mean = BS_MASS;
	if(tree=="ntKp") init_mean = BP_MASS;

	RooRealVar meanMC(Form("meanMC%d_%s",_count,pdf.Data()),"",init_mean,5.2,5.4) ;
	RooRealVar sigma1MC(Form("sigma1MC%d",_count),"",0.02,0.01,0.1) ;
	RooRealVar sigma2MC(Form("sigma2MC%d",_count),"",0.055,0.01,0.1) ;
	RooRealVar sigma3MC(Form("sigma3MC%d_%s",_count, pdf.Data()),"",0.0266,0.01,0.1) ;
	RooRealVar sigma4cbMC(Form("sigma4cbMC%d_%s",_count, pdf.Data()),"",0.0266,0.01,0.1) ;
	//RooRealVar sigma5cbMC(Form("sigma5cbMC%d_%s",_count, pdf.Data()),"",0.0266,0.01,0.1) ;
	RooRealVar alphaMC(Form("alphaMC%d_%s",_count,pdf.Data()),"",5.,0,50);
	//RooRealVar alphaMC1(Form("alphaMC1%d_%s",_count,pdf.Data()),"",5.,0,50);
	RooRealVar nMC(Form("nMC_%d_%s", _count, pdf.Data()),"",100,0,500);
	//RooRealVar nMC1(Form("nMC1_%d_%s", _count, pdf.Data()),"",100,0,500);

	RooRealVar* scale;
	scale = new RooRealVar("scale","scale",1,0,2);

	RooProduct scaled_sigma1MC("scaled_sigma1MC","scaled_sigma1MC", RooArgList(*scale,sigma1MC));
	RooProduct scaled_sigma2MC("scaled_sigma2MC","scaled_sigma2MC", RooArgList(*scale,sigma2MC));
	RooProduct scaled_sigma3MC("scaled_sigma3MC","scaled_sigma3MC", RooArgList(*scale,sigma3MC));
	RooProduct scaled_sigma4cbMC("scaled_sigma4cbMC","scaled_sigma4cbMC", RooArgList(*scale,sigma4cbMC));
	//RooProduct scaled_sigma5cbMC("scaled_sigma5cbMC","scaled_sigma5cbMC", RooArgList(*scale,sigma5cbMC));

	RooGaussian sig1MC(Form("sig1MC%d_%s",_count,pdf.Data()),"",*mass,meanMC,scaled_sigma1MC);  
	RooGaussian sig2MC(Form("sig2MC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma2MC);  
	RooGaussian sig3MC(Form("sig3MC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma3MC);  
	RooCBShape  CBMC(Form("CBMC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma4cbMC, alphaMC, nMC);
	//RooCBShape  CBMC1(Form("CBMC1%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma5cbMC, alphaMC1, nMC1);


	RooRealVar sig1fracMC(Form("sig1fracMC%d",_count),"",0.5,0.,1.);
	RooRealVar sig2fracMC(Form("sig2fracMC%d_%s",_count, pdf.Data()),"",0.,0.,1.);
	//RooRealVar sig3fracMC(Form("sig3fracMC%d_%s",_count, pdf.Data()),"",0.5,0.,1.);

	RooAddPdf* sigMC;
	RooRealVar nsigMC(Form("nsigMC%d",_count),"",1, 0, 1.2 * dsMC->sumEntries());

	//if((variation=="signal" && pdf=="2cb")) sigMC = new RooAddPdf(Form("sigMC%d_%s",_count,pdf.Data()),"",RooArgList(CBMC,CBMC1),sig1fracMC);
	if((variation=="" && pdf=="") || variation== "background" || (variation=="signal" && pdf=="fixed" )) sigMC = new RooAddPdf(Form("sigMC%d_%s",_count,pdf.Data()),"",RooArgList(sig1MC,sig2MC),sig1fracMC);
	if(variation=="signal" && pdf=="3gauss") sigMC = new RooAddPdf(Form("sigMC%d_%s",_count, pdf.Data()), "", RooArgList(sig1MC, sig2MC, sig3MC), RooArgList(sig1fracMC, sig2fracMC));
	if(variation=="signal" && pdf=="gauss_cb") sigMC = new RooAddPdf(Form("sigMC%d_%s",_count, pdf.Data()), "", RooArgList(sig1MC, CBMC), sig1fracMC);

	RooAddPdf* modelMC;
	//if(variation =="signal" && pdf=="1gauss") modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(sig1MC),RooArgList(nsigMC));

	if((variation=="signal" && (pdf=="gauss_cb"|| pdf=="3gauss"|| pdf=="fixed"))||variation=="background") modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(*sigMC),RooArgList(nsigMC)); //PUT IN CB
	if(variation =="" && pdf=="") modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(*sigMC),RooArgList(nsigMC));

//////////ROOFIT ROOFIT ROOFIT  MC MC MC MC 

	double SignalWidth = 0.2;
	mass->setRange("signal",init_mean-SignalWidth, init_mean+SignalWidth);    //focous the MC fit to the signal region to prevent statistical flutuations
	std::cout<<"sum Entries: "<<dsMC->sumEntries()<<std::endl;
	
	RooFitResult* fitResultMC;
	scale->setConstant();
	fitResultMC = modelMC->fitTo(*dsMC,Save(), Range("signal"));
	scale->setConstant(false);
	std::cout << "sigma1MC is " << sigma1MC.getVal() << std::endl;

	return fitResultMC;
}

//The function above is only for the extraction of parameters related to Monte Carlo, please use the below one for the actual fit.

RooFitResult *fit(TString variation, TString pdf,TString tree, TCanvas* c, TCanvas* cMC, RooDataSet* ds, RooDataSet* dsMC, RooDataHist* dh, RooDataHist* dhMC, RooRealVar* mass, RooPlot* &outframe, Double_t ptmin, Double_t ptmax, int isMC, TString npfit)
{
	//extern TString pdf_ext = pdf;	

	cout<<"total data: "<<ds->numEntries()<<endl;
	TH1* h = dh->createHistogram("Bmass");
	h->Sumw2(kFALSE);
	h->SetBinErrorOption(TH1::kPoisson);
	h->SetMarkerSize(1);
	h->SetMarkerStyle(20);
	h->SetLineColor(1);
	h->SetLineWidth(2);
	RooPlot* frameMC = mass->frame();
	frameMC->SetTitle("");
	if(tree=="ntKp")frameMC->SetXTitle("m_{J/#psiK^{#pm}} (GeV/c^{2})");
	if(tree=="ntphi")frameMC->SetXTitle("m_{J/#psiK^{+}K^{-}} (GeV/c^{2})");

	frameMC->SetYTitle("Events / (20 MeV/c^{2})");
	frameMC->GetXaxis()->CenterTitle();
	frameMC->GetXaxis()->SetTitleOffset(1.0);
	frameMC->GetYaxis()->SetTitleOffset(2.);
	frameMC->GetXaxis()->SetTitleSize(0.055);
	frameMC->GetYaxis()->SetTitleSize(0.035);
	frameMC->GetXaxis()->SetTitleFont(42);
	frameMC->GetYaxis()->SetTitleFont(42);
	frameMC->GetXaxis()->SetLabelFont(42);
	frameMC->GetYaxis()->SetLabelFont(42);
	frameMC->GetXaxis()->SetLabelSize(0.055);
	frameMC->GetYaxis()->SetLabelSize(0.035);
	frameMC->SetStats(0);
	frameMC->GetXaxis()->SetNdivisions(-50205);

	cMC->cd();
	double init_mean;
	if(tree=="ntphi") init_mean = BS_MASS;
	if(tree=="ntKp") init_mean = BP_MASS;

	RooRealVar meanMC(Form("meanMC%d_%s",_count,pdf.Data()),"",init_mean,5.2,5.4) ;
	RooRealVar sigma1MC(Form("sigma1MC%d",_count),"",0.02,0.01,0.1) ;
	RooRealVar sigma2MC(Form("sigma2MC%d",_count),"",0.055,0.01,0.1) ;
	RooRealVar sigma3MC(Form("sigma3MC%d_%s",_count, pdf.Data()),"",0.0266,0.01,0.1) ;
	RooRealVar sigma4cbMC(Form("sigma4cbMC%d_%s",_count, pdf.Data()),"",0.0266,0.01,0.1) ;
	//RooRealVar sigma5cbMC(Form("sigma5cbMC%d_%s",_count, pdf.Data()),"",0.0266,0.01,0.1) ;
	RooRealVar alphaMC(Form("alphaMC%d_%s",_count,pdf.Data()),"",5.,0,50);
	//RooRealVar alphaMC1(Form("alphaMC1%d_%s",_count,pdf.Data()),"",5.,0,50);
	RooRealVar nMC(Form("nMC_%d_%s", _count, pdf.Data()),"",100,0,500);
	//RooRealVar nMC1(Form("nMC1_%d_%s", _count, pdf.Data()),"",100,0,500);

	RooRealVar* scale;
	scale = new RooRealVar("scale","scale",1,0,2);

	RooProduct scaled_sigma1MC("scaled_sigma1MC","scaled_sigma1MC", RooArgList(*scale,sigma1MC));
	RooProduct scaled_sigma2MC("scaled_sigma2MC","scaled_sigma2MC", RooArgList(*scale,sigma2MC));
	RooProduct scaled_sigma3MC("scaled_sigma3MC","scaled_sigma3MC", RooArgList(*scale,sigma3MC));
	RooProduct scaled_sigma4cbMC("scaled_sigma4cbMC","scaled_sigma4cbMC", RooArgList(*scale,sigma4cbMC));
	//RooProduct scaled_sigma5cbMC("scaled_sigma5cbMC","scaled_sigma5cbMC", RooArgList(*scale,sigma5cbMC));

	RooGaussian sig1MC(Form("sig1MC%d_%s",_count,pdf.Data()),"",*mass,meanMC,scaled_sigma1MC);  
	RooGaussian sig2MC(Form("sig2MC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma2MC);  
	RooGaussian sig3MC(Form("sig3MC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma3MC);  
	RooCBShape  CBMC(Form("CBMC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma4cbMC, alphaMC, nMC);
	//RooCBShape  CBMC1(Form("CBMC1%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma5cbMC, alphaMC1, nMC1);


	RooRealVar sig1fracMC(Form("sig1fracMC%d_%s",_count, pdf.Data()),"",0.5,0.,1.);
	RooRealVar sig2fracMC(Form("sig2fracMC%d_%s",_count, pdf.Data()),"",0.,0.,1.);
	//RooRealVar sig3fracMC(Form("sig3fracMC%d_%s",_count, pdf.Data()),"",0.5,0.,1.);

	RooAddPdf* sigMC;
	RooRealVar nsigMC(Form("nsigMC%d",_count),"",1, 0, 1.2 * dsMC->sumEntries());

	//if((variation=="signal" && pdf=="2cb")) sigMC = new RooAddPdf(Form("sigMC%d_%s",_count,pdf.Data()),"",RooArgList(CBMC,CBMC1),sig1fracMC);
	if((variation=="" && pdf=="") || variation== "background" || (variation=="signal" && pdf=="fixed" )) sigMC = new RooAddPdf(Form("sigMC%d_%s",_count,pdf.Data()),"",RooArgList(sig1MC,sig2MC),sig1fracMC);
	if(variation=="signal" && pdf=="3gauss") sigMC = new RooAddPdf(Form("sigMC%d_%s",_count, pdf.Data()), "", RooArgList(sig1MC, sig2MC, sig3MC), RooArgList(sig1fracMC, sig2fracMC));
	if(variation=="signal" && pdf=="gauss_cb") sigMC = new RooAddPdf(Form("sigMC%d_%s",_count, pdf.Data()), "", RooArgList(sig1MC, CBMC), sig1fracMC);

	RooAddPdf* modelMC;
	//if(variation =="signal" && pdf=="1gauss") modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(sig1MC),RooArgList(nsigMC));

	if((variation=="signal" && (pdf=="gauss_cb"|| pdf=="3gauss"|| pdf=="fixed"))||variation=="background") modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(*sigMC),RooArgList(nsigMC));
	if(variation =="" && pdf=="") modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(*sigMC),RooArgList(nsigMC));

//////////ROOFIT ROOFIT ROOFIT  MC MC MC MC 

	double SignalWidth = 0.2;
	mass->setRange("signal",init_mean-SignalWidth, init_mean+SignalWidth);    //focous the MC fit to the signal region to prevent statistical flutuations
	std::cout<<"sum Entries: "<<dsMC->sumEntries()<<std::endl;
	
	RooFitResult * fitResultMC;
	scale->setConstant();
	fitResultMC = modelMC->fitTo(*dsMC,Save(), Range("signal"));
	scale->setConstant(false);
	std::cout << "sigma1MC is " << sigma1MC.getVal() << std::endl;
///////////ROOFIT ROOFIT ROOFIT MC MC MC MC

/*	std::cout<<"mean_MC= "<<meanMC.getVal()<<std::endl;
	std::cout<<"sigma1_MC= "<<sigma1MC.getVal()<<std::endl;
	std::cout<<"sigma2_MC= "<<sigma2MC.getVal()<<std::endl;
	std::cout<<"fraction1_MC= "<<sig1fracMC.getVal()<<std::endl;
	std::cout<<"fraction2_MC= "<<sig2fracMC.getVal()<<std::endl;
*/

	dsMC->plotOn(frameMC,Name(Form("dsMC_cut%d",_count)),Binning(nbinsmasshisto),MarkerSize(1.55),MarkerStyle(20),LineColor(1),LineWidth(4));
	if(pdf!="1gauss"){
		modelMC->plotOn(frameMC,Name(Form("sigMC%d_%s",_count, pdf.Data())),Components(*sigMC),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("L"),FillStyle(3002),FillColor(kOrange-3),LineStyle(7),LineColor(kOrange-3),LineWidth(4));
		modelMC->plotOn(frameMC,Name(Form("sigFMC%d_%s",_count, pdf.Data())),Components(*sigMC),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("F"),FillStyle(3002),FillColor(kOrange-3),LineStyle(7),LineColor(kOrange-3),LineWidth(4));
	} else {
		modelMC->plotOn(frameMC,Name(Form("sigMC%d_%s",_count, pdf.Data())),Components(sig1MC),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("L"),FillStyle(3002),FillColor(kOrange-3),LineStyle(7),LineColor(kOrange-3),LineWidth(4));
		modelMC->plotOn(frameMC,Name(Form("sigFMC%d_%s",_count, pdf.Data())),Components(sig1MC),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("F"),FillStyle(3002),FillColor(kOrange-3),LineStyle(7),LineColor(kOrange-3),LineWidth(4));
	}

	modelMC->plotOn(frameMC,Name(Form("modelMC%d_%s",_count, pdf.Data())),Normalization(1.0,RooAbsReal::RelativeExpected),Precision(1e-6),DrawOption("L"),LineColor(2),LineWidth(4));
	double x_1 = 0.58;
	double x_2 = 0.95;
	double y_2 = 0.92;
	double y_1;
	double y_space = 0.04;
	int nitems = 7;

	y_1 = y_2 - y_space*nitems;

	modelMC->paramOn(frameMC,Layout(x_1, x_2, y_1), Format("NEU",AutoPrecision(1)));
	frameMC->getAttText()->SetTextSize(0.02);
	frameMC->SetMaximum(nsigMC.getVal()*1.2);
	frameMC->Draw();
	cMC->RedrawAxis();
	//  cMC->SetLogy();
	cMC->SaveAs("mc_afterfit_test.pdf");
	//	RooPlot* frameMC = mass->frame();
	c->cd();

	RooPlot* frame = mass->frame();
	TPad *p1 = new TPad("p1","p1",0.,0.28,1.,0.99);
	// TPad *p1 = new TPad("p1","p1",0.05,0.05,0.99,0.99);
	p1->SetBorderMode(1); 
	p1->SetFrameBorderMode(0); 
	p1->SetBorderSize(2);
	p1->SetBottomMargin(0.10);
	p1->Draw(); 

	TPad *p2 = new TPad("p2","p2",0.,0.10,1.,0.32);// 0.26 
	p2->SetTopMargin(0.);    
	p2->SetBorderMode(0);
	p2->SetBorderSize(2); 
	p2->SetFrameBorderMode(0); 
	p2->SetTicks(1,1); 
	//  p2->SetBottomMargin(0.2);
	p2->Draw();


	double mass_peak = 0;
	if(tree=="ntphi") mass_peak = BS_MASS;
	if(tree=="ntKp") mass_peak = BP_MASS;

	double n_signal_initial = ds->sumEntries(TString::Format("abs(Bmass-%g)<0.05",mass_peak)) - ds->sumEntries(TString::Format("abs(Bmass-%g)<0.10&&abs(Bmass-%g)>0.05",mass_peak,mass_peak));
	if(n_signal_initial<0){n_signal_initial=1;}

	double n_combinatorial_initial = ds->sumEntries() - n_signal_initial;
	p1->cd();

	RooRealVar mean(Form("mean%d",_count),"",meanMC.getVal(),5.,6.) ;
	RooRealVar sigma1(Form("sigma1%d",_count),"",sigma1MC.getVal(),0.01,0.1) ;
	RooRealVar sigma2(Form("sigma2%d",_count),"",sigma2MC.getVal(),0.01,0.1) ;
	RooRealVar sigma3(Form("sigma3%d",_count),"",sigma3MC.getVal(),0.01,0.1) ;
	RooRealVar sigma4cb(Form("sigma4cb%d",_count),"",sigma4cbMC.getVal(),0.01,0.1) ;
	//RooRealVar sigma5cb(Form("sigma5cb%d",_count),"",sigma5cbMC.getVal(),0.01,0.1) ;
	RooRealVar alpha(Form("alpha%d_%s",_count,pdf.Data()),"",alphaMC.getVal(),0,50);
	//RooRealVar alpha1(Form("alpha1%d_%s",_count,pdf.Data()),"",alphaMC1.getVal(),0,50);
	RooRealVar n(Form("n_%d_%s", _count, pdf.Data()),"",nMC.getVal(),0,500);
	//RooRealVar n1(Form("n1_%d_%s", _count, pdf.Data()),"",nMC1.getVal(),0,500);
	
	RooProduct scaled_sigma1("scaled_sigma1","scaled_sigma1", RooArgList(*scale,sigma1));
	RooProduct scaled_sigma2("scaled_sigma2","scaled_sigma2", RooArgList(*scale,sigma2));
	RooProduct scaled_sigma3("scaled_sigma3","scaled_sigma3", RooArgList(*scale,sigma3));
	RooProduct scaled_sigma4cb("scaled_sigma4cb","scaled_sigma4cb", RooArgList(*scale,sigma4cb));
	//RooProduct scaled_sigma5cb("scaled_sigma5cb","scaled_sigma5cb", RooArgList(*scale,sigma5cb));
	RooGaussian sig1(Form("sig1%d",_count),"",*mass,mean,scaled_sigma1);  
	RooGaussian sig2(Form("sig2%d",_count),"",*mass,mean,scaled_sigma2);  
	RooGaussian sig3(Form("sig3%d",_count),"",*mass,mean,scaled_sigma3);  
	RooCBShape  CB(Form("CB%d_%s",_count, pdf.Data()),"",*mass,mean,scaled_sigma4cb, alpha, n);
	//RooCBShape  CB1(Form("CB1%d_%s",_count, pdf.Data()),"",*mass,mean,scaled_sigma5cb, alpha1, n1);
	RooRealVar c1(Form("c1%d",_count),"",1.,0.,5.);

	RooRealVar sig1frac(Form("sig1frac%d",_count),"",sig1fracMC.getVal(),0.,1.);
	RooRealVar sig2frac(Form("sig2frac%d",_count),"",sig2fracMC.getVal(),0.,1.);
	RooAddPdf* sig;
	if(variation=="signal" && pdf=="3gauss") sig = new RooAddPdf(Form("sig%d",_count), "", RooArgList(sig1, sig2, sig3), RooArgList(sig1frac, sig2frac));
	if(variation=="signal" && pdf=="gauss_cb") sig = new RooAddPdf(Form("sig%d",_count),"",RooArgList(sig1, CB), sig1frac);
	//if((variation=="signal" && pdf=="2cb") || variation== "background" || (variation=="signal" && pdf=="fixed")) sig = new RooAddPdf(Form("sig%d",_count),"",RooArgList(CB,CB1),sig1frac);
	if((variation=="" && pdf=="") || variation== "background" || (variation=="signal" && pdf=="fixed")) sig = new RooAddPdf(Form("sig%d",_count),"",RooArgList(sig1,sig2),sig1frac);
	
///////////////// BACKGROUND FUNCTIONS

	//RooChebychev bkg(Form("bkg%d",_count),"",*mass,RooArgSet(a0,a1,a2));
	RooRealVar a0(Form("a0%d",_count),"",1,-5,5);
	RooRealVar a1(Form("a1%d",_count),"",1,-5,5);
	RooRealVar a2(Form("a2%d",_count),"",1,-5e3,5e3);
	// RooRealVar a3(Form("a3%d",_count),"",1e0,-1e4,1e4);
	RooPolynomial bkg_1st(Form("bkg%d",_count),"",*mass,RooArgSet(a0));
	RooPolynomial bkg_2nd(Form("bkg%d",_count),"",*mass,RooArgSet(a0,a1));
	RooPolynomial bkg_3rd(Form("bkg%d",_count),"",*mass,RooArgSet(a0,a1,a2));
	RooRealVar lambda(Form("lambda%d", _count), "lambda",-0.5, -3., 1.);
	RooExponential bkg(Form("bkg%d",_count),"",*mass,lambda);
	RooGenericPdf peakbg(Form("peakbg%d",_count),"",Form("(%s)",npfit.Data()),RooArgSet(*mass));

///////////////// BACKGROUND FUNCTIONS

	RooRealVar nsig(Form("nsig%d",_count),"",n_signal_initial,-1,ds->sumEntries()*3);
	RooRealVar nbkg(Form("nbkg%d",_count),"",n_combinatorial_initial,0.,ds->sumEntries());
	RooRealVar npeakbg(Form("npeakbg%d",_count),"",1,0,1e5);
	RooAddPdf* model;

/////////////////Bs Bs Bs Bs Bs Bs Bs Bs

	if((variation=="" && pdf=="") || (pdf=="mass_range")) model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg),RooArgList(nsig,nbkg));
	if(variation=="background" && pdf=="1st") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg_1st),RooArgList(nsig,nbkg));
	if(variation=="background" && pdf=="2nd") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg_2nd),RooArgList(nsig,nbkg));
	if(variation=="background" && pdf=="3rd") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg_3rd),RooArgList(nsig,nbkg));
	if(variation=="signal" && pdf=="1gauss") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(sig1,bkg),RooArgList(nsig,nbkg));
	if(variation=="signal" && (pdf=="3gauss"|| pdf=="fixed"|| pdf=="gauss_cb")) model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig, bkg),RooArgList(nsig, nbkg));

/////////////////Bs Bs Bs Bs Bs Bs Bs Bs


/////////////////BP BP BP BP BP BP BP BP

        if(npfit != "1" && variation=="" && pdf=="") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(bkg,*sig,peakbg),RooArgList(nbkg,nsig,npeakbg));
	if(npfit != "1" && pdf=="mass_range"){ model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig,bkg),RooArgList(nsig,nbkg));}
        if(npfit != "1" && variation=="background" && pdf=="1st") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(bkg_1st,*sig,peakbg),RooArgList(nbkg,nsig,npeakbg));
	if(npfit != "1" && variation=="background" && pdf=="2nd") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(bkg_2nd,*sig,peakbg),RooArgList(nbkg,nsig,npeakbg));
	if(npfit != "1" && variation=="background" && pdf=="3rd") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(bkg_3rd,*sig,peakbg),RooArgList(nbkg,nsig,npeakbg));
	if(npfit != "1" && variation=="signal" && pdf=="1gauss") model = new RooAddPdf(Form("model%d",_count),"",RooArgList(bkg,sig1,peakbg),RooArgList(nbkg,nsig,npeakbg));
	if(npfit != "1" && (variation=="signal" && (pdf=="3gauss"|| pdf=="fixed"|| pdf=="gauss_cb"))) model = new RooAddPdf(Form("model%d",_count),"",RooArgList(*sig, bkg, peakbg),RooArgList(nsig, nbkg, npeakbg));

/////////////////BP BP BP BP BP BP BP BP
	sigma1.setConstant();
	if(pdf!="1gauss"){
		sigma2.setConstant();
		sig1frac.setConstant();
	}
	if(variation=="signal" && pdf=="3gauss"){
		sigma3.setConstant();
		sig2frac.setConstant();  
	}
	if(variation=="signal" && pdf=="gauss_cb"){
		sigma4cb.setConstant();
		n.setConstant();
		alpha.setConstant();
	}
	if(variation=="signal" && pdf=="fixed") mean.setConstant();
/*	if(variation=="signal" && pdf=="2cb"){
		sigma5cb.setConstant();
		n1.setConstant();
		alpha1.setConstant();
		sigma4cb.setConstant();
		n.setConstant();
		alpha.setConstant();

	}
*/
////// ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT

/*RooFitResult* fitResult;
        if(pdf =="mass_range"){
                mass->setRange("m_range", 5.19 , 6 );    //set a range to be used if pdf = mass_range
                fitResult = model->fitTo(*ds,Save(), Minos() , Extended(kTRUE),Range("m_range"));}
        else{
                mass->setRange("m_range", minhisto , maxhisto );
                fitResult = model->fitTo(*ds,Save(), Minos(),  RooFit::PrintLevel(0) , Extended(kTRUE));}
*/

mass->setRange("m_range", 5.19 , 6 );    //set a range to be used if pdf = mass_range
mass->setRange("all", minhisto, maxhisto);    //set a range to be used if pdf = mass_range
TString fitRange = (pdf == "mass_range")? "m_range" : "all";
RooFitResult* fitResult = model->fitTo(*ds,Save(), Minos(),Extended(kTRUE), Range(fitRange));

///// ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT

	w_val->import(*model);
	w_val->import(nsig);

	cout << "nbkg->getVal() = " << nbkg.getVal() << endl;
	cout << "npeakbg->getVal() = " << npeakbg.getVal() << endl;

	RooAbsReal* nll = model->createNLL(*ds);
	double log_likelihood= nll->getVal();

/*ds->plotOn(frame,Name(Form("ds_cut%d",_count)),Binning(nbinsmasshisto),MarkerSize(1),MarkerStyle(20),MarkerColor(1),LineColor(1),LineWidth(2),LineColor(1));
        if(npfit != "1"){
                model->plotOn(frame, Name(Form("peakbg%d",_count)) , Components(peakbg), RooFit::NormRange("m_range"), Precision(1e-6),DrawOption("L"),FillStyle(3005),FillColor(kGreen+4),LineStyle(1),LineColor(kGreen+4),LineWidth(3));
                model->plotOn(frame, Name(Form("peakbgF%d",_count)),Components(peakbg), RooFit::NormRange("m_range"), Precision(1e-6),DrawOption("F"),FillStyle(3002),FillColor(kGreen+4),LineStyle(7),LineColor(kGreen+4),LineWidth(3));
                        }
        model->plotOn(frame,Name(Form("bkg%d",_count)) , Components(bkg), RooFit::NormRange("m_range"),Precision(1e-6),DrawOption("L"),LineStyle(7),LineColor(4),LineWidth(3));
        model->plotOn(frame,Name(Form("model%d",_count)), RooFit:: NormRange("m_range"), Precision(1e-6),DrawOption("L"),LineColor(2),LineWidth(3));

        if(pdf!="1gauss"){
                model->plotOn(frame,Name(Form("sig%d",_count)), Components(*sig),RooFit::NormRange("m_range"), Precision(1e-6),DrawOption("L"),FillStyle(3002),FillColor(kOrange-3),LineStyle(7),LineColor(kOrange-3),LineWidth(3));
                model->plotOn(frame,Name(Form("sigF%d",_count)),Components(*sig),RooFit::NormRange("m_range"), Precision(1e-6),DrawOption("F"),FillStyle(3002),FillColor(kOrange-3),LineStyle(7),LineColor(kOrange-3),LineWidth(3));
        } else {
                model->plotOn(frame,Name(Form("sig%d",_count)), Components(sig1),NormRange("m_range"),Precision(1e-6),DrawOption("L"),FillStyle(3002),FillColor(kOrange-3),LineStyle(7),LineColor(kOrange-3),LineWidth(3));
                model->plotOn(frame,Name(Form("sigF%d",_count)),Components(sig1),NormRange("m_range"),Precision(1e-6),DrawOption("F"),FillStyle(3002),FillColor(kOrange-3),LineStyle(7),LineColor(kOrange-3),LineWidth(3));
        }
*/

 ds->plotOn(frame, Name(Form("ds_cut%d", _count)), Binning(nbinsmasshisto),
             MarkerSize(1), MarkerStyle(20), MarkerColor(1),
             LineColor(1), LineWidth(2));//draw an transparent hist
  if(npfit != "1") {
    TString option = (pdf == "mass_range")? "L" : "LF";
    RooCmdArg drawRange = (pdf == "mass_range")? Range(fitRange) : RooCmdArg();
    model->plotOn(frame,  Name(Form("peakbg%d", _count)) , Components(peakbg),
                  drawRange, Precision(1e-6), DrawOption(option), FillStyle(3005),
                  FillColor(kGreen+4), LineStyle(1), LineColor(kGreen+4), LineWidth(3));
  }

	model->plotOn(frame, Name(Form("bkg%d", _count)) ,  Components(bkg),
                Range(fitRange), Precision(1e-6), DrawOption("L"),
                LineStyle(7), LineColor(4), LineWidth(3));
	model->plotOn(frame, Name(Form("model%d", _count)),
                Range(fitRange),  Precision(1e-6), DrawOption("L"),
                LineColor(2), LineWidth(3));

	if(pdf!="1gauss") {
		model->plotOn(frame, Name(Form("sig%d", _count)),  Components(*sig),
                  Range(fitRange), Precision(1e-6), DrawOption("LF"),
                  FillStyle(3002), FillColor(kOrange-3),
                  LineStyle(7), LineColor(kOrange-3), LineWidth(3));
	} else {
		model->plotOn(frame, Name(Form("sig%d", _count)),  Components(sig1),
                  NormRange(fitRange), Precision(1e-6), DrawOption("L"), FillStyle(3002),
                  FillColor(kOrange-3), LineStyle(7), LineColor(kOrange-3), LineWidth(3));
		model->plotOn(frame, Name(Form("sigF%d", _count)), Components(sig1),
                  NormRange(fitRange), Precision(1e-6), DrawOption("F"), FillStyle(3002),
                  FillColor(kOrange-3), LineStyle(7), LineColor(kOrange-3), LineWidth(3));
	}




  if(drawLegend){model->paramOn(frame,Layout(1, 1, 1), Format("NEU",AutoPrecision(3)));}   //this one does not print parameters
  else{model->paramOn(frame,Layout(0.6, x_2, y_1), Format("NEU",AutoPrecision(3)));}

	frame->getAttText()->SetTextSize(0.00);
	frame->getAttFill()->SetFillStyle(0);
	frame->getAttLine()->SetLineWidth(0);
	frame->SetTitle("");
	frame->SetYTitle("Events / (20 MeV/c^{2})");

	frame->GetXaxis()->CenterTitle();
	frame->GetYaxis()->CenterTitle();
	frame->GetXaxis()->SetTitleOffset(1.1);
	frame->GetYaxis()->SetTitleOffset(1.4);
	frame->GetXaxis()->SetTitleSize(0.04);
	frame->GetYaxis()->SetTitleSize(0.05);
	frame->GetXaxis()->SetLabelFont(42);
	frame->GetYaxis()->SetTitleFont(42);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetLabelSize(0.045);
	frame->GetYaxis()->SetLabelSize(0.05);
	frame->SetStats(0);
	(frame->GetXaxis())->SetRangeUser(minhisto,maxhisto);
//	(frame->GetYaxis())->SetRangeUser(0,nsig.getVal()*0.8);

	frame->GetXaxis()->SetNdivisions(-50205);	
	frame->Draw();
	
if(tree == "ntKpi"){
if(ptmin==5 && ptmax==60){ (frame->GetYaxis())->SetRangeUser(0,12100);}
else if (ptmin == 5) { (frame->GetYaxis())->SetRangeUser(0,900);}
else if (ptmin == 7) { (frame->GetYaxis())->SetRangeUser(0,1200);}
else if (ptmin == 10) { (frame->GetYaxis())->SetRangeUser(0,5000);}
else if (ptmin == 15) { (frame->GetYaxis())->SetRangeUser(0,3000);}
else if (ptmin == 20) { (frame->GetYaxis())->SetRangeUser(0,2000);}
else if (ptmin == 30) { (frame->GetYaxis())->SetRangeUser(0,600);}
else if (ptmin == 50) { (frame->GetYaxis())->SetRangeUser(0,55);}
} else if (tree == "ntphi") {
if(ptmin==7 && ptmax==50){ (frame->GetYaxis())->SetRangeUser(0,1200);}
else if (ptmin == 7) { (frame->GetYaxis())->SetRangeUser(0,100);}
else if (ptmin == 10) { (frame->GetYaxis())->SetRangeUser(0,510);}
else if (ptmin == 15) { (frame->GetYaxis())->SetRangeUser(0,300);}
else if (ptmin == 20) { (frame->GetYaxis())->SetRangeUser(0,340);}}


	RooHist* pull_hist = frame->pullHist(Form("ds_cut%d",_count),Form("model%d",_count));
	//  RooHist* pull_hist = frame->pullHist("Data","Fit");

	RooPlot* pull_plot = mass->frame();
	(pull_plot->GetXaxis())->SetRangeUser(minhisto,maxhisto); //maxhisto
	RooRealVar x("x","",1e3,0,1e6);
	x.setVal(0.);
	RooGenericPdf* line_ref = new RooGenericPdf("ref_0", "ref_0", "x", RooArgList(x));
	line_ref->plotOn(pull_plot, LineStyle(7), LineColor(13), LineWidth(2));  

	pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"XP");
	pull_plot->SetTitle("");

	if(tree=="ntKp")pull_plot->SetXTitle("m_{J/#psiK^{#pm}} (GeV/c^{2})");
	if(tree=="ntphi")pull_plot->SetXTitle("m_{J/#psiK^{+}K^{-}} (GeV/c^{2})");

	pull_plot->SetYTitle("Pull");
	pull_plot->GetYaxis()->SetTitleFont(42);  
	pull_plot->GetYaxis()->SetTitleSize(0.17);
	pull_plot->GetYaxis()->CenterTitle();
	pull_plot->GetYaxis()->SetLabelOffset(0.01);
	pull_plot->GetYaxis()->SetLabelSize(0.15);
	pull_plot->GetYaxis()->SetNdivisions(305);

	pull_plot->GetXaxis()->SetTitleSize(0.20);
	pull_plot->GetXaxis()->SetTitleOffset(1.0);
	pull_plot->GetXaxis()->CenterTitle();
	pull_plot->GetXaxis()->SetLabelFont(42);
	pull_plot->GetXaxis()->SetLabelOffset(0.01);
	pull_plot->GetXaxis()->SetLabelSize(0.15);
	pull_plot->GetXaxis()->SetTickLength(0.15);
	pull_plot->GetXaxis()->SetNdivisions(-50205);

	std::cout<<"pull done"<<std::endl;
	/*
	   mass->setRange("signal",5.2,5.6);
	   RooAbsReal* sigIntegral = sig->createIntegral(*mass,NormSet(*mass),Range("signal"));
	   double sigIntegralErr = sigIntegral->getPropagatedError(*fitResult);
	   RooAbsReal* bkgIntegral = bkg.createIntegral(*mass,NormSet(*mass),Range("signal"));
	   double bkgIntegralErr = bkgIntegral->getPropagatedError(*fitResult);
	   cout<<"nsig: "<<nsig.getVal()<<endl;
	//cout<<"nsig: "<<model->getParameters(*mass)->getRealValue("nsig")<<endl;
	cout<<"nsig error: "<<nsig.getError()<<endl;
	cout<<"sig integral: "<<sigIntegral->getVal()<<endl;
	cout<<"sig integral error: "<<sigIntegralErr<<endl;
	cout<<"nbkg: "<<nbkg.getVal()<<endl;
	//cout<<"nbkg: "<<model->getParameters(*mass)->getRealValue("nbkg")<<endl;
	cout<<"nbkg error: "<<nbkg.getError()<<endl;
	cout<<"bkg integral: "<<bkgIntegral->getVal()<<endl;
	cout<<"bkg integral error: "<<bkgIntegralErr<<endl;*/
	//cout<<"#expected events: "<<model->expectedEvents(RooArgSet(*mass))<<endl;
	//cout<<"model integral: "<<model->getNormIntegral(RooArgSet(*mass))->getVal()<<endl;
	//fitResult->Print("v");

	//RooHist* datahist = new RooHist();
	//datahist = frame->getHist("ds");
	//TGraphAsymmErrors* datagraph = static_cast<TGraphAsymmErrors*>(datahist);
	//RooCurve* modelcurve = new RooCurve();
	//modelcurve = frame->getCurve("model");

	Double_t yield = nsig.getVal();
	Double_t yieldErr = nsig.getError();
	TH1D* fh = (TH1D*)h->Clone("fh");
	double nexpected = model->expectedEvents(RooArgSet(*mass));
	double dataArr[nbinsmasshisto]; double dataErrArr[nbinsmasshisto]; double fitArr[nbinsmasshisto];
	double val = 0;
	for(int i = 0; i < nbinsmasshisto; i++){
		dataArr[i] = h->GetBinContent(i+1);
		dataErrArr[i] = h->GetBinError(i+1);
		mass->setVal(h->GetBinCenter(i+1));
		val = model->getVal(RooArgSet(*mass))*1./nbinsmasshisto*nexpected;
		fitArr[i] = val;
		fh->SetBinContent(i+1, fitArr[i]);
		fh->SetBinError(i+1, sqrt(fitArr[i]));
	}
	cout<<"frame->chiSquare: "<<frame->chiSquare(Form("model%d",_count),Form("ds_cut%d",_count), fitResult->floatParsFinal().getSize())<<endl;
	cout<<"pars= "<<fitResult->floatParsFinal().getSize()<<endl;
	cout<<"chi2: "<<frame->chiSquare(Form("model%d",_count),Form("ds%d",_count))<<endl;
	//RooChi2Var chi2_lowstat("chi2_lowstat","chi2",*model,*dh);
	//cout<<chi2_lowstat.getVal()<<endl;

	double chiRoo = frame->chiSquare(Form("model%d",_count),Form("ds_cut%d",_count), fitResult->floatParsFinal().getSize());
	double chi2Std = 0;
	double chi2Neyman = 0;
	double chi2Peason = 0;
	double chi2BakerCousins = 0;
	chi2Cal(dataArr, dataErrArr, fitArr, nbinsmasshisto, chi2Std, chi2Neyman, chi2Peason, chi2BakerCousins);
	printf("chi2 Standard: %f\n",chi2Std);
	printf("chi2 Neyman: %f\n",chi2Neyman);
	printf("chi2 Peason: %f\n",chi2Peason);
	printf("chi2 Baker & Cousins: %f\n",chi2BakerCousins);

	TLegend *leg = new TLegend(0.62,0.55,0.89,0.75,NULL,"brNDC"); 
	leg = new TLegend(0.65/1.20+0.10,0.42,0.89+0.10,0.90,NULL,"brNDC");
	
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(h,"Data","pe");
	leg->AddEntry(frame->findObject(Form("model%d",_count)),"Fit","l");
	leg->AddEntry(frame->findObject(Form("sig%d",_count)),"Signal","f");
	leg->AddEntry(frame->findObject(Form("bkg%d",_count)),"Background","l");
	if(npfit != "1") leg -> AddEntry(frame->findObject(Form("peakbg%d",_count)),"B #rightarrow J/#psi X","f");
  	if (drawLegend) {leg -> Draw();}

	float posit = 0.58;
	int nDOF = nbinsmasshisto-(fitResult->floatParsFinal().getSize());
	float nChi2 = chi2BakerCousins/(nbinsmasshisto-(fitResult->floatParsFinal().getSize()));
	int nDigit_chi2BakerCousins = 2;
	int nDigit_nChi2 = 2;
	chi2BakerCousins = roundToNdigit(chi2BakerCousins);
	nChi2 = roundToNdigit(nChi2);
	nDigit_chi2BakerCousins = sigDigitAfterDecimal(chi2BakerCousins);
	nDigit_nChi2 = sigDigitAfterDecimal(nChi2);
	TLatex* texChi = new TLatex(0.65,0.64, Form("#chi^{2}/nDOF = %2f", chiRoo));
	texChi->SetNDC();
	texChi->SetTextAlign(12);
	texChi->SetTextSize(0.03);
	texChi->SetTextFont(42);

	nsig.setVal(0.);
	nsig.setConstant();
	RooFitResult* fitResult_nosig = model->fitTo(*ds,Save());
	RooAbsReal* nll_nosig = model->createNLL(*ds);
	double log_likelihood_nosig= nll_nosig->getVal();

        real_significance = sqrt(2*(-log_likelihood+log_likelihood_nosig));
	std::cout<<"REAL SIGNIFICANCE= "<<real_significance<<std::endl;
	std::cout<<"***********************************************************************"<<std::endl;

	std::cout<<"sigma1_MC error hi= "<<sigma1MC.getAsymErrorHi()<<std::endl;
	std::cout<<"sigma1_MC error = "<<sigma1MC.getError()<<std::endl;
	std::cout<<"sigma2_MC error hi= "<<sigma2MC.getAsymErrorHi()<<std::endl;
	std::cout<<"sigma2_MC error = "<<sigma2MC.getError()<<std::endl;
	//RooFitResult* fitResult = model->fitTo(*ds,Save(),Minos());
	//ds->plotOn(frame,Name(Form("ds%d",_count)),Binning(nbinsmasshisto),MarkerSize(1.55),MarkerStyle(20),LineColor(1),LineWidth(4));


	// double width = 0.08;
	// double BmassH = BS_MASS + width;
	// double BmassL = BS_MASS - width;
	// mass->setRange("signal",BmassL,BmassH);
	RooAbsReal *bkgIntegral = bkg.createIntegral(*mass,NormSet(*mass),Range("signal"));
	// bkgIntegralErr = bkgIntegral->getPropagatedError(*fitResult);	
	cout<<"bkg integral: "<<bkgIntegral->getVal()<<endl;
	//  cout<<"bkg integral error: "<<bkgIntegralErr<<endl;

	bkgd = nbkg.getVal()*bkgIntegral->getVal();
	Significance = yield/sqrt(bkgd+yield);

	Double_t Significance =  real_significance;
	//	Double_t Significance =  yield/TMath::Sqrt(bkgd+yield);
	int nDigit_Significance = 3;
	Significance = roundToNdigit(Significance);
	nDigit_Significance = sigDigitAfterDecimal(Significance);
	TLatex* texSig = new TLatex(0.65,0.67,Form("Significance = %.*f", nDigit_Significance, Significance));
	cout<<"Significance = "<<Significance<<endl;
	texSig->SetNDC();
	texSig->SetTextFont(42);
	texSig->SetTextSize(0.03);
	texSig->SetLineWidth(2);


	cMC->cd();
//	leg->Draw("same");

	p2->cd();
	pull_plot->GetYaxis()->SetTitleOffset(0.4);
	pull_plot->SetYTitle("Pull");
	pull_plot->Draw();

	p1->cd();
	outframe = frame;
	outputw->import(*model);

	cout << "------------------------------------------------------------------------------------------------" << endl;

	Double_t yieldPrintErr = nsig.getError();
	Double_t yieldPrintErrUp = nsig.getAsymErrorHi();
	Double_t yieldPrintErrDown = -1 * nsig.getAsymErrorLo();

	cout << "yield Error = " << yieldPrintErr << "     yield Error Up = " << yieldPrintErrUp << "    yieldPrintErrDown = " << yieldPrintErrDown << endl;

	cout << "------------------------------------------------------------------------------------------------" << endl;

	RooAbsReal* RangeBakground = bkg.createIntegral(*mass,NormSet(*mass),Range("signal")); 

	//cout << "RangeBakground = " << RangeBakground << endl;

	double Calback = RangeBakground->getVal() * nbkg.getVal();
	double StatSig = yield/sqrt(yield + Calback);
/*
	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.045);
	lat->DrawLatex(0.63,0.55,Form("S = %.1f",yield));	
	lat->DrawLatex(0.63,0.50,Form("B = %.1f",Calback));	
	lat->DrawLatex(0.63,0.45,Form("S/#sqrt{S+B} = %.1f",StatSig));	
	lat->Draw("SAME");
*/
	return fitResult;
}

void clean0(TH1D* h)
{
	for (int i=1;i<=h->GetNbinsX();i++)
	{
		if(h->GetBinContent(i)==0) h->SetBinError(i,1);
	}
}

double ErrorOnSigma(double width, double errwidth, double smear, double errsmearing){
	double squarederroronsigma=(1+smear)*(1+smear)*errwidth*errwidth+width*width*errsmearing*errsmearing;
	double erroronsigma=TMath::Sqrt(squarederroronsigma);
	return erroronsigma;}


void latex_table(std::string filename, int n_col, int n_lin, std::vector<std::string> col_name, std::vector<std::string> labels, 
		std::vector<std::vector<double> > numbers, std::string caption){
	
	std::ofstream file_check;
	std::ofstream file;

	//Begin Document
                                                                                    
	file.open(filename + ".tex");
	file_check.open(filename + "_check.tex");

	file_check << "\\documentclass{article}" << std::endl;
	//file << "\\usepackage[utf8]{inputenc}" << std::endl;     
	file_check << "\\usepackage{rotating}" << std::endl;                                                                                   
	// file_check << "\\usepackage{cancel}" << std::endl;
	file_check << "\\usepackage{geometry}" << std::endl;
	file_check << "\\usepackage{booktabs}" << std::endl;
	file_check << "\\geometry{a4paper, total={170mm,257mm}, left=20mm, top=20mm,}" << std::endl;

	file_check << "\\title{Bs/B+}" << std::endl;
	file_check << "\\author{Henrique Legoinha}" << std::endl;
	file_check << "\\date{June 2022}" << std::endl;
	file_check << "\\begin{document}" << std::endl;
	file_check << "\\maketitle" << std::endl;

	// Create table                                                                                                                                
	//file_check << "\\begin{table}[!h]" << std::endl;
	//file_check << "\\centering" << std::endl;                                                                                                         
	//setup table size                                                                                                                             
	std::string col="c";
	//std::setprecision(3);

	for(int i=1; i<n_col; i++)
		col+="|c";

	file_check << "\\begin{sidewaystable}"<< std::endl;
	file_check << "\\begin{tabular}{"+col+"}" << std::endl;
	file_check << "\\toprule" << std::endl;
	file << "\\begin{tabular}{"+col+"}" << std::endl;
	file << "\\toprule" << std::endl;

	for(int c=0; c<n_col-1; c++){
		file << col_name[c] << " & ";
		file_check << col_name[c] << " & ";
	}

	file << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;
	file_check << col_name[n_col-1] << " \\\\ \\midrule" << std::endl;


	for(int i=1; i<n_lin; i++)
	{
		file << labels[i-1] << " & ";
		file_check << labels[i-1] << " & ";

		for(int c=1; c<n_col-1; c++){
			file << std::setprecision(3)<< numbers[c-1][i-1]<< " \\% & ";
			file_check << std::setprecision(3) << numbers[c-1][i-1]<< " \\% & ";

		}

		file << std::setprecision(3)<< numbers[n_col-2][i-1]<< " \\% \\\\" << std::endl;
		file_check << std::setprecision(3)<< numbers[n_col-2][i-1]<< " \\% \\\\" << std::endl;
	}

	file << "\\bottomrule" << std::endl;
	file_check << "\\bottomrule" << std::endl;

	//End Table                                                                                                                                    
	file << "\\end{tabular}" << std::endl;
	file_check << "\\end{tabular}" << std::endl;
	file_check << "\\caption{"+caption+"}" << std::endl;
	file_check << "\\end{sidewaystable}"<< std::endl;

	//file_check << "\\end{table}" << std::endl;
	//End document                                                                                                                                 

	file_check << "\\end{document}" << std::endl;

	file.close();
	file_check.close();
	cout<<filename<<endl;
	system(("pdflatex " + filename + "_check.tex").c_str());
	system((filename + "_check.pdf").c_str());
}

/*
   void validate_fit(RooWorkspace* w, int particle, TString ModName, TString SigName)
   {

   std::cout << "Now Perform Check on Fit" << std::endl;



   RooRealVar Bmass = *(w->var("Bmass"));
//RooAbsPdf* model  = w->pdf("model");
RooAbsPdf* model  = w->pdf(ModName.Data());


vector<RooRealVar> params;
//	params.push_back(*(w->var("n_signal")));

std::cout << "Work Here Before Name" << std::endl;

params.push_back(*(w->var(SigName.Data())));

const int params_size = params.size();




std::cout << "Work Here Before ROOMC" << std::endl;


RooMCStudy* mcstudy = new RooMCStudy(*model, Bmass, Binned(kTRUE), Silence(), Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));

mcstudy->generateAndFit(5000);

std::cout << "Work Here Before VEc" << std::endl;


vector<RooPlot*> framesPull, framesParam;

for(int i = 0; i < params_size; ++i)
{
framesPull.push_back(mcstudy->plotPull(params.at(i),FrameBins(200),FrameRange(-5,5)));
framesPull[i]->SetTitle("");
framesParam.push_back(mcstudy->plotParam(params.at(i),FrameBins(50)));
framesParam[i]->SetTitle("");
}

vector<TGraph*> his;

for(int i = 0; i < params_size; ++i){
his.push_back(static_cast<TGraph*>(framesPull.at(i)->getObject(0)));
}

gStyle->SetOptFit(0111);

TCanvas* c_pull = new TCanvas("pulls", "pulls", 900, 800);

gPad->SetLeftMargin(0.15);

for(int i = 0; i < params_size; ++i){
c_pull->cd();
his[i]->SetTitle("");
his[i]->Draw();
c_pull->Update();
his[i]->Fit("gaus","","",-5,5);
his[i]->GetFunction("gaus")->SetLineColor(4);
his[i]->GetFunction("gaus")->SetLineWidth(5);
his[i]->GetXaxis()->SetTitle("Pull");
his[i]->GetYaxis()->SetTitle("Toy MCs");
his[i]->Draw("same");
}

TCanvas* c_params = new TCanvas("params", "params", 900, 800);

for(int i = 0; i < params_size; ++i){
	c_params->cd();
	framesParam.at(i)->GetYaxis()->SetTitleOffset(1.4);
	framesParam.at(i)->Draw();
}

if(particle == 0){
	c_pull->SaveAs("./results/Bu/pulls/pulls_poisson_Bu.pdf");
	c_pull->SaveAs("./results/Bu/pulls/pulls_poisson_Bu.gif");
	c_params->SaveAs("./results/Bu/pulls/pulls_params_poisson_Bu.pdf");
	c_params->SaveAs("./results/Bu/pulls/pulls_params_poisson_Bu.gif");
}else if(particle == 1){
	c_pull->SaveAs("./results/Bs/pulls/pulls_poisson_Bs.pdf");
	c_pull->SaveAs("./results/Bs/pulls/pulls_poisson_Bs.gif");
	c_params->SaveAs("./results/Bs/pulls/pulls_params_poisson_Bs.pdf");
	c_params->SaveAs("./results/Bs/pulls/pulls_params_poisson_Bs.gif");
}

}
*/



void validate_fit(RooWorkspace* w, TString tree, TString variable, int full, int CentMin, int CentMax, int ptMin, int ptMax)
{
	std::cout << "Now Perform Check on Fit" << std::endl;
	RooRealVar Bmass = *(w->var("Bmass"));
	RooAbsPdf* model  = w->pdf(Form("model%d",_count));
	//RooDataSet* data = (RooDataSet*) w->data("data");

	//model->fitTo(*data);

	vector<RooRealVar> params;
	params.push_back(*(w->var(Form("nsig%d",_count))));
	params.push_back(*(w->var(Form("mean%d",_count))));


	/*
	//Retrieving Parameters//

	RooRealVar* lambda =  (w->var(Form("lambda%d", _count)));
	RooRealVar* nbkg   =  (w->var(Form("nbkg%d", _count)));
	RooRealVar * mean  = (w->var(Form("mean%d", _count)));*/
	RooRealVar * sigma1  = (w->var(Form("sigma1%d", _count)));
	RooRealVar * sigma2  = (w->var(Form("sigma2%d", _count)));
	std::cout << "sigma1 is " << sigma1->getVal() << std::endl; 
	//	RooRealVar * c1  = (w->var(Form("c1%d", _count)));
/*	RooRealVar * sig1frac = (w->var(Form("sig1frac%d", _count)));




	lambda->setVal(lambda->getVal());
	nbkg->setVal(nbkg->getVal());
	mean->setVal(mean->getVal());*/
	sigma1->setVal(sigma1->getVal());
	sigma2->setVal(sigma2->getVal());
	//	c1->setVal(c1->getVal());
/*	sig1frac->setVal(sig1frac->getVal());
	*/




	/*
	// test removing bakcground
	RooRealVar* lambda =  (w->var(Form("lambda%d", _count)));
	RooRealVar* nbkg   =  (w->var(Form("nbkg%d", _count)));

	//	double step = 0.002;


	//	double nkgdValue = 51.15;
	//	double lambdaValue = -2.0236;


	cout << "BRO lambda = " << lambda->getVal() << endl;
	cout << "BRIS back = " << nbkg->getVal() << endl;

	double nkgdValue = nbkg->getVal();
	double lambdaValue = lambda->getVal();


	nbkg->setVal(50);
	lambda->setVal(lambdaValue);
	*/
	//	RooRealVar* nbkg   =  (w->var(Form("nbkg%d", _count)));
	//	nbkg->setVal(50);

	double n_signal_init = params[0].getVal();
	double n_signal_error_init = params[0].getError();

	cout << "N_signal initial value: " << n_signal_init << endl;
	cout << "N_signal initial Error value: " << n_signal_error_init << endl;



	int params_size = params.size();

	cout << "params_size " << params_size << endl;

	RooMCStudy* mcstudy = new RooMCStudy(*model, Bmass,  Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));

	mcstudy->generateAndFit(5000);

	cout << "DONE Generate and Fit " << endl;


	TString XName[2] = {"nsig","mean"};

	vector<RooPlot*> framesPull, framesParam, framesError;

	for(int i = 0; i < params_size; ++i)
	{
		framesPull.push_back(mcstudy->plotPull(params.at(i),FrameBins(50),FrameRange(-5,5)));
		framesPull[i]->SetTitle("");
		framesParam.push_back(mcstudy->plotParam(params.at(i),FrameBins(50)));
		framesParam[i]->SetTitle("");
		framesError.push_back(mcstudy->plotError(params.at(i),FrameBins(50)));
		framesError[i]->SetTitle("");
	}

	vector<TGraph*> h1;
	vector<TGraph*> h2;
	vector<TGraph*> h3;

	for(int i = 0; i < params_size; ++i){
		h1.push_back(static_cast<TGraph*>(framesPull.at(i)->getObject(0)));
		h2.push_back(static_cast<TGraph*>(framesParam.at(i)->getObject(0)));
		h3.push_back(static_cast<TGraph*>(framesError.at(i)->getObject(0)));
	}

	gStyle->SetOptFit(0111);

	TCanvas* c_pull = new TCanvas("pulls", "pulls", 900, 800);

	gPad->SetLeftMargin(0.15);
	//gPad->SetBottomMargin(0.0);


	for(int i = 0; i < params_size; ++i){
		c_pull->cd();
		h1[i]->SetTitle("");
		h1[i]->Draw();
		c_pull->Update();
		h1[i]->Fit("gaus","","",-5,5);
		h1[i]->GetFunction("gaus")->SetLineColor(4);
		h1[i]->GetFunction("gaus")->SetLineWidth(5);
		h1[i]->GetXaxis()->SetTitle("Pull");
		h1[i]->GetYaxis()->SetTitle("");
		h1[i]->Draw("APsame");
	}


	TCanvas* c_params = new TCanvas("params", "params", 900, 800);
	for(int i = 0; i < params_size; ++i){
		c_params->cd();
		framesParam.at(i)->GetYaxis()->SetTitle("Toy MCs");
		framesParam.at(i)->GetYaxis()->SetTitleOffset(1.4);
		framesParam.at(i)->Draw();
	}


	/*TCanvas* c_params = new TCanvas("params", "params", 900, 800);

	  for(int i = 0; i < params_size; ++i){
	  c_params->cd();
	  h2[i]->SetTitle("");
	  h2[i]->Draw();
	  c_params->Update();
	  if(particle == 0){h2[i]->Fit("gaus","","",900,1200);}
	  else if(particle == 1){h2[i]->Fit("gaus","","",60, 120);}
	  h2[i]->GetFunction("gaus")->SetLineColor(4);
	  h2[i]->GetFunction("gaus")->SetLineWidth(5);
	//h2[i]->GetXaxis()->SetTitle();
	h2[i]->GetYaxis()->SetTitle("Toy MCs");
	h2[i]->Draw("same");
	}*/



	TCanvas* c_errors = new TCanvas("errors", "errors", 900, 800);
	gPad->SetLeftMargin(0.15);
	for(int i = 0; i < params_size; ++i){
		c_errors->cd();
		framesError.at(i)->GetYaxis()->SetTitle("Toy MCs");
		framesError.at(i)->GetYaxis()->SetTitleOffset(1.4);
		framesError.at(i)->Draw();
	}

	for(int i = 0; i < params_size; ++i){

		n_signal_init = params[i].getVal();
		n_signal_error_init = params[i].getError();

		c_pull->cd();
		h1[i]->SetTitle("");
		h1[i]->Draw();
		c_pull->Update();
		h1[i]->Fit("gaus","","",-5,5);
		h1[i]->GetFunction("gaus")->SetLineColor(4);
		h1[i]->GetFunction("gaus")->SetLineWidth(5);
		h1[i]->GetXaxis()->SetTitle(Form("%s Pull",XName[i].Data()));
		h1[i]->GetYaxis()->SetTitle("");
		h1[i]->Draw("APsame");
		c_errors->cd();

		h3[i]->SetTitle("");
		h3[i]->Draw();
		c_errors->Update();
		h3[i]->Fit("gaus","","",n_signal_error_init*0.7,n_signal_error_init*1.3);
		h3[i]->GetFunction("gaus")->SetLineColor(4);
		h3[i]->GetFunction("gaus")->SetLineWidth(5);
		h3[i]->GetXaxis()->SetTitle(Form("%s Error",XName[i].Data()));
		h3[i]->GetYaxis()->SetTitle("");
		h3[i]->Draw("APsame");



		c_params->cd();

		h2[i]->SetTitle("");
		h2[i]->Draw();
		c_params->Update();
		h2[i]->Fit("gaus","","",n_signal_init * 0.7,n_signal_init * 1.3);
		h2[i]->GetFunction("gaus")->SetLineColor(4);
		h2[i]->GetFunction("gaus")->SetLineWidth(5);
		h2[i]->GetXaxis()->SetTitle(Form("%s Mean",XName[i].Data()));
		h2[i]->GetYaxis()->SetTitle("");
		h2[i]->Draw("APsame");


		if(tree=="ntKp"){
			//c_pull->SaveAs("./mcstudy/pulls_poisson_Bu.pdf");
			c_pull->SaveAs(Form("new/pull_signal_full%d_%s_%d_%s.png",full,variable.Data(),_count,tree.Data()));
			c_params->SaveAs(Form("new/params_signal_full%d_%s_%d_%s.png",full,variable.Data(),_count,tree.Data()));
			c_errors->SaveAs(Form("new/error_signal_full%d_%s_%d_%s.png",full,variable.Data(),_count,tree.Data()));
			//c_errors->SaveAs("./results/Bu/pulls/pulls_error.gif");
		}
		else if(tree=="ntphi"){
			c_pull->SaveAs(Form("newFIT/pull_signal_%s_%d_%d_%d_%d_%d_Bs.png",variable.Data(),CentMin,CentMax,ptMin,ptMax,i));
			//c_pull->SaveAs("./mcstudy/pulls_poisson_Bs.pdf");
			c_params->SaveAs(Form("newFIT/param_signal_%s_%d_%d_%d_%d_%d_Bs.png",variable.Data(),CentMin,CentMax,ptMin,ptMax,i));
			c_errors->SaveAs(Form("newFIT/error_signal_%s_%d_%d_%d_%d_%d_Bs.png",variable.Data(),CentMin,CentMax,ptMin,ptMax,i));
			//c_errors->SaveAs("./results/Bs/pulls/pulls_error.gif");
		}
	}
	/*
	   else if(tree=="ntphi"){
	   c_pull->SaveAs(Form("new/StepScan/pull/pull_signal_full%d_%s_%d_%d_Bs.png",full,variable.Data(),_count,NTrial));
//c_pull->SaveAs("./mcstudy/pulls_poisson_Bs.pdf");
c_params->SaveAs(Form("new/StepScan/param/param_signal_full%d_%s_%d_%d_Bs.png",full,variable.Data(),_count,NTrial));
c_errors->SaveAs(Form("new/StepScan/error/error_signal_full%d_%s_%d_%d_Bs.png",full,variable.Data(),_count,NTrial));
//c_errors->SaveAs("./results/Bs/pulls/pulls_error.gif");
}
*/
}


