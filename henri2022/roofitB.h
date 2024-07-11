#include "TAxis.h"
#include "uti.h"
#include "TSystem.h"
#include "RooCBShape.h"
#include "RooWorkspace.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooFormulaVar.h"
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

void plot_jpsifit(RooWorkspace& w, int nbin_hist, TString pdf, RooAbsPdf* model, RooDataSet* ds, TString plotName,  bool with_sig);
void fix_parameters(RooWorkspace& w, TString pdfName, bool release=false);
void fit_jpsinp (RooWorkspace& w, int nbin_hist, TString pdf, float bin_i, float bin_f, TString Var, TString bsbp);
template<typename... Targs>
void plot_mcfit(RooWorkspace& w, RooAbsPdf* model, RooDataSet* ds, TString plotName,  Targs... options);

using namespace RooFit;
using namespace std;

#define BS_MASS 5.36682
#define BP_MASS 5.27915

double minhisto=5.;
double maxhisto=6.;
int nbinsmasshisto=100;
Int_t _count=0;

RooFitResult *fit(TString variation, TString pdf,TString tree, TCanvas* c, TCanvas* cMC, RooDataSet* ds, RooDataSet* dsMC, RooRealVar* mass, float binmin, float binmax, RooWorkspace& w, TString which_var, TString bsbpbins ){
	
	//if (tree == "ntphi"){nbinsmasshisto = 50;} //100 bins is to much for bs case
	if ( which_var == "Bpt" && ( (int)binmin == 50 & (int)binmax == 60) ){nbinsmasshisto = 50;} //100 bins is to much for bp 50-60 mass bin case
	
	cMC->cd();
	RooPlot* frameMC = mass->frame();
	frameMC->SetTitle("");
	frameMC->GetYaxis()->SetTitle(TString::Format("Events / (%g MeV/c^{2})",(mass->getMax()-mass->getMin())/nbinsmasshisto*1000));
	frameMC->GetYaxis()->SetTitleOffset(2.);
	frameMC->GetYaxis()->SetTitleSize(0.035);
	frameMC->GetXaxis()->SetTitleFont(42);
	frameMC->GetYaxis()->SetTitleFont(42);
	frameMC->GetXaxis()->SetLabelFont(42);
	frameMC->GetYaxis()->SetLabelFont(42);
	frameMC->GetXaxis()->SetLabelSize(0.035);
	frameMC->GetYaxis()->SetLabelSize(0.035);
	frameMC->SetStats(0);
	frameMC->GetXaxis()->SetNdivisions(-50205);

	TPad *pMC1 = new TPad("pMC1","pMC1",0.,0.215,1.,1);
	// TPad *p1 = new TPad("p1","p1",0.05,0.05,0.99,0.99);
	pMC1->SetBorderMode(1); 
	pMC1->SetFrameBorderMode(0); 
	pMC1->SetBorderSize(2);
	pMC1->SetBottomMargin(0.10);
	pMC1->Draw(); 

	TPad *pMC2 = new TPad("pMC2","pMC2",0.,0.02,1.,0.24);
	pMC2->SetTopMargin(0.);    
	pMC2->SetBorderMode(0);
	pMC2->SetBorderSize(2); 
	pMC2->SetFrameBorderMode(0); 
	pMC2->SetTicks(1,1); 
	//pMC2->Draw();

	// give better initial values to the Bmesons mass
	double init_mean = BS_MASS;
	double m1 = 5.35;
	double m2 = 5.375;
	if (tree == "ntKp"){
		init_mean = BP_MASS;
		m1 = 5.27;
		m2 = 5.29;}
	// give better initial values to the Bmesons mass

	RooRealVar meanMC(Form("meanMC%d_%s",_count,pdf.Data()),"",init_mean,m1,m2) ;
	RooRealVar sigma1MC(Form("sigma1MC%d_%s",_count,pdf.Data()),"",0.05,0.005,0.11) ;
	RooRealVar sigma2MC(Form("sigma2MC%d_%s",_count, pdf.Data()),"",0.03,0.01,0.06) ;
	RooRealVar sigma3MC(Form("sigma3MC%d_%s",_count, pdf.Data()),"",0.01,0.005,0.2) ;
	RooRealVar sigma4cbMC(Form("sigma4cbMC%d_%s",_count, pdf.Data()),"",0.0266,0.01,0.5) ;
	RooRealVar sigma5cbMC(Form("sigma5cbMC%d_%s",_count, pdf.Data()),"",0.0266,0.01,0.5) ;
	RooRealVar alphaMC(Form("alphaMC%d_%s",_count,pdf.Data()),"",4.,0,15);
	RooRealVar nMC(Form("nMC_%d_%s", _count, pdf.Data()),"",10,-100,200);
	RooRealVar alphaMC2(Form("alphaMC2%d_%s",_count,pdf.Data()),"",4.,0,15);
	RooRealVar nMC2(Form("nMC2_%d_%s", _count, pdf.Data()),"",10,-100,200);
	RooRealVar* scale = new RooRealVar("scale","scale",1,0,2);
	RooProduct scaled_sigma1MC(Form("scaled_sigma1MC%d_%s",_count,pdf.Data()),"scaled_sigma1MC", RooArgList(*scale,sigma1MC));
	RooProduct scaled_sigma2MC(Form("scaled_sigma2MC%d_%s",_count,pdf.Data()),"scaled_sigma2MC", RooArgList(*scale,sigma2MC));
	RooProduct scaled_sigma3MC(Form("scaled_sigma3MC%d_%s",_count,pdf.Data()),"scaled_sigma3MC", RooArgList(*scale,sigma3MC));
	RooProduct scaled_sigma4cbMC(Form("scaled_sigma4cbMC%d_%s",_count,pdf.Data()),"scaled_sigma4cbMC", RooArgList(*scale,sigma4cbMC));
	RooProduct scaled_sigma5cbMC(Form("scaled_sigma5cbMC%d_%s",_count,pdf.Data()),"scaled_sigma5cbMC", RooArgList(*scale,sigma5cbMC));
	RooGaussian sig1MC(Form("sig1MC%d_%s",_count,pdf.Data()),"",*mass,meanMC,scaled_sigma1MC);  
	RooGaussian sig2MC(Form("sig2MC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma2MC);  
	RooGaussian sig3MC(Form("sig3MC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma3MC);  
	RooCBShape  CBMC(Form("CBMC%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma4cbMC, alphaMC, nMC);
	RooCBShape  CBMC2(Form("CBMC2%d_%s",_count, pdf.Data()),"",*mass,meanMC,scaled_sigma5cbMC, alphaMC2, nMC2);
	RooRealVar sig1fracMC(Form("sig1fracMC%d_%s",_count, pdf.Data()),"", 0.2, 0.001, 1);
	RooRealVar sig2fracMC(Form("sig2fracMC%d_%s",_count, pdf.Data()),"", 0.7, 0.001, 1);
	RooRealVar nsigMC(Form("nsigMC%d_%s",_count, pdf.Data()),"",dsMC->sumEntries(), 0.9*dsMC->sumEntries(), 1.2 * dsMC->sumEntries());
	RooAddPdf* sigMC ;
	if((variation=="" && pdf=="") || variation== "background" || (variation=="signal" && pdf=="fixed" )) sigMC = new RooAddPdf(Form("sigMC%d_%s",_count,pdf.Data()),"",RooArgList(sig1MC,sig2MC),sig1fracMC);
	if(variation=="signal" && pdf=="3gauss") sigMC = new RooAddPdf(Form("sigMC%d_%s",_count, pdf.Data()), "", RooArgList(sig1MC, sig2MC, sig3MC), RooArgList(sig1fracMC, sig2fracMC), true);
	if(variation=="signal" && pdf=="gauss_cb") sigMC = new RooAddPdf(Form("sigMC%d_%s",_count, pdf.Data()), "", RooArgList(sig1MC, CBMC), sig1fracMC);
	if(variation=="signal" && pdf=="2cb") sigMC = new RooAddPdf(Form("sigMC%d_%s",_count, pdf.Data()), "", RooArgList(CBMC, CBMC2), sig1fracMC);

	RooAddPdf* modelMC;
	if((variation=="signal" && (pdf=="gauss_cb"|| pdf=="3gauss"|| pdf=="fixed"|| pdf=="2cb"))||variation=="background") modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(*sigMC),RooArgList(nsigMC));
	if(variation =="" && pdf=="") modelMC = new RooAddPdf(Form("modelMC%d_%s",_count, pdf.Data()),"",RooArgList(*sigMC),RooArgList(nsigMC));

//////////ROOFIT ROOFIT ROOFIT  MC MC MC MC 

	double subt = 0.14;
	if (tree == "ntphi") subt = .10;
	mass->setRange("signal",init_mean-subt, init_mean+subt);  //focus the MC fit to the signal region to prevent statistical flutuations
	scale->setConstant();
	RooFitResult * fitResultMC = modelMC->fitTo(*dsMC,Save(), Range("signal"));
	w.import(*modelMC);
	RooRealVar * fMCchi2_params = new RooRealVar(Form("ndfMC_%d_%s",_count, pdf.Data()), "", fitResultMC->floatParsFinal().getSize() );
	w.import(*fMCchi2_params);

	cout << "Signal MC FIT Done " << endl;

///////////ROOFIT ROOFIT ROOFIT MC MC MC MC

//PLOT MC FIT
//PLOT MC 

	pMC1->cd();
	dsMC->plotOn(frameMC,Name(Form("dsMC_cut%d",_count)),Binning(nbinsmasshisto),MarkerSize(0.5),MarkerStyle(8),LineColor(1),LineWidth(1));
	modelMC->plotOn(frameMC,Name(Form("modelMC%d_%s",_count, pdf.Data())),Precision(1e-6),DrawOption("LF"), FillStyle(3002), FillColor(kOrange-3), LineStyle(7),LineColor(kOrange-3),LineWidth(1));
	modelMC->paramOn(frameMC,Layout(0.2, 0.5, 0.70), Format("NEU",AutoPrecision(2)));
	frameMC->getAttText()->SetTextSize(0.035);
	frameMC->getAttFill()->SetFillStyle(0);
	frameMC->getAttLine()->SetLineWidth(0);
    if(tree=="ntKp")frameMC->SetXTitle("m_{J/#psiK^{+}} [GeV/c^{2}]");
	if(tree=="ntphi")frameMC->SetXTitle("m_{J/#psiK^{+}K^{-}} [GeV/c^{2}]");
	frameMC->GetXaxis()->SetTitleOffset(1.2);
	frameMC->GetXaxis()->SetTitleSize(0.035);
	frameMC->GetXaxis()->SetTitleFont(42);
	frameMC->GetXaxis()->CenterTitle();
	frameMC->GetXaxis()->SetLabelSize(0.035);
	frameMC->GetXaxis()->SetRangeUser(init_mean-subt, init_mean+subt);
	frameMC->Draw();
	cMC->RedrawAxis();
	
	TLatex* texB_BINS = new TLatex(0.5,0.5,"");
	if ( which_var == "Bpt"){texB_BINS = new TLatex(0.22,0.8, Form("%d < p_{T} <%d GeV", (int) binmin, (int)binmax));}
	else if( which_var == "By"){texB_BINS = new TLatex(0.22,0.8, Form("%0.1f < |y| <%0.1f GeV",binmin, binmax));}
	else if( which_var == "nMult"){texB_BINS = new TLatex(0.22,0.8, Form("%d < Mult <%d GeV",(int)binmin, (int)binmax));}
	texB_BINS->SetNDC();
	texB_BINS->SetTextFont(42);
	texB_BINS->SetTextSize(0.03);
	texB_BINS->SetLineWidth(1);
	texB_BINS->Draw();
	
	TLatex* texB_N_vs_count = new TLatex(0.22,0.75, Form("Bin Events: %i", int (w.var(Form("Events_in_MC_%d",_count))->getVal()) ));
	texB_N_vs_count->SetNDC();
	texB_N_vs_count->SetTextFont(42);
	texB_N_vs_count->SetTextSize(0.03);
	texB_N_vs_count->SetLineWidth(1);
	texB_N_vs_count->Draw();
	
	TLatex* texB = new TLatex(0.5,0.5,"");
	if(tree=="ntphi"){ texB = new TLatex(0.21,0.85, "B^{0}_{s}");}
	if(tree=="ntKp"){ texB = new TLatex(0.21,0.85, "B^{+}");}
	texB->SetNDC();
	texB->SetTextFont(62);
	texB->SetTextSize(0.04);
	texB->SetLineWidth(2);
	texB->Draw();

	//cMC->SetLogy();
	TLegend *legMC = new TLegend(0.8,0.77,0.89,0.89,NULL,"brNDC"); 
	legMC->SetBorderSize(0);
	legMC->SetTextSize(0.025);     
	legMC->SetTextFont(42);
	legMC->SetFillStyle(0);
	legMC->AddEntry(frameMC->findObject(Form("dsMC_cut%d", _count)), " MC","lp");
	legMC->AddEntry(frameMC->findObject(Form("modelMC%d_%s",_count,pdf.Data()))," Sig. PDF","f");
  	legMC -> Draw();

//PLOT MC
//PULL MC

	pMC2->cd();
	RooRealVar x("x","",0,0,1e6);
	RooGenericPdf* line_ref = new RooGenericPdf("ref_0", "ref_0", "x", RooArgList(x));
	RooHist* pull_histMC = frameMC->pullHist(Form("dsMC_cut%d",_count),Form("modelMC%d_%s",_count, pdf.Data()));
	pull_histMC->SetMarkerSize(0.5);
	RooPlot* pull_plotMC = mass->frame();
	(pull_plotMC->GetXaxis())->SetRangeUser(init_mean-subt, init_mean+subt);
	line_ref->plotOn(pull_plotMC, LineStyle(7), LineColor(kOrange-3), LineWidth(1)); 
	pull_plotMC->addPlotable(static_cast<RooPlotable*>(pull_histMC),"XP");
	pull_plotMC->SetTitle("");
	pull_plotMC->SetMarkerStyle(6);
	if(tree=="ntKp")pull_plotMC->SetXTitle("m_{J/#psiK^{+}} [GeV/c^{2}]");
	if(tree=="ntphi")pull_plotMC->SetXTitle("m_{J/#psiK^{+}K^{-}} [GeV/c^{2}]");
	pull_plotMC->SetYTitle("Pull");
	pull_plotMC->GetYaxis()->SetTitleFont(42);  
	pull_plotMC->GetYaxis()->SetTitleSize(0.15);
	pull_plotMC->GetYaxis()->CenterTitle();
	pull_plotMC->GetYaxis()->SetLabelOffset(0.01);
	pull_plotMC->GetYaxis()->SetLabelSize(0.1);
	pull_plotMC->GetYaxis()->SetNdivisions(305);
	pull_plotMC->GetYaxis()->SetTitleOffset(0.4);
	pull_plotMC->GetXaxis()->SetTitleSize(0.16);
	pull_plotMC->GetXaxis()->SetTitleOffset(1.0);
	pull_plotMC->GetXaxis()->CenterTitle();
	pull_plotMC->GetXaxis()->SetLabelFont(42);
	pull_plotMC->GetXaxis()->SetLabelOffset(0.01);
	pull_plotMC->GetXaxis()->SetLabelSize(0.1);
	pull_plotMC->GetXaxis()->SetTickLength(0.16);
	pull_plotMC->GetXaxis()->SetNdivisions(-50205);
	//pull_plotMC->Draw();

//PULL MC
//PLOT MC FIT

// FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp
if(tree == "ntKp" && variation=="" && pdf==""){ 

	// DEFINE MODEL to fit the non prompt background
		//inclusive MC signal Model
		RooRealVar* meannp = 0;
		RooRealVar* sigma1np = 0;
		RooProduct* sigma2np;
		RooRealVar* ratio_sigma12np = 0;
		RooRealVar* cofs_b_np = 0;
		meannp = new RooRealVar(Form("meannp%d_%s",_count,""),"meannp",5.281,5.24,5.32);
		sigma1np = new RooRealVar(Form("sigma1np%d_%s",_count,""),"sigma1np",0.02,0.01,0.03);
		ratio_sigma12np = new RooRealVar(Form("ratio_sigma12np%d_%s",_count,""),"ratio_sigma12np", 2.4, 0.1, 10);
		sigma2np = new RooProduct(Form("sigma2np%d_%s",_count,""), "sigma2np", RooArgList(*sigma1np, *ratio_sigma12np));
		cofs_b_np = new RooRealVar(Form("cofs_b_np%d_%s",_count,""), "cofs_b_np", 0.5, 0., 1.);
		RooGaussian* signal1_b_np = new RooGaussian(Form("signal1_b_np%d_%s",_count,""),"signal_gauss1_b_np",*mass,*meannp,*sigma1np);
		RooGaussian* signal2_b_np = new RooGaussian(Form("signal2_b_np%d_%s",_count,""),"signal_gauss2_b_np",*mass,*meannp,*sigma2np); 
		RooAddPdf* signalnp = new RooAddPdf(Form("signalnp%d_%s",_count,""), "signalnp", RooArgList(*signal1_b_np,*signal2_b_np),*cofs_b_np);
		w.import(*signalnp);
		//inclusive MC signal Model

	// MC Part. Reconstructed Background Model
	RooRealVar* m_nonprompt_scale=0;
	RooRealVar* m_nonprompt_shift=0;
	m_nonprompt_scale = new RooRealVar(Form("m_nonprompt_scale%d_%s",_count,""), "m_nonprompt_scale",0.01, 0.005, 1);
    m_nonprompt_shift = new RooRealVar(Form("m_nonprompt_shift%d_%s",_count,""), "m_nonprompt_shift", 5.15, 5.1, 5.2);
	RooGenericPdf* erfc = new RooGenericPdf(Form("erfc%d_%s",_count,""), "erfc", Form("TMath::Erfc((Bmass-m_nonprompt_shift%d_%s)/m_nonprompt_scale%d_%s)",_count,"",_count,""), RooArgList(*mass, *m_nonprompt_scale, *m_nonprompt_shift));
	// MC Part. Reconstructed Background Model

	// MC Combinatorial Background Model
	//RooRealVar* alpha_np; 
	RooRealVar* alpha_np; 
	alpha_np = new RooRealVar(Form("alpha_np%d_%s", _count,""), "alpha_np",-0.6 , -10.0 , -0.1);
	RooExponential COMB_jpsi(Form("COMB_jpsi%d_%s",_count,""), "COMB_jpsi", *mass, *alpha_np);
		
	if ( (which_var == "Bpt") && ((int)binmin == 50 && (int)binmax == 60)){
		cout << "not enough data for bin " << (int)binmin << "_" << (int)binmax << " use last parameters instead" << endl;
		alpha_np->setVal(w.var(Form("alpha_np%d_%s", _count-1,""))->getVal());
		alpha_np->setConstant();
		m_nonprompt_scale->setVal(w.var(Form("m_nonprompt_scale%d_%s",_count-1,""))->getVal());
		m_nonprompt_scale->setConstant();
		m_nonprompt_shift->setVal(w.var(Form("m_nonprompt_shift%d_%s",_count-1,""))->getVal());
		m_nonprompt_shift->setConstant();
	} 
	// MC Combinatorial Background Model

	RooRealVar jpsinp_fraction(Form("jpsinp_fraction%d_%s",_count,""), "fraction", 0.35, 0.05, 1);
	RooAddPdf* m_jpsinp_cont = new RooAddPdf(Form("m_jpsinp_cont%d_%s",_count,""), "model for jpsi nonprompt bg", RooArgList(COMB_jpsi, *erfc), RooArgList(jpsinp_fraction));
	w.import(*m_jpsinp_cont);
	// DEFINE MODEL to fit the non prompt background

	fit_jpsinp(w,  nbinsmasshisto, pdf, binmin, binmax, which_var.Data(), bsbpbins);
	}

// FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp FIT MCnp

	//FIT THE DATA FIT THE DATA FIT THE DATA
	scale->setConstant(false);	
	RooAddPdf* model;

	c->cd();
	RooPlot* frame = mass->frame("");
	TPad *p1 = new TPad("p1","p1",0.,0.215,1.,1);
	p1->SetBorderMode(1); 
	p1->SetFrameBorderMode(0); 
	p1->SetBorderSize(2);
	p1->SetBottomMargin(0.10);
	p1->Draw(); 

	TPad *p2 = new TPad("p2","p2",0.,0.,1.,0.24); 
	p2->SetTopMargin(0.);    
	p2->SetBorderMode(0);
	p2->SetBorderSize(2); 
	p2->SetFrameBorderMode(0); 
	p2->SetTicks(1,1); 
	p2->Draw();

	double n_signal_initial = ds->sumEntries(TString::Format("abs(Bmass-%g)<0.05",init_mean)) - ds->sumEntries(TString::Format("abs(Bmass-%g)<0.10&&abs(Bmass-%g)>0.05",init_mean,init_mean));
	if(n_signal_initial<0){n_signal_initial=1;}
	cout << "n_signal_initial " << n_signal_initial << endl; 

///////////////// SIGNAL FUNCTIONS

	RooRealVar nsig(Form("nsig%d_%s",_count,pdf.Data()),"",n_signal_initial,-1,ds->sumEntries()*3);
	RooRealVar mean(Form("mean%d_%s",_count,pdf.Data()),"",meanMC.getVal(),m1,m2) ;
	RooRealVar sigma1(Form("sigma1%d_%s",_count,pdf.Data()),"",sigma1MC.getVal(),0.01,0.1) ;
	RooRealVar sigma2(Form("sigma2%d_%s",_count,pdf.Data()),"",sigma2MC.getVal(),0.01,0.1) ;
	RooRealVar sigma3(Form("sigma3%d_%s",_count,pdf.Data()),"",sigma3MC.getVal(),0.01,0.1) ;
	RooRealVar sigma4cb(Form("sigma4cb%d_%s",_count,pdf.Data()),"",sigma4cbMC.getVal(),0.01,0.1) ;
	RooRealVar sigma5cb(Form("sigma5cb%d_%s",_count,pdf.Data()),"",sigma5cbMC.getVal(),0.01,0.1) ;
	RooRealVar alpha(Form("alpha%d_%s",_count,pdf.Data()),"",alphaMC.getVal(),0,50);
	RooRealVar n(Form("n_%d_%s", _count, pdf.Data()),"",nMC.getVal(),0,500);
	RooRealVar alpha2(Form("alpha2%d_%s",_count,pdf.Data()),"",alphaMC2.getVal(),0,50);
	RooRealVar n2(Form("n2_%d_%s", _count, pdf.Data()),"",nMC2.getVal(),0,500);
	RooProduct scaled_sigma1(Form("scaled_sigma1%d_%s",_count,pdf.Data()),"scaled_sigma1", RooArgList(*scale,sigma1));
	RooProduct scaled_sigma2(Form("scaled_sigma2%d_%s",_count,pdf.Data()),"scaled_sigma2", RooArgList(*scale,sigma2));
	RooProduct scaled_sigma3(Form("scaled_sigma3%d_%s",_count,pdf.Data()),"scaled_sigma3", RooArgList(*scale,sigma3));
	RooProduct scaled_sigmacb(Form("scaled_sigmacb%d_%s",_count,pdf.Data()),"scaled_sigmacb", RooArgList(*scale,sigma4cb));
	RooProduct scaled_sigmacb2(Form("scaled_sigmacb2%d_%s",_count,pdf.Data()),"scaled_sigmacb2", RooArgList(*scale,sigma5cb));
	RooGaussian sig1(Form("sig1%d_%s",_count,pdf.Data()),"",*mass,mean,scaled_sigma1);  
	RooGaussian sig2(Form("sig2%d_%s",_count,pdf.Data()),"",*mass,mean,scaled_sigma2);  
	RooGaussian sig3(Form("sig3%d_%s",_count,pdf.Data()),"",*mass,mean,scaled_sigma3);  
	RooCBShape  CB(Form("CB%d_%s",_count, pdf.Data()),"",*mass,mean,scaled_sigmacb, alpha, n);
	RooCBShape  CB2(Form("CB2%d_%s",_count, pdf.Data()),"",*mass,mean,scaled_sigmacb2, alpha2, n2);
	RooRealVar c1(Form("c1%d_%s",_count,pdf.Data()),"",1.,0.,5.);
	RooRealVar sig1frac(Form("sig1frac%d_%s",_count,pdf.Data()),"",sig1fracMC.getVal(),0.,1.);
	RooRealVar sig2frac(Form("sig2frac%d_%s",_count,pdf.Data()),"",sig2fracMC.getVal(),0.,1.);

	RooAddPdf* sig;
	if(variation=="signal" && pdf=="3gauss") sig = new RooAddPdf(Form("sig%d_%s",_count,pdf.Data()), "", RooArgList(sig1, sig2, sig3), RooArgList(sig1frac, sig2frac), true);
	if(variation=="signal" && pdf=="gauss_cb") sig = new RooAddPdf(Form("sig%d_%s",_count,pdf.Data()),"",RooArgList(sig1, CB), sig1frac);
	if((variation=="" && pdf=="") || variation== "background" || (variation=="signal" && pdf=="fixed")) sig = new RooAddPdf(Form("sig%d_%s",_count,pdf.Data()),"",RooArgList(sig1,sig2),sig1frac);
	if(variation=="signal" && pdf=="2cb") sig = new RooAddPdf(Form("sig%d_%s",_count,pdf.Data()),"",RooArgList(CB, CB2), sig1frac);	
///////////////// SIGNAL FUNCTIONS


/////////////////  BACKGROUND FUNCTIONS

	RooRealVar nbkg(Form("nbkg%d_%s",_count,pdf.Data()),"",ds->sumEntries() - n_signal_initial,0.,ds->sumEntries());
	RooRealVar a0(Form("a0%d_%s",_count,pdf.Data()),"",-0.1,-5,5);
	RooRealVar a1(Form("a1%d_%s",_count,pdf.Data()),"",0.01,-5,5);
	RooRealVar a2(Form("a2%d_%s",_count,pdf.Data()),"",1,-5e3,5e3);
	RooPolynomial bkg_1st(Form("bkg%d_%s",_count,pdf.Data()),"",*mass,RooArgSet(a0));
	RooPolynomial bkg_2nd(Form("bkg%d_%s",_count,pdf.Data()),"",*mass,RooArgSet(a0,a1));
	RooPolynomial bkg_3rd(Form("bkg%d_%s",_count,pdf.Data()),"",*mass,RooArgSet(a0,a1,a2));
	RooRealVar lambda(Form("lambda%d_%s", _count,pdf.Data()), "lambda",-1.5, -5., 1.);
	RooExponential bkg(Form("bkg%d_%s",_count,pdf.Data()),"",*mass,lambda);
	RooRealVar nbkg_part_r(Form("nbkg_part_r%d_%s",_count,pdf.Data()),"",100,0,1e5);

	// B+ PEAKING AND PART. RECONSTRUCTED BACKGROUNDS
	RooProduct *nbkg_peaking = 0;
	RooAbsPdf* jpsipi ;   
	RooAbsPdf* erfc ;   
	if(tree== "ntKp"){ 
		jpsipi = w.pdf("jpsipi");   
		erfc = w.pdf(Form("erfc%d_%s",_count,""));   
		RooRealVar* jpsipi_to_signal_ratio;
		if ( variation == "background" && pdf == "jpsi_sig"){
			jpsipi_to_signal_ratio = new RooRealVar("jpsipi_to_signal_ratio","jpsipi_to_signal_ratio", 0.0517);    // from average of pT bins
			jpsipi_to_signal_ratio->setConstant();}      
		else { jpsipi_to_signal_ratio = w.var("jpsipi_to_signal_ratio"); }
		nbkg_peaking = new RooProduct(Form("nbkg_peaking%d_%s",_count,pdf.Data()), "number of jpsi pi with fixed ratio to n_signal", RooArgList(nsig, *jpsipi_to_signal_ratio));
	}
	// B+ PEAKING AND PART. RECONSTRUCTED BACKGROUNDS
/////////////////  BACKGROUND FUNCTIONS

//////////////// MODEL MODEL MODEL MODEL
/////////////////Bs Bs Bs Bs Bs Bs Bs Bs
if(tree == "ntphi"){
	if((variation=="" && pdf=="") || (pdf=="mass_range") ||(pdf=="jpsi_sig")) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg),RooArgList(nsig,nbkg));
	if(variation=="background" && pdf=="1st") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg_1st),RooArgList(nsig,nbkg));
	if(variation=="background" && pdf=="2nd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg_2nd),RooArgList(nsig,nbkg));
	if(variation=="background" && pdf=="3rd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg_3rd),RooArgList(nsig,nbkg));
	if(variation=="signal" && pdf=="1gauss") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(sig1,bkg),RooArgList(nsig,nbkg));
	if(variation=="signal" && (pdf=="3gauss"|| pdf=="fixed"|| pdf=="gauss_cb" || pdf=="2cb")) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig, bkg),RooArgList(nsig, nbkg));
}

/////////////////BP BP BP BP BP BP BP BP
if(tree == "ntKp"){
    if((pdf=="jpsi_sig" || pdf=="")) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(bkg,*sig, *erfc,*jpsipi),RooArgList(nbkg,nsig,nbkg_part_r,*nbkg_peaking));
	if(pdf=="mass_range"){ model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig,bkg,*jpsipi),RooArgList(nsig,nbkg,*nbkg_peaking));}
    if(variation=="background" && pdf=="1st") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(bkg_1st,*sig,*erfc,*jpsipi),RooArgList(nbkg,nsig,nbkg_part_r,*nbkg_peaking));
	if(variation=="background" && pdf=="2nd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(bkg_2nd,*sig,*erfc,*jpsipi),RooArgList(nbkg,nsig,nbkg_part_r,*nbkg_peaking));
	if(variation=="background" && pdf=="3rd") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(bkg_3rd,*sig,*erfc,*jpsipi),RooArgList(nbkg,nsig,nbkg_part_r,*nbkg_peaking));
	if(variation=="signal" && pdf=="1gauss") model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(bkg,sig1,*erfc,*jpsipi),RooArgList(nbkg,nsig,nbkg_part_r,*nbkg_peaking));
	if((variation=="signal" && (pdf=="3gauss"|| pdf=="fixed"|| pdf=="gauss_cb" ))) model = new RooAddPdf(Form("model%d_%s",_count,pdf.Data()),"",RooArgList(*sig, bkg, *erfc, *jpsipi),RooArgList(nsig, nbkg, nbkg_part_r,*nbkg_peaking));
}
//////////////// MODEL MODEL MODEL MODEL

//////////////// SET PARAMETERS FROM MC FITS

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
	if(variation=="signal" && pdf=="2cb"){
		sigma4cb.setConstant();
		n.setConstant();
		alpha.setConstant();
		sigma5cb.setConstant();
		n2.setConstant();
		alpha2.setConstant();	
	}

	if(variation=="signal" && pdf=="fixed") mean.setConstant();

//////////////// SET PARAMETERS FROM MC FITS

////// ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT

  	TString fitRange = (pdf == "mass_range") ? "m_range" : "all";
	RooFitResult* fitResult = model->fitTo(*ds,Save(), Minos(), Extended(kTRUE), Range(fitRange));
	fitResult->Print("v"); 
	w.import(*model);

////// ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT ROOFIT

  p1->cd();
  ds->plotOn(frame, Name(Form("ds_cut%d", _count)), Binning(nbinsmasshisto), MarkerSize(0.5), MarkerStyle(8), MarkerColor(1), LineColor(1), LineWidth(1)); 
  model->plotOn(frame, Name(Form("model%d_%s", _count, pdf.Data())), Range(fitRange), Precision(1e-6), DrawOption("L"), LineColor(2), LineWidth(1));
  model->plotOn(frame, Name(Form("sig%d_%s", _count, pdf.Data())),  Components(*sig), Range(fitRange), Precision(1e-6), DrawOption("LF"), FillStyle(3002), FillColor(kOrange-3), LineStyle(7), LineColor(kOrange-3), LineWidth(1)); 
  if(tree== "ntKp")	{
	//TString option = (pdf == "mass_range")? "L" : "LF";
    //RooCmdArg drawRange = (pdf == "mass_range")? Range(fitRange) : RooCmdArg();
    model->plotOn(frame, RooFit::Name(Form("erfc%d_%s",_count,"")) , Components(*erfc), Range(fitRange),  NormRange(fitRange), LineColor(kGreen+3), LineStyle(9), LineWidth(2), DrawOption("L"));
	model->plotOn(frame, RooFit::Name("B->J/#psi #pi"), Components(*jpsipi), NormRange(fitRange), DrawOption("LF"), FillColor(kMagenta+1), LineStyle(1), LineColor(kMagenta+1), LineWidth(1)); 
					}
   	model->plotOn(frame, Name(Form("bkg%d_%s",_count,pdf.Data())) ,  Components(bkg), Range(fitRange), Precision(1e-6),  DrawOption("L"), LineStyle(7), LineColor(4), LineWidth(1));

	model->paramOn(frame,Layout(0.2, 0.45, 0.5), Format("NEU",AutoPrecision(2)));
	//frame->getAttText()->SetTextSize(0.035);
	//frame->getAttFill()->SetFillStyle(0);
	//frame->getAttLine()->SetLineWidth(0);
	frame->SetTitle("");
	frame->SetXTitle("");
	frame->GetYaxis()->SetTitle(TString::Format("Events / (%g MeV/c^{2})",(mass->getMax()-mass->getMin())/nbinsmasshisto*1000));
	frame->GetYaxis()->SetTitleOffset(1.7);
	frame->GetYaxis()->SetTitleSize(0.04);
	frame->GetYaxis()->SetTitleFont(42);
	frame->GetYaxis()->SetLabelFont(42);
	frame->GetYaxis()->SetLabelSize(0.04);
	frame->SetStats(0);
	double plot_min = 5.25;  
	if(tree=="ntKp") plot_min = 5.05;    
	frame->GetXaxis()->SetRangeUser(plot_min,5.5);  
	frame->GetXaxis()->SetNdivisions(-50205);
	frame->GetXaxis()->SetTitleSize(0.04);
	frame->GetXaxis()->SetTitleFont(42);
	frame->GetXaxis()->SetLabelFont(42);
	frame->GetXaxis()->SetLabelSize(0.04);
	frame->GetXaxis()->SetLabelOffset(0.01);
	frame->Draw();

	TLegend *leg = new TLegend(0.75,0.58,0.9,0.89,NULL,"brNDC"); 
	if (tree == "ntphi"){leg = new TLegend(0.75,0.65,0.9,0.89,NULL,"brNDC");}
	leg->SetBorderSize(0);
	leg->SetTextSize(0.04);
	leg->SetTextFont(42);
	leg->SetFillStyle(0);
	leg->AddEntry(frame->findObject(Form("ds_cut%d", _count)), " Data","LEP");
	leg->AddEntry(frame->findObject(Form("model%d_%s",_count,pdf.Data()))," Model","l");
	leg->AddEntry(frame->findObject(Form("sig%d_%s",_count,pdf.Data()))," Signal","f");
	leg->AddEntry(frame->findObject(Form("bkg%d_%s",_count,pdf.Data()))," Comb. Bkg.","l");
	if(tree== "ntKp"){
		leg -> AddEntry(frame->findObject("B->J/#psi #pi")," B^{+} #rightarrow J/#psi #pi^{+}","f");
		leg -> AddEntry(frame->findObject(Form("erfc%d_%s",_count,pdf.Data()))," B #rightarrow J/#psi X","l");
	}
	leg -> Draw();
	
	p2->cd();
	RooHist* pull_hist = frame->pullHist(Form("ds_cut%d",_count),Form("model%d_%s",_count,pdf.Data()));
	pull_hist->SetMarkerSize(0.5);
	RooPlot* pull_plot = mass->frame();
	(pull_plot->GetXaxis())->SetRangeUser(minhisto, maxhisto); //maxhisto
	//(pull_plot->GetXaxis())->SetRangeUser(plot_min,5.5); //maxhisto
	line_ref->plotOn(pull_plot, LineStyle(1), LineColor(2), LineWidth(1));  
	pull_plot->addPlotable(static_cast<RooPlotable*>(pull_hist),"XP");
	pull_plot->SetTitle("");
	if(tree=="ntKp")pull_plot->SetXTitle("m_{J/#psiK^{+}} [GeV/c^{2}]");
	if(tree=="ntphi")pull_plot->SetXTitle("m_{J/#psiK^{+}K^{-}} [GeV/c^{2}]");
	pull_plot->SetYTitle("Pull");
	pull_plot->GetYaxis()->SetTitleFont(42);  
	pull_plot->GetYaxis()->SetTitleSize(0.15);
	pull_plot->GetYaxis()->CenterTitle();
	pull_plot->GetYaxis()->SetLabelOffset(0.01);
	pull_plot->GetYaxis()->SetLabelSize(0.14);
	pull_plot->GetYaxis()->SetNdivisions(305);
	pull_plot->GetYaxis()->SetTitleOffset(0.45);
	pull_plot->GetXaxis()->SetTitleSize(0.15);
	pull_plot->GetXaxis()->SetTitleOffset(1.2);
	pull_plot->GetXaxis()->CenterTitle();
	pull_plot->GetXaxis()->SetLabelFont(42);
	pull_plot->GetXaxis()->SetLabelOffset(0.01);
	pull_plot->GetXaxis()->SetLabelSize(0.14);
	pull_plot->GetXaxis()->SetTickLength(0.16);
	pull_plot->GetXaxis()->SetNdivisions(-50205);
	pull_plot->Draw();

/*	Double_t yield = nsig.getVal();
	Double_t yieldErr = nsig.getError();
	TH1D* fh = (TH1D*)h->Clone("fh");
	double nexpected = model->expectedEvents(RooArgSet(*mass));
	double dataArr[nbinsmasshisto]; 
	double dataErrArr[nbinsmasshisto]; 
	double fitArr[nbinsmasshisto];
	double val = 0;
	for(int i = 0; i < nbinsmasshisto; i++){
		dataArr[i] = h->GetBinContent(i+1);
		dataErrArr[i] = h->GetBinError(i+1);
		mass->setVal(h->GetBinCenter(i+1));
		val = model->getVal(RooArgSet(*mass))*1./nbinsmasshisto*nexpected;
		fitArr[i] = val;
		fh->SetBinContent(i+1, fitArr[i]);
		fh->SetBinError(i+1, sqrt(fitArr[i])); 
	} */
	
	////////////////////////////////// SIGNIFICANCE
	/*float bkgd;
	double Significance;
	double real_significance;
	nsig.setVal(0.);
	nsig.setConstant();
	RooFitResult* fitResult_nosig = model->fitTo(*ds,Save());
	RooAbsReal* nll_nosig = model->createNLL(*ds);
	double log_likelihood_nosig= nll_nosig->getVal();
	RooAbsReal* nll = model->createNLL(*ds);
	double log_likelihood= nll->getVal();
    real_significance = sqrt(2*(-log_likelihood+log_likelihood_nosig));
	std::cout<<"REAL SIGNIFICANCE= "<<real_significance<<std::endl;
	std::cout<<"***********************************************************************"<<std::endl;
	RooAbsReal *bkgIntegral = bkg.createIntegral(*mass,NormSet(*mass),Range("signal"));
	// bkgIntegralErr = bkgIntegral->getPropagatedError(*fitResult);	
	cout<<"bkg integral: "<<bkgIntegral->getVal()<<endl;
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
	RooAbsReal* RangeBakground = bkg.createIntegral(*mass,NormSet(*mass),Range("signal")); 
	double Calback = RangeBakground->getVal() * nbkg.getVal();
	double StatSig = yield/sqrt(yield + Calback);
	TLatex *lat = new TLatex();
	lat->SetNDC();
	lat->SetTextSize(0.045);
	lat->DrawLatex(0.63,0.55,Form("S = %.1f",yield));	
	lat->DrawLatex(0.63,0.50,Form("B = %.1f",Calback));	
	lat->DrawLatex(0.63,0.45,Form("S/#sqrt{S+B} = %.1f",StatSig));	
	lat->Draw("SAME"); */
	////////////////////////////////// SIGNIFICANCE

	Double_t yieldPrintErr = nsig.getError();
	Double_t yieldPrintErrUp = nsig.getAsymErrorHi();
	Double_t yieldPrintErrDown = nsig.getAsymErrorLo();

	cout << "----------------------------------------------------------------------------------------------------------------" << endl;
	cout << "Signal Yield = " << nsig.getVal() << "     yield Error = " << yieldPrintErr << "     yield Error Up = " << yieldPrintErrUp << "     yieldPrintErrDown = " << yieldPrintErrDown << endl;
	cout << "----------------------------------------------------------------------------------------------------------------" << endl;

	p1->cd();
	return fitResult;

} // END OF MAIN FITTING FUNCTION











// FIT TO BINED PAR_R_BKG
void fit_jpsinp(RooWorkspace& w, int nbin_hist, TString pdf, float bin_i, float bin_f, TString Var, TString bsbp) {

	RooDataSet* d_s = (RooDataSet*) w.data("jpsinp");
	cout << Var.Data() << " bins study of PART_R_BKG" << endl;  

	if(Var == "Bpt" ){
	  	d_s = (RooDataSet*) d_s->reduce(Form("(Bpt > %d && Bpt < %d)", (int) bin_i , (int) bin_f) );

	} else if(Var == "By") { 		
	  	d_s = (RooDataSet*) d_s->reduce(Form("abs(By)>%f && abs(By)< %f", bin_i , bin_f) );
		
	} else if(Var == "nMult"){ 
		cout << "for now MULTIPLICITY uses the INCLUSIVE bin";

	}
  	
	// Get rid of B+ at gen level
	RooDataSet* ds_cont = (RooDataSet*) d_s->reduce("Bgen != 23333 && Bgen != 23335 && Bgen > 5000");
	// Signal MC inclusive
	RooDataSet* ds_sig = (RooDataSet*) d_s->reduce("Bgen == 23333");

	// Create the necessary folders and define paths
	TString path_to_save = Form("./results/BP/PAR_R_bkg/%s", bsbp.Data() );
  	gSystem->mkdir(path_to_save, true);

	TString jpsi_fit_plot = "" ;
	TString incSIG_fit_plot = "" ;
	TString jpsi_plot_with_sig = "" ;
	if(Var == "Bpt" ){
		jpsi_fit_plot = path_to_save + TString::Format("/np_fit_%i-%i_%s.pdf", (int) bin_i, (int) bin_f, Var.Data() );
		incSIG_fit_plot = path_to_save + TString::Format("/InclusiveMC_Signal_fit_%i-%i_%s.pdf",(int) bin_i,(int) bin_f, Var.Data() );
  		jpsi_plot_with_sig = path_to_save + TString::Format("/np_fit_signal_%i-%i_%s.pdf",(int) bin_i,(int) bin_f, Var.Data() );
	} else if (Var == "By" ){
		jpsi_fit_plot = path_to_save + TString::Format("/np_fit_%0.1f-%0.1f_%s.pdf", bin_i, bin_f, Var.Data() );
		incSIG_fit_plot = path_to_save + TString::Format("/InclusiveMC_Signal_fit_%0.1f-%0.1f_%s.pdf", bin_i, bin_f, Var.Data() );
  		jpsi_plot_with_sig = path_to_save + TString::Format("/np_fit_signal_%0.1f-%0.1f_%s.pdf", bin_i, bin_f, Var.Data() );
	} else if(Var == "nMult" ){
		jpsi_fit_plot = path_to_save + TString::Format("/np_fit_%s.pdf", Var.Data());
		incSIG_fit_plot = path_to_save + TString::Format("/InclusiveMC_Signal_fit_%s.pdf", Var.Data());
  		jpsi_plot_with_sig = path_to_save + TString::Format("/np_fit_signal_%s.pdf", Var.Data());
	}

	//[START] FIX SHAPE (NP background)
	RooRealVar Bmass = *(w.var("Bmass"));
	RooAbsPdf* m_jpsinp_cont = w.pdf(Form("m_jpsinp_cont%d_%s",_count, pdf.Data()));
	// FIT
	auto cont_result = m_jpsinp_cont->fitTo(*ds_cont, Save());
	// FIT
	plot_jpsifit(w, nbin_hist, pdf, m_jpsinp_cont, ds_cont, jpsi_fit_plot, false);
	fix_parameters(w, Form("m_jpsinp_cont%d_%s",_count,pdf.Data()));
	//[END] FIX SHAPE (NP background) 

	//[START] FIT inclusiveMC SIGNAL 
	RooAbsPdf* signalnp = w.pdf(Form("signalnp%d_%s",_count,pdf.Data()));  
  	RooRealVar n_signal_np(Form("n_signal_np%d_%s",_count,pdf.Data()), "n_signal_np", 1000, 0., 150000); 
	RooExtendPdf signal_ext(Form("signal_ext%d_%s",_count,pdf.Data()), "extended signal pdf", *signalnp, n_signal_np);
	Bmass.setRange("bmc", 5.15, 5.4);
	// FIT
	auto signal_result = signal_ext.fitTo(*ds_sig, Range("bmc"), Save(), Extended());
	// FIT
	plot_mcfit(w, &signal_ext, ds_sig, incSIG_fit_plot,  Range("bmc"), LineColor(kRed),LineStyle(1), LineWidth(2));		
	fix_parameters(w, Form("signalnp%d_%s",_count,pdf.Data()));	
	//[END] FIT inclusiveMC SIGNAL 

	//FIX RATIO
		// Fix the ratio of jpsipi to signal
		RooRealVar jpsipi_to_signal_ratio("jpsipi_to_signal_ratio", "jpsipi_to_signal_ratio",0.05, 0, 1);
		jpsipi_to_signal_ratio.setVal(0.0384);   		// from PDG 	
  		jpsipi_to_signal_ratio.setConstant();
		w.import(jpsipi_to_signal_ratio);
	//FIX RATIO

	// TEST THE FIT TEplotNameST THE FIT
	// Import from WorkSpace and define variables
	RooAbsPdf* erfc = w.pdf(Form("erfc%d_%s",_count,pdf.Data()));
	RooAbsPdf* COMB_jpsi = w.pdf(Form("COMB_jpsi%d_%s",_count,pdf.Data()));  
	RooAbsPdf* jpsipi = w.pdf("jpsipi");    
  	RooRealVar* sigma1_np = w.var(Form("sigma1np%d_%s",_count,pdf.Data()));      
	RooRealVar n_cont(Form("n_cont_np%d_%s",_count,pdf.Data()), "n_cont_np", 1000, 0., (d_s->sumEntries())*2);
	RooRealVar n_erfc(Form("n_nonprompt%d_%s",_count,pdf.Data()), "n_nonprompt", 1000, 0., (d_s->sumEntries())*2);
	RooProduct* n_jpsipi = new RooProduct(Form("n_jpsipi_by_signal%d_%s",_count,pdf.Data()), "number of jpsis pi with fixed ratio to n_signal", RooArgList(n_signal_np, jpsipi_to_signal_ratio));
	RooRealVar* alpha_comb_m =  w.var(Form("alpha_np%d_%s", _count,""));

	// Unfix the signal and back to better describe each pT bin peak
  	sigma1_np->setConstant(false);
	alpha_comb_m->setConstant(false);

	// BUILD TOTAL PDF
	RooAddPdf* model_inclusive = new RooAddPdf(Form("model_inclusive%d_%s",_count,pdf.Data()), "NP with B+", RooArgList(*signalnp, *jpsipi, *erfc, *COMB_jpsi), RooArgList(n_signal_np, *n_jpsipi, n_erfc, n_cont));
    // FITFITFIT JUST TO CHECK 
    model_inclusive->fitTo(*d_s, Save(), Extended(), NumCPU(4));
	// Plot
    plot_jpsifit(w, nbin_hist, pdf, model_inclusive, d_s, jpsi_plot_with_sig, true);
	// TEST THE FIT TEST THE FIT
}
// FIT TO BINED PAR_R_BKG











void plot_jpsifit(RooWorkspace& w, int nbin_hist, TString pdf, RooAbsPdf* model, RooDataSet* ds, TString plotName, bool with_sig) {
  
  TCanvas* can_np= new TCanvas("can_mc","",600,600);
  can_np->cd();
  TPad *p1 = new TPad("p1","p1",0.,0.2,1.,0.99);
  p1->SetBorderMode(1); 
  p1->SetFrameBorderMode(0); 
  p1->SetBorderSize(2);
  p1->SetBottomMargin(0.10);
  p1->Draw(); 
  p1->cd();
  
  RooRealVar Bmass = *(w.var("Bmass"));
  Bmass.setRange("bmass", 5.0, 6.0);
  RooPlot* massframe = Bmass.frame(Title(" "));
  massframe->GetYaxis()->SetTitle(TString::Format("Events / (%g MeV/c^{2})",(Bmass.getMax()-Bmass.getMin())/100*1000));
  massframe->GetXaxis()->SetTitle("m_{J/#psi #pi^{+}} [GeV/c^{2}]");
  massframe->GetXaxis()->CenterTitle();
  massframe->GetYaxis()->SetTitleOffset(1.5);
  massframe->GetXaxis()->SetTitleOffset(1.2);
  massframe->GetYaxis()->SetTitleSize(0.035);
  massframe->GetYaxis()->SetTitleFont(42);
  massframe->GetXaxis()->SetLabelFont(42);
  massframe->GetYaxis()->SetLabelFont(42);
  massframe->GetXaxis()->SetLabelSize(0.035);
  massframe->GetYaxis()->SetLabelSize(0.035);	
  massframe->GetXaxis()->SetRangeUser(5.0,5.6);
  ds->plotOn(massframe, RooFit::Name("NP"),Binning(nbin_hist), MarkerSize(0.5),MarkerStyle(8),LineColor(1),LineWidth(1));
  model->plotOn(massframe, RooFit::Name("NP Fit"), NormRange("bmass"),LineColor(kRed), LineStyle(1), LineWidth(1));
  model->plotOn(massframe, RooFit::Name("par"),Components(Form("erfc%d_%s",_count,pdf.Data())), NormRange("bmass"), LineColor(kGreen+3), LineStyle(9), LineWidth(2), DrawOption("L"));
  model->plotOn(massframe, RooFit::Name("COMB_jpsi"),Components(Form("COMB_jpsi%d_%s",_count,pdf.Data())), NormRange("bmass"),LineColor(kBlue), LineWidth(1), LineStyle(kDashed));
  if (with_sig) {
   model->plotOn(massframe, RooFit::Name("signal"),Components(Form("signalnp%d_%s",_count,pdf.Data())), NormRange("bmass"), LineColor(kOrange-3), LineStyle(1), LineWidth(1), FillStyle(3002),FillColor(kOrange-3), VLines(), DrawOption("LF"));
   model->plotOn(massframe, RooFit::Name("B->J/#psi #pi"),Components("jpsipi"), NormRange("bmass"),DrawOption("LF"), FillColor(kMagenta+1), LineStyle(1), LineColor(kMagenta+1), LineWidth(1)); 
  }

  model->paramOn(massframe,  Layout(0.6, 0.95, 0.65),Format("NEU", AutoPrecision(1)));
  massframe->getAttText()->SetTextSize(0.025);
  massframe->getAttFill()->SetFillStyle(0);
  massframe->getAttLine()->SetLineWidth(0);
  massframe->Draw();

  TLatex txt;
  TLegend *leg = new TLegend(0.75,0.71,0.89,0.89,NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.025);     
  leg->SetTextFont(42);
  leg->SetFillStyle(0);
  leg->AddEntry(massframe->findObject("NP"), " MC", "ple");
  leg->AddEntry(massframe->findObject("par"), " Part. Rec. Bkg.", "l");
  leg->AddEntry(massframe->findObject("NP Fit"),"Part. + Comb Bkg.","l");
  // compare yields with gen particles
  if (with_sig) {
    leg->AddEntry(massframe->findObject("signal"), "signal", "f");
    leg->AddEntry(massframe->findObject("COMB_jpsi"), " Comb. Bkg.", "l");
    leg->AddEntry(massframe->findObject("B->J/#psi #pi"), "B->J/#psi #pi", "f");
  				}
  leg->Draw();
  can_np->SaveAs(plotName);
				}


template<typename... Targs>
void plot_mcfit(RooWorkspace& w, RooAbsPdf* model, RooDataSet* ds, TString plotName, Targs... options) {
  
  TCanvas* can_mc= new TCanvas("can_mc","",700,700);
  can_mc->cd();
  TPad *p1 = new TPad("p1","p1",0.,0.215,1.,1);
  p1->SetBorderMode(1); 
  p1->SetFrameBorderMode(0); 
  p1->SetBorderSize(2);
  p1->SetBottomMargin(0.10);
  p1->Draw(); 
  p1->cd();

  RooRealVar Bmass = *(w.var("Bmass"));
  RooPlot* massframe = Bmass.frame(Title(" "));
  massframe->GetYaxis()->SetTitle(TString::Format("Events / (%g MeV/c^{2})",(Bmass.getMax()-Bmass.getMin())/100*1000));
  massframe->GetXaxis()->SetTitle("m_{J/#psi #pi^{+}} [GeV/c^{2}]");
  massframe->GetXaxis()->CenterTitle();
  massframe->GetYaxis()->SetTitleOffset(1.5);
  massframe->GetXaxis()->SetTitleOffset(1.2);
  massframe->GetYaxis()->SetTitleSize(0.035);
  massframe->GetYaxis()->SetTitleFont(42);
  massframe->GetXaxis()->SetLabelFont(42);
  massframe->GetYaxis()->SetLabelFont(42);
  massframe->GetXaxis()->SetLabelSize(0.035);
  massframe->GetYaxis()->SetLabelSize(0.035);	
  massframe->GetXaxis()->SetRangeUser(5.15,5.45);
  if (plotName == "./results/BP/InclusiveMC_JPsipi_fit.pdf") {massframe->GetXaxis()->SetRangeUser(5,6);}
  else if (plotName == "./results/BP/InclusiveMC_JPsipi_fit_BsBPBINS.pdf") {massframe->GetXaxis()->SetRangeUser(5,6);}

  ds->plotOn(massframe, Name("dsjpp") , MarkerSize(0.5), MarkerStyle(8), LineColor(1), LineWidth(1));
  model->plotOn(massframe, RooFit::Name("MCFit"), options...);
  model->paramOn(massframe, Layout(0.62, 0.89, 0.75), "", Format("NEU", AutoPrecision(1)) ) ;
  massframe->getAttText()->SetTextSize(0.025);
  massframe->getAttFill()->SetFillStyle(0);
  massframe->getAttLine()->SetLineWidth(0);
  massframe->Draw();
  TLegend *leg_jpp = new TLegend(0.75,0.8,0.89,0.89,NULL,"brNDC");
  leg_jpp->SetBorderSize(0);
  leg_jpp->SetTextSize(0.025);     
  leg_jpp->SetTextFont(42);
  leg_jpp->SetFillStyle(0);
  leg_jpp->AddEntry(massframe->findObject("dsjpp"), " MC","p");
  leg_jpp->AddEntry(massframe->findObject("MCFit")," Peaking Bkg. PDF","f");
  leg_jpp->Draw();
  can_mc->SaveAs(plotName);
}


/* Fix or release the parameters for a given PDF */
void fix_parameters(RooWorkspace& w, TString pdfName, bool release) {
  RooAbsPdf* pdf = w.pdf(pdfName);
  RooAbsData* ds = w.data("jpsinp");
  RooArgSet* par_set = pdf->getParameters(*ds);
  auto itr = par_set->createIterator();
  bool toFix = ! release;
  std::string fix_or_float = (toFix)? "fix " : "float ";
  std::cout << fix_or_float << "parameters:";
  for (auto i = 0; i < par_set->getSize(); ++i) {
    RooRealVar* var = (RooRealVar*) itr->Next();
    var->setConstant(true);
    TString name = var->GetName();
    std::cout << name << ", ";
    var->setConstant(toFix);
  }
  std::cout << "\n";
}

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

	file_check << "\\begin{document}" << std::endl;
	// Create table                                                                                                                                                                                                                                                
	std::string col="c";
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
			file << /*std::setprecision(3)<<*/  numbers[c-1][i-1]<< " \\% & ";
			file_check << /*std::setprecision(3)<<*/  numbers[c-1][i-1]<< " \\% & ";
									}
		file << /*std::setprecision(3)<<*/  numbers[n_col-2][i-1]<< " \\% \\\\" << std::endl;
		file_check << /*std::setprecision(3)<<*/  numbers[n_col-2][i-1]<< " \\% \\\\" << std::endl; 
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
	system(("pdflatex " + filename+ "_check.tex").c_str());
	//system(("open " + filename + "_check.pdf").c_str());
}

void validate_fit(RooWorkspace* w, TString pdf, TString tree, TString variable, int full, float binmin, float binmax, string Path)
{
	std::cout << "Performing Check on Fit" << std::endl;
	RooRealVar Bmass = *(w->var("Bmass"));
	RooAbsPdf* model  = w->pdf(Form("model%d_%s",_count,pdf.Data()));
	//RooDataSet* data = (RooDataSet*) w->data("data");

	//model->fitTo(*data);
	vector<RooRealVar> params;
	params.push_back(*(w->var(Form("nsig%d_%s",_count,pdf.Data()))));
	//params.push_back(*(w->var(Form("mean%d_%s",_count,pdf.Data()))));

	/*
	//Retrieving Parameters//

	RooRealVar* lambda =  (w->var(Form("lambda%d", _count)));
	RooRealVar* nbkg   =  (w->var(Form("nbkg%d", _count)));
	RooRealVar * mean  = (w->var(Form("mean%d", _count)));
	RooRealVar * sigma1  = (w->var(Form("sigma1%d", _count)));
	RooRealVar * sigma2  = (w->var(Form("sigma2%d", _count)));
	//	RooRealVar * c1  = (w->var(Form("c1%d", _count)));
	RooRealVar * sig1frac = (w->var(Form("sig1frac%d", _count)));

	lambda->setVal(lambda->getVal());
	nbkg->setVal(nbkg->getVal());
	mean->setVal(mean->getVal());
	sigma1->setVal(sigma1->getVal());
	sigma2->setVal(sigma2->getVal());
	//	c1->setVal(c1->getVal());
	sig1frac->setVal(sig1frac->getVal());
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
	int params_size = params.size();

	RooMCStudy* mcstudy = new RooMCStudy(*model, Bmass,  Extended(), FitOptions(Save(kTRUE), PrintEvalErrors(0)));
	mcstudy->generateAndFit(5000);

	cout << "DONE Generate and Fit " << endl;

	TString XName[2] = {"P(Y)","mean"};
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


	TCanvas* c_pull = new TCanvas("pulls", "pulls", 700, 700);
	TCanvas* c_params = new TCanvas("params", "params", 700, 700);
	TCanvas* c_errors = new TCanvas("errors", "errors", 700, 700);

	gStyle->SetOptFit(0111);
	gPad->SetLeftMargin(0.15);
	gStyle->SetStatX(0.95);		//Stat box x position (top right hand corner)	
	gStyle->SetStatY(0.9); 		//Stat box y position 	
	gStyle->SetStatW(0.1);	 		//Stat box width as fraction of pad size	0.05	
	gStyle->SetStatFont(62);  		//Stat box font
	gStyle->SetStatFontSize(0);
	gStyle->SetStatBorderSize(0);
	gStyle->SetStatColor(0);
	gStyle->SetStatStyle(0);		//Stat box fill style hollow

	for(int i = 0; i < params_size; ++i){

		n_signal_init = params[i].getVal();
		n_signal_error_init = params[i].getError();

		c_pull->cd();
		h1[i]->SetTitle("");
		h1[i]->Fit("gaus","","",-5,5);
		
		c_pull->Update();
		h1[i]->GetFunction("gaus")->SetLineColor(kCyan+1);
		h1[i]->GetFunction("gaus")->SetLineWidth(1);
		h1[i]->GetFunction("gaus")->SetFillStyle(3019);
		h1[i]->GetFunction("gaus")->SetFillColor(kCyan+1);
		h1[i]->GetXaxis()->SetTitle(Form("%s",XName[i].Data()));
		h1[i]->GetYaxis()->SetTitle("Toy MCs");
		h1[i]->GetXaxis()->SetRangeUser(-5, 5);
		h1[i]->Draw("APsame");
	
		c_errors->cd();
		h3[i]->SetTitle("");
		h3[i]->Fit("gaus","","",n_signal_error_init*0.5,n_signal_error_init*1.5);
		
		c_errors->Update();
		h3[i]->GetFunction("gaus")->SetLineColor(4);
		h3[i]->GetFunction("gaus")->SetLineWidth(1);
		h3[i]->GetXaxis()->SetTitle(Form("%s Error",XName[i].Data()));
		h3[i]->GetYaxis()->SetTitle("");
		h3[i]->Draw("APsame");

		c_params->cd();
		h2[i]->SetTitle("");
		h2[i]->Fit("gaus","","",n_signal_init * 0.5,n_signal_init * 1.5);

		c_params->Update();
		h2[i]->GetFunction("gaus")->SetLineColor(4);
		h2[i]->GetFunction("gaus")->SetLineWidth(1);
		h2[i]->GetXaxis()->SetTitle(Form("%s Mean",XName[i].Data()));
		h2[i]->GetYaxis()->SetTitle("");
		h2[i]->Draw("APsame");

		c_pull->SaveAs(Form("%s/pull_signal_%s_%0.1f_%0.1f_%d_%s.pdf",Path.c_str(),variable.Data(),binmin,binmax,i,tree.Data()));
		c_params->SaveAs(Form("%s/param_signal_%s_%0.1f_%0.1f_%d_%s.pdf",Path.c_str(),variable.Data(),binmin,binmax,i,tree.Data()));
		c_errors->SaveAs(Form("%s/error_signal_%s_%0.1f_%0.1f_%d_%s.pdf",Path.c_str(),variable.Data(),binmin,binmax,i,tree.Data()));
		
	}
}


