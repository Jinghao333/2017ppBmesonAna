#include "roofitB.h"
#include "CMS_lumi.C"
#include "parameter.h"                             //UPDATED
#include "parametersNew.h"
#include "TSystem.h"
#include <string>
#include <sstream>
#include <TGraph.h>
#include "TMultiGraph.h"
#include "TGraphErrors.h"

int syst=1;

TTree* makeTTree(TTree* intree, TString treeTitle) 
{
	TTree* outtree = new TTree(treeTitle.Data(),treeTitle.Data());
	float mass;
	outtree->Branch("mass",&mass,"mass/F") ;
	int Bsize;
	float Bmass[20000];
	intree->SetBranchAddress("Bsize",&Bsize);
	intree->SetBranchAddress("Bmass",Bmass);
	int nentries = intree->GetEntries();
	for(int n=0; n<nentries; n++){
		intree->GetEntry(n);
		for(int b=0; b<Bsize; b++){
			mass = Bmass[b];
			outtree->Fill();
		}
	}
	return outtree ;
}

//void roofitB(int usePbPb = 0, int fitOnSaved = 0, TString inputmc = "", TString varExp = "", TString trgselection = "",  TString cut = "", TString cutmcgen = "", int isMC = 0, Double_t luminosity = 1., int doweight = 1, TString collsyst = "", TString outputfile = "", TString outplotf = "", TString npfit = "", int doDataCor = 0, Float_t centmin = 0., Float_t centmax = 100.)

void roofitB(int doubly = 0, TString tree = "ntphi", int full = 1, int usePbPb = 0, int fitOnSaved = 0, TString inputdata = "", TString inputmc = "", TString varExp = "", TString trgselection = "",  TString cut = "", TString cutmcgen = "", int isMC = 0, Double_t luminosity = 1., int doweight = 1, TString outputfile = "", TString outplotf = "", TString npfit = "", int doDataCor = 0){ 

	double MyBackground;
	double	yieldRec;
	int _nBins;

// Bpt Bpt Bpt Bpt	
	if(varExp == "Bpt"){ 
	if(full == 1){ _nBins = 1 ;}
	else if(full == 0) {
		if(tree=="ntphi"){ _nBins = nptBins;}
//		if(tree=="ntphi"){ _nBins = nptBins_test;}

		else if(tree=="ntKp"){ _nBins = nptBinsBP;}}
	} else if(varExp == "By"){
		if(full == 1){_nBins = 1;}
		else if(full == 0){_nBins = nyBins_both;}
	}
		else if(varExp == "nMult"){
		if(full == 1){_nBins = 1;}
		else if(full == 0){_nBins = nmBins_both;}
	}

	cout << "number of bins" << _nBins << endl;	
	double _ptBins[_nBins+1];

	if(varExp == "Bpt"){ 
	if(full == 1){
                if(tree=="ntphi"){for( int c=0; c<_nBins+1; c++){_ptBins[c]=ptBins_full[c];}}
								else if(tree=="ntKp"){for( int c=0; c<_nBins+1; c++){_ptBins[c]=ptBins_fullBP[c];}}
	}else if(full == 0) {
		if(tree=="ntphi"){for( int c=0; c<_nBins+1; c++){_ptBins[c]=ptbinsvec[c];}}
//		if(tree=="ntphi"){for( int c=0; c<_nBins+1; c++){_ptBins[c]=ptbinsvec_test[c];}}
		else if(tree=="ntKp"){for( int c=0; c<_nBins+1; c++){_ptBins[c]=ptbinsvecBP[c];}}
	                     }
	} else if(varExp == "By"){
		if(full == 0){
			for(int c=0; c<_nBins+1; c++){_ptBins[c]=ybinsvec[c];}
		}
		else if(full == 1){
			for(int c=0; c<_nBins+1; c++){_ptBins[c]=ybinsvec_full[c];}
		}
	}
		else if(varExp == "nMult"){
		if(full == 0){
			for(int c=0; c<_nBins+1; c++){_ptBins[c]=nmbinsvec[c];}
		}
		else if(full == 1){
			for(int c=0; c<_nBins+1; c++){_ptBins[c]=nmbinsvec_full[c];}
		}
	}

std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@"<<_ptBins<<std::endl;
cout << endl << endl;
std::cout<<"Variable "<<varExp<<std::endl;
cout << tree << " BINS: ";
for(int t; t< sizeof(_ptBins)/sizeof(_ptBins[0]);t++){cout <<"__"<<  _ptBins[t]<<"__";}
cout << endl << endl;

/*	if(varExp=="nMult"){
		if(full==0) {
			nBins_check = nBins_Mult;
			_ptBins = MultBins;
			_nBins = nBins_Mult;
			_ptBins = MultBins;
		}
		else if(full==1){
			nBins_check = nBins_full;
			_ptBins = nMults_full;
		}
	}
	else if(varExp=="By"){
		if(full==0) {
			nBins_check = nBinsY;
			_ptBins = ptBinsY;
		}
		else if(full==1){
			nBins_check = nBins_full;
			_ptBins = hiBins_full;
		}

	}*/

	/*	if(varExp == "By"){
		_nBins = nBinsY;
		_ptBins = ptBinsY;
		}*/
	

	if (!(usePbPb==1||usePbPb==0)) std::cout<<"ERROR!!, you are using a non valid isPbPb option"<<std::endl;
	bool isPbPb=(bool)(usePbPb);

	std::cout<<"DEBUG 3"<<std::endl;
	if(!isPbPb)
	{
		std::cout<<"DEBUG 4 !PbPb"<<std::endl;
		seldata = Form("%s",trgselection.Data());
		std::cout<<"DEBUG 5"<<std::endl;
		selmc = Form("%s",cut.Data());
		std::cout<<"DEBUG 6"<<std::endl;
		selmcgen = Form("%s",cutmcgen.Data());
		std::cout<<"DEBUG 7"<<std::endl;
	}
	else
	{
		std::cout<<"DEBUG 4 else"<<std::endl;
		seldata = Form("%s&&%s",trgselection.Data(),cut.Data());
		std::cout<<"DEBUG 5"<<std::endl;
		selmc = Form("%s",cut.Data());
		std::cout<<"DEBUG 6"<<std::endl;
		selmcgen = Form("%s",cutmcgen.Data());
		std::cout<<"DEBUG 7"<<std::endl;
		std::cout<<"cut = "<<cut<<std::endl;
		std::cout<<"seldata= "<<seldata<<std::endl;
		std::cout<<"selmc= "<<selmc<<std::endl;
	}
	//return;
	gStyle->SetTextSize(0.05);
	gStyle->SetTextFont(42);
	gStyle->SetPadRightMargin(cRightMargin);
	gStyle->SetPadLeftMargin(cLeftMargin);
	gStyle->SetPadTopMargin(cTopMargin);
	gStyle->SetPadBottomMargin(cBottomMargin);
	cout << "cBottomMargin = " << cBottomMargin << endl;
	gStyle->SetPadBottomMargin(0.45);
	gStyle->SetTitleX(.0f);

	TFile* inf = new TFile(inputdata.Data());
	TTree* skimtree_new = (TTree*)inf->Get(tree);
	//TTree* skimtree_new =skimtree->CloneTree();
	//  inf->Close();
	TFile* infMC = new TFile(inputmc.Data());
	std::cout<<"MC file: "<<inputmc.Data()<<std::endl;
	//  infMC->cd();
	std::cout<<"cd"<<std::endl;
	TTree* skimtreeMC_new = (TTree*)infMC->Get(tree);
	//TTree* skimtreeMC_new =skimtreeMC->CloneTree();
	//  fout->cd();
	std::cout<<"depois do get"<<std::endl;
	std::cout<<"clonou"<<std::endl;
	//  inf->Close();
	//  fout->Close();

	TH1D* h;
	TH1D* hMC;
	TH1D* hpull;

	//TTree* ntGen = new TTree();
	/*TTree* skimtree= new TTree();
	  TTree* skimtreeMC = new TTree();
	  TTree* skimtree_new= new TTree();
	  TTree* skimtreeMC_new = new TTree();*/

	//if(fitOnSaved == 0){

//	std::cout<<"Creating trees "<<std::endl;

	/*skimtree->SetObject("ntphi_data", "ntphi_data");  
	  skimtree_new->SetObject("ntphi_data", "ntphi_data");  */
	//nt->AddFriend("hltanalysis/HltTree");
	//nt->AddFriend("hiEvtAnalyzer/HiTree");
	//nt->AddFriend("skimanalysis/HltTree");
	/*
	   nt->AddFriend("BDT_pt_15_20");
	   nt->AddFriend("BDT_pt_7_15");
	   nt->AddFriend("BDT_pt_5_7");
	   nt->AddFriend("BDT_pt_20_50");	
	   */
	//nt->AddFriend("BDT");	

	/*
	   ntGen = (TTree*)infMC->Get("ntGen");
	   ntGen->AddFriend("ntHlt");
	   ntGen->AddFriend("ntHi");
	   cout << "Pass 2" << endl;
	   ntMC = (TTree*)infMC->Get("ntphi");
	   ntMC->AddFriend("ntHlt");
	   ntMC->AddFriend("ntHi");
	   ntMC->AddFriend("ntSkim");
	   ntMC->AddFriend("mvaTree");
	   ntMC->AddFriend("ntGen");
	   */
	/*	ntGen = (TTree*)infMC->Get("Bfinder/ntGen");
		ntGen->AddFriend("hltanalysis/HltTree");
		ntGen->AddFriend("hiEvtAnalyzer/HiTree");*/
	/*  infMC->cd();
		std::cout<<"cd"<<std::endl;
		skimtreeMC = (TTree*)infMC->Get("ntphi");
		std::cout<<"depois do get"<<std::endl;
		skimtreeMC_new =skimtreeMC->CloneTree();
		std::cout<<"clonou"<<std::endl;
		TFile* fout = new TFile("/afs/cern.ch/work/j/jusilva/CMSSW_7_5_8_patch5/src/UserCode/Bmeson2018Ana/Bs/CrossSection/test_trees/test_mc_new.root", "RECREATE");
		fout->cd();
		skimtreeMC->Write();
		fout->Close();*/
	/*		ntMC->AddFriend("hltanalysis/HltTree");
			ntMC->AddFriend("hiEvtAnalyzer/HiTree");
			ntMC->AddFriend("skimanalysis/HltTree");*/
	/*
	   ntMC->AddFriend("BDT_pt_15_20");
	   ntMC->AddFriend("BDT_pt_7_15");
	   ntMC->AddFriend("BDT_pt_5_7");
	   ntMC->AddFriend("BDT_pt_20_50");	
	   */
	//		ntMC->AddFriend("BDT");
	//	ntMC->AddFriend("Bfinder/ntGen");	
	//}



	TH1D* hPt = new TH1D("hPt","",_nBins,_ptBins);

	TH1D* hPtMC = new TH1D("hPtMC","",_nBins,_ptBins);
	//TH1D* hPtGen = new TH1D("hPtGen","",_nBins,_ptBins);
	TH1D* hMean = new TH1D("hMean","",_nBins,_ptBins);                       
	TH1D* hSigmaGaus1 = new TH1D("hSigmaGaus1","",_nBins,_ptBins); 
	TH1D* hSigmaGaus2 = new TH1D("hSigmaGaus2","",_nBins,_ptBins); 
	std::cout<<"Histograms"<<std::endl;
	//RooWorkspace* w = new RooWorkspace("w");
	RooRealVar* mass = new RooRealVar("Bmass","Bmass",minhisto,maxhisto);
	RooRealVar* pt = new RooRealVar("Bpt","Bpt",0,100);
	RooRealVar* y = new RooRealVar("By","By",-2.4, 2.4);
//	RooRealVar* HiBin = new RooRealVar("hiBinNew","hiBinNew",0.,90.*2);
	RooRealVar* nMult = new RooRealVar("nMult","nMult",0,200);
	RooRealVar* trackSelection = new RooRealVar("track", "track", 0, 5);
//	RooRealVar* w = new RooRealVar("Pthatweight","Pthatweight",0.,12.) ;
	RooDataSet* ds = new RooDataSet();
	RooDataSet* dsMC = new RooDataSet();   
	std::cout<<"Created dataset"<<std::endl;
	RooDataHist* dh = new RooDataHist();   
	RooDataHist* dhMC = new RooDataHist();   
	std::cout<<"Created roodatahists"<<std::endl;
	RooPlot* frame = new RooPlot();
	std::cout<<"Created d"<<std::endl;
	RooHist* datahist = new RooHist();
	RooCurve* modelcurve = new RooCurve();
	std::cout<<"Created roocurve"<<std::endl;

	//weightgen = weightgen_pp;
	weightmc  = weightmc_pp;
	if(usePbPb){
		weightgen = weightgen_PbPb;
		//weightmc = weightmc_PbPb;
		if(doweight == 0 ) weightmc = "1"; 
		if(doweight == 1) weightmc = "pthatweightNew";
	    weightmc = "1"; 
		
		std::cout<<"weights defined"<<std::endl;	

	}

	TString _prefix = "";
	TString _isMC = "data";
	if(isMC) _isMC = "mcAsData";
	TString _isPbPb = "pp";
	//if(usePbPb) _isPbPb ="PbPb"
	TString _postfix = "";
	if(weightdata!="1") _postfix = "_EFFCOR";


	std::cout<<"Passed AddFriends"<<std::endl;
	TString outputf;
	outputf = Form("%s",outputfile.Data());
	std::cout<<"Formed outputf "<<outputf<<std::endl;
	TFile* outf = new TFile(outputf.Data(),"recreate");
	outf->cd();
	

	ds = new RooDataSet(Form("ds%d",_count),"",skimtree_new,
                      RooArgSet(*mass, *pt, *y, *nMult, *trackSelection));
	dsMC = new RooDataSet(Form("dsMC%d",_count),"",skimtreeMC_new,
                        RooArgSet(*mass, *pt, *y, *nMult, *trackSelection));

	//    TString ptbinning;

	std::vector<std::string> background = {"1st", "2nd","mass_range"};
	std::vector<std::string> signal = {"3gauss", "fixed", "gauss_cb"};
	
	std::vector<std::vector<double>> background_syst;
	std::vector<std::vector<double>> signal_syst;
	std::vector<std::vector<double>> general_syst;
	std::vector<double> nominal_yields;
	std::vector<std::vector<double>> back_syst_rel_values;
	std::vector<std::vector<double>> sig_syst_rel_values;


	double yield_vec[_nBins];
	double yield_vec_err_low[_nBins];
	double yield_vec_err_high[_nBins];
	double scale_vec[_nBins];
	double scale_vec_err_low[_nBins];
	double scale_vec_err_high[_nBins];
	double resol_vec[_nBins];
	double resol_vec_err_low[_nBins];
	double resol_vec_err_high[_nBins];
	
	std::vector<std::vector<double>> stat_error;
	double var_mean[_nBins];
	double hori_low[_nBins];
	double hori_high[_nBins];
	for(int i=0;i<_nBins;i++)
	{
		_count++;
		TCanvas* c= new TCanvas(Form("c%d",_count),"",600,550);
		TCanvas* cMC= new TCanvas(Form("cMC%d",_count),"",600,600);
		//		c->cd();
		//	if(fitOnSaved == 0){
		// tree with selecitons applied, create RooDataSet
		/*	if(isMC==1) {
			skimtree = (TTree*)nt->CopyTree(Form("%s*(%s&&%s>%f&&%s<%f)*(1/%s)&&BgenNew==23333",weightmc.Data(),seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data()), "skimtree");
			}
			else{
			skimtree = (TTree*)nt->CopyTree(   Form("(%s&&%s>%f&&%s<%f)*(1/%s)",                seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data()), "skimtree");		
			std::cout<<"copy tree for data"<<std::endl;
			}
			fout->cd();*/ 
		//skimtree->Write(); 
		//fout->Close();
		//      skimtreeMC = (TTree*)ntMC->CopyTree(Form("%s*(%s&&%s>%f&&%s<%f)",weightmc.Data(),Form("%s&&BgenNew==23333",selmc.Data()),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1]), "skimtreeMC");
		//			ds = new RooDataSet(Form("ds%d",_count),"",skimtree_new, RooArgSet(*mass, *pt));
		//      return;
		//SIMAO 
		RooDataSet* ds_cut;
if(doubly==0) {if(varExp == "Bpt"){
							      ds_cut =
							        new RooDataSet(Form("ds_cut%d",_count),"",ds,
							                       RooArgSet(*mass, *pt, *y, *trackSelection),
							                       Form("(Bpt>%f && Bpt < %f)&&((Bpt < 10 &&  abs(By) > 1.5 ) || (Bpt > 10))",_ptBins[i] , _ptBins[i+1]));}
						        else if(varExp == "By"){
											       ds_cut = new RooDataSet(Form("ds_cut%d", _count),"", ds, RooArgSet(*mass, *pt, *y, *trackSelection), Form("By>%f && By< %f", _ptBins[i], _ptBins[i+1]));}
										else if(varExp == "nMult"){
											       ds_cut = new RooDataSet(Form("ds_cut%d", _count),"", ds, RooArgSet(*mass, *pt, *y, *trackSelection), Form("nMult>%f && nMult< %f", _ptBins[i], _ptBins[i+1]));}
  }
		if(doubly==1)ds_cut = new RooDataSet(Form("ds_cut%d",_count),"",ds, RooArgSet(*mass, *pt, *y, *nMult), Form("%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto)); 
		if(doubly==2)ds_cut = new RooDataSet(Form("ds_cut%d",_count),"",ds, RooArgSet(*mass, *pt, *y), Form("%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto)); 





		//  RooRealVar* w = (RooRealVar*) data->addColumn(wFunc) ;
		//      RooDataSet wdata(data->GetName(),data->GetTitle(),data,*data->get(),0,w->GetName()) ;
		//dsMC = new RooDataSet(Form("dsMC%d",_count),"",skimtreeMC_new,RooArgSet(*mass, *pt));


		RooDataSet* dsMC_cut;
if(doubly==0) {if(varExp == "Bpt"){
      dsMC_cut =
        new RooDataSet(Form("dsMC_cut%d",_count), "", dsMC,
                       RooArgSet(*mass, *pt, *y, *trackSelection), Form("(%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f)&&((Bpt < 10 &&  abs(By) > 1.5 ) || (Bpt > 10))",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto));}
        						else if(varExp=="By" || varExp=="nMult"){
        								       dsMC_cut = new RooDataSet(Form("dsMC_cut%d", _count),"", dsMC,  RooArgSet(*mass, *pt, *y, *trackSelection), Form("%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto));}
        					/*	else if(varExp=="nMult"){
        								       dsMC_cut = new RooDataSet(Form("dsMC_cut%d", _count),"", dsMC,  RooArgSet(*mass, *pt, *y, *trackSelection), Form("%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto));}	*/       
    }
		if(doubly==1) dsMC_cut = new RooDataSet(Form("dsMC_cut%d",_count),"",dsMC, RooArgSet(*mass, *pt, *y, *nMult), Form("%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto), "1"); 
		if(doubly==2) dsMC_cut = new RooDataSet(Form("dsMC_cut%d",_count),"",dsMC, RooArgSet(*mass, *pt,  *y), Form("%s>=%f&&%s<=%f&&Bmass>%f&&Bmass<%f",varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],minhisto, maxhisto), "1"); 


    // Apply track selection cut
    std::cout << "data entries: " << ds_cut->sumEntries() << "\n";
    std::cout << "MC entries: " << dsMC_cut->sumEntries() << "\n";
    ds_cut = (RooDataSet*) ds_cut->reduce(seldata);
    dsMC_cut = (RooDataSet*) dsMC_cut->reduce(selmc);
    std::cout << "data entries: " << ds_cut->sumEntries() << "\n";
    std::cout << "MC entries: " << dsMC_cut->sumEntries() << "\n";

		std::cout<<"Really created datasets"<<std::endl;
				// create RooDataHist
		
		h = new TH1D(Form("h%d",_count),"",nbinsmasshisto,minhisto,maxhisto);
		hMC = new TH1D(Form("hMC%d",_count),"",nbinsmasshisto,minhisto,maxhisto);
		if(isMC==1) skimtree_new->Project(Form("h%d",_count),"Bmass",Form("%s*(%s&&%s>%f&&%s<%f && BgenNew == 23333)*(1/%s)",weightmc.Data(),seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data()));
		else        skimtree_new->Project(Form("h%d",_count),"Bmass",   Form("(%s&&%s>%f&&%s<%f)*(1/%s)",                seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data()));
//		std::string project_str=Form("(%s&&%s>%f&&%s<%f)*(1/%s)",                seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data());
//		std::cout<<"string= "<<project_str<<std::endl;
//		std::string project_str_mc= Form("%s*(%s&&%s>%f&&%s<%f && BgenNew == 23333)*(1/%s)",weightmc.Data(),seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data());
//		std::cout<<"string isMC= "<<project_str_mc<<std::endl;
		skimtreeMC_new->Project(Form("hMC%d",_count),"Bmass",Form("%s*(%s&&%s>%f&&%s<%f)",weightmc.Data(),Form("%s&&BgenNew==23333",selmc.Data()),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1]));
		dh = new RooDataHist(Form("dh%d",_count),"",*mass,Import(*h));
		dhMC = new RooDataHist(Form("dhMC%d",_count),"",*mass,Import(*hMC));
		std::cout<<"RooDataHist "<<std::endl;
		h->SetAxisRange(0,h->GetMaximum()*1.4,"Y");
		outputw->import(*ds);
		outputw->import(*dsMC);
		outputw->import(*dh);
		outputw->import(*dhMC);
		//}
		/*if(fitOnSaved == 1){
		  inputw = (RooWorkspace*)inf->Get("w");
		  ds = (RooDataSet*)inputw->data(Form("ds%d",_count));
		  dsMC = (RooDataSet*)inputw->data(Form("dsMC%d",_count));
		  dh = (RooDataHist*)inputw->data(Form("dh%d",_count));
		  dhMC = (RooDataHist*)inputw->data(Form("dhMC%d",_count));
		  }*/


		
		RooFitResult* f = fit("", "", tree, c, cMC, ds_cut, dsMC_cut, dh, dhMC, mass, frame, _ptBins[i], _ptBins[i+1], isMC, npfit);
		RooFitResult* f_MC = fitMC("", "", tree, c, cMC, ds_cut, dsMC_cut, dh, dhMC, mass, frame, _ptBins[i], _ptBins[i+1], isMC, npfit);


	//	std::cout << "Now Finall We Validate Our Fits" << std::endl;
	//	validate_fit(w_val, tree, varExp, full,centmin, centmax, _ptBins[i], _ptBins[i+1]);
		//std::cout << "Now Finall We Scan Our Significance" << std::endl;

		//scan_significance(w_val, tree, varExp, full,centmin, centmax, _ptBins[i], _ptBins[i+1]);

		/*
		for(int q= 0; q < 100; q++){
			validate_fit(w_val, tree, varExp, full,q);
		}*/
	
	
		//datahist = frame->getHist("ds");
		//TGraphAsymmErrors* datagraph = static_cast<TGraphAsymmErrors*>(datahist);


		//DF TRY
		modelcurve = frame->getCurve(Form("model%d",_count));
		std::cout << "The fit results are" << f->floatParsFinal() << std::endl;
		RooRealVar* fitYield = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index(Form("nsig%d",_count))));

		double yield = fitYield->getVal();
		yieldRec = fitYield->getVal();


		RooRealVar* BackGround = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index(Form("nbkg%d",_count))));
		MyBackground = BackGround->getVal();

		RooRealVar* width_scale = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index("scale")));
		double Myscale = width_scale->getVal();
		double Myscale_err = width_scale->getError();
		scale_vec[i] = Myscale;
		scale_vec_err_low[i] = Myscale_err;
		scale_vec_err_high[i] = Myscale_err;

		nominal_yields.push_back(yield);
		double yieldErr = fitYield->getError();
		printf("yield: %f, yieldErr: %f\n", yield, yieldErr);
		yieldErr = yieldErr*_ErrCor;

		yield_vec[i]=yield;
		yield_vec_err_low[i]=yieldErr;
		yield_vec_err_high[i]=yieldErr;

		yield_vec[i]=yield_vec[i]/(_ptBins[i+1]-_ptBins[i]);
		yield_vec_err_low[i]=yield_vec_err_low[i]/(_ptBins[i+1]-_ptBins[i]);
		yield_vec_err_high[i]=yield_vec_err_high[i]/(_ptBins[i+1]-_ptBins[i]);
		var_mean[i] = (_ptBins[i+1]+_ptBins[i]) / 2;
		hori_low[i] = var_mean[i]-_ptBins[i];
		hori_high[i] = _ptBins[i+1]-var_mean[i];

//Resolution MC
		//TString nominal = "";
		std::cout << "sec1sec1sec1" << endl;
		RooRealVar* sigma1MC = static_cast<RooRealVar*>(f_MC->floatParsFinal().at(f_MC->floatParsFinal().index(Form("sigma1MC%d", _count))));
		//RooRealVar* sigma1 = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index("sigma1")));
		
		//RooRealVar* test_var = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index("nbkg1")));
		//std::cout <<"sigma1 is " << sigmaMC1->getVal() << std::endl;
		std::cout << "sec2sec2sec2" << endl;	
		double Mysigma1MC = sigma1MC->getVal();
		std::cout << "sec3sec3sec3" << endl;	
		double Mysigma1MC_err = sigma1MC->getError();
		std::cout << "sec4sec4sec4" << endl;	

		RooRealVar* sigma2MC = static_cast<RooRealVar*>(f_MC->floatParsFinal().at(f_MC->floatParsFinal().index(Form("sigma2MC%d",_count))));
		double Mysigma2MC = sigma2MC->getVal();
		double Mysigma2MC_err = sigma2MC->getError();

		RooRealVar* weightMC = static_cast<RooRealVar*>(f_MC->floatParsFinal().at(f_MC->floatParsFinal().index(Form("sig1fracMC%d",_count))));
		double MyweightMC  = weightMC->getVal();
		double MyweightMC_err = weightMC->getError();
		
		double resol = sqrt(MyweightMC * pow(Mysigma1MC, 2) + (1 - MyweightMC) * pow(Mysigma2MC, 2));
		double resol_err = (MyweightMC * Mysigma1MC * Mysigma2MC_err + (1-MyweightMC) * Mysigma2MC * Mysigma2MC_err) / resol;

		resol_vec[i] = resol;
		resol_vec_err_low[i] = resol_err;
		resol_vec_err_high[i] = resol_err;
//Resolution MC
 
		std::vector<double> aa;
		double a=yield_vec_err_low[i]/yield_vec[i]*100;
		aa.push_back(a);
		












		
		if(fitOnSaved == 0){
			TH1D* htest = new TH1D(Form("htest%d",_count),"",nbinsmasshisto,minhisto,maxhisto);
			TString sideband = "(abs(Bmass-5.367)>0.2&&abs(Bmass-5.367)<0.3";
			//    	skimtree_new->Project(Form("htest%d",_count),"Bmass",Form("%s&&%s&&%s>%f&&%s<%f)*(1/%s)",sideband.Data(),seldata.Data(),varExp.Data(),_ptBins[i],varExp.Data(),_ptBins[i+1],weightdata.Data()));
			std::cout<<"yield bkg sideband: "<<htest->GetEntries()<<std::endl;
		}

		if(varExp!="nMult"){
			hPt->SetBinContent(i+1,yield/(_ptBins[i+1]-_ptBins[i]));
			hPt->SetBinError(i+1,yieldErr/(_ptBins[i+1]-_ptBins[i]));
		}

		else{
			hPt->SetBinContent(i+1,yield/(_ptBins[1] - _ptBins[0]));
			hPt->SetBinError(i+1,yieldErr/(_ptBins[1] - _ptBins[0]));
		}

		if(f->floatParsFinal().index(Form("nsig%d",_count)) != -1){
			RooRealVar* fitMean = static_cast<RooRealVar*>(f->floatParsFinal().at(f->floatParsFinal().index(Form("mean%d",_count))));
			hMean->SetBinContent(i+1,fitMean->getVal());
			hMean->SetBinError(i+1,fitMean->getError());  
		}



		TLatex* tex_pt;
		TLatex* tex_nMult;
	  TLatex* tex_y;
		TLatex* tex_y1;
		TLatex* tex_y11;
		TLatex* tex_y2;
		

	if(varExp=="Bpt"){
      //for the paper run these
      if (drawLegend) {
        tex_pt = new TLatex(0.55,0.4,Form("%d < p_{T} < %d GeV/c",(int)_ptBins[i],(int)_ptBins[i+1]));
        tex_y = new TLatex(0.55,0.34,"p_{T} > 10 GeV/c : |y| < 2.4");
        tex_y2 = new TLatex(0.55,0.28,"p_{T} < 10 GeV/c : 1.5 < |y| < 2.4");
        tex_y1 = new TLatex(0.55,0.34,"1.5 < |y| < 2.4");
        tex_y11 =new TLatex(0.55,0.34,"|y| < 2.4");
      } else {
        //fr the AN run these
        tex_pt = new TLatex(0.55,0.8,Form("%d < p_{T} < %d GeV/c",(int)_ptBins[i],(int)_ptBins[i+1]));
        tex_y = new TLatex(0.55,0.74,"p_{T} > 10 GeV/c : |y| < 2.4");
        tex_y2 = new TLatex(0.55,0.68,"p_{T} < 10 GeV/c : 1.5 < |y| < 2.4");
        tex_y1 = new TLatex(0.55,0.74,"1.5 < |y| < 2.4");
        tex_y11 =new TLatex(0.55,0.74,"|y| < 2.4");
      }
		}

if(varExp=="By"){
      //for the paper run these
      if (drawLegend) {    
        tex_y = new TLatex(0.55,0.4,Form("%.1f < y < %.1f ",_ptBins[i],_ptBins[i+1]));
      } else {
        //fr the AN run these
        tex_y = new TLatex(0.55,0.8,Form("%.1f < y < %.1f ",_ptBins[i],_ptBins[i+1]));
      }
		}
if(varExp=="nMult"){ 
	//for the paper run these
      if (drawLegend) {
      	tex_nMult = new TLatex(0.55,0.4,Form("%d < nTrks < %d",(int)_ptBins[i],(int)_ptBins[i+1]));
      } else {
        //fr the AN run these
        tex_nMult = new TLatex(0.55,0.8,Form("%d < nTrks < %d",(int)_ptBins[i],(int)_ptBins[i+1]));
      }
	}
/*	if(varExp=="By"){
	tex_pt = new TLatex(0.21,0.75,"5 < p_{T} < 50 GeV/c");
	tex_y = new TLatex(0.21,0.69,Form("%.1f < |y| < %.1f",_ptBins[i],_ptBins[i+1]));
		}*/
		
if (varExp=="Bpt"){
	

		tex_pt->SetNDC();
		tex_pt->SetTextFont(42);
		tex_pt->SetTextSize(0.045);
		tex_pt->SetLineWidth(2);
		tex_pt->Draw();

		tex_y->SetNDC();
		tex_y->SetTextFont(42);
		tex_y->SetTextSize(0.045);
		tex_y->SetLineWidth(2);

		tex_y1->SetNDC();
		tex_y1->SetTextFont(42);
		tex_y1->SetTextSize(0.045);
		tex_y1->SetLineWidth(2);
		tex_y11->SetNDC();
		tex_y11->SetTextFont(42);
		tex_y11->SetTextSize(0.045);
		tex_y11->SetLineWidth(2);
	
		tex_y2->SetNDC();
		tex_y2->SetTextFont(42);
		tex_y2->SetTextSize(0.045);
		tex_y2->SetLineWidth(2);


		if(_ptBins[i] >= 10){tex_y11->Draw();}
		else if(_ptBins[i+1]==50){tex_y->Draw();
					    tex_y2->Draw();}
		else{tex_y1->Draw();}
	} 
	else if (varExp=="By") {

		tex_y->SetNDC();
		tex_y->SetTextFont(42);
		tex_y->SetTextSize(0.045);
		tex_y->SetLineWidth(2);
		tex_y->Draw();}
	else if (varExp=="nMult") {
		tex_nMult->SetNDC();
		tex_nMult->SetTextFont(42);
		tex_nMult->SetTextSize(0.045);
		tex_nMult->SetLineWidth(2);
		tex_nMult->Draw();
			}
	CMS_lumi(c,19011,0);
	c->Update();
	TLatex* texB;
	if(tree=="ntphi"){ texB = new TLatex(0.21,0.82, "B^{0}_{s}");}
	if(tree=="ntKp"){ texB = new TLatex(0.21,0.8, "B^{#pm}");}
	texB->SetNDC();
	texB->SetTextFont(62);
	texB->SetTextSize(0.065);
	texB->SetLineWidth(2);
	texB->Draw();

/*
		TLatex *lat = new TLatex();
		lat->SetNDC();
		lat->SetTextSize(0.035);
	
		lat->DrawLatex(0.64,0.85,Form("S = %.1f", yield));
		lat->DrawLatex(0.64,0.80,Form("S_err = %.1f", yieldErr));		
		lat->DrawLatex(0.64,0.75,Form("B = %.1f", bkgd));
		//lat->DrawLatex(0.48,0.70,Form("Significance: S/#sqrt{S+B} = %.1f", Significance));
		lat->DrawLatex(0.64,0.70,Form("Significance: %.1f", real_significance));
*/
		gSystem->mkdir("filesbp",true); 
		gSystem->mkdir("filesbs",true); 
		gSystem->mkdir(Form("%s", outplotf.Data()),true); 
		//c->SaveAs(Form("%s%s/%s_%s_%s_%d_%d_cutY%d_",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),varExp.Data(),(int)_ptBins[i],(int)_ptBins[i+1], doubly)+tree+".pdf");
		c->SaveAs(Form("%s%s/%s_%s_%s_%d_%d_cutY%d_",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),varExp.Data(),(int)_ptBins[i],(int)_ptBins[i+1], doubly)+tree+".png");
		//c->SaveAs(Form("%s%s/%s_%s_%d%s_%s_%d%d_cutY%d",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),_postfix.Data(),varExp.Data(),(int)_ptBins[i],(int)_ptBins[i+1], doubly)+tree+".C");
		cMC->SaveAs(Form("%s%s/%s_%s_%s_%d_%d_cutY%d_",outplotf.Data(),_prefix.Data(),"mc",_isPbPb.Data(),varExp.Data(), (int)_ptBins[i], (int)_ptBins[i+1], doubly)+tree+".pdf");
		// cMC->SaveAs(Form("%s%s/%s_%s_%s_%d_%d_cutY%d_",outplotf.Data(),_prefix.Data(),"mc",_isPbPb.Data(),varExp.Data(), (int)_ptBins[i], (int)_ptBins[i+1], doubly)+tree+".png");



	/*	    TH1* h = dh->createHistogram("Bmass");
				h->GetEntries();
				h->Sumw2(kFALSE);
				h->SetBinErrorOption(TH1::kPoisson);roofitB.C
				TCanvas* cpull= new TCanvas(Form("cpull%d",_count),"",600,600);
				cpull->cd();
				TGraphAsymmErrors* pullgraph = new TGraphAsymmErrors();
				pullgraph->SetTitle("");
				pullgraph->SetMaximum(5);
				pullgraph->SetMinimum(-5);
				pullgraph->SetMarkerSize(1.55); pullgraph->SetMarkerStyle(20); pullgraph->SetLineColor(1); pullgraph->SetLineWidth(4);
				double x; double xfit; double y; double yfit;
				double binerr;
				for(int b = 0; b < h->GetNbinsX(); b++){
				modelcurve->GetPoint(modelcurve->findPoint(h->GetBinCenter(b+1)),xfit,yfit);
				binerr = h->GetBinContent(b+1) > yfit ? h->GetBinErrorLow(b+1) : h->GetBinErrorUp(b+1);
				pullgraph->SetPoint(b,h->GetBinCenter(b+1),(h->GetBinContent(b+1)-yfit)/binerr);
		//pullgraph->SetPointEYlow(b,h->GetBinErrorLow(b+1)/binerr);
		//pullgraph->SetPointEYhigh(b,h->GetBinErrorUp(b+1)/binerr);
		pullgraph->SetPointEYlow(b,1);
		pullgraph->SetPointEYhigh(b,1);
		}
		TLine* line = new TLine(5., 0., 6., 0.);
		line->SetLineStyle(9);
		line->SetLineWidth(4);
		line->SetLineColor(kGreen+1);
		line->Draw();
		pullgraph->Draw();
		cpull->SaveAs(Form("%s%s/%s_%s_%d%s_pull.pdf",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),_count,_postfix.Data()));*/

		RooCurve* modelcurve_back = new RooCurve();
		std::vector<double> back_variation; 
		std::vector<double> back_err;  
		RooCurve* modelcurve_signal = new RooCurve();
		std::vector<double> signal_variation; 
		std::vector<double> signal_err;
		std::vector<double> general_err;  
		double max_signal=0.;
		double max_back=0.;
		double full_err=0;
	
		if(syst==1){

				for(int j=0; j<background.size(); j++){
				RooFitResult* f_back = fit("background", background[j], tree, c, cMC, ds_cut, dsMC_cut, dh, dhMC, mass, frame, _ptBins[i], _ptBins[i+1], isMC, npfit);
				
				texB->Draw();
				
				if (varExp=="Bpt"){
					tex_pt->Draw();
					if(_ptBins[i] >= 10){tex_y11->Draw();}
					else if(_ptBins[i+1]==50){
							tex_y->Draw();
							tex_y2->Draw();}
					else{tex_y1->Draw();}
				} else if (varExp=="By"){tex_y->Draw();}
					else if (varExp=="nMult"){tex_nMult->Draw();}
				CMS_lumi(c,19011,0);
				c->Update();
				
				c->SaveAs(Form("%s/%s_%s_%s_%d_%d_%s_cutY%d_", outplotf.Data(), _isMC.Data(), _isPbPb.Data(), varExp.Data(),(int)_ptBins[i],(int)_ptBins[i+1],background[j].c_str(), doubly)+tree+".png");
				//c->SaveAs(Form("%s/%s_%s_%s_%d_%d_%s_cutY%d_", outplotf.Data(), _isMC.Data(), _isPbPb.Data(), varExp.Data(),(int)_ptBins[i],(int)_ptBins[i+1],background[j].c_str(), doubly)+tree+".pdf");
			
				modelcurve_back = frame->getCurve(Form("model%d",_count));
				RooRealVar* fitYield_back = static_cast<RooRealVar*>(f_back->floatParsFinal().at(f_back->floatParsFinal().index(Form("nsig%d",_count))));
				back_variation.push_back(fitYield_back->getVal());
				back_err.push_back(abs(((yield-fitYield_back->getVal())/yield)*100));
				if(abs(((yield-fitYield_back->getVal())/yield)*100)>max_back) max_back=abs(((yield-fitYield_back->getVal())/yield)*100);
			}

			general_err.push_back(max_back);

			for(int j=0; j<signal.size(); j++){
				RooFitResult* f_signal = fit("signal", signal[j], tree, c, cMC, ds_cut, dsMC_cut, dh, dhMC, mass, frame, _ptBins[i], _ptBins[i+1], isMC, npfit);
				texB->Draw();

				if (varExp=="Bpt"){
					tex_pt->Draw();
					if(_ptBins[i] >= 10){tex_y11->Draw();}
					else if(_ptBins[i+1]==50){
							tex_y->Draw();
							tex_y2->Draw();}
					else{tex_y1->Draw();}
				} else if (varExp=="By"){tex_y->Draw();}
					else if (varExp=="nMult"){tex_nMult->Draw();}

				CMS_lumi(c,19011,0);
				c->Update();

				//tex_y->Draw();
			//	c->SaveAs(Form("%s%s/%s_%s_%d%s_%s_%d_%d_%s_cutY%d",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),_count,_postfix.Data(),varExp.Data(),(int)_ptBins[i],(int)_ptBins[i+1],signal[j].c_str(), doubly)+tree+".pdf");
				cMC->SaveAs(Form("%s%s/%s_%s_%s_%d_%d_%s_cutY%d_",outplotf.Data(),_prefix.Data(),"mc",_isPbPb.Data(),varExp.Data(), (int)_ptBins[i], (int)_ptBins[i+1],signal[j].c_str(), doubly)+tree+".pdf");
				c->SaveAs(Form("%s%s/%s_%s_%s_%d_%d_%s_cutY%d_",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),varExp.Data(),(int)_ptBins[i],(int)_ptBins[i+1],signal[j].c_str(), doubly)+tree+".png");
				//c->SaveAs(Form("%s%s/%s_%s_%s_%d_%d_%s_cutY%d_",outplotf.Data(),_prefix.Data(),_isMC.Data(),_isPbPb.Data(),varExp.Data(),(int)_ptBins[i],(int)_ptBins[i+1],signal[j].c_str(), doubly)+tree+".pdf");
				modelcurve_signal = frame->getCurve(Form("model%d",_count));
				RooRealVar* fitYield_signal = static_cast<RooRealVar*>(f_signal->floatParsFinal().at(f_signal->floatParsFinal().index(Form("nsig%d",_count))));
				signal_variation.push_back(fitYield_signal->getVal());
				signal_err.push_back(abs(((yield-fitYield_signal->getVal())/yield)*100));
				if(abs(((yield-fitYield_signal->getVal())/yield)*100)>max_signal) max_signal=abs(((yield-fitYield_signal->getVal())/yield)*100);
			}

			general_err.push_back(max_signal);
			full_err=sqrt(max_back*max_back+max_signal*max_signal);
			general_err.push_back(full_err);

		}


		if(syst==1){
			stat_error.push_back(aa);
			background_syst.push_back(back_variation);
			signal_syst.push_back(signal_variation);
			back_syst_rel_values.push_back(back_err);
			sig_syst_rel_values.push_back(signal_err);
			general_syst.push_back(general_err);

		}
	
std::cout << "The size of general_err is " << general_err.size() << std::endl;		

	//validate_fit(outputw,1,Form("model%d",_count),Form("nsig%d",_count));
	}


	hMean->Write();
	hPt->Write();
	outf->Close();	

	string Path;
	if(tree == "ntphi"){ Path = "./filesbs/";}
	else if (tree == "ntKp"){Path = "./filesbp/";}
	std::ofstream myfile;
	myfile.open (Path + "systematics_" + tree.Data() + ".txt");

	if(syst==1){ 
		for(int i=0; i<_nBins; i++){
			std::cout<<"pt bin = "<<_ptBins[i]<<" "<<_ptBins[i+1]<<std::endl;
			std::cout<<"nominal yield = "<<nominal_yields[i]<<std::endl;
			std::cout<<"count= "<<i<<std::endl;

			for(int j=0; j<background.size(); j++){
				std::cout<<" back sys in bin "<<i<<" ; with pdf "<< background[j] << " ="<<background_syst[i][j]<<std::endl;
				std::cout<<" back sys in bin "<<i<<" ; with pdf "<< background[j] << " ="<<back_syst_rel_values[i][j]<<" % "<<std::endl;
				myfile<<" back sys in bin "<<i<<" ; with pdf "<< background[j] << " ="<<back_syst_rel_values[i][j]<<" % "<<std::endl;

			}
			for(int j=0; j<signal.size(); j++){
				std::cout<<" signal sys in bin "<<i<<" ; with pdf "<< signal[j] << " ="<<signal_syst[i][j]<<std::endl;
				std::cout<<" signal sys in bin "<<i<<" ; with pdf "<< signal[j] << " ="<<sig_syst_rel_values[i][j]<<" % "<<std::endl;
				myfile<<" signal sys in bin "<<i<<" ; with pdf "<< signal[j] << " ="<<sig_syst_rel_values[i][j]<<" % "<<std::endl;

			}
		}

	}


	myfile.close();

	/*if(fitOnSaved == 1){
	  outf->Close();	
	  return;
	  }
	  ntMC->Project("hPtMC",varExp.Data(),TCut(weightmc)*(TCut(selmc.Data())&&"(BgenNew==23333)"));
	  divideBinWidth(hPtMC);
	//	ntGen->Project("hPtGen","Gpt",TCut(weightgen)*(TCut(selmcgen.Data())));
	//	divideBinWidth(hPtGen);*/



	TCanvas* cPt =  new TCanvas("cPt","",600,600);
	cPt->SetLogy();
	hPt->SetXTitle("B_{s} p_{T} (GeV/c)");
	hPt->SetYTitle("Uncorrected dN(B_{s})/dp_{T}");
	hPt->Draw();


	std::vector<std::string> labels_back = {"Linear", "2nd Poly", "mass range" };
	std::vector<std::string> col_name_back;
	std::vector<std::string> labels_signal = {"Triple Gaussian", "Fixed Mean", "CB+Gaussian"};
	std::vector<std::string> labels_general = {"Background", "Signal", "Total"};
	std::vector<std::string> labels_general_stat = {"Statistical error"};
	std::vector<std::string> col_name_general;
	std::vector<std::string> col_name_signal;
	std::vector<std::string> col_name_general_stat;

	string name;
	int m=_nBins;
	
	col_name_back.push_back("Background Model");
	col_name_signal.push_back("Signal Model");
	col_name_general.push_back("Systematic Source");
	col_name_general_stat.push_back(" ");

	if(varExp=="Bpt"){name="$<p_T<$";} else if(varExp=="By"){name="$<y<$";} else if(varExp=="nMult"){name="$<nTrks<$";}
	for(int i=0;i<m;i++){
		std::ostringstream collabel;
		std::ostringstream siglabel;
		std::ostringstream genlabel;
		std::ostringstream statlabel;
		collabel<<_ptBins[i]<<name<<_ptBins[i+1];
		siglabel<<_ptBins[i]<<name<<_ptBins[i+1];
		genlabel<<_ptBins[i]<<name<<_ptBins[i+1];
		statlabel<<_ptBins[i]<<name<<_ptBins[i+1];
		std::string label1 = collabel.str();
		col_name_back.push_back(label1);
		col_name_signal.push_back(label1);
		col_name_general.push_back(label1);
		col_name_general_stat.push_back(label1);

	
	}
	



//for(int t; t< sizeof(_ptBins)/sizeof(_ptBins[0]);t++){cout <<"__"<<  _ptBins[t]<<"__";}
//	cout<<stat_error.size()<<general_syst.size()<<endl;
	//stat_error.push_back(aa);
	if(syst==1 && full==0){
		gSystem->mkdir("tabels",true); 
		latex_table(Path + "background_systematics_table_"+std::string (varExp.Data())+"_"+std::string (tree.Data()), _nBins+1,  (int)(1+background.size()),  col_name_back,labels_back,back_syst_rel_values, "Background PDF Systematic Errors");
		latex_table(Path + "signal_systematics_table_"+std::string (varExp.Data())+"_"+std::string (tree.Data()), _nBins+1, (int)(1+signal.size()),    col_name_signal, labels_signal,sig_syst_rel_values, "Signal PDF Systematic Errors");
		latex_table(Path + "general_systematics_table_"+std::string (varExp.Data())+"_"+std::string (tree.Data()),  _nBins+1, 4 , col_name_general, labels_general, general_syst, "Overall PDF Variation Systematic Errors");	
		latex_table(Path + "Statistical_error_table_"+std::string (varExp.Data())+"_"+std::string (tree.Data()),  _nBins+1, 2 , col_name_general_stat, labels_general_stat, stat_error, "Statistical error");	
		

		double zero[_nBins];
			for (int i=0;i<_nBins;i++){zero[i]=0.;}

		TCanvas* c_back= new TCanvas();
		TLegend* legback=new TLegend(0.7,0.7,0.9,0.9);
		TMultiGraph* m_back= new TMultiGraph();
		const char* backlabel[3]={"Linear", "2nd Poly", "mass range"};
		for (int j=0;j<(int)(background.size());j++){
			Double_t x[_nBins],y[_nBins];
			for (int i=0;i<_nBins;i++){
				x[i]=(_ptBins[i]+_ptBins[i+1])/2;
				y[i]=	back_syst_rel_values[i][j];
		}
			TGraph *g_back= new TGraphAsymmErrors (_nBins,x,y,zero,zero,zero,zero);
			g_back->SetMarkerColor(j+1);
			g_back->SetMarkerStyle(21);
			m_back->Add(g_back);
			legback->AddEntry(g_back,backlabel[j],"p");
	}
		m_back->GetXaxis()->SetTitle(varExp.Data());
		m_back->GetYaxis()->SetTitle("Systematic Uncertainty(%)");
		m_back->Draw("A*");
		legback->Draw();
		c_back->SaveAs(Form("./tabels/background_systematics_plot_%s_%s.png",tree.Data(),varExp.Data())); 

	TCanvas* c_sig= new TCanvas();
	TLegend* legsig=new TLegend(0.7,0.7,0.9,0.9);
	TMultiGraph* m_sig= new TMultiGraph();
	const char* siglabel[4]={"Triple Gaussian", "Fixed Mean", "CB+Gaussian"};
	for (int j=0;j<(int)(signal.size());j++){
		Double_t x[_nBins],y[_nBins];
		for (int i=0;i<_nBins;i++){
			x[i]=(_ptBins[i]+_ptBins[i+1])/2;
			y[i]=	sig_syst_rel_values[i][j];
	}
		TGraph *g_sig= new TGraphAsymmErrors (_nBins,x,y,zero,zero,zero,zero);
		g_sig->SetMarkerColor(j+1);
		g_sig->SetMarkerStyle(21);
		m_sig->Add(g_sig);
		legsig->AddEntry(g_sig,siglabel[j],"p");
}
	m_sig->GetXaxis()->SetTitle(varExp.Data());
	m_sig->GetYaxis()->SetTitle("Systematic Uncertainty(%)");
	m_sig->Draw("A*");
	legsig->Draw();
	c_sig->SaveAs(Form("./tabels/signal_systematics_plot_%s_%s.png",tree.Data(),varExp.Data())); 

	

	TCanvas* c_gen= new TCanvas();
	TLegend* legen=new TLegend(0.7,0.7,0.9,0.9);
	TMultiGraph* m_gen= new TMultiGraph();
	const char* genlabel[3]={"Background", "Signal", "Total"};
	for (int j=0;j<(int)(labels_general.size());j++){
		Double_t x[_nBins],y[_nBins];
		for (int i=0;i<_nBins;i++){
			x[i]=(_ptBins[i]+_ptBins[i+1])/2;
			y[i]=	general_syst[i][j];
	}
		TGraph *g_gen= new TGraphAsymmErrors (_nBins,x,y,zero,zero,zero,zero);
		g_gen->SetMarkerColor(j+1);
		g_gen->SetMarkerStyle(21);
		m_gen->Add(g_gen);
		legen->AddEntry(g_gen,genlabel[j],"p");
}
	Double_t x[_nBins],y[_nBins];
	for (int i=0;i<_nBins;i++){
			x[i]=(_ptBins[i]+_ptBins[i+1])/2;
			y[i]=	stat_error[i][0];
	}
		TGraph *g_gen= new TGraphAsymmErrors (_nBins,x,y,zero,zero,zero,zero);
		g_gen->SetMarkerColor(9);
		g_gen->SetMarkerStyle(21);
		m_gen->Add(g_gen);
		legen->AddEntry(g_gen,"Statistical","p");



	m_gen->GetXaxis()->SetTitle(varExp.Data());
	m_gen->GetYaxis()->SetTitle("Total Uncertainty(%)");
	m_gen->Draw("A*");
	legen->Draw();
	c_gen->SaveAs(Form("./tabels/general_systematics_plot_%s_%s.png",tree.Data(),varExp.Data())); 

	}



	cout << "Final Background = " << MyBackground << endl;

	cout << "Final Yield = " << yieldRec << endl;

// Differential plot part starts
	 gSystem->mkdir("Graphs",true); 
	 TCanvas c_diff;
	 TMultiGraph* mg = new TMultiGraph();

	 TGraphAsymmErrors* gr_staterr = new TGraphAsymmErrors(_nBins,var_mean,yield_vec,hori_low,hori_high,yield_vec_err_low,yield_vec_err_high);
	 gr_staterr->SetLineColor(1); 
	
	 if(varExp == "By"){
		 mg->GetXaxis()->SetTitle("Rapidity (y)");
		 mg->GetYaxis()->SetTitle("dN_{S}/dy");
		 mg->GetXaxis()->SetLimits(-2.4 ,2.4);
	 }
	 if(varExp == "Bpt"){
		 mg->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		 mg->GetYaxis()->SetTitle("dN_{S}/dp_{T}");
		if (tree == "ntKp"){ mg->GetXaxis()->SetLimits(0 ,80); }
		if (tree == "ntphi"){ mg->GetXaxis()->SetLimits(0 ,60); }
	 }
	 if(varExp == "nMult"){
		 mg->GetXaxis()->SetTitle("Multiplicity (Mult)");
		 mg->GetYaxis()->SetTitle("dN_{S}/dMult");
		 mg->GetXaxis()->SetLimits(0, 110);
	 }

	 mg->Add(gr_staterr);
	 mg->Draw("ap");
	 //mg->SetTitle("Differential Signal Yield");  
	 
         TLegend *leg_d = new TLegend(0.7,0.7,0.9,0.9);
	 leg_d->AddEntry(gr_staterr, "Statistical Uncertainty", "e");
	 //leg_d->AddEntry(grs, "Systematic Uncertainty", "e");
	 leg_d->SetBorderSize(0);
	 leg_d->SetFillStyle(0);
	 leg_d->SetTextSize(0);
	 leg_d->Draw();

	 const char* pathc =Form("./Graphs/raw_yield_%s_%s.png",tree.Data(),varExp.Data()); 
	 c_diff.SaveAs(pathc);
// Differential plot part ends

// Parameters vs variables part starts
	 	 TCanvas c_par;
	 TMultiGraph* mg_par = new TMultiGraph();

	 TGraphAsymmErrors* gr_scale = new TGraphAsymmErrors(_nBins,var_mean,scale_vec,hori_low,hori_high,scale_vec_err_low,scale_vec_err_high);
	 gr_scale->SetLineColor(1); 
	
	 if(varExp == "By"){
		 mg_par->GetXaxis()->SetTitle("Rapidity (y)");
		 mg_par->GetYaxis()->SetTitle("Mass resolution scale factor");
		 mg_par->GetXaxis()->SetLimits(-2.4 ,2.4);
		 mg_par->GetYaxis()->SetLimits(0, 2.0);

	 }
	 if(varExp == "Bpt"){
		 mg_par->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		 mg_par->GetYaxis()->SetTitle("Mass resolution scale factor");
		 if (tree == "ntKp"){ mg_par->GetXaxis()->SetLimits(0 ,80); }
		 if (tree == "ntphi"){ mg_par->GetXaxis()->SetLimits(0 ,60); }
		 mg_par->GetYaxis()->SetLimits(0, 2.0);
	 }
	 if(varExp == "nMult"){
		 mg_par->GetXaxis()->SetTitle("Multiplicity (Mult)");
		 mg_par->GetYaxis()->SetTitle("Mass resolution scale factor");
		 mg_par->GetXaxis()->SetLimits(0, 110);
		 mg_par->GetYaxis()->SetLimits(0, 2.0);
	 }

	 mg_par->Add(gr_scale);
	 mg_par->Draw("ap");
	 //mg->SetTitle("Differential Signal Yield");  
/*	 
         TLegend *leg_par = new TLegend(0.7,0.7,0.9,0.9);
	 leg_par->AddEntry(gr_scale, "Scale", "e");
	 //leg_d->AddEntry(grs, "Systematic Uncertainty", "e");
	 leg_par->SetBorderSize(0);
	 leg_par->SetFillStyle(0);
	 leg_par->SetTextSize(0);
	 leg_par->Draw();
*/
	 const char* pathc_par =Form("./Graphs/parameters_variation_%s_%s.png",tree.Data(),varExp.Data()); 
	 c_par.SaveAs(pathc_par);
//Parameters vs variables part ends

//Resolution plot part starts
	 	 TCanvas c_resol;
	 TMultiGraph* mg_resol = new TMultiGraph();

	 TGraphAsymmErrors* gr_resol = new TGraphAsymmErrors(_nBins,var_mean,resol_vec,hori_low,hori_high,resol_vec_err_low,resol_vec_err_high);
	 gr_resol->SetLineColor(1); 
	
	 if(varExp == "By"){
		 mg_resol->GetXaxis()->SetTitle("Rapidity (y)");
		 mg_resol->GetYaxis()->SetTitle("Resolution(MC)");
		 mg_resol->GetXaxis()->SetLimits(-2.4 ,2.4);
		 //mg_par->GetYaxis()->SetLimits(0, 2.0);

	 }
	 if(varExp == "Bpt"){
		 mg_resol->GetXaxis()->SetTitle("Transverse Momentum (p_{T})");
		 mg_resol->GetYaxis()->SetTitle("Resolution(MC)");
		 if (tree == "ntKp"){ mg_resol->GetXaxis()->SetLimits(0 ,80); }
		 if (tree == "ntphi"){ mg_resol->GetXaxis()->SetLimits(0 ,60); }
		 //mg_par->GetYaxis()->SetLimits(0, 2.0);
	 }
	 if(varExp == "nMult"){
		 mg_resol->GetXaxis()->SetTitle("Multiplicity (Mult)");
		 mg_resol->GetYaxis()->SetTitle("Resolution(MC)");
		 mg_resol->GetXaxis()->SetLimits(0, 110);
		 //mg_par->GetYaxis()->SetLimits(0, 2.0);
	 }

	 mg_resol->Add(gr_resol);
	 mg_resol->Draw("ap");

	 const char* pathc_resol =Form("./Graphs/resolution_%s_%s.png",tree.Data(),varExp.Data()); 
	 c_resol.SaveAs(pathc_resol);

//Resolution plot part ends

}
