#include <string.h>
#include "TChain.h"
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <TString.h>
#include <TH1I.h>
#include <TFile.h>
#include <vector>
#include <TH1F.h>
#include <TList.h>
#include <TH1D.h>
#include <math.h>
#include <TAxis.h>
#include <map>
 
 
TH1D * hist_shift(TH1D * h, int shift = 1, int xmin = 0, int xmax = 100, int void_vals = 0);
 
int get_efficiency(){ 
	
	TFile *f = new TFile("/home/csc2168/particle/stilbene/histos3/full_corr_ind_result.root");
	//TFile *f2 = new TFile("/home/csc2168/particle/stilbene/histos3/result_ind_early.root");
	TFile *f3 = new TFile("/home/csc2168/particle/stilbene/histos3/fission_result.root");
	
	TH1D *h0 = (TH1D*)f->Get("hit_num_ch2");
    TH1D *h1 = (TH1D*)f3->Get("tof_fission");
	TH1D *h2 = (TH1D*)f->Get("tof_stilbene");
	TH1D *h3 = (TH1D*)f->Get("live_time");
	TH1D h4;
	TH1D *h5 = (TH1D*)f->Get("tof_stilbene");
	
	
//##### Perform Fits

//Fit Fission

	TF1 *func1 = new TF1("func1", "gaus", -424, -418); //gamma
//	TF1 *func2 = new TF1("func2", "expo", -408, -350); //neutron
	//TF1 *func3 = new TF1("func3", "gaus(0) + gaus(3)", -550, 450);
	h1->Fit("func1", "R");
	h1->Fit("func2", "R+");
	
	TF1 *func4 = new TF1("func4", "gaus", -525, -510); //gamma
	//TF1 *func5 = new TF1("func5", "gaus", -510, -505); //neutron
	//TF1 *func6 = new TF1("func6", "gaus(0) + gaus(3)", -550, 450);
	h5->Fit("func4", "R");
	//h5->Fit("func5", "R+");
	
	
	printf("%f\n", (double)h1->GetMaximumBin());
	
//##### Live Time

	 double trigs = h0->GetEntries();
	 printf("%f\n", trigs);
     h3->Scale(1/trigs);
	 h2->Divide(h3);
	


	
	//###################################################
	//Save all histos
	
	
	TList *l = new TList();
	//l->Add(h1);
	//l->Add(h2);
	l->Add(h5);
	
	TIter iter(l);
	TFile* f4 = new TFile(Form("stilbene_fit.root"),"recreate");
	TH1D* h;

	while((h = (TH1D*)iter.Next())){
		
		h->Write();
	}
	
	f->Close();
	//f2->Close();
	f3->Close();  
	f4->Close();
	return 0;
}

TH1D *hist_shift(TH1D * h, int shift =1, int xmin = 0, int xmax = 100, int void_vals = 0){ //void vals is what to set the eliminated bins
	int nbins = h->GetNbinsX();

	TH1D *h2 = new TH1D(h->GetName(), h->GetTitle(), nbins, xmin, xmax);
	if(shift > 0){
		for(int i = 1; i < shift; i++){
			h2->SetBinContent(i, void_vals);
		}
		for(int i = shift; i< nbins; i++){
			h2->SetBinContent(i, h->GetBinContent(i - shift));
		}
	}
	else if(shift < 0){
		shift = -1*shift;
		for(int i = (nbins - shift); i < nbins; i++){
			h2->SetBinContent(i, void_vals);
		}
		for(int i = 1; i< (nbins - shift); i++){
			h2->SetBinContent(i, h->GetBinContent(i + shift));
		}
	}
	else{
		return h2;
	}
		return h2;
}








