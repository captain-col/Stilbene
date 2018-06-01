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
 
int correct_stil(){ 
	
	TFile *f = new TFile("/home/csc2168/particle/stilbene/histos3/neutrons_fission.root");
	TFile *f2 = new TFile("/home/csc2168/particle/stilbene/histos3/neutrons_stilbene.root");
	
	//TH1D *h0 = (TH1D*)f->Get("hit_num_ch2");
    TH1D *h1 = (TH1D*)f->Get("tof_fission");
	TH1D *h2 = (TH1D*)f2->Get("tof_stilbene");

    printf("Fission Integral: %f\n" , h1->Integral(280, 1000));
	printf("Stilbene Integral: %f\n", h2->Integral(280, 1000));	
	

	
//##### Fission Efficiency

	 double fission_eff = 2.52*pow(10, -8);
     h1->Scale(1/fission_eff);
	 
	 h2->Divide(h1);
	
	

	
	//###################################################
	//Save all histos
	
	
	TList *l = new TList();
	l->Add(h1);
	l->Add(h2);

	
	TIter iter(l);
	TFile* f4 = new TFile(Form("final_stilbene.root"),"recreate");
	TH1D* h;

	while((h = (TH1D*)iter.Next())){
		
		h->Write();
	}
	
	f->Close();
	f2->Close();
	//f3->Close(); 
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