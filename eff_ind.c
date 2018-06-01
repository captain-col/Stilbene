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
  

int get_efficiency(int file_begin, int file_end, int shift = 0){ 
	
	char filename[100];
	sprintf(filename, "/home/csc2168/particle/stilbene/histos3/histos_%d_%d.root", file_begin, file_end);
	
	TFile *f;
	TFile *f2;
	TH1D *h1;
	TH1D *h2;
	TH1D *h3;
	
	
	if(shift == 0){
		f = new TFile(filename);	

		h1 = (TH1D*)f->Get("hit_num_ch2");
		h2 = (TH1D*)f->Get("tof_stilbene");
		h3 = (TH1D*)f->Get("live_time");
	}
	else{
		char filename2[100];
		sprintf(filename2, "/home/csc2168/particle/stilbene/histos3/histos_%d_%d_shift.root", file_begin, file_end);
		
		f = new TFile(filename);	
		h1 = (TH1D*)f->Get("hit_num_ch2");
		
		f2 = new TFile(filename2);	
		h2 = (TH1D*)f2->Get("tof_stilbene");
		h3 = (TH1D*)f2->Get("live_time");
		
	}
	
	
//##### Perform Fits
	
	TF1 *func4 = new TF1("func4", "gaus", -525, -510); //gamma
	TF1 *func5 = new TF1("func5", "gaus", -505, -495); //neutron
	//TF1 *func6 = new TF1("func6", "gaus(0) + gaus(3)", -550, 450);
	//h5->Fit("func4", "R");
//	h5->Fit("func5", "R+");
	
	
	printf("Max of Stilbene %d\n", h2->GetMaximumBin() -800);
	printf("Min of Live Time %d\n", h3->GetMinimumBin() -800);
	
//##### Live Time

	 double trigs = h1->GetEntries();
	 printf("%f\n", trigs);
	 h3->Scale(1/trigs);
	 
	 h2->Divide(h3);
	 h3->SetTitle("Dead Time Efficiency");
	 h2->SetTitle("Dead Time Corrected Stilbene");
	 


	
	//###################################################
	//Save all histos
	
	
	TList *l = new TList();
	l->Add(h2);
	l->Add(h3);
	
	char outfile[100];
	sprintf(outfile, "corr_eff_res_%d_%d.root", file_begin, file_end);
	TIter iter(l);
	TFile* f4 = new TFile(outfile,"recreate");
	TH1D* h;

	while((h = (TH1D*)iter.Next())){
		h->Write();
	}
	
	f->Close();
	f4->Close();
	if(shift != 0){f2->Close();}
	return 0;
}