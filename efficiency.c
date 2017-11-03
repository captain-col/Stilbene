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
 

int get_efficiency(){ 
	
	TFile *f = new TFile("/home/csc2168/particle/stilbene/histos/result.root");
	TFile *f2 = new TFile("/home/csc2168/particle/stilbene/histos/result.root");
	TFile *f3 = new TFile("/home/csc2168/particle/stilbene/histos/fission_result.root");
	
    TH1D *h1 = (TH1D*)f3->Get("tof_fission");
	TH1D *h2 = (TH1D*)f->Get("tof_stilbene");
	TH1D *h3 = (TH1D*)f->Get("live_time");
	TH1D h4;
	TH1D *h5 = (TH1D*)f2->Get("tof_stilbene");
	
	
//##### Perform Fits

//Fit Fission

	TF1 *func1 = new TF1("func1", "gaus", -422, -418); //gamma
	TF1 *func2 = new TF1("func2", "gaus", -410, -360); //neutron
	//TF1 *func3 = new TF1("func3", "gaus(0) + gaus(3)", -550, 450);
	h1->Fit("func1", "R");
	h1->Fit("func2", "R+");
	
	TF1 *func4 = new TF1("func4", "gaus", -530, -508); //gamma
	TF1 *func5 = new TF1("func5", "gaus", -505, -492); //neutron
	//TF1 *func6 = new TF1("func6", "gaus(0) + gaus(3)", -550, 450);
	h5->Fit("func4", "R");
	h5->Fit("func5", "R+");
	
	
//##### Live Time

	 h2->Divide(h3);
	


	
	//###################################################
	//Save all histos
	
	
	TList *l = new TList();
	l->Add(h1);
	l->Add(h2);
	l->Add(h5);
	
	TIter iter(l);
	TFile* f4 = new TFile(Form("eff_res.root"),"recreate");
	TH1D* h;

	while((h = (TH1D*)iter.Next())){
		
		h->Write();
	}
	
	f->Close();
	f2->Close();
	f3->Close();
	f4->Close();
	return 0;
}




