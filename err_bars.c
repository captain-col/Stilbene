double *prop_err(double *err, double *err2, int nbins)



int err_bars(){
	
	int nbins = 1800;

	TFile *f = new TFile("full_result.root");
	TH1D *h = (TH1D *)f->Get("e_live_time");
	TH1D *h2 = (TH1D *)f->Get("e_tof_stilbene");
	Th1D *h3 = (TH1D *)f->Get("e_tof_fission");
	
	Double_t errs[nbins];
	Double_t errs2[nbins];
	Double_t errs3[nbins];

	for(int i = 1; i<=nbins; i++){
		errs[i] = h->GetBinError(i);
		errs2[i] = h2->GetBinError(i);
		errs3[i] = h3->GetBinError(i);
	}

	return 0;
}

double *prop_err(double *err, double *err2, int nbins){
	double tot[nbins];
	for(int i = 0; i<nbins; i++){
		tot[i] = sqrt(err[i] + err2[i]);
	}
	return tot;
}
