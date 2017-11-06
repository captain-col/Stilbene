#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <TString.h>
#include <TH1I.h>
#include <TFile.h>
#include <vector>
#include <TH1F.h>
#include <TList.h>
#include <TH2D.h>
#include <math.h>
#include <TAxis.h>
#include "/home/csc2168/particle/stilbene/MuonBuffer_test.hh"
#include <map>


event_t event_buffer;
int analyzeFlux(int filenumberstart, int filenumberend, bool bl);
void FillHistos(event_t& ev);
int tdc_ev_count = -1;
int b = 1;
int nbins = 50;
int binl = 1800/nbins;
long startRun;
long startCurrentRun;
char path[70] = "/home/csc2168/particle/stilbene/";
double lsb = 0.025;
double TOF2Energy(double TOF);
double scintThresh = 50.;

long timeBin = 10000; //1 s in tenths of ms
int nBins =0;
int evCounter;
int hitCounter;
int tdc_ev_nb = 0;
double time_gamma_bacon = 744.6;
//double time_gamma = 774.;
//double time_gamma = 780.;
double time_gamma = -500;
long endRun;
int above, below;
double runTime, deltaT;

TH1I* h1; // hits on channel 1 stilbene
TH1I* h2; // hits on channel 2 trigger
TH1I* h3; // hits on channel 3 fission
TH1D* h4; // channel 1 distribution
TH1D* h5; // channel 2 distribution
TH1D* h6; // channel 3 distribution
TH1D* h7; // trigger-subtracted TOF of fission chamber
TH1D* h8; // trigger-subtracted TOF of stilbene
TH1D* h9; // bunch_id distribution (macro pulse)
TH1D* h10; // Kinetic Energy of Neutrons
TH1D* h11; // Delta same event hits Ch2
TH1D* h12; // Delta same event hits Ch3
TH1D* h13; // Delta_T
TH1D* h14; // Live Time
TH1D* h15; // Time difference between stilbene hits

TAxis* live_axis;// live time axis


int analyzeFlux(int filenumberstart, int filenumberend, bool bl = false)
{	
	b = 1;
	evCounter = 0;
	char filename[90];
	tdc_ev_nb = 0;
	above = 0;
	below = 0; 
	runTime = 0.;  
	deltaT = 0.;
	h1 = new TH1I("hit_num_ch1","Number of hits on channel 1 Stilbene",10,0,10);
	h2 = new TH1I("hit_num_ch2","Number of hits on channel 2 Trigger",10,0,10);
	h3 = new TH1I("hit_num_ch3","Number of hits on channel 3 Fission",10,0,10);
	h4 = new TH1D("hit_ch1","Hit distribution for channel 1 Stilbene",1800, 0,1800);
	h4->GetXaxis()->SetName("ns");
	h4->GetXaxis()->SetTitle("ns");
	h5 = new TH1D("hit_ch2","Hit distrubtion on channel 2 Trigger",1800,0,1800);
	h5->GetXaxis()->SetName("ns");
	h5->GetXaxis()->SetTitle("ns");	
	h6 = new TH1D("hit_ch3","Hit distribution for channel 3 Fission",1800,0,1800);
	h6->GetXaxis()->SetName("ns");
	h6->GetXaxis()->SetTitle("ns");
	h7 = new TH1D("tof_stilbene","TOF Distribution for Stilbene",1800,-800,1000);
	h8 = new TH1D("tof_fission","TOF Time Distribution for Fission Chamber",1800,-800,1000);
	h9 = new TH1D("bunch_id", "Distribution of Bunch ID's of Triggers", 1000, 0, 5000);
	h10 = new TH1D("Ekin","Neutron kinetic energy",400,0.,800);
	h11 = new TH1D("delta_ch1","Time difference for same event hits Stilbene",100,0,200);
	h12 = new TH1D("delta_ch3","Time difference for same event hits Fission",100,0,200);
	h13 = new TH1D("deltaT","Number of Triggers over Absolute Time Elapsed", 500,0,2500);
	h13->GetXaxis()->SetTitle("Ms Elapsed Since File Beginning");
	h14 = new TH1D("live_time", "Live Time", 1800, -800, 1000); // 180 ns dead time
	live_axis = h14->GetXaxis();
	h15 = new TH1D("stil_hit_dif", "Time Difference Between Stilbene Hits", 1800, 0, 1800);
	
	TList *l = new TList();
	l->Add(h1);
	l->Add(h2);
	l->Add(h3);
	l->Add(h4);
	l->Add(h5);
	l->Add(h6);
	l->Add(h7);
	l->Add(h8);
	l->Add(h9);
	l->Add(h10);
	l->Add(h11);
	l->Add(h12);
	l->Add(h13);
	l->Add(h14);
	l->Add(h15);

	char BLFName[100] = "/home/csc2168/particle/stilbene/errorEventFile.txt";
	FILE* blackListFile = fopen(BLFName,"r");
	std::map<int,int> BLmap;
	std::map<int,int>::iterator blit;
	int maxEvent = 100001;
	int filenb, evnb;
	while(fscanf(blackListFile,"%d %d\n",&filenb, &evnb)!=EOF){
	BLmap.insert ( std::pair<int,int>(filenb,evnb) );
	} 
	for(int filenumber = filenumberstart;filenumber<= filenumberend; filenumber++){ 
		sprintf(filename,"%s/outFile_%d.dat",path,filenumber);
		FILE*  outFile = fopen(filename,"r");
		if(!outFile){
			printf(" file %s not found. Exiting...\n",filename);
			exit(0);
		}
		if(bl == true && BLmap.find(filenumber) != BLmap.end()){
		maxEvent = BLmap[filenumber];
		}
		else maxEvent = 100001;
		printf("Opening file %s\n",filename);
		int evCt = 0;
		while(fread(&event_buffer,sizeof(event_buffer),1,outFile) == 1){
			if(evCt < maxEvent){	
			FillHistos(event_buffer);
			}
			if(evCt%1000 == 0){ 
				printf("Event number %d, runTime %f s\n",evCt, deltaT);
			}
			evCt++;
		}
		fclose(outFile);
	}

	runTime += deltaT;
	printf("events processed: %d, duration %f s \n",evCounter, runTime);
	printf("Number of neutrons with energy <%lg MeV %d, > 800 MeV %d\n", scintThresh, below, above);
	TFile* f = new TFile(Form("histos_%d_%d.root",filenumberstart, filenumberend),"recreate");
	TIter iter(l);
	TH1D* h;
	while ((h = (TH1D*)iter.Next())){
		h->Write();
	}
	f->Close();
	return 0;
}

void FillHistos(event_t &ev){
	int bin_num = 0;
	int hit_num = 0;

	tdc_ev_nb = ev.eventNumberTDC;
	if(b == 1 || ev.startTimeSec != startCurrentRun){
		startRun = ev.clockTimeTenthsMs;
		startCurrentRun = ev.startTimeSec;
		b = 0;
		runTime += deltaT;
		printf("start time for this run: %ld \n", startCurrentRun);
	}
	else{
		deltaT = (ev.clockTimeTenthsMs - startRun)/10.;
	}
	evCounter++;

	int n_hits_ch1 = 0;
	int n_hits_ch2 = 0;
	int n_hits_ch3 = 0;

	h13->Fill(deltaT);

	std::vector<double> times1; 
	std::vector<double> times2; 
	std::vector<double> times3; 

	for(int j = 0; j< ev.n_hits_tdc; j++){
		if((ev.time[j]).channel == 1){
			times1.push_back((ev.time[j]).value);
			n_hits_ch1++;
		}
		else if ((ev.time[j]).channel == 2){
			times2.push_back((ev.time[j]).value);
			n_hits_ch2++;
			h9->Fill(ev.time[j].bunch_id);
		}
		else if ((ev.time[j]).channel == 3){
			times3.push_back((ev.time[j]).value);
			n_hits_ch3++;
		}
	}

	if(n_hits_ch1 >0 && n_hits_ch2 >0 && n_hits_ch3 > 0){

		for(int i =0; i < n_hits_ch1; i++){
			
			double time_neutron = (times1[i]-times2[0])*lsb;	
			double deltaTOF = time_gamma - time_neutron;
			double Ekin = TOF2Energy(deltaTOF);	
			if(Ekin <=800. && Ekin >=scintThresh){
				h10->Fill(Ekin);
			}
			else if (Ekin > 800.){
				above++;
			}
			else if(Ekin < scintThresh){
				below++;
			}
		}
		
		std::vector<double> delta_t;
		std::vector<double> delta_t2;
		
		for(int i =0; i < n_hits_ch2; i++){
            h5->Fill(times2[i]*lsb);
        }
		for(int i =0; i < n_hits_ch1; i++){
			delta_t.push_back((times1[i]-times2[0])*lsb);
			h4->Fill(times1[i]*lsb);
			h7->Fill((times1[i]-times2[0])*lsb);
		}
		for(int i =0; i < n_hits_ch3; i++){
			delta_t2.push_back((times3[i]-times2[0])*lsb);
			h6->Fill(times3[i]*lsb);
			h8->Fill((times3[i]-times2[0])*lsb);
		}
		for(int k = 0; k < delta_t.size() - 1; k++){
			h11->Fill(- delta_t[k+1] + delta_t[k]);
		}
		for(int k = 0; k < delta_t2.size() - 1; k++){
			h12->Fill(- delta_t2[k+1] + delta_t2[k]);
		}
	}
	
	else if(n_hits_ch1 >0 && n_hits_ch2 >0){
		  
		for(int i =0; i < n_hits_ch1; i++){
			
			double time_neutron = (times1[i]-times2[0])*lsb;	
			double deltaTOF = time_gamma - time_neutron;
			double Ekin = TOF2Energy(deltaTOF);	
			if(Ekin <=800. && Ekin >=scintThresh){
				h10->Fill(Ekin);
			}
			else if (Ekin > 800.){
				above++;
			}
			else if(Ekin < scintThresh){
				below++;
			}
		}
		std::vector<double> delta_t;
		for(int i =0; i < n_hits_ch2; i++){
            h5->Fill(times2[i]*lsb);
        }
		for(int i =0; i < n_hits_ch1; i++){
			delta_t.push_back((times1[i]-times2[0])*lsb);
			h7->Fill((times1[i]-times2[0])*lsb);
			h4->Fill(times1[i]*lsb);
		}
		for(int k = 0; k < delta_t.size() - 1; k++){
			h11->Fill(- delta_t[k+1] + delta_t[k]);
		}
	
	}
	
	else if(n_hits_ch3 >0 && n_hits_ch2 >0){
		  
		std::vector<double> delta_t2;
		for(int i =0; i < n_hits_ch2; i++){
            h5->Fill(times2[i]*lsb);
        }
		for(int i =0; i < n_hits_ch3; i++){
			delta_t2.push_back((times3[i]-times2[0])*lsb);
			h6->Fill(times3[i]*lsb);
			h8->Fill((times3[i]-times2[0])*lsb);
		}
		for(int k = 0; k < delta_t2.size() - 1; k++){
			h12->Fill(- delta_t2[k+1] + delta_t2[k]);
		}
	}
	
	else if(n_hits_ch2 > 0){
	for(int i =0; i < n_hits_ch2; i++){
          h5->Fill(times2[i]*lsb);
      }
	}
	
	//dead time correction
	if(n_hits_ch1 > 0){
		hit_num = 0;
		bin_num = live_axis->FindBin((times1[hit_num] - times2[0])*lsb);
		for(int i = 0; i < nbins; i++){
			if(i == bin_num){
				h14->Fill(i*binl - 800);
				i = i + 5; //5 * 36ns bins = 180ns
				hit_num++;
				if(hit_num < n_hits_ch1){
					bin_num = live_axis->FindBin((times1[hit_num] - times2[0])*lsb);
					if(bin_num < i){
						printf("Error! Second hit within dead time!\n");
					}
				}
			}
			else{
				h14->Fill(i*binl - 800);
			}
		}
	}
	else{
		for(int i = 0; i < nbins; i++){
			h14->Fill(i*binl - 800);
		}
	}
	
	if(n_hits_ch1 > 1){
		for(int i = 0; i < n_hits_ch1 -1; i++){
			h15->Fill((times1[i+1] - times1[i])*lsb);
		}
	}

	times1.clear();
	times2.clear();
	times3.clear();

	h1->Fill(n_hits_ch1);
	h2->Fill(n_hits_ch2);
	h3->Fill(n_hits_ch3);

}


double TOF2Energy(double deltaTOF){
	deltaTOF=deltaTOF/1.E9;  //in s
	double c = 299792458; //in m/s
	double l = 15; //target to detector distance in m
	double m = 939.5/(c*c);  //neutron mass in MeV
	double TOFgamma = l/c;
	double result = 900;
	if(deltaTOF >=0){
		double TOF = deltaTOF+TOFgamma;
		double restE = m*c*c;
		//	printf("RestE = %f, TOF = %e, tc - l = %f\n",restE, TOF, TOF*c - l);
		result = ((TOF*c/sqrt((TOF*c*TOF*c)-(l*l)))-1)*restE;
	}

	return result;
}


