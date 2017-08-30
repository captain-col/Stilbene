
#include <stdio.h>
#include <sys/times.h>
#include <sys/types.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TFile.h>
#include <vector>
#include <TList.h>
#include <TH2D.h>
#include "MuonBuffer_test.hh"
#include <map>

event_t event_buffer;
void FillHistos(event_t& ev);
int tdc_ev_count = -1;
int qdc_ev_count = -1;
time_t startRun;
clock_t startCurrentRun;

bool begin;
time_t timeBin = 10000; //1 s in tenths of ms
int nBins =0;
int evCounter;
int oldMisalignment = 0;
int qdc_ev_nb = 0;
int tdc_ev_nb = 0;
double time_gamma_bacon = 744.6;
//double time_gamma = 774.;
double time_gamma = 780.;
time_t endRun;
int above, below;
time_t runTime, deltaT;

TH1D* h1;
TH1I* h2;
TH1D* h3;
TH1D* h4;
TH1I* h5;
TH1D* h6;
TH1D* h7;
TH1D* h8;
TH1D* h9;
TH1D* h10;
TH1I* h12;
TH2D* h13;
TH1D* h14;
TH1D* h15;
TH1D* h16;
TH1D* h17;
TH1D* h18;
TH1D* h19;
TH1D* h20;
TH1D* h21;
TH1D* h22;
TH1I* h23;
TH1D* h24;
TH2D* h25;
TH2D* h26;
char path[70] = "/Users/elena/Documents/LDRD/Neutrons/data/2013-11-27_Day";
double lsb = 0.025;
double TOF2Energy(double TOF);
double scintThresh = 50.;
int MakePlots(int filenumberstart, int filenumberend, bool bl = false)
{
	begin = true;		
	evCounter = 0;
	char filename[90];
	qdc_ev_nb = 0;
	tdc_ev_nb = 0;
	oldMisalignment = 0;
	above = 0;
	below = 0;
	runTime = 0.;
	deltaT = 0.;
	h1 = new TH1D("event_rate","Event rate vs time",10000,0,10000);
	h1->GetYaxis()->SetTitle("Hz");
	h1->GetXaxis()->SetTitle("time(s)");
	h2 = new TH1I("hit_distrib1","Probability distribution for the number of hits on PMT1",1000,0,100);
	h5 = new TH1I("hit_distrib2","Probability distribution for the number of hits on PMT2",1000,0,100);
	h12 = new TH1I("hit_distrib_PO","Probability distribution for the number of hits on Pick off",1000,0,100);
	h23 = new TH1I("hit_distrib_CU","Probability distribution for the number of hits on coincidence unit",1000,0,100);
	h3 = new TH1D("hit_time_ave","Hit time distribution (average of PMTS)",1000,0,200);
	h4 = new TH1D("hit_time_diff","Hit time difference distribution",1000,-50,50);
	h4->GetXaxis()->SetName("ns");
	h4->GetXaxis()->SetTitle("ns");
	h3->GetXaxis()->SetName("ns");
	h3->GetXaxis()->SetTitle("ns");
	h6 = new TH1D("Q_ch1","Charge on ch1",1920,0,3841);
	h7 = new TH1D("Q_ch2","Charge on ch2",1920,0,3841);
	h16 = new TH1D("Q_ch1_bacon","Charge on ch1",1920,0,3841);
	h17 = new TH1D("Q_ch2_bacon","Charge on ch2",1920,0,3841);
	h19 = new TH1D("Q_ch1_bacon_neu","Charge on ch1 in neutron time window",1920,0,3841);
	h20 = new TH1D("Q_ch2_bacon_neu","Charge on ch2 in neutron time window",1920,0,3841);
	h15 = new TH1D("Q_ch1_neu","Charge on ch1 in neutron time window",1920,0,3841);
	h14 = new TH1D("Q_ch2_neu","Charge on ch2 in neutron time window",1920,0,3841);
	h8 = new TH1D("t_ch1","Time on PMT1",1000,0,1000);
	h9 = new TH1D("t_ch2","Time on PMT2",1000,0,1000);
	h10 = new TH1D("t_ch0","Time of PO coil wrt coincindence signal",1000,0,1000);
	h24 = new TH1D("delta_times","Time differences between same events hits",1000,0,1000);
	h18 = new TH1D("t_ch0_bacon","Time of PO coil wrt coincindence signal",1000,0,1000);
	h10->GetXaxis()->SetTitle("ns");
	h13 = new TH2D("Q_t","Q vs time",500,0,1000,1920,0,3841);
	h25 = new TH2D("Q_t_bacon","Q vs time for bacon",500,0,1000,1920,0,3841);
	h26 = new TH2D("Q_E_bacon","Q vs energy for bacon",400/2,0.,800.,1920/8,0,3841);
	h26->GetXaxis()->SetTitle("E_{kin}(MeV)");
	h26->GetYaxis()->SetTitle("E_{dep}(ADC)");
	h21 = new TH1D("Ekin","Neutron kinetic energy",400,0.,800);
	h22 = new TH1D("Ekin_bacon","Neutron kinetic energy",200,0.,800);
	h21->GetXaxis()->SetTitle("MeV");
	h22->GetXaxis()->SetTitle("MeV");
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
	l->Add(h12);
	l->Add(h13);
	l->Add(h14);
	l->Add(h15);
	l->Add(h16);
	l->Add(h17);
	l->Add(h18);
	l->Add(h19);
	l->Add(h20);
	l->Add(h21);
	l->Add(h22);
	l->Add(h23);
	l->Add(h24);
	l->Add(h25);
	l->Add(h26);

	char BLFName[100] = "/Users/elena/Documents/LDRD/Neutrons/data/codeNovember/errorEventFile.txt";
	FILE* blackListFile = fopen(BLFName,"r");
	std::map<int,int> BLmap;
	std::map<int,int>::iterator blit;
	int maxEvent = 5001;
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
		else maxEvent = 5001;
		printf("Opening file %s\n",filename);
		int evCt = 0;
		while(fread(&event_buffer,sizeof(event_buffer),1,outFile) == 1){
			if(evCt < maxEvent){
			FillHistos(event_buffer);
			}
		/*	else{
			printf("Skip file %d, event %d\n",filenumber,evCt);
			}
			*/
			if(evCt%1000 == 0){ 
				printf("Event number %d, runTime %d s\n",evCt, deltaT);
			}
			evCt++;
		}
		fclose(outFile);
	}
	runTime += deltaT;
	printf("events processed: %d, duration %d s \n",evCounter, runTime);
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
	qdc_ev_nb = ev.eventNumberQDC;
	tdc_ev_nb = ev.eventNumberTDC;
	int misalignment = qdc_ev_nb - tdc_ev_nb;
	if(misalignment != oldMisalignment){
		printf("qdc_ev_nb = %d, tdc_ev_nb = %d at event %d\n",qdc_ev_nb, tdc_ev_nb,ev.eventNumberQDC);
	}
	if(misalignment != 0)
		printf("misalignment = %d\n",misalignment);
	oldMisalignment = misalignment;
	if(begin == true || ev.startTimeSec != startCurrentRun){
		startRun = ev.clockTimeTenthsMs;
		startCurrentRun = ev.startTimeSec;
		begin = false;
		runTime += deltaT;
		printf("start time for this run: %d \n", startCurrentRun);
	}
	else{
		deltaT = (ev.clockTimeTenthsMs - startRun)/10000.;
	}
	evCounter++;

	int n_hits_ch3 = 0;
	int n_hits_ch2 = 0;
	int n_hits_ch1 = 0;
	int n_hits_ch4 = 0;
	int n_hits_ch6 = 0;
	int n_hits_ch7 = 0;
	int n_hits_ch8 = 0;
	std::vector<double> times1; 
	std::vector<double> times2; 
	std::vector<double> times3; 
	std::vector<double> times4; 
	std::vector<double> times6; 
	std::vector<double> times7; 
	std::vector<double> times8; 
	for(int j = 0; j< ev.n_hits_tdc; j++){

		if((ev.time[j]).channel == 2){
			n_hits_ch2++;
			times2.push_back((ev.time[j]).value);

		}

		else if ((ev.time[j]).channel == 3){
			n_hits_ch3++;
			times3.push_back((ev.time[j]).value);

		}

		else if ((ev.time[j]).channel == 1){
			times1.push_back((ev.time[j]).value);
			n_hits_ch1++;

		}

		else if ((ev.time[j]).channel == 4){
			times4.push_back((ev.time[j]).value);

			n_hits_ch4++;
		}
		else if ((ev.time[j]).channel == 6){
			times6.push_back((ev.time[j]).value);

			n_hits_ch6++;
		}
		else if ((ev.time[j]).channel == 7){
			times7.push_back((ev.time[j]).value);

			n_hits_ch7++;
		}
		else if ((ev.time[j]).channel == 8){
			times8.push_back((ev.time[j]).value);

			n_hits_ch8++;
		}
	}
	double totCharge = 0;
	double charge1 =0;
	double charge2 = 0;
	double charge3 = 0;
	double charge4 = 0;


	for(int i = 0; i <ev.n_chans_qdc; i++){
		if((ev.charge[i]).channel == 0){
			h6->Fill((ev.charge[i]).value);
			totCharge+=(ev.charge[i]).value;
			charge1 =(ev.charge[i]).value;
		}
		else if((ev.charge[i]).channel == 1){
			h7->Fill((ev.charge[i]).value);
			totCharge+= (ev.charge[i]).value;
			charge2 =(ev.charge[i]).value;
		}
		else if((ev.charge[i]).channel == 2){
			h16->Fill((ev.charge[i]).value);
			totCharge+= (ev.charge[i]).value;
			charge3=(ev.charge[i]).value;
		}
		else if((ev.charge[i]).channel == 3){
			h17->Fill((ev.charge[i]).value);
			totCharge+= (ev.charge[i]).value;
			charge4=(ev.charge[i]).value;
		}
	}
	if(n_hits_ch1 >0){
		for(int i = 0; i <ev.n_chans_qdc; i++){
			if((ev.charge[i]).channel == 0){
				h15->Fill((ev.charge[i]).value);    
			}
			else if((ev.charge[i]).channel == 1){
				h14->Fill((ev.charge[i]).value);
			}
			else if((ev.charge[i]).channel == 2){
				h19->Fill((ev.charge[i]).value);
			}
			else if((ev.charge[i]).channel == 3){
				h20->Fill((ev.charge[i]).value);
			}


		}
	}

	if(n_hits_ch1 >0 && n_hits_ch2 >0 && n_hits_ch3 >0 && n_hits_ch4>0){
		for(int i =0; i < 1; i++){

			double time_neutron = (times1[0]-times4[i])*lsb;
			h10->Fill(time_neutron);	
			double deltaTOF = time_gamma - time_neutron;
			double Ekin = TOF2Energy(deltaTOF);	
			if(Ekin <=800. && Ekin >=scintThresh){
				h21->Fill(Ekin);
			}
			else if (Ekin > 800.){
				above++;
			}
			else if(Ekin < scintThresh){
				below++;
			}
			h23->Fill(n_hits_ch4);

			h13->Fill((times1[0]-times4[i])*lsb,charge1);
		}
		std::vector<double> delta_t;
		for(int i =0; i < n_hits_ch4; i++){
			delta_t.push_back((times1[0]-times4[i])*lsb);
		}
		for(int k = 0; k < delta_t.size() - 1; k++){
			h24->Fill(- delta_t[k+1] + delta_t[k]);
		}
		h8->Fill((times1[0]-times2[0])*lsb);
		h9->Fill((times1[0]-times3[0])*lsb);

	}

	if(n_hits_ch1 >0 && n_hits_ch6 >0 && n_hits_ch7 >0 && n_hits_ch8>0){
		for(int i =0; i < 1; i++){

			double time_neutron = (times1[0]-times8[i])*lsb;
			h18->Fill(time_neutron);
			double deltaTOF = time_gamma_bacon - time_neutron;
			double Ekin = TOF2Energy(deltaTOF);	
			if(Ekin <=800. && Ekin >=scintThresh){
				h22->Fill(Ekin);
				h26->Fill(Ekin,charge3);
			}
			else if (Ekin > 800.){
				above++;
			}
			else if(Ekin < scintThresh){
				below++;
			}
			h25->Fill((times1[0]-times8[i])*lsb,charge3);

		}
	}
	times1.clear();
	times2.clear();
	times3.clear();
	times4.clear();

	times6.clear();
	times7.clear();
	times8.clear();

	h2->Fill(n_hits_ch2);
	h5->Fill(n_hits_ch3);
	h12->Fill(n_hits_ch1);



}
double TOF2Energy(double deltaTOF){
	deltaTOF=deltaTOF/1.E9;  //in s
	double c = 299792458; //in m/s
	double l = 25; //target to detector distance in m
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

