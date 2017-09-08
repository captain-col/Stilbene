#ifndef __QMUONBUFFER_HH_
#define __QMUONBUFFER_HH_
#include <stdio.h>
#include <sys/time.h>
#include <sys/shm.h>
#include <sys/ipc.h>
#include <sys/times.h>
#define NCHANNELSQDC 0
#define NCHANNELSTDC 8
//maximum number of hits per event on TDC
#define MAXHITTDC 32*(NCHANNELSTDC-1)
#define MAXEVENTSPERFILE 5000 //was 5000
#define NEVENTSSHM 1000  //was 1000
#define READSIZE  2+4*32*NCHANNELSTDC

int n_active_channels_qdc = 0;
int n_active_channels_tdc = 8;

struct datum_qdc_t{
  short channel;
  int value;
  short underThrOverFlow; // 0 -> good datum
  // 1 -> overflow
  // 2 -> underthreshold
  // 4 -> invalid datum (trigger not accepted)
  short validDatum;   // 1 -> valid datum
  //0 -> not valid datum written in the MEB
  int ev_count;
};


struct datum_tdc_t{
  int value;
  short channel;
  //incremented by 1 every 32 (or 16 words) if I remember correctly
  short bunch_id;  
  //error flag for the hit (from TDC board)                                                                                 
  short error_flag;
  
};

struct event_t{
  //start time of the run
  clock_t startTimeSec;
  //lasts for ~ 4 days
  //time elapsed since the start of the run in tenths of ms
  time_t clockTimeTenthsMs;
  //pointer to the TDC data
  datum_tdc_t time[MAXHITTDC];
  //pointer to the QDC data
 // datum_qdc_t charge[NCHANNELSQDC];
  //number of hits on TDC for this event
  size_t n_hits_tdc;
  //error flag for the event (from TDC board)
  size_t error_flag_tdc;
  //event number as read from the QDC board
//  size_t eventNumberQDC;
  //event number as read from the TDC board
  size_t eventNumberTDC;
  //number of memorised channels
//  size_t n_chans_qdc;
 
};
//structure of data in the shared memory
//irrelevant for data analysis
struct shmbuffer_t{
	unsigned long long p_read;
	unsigned long long p_write;
	event_t data[NEVENTSSHM];
};

//#define READSIZE sizeof(event_t)
void DumpEvent(event_t&);
void DumpEventNb(event_t&);
//used for debugging purposes
/*void DumpEvent(event_t &ev){
 printf("ev = %d\n",ev.eventNumberQDC);
  printf("n_active_channels_qdc = %d\n",n_active_channels_qdc);
  printf("StartTimeSec = %lu, in hex = %x, size = %lu\n",ev.startTimeSec,ev.startTimeSec, sizeof(clock_t));
  printf("ClockTimeTenthsMs = %lu, in hex = %x, size = %lu\n",ev.clockTimeTenthsMs,ev.clockTimeTenthsMs, sizeof(time_t));
  printf("Printing time info ***************************************\n");
  for(int j = 0; j< ev.n_hits_tdc; j++){
    printf("Value = %d, in hex = %x, size = %lu, j=%d \n",(ev.time[j]).value, (ev.time[j]).value  , sizeof(int), j);
    printf("Channel = %d, in hex = %x, size = %lu, j=%d\n",(ev.time[j]).channel, (ev.time[j]).channel  , sizeof(int),j);
    printf("Bunch = %d, in hex = %x, size = %lu, j=%d\n",(ev.time[j]).bunch_id, (ev.time[j]).bunch_id, sizeof(int),j);
    printf("Error = %d, in hex = %x, size = %lu, j=%d\n",(ev.time[j]).error_flag, (ev.time[j]).error_flag, sizeof(int),j);
   
  }*/
  
  /*printf("Printing charge info ***************************************\n");
  for(int i = 0; i < n_active_channels_qdc; i++){
    printf("Q.ch = %hu, in hex = %x, size = %lu, i = %d\n",(ev.charge[i]).channel,(ev.charge[i]).channel, sizeof(short), i);
    printf("Q.underThrOverFlow = %hu, in hex = %x, size = %lu, i = %d\n",(ev.charge[i]).underThrOverFlow,(ev.charge[i]).underThrOverFlow, 
	   sizeof(short), i);
    printf("Q.validDatum = %hu, in hex = %x, size = %lu, i = %d\n",(ev.charge[i]).validDatum,(ev.charge[i]).validDatum, 
	   sizeof(short), i);
  }
  
 printf("EventNumberQDC = %lu,in hex = %x, size = %lu\n",ev.eventNumberQDC,ev.eventNumberQDC , sizeof(size_t));
printf("hits_on_tdc = %lu,in hex = %x, size = %lu\n",ev.n_hits_tdc,ev.n_hits_tdc, sizeof(size_t));
printf("tdc_error_flag = %lu,in hex = %x, size = %lu\n",ev.error_flag_tdc,ev.error_flag_tdc, sizeof(size_t));
printf("******************************************************************\n \n");
}*/

void DumpEventNb(event_t &ev){
  if((ev.eventNumberTDC)%1000 == 0)
    printf("EventNumberTDC = %lu\n",ev.eventNumberTDC);
}
#endif

