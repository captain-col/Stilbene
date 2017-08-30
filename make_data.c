#include <stdlib.h>
#include <stdio.h>
#include "/home/daquser/Desktop/DAQ2013/MuonBuffer_test.hh"


int main(int argc, char* argv[])
{
  if(argc !=4){
    printf("Error. Need beginning and ending file numbers or hard drive #\n");
    exit(1);
  }

  int begin_filen =  atoi(argv[1]);
  int end_filen = atoi(argv[2]);
  int hard_drive = atoi(argv[3]);
  
  int counter;
  FILE *ptr_myfile;
  FILE *data_file;
  struct event_t my_record;
  short nc_flag = 0;
  int last = 0;
  char filename[50]; 
  
  data_file = fopen("data_file.txt", "wb");
  
  for(int filen = begin_filen; filen <= end_filen; filen++){
 
    if(hard_drive == 1){
      sprintf(filename, "/media/X/2017_WNRstilbene/outFile_%d.dat", filen);
    }
    else{
      sprintf(filename, "/media/LACIE SHARE/2017_WNRstilbene/outFile_%d.dat", filen); 
    }
    

    ptr_myfile = fopen(filename, "rb");
    if (!ptr_myfile)
      {
	printf("Unable to open file! Correct hard drive (command line parameter 3)? \n");
	return 1;
      }
  for( counter=1; counter <= MAXEVENTSPERFILE; counter++)
    {
      fread(&my_record,sizeof(struct event_t),1,ptr_myfile);
      //   printf("Nhits for %d is %d\n", counter, my_record.n_hits_tdc);
       for(int i = 0; i< my_record.n_hits_tdc; i++){
	  fprintf(data_file, "%d ", my_record.time[i].channel);
	  fprintf(data_file, "%d\n", my_record.time[i].value);
	  if(i == my_record.n_hits_tdc -1){
	    fprintf(data_file, "%d\n", -3);
	    fprintf(data_file, "%d\n", -3); //2 for parity
	  }
      }
     }
       fclose(ptr_myfile);
  }	
  fclose(data_file);
  return 0;
}
