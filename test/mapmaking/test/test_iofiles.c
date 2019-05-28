//iofiles
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <mpi.h>
//#include <time.h>
#include <string.h>
//#include <midapack.h>


int main (int argc, char const *argv[]) {

  int i;

  int block_size = pow(2,20);//+200000;//163936;
  int part_id   =  0;  // could be rank of processor

    printf(" block_size = %d\n", block_size );


  int *point_data;
  point_data  = (int *) malloc(block_size*sizeof(int));
  double *signal;
  signal = (double *) malloc(block_size*sizeof(double));
 
//  unsigned int point_data[block_size];
//  double signal[block_size];


  ioReadfile(block_size, part_id, point_data, signal);


  for( i = 0; i < 10; ++i) {
    printf(" point_data[%d] = %u\n", i, point_data[i] );
    printf("signal[%d] = %lf\n", i, signal[i] );
  }


  free(point_data);
  free(signal);

        return 0;
}


int ioReadfile( int block_size, int part_id, unsigned int *point_data, double *signal)
{

        int i;

        char p_vectorFile[256];
        char *p_vectorFileNameBasis = "point_data_0_";

        char s_vectorFile[256];
        char *s_vectorFileNameBasis = "signal_";

//        unsigned int point_data[block_size];
//        double signal[block_size];

//  unsigned int *point_data;
//  point_data  = (unsigned int *) malloc(block_size*sizeof(unsigned int));

//  double *signal;
//  signal  = (double *) malloc(block_size*sizeof(double));


        FILE *fp;

        sprintf(p_vectorFile, "%s%01d.dat", p_vectorFileNameBasis, part_id);
        sprintf(s_vectorFile, "%s%01d.dat", s_vectorFileNameBasis, part_id); // s_vectorFile is now = signal_0.dat

        printf(" Pointing file name: %s\n", p_vectorFile);
        printf("   Signal file name: %s\n", s_vectorFile);


        fp=fopen(p_vectorFile, "rb");
        fread(point_data, sizeof(unsigned int), block_size, fp);

        fp=fopen(s_vectorFile, "rb");
        fread(signal, sizeof(double), block_size, fp);


//        for( i = 0; i < 10; ++i) {
//                point_data_out[i] = point_data[i] ;
        //        signal_out[i] = signal[i] ;
//        }

        return 0;
}



