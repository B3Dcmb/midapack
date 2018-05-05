/** @file   sky_scan.c
    @brief  write the direction (cartesian coordinnates) of observation into files 

    Typically the pointing of observation instrument (a telescope in orbit or on earth) is defined by many kinematic parameters.
    To make it easy, I just defined poiting stream with a set of three rotations, and an angle. As telescope was a giroscope. 
    Rotations represents 3 successives motions which are usually : 
     - sun revolution (1 year),
     - earth or L2 revolution (24 or 1hours),
     - spacecraft spin.
    Actually, observation mission can have other specific rotation.
    Anyway, with 3 rotations you should be able to simulates whatever sky scan strategy. 
    Last parameters is the telescope direction on the spacecraft. 
    <br>
    Input : few parameters can be passed at run time as rotations, sampling period, and so on...
    There are also default values. To simulate realistic mission data, you can use the following default values :   
    -w1 = 132.0 rph (2.2 min per period) equal satellite spin.
    -w2 = 1.0 rph is revolution speed in L2 orbit,
    -w3 = 0.000114155 rph (1 per year) sun revolution, 
    <br>
    Output : pointing stream is a list of successive pointing directions. You can specify a filename to write those stream.
    Number of data per file and number of files can be passed as options.
    @author Pierre Cargemel
    @date   December 2011 
 
    Bugs fixed (relevant to long timestreams) - rs
    @date February 2015*/

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <chealpix.h>
#include "midapack.h"


//extern void __mapoutput_MOD_map2gif(int* nside, double* map, int* filenamelen,char* output_file);

extern char *optarg;

void usage(){
  fprintf(stderr, "\n*#  Usage: sky_stream -f [filename]");
  fprintf(stderr, "\n*#  -a [value] opening angle (rad)");
  fprintf(stderr, "\n*#  -b [value] angle swept by ring center (rad)");
  fprintf(stderr, "\n*#  -n [value] number of samples per process");
  fprintf(stderr, "\n*#  -m [value] healpix resolution");
  fprintf(stderr, "\n*#  -h print hitmap");
  fprintf(stderr, "\n*#  -i print matrix info");
  fprintf(stderr, "\n*#  -p print partition map");
  fprintf(stderr, "\n*#  -o [filename]");
}



int main(int argc, char *argv[]){

  double timer[6], tmax[6];
  Mat A;                                  //matrix
  double *y1;                               //resulting values vectors
  double *a1, *x1;                                 //parameter vectors

  int i, j, err, iter, niter=25;
  char ch;
  char *bn;
  char *endp;
  char fn[100];
  FILE *stream;
  int filenamelen;
  int info = 0;
  int part = 0;
  int hit = 0;
  int rank, size;             
  double Wz0, Wu, Wz;			//vitesse de rotation autour de z0, u, et w en révolutuion par heure
  double psi, teta, phi;	 	//angles d'Euler
  int    nt = 50000;			//nombre de mesures par proc
  double t  = 0.0;     			//temps initiale sur le proc
  int nb=1;    			        //nombre de fichiers à écrire 
  double x0[3], u, v, w, x[3];		//coordonnées
  double R0_v[3][3], Rv_w[3][3], Rw_[3][3];   //matrices de rotation de R0 dans Rv, Rv dans Rw et Rw dans R.
  double N;            			//norme
  double max_sweep;			//arc balayé
  int nside = 64;                                                             
  long npix, ipring;
  int* pix;
  int once, twice;
  double *hitmap, *ghitmap;
  double *partmap, *gpartmap;
  double *cover, *gcover;
  char coordsys[1];      
  coordsys[0]='C';      
  int flag=2;
  
  double r=5000;
  double alpha=30.0;	//"angle d'ouverture"
  double beta=360.0;	//"angle balayage à l équateur"


  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  //-----------------------------------options-----------------------
  while ((ch = getopt( argc, argv, "a:b:r:n:m:o:f:ihp")) != EOF) {
    switch(ch) {
      case 'a':
         alpha = strtod(optarg, &endp);
        if (errno == ERANGE){ 
          usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      break;
      case 'b':
         beta = strtod(optarg, &endp);
        if (errno == ERANGE){
          usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      break;
      case 'r':
         r = strtod(optarg, &endp);
        if (errno == ERANGE){
          usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      break;
      case 'n':
        nt = strtoul(optarg, 0, 10);
        if ((errno == ERANGE) || (errno != 0 && nt <= 0)) {
          usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      break;
      case 'm':
        nside = strtoul(optarg, 0, 10);
        if ((errno == ERANGE) || (errno != 0 && nside <= 0)) {
          usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      break;
      case 'f':
 	flag = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && flag < 0)) {
	  usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
	}
      break;
      case 'o':
        bn = strdup(optarg);
      break;
      case 'h':
        hit=1;
      break;
      case 'i':
        info=1;
      break;
      case 'p':
        part=1;
      break;
    }
  }





  //----------------------------------------------scanning----------------------
  x0[0]= cos(2*M_PI*alpha/360.0);	//pointage instrument   
  x0[1]= sin(2*M_PI*alpha/360.0);	//
  x0[2]= 0.0;		//

  Wz0 = 0.0;			//vitesse  
  Wu = 1.0/r;			//vitesse
  Wz = beta/(360.0*size*nt);	//vitesse 
  npix = nside2npix(nside);

  if(rank==0){
    printf("\n========Parameters============");
    printf("\nalpha beta (degrees) :  %.2lf %.2lf", alpha, beta); 
    printf("\ntotal samples        :  %ld", nt*size);
    printf("\nper ring per process :  %.2lf %d", r, nt); 
    printf("\nnside npix           :  %d %d", nside, npix);
    printf("\nfilename             :  %s", bn);
    printf("\nrank size            :  %d %d", rank, size);
    printf("\nflag                 :  %d\n", flag);
  }
  
  timer[0] = MPI_Wtime();
 
  teta = 1.0*rank*nt-(int)((1.0*rank*nt)*Wu)*r;  // rs 02/04/2015 - remove irrelevant multiples of 2*PI                                         
  teta *= 2*M_PI*Wu;

  phi = 2*M_PI*(rank*beta/360.0)/size;   // rs 02/02/2015 - this is the global rotation around z axis for proc rank 

  pix= (int*) malloc((int64_t)(nt)*sizeof(int));  
  for(i=0; i< nt; i++){

    psi = 2*M_PI*Wz0*t;                                    //compute angular rotation   - N.B. psi not corrected for as is 0 here - rs 02/04/2015                 
    teta += 2*M_PI*Wu;
    phi += 2*M_PI*Wz;

    //matrices de rotation suivant angles d'euler
    //rotation du repère Rv(u,v,z0) / R0(x0,y0,z0 autour de z0 d'angle psi(t)
    R0_v[0][0] = cos(psi);     R0_v[0][1] = sin(psi);    R0_v[0][2] = 0.0;
    R0_v[1][0] = -sin(psi);    R0_v[1][1] = cos(psi);    R0_v[1][2] = 0.0; 
    R0_v[2][0] = 0.0;          R0_v[2][1] = 0.0;         R0_v[2][2] = 1.0;

    u = R0_v[0][0]*x0[0] + R0_v[0][1]*x0[1];
    v = R0_v[1][0]*x0[0] + R0_v[1][1]*x0[1];

    //rotation du repère Rw(u,w,z) / R0(u,v,z0) autour de u d'angle teta(t)
    Rv_w[0][0] = 1.0;          Rv_w[0][1] = 0.0;         Rv_w[0][2] = 0.0;
    Rv_w[1][0] = 0.0;          Rv_w[1][1] = cos(teta);   Rv_w[1][2] = sin(teta);
    Rv_w[2][0] = 0.0;          Rv_w[2][1] = -sin(teta);  Rv_w[2][2] = cos(teta);

    w    = Rv_w[1][1]*v + Rv_w[1][2]*x0[2];
    x[2] = Rv_w[2][1]*v + Rv_w[2][2]*x0[2];
      
    //rotation du repère R(x,y,z) / Rw(u,w,z) autour de z d'angle phi(t)
    Rw_[0][0] = cos(phi);      Rw_[0][1] = sin(phi);    Rw_[0][2] = 0.0;
    Rw_[1][0] = -sin(phi);      Rw_[1][1] = cos(phi);     Rw_[1][2] = 0.0; 
    Rw_[2][0] = 0.0;           Rw_[2][1] = 0.0;          Rw_[2][2] = 1.0;

    x[0] = Rw_[0][0]*u + Rw_[0][1]*w;
    x[1] = Rw_[1][0]*u + Rw_[1][1]*w;
    
    vec2pix_ring(nside, x, &ipring);
    pix[i]= (int) ipring;
  } 

  timer[0] = MPI_Wtime() - timer[0];



  //------------------------------------------matrix operations----------
  
  MatSetIndices(&A, nt, 1, pix);

  a1 = (double *) malloc((int64_t)(nt)* sizeof(double));          //< allocate vector values  according partitioning ant init to 0.0
  for(i=0; i<nt; i++)   
    a1[i]=1.0;
  
  MatSetValues(&A, nt, 1, a1);
  
  timer[1] = MPI_Wtime();
  MatLocalShape(&A, 3);
  timer[1] = MPI_Wtime() - timer[1];
  
  x1 = (double *) malloc(A.lcount *sizeof(double));          //< allocate full vector values  on each proc and init to 0.0 
  y1 = (double *) malloc(nt* sizeof(double));          //< allocate vector values  according partitioning ant init to 0.0

  for(iter=0;iter<niter;iter++){

    if(rank==0)
      printf("iter %d\n",iter);
  
    for(j=0; j<A.lcount; j++)   
      x1[j]=1.0;
  
    MPI_Barrier(MPI_COMM_WORLD);
  
    timer[2] = MPI_Wtime();
    MatComShape(&A, flag, MPI_COMM_WORLD);
    timer[2] = MPI_Wtime() - timer[2];
  
    timer[3] = MPI_Wtime();
    MatVecProd(&A, x1, y1, 0);
    timer[3] = MPI_Wtime() - timer[3];
  
    MPI_Barrier(MPI_COMM_WORLD);                               
  
    timer[4] = MPI_Wtime();
    for(i=0; i < A.lcount; i++)				//refresh vector
      x1[i]=0.0;						//
  
    int e=0;
    for(i=0; i< A.m*A.nnz; i+=A.nnz){			//local transform reduce          
      for(j=0; j< A.nnz; j++){
         x1[A.indices[i+j]] += A.values[i+j] * y1[e];	//
      }							//
      e++;						//
    }
  							//
    timer[4] = MPI_Wtime() - timer[4];
  
    timer[5] = MPI_Wtime();
    greedyreduce(&A, x1);					//global reduce
    timer[5] = MPI_Wtime() - timer[5];

    if(iter<niter-1)
      MatReset(&A); 


    MPI_Reduce(timer, tmax, 6, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if(rank==0){
      printf("\n========Timings==============");
      printf("\npre %lf", tmax[0]);
      printf("\nlocal %lf", tmax[1]);
      printf("\ncomm %lf", tmax[2]);
      printf("\nmatvec %lf", tmax[3]);
      printf("\nlocalreduce %lf", tmax[4]);
      printf("\ngreedyreduce %lf\n", tmax[5]);
    }

  }
  //TrMatVecProd(&A, y1, x1, 0);
  free(y1);




#if 0
  //--------------------------------------------------statistics
  if(rank==0)
    printf("\n========Stats==============");
  if(rank==0)
    printf("\ninfo: ");
  timer[0] = MPI_Wtime();
  MatInfo(&A,info, bn);					//write communication matrix
  timer[0] = MPI_Wtime() - timer[0];
  

  if(hit==1){						//build hitmap
    if(rank==0)
      printf("\nhit: ");
    timer[1] = MPI_Wtime();
    hitmap  = (double *) calloc(npix, sizeof(double));	//gather map on proc 0
    ghitmap  = (double *) calloc(npix, sizeof(double));	//gather map on proc 0
    for(i=0; i<A.lcount; i++)				//
      hitmap[A.lindices[i]]= x1[i];		    	//
    MPI_Reduce(hitmap, ghitmap, npix, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    timer[1] = MPI_Wtime() - timer[1];

    timer[2] = MPI_Wtime();
    /*if(rank==0){					//write map 
      sprintf(fn, "hit_%s_%d.gif", bn, size);		//
      filenamelen = strlen(fn);
      __mapoutput_MOD_map2gif((int*)(&nside), ghitmap, &filenamelen, fn);
    }*/
    free(hitmap);					//
    free(ghitmap);					//
    timer[2] = MPI_Wtime() - timer[2];
  }

  
  for(j=0; j < A.lcount; j++)				//cover stats  
    x1[j]=1.0;						//
    MPI_Barrier(MPI_COMM_WORLD);                               
    timer[3] = MPI_Wtime();				//
    greedyreduce(&A, x1);				//
    timer[3] = MPI_Wtime() - timer[3];			//
 
    timer[4] = MPI_Wtime();				//
    cover  = (double *) calloc(size, sizeof(double));	//compute cover
    gcover  = (double *) calloc(size, sizeof(double));	//compute cover
    for(i=1; i<size; i++){				//
      for(j=0; j<A.lcount; j++){			//
        if(x1[j]-i==0.0)				//
          cover[i] += 1.0;				//
      }
    } 
    MPI_Reduce(cover, gcover, size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(rank==0){					//write map 
      printf("\ncover ", gcover[i]/i);			//
      once=0.0;						//
      for(i=1; i<size; i++){
        once += gcover[i]/i;				//
        printf(" %.0lf", gcover[i]/i);			//
      }	
      printf("\nonce %d", once);
      twice= once - gcover[1];				//
      printf("\ntwice %d", twice);	//
    }							//
    free(cover);					//
    free(gcover);					//
    MPI_Bcast(&twice, 1, MPI_INT, 0, MPI_COMM_WORLD);
    timer[4] = MPI_Wtime() - timer[4];			//

    MPI_Barrier(MPI_COMM_WORLD);                               
    partmap  = (double *) calloc(twice, sizeof(double));	//gather map on proc 0
    gpartmap  = (double *) calloc(twice, sizeof(double));//	
    timer[5] = MPI_Wtime();				//
    MPI_Allreduce(partmap, gpartmap, twice, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    timer[5] = MPI_Wtime() - timer[5];			//

    free(partmap);					//
    free(gpartmap);					//
  /*if(part==1){						//partition map
    if(rank==0)
      printf("\npart: ");
    timer[5] = MPI_Wtime();
    partmap  = (double *) calloc(npix, sizeof(double));	//gather map on proc 0
    gpartmap  = (double *) calloc(npix, sizeof(double));//	
    for(i=0; i<A.lcount; i++)				//
      partmap[A.lindices[i]]= x1[i];			//
    MPI_Reduce(partmap, gpartmap, npix, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
                    
    if(rank==0){					//write map 
      sprintf(fn, "part_%s_%d.gif", bn, size);		//
      filenamelen = strlen(fn);
      __mapoutput_MOD_map2gif((int*)(&nside), gpartmap, &filenamelen, fn);
    }							//
    timer[5] = MPI_Wtime() - timer[5];
  }*/ 


  MPI_Reduce(timer, tmax, 6, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if(rank==0){
    if(info==1)
      printf("\ninfo %lf", tmax[0]);
    if(hit==1){
      printf("\nhit1 %lf", tmax[1]);
      printf("\nhit2 %lf", tmax[2]);
    }
    printf("\ngreedyreduce %lf", tmax[3]);
    printf("\ncover %lf", tmax[4]);
    printf("\nallreduce %lf\n", tmax[5]);
  }
#endif

  MatFree(&A); 
  free(pix);
  free(x1);						//
  free(a1);
  MPI_Finalize();     
  return 0;
}


