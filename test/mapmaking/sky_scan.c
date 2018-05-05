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
    @date   December 2011 */

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <math.h>
#include <mpi.h>
#include <chealpix.h>
#include <midapack.h>


extern void __mapoutput_MOD_map2gif(int* nside, double* map, int* filenamelen,char* output_file);

extern char *optarg;

/** @brief display common available options

    Usage: sky_scan -f [filename] -p [val] -q [val] -r [val] ... */
void usage(){
  fprintf(stderr, "\n*#  Usage: sky_stream -f [filename]");
  fprintf(stderr, "\n*#  -p|q|r [value] rotations");
  fprintf(stderr, "\n*#  -x|y|z [value] pointing direction");
  fprintf(stderr, "\n*#  -d [value] time between successive  measures, hour");
  fprintf(stderr, "\n*#  -n [value] number of samples by file");
  fprintf(stderr, "\n*#  -m [value] ihealpix resolution");
  fprintf(stderr, "\n*#  -h [value] print usage");
}



int main(int argc, char *argv[]){

  double t0, t1, max_t1, sum_t1, min_t1;
  Mat At;                                  //matrix
  int *y1_i;                                 //resulting indices vectors
  double *y1_v;                               //resulting values vectors
  double *x1;                                 //parameter vectors

  int i, err;
  char ch;
  char *bn;
  char *endp;
  char fn[100];
  FILE *stream;

  int rank, size;             
  double Wz0, Wu, Wz;			//vitesse de rotation autour de z0, u, et w en révolutuion par heure
  double Ad, Wd;			//fluctuations
  double psi, teta, phi;	 	//angles d'Euler
  int    nt = 50000;			//nombre de mesures par proc
  double t  = 0.0;     			//temps initiale sur le proc
  double dt = 0.00001; 			//intervalle de temps entre 2 mesures
  int chunk;				
  int nb=1;    			        //nombre de fichiers à écrire 
  double x0[3], u, v, w, x[3];		//coordonnées
  double R0_v[3][3], Rv_w[3][3], Rw_[3][3];   //matrices de rotation de R0 dans Rv, Rv dans Rw et Rw dans R.
  double N;            			//norme
  double max_sweep;			//arc balayé
  double T;				//période
  int nside = 1;                                                             
  long npix, np, ns, ipring;
  int* pix;
  double *map;
  char coordsys[1];      
  coordsys[0]='C';      
  int flag=1;
  x0[0]= 0.5;		//pointage instrument   
  x0[1]= 0.86;		//
  x0[2]= 0.0;		//

  Ad = 0.1;		//amplitude fluctuation 
  Wd = 0.0023;		//frequence fluctuation

  Wz0 = 0.0;		//vitesse tanguage 
  Wu = 60.0;		// orbite L2 + spin sat (60/h + 1/h)
  Wz = 0.00014155;	//vitesse rotation autour du soleil 1/an


  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /*option*/
  while ((ch = getopt( argc, argv, "f:p:q:r:d:n:m:h:o:")) != EOF) {
    switch(ch) {
      case 'p':
         Wz0 = strtod(optarg, &endp);
        if ((errno == ERANGE) || (errno != 0 && Wz0 <= 0)){
          usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      break;
      case 'q':
         Wu = strtod(optarg, &endp);
        if ((errno == ERANGE) || (errno != 0 && Wu <= 0)){
          usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      break;
      case 'r':
         Wz = strtod(optarg, &endp);
        if ((errno == ERANGE) || (errno != 0 && Wz <= 0)) {
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
      case 'd':
        dt = strtod(optarg, &endp);
        if ((errno == ERANGE) || (errno != 0 && dt <= 0)) {
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
        usage();
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 0;
      break;
    }
  }
  npix = nside2npix(nside);

  if(rank==0){
    printf("\n========sky scan====================\n", nside, npix);
    printf("x0        (x,y,z)  : \t ( %lf %lf %lf )\n", x0[0], x0[1], x0[2]); 
    printf("rotations (p,q,r)  : \t %lf, %lf and %lf rph\n", Wz0, Wu, Wz); 
    printf("nb parts           : \t %d ", size);
    printf("sampling pulse(dt) : \t %lf hours (%lf sec)\n", dt, dt*3600 );
    printf("nb samples         : \t %d per part ( global %d )\n", nt, size*nt);
    printf("observation period : \t %lf hour per chunk (whole %lf days)\n", nt*dt, size*dt*nt/24);
    printf("nside              : \t %d \n", nside);
    printf("npix               : \t %d \n", npix);
    printf("filename           : \t %s \n", bn);
  }
  
  t=rank * nt * dt;                                      //set init time on each proc 
  pix= (int*) malloc(nt*sizeof(int));  
  for(i=0; i< nt; i++){

    psi = 2*M_PI*Wz0 + Ad*sin(2*M_PI*Wd*t);                                    //compute angular rotation
    teta = 2*M_PI*Wu*t;
    phi = 2*M_PI*Wz*t ;+ Ad*cos(2*M_PI*(Wd+0.01)*t+0.02);

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
    t+=dt;
   //printf("\t%ld", ipring);
    //fprintf(stream, "%lf %lf %lf \n ", x[0], x[1], x[2]);
  } 
  
  MatAllocate(&At, nt, 1, MPI_COMM_WORLD);             //init distribuated matrix
  

  x1    = (double *) malloc(nt* sizeof(double));          //< allocate vector values  according partitioning ant init to 0.0
  
  for(i=0; i<nt; i++){   
      x1[i]=(1.0);
  }

  
  MatSetIndices(&At, nt, 0, 1, pix);
  MatSetValues(&At, nt, 0, 1, x1);

  free(pix);

  MPI_Barrier(MPI_COMM_WORLD);       
  t0 = MPI_Wtime();
  MatLocalShape(&At,0,0);
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  t1=t1-t0;
  MPI_Reduce(&t1, &max_t1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&t1, &min_t1, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&t1, &sum_t1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(rank==0)
    printf("\nLocalShape :\t%lf\t%lf\t%lf", max_t1, sum_t1/size, min_t1);
  

  y1_i  = (int *) calloc(At.lcount, sizeof(int));                //< allocate full vector indices on each proc and init to 0 
  y1_v  = (double *) calloc(At.lcount, sizeof(double));          //< allocate full vector values  on each proc and init to 0.0 

/*  MPI_Barrier(MPI_COMM_WORLD);                              // Pre  

  t0 = MPI_Wtime();
  MatVecProd(&At, x1, y1_v);
  t1 = MPI_Wtime();
  if(rank==0)
    printf("\n$ MapMatVecProd execution time :\t%lf sec \n", t1-t0);
  */

  MPI_Barrier(MPI_COMM_WORLD);                              // Pre  
  t0 = MPI_Wtime();
  MatComShape(&At, flag);
  MPI_Barrier(MPI_COMM_WORLD);
  t1 = MPI_Wtime();
  t1=t1-t0;
  MPI_Reduce(&t1, &max_t1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&t1, &min_t1, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&t1, &sum_t1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(rank==0)
    printf("\nComShape :\t%lf\t%lf\t%lf", max_t1, sum_t1/size, min_t1);

  MatInfo(&At,0, "scan_mat");   
  

  MPI_Barrier(MPI_COMM_WORLD);                               

  t0 = MPI_Wtime();
  TrMatVecProd(&At, x1, y1_v, 0);
  t1 = MPI_Wtime();
  t1=t1-t0;
  MPI_Reduce(&t1, &max_t1, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&t1, &min_t1, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&t1, &sum_t1, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(rank==0) 
    printf("\nProd :\t%lf\t%lf\t%lf", max_t1, sum_t1/size, min_t1);

  map    = (double *) calloc(npix, sizeof(double));   
  for(i=0; i<At.lcount; i++){
    map[At.lindices[i]]= y1_v[i];
  }
  //write_healpix_map(map, nside, bn, '1', coordsys);
  sprintf(fn, "%s_%d.gif", bn, rank);
     
  int filenamelen = strlen(fn);
  //__mapoutput_MOD_map2gif((int*)(&nside), map, &filenamelen, fn);
  __mapoutput_MOD_map2gif((int*)(&nside), map, &filenamelen, fn);
  printf("\n:: --- Returned from Fortran routine\n\n");

  MatFree(&At); 
//  free(y1_i);
  free(y1_v);
  free(x1);
//  free(map); 
  MPI_Finalize();       
  return 0;
}


