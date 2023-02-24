/** @file   test_pcg_mapmat.c
  
    @author Pierre Cargemel
    @date   December 2011 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "midapack.h"


extern char *optarg;

void usage(){ 
    printf("\n####################################################################");
    printf("\n#                                                                  #");
    printf("\n#  cg_mapmat est un exemple d'utilisation de la librairie midapack #");
    printf("\n#  dont il illustre les fonctionnalités (produit matrice-vecteur,  #");
    printf("\n#  transposée-vecteur). Le programme réalise une carte aux moindres#");
    printf("\n#  carrés. Pour cela il résoud, par méthode du gradient conjugué,  #");
    printf("\n#  le système :                                                    #");
    printf("\n#                                                                  #");
    printf("\n#                      A^t A x = A^t b                             #");
    printf("\n#                                                                  #");
    printf("\n#  avec A une matrice creuse rectangulaire, dont les Aij sont non- #");
    printf("\n#  nuls si l'éléments j de la carte a été observé par la mesure i. #");
    printf("\n#  b est un vecteur contenant les valeurs des différentes mesures. #");
    printf("\n#                                                                  #");
    printf("\n#  Usage : cg_mapmat [option][valeur]                              #");
    printf("\n#  -M [int] spécifie le nombre globales de lignes de A(nb mesures) #");
    printf("\n#  -Z [int] nombre éléments non-nuls par ligne                     #");
    printf("\n#  -f [int] option algorithme de communication [0:4]               #");
    printf("\n#  -s [int] option algorithme de trie [0 5]                        #");
    printf("\n#  -p [int] option mutli-thread [0,1]                              #");
    printf("\n#  -k [int] nombre maximum d'itération du gradient conjugué        #");
    printf("\n#  -e [double] tolerence residu final du gradient conjugué         #");
    printf("\n#  -o [filename] nom de fichier sortie                             #");
    printf("\n#  -i option affichant les information sur la matrice A            #"); 
    printf("\n#                                                                  #");
    printf("\n####################################################################\n");
}

int main(int argc, char *argv[]){
  int		M, N, Nnz;		//nbre globales de lignes, de colonnes, de valeurs non-nulles par colonne
  int		m, n;			//nbre locales de lignes, de colonnes 
  int 		gif;			//indice globale de la première ligne locale
  int		err, i, j, k;
  int           K;	                //nbre maximum d'iteration du CG
  double 	tol;			//tolerence du résidu final du CG
  Mat	A;			        //structure de matrice midapack
  int 		*indices;
  double 	*values;
  int 		flag, sort, omp ;	//option pour le schéma de communication, le trie, le parallélisme
  double	*b, *z, *Ag, *Ad, *NAd; 	 	//vecteurs du domaine temporelle 
  double	*x, *g, *d, *AtNAd;	//vecteurs du domaine de la carte
  double        alpha, gamma, resold, resnew;
  double 	localreduce;
  double	st, t;		 	//timer, start time
  char 		ch;
  char		*endp;
  char		*filename;		//nom de fichier de sortie
  int 		output, timer, info;          
  int 		rank, size;	        //numéro et nombre de processus
  
  MPI_Init(&argc, &argv);               //démarage MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); //
  MPI_Comm_size(MPI_COMM_WORLD, &size); //
 
  flag=1;				//valeur par défaut
  sort=1;                               //
  omp=0;                                //
  output=0;				//
  timer=0;				//
  info=0; 				//
  err=0;				//
  Nnz=1;				//
  K=100;				//
  tol=0.001;				//

  if(argc < 2){
      usage();
  }

  while((ch = getopt( argc, argv, "M:Z:f:s:k:e:p:o:i" )) != EOF){
    switch(ch) {
      case 'M':
	M = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && (M <= 0 || M >= 2147483648))){
          err=1;
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'Z':
        Nnz = strtol(optarg, 0, 10);
        if ((errno == ERANGE) || (errno != 0 && (Nnz <= 0 || Nnz >= 2147483648))){
	  err=1;
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'f':
 	flag = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && (flag < 0 || flag > 5))){
	  err=1;
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 's':
 	sort = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && (sort < 0 || sort > 5))){
	  err=1;
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'k':
 	K = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && K < 0)){
	  err=1;
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'e':
        tol = strtod(optarg, &endp);
	if ((errno == ERANGE) || (errno != 0 && tol < 0)){
	  err=1;
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      case 'p':
 	omp = strtoul(optarg, 0, 10);
	if ((errno == ERANGE) || (errno != 0 && (omp < 0 || omp > 3))){
	  err=1;
          usage();
          MPI_Abort(MPI_COMM_WORLD, err);
        }
      break;
      case 'o':
        filename = strdup(optarg);
        output=1;
      break;
      case 'i':
        info=1;
      break;
    }
  }

  N=M;

  partition(&gif, &m, M, rank, size);                     

  MatAllocate(&A, m, Nnz, MPI_COMM_WORLD);	//crée une matrice avec m lignes locales et Nnz éléments non-nuls par ligne

  indices  = (int *) malloc(m * sizeof(int));	
  values  = (double *) malloc(m * sizeof(double));
  for(j=0; j<Nnz; j++){
    for(i=0; i<m; i++){
      if(j==0){                    
        values[i] =  ((double)(gif+i+1.0))/M;
        indices[i] = (gif+i)%N;
      }
      else{
        values[i] =  0.0;
        indices[i] = (gif+i+j)%N;
      }  
    }
    MatSetIndices(&A, m, j, Nnz, indices);
    MatSetValues(&A, m, j, Nnz, values);
  }
  free(indices);
  free(values);
  
  if(output==1)
    MatSave(&A, filename);  			//écrit la matrice dans des fichiers

  MatLocalShape(&A, sort, omp);			//organise la structure de données de la matrice

  x   = (double *) malloc(A.lcount*sizeof(double));//crée une carte initiale 
  for(j=0; j<A.lcount; j++)
    x[j] = 0.0;

  srand(gif);					//initialise le générateur
  b   = (double *) malloc(m*sizeof(double));    //crée un vecteur de mesure
  for(i=0; i<m; i++)
    b[i] = 1.0 ;//+ (2*((double) rand()) / RAND_MAX -1);
  
  err = MatComShape(&A, flag);         		//crée un schéma de communication des données de la matrice

  if(info==1)
    MatInfo(&A, 0, "lhs");

  /** Algorithme du Gradient Conjugué légèrement modifier par rapport à l'algorithme classique.
  On effectue succesivement les produits A^t puis A. Pour des raisons de simplicité,
  les produits scalaires ne s'effectue pas dans le domaine carte mais dans le domaine temporelle.
  En particulier on utilise : < A^t y, x > = < y , Ax > .

  Utilisation mémoire : pour l'algorithme on a besoin de 5 vecteurs supplémentaire :
  - 2 dans le domaine de la carte(gradient, direction),
  - 3 dans le domaine temporelle. 
  
  Complexité : à chaque itération de l'algorithme du gradient conjugué on a :
  - 4 produits scalaires dans le domaine temporelle, 
  - 3 axpy dans le domaine de la carte,
  - 3 multiplication par A,
  - 1 multiplication par A^t.
  **/
  st=MPI_Wtime();			//initialisation
  g  = (double *) malloc(A.lcount*sizeof(double));	//vecteurs domaine carte 
  d  = (double *) malloc(A.lcount*sizeof(double));	// 
  AtNAd  = (double *) malloc(A.lcount*sizeof(double));	// 

  z = (double *) malloc(m*sizeof(double));    	//vecteurs domaine temps 
  Ad = (double *) malloc(m*sizeof(double)); 	//
  NAd = (double *) malloc(m*sizeof(double));	//
  Ag = (double *) malloc(m*sizeof(double));	//
  
  MatVecProd(&A, x, z, 0);		//z = N (Ax - b)
  for(i=0; i<m; i++)			//...
    z[i] = z[i]-b[i];			//...
  //mpi stbmm(toeplitz_N, z, z)		//...

  TrMatVecProd(&A, z, g, 0);		//g = At N (Ax - b) = Atz
  for(j=0; j<A.lcount; j++)             //d = g 
    d[j]=g[j];

  MatVecProd(&A, g, Ag, 0);		//Ag = A g = A Atz
  for(i=0; i<m; i++)			//Ad = Ag         
    Ad[i]=Ag[i];

  resnew=0.0;				//resnew = ||g|| = <Ag, z>
  localreduce=0.0;			//
  for(i=0; i<m; i++)			//         
    localreduce+=Ag[i]*z[i];  		//
  MPI_Allreduce(&localreduce, &resnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

  if(rank==0){
    printf("\ngradient conjugué : résidu initial %lf max iterations %d tolerence %lf", resnew, K, tol);
    printf("\niter\tresidu   \ttemps");
  }

  for(k=0; k<K ; k++){                  //début boucle

    //mpi stbmm(toeplitz_N, Ad, NAd)	//NAd= N Ad
  
    gamma =0.0;
    localreduce=0.0;			//gamma = <d, At N Ad> = <Ad, N Ad>
    for(i=0; i<m; i++)                  //         
      localreduce+=Ad[i]*NAd[i];	//
    MPI_Allreduce(&localreduce, &gamma, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);//printf("\ngamma %lf ", gamma);


    alpha=resnew/gamma;                 //x = x - alpha d
    for(j=0; j<A.lcount; j++)           // 
      x[j] = x[j] - alpha*d[j];		//
      
    for(i=0; i<m; i++)			//z = z - alpha N Ad
      z[i] = z[i] - alpha* NAd[i];	//

    TrMatVecProd(&A, NAd, AtNAd, 0);	//AtNAd = At N Ad

    for(j=0; j<A.lcount; j++)           //g = g - alpha At N Ad
      g[j] = x[j] - alpha*AtNAd[j];	//
  
    MatVecProd(&A, g, Ag, 0);		//Ag = A g = A Atz

    resold=resnew;              	//resnew = <g, g> = <Ag, z>
    resnew=0.0;				//
    localreduce=0.0;			//
    for(i=0; i<m; i++)			//         
      localreduce+=Ag[i]*z[i];		//
    MPI_Allreduce(&localreduce, &resnew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD); 

    t=MPI_Wtime();
    if(rank==0)
      printf("\n %d\t %lf\t %lf", k, resnew, t-st); 
    if(resnew<tol){
      if(rank ==0)
        printf("\n--> convergé (%lf < %lf)\n", resnew, tol);
      break;
    }
    st=MPI_Wtime();
 
    for(j=0; j<A.lcount; j++)     	//d = -g + (resnew/resold) d
      d[j]= -g[j] + d[j]*resnew/resold; // 

    MatVecProd(&A, d, Ad, 0);		//Ad = A d

  }
  free(z);
  free(g);
  free(Ag);
  free(d);
  free(Ad);
  free(NAd);
  free(AtNAd);
  /**
  Fin algorithme gradient conjugué
  **/ 
  printf("\n");

  MatFree(&A);                                                //free memory  
  free(b);
  free(x);
  MPI_Finalize();
  
  return 0;
}


