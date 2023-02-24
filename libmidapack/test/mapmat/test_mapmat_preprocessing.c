/** @file   test_mapmat_preprocessing.c
    @brief  Building of a PT-SCOTCH distribuated graph associate with a MAPMAT 

    Just Creates a matrix, builds and destroys a distribuated graph.
    Atfer building, the distributed graph is saved into files.
    You can specify the filename and any parameters describing the graph.
    For more information, see usage.  
    @author Pierre Cargemel
    @date   December 2011 */
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "mapmat.h"


extern char *optarg;

/** @brief display common available options

    Usage: test_mapmat_preprocessing -f [filename] -v [val] -w [val] -c [val] -m [val] 
    file should be set to the .dgr format. Each processor open a file named filename_%p_%r.dgr (where %p is number of processes and %r is rank of processor).
    So the given filename is just the root of the full name.
    Exemple : to save file graph_2_0.dgr and graph_2_1.dgr, just execute mpirun -n 2 test_mapmat_preprocessing -f graph*/
void usage(){
  fprintf(stderr, "\n*#  Usage: test_mapmat_preprocessing -f [filename]");
  fprintf(stderr, "\n*#  -M [number of rows]");
  fprintf(stderr, "\n*#  -N [number of columns]");
  fprintf(stderr, "\n*#  -Z [number non-zero values per columns]");
  fprintf(stderr, "\n*#  -B [number of blocks]");
  fprintf(stderr, "\n*#  Exemple : to load file graph_2_0.dgr and graph_2_1.dgr, just execute mpirun -n 2 test_mapat_prprocessing -f graph");
}

/** @brief check and return error code

    */
int check(int err){
    if(err==1){
      printf("..............test failed\n");
      return 1;
    }
    printf("......successfully tested\n");
    return 0;
}


int main(int argc, char *argv[]){

  int i, k, err, provided, rank, size;
  int connect, com, vertices, weigh;
  FILE *stream;
  char *bn;
  char fn[100];
  char ch;

  UINT rows, cols;
  int Nnz, Nbls;
  double t0, t1;
  MAPMAT A;

  SCOTCH_Dgraph dgraph;                


#ifdef SCOTCH_PTHREAD
  int                 thrdlvlreqval;
  int                 thrdlvlproval;
#endif /* SCOTCH_PTHREAD */

#ifdef SCOTCH_PTHREAD
  thrdlvlreqval = MPI_THREAD_MULTIPLE;
  if (MPI_Init_thread (&argc, &argv, thrdlvlreqval, &thrdlvlproval) != MPI_SUCCESS)
    printf("main: Cannot initialize (1)");
  if (thrdlvlreqval > thrdlvlproval)
    printf("main: MPI implementation is not thread-safe: recompile without SCOTCH_PTHREAD");
#else /* SCOTCH_PTHREAD */
  if (MPI_Init (&argc, &argv) != MPI_SUCCESS)
    printf("main: Cannot initialize (2)");
#endif /* SCOTCH_PTHREAD */


  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  /*option*/
  while ((ch = getopt( argc, argv, "M:N:Z:B:f:" )) != EOF) {
    switch(ch) {
      case 'N':
        cols = strtoul(optarg, 0, 10);
        if ((errno == ERANGE) || (errno != 0 && cols <= 0)) {
          usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      break;
      case 'M':
        rows = strtoul(optarg, 0, 10);
        if ((errno == ERANGE) || (errno != 0 && rows <= 0)) {
          usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      break;
      case 'Z':
        Nnz = strtol(optarg, 0, 10);
        if ((errno == ERANGE) || (errno != 0 && Nnz <= 0)) {
        usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      break;
      case 'B':
       Nbls = strtol(optarg, 0, 10);
        if ((errno == ERANGE) || (errno != 0 && Nbls <= 0)) {
        usage();
          MPI_Abort(MPI_COMM_WORLD, 1);
        }
      break;
      case 'f':
        bn = strdup(optarg);
      break;
    }
  }

  sprintf(fn, "%s_%d_%d.dgr", bn, size, rank );
  printf ("(%d) open %s \n",rank, fn);
  
  err = MapMatGen(&A, MPI_COMM_WORLD, rows, cols, Nnz, Nbls);
  
  err = SCOTCH_dgraphInit(&dgraph, MPI_COMM_WORLD);

  err = MapMatBuildGraph(&A , &dgraph);
  
  /*check*/
  err    = SCOTCH_dgraphCheck(&dgraph);
  if(err!=0){
    printf("error : invalid graph \n");
    MPI_Finalize();
    return 1;
  }
  //open file
  stream=fopen(fn,"w");
 
  err = SCOTCH_dgraphSave(&dgraph, stream);
  if(err!=0)
      printf("cannot save dgraph\n");
  
  fclose(stream);

  SCOTCH_dgraphExit(&dgraph);
  
  MPI_Finalize();

  return check(err);
}
