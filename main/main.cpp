#include <omp.h>

#include "tmLQCD.h"

#include "distillery.h"
#include "input_parms.h"

#include "mpi.h"

int main(int argc, char *argv[]){

  //MPI initialisation stuff
  //MPI_Init(&argc, &argv);
  int mpi_thread_provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &mpi_thread_provided);
  int numprocs = 0, myid = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  
  Eigen::initParallel();
  Eigen::setNbThreads(12);

  // initialisation of the twisted mass stuff - MUST BE the first thing to do
  int verbose = 1; // set to 1 to make tmLQCD more verbose
  tmLQCD_invert_init(argc, argv, verbose, myid);
  MPI_Barrier(MPI_COMM_WORLD);

  // initialisation of distillery
  LapH::input_parameter param;
  param.parse_input_file(argc, argv);
  if(myid == 0)
    std::cout << "processing config: " << param.config << "\n" << std::endl;
  MPI_Barrier(MPI_COMM_WORLD);

  tmLQCD_read_gauge(param.config);
  MPI_Barrier(MPI_COMM_WORLD);

  LapH::distillery dis;
  dis.initialise(param);
  MPI_Barrier(MPI_COMM_WORLD);

  // preparing source creation -------------------------------------------------
  size_t nb_of_inversions =  param.dilution_size_so[2];

  int length = 3*4*param.Lt*param.Ls*param.Ls*param.Ls/numprocs;
  std::complex<double>** sources = new std::complex<double>*[nb_of_inversions];
  for(size_t i = 0; i < nb_of_inversions; ++i)
        sources[i] = new std::complex<double>[length];

  std::complex<double>* propagator = new std::complex<double>[length];

  // loop over random vectors
  for(size_t rnd_id = 0; rnd_id < param.nb_rnd; ++rnd_id) {
    // loop over all inversions
    for(size_t dil_t = 0; dil_t < param.dilution_size_so[0]; ++dil_t){
      for(size_t dil_e = 0; dil_e < param.dilution_size_so[1]; ++dil_e){        
        
        dis.create_source(dil_t,dil_e,sources);
        
        for(size_t dil_d = 0; dil_d < param.dilution_size_so[2]; ++dil_d){        

          if(myid == 0) 
            std::cout << "\t\nDoing inversions at: t = " << dil_t << "\t e = " << dil_e 
                      << "\t d = " << dil_d << "\n" << std::endl;
  
          // tmLQCD can also write the propagator, if requested
          unsigned int op_id = 0;
          // the 1 at the end indicates that we would like the gauge field to persist in device memory
          // as a result, need to have one job per gauge configuration
          invert_quda_direct((double*) propagator, (double*) sources[dil_d], op_id, 1);
          MPI_Barrier(MPI_COMM_WORLD);
          memcpy((void*)sources[dil_d], (void*)propagator, 2*sizeof(double)*length);
          MPI_Barrier(MPI_COMM_WORLD);

        }

        dis.add_to_perambulator(dil_t,dil_e,sources);

      }
    } // end of loop over inversions

    // constructing the perambulator
    MPI_Barrier(MPI_COMM_WORLD);
    // creating new random vector and writing perambulator to disk -------------------
    dis.write_perambulator_to_disk(rnd_id); // ---------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    if(rnd_id < param.nb_rnd - 1)
      dis.reset_perambulator_and_randomvector(rnd_id+1);
    MPI_Barrier(MPI_COMM_WORLD);
  } // end of loop over random vectors

  tmLQCD_finalise();
  dis.clean();

  MPI_Finalize();
	return 0;
}

