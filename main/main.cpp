#include "tmLQCD.h"

#include "distillery.h"
#include "input_parms.h"

#include "mpi.h"

#include "unistd.h"

#include <cstdlib>

int main(int argc, char *argv[]){

  //MPI initialisation stuff
  int mpi_thread_provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &mpi_thread_provided);
  int numprocs = 0, myid = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Group world_group;
  MPI_Comm_group( MPI_COMM_WORLD, &world_group );
  MPI_Comm mpi_comm_world_2;
  MPI_Comm_create( MPI_COMM_WORLD, world_group, &mpi_comm_world_2 );

  Eigen::initParallel();
  Eigen::setNbThreads(1);

  // initialisation of the twisted mass stuff - MUST BE the first thing to do
  int verbose = 1; // set to 1 to make tmLQCD more verbose
  tmLQCD_invert_init(argc, argv, verbose, myid);
  MPI_Barrier(MPI_COMM_WORLD);

  // initialisation of distillery
  LapH::input_parameter param;
  param.parse_input_file(argc, argv);
  if( param.nb_rnd > 1 ){
    if( myid == 0 ){
      std::cout 
        << "This version of peram_gen is a complete and utter hack and can only " << std::endl
        << "do one random vector at a time. This is the result of a number of " << std::endl
        << "memory conserving modifications. If you're not on dense GPU nodes " << std::endl
        << "with far too little memory, you're better off using the master or zgemm " << std::endl
        << "branches. peram_gen will terminate now" << std::endl;
      fflush(stdout);
    }
    tmLQCD_finalise();
    MPI_Finalize();
    exit(123);
  }

  
  if(myid == 0) {
    std::cout << "processing config: " << param.config << "\n" << std::endl;
  }
  MPI_Barrier(MPI_COMM_WORLD);

  int gauge_read = tmLQCD_read_gauge(param.config);
  MPI_Barrier(MPI_COMM_WORLD);
  
  if( gauge_read < 0 ){
    tmLQCD_finalise();
    printf("There was an error in tmLQCD_read_gauge!\n");
    MPI_Finalise();
    exit(222);
  }

  LapH::distillery dis;
  dis.initialise(param);
  MPI_Barrier(MPI_COMM_WORLD);

  // preparing source creation -------------------------------------------------
  size_t nb_of_inversions =  param.dilution_size_so[2];
    std::complex<double>** sources = new std::complex<double>*[nb_of_inversions];
    std::complex<double>** propagators_t0 = new std::complex<double>*[nb_of_inversions];
    std::complex<double>** propagators_t1 = new std::complex<double>*[nb_of_inversions];

  int length = 3*4*param.Lt*param.Ls*param.Ls*param.Ls/numprocs;

  // loop over random vectors
  for(size_t rnd_id = 0; rnd_id < param.nb_rnd; ++rnd_id) {
    // take some memory allocation overhead to reduce overall memory usage
    // these will be deleted before the perambulator is written (which
    // requires significant temporary storage)
    for(size_t i = 0; i < nb_of_inversions; ++i){
      sources[i] = new std::complex<double>[length];
      propagators_t1[i] = new std::complex<double>[length];
      propagators_t0[i] = new std::complex<double>[length];
    }

    std::complex<double>** propagators_temp;
    omp_set_num_threads(param.peram_gen_omp_num_threads);
    #pragma omp parallel
    {
      int num_threads = omp_get_num_threads();
      int thread_id   = omp_get_thread_num();
      
      // loop over all inversions
      for(size_t dil_t = 0; dil_t < param.dilution_size_so[0]; ++dil_t){
        for(size_t dil_e = 0; dil_e < param.dilution_size_so[1]; ++dil_e){
          // the first thread generates the sources and does the inversions
          // the second thread moves on to the barrier and waits      
          if( num_threads == 1 || thread_id == 0 ){

            dis.create_source(dil_t,dil_e,sources);
            
            for(size_t dil_d = 0; dil_d < param.dilution_size_so[2]; ++dil_d){        

              if(myid == 0) 
                std::cout << "\t\nDoing inversions at: t = " << dil_t << "\t e = " << dil_e 
                          << "\t d = " << dil_d << "\n" << std::endl;
  
              // tmLQCD can also write the propagator, if requested
              unsigned int op_id = 0;
              unsigned int write_prop = 0;
#ifndef PG_QUDA_DIRECT
              {  // open a block for convenience 
#else
              if(param.quda_direct){
                // use direct passthrough to QUDA, saving a lot of time when memory usage is high
                // the last parameter indicates that we would like the gauge field to be kept resident
                // in device memory between inversions, this breaks the processing of multiple configurations
                // FIXME/TODO: support multiple configurations
                invert_quda_direct((double*) propagators_t0[dil_d], (double*) sources[dil_d], op_id);
              } else {
#endif
                // use standard inverter interface
                tmLQCD_invert((double *) propagators_t0[dil_d], (double *) sources[dil_d], op_id, write_prop);
              } // close either the if statement or the block (see #ifndef PG_QUDA_DIRECT)
            }
          }
          
          // this barrier is necessary to prevent concurrent access to propagators_t0 and propagators_t1
          #pragma omp barrier
          if( num_threads == 1 || thread_id == 0){
            // exchange pointers 
            propagators_temp = propagators_t1;
            propagators_t1 = propagators_t0;
            propagators_t0 = propagators_temp;
          }
          #pragma omp barrier
          
          // the second thread will enter add_to_perambulator while the first thread will move on to the next
          // iteration, generate the next set of sources and do the inversion
          if( num_threads == 1 || thread_id == 1 ){
            dis.add_to_perambulator(dil_t,dil_e,propagators_t1);
          }
        }
      } // end of loop over inversions
    } // OpenMP parallel section closing brace

    // ----------- this is a total hack to use less memory
    for(size_t i = 0; i < nb_of_inversions; ++i){
      delete[] sources[i];
      delete[] propagators_t1[i];
      delete[] propagators_t0[i];
    }
    printf("HACK finalize tmLQCD\n");
    fflush(stdout);
    tmLQCD_finalise();
    sleep(5);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("HACK clean distillery\n");
    fflush(stdout);
    dis.hack_clean();
    sleep(5);
    // ------------ this is a total hack to use less memory

    // constructing the perambulator
    MPI_Barrier(MPI_COMM_WORLD);
    // creating new random vector and writing perambulator to disk -------------------
    dis.write_perambulator_to_disk(rnd_id); // ---------------------------------------
    MPI_Barrier(MPI_COMM_WORLD);
    // memory reduction hack implies that this doesn't work anymore
    //if(rnd_id < param.nb_rnd - 1)
    //  dis.reset_perambulator_and_randomvector(rnd_id+1);
    //MPI_Barrier(MPI_COMM_WORLD);
  } // end of loop over random vectors

  dis.clean();

  MPI_Finalize();
	return 0;
}

