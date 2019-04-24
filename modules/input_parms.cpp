#include "input_parms.h"

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// deep copy constructor
LapH::input_parameter::input_parameter(const input_parameter& other){

  for (size_t i = 0; i < 3; i++){
    dilution_size_so[i] = other.dilution_size_so[i];
    dilution_type_so[i] = other.dilution_type_so[i];
  }
  nb_of_sinks = other.nb_of_sinks;
  nb_of_sink_rnd_vec = other.nb_of_sink_rnd_vec;
  seed_si.resize(0);
  std::copy(other.seed_si.begin(), other.seed_si.end(), 
            std::back_inserter(seed_si));
  for(size_t nbs = 0; nbs < dilution_size_si.size(); nbs++){
    dilution_size_si[nbs].resize(0);
    dilution_type_si[nbs].resize(0);
  }
  dilution_size_si.resize(0);
  dilution_type_si.resize(0);
  std::copy(other.dilution_size_si.begin(), other.dilution_size_si.end(), 
            std::back_inserter(dilution_size_si));
  std::copy(other.dilution_type_si.begin(), other.dilution_type_si.end(), 
            std::back_inserter(dilution_type_si));

  config = other.config;
  nb_config = other.nb_config;
  delta_config = other.delta_config;
  Ls = other.Ls;
  Lt = other.Lt;
  nb_ev = other.nb_ev;
  verbose = other.verbose;
  endianness = other.endianness;

  nb_rnd = other.nb_rnd;
  if (rnd_id != NULL)
    delete[] rnd_id;
  rnd_id = new int[nb_rnd];
  if (seed != NULL)
    delete[] seed;
  seed = new int[nb_rnd];
  for (size_t i = 0; i < nb_rnd; i++){
    rnd_id[i] = other.rnd_id[i];
    seed[i] = other.seed[i];
  } 
  
  quarktype = other.quarktype;
  outpath = other.outpath;
  inpath_ev = other.inpath_ev;
  peram_file_name = other.peram_file_name;
  rnd_vec_file_name = other.rnd_vec_file_name;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// deep assigment operator
LapH::input_parameter& LapH::input_parameter::operator=
                                                 (const input_parameter& other){
  if(this != &other){
    for (size_t i = 0; i < 3; i++){
      dilution_size_so[i] = other.dilution_size_so[i];
      dilution_type_so[i] = other.dilution_type_so[i];
    }
    nb_of_sinks = other.nb_of_sinks;
    nb_of_sink_rnd_vec = other.nb_of_sink_rnd_vec;
    seed_si.resize(0);
    std::copy(other.seed_si.begin(), other.seed_si.end(), 
              std::back_inserter(seed_si));
    for(size_t nbs = 0; nbs < dilution_size_si.size(); nbs++){
      dilution_size_si[nbs].resize(0);
      dilution_type_si[nbs].resize(0);
    }
    dilution_size_si.resize(0);
    dilution_type_si.resize(0);
    std::copy(other.dilution_size_si.begin(), other.dilution_size_si.end(), 
              std::back_inserter(dilution_size_si));
    std::copy(other.dilution_type_si.begin(), other.dilution_type_si.end(), 
              std::back_inserter(dilution_type_si));
    
    config = other.config;
    nb_config = other.nb_config;
    delta_config = other.delta_config;
    Ls = other.Ls;
    Lt = other.Lt;
    nb_ev = other.nb_ev;
    verbose = other.verbose;
    endianness = other.endianness;

    peram_gen_omp_num_threads = other.peram_gen_omp_num_threads;
    eigen_omp_num_threads = other.eigen_omp_num_threads;
    evec_read_omp_num_threads = other.evec_read_omp_num_threads;

    use_zgemm = other.use_zgemm;

#ifdef PG_QUDA_DIRECT
    quda_direct = other.quda_direct;
#endif

    create_rnd_vec_subdirs = other.create_rnd_vec_subdirs;
    
    nb_rnd = other.nb_rnd;
    if (rnd_id != NULL)
      delete[] rnd_id;
    rnd_id = new int[nb_rnd];
    if (seed != NULL)
      delete[] seed;
    seed = new int[nb_rnd];
    for (size_t i = 0; i < nb_rnd; i++){
      rnd_id[i] = other.rnd_id[i];
      seed[i] = other.seed[i];
    } 

    quarktype = other.quarktype;
    outpath = other.outpath;
    inpath_ev = other.inpath_ev;
    peram_file_name = other.peram_file_name;
    rnd_vec_file_name = other.rnd_vec_file_name;
  }
  return *this;
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// checking that the dilution parameters are correctly set in infile
// input:  type       -> gives the type of the dilution scheme
//         max_size   -> maximal size of space
// output: nb_dil_vec -> sets (in case of full and no dilution) and checks the 
//                       number of dilution vectors
static void check_dilution_input(const std::string type, const size_t max_size, 
                                 size_t& nb_dil_vec){
  // check for type
  if( (type.compare("F") != 0) && 
      (type.compare("I") != 0) && 
      (type.compare("B") != 0) && 
      (type.compare("N") != 0) ) {
        std::cerr << "Inversion type has to be one of \"F\", \"I\"," \
                     " \"B\" or \"N\"." << std::endl;
        std::cerr << "Aborting..." << std::endl;
        MPI_Finalize();
        std::exit(0);
  }
  // check and set number of inversions in corresponding space
  if (type.compare("F") == 0)
     nb_dil_vec = max_size;
  else if (type.compare("N") == 0 )
     nb_dil_vec = 1;
  else {
    if(max_size % nb_dil_vec != 0) {
      std::cerr << "check number of inversions" << std::endl;
      MPI_Finalize();
      std::exit(0);
    }
  }

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::input_parameter::check_input_parameters(){

  if(peram_gen_omp_num_threads > 2 || peram_gen_omp_num_threads < 1){
    std::cout << "The number of OMP threads for this version of peram_gen can only be 1 or 2! The selected value was " <<
                 peram_gen_omp_num_threads << std::endl;
    exit(1);
  }

  if(config > 100000) {
    std::cout << "Please check whether configuration number is correct: " 
              << config << std::endl;
    MPI_Finalize();
    std::exit(0);
  }
  if(Ls > 500) {
    std::cout << "Please check whether Ls is correct: " <<  Ls << std::endl;
    MPI_Finalize();
    std::exit(0);
  }
  if(Lt > 500) {
    std::cout << "Please check whether Lt is correct: " <<  Lt << std::endl;
    MPI_Finalize();
    std::exit(0);
  }
  if(nb_ev > 5000) {
    std::cout << "Please check whether number of eigenvectors is correct: " 
              <<  nb_ev << std::endl;
    MPI_Finalize();
    std::exit(0);
  }
  if(nb_rnd > 1000) {
    std::cout << "Please check whether starting number of randomvector" \
                 " is correct: " <<  nb_rnd << std::endl;
    MPI_Finalize();
    std::exit(0);
  }
  if((quarktype.compare("u") != 0) && 
     (quarktype.compare("d") != 0) && 
     (quarktype.compare("s") != 0) && 
     (quarktype.compare("c") != 0) ) {
       std::cerr << "Quarktype has to be one of \"u\", \"d\", \"s\" or \"c\"." 
                 << std::endl;
       std::cerr << "Aborting..." << std::endl;
       MPI_Finalize();
       std::exit(0);
  }
  if(endianness != "little" && endianness != "big"){
    std::cerr << "\n\nEndianness in infile must be 'little' or 'big'!\n" 
              << "Aborting...\n" << std::endl;
       MPI_Finalize();
       std::exit(0);
  }
     
  check_dilution_input(dilution_type_so[0], Lt, dilution_size_so[0]); 
  check_dilution_input(dilution_type_so[1], nb_ev, dilution_size_so[1]); 
  check_dilution_input(dilution_type_so[2], 4, dilution_size_so[2]);
  for(size_t nbs = 0; nbs < nb_of_sinks; nbs++){
    check_dilution_input(dilution_type_si[nbs][0], Lt, 
                         dilution_size_si[nbs][0]); 
    check_dilution_input(dilution_type_si[nbs][1], Ls*Ls*Ls, 
                         dilution_size_si[nbs][1]); 
    check_dilution_input(dilution_type_si[nbs][2], 4, 
                         dilution_size_si[nbs][2]); 
    check_dilution_input(dilution_type_si[nbs][3], 3, 
                         dilution_size_si[nbs][3]); 
  } 
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::input_parameter::check_and_create_filenames(){

}

/**
 * @brief check if the number of elements read by fscanf corresponds to expectations
 *
 * @param nelem_actual number of elements actually returned by fscanf
 * @param nelem_expected number of elements that was expected
 * @param parname string constant name of the parameter
 */
void LapH::input_parameter::verify_fscanf(int const nelem_actual, int const nelem_expected, char const * const parname){
  int mpi_rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  if( nelem_actual != nelem_expected ){
    char msg[200];
    snprintf(msg, 200, "Error parsing %s!\n", parname);
    if( mpi_rank == 0 ) std::cout << msg << std::endl;
    exit(-10);
  }
}

// -----------------------------------------------------------------------------
// TODO: rework this parser completely and have some streamlined input
// and output routines
// -----------------------------------------------------------------------------
void LapH::input_parameter::parse_input_file(int argc, char *argv[]) {

  int opt = -1;
  int reader = 0;
  char infilename[200];
  char readin[256];
  FILE* infile = NULL;

  // search for command line option and put filename in "infilename"
  for(int i = 0; i < argc; ++i) {
    if(std::strcmp(argv[i], "-LapHsin") == 0) {
      opt = i+1;
      break;
    }
  }
  int myid = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  if(opt < 0) {
    if(myid==0){
      std::cout << "No input file specified, trying infile.in" << std::endl;
    }
    sprintf(infilename, "infile.in");
  } else {
    sprintf(infilename, "%s", argv[opt]);
    if(myid==0){
      std::cout << "Trying input file " << infilename << std::endl;
    }
  }
  // open file for reading
  if ((infile = fopen(infilename, "r")) == NULL ) {
    std::cerr << "Could not open file " << infilename << std::endl;
    std::cerr << "Aborting..." << std::endl;
    exit(-10);
  }
  // scan infile and check arguments ----------------------------------------
  // configs
  peram_gen_omp_num_threads = 1;
  reader = fscanf(infile, "omp_num_threads = %u\n", &peram_gen_omp_num_threads);
  verify_fscanf(reader, 1, "omp_num_threads");
  if(myid==0) std::cout << "omp_num_threads = " << peram_gen_omp_num_threads << std::endl;

  evec_read_omp_num_threads = 1;
  reader = fscanf(infile, "evec_read_omp_num_threads = %u\n", &evec_read_omp_num_threads);
  verify_fscanf(reader, 1, "evec_read_omp_num_threads");
  if(myid==0) std::cout << "evec_read_omp_num_threads = " << evec_read_omp_num_threads << std::endl;

#ifdef PG_QUDA_DIRECT
  {
    quda_direct = false;
    size_t temp;
    reader = fscanf(infile, "quda_direct = %zu\n", &temp);
    verify_fscanf(reader, 1, "quda_direct");
    if(myid==0) std::cout << "quda_direct = " << quda_direct << std::endl;
    quda_direct = (bool) temp;
  }
#endif

  reader = fscanf(infile, "config = %zu \n", &config);
  verify_fscanf(reader, 1, "config");
  if(myid==0) std::cout << "config = " << config << std::endl;
  
  reader = fscanf(infile, "total number of configs = %zu \n", &nb_config);
  verify_fscanf(reader, 1, "total number of configs");
  if(myid==0) std::cout << "total number of configs = " << nb_config << std::endl;
  
  reader = fscanf(infile, "distance between configs = %zu \n", &delta_config);
  verify_fscanf(reader, 1, "distance between configs");
  if(myid==0) std::cout << "distance betweeen configs = " << delta_config << std::endl;
  
  // spatial extent
  reader = fscanf(infile, "Ls = %zu \n", &Ls);
  verify_fscanf(reader, 1, "Ls");
  if(myid==0) std::cout << "Ls = " << Ls << std::endl;
  
  // temporal extent
  reader = fscanf(infile, "Lt = %zu \n", &Lt);
  verify_fscanf(reader, 1, "Lt");
  if(myid==0) std::cout << "Lt = " << Lt << std::endl;
  
  // number of eigenvectors
  reader = fscanf(infile, "nb_ev = %zu \n", &nb_ev);
  verify_fscanf(reader, 1, "nb_ev");
  if(myid==0) std::cout << "nb_ev = " << nb_ev << std::endl;
  
  // starting number of randomvectors
  reader = fscanf(infile, "nb_rnd = %zu \n", &nb_rnd);
  verify_fscanf(reader, 1, "nb_rnd");
  if(myid==0) std::cout << "nb_rnd = " << nb_rnd << std::endl;
  
  // seed and id for randomvectors
  rnd_id = new int[nb_rnd];
  seed = new int[nb_rnd];
  for (size_t i = 0; i < nb_rnd; i++){ 
    reader = fscanf(infile, "id %i seed %i\n", rnd_id+i, seed+i);
    char parstring[100];
    snprintf(parstring, 100, "id %lu seed", i);
    verify_fscanf(reader, 2, parstring); 
    if(myid==0) std::cout << "seed for rnd_vec " << i << " = " << seed[i] << std::endl;
  }
  
  reader = fscanf(infile, "verbose = %zu\n", &( verbose));
  verify_fscanf(reader, 1, "verbose");
  if(myid==0) std::cout << "verbose = " << verbose << std::endl;
  
  reader = fscanf(infile, "endianness = %255s\n", readin);
  verify_fscanf(reader, 1, "endianness");
  endianness.assign(readin);
  if(myid==0) std::cout << "endianness = " << readin << std::endl;


  {
    unsigned int temp;
    reader = fscanf(infile, "use_zgemm = %u\n", &(temp));
    verify_fscanf(reader, 1, "use_zgemm");
    use_zgemm = (bool)temp;
    if(myid==0) std::cout << "use_zgemm = " << use_zgemm << std::endl;
  }

  {
    unsigned int temp;
    reader = fscanf(infile, "hack_clean = %u\n", &temp);
    verify_fscanf(reader, 1, "hack_clean");
    hack_clean = (bool)temp;
    if(myid==0) std::cout << "hack_clean = " << hack_clean << std::endl;
  }

  {
    unsigned int temp;
    reader = fscanf(infile, "create_rnd_vec_subdirs = %u\n", &temp);
    verify_fscanf(reader, 1, "hack_clean");
    create_rnd_vec_subdirs = temp;
    if(myid==0) std::cout << "create_rnd_vec_subdirs = " << create_rnd_vec_subdirs << std::endl;
  }

  // quarktype
  reader = fscanf(infile, "quarktype = %255s\n", readin);
  verify_fscanf(reader, 1, "quarktype");
  if(myid==0) std::cout << "quarktype = " << readin << std::endl;
  quarktype.assign(readin);

  // SOURCE --------------------------------------------------------------------
  // type and number of inversions for the source in time ----------------------
  reader = fscanf(infile, "inversion_source_type_t = %255s\n", readin);
  verify_fscanf(reader, 1, "inversion_source_type_t");
  if(myid==0) std::cout << "inversion_source_type_t = " << readin << std::endl;
  dilution_type_so[0].assign(readin);
  
  reader = fscanf(infile, "inversion_source_number_t = %zu\n", &(dilution_size_so[0]));
  verify_fscanf(reader, 1, "inversion_source_number_t");
  if(myid==0) std::cout << "inversion_source_number_t = " << dilution_size_so[0] << std::endl;
  
  // type and number of inversions for the source in eigenvector space ---------
  reader = fscanf(infile, "inversion_source_type_v = %255s\n", readin);
  verify_fscanf(reader, 1, "inversion_source_type_v");
  if(myid==0) std::cout << "inversion_source_type_v = " << readin << std::endl;
  dilution_type_so[1].assign(readin);
  
  reader = fscanf(infile, "inversion_source_number_v = %zu\n", &(dilution_size_so[1]));
  verify_fscanf(reader, 1, "inversion_source_number_v");
  if(myid==0) std::cout << "inversion_source_number_v = " << dilution_size_so[1] << std::endl;
  // type and number of inversions for the soure in Dirac space ----------------
  
  reader = fscanf(infile, "inversion_source_type_d = %255s\n", readin);
  verify_fscanf(reader, 1, "inversion_source_type_d");
  if(myid==0) std::cout << "inversion_source_type_d = " << readin << std::endl;
  dilution_type_so[2].assign(readin);
  
  reader = fscanf(infile, "inversion_source_number_d = %zu\n", &(dilution_size_so[2]));
  verify_fscanf(reader, 1, "inversion_source_number_d");
  if(myid==0) std::cout << "inversion_source_number_d = " << dilution_size_so[2] << std::endl;

  // SINK ----------------------------------------------------------------------
  reader = fscanf(infile, "nb_of_sinks = %zu\n", &nb_of_sinks);
  verify_fscanf(reader, 1, "nb_of_sinks");
  if(myid==0) std::cout << "nb_of_sinks = " << nb_of_sinks << std::endl;

  dilution_size_si.resize(nb_of_sinks);
  dilution_type_si.resize(nb_of_sinks);
  seed_si.resize(nb_of_sinks);
  nb_of_sink_rnd_vec.resize(nb_of_sinks);
 
  // BaKo: surely this is still not working correctly...
  for(size_t nbs = 0; nbs < nb_of_sinks; nbs++){
    dilution_size_si[nbs].resize(4);
    dilution_type_si[nbs].resize(4);
    reader = fscanf(infile, "nb_of_sink_rnd_vec = %zu\n", &(nb_of_sink_rnd_vec[nbs]));
    verify_fscanf(reader, 1, "nb_of_sink_rnd_vec");
    seed_si[nbs].resize(nb_of_sink_rnd_vec[nbs]);
    if(myid==0) std::cout << "nb_of_sink_rnd_vec = " << nb_of_sink_rnd_vec[nbs] << std::endl;


    // seed for sink random vector
    for(size_t nb_rnd = 0; nb_rnd < nb_of_sink_rnd_vec[nbs]; nb_rnd++){
      reader = fscanf(infile, "seed = %i\n", &(seed_si[nbs][nb_rnd]));
      verify_fscanf(reader, 1, "seed");
      if(myid==0) std::cout << "seed = " << seed_si[nbs][nb_rnd] << std::endl;
    }

    // type and number of dilution vectors for the sink in time ----------------
    reader = fscanf(infile, "inversion_sink_type_t = %255s\n", readin);
    dilution_type_si[nbs][0].assign(readin);
    verify_fscanf(reader, 1, "inversion_sink_type_t");
    if(myid==0) std::cout << "inversion_sink_type_t = " << readin << std::endl;
    reader = fscanf(infile, "inversion_sink_number_t = %zu\n", &(dilution_size_si[nbs][0]));
    verify_fscanf(reader, 1, "inversion_sink_number_t");
    if(myid==0) std::cout << "inversion_sink_number_t = " << dilution_size_si[nbs][0] << std::endl;
    

    // type and number of dilution vectors for the sink in space ---------------
    reader = fscanf(infile, "inversion_sink_type_s = %255s\n", readin);
    verify_fscanf(reader, 1, "inversion_sink_type_s");
    dilution_type_si[nbs][1].assign(readin);
    if(myid==0) std::cout << "inversion_sink_type_s = " << readin << std::endl;
    reader = fscanf(infile, "inversion_sink_number_s = %zu\n", &(dilution_size_si[nbs][1]));
    verify_fscanf(reader, 1, "inversion_sink_number_s");
    if(myid==0) std::cout << "inversion_sink_number_s = " << dilution_size_si[nbs][1] << std::endl;
    
    // type and number of dilution vectors for the sink in Dirac space ---------
    reader = fscanf(infile, "inversion_sink_type_d = %255s\n", readin);
    verify_fscanf(reader, 1, "inversion_sink_type_d");
    dilution_type_si[nbs][2].assign(readin);
    if(myid==0) std::cout << "inversion_sink_type_d = " << readin << std::endl;
    reader = fscanf(infile, "inversion_sink_number_d = %zu\n", &(dilution_size_si[nbs][2]));
    verify_fscanf(reader, 1, "inversion_sink_number_d");
    if(myid==0) std::cout << "inversion_sink_number_d = " << dilution_size_si[nbs][2] << std::endl;
    
    // type and number of dilution vectors for the sink in colour space --------
    reader = fscanf(infile, "inversion_sink_type_c = %255s\n", readin);
    verify_fscanf(reader, 1, "inversion_sink_type_c");
    dilution_type_si[nbs][3].assign(readin);
    if(myid==0) std::cout << "inversion_sink_type_c = " << readin << std::endl;
    reader = fscanf(infile, "inversion_sink_number_c = %zu\n", &(dilution_size_si[nbs][3]));
    verify_fscanf(reader, 1, "inversion_sink_number_c");
    if(myid==0) std::cout << "inversion_sink_number_c = " << dilution_size_si[nbs][3] << std::endl;
  }

  // output path for the perambulators and the randomvectors
  reader = fscanf(infile, "outpath = %255s\n", readin);
  verify_fscanf(reader, 1, "outpath");
  outpath.assign(readin);
  if(myid==0) std::cout << "outpath = " << readin << std::endl;

  // input path for the eigensystems
  reader = fscanf(infile, "inpath_ev = %255s\n", readin);
  verify_fscanf(reader, 1, "inpath_ev");
  inpath_ev.assign(readin);
  if(myid==0) std::cout << "inpath_ev = " << readin << std::endl;
 
  // close input file
  fclose(infile);
  
  // checking input parameters
  check_input_parameters();
  check_and_create_filenames();
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
void LapH::input_parameter::print_options() {
  
  int myid = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if( myid == 0 ){
    std::cout << "omp_num_threads = " << peram_gen_omp_num_threads << std::endl;
    std::cout << "evec_read_omp_num_threads = " << evec_read_omp_num_threads << std::endl;
    std::cout << "quda_direct = " << quda_direct << std::endl;
    std::cout << "config = " <<  config << std::endl;
    std::cout << "Ls = " <<  Ls << ", Lt = " <<  Lt << std::endl;
    std::cout << "nb_ev = " <<  nb_ev << ", nb_rnd = " <<  nb_rnd << std::endl;
    for(size_t i = 0; i < nb_rnd; i++)
      std::cout << "rnd_id = " << rnd_id[i] << " seed = " <<  seed[i] << std::endl;
    std::cout << "verbose = " <<  verbose << std::endl;
    std::cout << "use_zgemm = " << use_zgemm << std::endl;
    std::cout << "endianness = " << endianness << std::endl;
    std::cout << "quarktype = " <<  quarktype << std::endl;
    std::cout << "inversion source time: " <<  dilution_type_so[0] << " " 
              <<  dilution_size_so[0] << std::endl;
    std::cout << "inversion source eigenspace: " <<  dilution_type_so[1] << " " 
              <<  dilution_size_so[1] << std::endl;
    std::cout << "inversion source Dirac: " <<  dilution_type_so[2] << " " 
              <<  dilution_size_so[2] << std::endl;
    std::cout << "number of sinks: " << nb_of_sinks << std::endl;
    for(size_t nbs = 0; nbs < nb_of_sinks; nbs++){
      std::cout << "sink number: " << nbs << std::endl;
      std::cout << "\tinversion sink time: " <<  dilution_type_si[nbs][0] << " " 
                <<  dilution_size_si[nbs][0] << std::endl;
      std::cout << "\tinversion sink spatial space: " <<  dilution_type_si[nbs][1] 
                << " " <<  dilution_size_si[nbs][1] << std::endl;
      std::cout << "\tinversion sink Dirac: " <<  dilution_type_si[nbs][2] << " " 
                <<  dilution_size_si[nbs][2] << std::endl;
      std::cout << "\tinversion sink color: " <<  dilution_type_si[nbs][3] << " " 
                <<  dilution_size_si[nbs][3] << "\n" << std::endl;
    }
    std::cout << "output path: " <<  outpath << std::endl;
    std::cout << "input path ev: " <<  inpath_ev << "\n\n" <<  std::endl;
  }

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------



