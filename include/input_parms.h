#ifndef INPUT_PARMS_H_
#define INPUT_PARMS_H_

#include <cstring> 
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <fstream>
#include <string>
#include <vector>

#include "mpi.h"

#include "macros.h"

namespace LapH {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
class input_parameter{
public:
  // ---------------------- METHODS ------------------------------------------
  // standard constructor, does not do anything
  input_parameter() {
    rnd_id = NULL;
    seed = NULL;
    use_zgemm = false;
  }; 
  // copy constructor
  input_parameter(const input_parameter& other);
  // assignment operator
  input_parameter& operator=(const input_parameter& other);
  // standard destructor, does not do anything
  ~input_parameter() {
    delete[] rnd_id;
    delete[] seed;
  }; 

public:
  // --------------------- INTERFACE -----------------------------------------
  void parse_input_file(int argc, char *argv[]);
  void print_options();
  // -------------------- END INTERFACE --------------------------------------

private:
  // checks the input parameters for consistency
  void check_input_parameters();
  // creates the filenames for the perambulator and the random vectors and
  // checks if these files already exist. If they exist the program is
  // terminated to avoid overwriting these files!
  void check_and_create_filenames();

  /**
   * @brief check if the number of elements read by fscanf corresponds to expectations
   *
   * @param nelem_actual number of elements actually returned by fscanf
   * @param nelem_expected number of elements that was expected
   * @param parname string constant name of the parameter
   */
  void verify_fscanf(int const nelem_actual, int const nelem_expected, char const * const parname);

public:
  // -------------------------- DATA -----------------------------------------
  size_t dilution_size_so[3];      // SOURCE dilution sizes and types
  std::string dilution_type_so[3]; // 0 = time, 1 = eigevector space, 2 = Dirac
  size_t nb_of_sinks;              // number of sinks
  std::vector<size_t> 
               nb_of_sink_rnd_vec; // number of sink random vectors
  std::vector<
        std::vector<size_t> > 
                 dilution_size_si; // SINK dilution sizes and types
  std::vector<
        std::vector<std::string> >
                 dilution_type_si; // This is a vector now! Second dimension is:
                                   // 0 = time, 1 = space, 2 = Dirac, 3 = color
                                   // NOTE: Please use capital Latters I, B, F!!
                                   // All of these numbers are NOT the blocksize
                                   // or interlace sizes, but the number of 
                                   // inversions in the corresponding space 
  size_t config;                   // configuration number
  size_t nb_config;                // total number of configs to be processed
  size_t delta_config;             // distance between configurations
  size_t Ls;                       // spatial extend
  size_t Lt;                       // temporal extend
  size_t nb_ev;                    // number of eigenvectors
  size_t verbose;                  // displaying more informations

  bool use_zgemm;                  // for add_to_perambulator:  (true) -> use zgemm for multiplication with Vdagger
                                   //                          (false) -> use eigen multiplication

  bool hack_clean;                 // when running in a reduced memory environment, perambulator writing
                                   // may fail due to a lack of memory for the write buffer for the perambulator
                                   // as a workaround for this, we offer to free all memory from the distillery
                                   // at the cost of being able to do just a single random vector per job
                                   // this cost is acceptable for lattices up to 32c64 as eigenvector reading
                                   // only becomes a truly significant chunk of total time beyond that point

  // when tmLQCD is compiled without QUDA support, invert_quda_direct is not going to be defined
  // so we use the PG_QUDA_DIRECT preprocessor constant to enable and disable support for this
#ifdef PG_QUDA_DIRECT
  bool quda_direct;
#endif

  unsigned int peram_gen_omp_num_threads; // number of OpenMP threads used by peram_gen (1 or 2)
  unsigned int eigen_omp_num_threads;     // number of OpenMP threads used by eigen
  unsigned int evec_read_omp_num_threads; // number of OpenMP threads used for reading eigenvectors 
  std::string endianness;           // needed for endianess change 

  size_t nb_rnd;                    // random vector id of vector in process
  int* rnd_id;                      // id of random vectors
  int* seed;                        // seeds for random vector generation
  std::vector<std::vector<int> >
                          seed_si;  // seeds for sink random vector generation

  std::string quarktype;            // quark type: u,d,s,c -> for naming outfiles
  std::string outpath;              // path to write everything
  std::string inpath_ev;            // path to read eigenvectors
  std::string peram_file_name;      // perambulator file name
  std::string rnd_vec_file_name;    // random vektor file name
  // ------------------------- END DATA ---------------------------------------

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

} // namespace ends here

#endif /* INPUT_PARMS_H_ */

