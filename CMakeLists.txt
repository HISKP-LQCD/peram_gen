cmake_minimum_required(VERSION 3.8)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
message(STATUS "CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
project(peram_gen LANGUAGES C CXX Fortran)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CNAME_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED OIN)

set(CMAKE_C_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Wall")

set(CMAKE_C_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -ggdb")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -g -ggdb")

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(external_libraries)

find_package(Eigen3 REQUIRED NO_MODULE)
if(EIGEN3_FOUND)
  message(STATUS "Eigen3 found")
  include_directories(${EIGEN3_INCLUDE_DIR})
else()
  message(FATAL_ERROR "Eigen3 is required to build peram_gen!")
endif()

# tmLQCD links against LAPACK, so we have to do this too
set(LAPACK_HOME CACHE PATH "Path to custom LAPACK library")
if(NOT LAPACK_HOME)
  find_package(LAPACK REQUIRED)
  if( LAPACK_FOUND )
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
    list(APPEND external_libraries "${LAPACK_LIBRARIES}" )
    message(STATUS "Using LAPACK libraries: ${LAPACK_LIBRARIES}")
  else()
    message( FATAL_ERROR "Unable to find LAPACK!" )
  endif()
else()
  link_directories("${LAPACK_HOME}/lib")
  set(LAPACK_LIBRARIES CACHE STRING "libraries for custom LAPACK (e.g. '-lsci_gnu_mp' or '-llapack')")
  if(LAPACK_LIBRARIES)
    list(APPEND external_libraries "${LAPACK_LIBRARIES}")
    message(STATUS "Using LAPACK libraries: ${LAPACK_LIBRARIES}")
  else()
    message(FATAL_ERROR "When using custom LAPACK libraries, -DLAPACK_HOME and -DLAPACK_LIBRARIES must be provided")
  endif()
endif()

# tmLQCD links against BLAS, so we have to do this too
set(BLAS_HOME CACHE PATH "Path to custom BLAS library")
if(NOT BLAS_HOME)
  find_package(BLAS REQUIRED)
  if( BLAS_FOUND )
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${BLAS_LINKER_FLAGS}")
    list(APPEND external_libraries "${BLAS_LIBRARIES}" )
    message(STATUS "Using BLAS libraries: ${BLAS_LIBRARIES}")
  else()
    message( FATAL_ERROR "Unable to find BLAS!" )
  endif()
else()
  link_directories("${BLAS_HOME}/lib")
  set(BLAS_LIBRARIES CACHE STRING "libraries for custom LAPACK (e.g. '-lsci_gnu_61_mp' or '-lblas')")
  if(BLAS_LIBRARIES)
    list(APPEND external_libraries "${BLAS_LIBRARIES}")
    message(STATUS "Using BLAS libraries: ${BLAS_LIBRARIES}")
  else()
    message(FATAL_ERROR "When using custom BLAS libraries, -DBLAS_HOME and -DBLAS_LIBRARIES must be provided")
  endif()
endif()

find_package(OpenMP REQUIRED)
if ( OPENMP_FOUND )
  add_definitions(-DHAVE_OPENMP)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  list(APPEND external_libraries "${OPENMP_LIBRARIES}" )
else()
  message( FATAL_ERROR "Unable to find OpenMP!" )
endif()

find_package(Boost
             REQUIRED
             program_options
             filesystem)
if( Boost_FOUND )
  link_directories( ${Boost_LIBRARY_DIRS} )
  include_directories( ${Boost_INCLUDE_DIRS} )
  link_libraries( ${Boost_LIBRARIES} )
else()
  message( FATAL_ERROR "Unable to find Boost (or one of the requested components, see CMakeLists.txt!")
endif()

set( MPI_HOME CACHE PATH "path to MPI library installation directory" )
if( MPI_HOME )
  add_definitions(-DHAVE_MPI)
  add_definitions(-DOMPI_SKIP_MPICXX)
  include_directories(${MPI_HOME}/include)
  link_directories(${MPI_HOME}/lib)
  list(APPEND external_libraries "-lmpi")
else()
  find_package(MPI)
  if( MPI_C_FOUND )
    add_definitions(-DHAVE_MPI)
    # we never want to have the C++ MPI bindings!
    add_definitions(-DOMPI_SKIP_MPICXX)
    include_directories(${MPI_C_INCLUDE_PATH})
    list(APPEND external_libraries "${MPI_C_LIBRARIES}")
    message(STATUS "MPI found")
    message(STATUS "  MPI_C_INCLUDE_PATH=${MPI_C_INCLUDE_PATH}")
    message(STATUS "  MPI_C_LIBRARIES=${MPI_C_LIBRARIES}")
  else()
    message( FATAL_ERROR "MPI could not be found automatically and MPI_HOME has not been specified!" )
  endif()
endif()

set( TMLQCD_SRC CACHE PATH "path to tmLQCD source directory" )
if( TMLQCD_SRC )
  include_directories( ${TMLQCD_SRC}/include )
  include_directories( ${TMLQCD_SRC} )
else()
	message( FATAL_ERROR "The tmLQCD source directory must be specified!" )
endif()

set( TMLQCD_BUILD CACHE PATH "path to tmLQCD build directory (must be build directory for source in TMLQCD_SRC!)" )
if( TMLQCD_BUILD )
  link_directories(${TMLQCD_BUILD}/lib)
  include_directories(${TMLQCD_BUILD})
  include_directories(${TMLQCD_BUILD}/include)
  list(APPEND external_libraries 
    "-lwrapper" "-lhmc" "-lmonomial"
    "-loperator" "-lmeas" "-linit" 
    "-lhmc" "-lxchange" "-lrational"
    "-lsolver" "-lio" "-llinalg" "-lm")
else()
  message( FATAL_ERROR "The tmLQCD Build directory must be specified!" )
endif()

set( LIME_HOME CACHE PATH "path to the c-lime library installation" )
if( LIME_HOME )
  link_directories(${LIME_HOME}/lib)
  include_directories(${LIME_HOME}/include)
  list(APPEND external_libraries "-llime")
else()
  message(FATAL_ERROR "LIME_HOME must be defined!" )
endif()

set( LEMON_HOME CACHE PATH "path to the LEMON library installation (optional, depends on how tmLQCD was configured)" )
if( LEMON_HOME )
  add_definitions(-DHAVE_LIBLEMON)
  link_directories(${LEMON_HOME}/lib)
  include_directories(${LEMON_HOME}/include)
  list(APPEND external_libraries "-llemon")
  # some of the conversions done in LEMON are rather dubious at best
  # let's use -fpermissive to allow them to compile anyway
  #string(APPEND CMAKE_CXX_FLAGS " -fpermissive")
  message( STATUS "LEMON support enabled" )
else()
  message( STATUS "LEMON support disabled" )
endif()

set( DD_ALPHA_AMG_HOME CACHE PATH "path to the DDalphaAMG library installation (optional, depends on how tmLQCD was configured)" )
if( DD_ALPHA_AMG_HOME )
  link_directories(${DD_ALPHA_AMG_HOME}/lib)
  include_directories(${DD_ALPHA_AMG_HOME}/include)
  list(APPEND external_libraries "-lDDalphaAMG")
  message( STATUS "DDalphaAMG support enabled" )
endif()

set( QUDA_HOME CACHE PATH "path to QUDA installation (optional, depends on whether tmLQCD uses QUDA sovers)" ) 
if( QUDA_HOME )
  link_directories(${QUDA_HOME}/lib)
  include_directories(${QUDA_HOME}/include)
  list(APPEND external_libraries "-lquda")

  # enable CUDA language support
  enable_language(CUDA)
  # we still use find_package because we want the directories
  find_package(CUDA REQUIRED)
  if( CUDA_TOOLKIT_ROOT_DIR )
    link_directories( "${CUDA_TOOLKIT_ROOT_DIR}/lib" "${CUDA_TOOLKIT_ROOT_DIR}/lib64" )
    include_directories( ${CUDA_INCLUDE_DIRS} )
    # QUDA requires -lcuda for cuMemAlloc in addition to -lcudart
    list(APPEND external_libraries "-lcuda" "-lcudart" "-lcublas")
    message(STATUS "CUDA found")
    message(STATUS "  CUDA_SDK_ROOT_DIR=${CUDA_SDK_ROOT_DIR}")
    message(STATUS "  CUDA_TOOLKIT_ROOT_DIR=${CUDA_TOOLKIT_ROOT_DIR}")
    message(STATUS "  CUDA_INCLUDE_DIRS=${CUDA_INCLUDE_DIRS}")
  else()
    message( FATAL_ERROR "Unable to find CUDA!")
  endif()
  message( STATUS "compiling against tmLQCD with QUDA" )
endif()

set( SUPPORT_QUDA_DIRECT FALSE CACHE BOOL "whether to support the quda_direct pathway for inversions (can be toggled further at runtime)" )
if( SUPPORT_QUDA_DIRECT )
  if( NOT QUDA_HOME )
    message(FATAL_ERROR "The direct GPU solver is only available when tmLQCD has been compiled with QUDA support. See QUDA_HOME.")
  endif()
  add_definitions(-DPG_QUDA_DIRECT)
endif()

include(GetGitRevisionDescription)
get_git_head_revision(git_refspec git_sha1)
git_local_changes(git_changes)
message(STATUS "git sha1: ${git_sha1}")
message(STATUS "git refspec: ${git_refspec}")
message(STATUS "git state: ${git_changes}")
configure_file("${CMAKE_SOURCE_DIR}/modules/gitversion.cpp.in"
  "${CMAKE_CURRENT_BINARY_DIR}/modules/gitversion.cpp"
  @ONLY)
list(APPEND SOURCES "${CMAKE_CURRENT_BINARY_DIR}/modules/gitversion.cpp")

set( QMP_HOME CACHE PATH "path to the QMP library required by QPhiX (optional, depends on how tmLQCD was configured)" )
if( QMP_HOME )
  link_directories(${QMP_HOME}/lib)
  include_directories(${QMP_HOME}/include)
  list(APPEND external_libraries "-lqmp")
endif()

set( QPHIX_HOME CACHE PATH "path to the QPhiX library installation (optional, depends on how tmLQCD was configured)" )
if( QPHIX_HOME )
  ## this is a really ugly hack to make sure that all symbols are found, bit it works!
  ## it can only be done as the very last step because it changes the linker flags
  ## and makes many tests break
  ## the way the QPhiX interface is linked into tmLQCD should be adjusted for this
  ## to become unnecessary
  set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${TMLQCD_BUILD}/qphix_interface.o" )
  link_directories(${QPHIX_HOME}/lib)
  include_directories(${QPHIX_HOME}/include)
  list(APPEND external_libraries "-lqphix_codegen" "-lqphix_solver")
  message( STATUS "compiling against tmLQCD with QPhiX" )
endif()

include_directories(${CMAKE_SOURCE_DIR}/include)

add_library(pgen
  STATIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules/gitversion.cpp
  modules/distillery.cpp
  modules/random_vector.cpp
  modules/ranlxd.cpp
  modules/ranlxs.cpp
  modules/input_parms.cpp
)

add_executable(peram_gen main/main.cpp)
target_link_libraries(peram_gen
  pgen
	${external_libraries})

