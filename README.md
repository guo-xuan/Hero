# HERO
HERO is a multi-threaded and multiprocess distributed sequencing read correction tool. It calculates the alignments between possible overlapped reads and calls the consensus. The advantage of HERO is that it can correct sequencing errors according to more reliable alignment than K-mers and also can estimate the sequencing coverage.

### Currrent Version
* v1.0

### Basic Dependencies

1. GNU GCC with C++11 support i.e. gcc4.9+ or above
2. MPI Library with MPI-3 support i.e. OpenMPI 1.8 and above or cray-mpich/7.4.0 and above. By default the mpic++ wrapper is needed. If you are on a Cray cluster and the wrapper is "CC". You will need to edit the compiler.mk file. Uncomment the line "CC := CC" and comment out "CC := mpic++".  

### Installation Steps
1. Download the tarball with compiled executables for Linux with GCC 4.9 and above from  [https://github.com/abiswas-odu/Disco/releases](https://github.com/abiswas-odu/Disco/releases). The code has been tested only on Linux and compiled with GCC4.9 and opemnpi 1.8.4.
2. If you decide to download the source code, use the following commands to build:
  1. OpenMP version "make openmp". This is also the default make option.  
  2. MPI distributed version "make mpi-dist-comp" 
If compiled successfully, the required executables will be built. 

### Running The HERO

There are two basic versions of the HERO, one for running on a single machine and another for running with MPI on a cluster.  