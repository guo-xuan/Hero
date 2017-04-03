# HERO
HERE is a multi threaded and multiprocess distributed sequencing read correction tool.

### Currrent Version
* v1.0

### Setup and Installation

#### Basic Dependencies

1. GNU GCC with C++11 support i.e. gcc4.9+ or above
2. MPI Library with MPI-3 support i.e. OpenMPI 1.8 and above or cray-mpich/7.4.0 and above. By default the mpic++ wrapper is needed. If you are on a Cray cluster and the wrapper is "CC". You will need to edit the compiler.mk file. Uncomment the line "CC := CC" and comment out "CC := mpic++".  
