# Makefile.inc example. Works with the GNU compilers & OpenMPI.

# Flags to be passed to compilers (C++, C, Fortran, mex)
# CXXFLAGS is mandatory. Flags for the C++ compiler.
# CFLAGS is optional (used only for the C and Fortran interface).
# FFLAGS is optional (used only for the Fortran interface, or with -DHQR). 
# MEXFLAGS is optional (used only for the Octave/Matlab interface).
CXXFLAGS  = -O3 #-DHQR #-DRANDGEN -DRANGENNORMAL -std=c++11
CFLAGS    = -O3
FFLAGS    = -O3
MEXFLAGS  = --mex

# Compilers  (C++, C, Fortran, mex)
# CXX is mandatory. C++ compiler.
# CC is optional (used only for the C and Fortran interface).
# FC is optional (used only for the Fortran interface, or with -DHQR). 
# MEXis optional (used only for the Octave/Matlab interface).
CXX       = mpic++
CC        = mpicc
FC        = mpifort
MEX       = mkoctfile

# Libraries:
# LIB should provide ScaLAPACK, BLACS, LAPACK, BLAS
# LIBCXX is used for the Fortran interface, to provide C++ libraries to the linker
# LIBMEX is used for the mex interface
LIB       = -lscalapack -lblas -lm
LIBCXX    = -lstdc++ -lmpi_cxx
LIBMEX    = -L/usr/lib/openmpi -lmpi_cxx -lscalapack -lblas
