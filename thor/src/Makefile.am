AUTOMAKE_OPTIONS = subdir-objects

lib_LIBRARIES = libthor.a

libthor_a_SOURCES = \
	thor_lib.cc \
	../../TPSA/TPSALib.f ../../TPSA/LieLib.f

GSL_INCLUDE = -I/usr/local/opt/gsl/include

# C Preprocessor flags.
AM_CPPFLAGS = $(GSL_INCLUDE) \
	      -I../inc -I$(NUM_REC)/inc

# FORTRAN-90 flags.
FFLAGS = \
	-g -O2 -Wall -fbounds-check -fdefault-integer-8 -mcmodel=medium \
	-w -Wno-argument-mismatch -fmax-stack-var-size=10

#AM_CXXFLAGS = -g -O0 -Wall -Wno-non-template-friend -fPIC
AM_CXXFLAGS = -g -O0 -Wall -Wno-error=all -fPIC
