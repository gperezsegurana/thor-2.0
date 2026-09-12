thor-2.0

Author: Johan Bengtsson

Self-Consistent Symplectic Integrator for charged particle beam dynamics,
based on TPSA (Truncated Power Series Algebra), aka PTC (Polymorphic Tracking
Code), originated 1994; by implementing a transparent polymorphic number object
with reference counting for FP/TPSA in C++.

The symplectic integrator for realistic modeling of magnetic lattices for
ring-based synchrotrons was initially implemented in Pascal, by the author,
with care taken for the software architecture and resulting records/modules
(-> "objects") to reflect the structure of the mathematical objects describing
the underlying beam dynamics model.


Requirements:

   GNU C/C++ and FORTRAN-95 compilers: gcc and gfortran.
   GNU Scientific Library GSL.
   Boost headers.
   GNU autoconf/automake environment and libtool.
   "Numerical Recipes in C": http://www.nr.com.

On macOS:

   brew install autoconf automake libtool gcc gsl boost

On Debian/Ubuntu:

   sudo apt install build-essential gfortran autoconf automake libtool libgsl-dev libboost-dev

Build:

   export THOR_LIB=$PWD
   export NUM_REC=/path/to/numrec
   ./make_thor-2.0.sh

The library build installs into ``thor/lib``. Build the project programs in a
second step:

   cd projects
   ./make_projects.sh

Both scripts accept extra ``configure`` arguments. If GSL is not on ``PATH``,
pass its prefix, for example ``./make_thor-2.0.sh --with-gsl=/opt/homebrew/opt/gsl``.
Homebrew Boost is found automatically. For another Boost installation, pass
its prefix to the projects build, e.g. ``./make_projects.sh --with-boost=DIR``.
The projects configure step checks ``THOR_LIB`` and ``NUM_REC`` before building
and reports the GSL installation it found.
