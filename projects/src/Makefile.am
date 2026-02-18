AUTOMAKE_OPTIONS = subdir-objects

#SUBDIRS = Powell

bin_PROGRAMS = \
	thor H_min max_iv_match max_iv_uc max_iv_dnu vlad vlad_2 vlad_3 vlad_4 \
	vlad_5 sls_2 diamond_ii delta delta_lo delta_lo2 quad2d quad3d \
	tune_conf tune_conf_2 tune_conf_3 tbs get_dnu mpole_scheme \
	poincare_map chk_drv_terms sxt_resp_mat max_4u max_4u_2 max_4u_3 \
	max_4u_4 max_4u_conj_grad_3 \
	gustavo gustavo_s
#               sext_resp_mat

thor_SOURCES               = thor.cc
H_min_SOURCES              = H_min.cc
max_iv_match_SOURCES       = max_iv_match.cc
max_iv_uc_SOURCES          = max_iv_uc.cc
max_iv_dnu_SOURCES         = max_iv_dnu.cc
vlad_SOURCES               = vlad.cc
vlad_2_SOURCES             = vlad_2.cc
vlad_3_SOURCES             = vlad_3.cc
vlad_4_SOURCES             = vlad_4.cc
vlad_5_SOURCES             = vlad_5.cc
sls_2_SOURCES              = sls_2.cc
diamond_ii_SOURCES         = diamond_ii.cc
delta_SOURCES              = delta.cc
delta_lo_SOURCES           = delta_lo.cc
delta_lo2_SOURCES          = delta_lo2.cc
quad2d_SOURCES             = quad2d.cc
quad3d_SOURCES             = quad3d.cc
tune_conf_SOURCES          = tune_conf.cc
tune_conf_2_SOURCES        = tune_conf_2.cc
tune_conf_3_SOURCES        = tune_conf_3.cc
tbs_SOURCES                = tbs.cc
get_dnu_SOURCES            = get_dnu.cc
mpole_scheme_SOURCES       = mpole_scheme.cc param_type.cc
poincare_map_SOURCES       = poincare_map.cc
chk_drv_terms_SOURCES      = chk_drv_terms.cc
sxt_resp_mat_SOURCES       = sxt_resp_mat.cc
max_4u_SOURCES             = max_4u.cc param_type.cc
max_4u_2_SOURCES           = max_4u_2.cc param_type.cc
max_4u_3_SOURCES           = max_4u_3.cc param_type.cc
max_4u_4_SOURCES           = max_4u_4.cc param_type.cc
max_4u_conj_grad_3_SOURCES = max_4u_3.cc param_type.cc
gustavo_SOURCES            = gustavo.cc param_type.cc
gustavo_s_SOURCES            = gustavo_s.cc param_type.cc
# sext_resp_mat_SOURCES      = sext_resp_mat.c

LOCAL = /usr/local

LIBS =	$(THOR_LIB)/thor/lib/libthor.a \
	$(THOR_LIB)/thor/lib/libTPSALib.a $(THOR_LIB)/thor/lib/libLieLib.a \
	$(NUM_REC)/lib/libnum_rec.a \
	-L$(LOCAL)/opt/gsl/lib -lgsl -lgslcblas \
	-L$(LOCAL)/opt/gfortran/lib/gcc/current -lgfortran

# C Preprocessor flags.
AM_CPPFLAGS = \
	-I$(THOR_LIB)/thor/inc \
	-I$(NUM_REC)/inc \
	-I$(LOCAL)/opt/boost/include

# FORTRAN 77 flags (The package maintainer’s variant of FFLAGS).
AM_FFLAGS = \
	-g -Og -fbacktrace -Wall -fbounds-check-fmax

# C++ flags.
AM_CXXFLAGS = -g -Og -Wall -fPIC -std=c++11

# Allocate on heap & set stack size to 1 GB.
#AM_LDFLAGS = -F/System/Library/Frameworks -framework CoreServices

UNLIMIT_RESOURCES = datasize memorysize stacksize
# System limit:
#   ulimit -a
# Set stack size for shell to hard limit (can only be changed once):
#   ulimit -s hard
#   ulimit -s
