###############################################################################
#####                                                                     #####
#####                     OPTIONS FOR THE HOST = perlmutter               #####
#####                                                     March 11, 2013  #####
###############################################################################

#####  LIBRARIES PATHS
###############################################################################
###MKLPATH = 
###FFTPATH = 
###PETSC_DIR = 
###SLEPC_DIR = 


#####  LIBRARIES AND INCLUDE FILES
###############################################################################
# MKLLIBS = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
MKLLIBS = 


ifeq ($(PRECISION),double)
  FFTLIBS +=  -Wl,-rpath,$(FFTW_ROOT)/lib -L$(FFTW_ROOT)/lib -lfftw3_mpi -lfftw3
else
  FFTLIBS +=  -Wl,-rpath,$(FFTW_ROOT)/lib -L$(FFTW_ROOT)/lib -lfftw3_mpif -lfftw3f
endif

LIBS = $(MKLLIBS) $(FFTLIBS)
INCPATHS = -I$(FFTW_ROOT)/include

ifeq ($(SLEPC),yes)

   include $(SLEPC_DIR)/conf/slepc_common
   INCPATHS += $(SLEPC_INCLUDE)
   INCPATHS += -I${PETSC_DIR}/include -I${PETSC_DIR} -I${PETSC_DIR}/${PETSC_ARCH}/include
   LIBS += $(SLEPC_LIB)
   PREPROC += -DWITHSLEPC

endif

#####  PREPROCESSING OPTIONS
###############################################################################
PREPROC =

#####  COMPILERS AND LINKER TYPE
###############################################################################
FC = mpifort
LD = $(FC) 


#####  COMPILING OPTIONS
###############################################################################
ifeq ($(CMPLTYPE),optim)
  FFLAGS   = -O0
endif

ifeq ($(CMPLTYPE),debug)
  FFLAGS = -g -O0 -fdebug-info-for-profiling -Weverything -gdwarf-5
  LDLAGS += -g -O0 -fdebug-info-for-profiling -Weverything -gdwarf-5
endif

# use PrgEnv-aocc
ifeq ($(PRECISION),double)
  FFLAGS += -r8
  LDLAGS += -r8
endif


###############################################################################
###############################################################################
