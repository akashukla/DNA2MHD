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
  FFTLIBS +=  -L$(P3DFFT_AOCC)/lib -L$(FFTW_ROOT)/lib -lp3dfft -lfftw3_mpi -lfftw3
else # no single precision option
  -L$(P3DFFT_AOCC)/lib -L$(FFTW_ROOT)/lib -lp3dfft -lfftw3_mpi -lfftw3
endif

LIBS = $(MKLLIBS) $(FFTLIBS)
INCPATHS =  -I$(P3DFFT_AOCC)/include -I$(FFTW_ROOT)/include

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
FC = mpif90
LD = $(FC) 


#####  COMPILING OPTIONS
###############################################################################
ifeq ($(CMPLTYPE),optim)
  FFLAGS   = -O2
LDLAGS += -O2
endif

ifeq ($(CMPLTYPE),debug)
  FFLAGS = -g -O2 -fdebug-info-for-profiling -Weverything -gdwarf-5
  LDLAGS += -g -O2 -fdebug-info-for-profiling -Weverything -gdwarf-5
endif

# use PrgEnv-aocc
ifeq ($(PRECISION),double)
  FFLAGS += -r8
  LDLAGS += -r8
endif


###############################################################################
###############################################################################
