
###############################################################################
#####                                                                     #####
#####                     OPTIONS FOR THE HOST = stampede                 #####
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
MKLLIBS = -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
		  
ifeq ($(PRECISION),double)
  FFTLIBS +=  -Wl,-rpath,$(TACC_FFTW3_LIB) -L$(TACC_FFTW3_LIB) -lfftw3_mpi -lfftw3
else
  FFTLIBS +=  -Wl,-rpath,$(TACC_FFTW3_LIB) -L$(TACC_FFTW3_LIB) -lfftw3_mpif -lfftw3f
endif

LIBS = $(MKLLIBS) $(FFTLIBS)
INCPATHS = -I$(TACC_FFTW3_INC)

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
endif

ifeq ($(CMPLTYPE),debug)
  OPT = -O0
  FFLAGS += -g 
endif

ifeq ($(PRECISION),double)
  FFLAGS += -r8 
  LDLAGS += -r8
endif


###############################################################################
###############################################################################
