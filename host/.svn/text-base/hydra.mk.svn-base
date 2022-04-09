###############################################################################
#####                                                                     #####
#####                     OPTIONS FOR THE HOST = hydra                    #####
#####                                                         17/05/2013  #####
###############################################################################


#####  LIBRARIES PATHS
###############################################################################


#####  LIBRARIES AND INCLUDE FILES
###############################################################################
		  
INCPATHS = -I$(FFTW_HOME)/include
LIBS = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f

ifeq ($(SLEPC),yes)
   include $(SLEPC_DIR)/conf/slepc_common
   INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
#   INCPATHS += $(SLEPC_INCLUDE) $(PETSC_INC) -I$(PETSC_DIR)/$(PETSC_ARCH)/include
#  compiling with SLEPc requires libsvml to avoid unresolved dependencies
   LIBS += -L$(IFORT_BASE)/compiler/lib/intel64 -lsvml
   LIBS += $(SLEPC_LIB)
   PREPROC += -DWITHSLEPC
endif

INCPATHS += -I$(MKLROOT)/include/intel64/lp64 -I$(MKLROOT)/include
#LIBS += -L$(MKLROOT)/lib/intel64  -L$(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core  -lpthread -lm
LIBS += $(MKLROOT)/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a $(MKLROOT)/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm
 
#####  COMPILERS AND LINKER TYPE
###############################################################################
FC = mpiifort
LD = $(FC) 



#####  COMPILING OPTIONS
###############################################################################
ifeq ($(CMPLTYPE),optim)
  FFLAGS   = -stand f03 -diag-disable 7416,7025 -O3 -march=corei7-avx -xAVX -ip -g
  LDFLAGS  = $(FFLAGS)
endif

ifeq ($(CMPLTYPE),debug)
  FFLAGS   = 
  LDFLAGS += $(FFLAGS)
endif

ifeq ($(PRECISION),double)
  FFLAGS += -r8 
  LDLAGS += -r8
endif


###############################################################################
###############################################################################
