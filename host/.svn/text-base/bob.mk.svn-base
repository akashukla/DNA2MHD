###############################################################################
#####                                                                     #####
#####                     OPTIONS FOR THE HOST = bob                      #####
#####                                                         10/12/2012  #####
###############################################################################


#####  LIBRARIES PATHS
###############################################################################
MKLPATH = $(MKL_HOME)/lib/intel64
FFTPATH = $(FFTW_HOME)
PETSC_DIR = /afs/ipp-garching.mpg.de/home/t/tbg/soft/@sys/petsc
SLEPC_DIR = /afs/ipp-garching.mpg.de/home/t/tbg/soft/@sys/slepc


#####  LIBRARIES AND INCLUDE FILES
###############################################################################
MKLLIBS  = -L$(MKLPATH)
MKLLIBS += -Wl,--start-group $(MKLPATH)/libmkl_scalapack_lp64.a \
		$(MKLPATH)/libmkl_blacs_intelmpi_lp64.a $(MKLPATH)/libmkl_intel_lp64.a \
		$(MKLPATH)/libmkl_sequential.a $(MKLPATH)/libmkl_core.a -Wl,--end-group \
		-lpthread 
		  
ifeq ($(PRECISION),double)
  FFTLIBS = $(FFTPATH)/lib/libfftw3.a
else
  FFTLIBS = $(FFTPATH)/lib/libfftw3f.a
endif

LIBS = $(MKLLIBS) $(FFTLIBS)
INCPATHS = -I$(FFTPATH)/include

ifeq ($(SLEPC),yes)
   PETSC_ARCH = 
   include $(SLEPC_DIR)/conf/slepc_common
   LIBS += $(SLEPC_LIB)
   INCPATHS +=$(PETSC_FC_INCLUDES)  $(SLEPC_INCLUDE)
   PREPROC += -DWITHSLEPC
endif
 
 
#####  PREPROCESSING OPTIONS
###############################################################################
PREPROC =
 
 
#####  COMPILERS AND LINKER TYPE
###############################################################################
FC = mpiifort
LD = $(FC) 


#####  COMPILING OPTIONS
###############################################################################
ifeq ($(CMPLTYPE),optim)
  FFLAGS   = -stand f03 -diag-disable 7416,7025 -O3 -march=core2 -xSSSE3 -ip
  FFLAGS  += -vec_report0  
  FFLAGS  += $(MKLFLAGS)
  LDFLAGS  = $(FFLAGS)
endif

ifeq ($(CMPLTYPE),debug)
  FFLAGS   = -O0 -g -stand f03 -C -traceback -debug inline-debug-info -warn all  
  FFLAGS  += -check all -ftrapuv -fstack-security-check -fpe0 
  LDFLAGS += $(FFLAGS)
  INCLUDE_SYMBOLS=yes
endif

ifeq ($(PRECISION),double)
  FFLAGS += -r8 
  LDLAGS += -r8
endif


###############################################################################
###############################################################################
