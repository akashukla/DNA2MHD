
###############################################################################
#####                                                                     #####
#####                     OPTIONS FOR THE HOST = stampede                 #####
#####                                                     March 11, 2013  #####
###############################################################################


#####  LIBRARIES PATHS
###############################################################################
###MKLPATH = 
FFTPATH = /usr/lib
###PETSC_DIR = 
###SLEPC_DIR = 


#####  LIBRARIES AND INCLUDE FILES
###############################################################################
ifeq ($(PRECISION),double)
  FFTLIBS +=  -lfftw3_mpi -lfftw3 -lm -L$(FFTPATH)
else
  FFTLIBS +=  -lfftw3_mpif -lfftw3f -lm  -L$(FFTPATH)
endif

LIBS = $(FFTLIBS)
INCPATHS = -I/usr/include

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
  # FFLAGS += -r8 
  # LDLAGS += -r8
endif


###############################################################################
###############################################################################
