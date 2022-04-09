
###############################################################################
#####                                                                     #####
#####                     OPTIONS FOR THE HOST = edison                   #####
#####                                                     May 10, 2013    #####
###############################################################################


#####  LIBRARIES PATHS
###############################################################################
###MKLPATH = 
###FFTPATH = 
###PETSC_DIR = 
###SLEPC_DIR = 


#####  LIBRARIES AND INCLUDE FILES
###############################################################################
		  
ifeq ($(PRECISION),double)
  INCPATHS += -I$(FFTW_INC)
  FFTLIBS += -L$(FFTW_DIR) -lfftw3 -lfftw3f
else
  FFTLIBS +=  
endif

LIBS = $(FFTLIBS)

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
FC = ftn
LD = $(FC) 


#####  COMPILING OPTIONS
###############################################################################
ifeq ($(CMPLTYPE),optim)
  FFLAGS   = 
endif

ifeq ($(CMPLTYPE),debug)
  OPT = 
  FFLAGS +=  
endif

ifeq ($(PRECISION),double)
  FFLAGS += -fdefault-real-8
  LDLAGS += -fdefault-real-8
endif


###############################################################################
###############################################################################
