
###############################################################################
#####                                                                     #####
#####                     OPTIONS FOR THE HOST = hopper                   #####
#####                                                     May 15, 2013    #####
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
  INCPATHS += $(FFTW_INCLUDE_OPTS)
  FFTLIBS += $(FFTW_POST_LINK_OPTS)
else
  FFTLIBS +=  
endif

LIBS = $(FFTLIBS)

ifeq ($(SLEPC),yes)


   include $(SLEPC_DIR)/conf/slepc_common

   INCPATHS +=$(PETSC_FC_INCLUDES) $(SLEPC_INCLUDE)
   LIBS += $(PETSC_LIB) $(SLEPC_LIB)
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
  OPTLEVEL = 4
  FFLAGS = -Mipa=fast,inline:3
  FFLAGS += -O$(OPTLEVEL) -fastsse
endif

ifeq ($(CMPLTYPE),debug)
  FLAGS = -g -C -O0 -Mbounds -Mchkfpstk -Mchkptr -traceback
  FFLAGS +=  
endif

ifeq ($(PRECISION),double)
  FFLAGS += -r8
  LDLAGS += -r8
endif


###############################################################################
###############################################################################
