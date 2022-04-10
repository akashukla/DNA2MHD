###############################################################################
#####                                                                     #####
#####                          DNA MAIN MAKEFILE                          #####
#####                                                         29/12/2012  #####
###############################################################################

default	: all 


###############################################################################
#####                            USER CHOICES                             #####  
###############################################################################

#####  CHOSSE THE DESIRED COMPUTER HOST 
###############################################################################
#  HOST = ls5
  HOST = stampede2
#  HOST = edison
#  HOST = hopper
# HOST = bob
# HOST = hydra
# HOST = btmac
# HOST = vip

#####  CHOSSE THE DESIRED COMPILATION TYPE 
###############################################################################
#  CMPLTYPE = debug
  CMPLTYPE = optim

#####  CHOSSE THE DESIRED PRECISION TYPE 
###############################################################################
  PRECISION = double
# PRECISION = single

#####  CHOSSE THE USE OF A SLEPC LIBRARY 
###############################################################################
#  SLEPC = yes
 SLEPC = no  

#####  CHOSSE THE EXECUTABLE NAME 
###############################################################################
  EXEC = dna
# EXEC = rna


###############################################################################
#####                          GLOBAL PARAMETERS                          #####
###############################################################################

#####  SPECIFYING THE DIRECTORIES RELATIVE TO THE MAKEFILE LOCATION
###############################################################################
CHDIR_SHELL := $(SHELL)
define chdir
   $(eval _D=$(firstword $(1) $(@D)))
   $(info $(MAKE): cd $(_D)) $(eval SHELL = cd $(_D); $(CHDIR_SHELL))
endef

BASEDIR = $(dir $(CURDIR)/)
HOSTDIR = $(BASEDIR)host
SRCDIR  = $(BASEDIR)src
OBJDIR  = $(BASEDIR)obj
BINDIR  = $(BASEDIR)bin


#####  REDING THE LIBRARYS ANS COMPILERS OPTIONS DEPENDING ON THE HOST  
###############################################################################
include $(HOSTDIR)/$(HOST).mk


##### THE FILES                                                
###############################################################################
F90SRC = cc_comm.f90 \
		 cc_get_rhs_lin.f90 \
		 cc_get_rhs_nl.f90 \
		 cc_init.f90 \
		 cc_initial_condition.f90 \
		 cc_par_io.f90 \
		 cc_main.f90 \
		 cc_par_mod.f90 \
		 cc_time_advance.f90 \
		 cc_calc_dt.f90 \
#		 ee_diagnostics.f90 \
		 ee_eigen_direct.f90 \
		 ee_performance.f90 \
		 ee_triple_transfers.f90 \
		 ee_mtrandom.f90 \


#ifeq ($(SLEPC),yes)
#  F90SRC2 = ee_eigen_iterative.F90 \
#  		  ee_petsc_aux.F90 \
#  		  ee_slepc_aux.F90
#endif
		  
F90OBJ  = $(F90SRC:.f90=.o)
F90OBJ2 = $(F90SRC2:.F90=.o)
OBJLIST = $(F90OBJ) $(F90OBJ2)


##### THE DEPENDENCIES                                                
###############################################################################
#ifeq ($(SLEPC),yes)
#  ee_petsc_aux.o:	cc_par_mod.o cc_field_solver.o cc_get_rhs_lin.o 
#  ee_slepc_aux.o:	cc_par_mod.o ee_petsc_aux.o 
#  ee_eigen_iterative.o:	cc_par_mod.o ee_petsc_aux.o ee_slepc_aux.o cc_get_rhs_lin.o \
#				cc_initial_condition.o cc_field_solver.o
#  cc_calc_dt.o:	cc_par_mod.o ee_eigen_iterative.o  cc_get_rhs_nl.o ee_eigen_direct.o
#else
cc_calc_dt.o:	cc_par_mod.o  cc_get_rhs_nl.o #ee_eigen_direct.o
#endif

cc_comm.o:	cc_par_mod.o
#cc_field_solver.o:	cc_par_mod.o cc_flr.o cc_aux_func.o cc_hk.o
#cc_flr.o:	cc_par_mod.o cc_aux_func.o
#cc_hk.o:	cc_par_mod.o cc_aux_func.o cc_comm.o
cc_get_rhs_lin.o:	cc_par_mod.o #cc_flr.o cc_hk.o
cc_get_rhs_nl.o:	cc_par_mod.o #cc_field_solver.o
cc_init.o:	cc_par_mod.o cc_calc_dt.o #cc_calc_dt.o #ee_diagnostics.o \
                    #cc_hk.o cc_gaussquadrature.o #cc_field_solver.o cc_flr.o 
cc_initial_condition.o:	cc_par_mod.o cc_par_io.o ee_mtrandom.o
cc_par_io.o:	cc_par_mod.o #cc_gaussquadrature.o
#cc_gaussquadrature.o:	cc_par_mod.o
cc_time_advance.o:	cc_par_mod.o cc_get_rhs_lin.o cc_get_rhs_nl.o cc_calc_dt.o \
					#cc_calc_dt.o #cc_field_solver.o ee_diagnostics.o 
cc_main.o:			ee_performance.o cc_par_mod.o cc_comm.o cc_init.o \
				cc_par_io.o cc_time_advance.o  \
				cc_get_rhs_nl.o ee_triple_transfers.o  cc_calc_dt.o #cc_hk.o cc_flr.o  ee_diagnostics.o `

#ee_diagnostics.o:	cc_par_mod.o cc_get_rhs_lin.o cc_get_rhs_nl.o \
				cc_par_io.o #cc_flr.o cc_hk.o ee_Gyro_LES.o #cc_field_solver.o 
#ee_eigen_direct.o:	cc_par_mod.o cc_get_rhs_lin.o #cc_field_solver.o
#ee_performance.o:	cc_par_mod.o cc_get_rhs_lin.o cc_get_rhs_nl.o \
				cc_time_advance.o #cc_field_solver.o 
#ee_triple_transfers.o:	cc_par_mod.o cc_get_rhs_nl.o #cc_field_solver.o 			
#ee_Gyro_LES.o:	cc_field_solver.o cc_par_mod.o cc_get_rhs_nl.o \
				cc_par_io.o cc_flr.o


#####  COMPILING THE CODE 
###############################################################################
all: $(EXEC) 

$(EXEC): directory $(OBJLIST) 
	$(LD) $(LDFLAGS) $(PREPROC) -o $(BINDIR)/$(EXEC) $(OBJLIST) $(LIBS)
	@echo !!!!!!!!!!!!!!!!!!!!!! SUCCESS.
	
directory:
	@echo  !!!!!!!!!!!!!!!!!!!!!! INITIAL DIRECTORY: $(shell pwd)
	test -d $(BINDIR) || mkdir -p $(BINDIR)
	test -d $(OBJDIR) || mkdir -p $(OBJDIR)
	$(call chdir,$(OBJDIR))
	@echo  !!!!!!!!!!!!!!!!!!!!!! COMPILING DIRECTORY: $(shell pwd)	

%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(PREPROC) $(INCPATHS) -c -o $@ $< 
	
%.o: $(SRCDIR)/%.F90
	$(FC) $(FFLAGS) $(PREPROC) $(INCPATHS) -c -o $@ $< 


#####  MAKING THE CLEANING TOOLS 
###############################################################################
clean:: cl

cl: 
	rm -f $(OBJDIR)/*.o
	rm -f $(OBJDIR)/*.mod
	rm -f $(BINDIR)/$(EXEC)
	@echo !!!!!!!!!!!!!!!!!!!!!!  CLEANING DONE.


###############################################################################
###############################################################################
