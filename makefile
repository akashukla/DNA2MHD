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
#  HOST = eh
 HOST = perl
#  HOST = stampede2
#  HOST = edison
#  HOST = hopper
# HOST = bob
# HOST = hydra
# HOST = btmac
# HOST = vip

#####  CHOSSE THE DESIRED COMPILATION TYPE 
###############################################################################

# CMPLTYPE = debug
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
   $(info $(MAKE): cd $(_D)) cd $(_D)
endef

BASEDIR = $(dir $(CURDIR)/)
HOSTDIR = $(BASEDIR)host
SRCDIR  = $(BASEDIR)src
OBJDIR  = $(BASEDIR)obj
BINDIR  = $(BASEDIR)bin2


#####  REDING THE LIBRARYS ANS COMPILERS OPTIONS DEPENDING ON THE HOST  
###############################################################################
include $(HOSTDIR)/$(HOST).mk


##### THE FILES                                                
###############################################################################
F90SRC = cc_comm.f90 \
		 cc_get_rhs_lin.f90 \
		 cc_get_rhs_nl_mpi.f90 \
		 cc_init.f90 \
		 cc_initial_condition.f90 \
		 cc_par_io.f90 \
		 cc_main.f90 \
		 cc_par_mod.f90 \
		 cc_time_advance.f90 \
		 ee_diagnostics.f90 \
		 cc_random.f90\
		 ee_mtrandom.f90
#		 cc_calc_dt.f90 \
#		 ee_diagnostics.f90 \
		 ee_eigen_direct.f90 \
		 ee_performance.f90 \
		 ee_triple_transfers.f90 \
		 ee_mtrandom.f90 \


#ifeq ($(SLEPC),yes)
#  F90SRC2 = ee_eigen_iterative.F90 \
#/work2/04943/akshukla/stampede2/DNA2MHD/src/cc_calc_dt.f90  		  ee_petsc_aux.F90 \
#  		  ee_slepc_aux.F90
#endif

F90OBJ  = $(F90SRC:.f90=.o)
F90OBJ2 = $(F90SRC2:.F90=.o)
OBJLIST = $(OBJDIR)/cc_comm.o \
                 $(OBJDIR)/cc_get_rhs_lin.o \
                 $(OBJDIR)/cc_get_rhs_nl_mpi.o \
                 $(OBJDIR)/cc_init.o \
                 $(OBJDIR)/cc_initial_condition.o \
                 $(OBJDIR)/cc_par_io.o \
                 $(OBJDIR)/cc_main.o \
                 $(OBJDIR)/cc_par_mod.o \
                 $(OBJDIR)/cc_time_advance.o \
                 $(OBJDIR)/ee_diagnostics.o \
                 $(OBJDIR)/cc_random.o \
                 $(OBJDIR)/ee_mtrandom.o

##### THE DEPENDENCIES                                                
###############################################################################
#ifeq ($(SLEPC),yes)
#  ee_petsc_aux.o:	cc_par_mod.o cc_field_solver.o cc_get_rhs_lin.o 
#  ee_slepc_aux.o:	cc_par_mod.o ee_petsc_aux.o 
#  ee_eigen_iterative.o:	cc_par_mod.o ee_petsc_aux.o ee_slepc_aux.o cc_get_rhs_lin.o \
#				cc_initial_condition.o cc_field_solver.o
#  cc_calc_dt.o:	cc_par_mod.o ee_eigen_iterative.o  cc_get_rhs_nl_mpi.o ee_eigen_direct.o
#else
#cc_calc_dt.o:	cc_par_mod.o  cc_get_rhs_nl_mpi.o #ee_eigen_direct.o
#endif

$(OBJDIR)/cc_comm.o:	$(OBJDIR)/cc_par_mod.o
#cc_field_solver.o:	cc_par_mod.o cc_flr.o cc_aux_func.o cc_hk.o
#cc_flr.o:	cc_par_mod.o cc_aux_func.o
#cc_hk.o:	cc_par_mod.o cc_aux_func.o cc_comm.o
$(OBJDIR)/cc_get_rhs_lin.o:	$(OBJDIR)/cc_par_mod.o $(OBJDIR)/cc_random.o #cc_flr.o cc_hk.o
$(OBJDIR)/cc_get_rhs_nl_mpi.o:	$(OBJDIR)/cc_par_mod.o #cc_field_solver.o
$(OBJDIR)/cc_init.o:	$(OBJDIR)/cc_par_mod.o $(OBJDIR)/ee_diagnostics.o #cc_calc_dt.o #ee_diagnostics.o \
                    #cc_hk.o cc_gaussquadrature.o #cc_field_solver.o cc_flr.o 
$(OBJDIR)/cc_initial_condition.o:	$(OBJDIR)/cc_par_mod.o $(OBJDIR)/cc_par_io.o $(OBJDIR)/ee_mtrandom.o $(OBJDIR)/ee_diagnostics.o
$(OBJDIR)/cc_par_io.o:	$(OBJDIR)/cc_par_mod.o #cc_gaussquadrature.o
#cc_gaussquadrature.o:	cc_par_mod.o
$(OBJDIR)/cc_time_advance.o:	$(OBJDIR)/cc_par_mod.o $(OBJDIR)/cc_get_rhs_lin.o $(OBJDIR)/cc_get_rhs_nl_mpi.o $(OBJDIR)/ee_diagnostics.o 

$(OBJDIR)/cc_main.o:			$(OBJDIR)/cc_par_mod.o $(OBJDIR)/cc_comm.o $(OBJDIR)/cc_init.o \
				$(OBJDIR)/cc_par_io.o $(OBJDIR)/cc_time_advance.o  \
				$(OBJDIR)/cc_get_rhs_nl_mpi.o $(OBJDIR)/ee_diagnostics.o  #cc_hk.o cc_flr.o  ee_diagnostics.o `ee_performance.o  ee_triple_transfers.o

$(OBJDIR)/ee_diagnostics.o:	$(OBJDIR)/cc_par_mod.o $(OBJDIR)/cc_get_rhs_lin.o $(OBJDIR)/cc_get_rhs_nl_mpi.o \
				$(OBJDIR)/cc_par_io.o #cc_flr.o cc_hk.o ee_Gyro_LES.o #cc_field_solver.o 
#ee_eigen_direct.o:	cc_par_mod.o cc_get_rhs_lin.o #cc_field_solver.o
#ee_performance.o:	cc_par_mod.o cc_get_rhs_lin.o cc_get_rhs_nl_mpi.o \
				cc_time_advance.o #cc_field_solver.o 
#ee_triple_transfers.o:	cc_par_mod.o cc_get_rhs_nl_mpi.o #cc_field_solver.o 			
#ee_Gyro_LES.o:	cc_field_solver.o cc_par_mod.o cc_get_rhs_nl_mpi.o \
				cc_par_io.o cc_flr.o


#####  COMPILING THE CODE 
###############################################################################
all: $(EXEC) 

$(EXEC): directory $(OBJLIST) 
	$(LD) $(LDLAGS) $(PREPROC) -o $(BINDIR)/$(EXEC) $(OBJLIST) $(LIBS)
	@echo !!!!!!!!!!!!!!!!!!!!!! SUCCESS.

directory:
	@echo  !!!!!!!!!!!!!!!!!!!!!! INITIAL DIRECTORY: $(shell pwd)
	test -d $(BINDIR) || mkdir -p $(BINDIR)
	test -d $(OBJDIR) || mkdir -p $(OBJDIR)
	@echo  !!!!!!!!!!!!!!!!!!!!!! COMPILING DIRECTORY: $(shell pwd)	

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	cd $(OBJDIR) && $(FC) $(FFLAGS) $(PREPROC) $(INCPATHS) -c -o $@ $<

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
