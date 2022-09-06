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
   $(info $(MAKE): cd $(_D)) $(cd $(_D))
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
F90SRC =	cc_comm.f90\
		cc_get_rhs_lin.f90\
	 	cc_get_rhs_nl.f90\
		cc_init.f90\
		cc_initial_condition.f90\
		cc_par_io.f90\
		cc_par_mod.f90\
		cc_time_advance.f90\
		cc_main.f90\
		ee_diagnostics.f90

F90OBJ  = $(F90SRC:.f90=.o)
OBJLIST = $(addprefix $(OBJDIR)/, $(F90OBJ))


##### THE DEPENDENCIES                                                
###############################################################################
#ifeq ($(SLEPC),yes)
#  ee_petsc_aux.o:	cc_par_mod.o cc_field_solver.o cc_get_rhs_lin.o 
#  ee_slepc_aux.o:	cc_par_mod.o ee_petsc_aux.o 
#  ee_eigen_iterative.o:	cc_par_mod.o ee_petsc_aux.o ee_slepc_aux.o cc_get_rhs_lin.o \
#				cc_initial_condition.o cc_field_solver.o
#  cc_calc_dt.o:	cc_par_mod.o ee_eigen_iterative.o  cc_get_rhs_nl.o ee_eigen_direct.o
#else
#cc_calc_dt.o:	cc_par_mod.o  cc_get_rhs_nl.o #ee_eigen_direct.o
#endif

$(OBJDIR)/cc_comm.o:	$(OBJDIR)/cc_par_mod.o
$(OBJDIR)/cc_get_rhs_lin.o:	$(OBJDIR)/cc_par_mod.o
$(OBJDIR)/cc_get_rhs_nl.o:	$(OBJDIR)/cc_par_mod.o
$(OBJDIR)/cc_init.o:	$(OBJDIR)/cc_par_mod.o\
			$(OBJDIR)/ee_diagnostics.o
$(OBJDIR)/cc_initial_condition.o:	$(OBJDIR)/cc_par_mod.o\
					$(OBJDIR)/cc_par_io.o\
					$(OBJDIR)/ee_mtrandom.o
$(OBJDIR)/cc_par_io.o:	$(OBJDIR)/cc_par_mod.o
$(OBJDIR)/cc_time_advance.o:	$(OBJDIR)/cc_par_mod.o\
				$(OBJDIR)/cc_get_rhs_lin.o\
				$(OBJDIR)/cc_get_rhs_nl.o\
				$(OBJDIR)/ee_diagnostics.o
$(OBJDIR)/cc_main.o:	$(OBJDIR)/cc_par_mod.o\
			$(OBJDIR)/cc_comm.o\
			$(OBJDIR)/cc_init.o\
			$(OBJDIR)/cc_par_io.o\
			$(OBJDIR)/cc_time_advance.o\
			$(OBJDIR)/cc_get_rhs_nl.o\
			$(OBJDIR)/ee_diagnostics.o

$(OBJDIR)/ee_diagnostics.o:	$(OBJDIR)/cc_par_mod.o\
				$(OBJDIR)/cc_get_rhs_lin.o\
				$(OBJDIR)/cc_get_rhs_nl.o\
				$(OBJDIR)/cc_par_io.o


#####  COMPILING THE CODE 
###############################################################################
all: $(EXEC) MAKE

$(EXEC): directory $(OBJLIST)
	$(LD) $(LDFLAGS) $(PREPROC) -o $(BINDIR)/$(EXEC) $(OBJLIST) $(LIBS)
	@echo !!!!!!!!!!!!!!!!!!!!!! SUCCESS.

directory:
	@echo $(SRCDIR) $(OBJDIR) $(OBJLIST)
	test -d $(BINDIR) || mkdir -p $(BINDIR)
	test -d $(OBJDIR) || mkdir -p $(OBJDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) $(PREPROC) $(INCPATHS) -c -o $@ $< 

$(OBJDIR)/%.o: $(SRCDIR)/%.F90
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
