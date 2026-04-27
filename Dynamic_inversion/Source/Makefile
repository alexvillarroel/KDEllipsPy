#
#	Make file for program NA.f
#

# setup NA macros for compilers and archiver

include makena.macros	# Edit this file to change location or
                        # type of compilers (e.g. to gcc, g77)


MPI = -DMPI=1           # MPI option (only used if TYPE=mpi)


#---- start user specific make --------------------------------------------

# edit the following to compile user routines and main program

DIR = fd3d_subs

USER_SUBS = fd3d

fd3d: $(RFIOBJ)
	cd fd3d_subs; make -f Makefile


# user supplied subroutines called by user interface routines

# must be edited for user problem

RFIOBJ = $(DIR)/fd3d_theo.o  $(DIR)/deriv.o $(DIR)/fsource.o \
      $(DIR)/fft2cd.o $(DIR)/four1.o $(DIR)/XAPIIR.o \
      $(DIR)/forward_modelling.o $(DIR)/calcmisfit.o $(DIR)/display_model.o \
      $(DIR)/display_final.o $(DIR)/ll2km.o $(DIR)/convm_dynsource.o \
      $(DIR)/output_summary.o $(DIR)/lofw.o \
      $(DIR)/mkpeak.o $(DIR)/mkstress.o
 

SYNTHOBJ = $(DIR)/fd3d_theo.o  $(DIR)/deriv.o $(DIR)/fsource.o \
      $(DIR)/forward_modelling.o $(DIR)/calcmisfit.o \
      $(DIR)/fft2cd.o $(DIR)/four1.o $(DIR)/XAPIIR.o \
      $(DIR)/ll2km.o  $(DIR)/convm_dynsource.o $(DIR)/lofw.o\
      $(DIR)/mkpeak.o $(DIR)/mkstress.o


USERFILES = $(RFIOBJ)

# compile example program
# edit the following to link in user routines

fd3d_na: fd3d
	$(F77) $(FFLAGS) -o fd3d_na fd3d_na.f \
	na.o $(USERFILES) 

fd3d_na_mpi: fd3d na_mpi
	$(MF77) $(MPI) $(FFLAGS) -o fd3d_na_mpi fd3d_na.f \
	na.o $(USERFILES) 


fd3d_fixed: fd3d
	$(F77) $(FFLAGS) -o fd3d_fixed fd3d_fixed.f \
	NA_src/na.o $(USERFILES) 


fd3d_direct: fd3d_direct.f $(SYNTHOBJ)
	$(F77) $(FFLAGS) -o fd3d_direct fd3d_direct.f $(SYNTHOBJ)

#	NA_src/na.o $(USERFILES) 
all:: $(USER_SUBS) na_mpi fd3d_na_mpi mes_mpi

# remove user object files and libraries

cleanuser:
	cd fd3d_subs; \rm -f *.o 
#
#---- end user specific make ----------------------------------------------
#
# do not change below
#
na:     
	$(F77) $(FFLAGS) -c NA_src/na.f
	mv na.o NA_src

na_mpi:
	$(MF77) $(MPI) $(FFLAGS) -c NA_src/na.F


# remove all object files and libraries 

distclean: clean

clean: cleanuser
#	\rm -f ../bin/*
	\rm *.o
	\rm -f lib/*

mes:
	echo "make all completed"


mes_mpi:
	echo
	echo "make all completed (MPI version) "
	echo
