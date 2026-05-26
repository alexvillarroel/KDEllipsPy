# # Make file for program NA.f #

# setup NA macros for compilers and archiver


#include makena.macros	# Edit this file to change location or
                        # type of compilers (e.g. to gcc, g77)

F77 = ifx -xhost -ipo -O4


MPI = -DNA_MPI=1           # MPI option (only used if TYPE=mpi)

#---- start user specific make  (Edit below) ----------------------------------

# edit the three variables below to compile user routines and main program

                         # Enter here the directory name containing all 
			 # user code other than main program
USER_DIR = kine_subs
                         # Enter here the list of object files or libraries 
			 # to be linked in with user code
USER_OBJ = \
         $(USER_DIR)/forward_modelling_ell.o \
         $(USER_DIR)/convm_new_src.o \
         $(USER_DIR)/calcmisfit.o \
         $(USER_DIR)/display_model.o \
         $(USER_DIR)/display_final.o \
         $(USER_DIR)/lofw.o \
         $(USER_DIR)/output_summary.o \
         $(USER_DIR)/read_vmodelNA.o

                         # Enter here name of user's main routine
			 # Note that the main program name is
			 # assumed to be `USER_MAIN.F' 
USER_MAIN = kine_na_ifx

user : na
	cd $(USER_DIR); make -f Makefile

all: user kine_na  kine_direct 

serial:: user na user_main 

mpi:: user na_mpi user_main_mpi utl mes_mpi

kine_na:
	$(F77) $(FFLAGS) -o ../bin/kine_na kine_na_new.f \
	lib/na.o  $(USER_OBJ)

kine_direct: kine_na_direct.f
	$(F77) $(FFLAGS) -o ../bin/kine_direct kine_na_direct.f \
	$(USER_OBJ)


# remove user object files and libraries

cleanuser:
	cd $(USER_DIR); \rm -fr *.o
	\rm -f ../bin/$(USER_MAIN)

#---- end user specific make ----------------------------------------------
#
# do not change below
#
na:     
	$(F77) $(FFLAGS) $(F77_COMPAT) -c NA_src/na.F
	mv na.o lib

na_mpi:     
	$(MF77) $(MPI) $(FFLAGS) $(F77_COMPAT) -c NA_src/na.F
	mv na.o lib

utl:: 
	cd utl/splot; make -f Makefile 
	cd utl/nad; make -f Makefile 
	cd utl/nad/nad_convert;  make -f Makefile 





# remove all object files and libraries 

distclean: clean

clean: cleanuser
	-rm -f ../bin/*
	-rm -f lib/*
	-rm -f utl/splot/Lib/*
	-rm -f utl/splot/subs/*.o
	-rm -f utl/splot/*.o
	-rm -f utl/nad/*.o
	-rm -f utl/nad/nad_convert/*.o
	-rm -f ../data/na.nad
	-rm -f ../data/na.sum
	-rm -f ../data/sobol.coeff
	-rm -f ../data/splot/*.ps
	-rm -f ../data/rfi_plot/*.ps
	-rm -f ../data/rfi_files/NA_ORF/rec_func_demo
	-rm -f ../data/rfi_files/NA_SRF/rec_func_demo
	-rm -f ./kine_subs/rfi_plot/graphic/xplot/lib/*.a
	-rm -f ./kine_subs/rfi_plot/graphic/lib/*.a
	-rm -f ./kine_subs/rfi_plot/graphic/xplot/*.o
	-rm -f ./kine_subs/rfi_plot/graphic/*.o
	-rm -f ./kine_subs/rfi_plot/*.o
