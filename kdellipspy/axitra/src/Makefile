################## PARAMETERS TO MODIFY *********************************

CC=gcc		#	C Compiler
CFLAGS=-O3 -std=gnu89	# 	and its options
#
# GFORTRAN compiler
#
 FC=gfortran
 F2PY=-lgomp --opt="-Ofast -march=native" --f90flags="-fopenmp -cpp"
 FFLAGS = -Ofast -cpp -fopenmp -fPIC

#
# Intel ifort compiler (>V15)
#
#FC=ifort
#F2PY=-L/opt/intel/compilers_and_libraries_2017.4.181/mac/compiler/lib -liomp5 --f90flags="-qopenmp -fpp"
# debug with openmp
#FFLAGS = -g -qopenmp -fpp -assume byterecl 
# debug without openmp
#FFLAGS = -g -fpp -assume byterecl 
#
# run without openmp
#FFLAGS = -O3 -fpp -assume byterecl 
# run with openmp
#FFLAGS = -O3 -fpp -qopenmp -assume byterecl 


################## END OF MODIFICATIONS *************************
################## Nothing to be changed below #################
.SUFFIXES: .o .f90
.SUFFIXES: .o .f

.c.o:
	$(CC) $(CFLAGS) -c $<

.f90.o:	$(INC)
	$(FC) $(FFLAGS) -c $<

all:	axitra convms

OBJ1=parameter.o dimension1.o dimension2.o allocateArray.o initdata.o reflect0.o reflect1.o reflect2.o reflect3.o \
reflect4.o reflect5.o ff0ad.o dll2km.o axitra.o
OBJ2=parameter.o dll2km.o dimension1.o dimension2.o  fft2cd.o fsource.o convm.o
#######################################
#			AXITRA
#######################################

axitra:	$(INC) $(OBJ1)
	$(FC) $(FFLAGS) -o $@ $(OBJ1) 

######################################
#			CONVM
######################################

convms: $(OBJ2) parameter.o writesac.o
	$(FC) $(FFLAGS) -o $@ writesac.o $(OBJ2)

######################################
#			CLEAN, TAR, ...
######################################
clean: 
	$(RM) *.o *.mod *.a *.so axitra convm convms 2>/dev/null

all: axitra convms
	
tar:	
	tar cvf axitra.tar * ;gzip axitra.tar\

######################################
#			PYTHON
######################################
PYTHON?=python
python: parameter.o dimension1.o fft2cd.o fsource.o dll2km.o axitra2py
	@echo "dict(real=dict(fd='double', fs='float'), complex=dict(fd='complex_double', fs='complex_float'))" > .f2py_f2cmap
	$(PYTHON) -m numpy.f2py convmPy.F90 parameter.f90 fft2cd.f90 fsource.f90 dll2km.f90 -m convmPy only: moment_conv :
	$(FC) $(FFLAGS) -c convmPy-f2pywrappers2.f90
	$(FC) $(FFLAGS) -c convmPy.F90
	$(CC) $(CFLAGS) -fPIC -c convmPymodule.c \
		-I`$(PYTHON) -c "import numpy; print(numpy.get_include())"` \
		-I`$(PYTHON) -c "import sysconfig; print(sysconfig.get_path('include'))"` \
		-I`$(PYTHON) -c "import numpy.f2py; print(numpy.f2py.get_include())"`
	$(CC) $(CFLAGS) -fPIC -c `$(PYTHON) -c "import numpy.f2py; print(numpy.f2py.get_include())"`/fortranobject.c \
		-I`$(PYTHON) -c "import numpy; print(numpy.get_include())"` \
		-I`$(PYTHON) -c "import sysconfig; print(sysconfig.get_path('include'))"` \
		-I`$(PYTHON) -c "import numpy.f2py; print(numpy.f2py.get_include())"`
	$(FC) $(FFLAGS) -shared -o convmPy`$(PYTHON) -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"` \
		convmPymodule.o fortranobject.o convmPy-f2pywrappers2.o convmPy.o \
		parameter.o dimension1.o fft2cd.o fsource.o dll2km.o

axitra2py: dimension1.o dimension2.o parameter.o allocateArray.o initdata.o dll2km.o axitra2py.o 
	$(FC) $(FFLAGS) -o $@ axitra2py.o dimension1.o dimension2.o allocateArray.o initdata.o dll2km.o parameter.o
