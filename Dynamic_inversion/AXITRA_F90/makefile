
# sur tarapaca
#FF = ifort  -i4 -assume byterecl -auto  -FR  -O3 -ipo -pad -opt-report-levelmin -xP -vec-report0 
# test dgf
##FF = ifort -O3 -pad -ipo -xP -parallel
FF = ifort -xHOST -pad -ipo -xP -parallel 
# AUTRES
#FF= /usr/local/lf9560/bin/lf95 -g --chk a,e,s,u,x --chkglobal --trap
#FF= /usr/local/lf9560/bin/lf95 --o1 --ntrap
#FF= /usr/local/intel/fc/9.1/bin/ifort -O3 -pad -ipo -xB -assume byterecl -vec-report0 -opt-report-levelmin
#FF = /usr/FortranIntel/fce/9.1.037/bin/ifort -i4 -assume byterecl -auto  -FR  -mcmodel=large -O3 -ipo -pad -opt-report-levelmin -xP -vec-report0 
#FF= /usr/local/intel/fc/9.1/bin/ifort -warn all 
#FF= /usr/mpich-1.2.7p1-ifort//bin/mpif90    -i4 -assume byterecl -auto  -FR  -mcmodel=large      -O3 -ipo -pad -opt-report-levelmin   -xP -vec-report0 

axitra:
	$(FF)$ -c axitra_module.f90
	$(FF)$ -c axitra.f90
	$(FF)$ -c initdata.f90
	$(FF)$ -c reflect0.f90
	$(FF)$ -c reflect1.f90
	$(FF)$ -c reflect2.f90
	$(FF)$ -c reflect3.f90
	$(FF)$ -c reflect4.f90
	$(FF)$ -c reflect5.f90
	$(FF)$ -c ff01ad.f90
	$(FF)$ -c ff02ad.f90
	$(FF)$ -o axitra.out axitra_module.o axitra.o initdata.o reflect0.o reflect1.o reflect2.o reflect3.o reflect4.o reflect5.o ff01ad.o ff02ad.o 
	/bin/rm -f *.o *.mod > /dev/null

convm:
	$(FF) -c convm.f90
	$(FF) -c fsource.f90
	$(FF) -c fft2cd.f90
	$(FF) -o convm.out convm.o fsource.o fft2cd.o 
	/bin/rm -f *.o *.mod > /dev/null
sismo:
	$(FF) -c sismo.f90
	$(FF) -o sismo.out sismo.o
	/bin/rm -f *.o *.mod > /dev/null
