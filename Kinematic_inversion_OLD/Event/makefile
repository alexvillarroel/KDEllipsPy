FFLAGS= 
#FF= f77 -fast -O4 -autopar -xarch=v8plusa -fsimple=2 -xvector -xprefetch
#FF= gfortran -frecord-marker=4 -O4
FF=ifx -fast
#LIBSAC= sacio.a
#LIBSAC= /usr/local/SAC2000/lib/sacio.a
#LIBSAC= /net/quake/local/lib/libsac.a
#LIBSAC= sac.a

all: axitra

axitra08:axitra08.o initdata.o reflect0.o reflect1.o reflect2.o reflect3.o reflect4.o reflect5.o ff0ad.o
	$(FF) -O axitra08.o initdata.o reflect0.o reflect1.o reflect2.o reflect3.o reflect4.o reflect5.o ff0ad.o -o axitra08

axitra:axitra.o initdata.o reflect0.o reflect1.o reflect2.o reflect3.o reflect4.o reflect5.o ff0ad.o
	$(FF) -O axitra.o initdata.o reflect0.o reflect1.o reflect2.o reflect3.o reflect4.o reflect5.o ff0ad.o -o axitra

convm:convm.o fsource.o ll2km.o fft2cd.o XAPIIR.o
	$(FF) -O convm.o fsource.o ll2km.o fft2cd.o XAPIIR.o -o convm 

faille: faille.o
	$(FF) -O3 faille.o -o faille
clean:
	rm -f axitra axitra08
	rm -f *.o

SUFFIXES: .f .o
.f.o:
	$(FF)  -c  ${FFLAGS}  $*.f

