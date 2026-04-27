#
#	Make discription file for Receiver function example subroutines
#

#FFLAGS = -O4
#F77 = /opt/SUNWspro/bin/f77

DEMO = ../

include $(DEMO)/makena.macros   # Edit this file to change location or
                        # type of compilers (e.g. to gcc, g77)

DIR = .
RFI = $(DIR)/fd3d_theo.o  $(DIR)/deriv.o $(DIR)/fsource.o \
      $(DIR)/fft2cd.o $(DIR)/four1.o $(DIR)/XAPIIR.o $(DIR)/ll2km.o \
      $(DIR)/mkpeak.o $(DIR)/mkstress.o $(DIR)/convm_dynsource.o\
      $(DIR)/forward_modelling.o $(DIR)/calcmisfit.o $(DIR)/display_model.o \
      $(DIR)/display_final.o $(DIR)/lofw.o  \
      $(DIR)/output_summary.o  $(DIR)/mkpeak.o $(DIR)/mkstress.o
      

# make fd3d_subs 

fd3d_subs: $(RFI)

.f.o:
	$(F77) $(FFLAGS) -c $*.f
.c.o:
	$(CC) $(CFLAGS) -c $*.c
