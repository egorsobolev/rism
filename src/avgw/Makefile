include ../config.mk

SRC=$(wildcard *.c)
OBJMPI = $(addprefix .o.mpi/,$(SRC:.c=.o))
OBJSRL = $(addprefix .o.serial/,$(SRC:.c=.o))

all:    serial mpi

serial: mkavgw

mpi:    mkavgw.MPI

.o.serial/%.o:  %.c
	$(CC) $(CFLAGS) -I $(RISMROOT)/include -c -o $@ $<

.o.mpi/%.o:     %.c
	$(MPICC) -DMPI $(CFLAGS) $(MPICFLAGS) -I $(RISMROOT)/include -c -o $@ $<

$(OBJMPI):      | .o.mpi

$(OBJSRL):      | .o.serial

.o.mpi:
	mkdir .o.mpi

.o.serial:
	mkdir .o.serial

mkavgw:	librism.a ../../lib/libcpocketfft.a $(OBJSRL)
	$(CC) $(LFLAGS) -L../../lib $(OBJSRL) -lrism -largtable2 \
         -lblas -lm -lcpocketfft -lstdc++ \
         -o $(RISMROOT)/bin/$@

mkavgw.MPI:	librism-mpi.a ../../lib/libcpocketfft.a $(OBJMPI)
	$(MPICC) $(LFLAGS) $(MPILFLAGS) -L../../lib $(OBJMPI) \
         -lrism-mpi -largtable2 -lblas -lm -lcpocketfft -lstdc++ \
         -o $(RISMROOT)/bin/$@

librism.a:
	cd ../rism && make serial

librism-mpi.a:
	cd ../rism && make mpi

../../lib/libcpocketfft.a:
	{ \
	cd ../../cpocketfft/src ; \
	make ; \
	cp libcpocketfft.a ../../lib ; \
	cp fft.h ../../include ; \
	}

.PHONY: clean
clean:
	-rm -r $(OBJMPI) $(OBJSRL) $(RISMROOT)/bin/mkavgw $(RISMROOT)/bin/mkavgw.MPI
	-rmdir .o.mpi .o.serial
	-rm -rf ../../lib/libcpocketfft.a ../../include/fft.h
