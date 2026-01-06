include ../config.mk

SRC=$(wildcard *.c)
OBJMPI = $(addprefix .o.mpi/,$(SRC:.c=.o))
OBJSRL = $(addprefix .o.serial/,$(SRC:.c=.o))

all:    serial mpi

serial: awrism

mpi:    awrism.MPI

.o.serial/%.o:  %.c
	$(CC) $(CFLAGS) -I ../../include -c -o $@ $<

.o.mpi/%.o:     %.c
	$(MPICC) -DMPI $(CFLAGS) $(MPICFLAGS) -I ../../include -c -o $@ $<

$(OBJMPI):      | .o.mpi

$(OBJSRL):      | .o.serial

.o.mpi:
	mkdir .o.mpi

.o.serial:
	mkdir .o.serial

awrism:  librism.a ../../lib/libcpocketfft.a $(OBJSRL)
	$(CC) $(LFLAGS) -L../../lib $(OBJSRL) -lrism -largtable2 \
         -lblas -lm -lcpocketfft -lstdc++ \
         -o ../../bin/$@

awrism.MPI:      librism-mpi.a ../../lib/libcpocketfft.a $(OBJMPI)
	$(MPICC) $(LFLAGS) $(MPILFLAGS) -L../../lib $(OBJMPI) \
         -lrism-mpi -largtable2 -lblas -lm -lcpocketfft -lstdc++ \
         -o ../../bin/$@

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
	-rm -rf .o.mpi .o.serial ../../bin/awrism ../../bin/awrism.MPI
	-rm -rf ../../lib/libcpocketfft.a ../../include/fft.h
