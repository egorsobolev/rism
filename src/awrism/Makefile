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

awrism:  $(OBJSRL) librism.a
	$(CC) $(LFLAGS) -L../../lib $(OBJSRL) -lrism -largtable2 \
         -lfftw3f -lfftw3 -lblas -lm -o ../../bin/$@

awrism.MPI:      $(OBJMPI) librism-mpi.a
	$(MPICC) $(LFLAGS) $(MPILFLAGS) -L../../lib $(OBJMPI) \
         -lrism-mpi -largtable2 -lfftw3f -lfftw3 -lblas -lm \
         -o ../../bin/$@

librism.a:
	cd ../rism && make serial

librism-mpi.a:
	cd ../rism && make mpi

.PHONY: clean
clean:
	-rm -rf .o.mpi .o.serial ../../bin/awrism ../../bin/awrism.MPI
