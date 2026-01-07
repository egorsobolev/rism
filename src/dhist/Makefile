include ../config.mk

SRC=$(wildcard *.c)
OBJMPI = $(addprefix .o.mpi/,$(SRC:.c=.o))
OBJSRL = $(addprefix .o.serial/,$(SRC:.c=.o))

all:    serial mpi

serial: mkhist

mpi:    mkhist.MPI

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

mkhist:  $(OBJSRL) librism.a
	$(CC) $(LFLAGS) -L../../lib $(OBJSRL) -lrism -largtable2 \
		-lm -o ../../bin/$@

mkhist.MPI:      $(OBJMPI) librism-mpi.a
	$(MPICC) $(LFLAGS) $(MPILFLAGS) -L../../lib $(OBJMPI) \
		-lrism-mpi -largtable2 -lm \
		-o ../../bin/$@

librism.a:
	cd ../rism && make serial

librism-mpi.a:
	cd ../rism && make mpi

.PHONY: clean
clean:
	-rm -rf .o.mpi .o.serial ../../bin/mkhist ../../bin/mkhist.MPI
