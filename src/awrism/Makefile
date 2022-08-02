include ../config.mk

SRC=$(wildcard *.c)
OBJMPI = $(addprefix .o.mpi/,$(SRC:.c=.o))
OBJSRL = $(addprefix .o.serial/,$(SRC:.c=.o))

all:    serial mpi

serial: awrism

mpi:    awrism.MPI

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

awrism:  $(OBJSRL) librism.a
	$(CC) $(LFLAGS) -L$(RISMROOT)/lib $(OBJSRL) -lrism -largtable2 \
         -lfftw3f -lfftw3 -lblas -lm -o $(RISMROOT)/bin/$@

awrism.MPI:      $(OBJMPI) librism-mpi.a
	$(MPICC) $(LFLAGS) $(MPILFLAGS) -L$(RISMROOT)/lib $(OBJMPI) \
         -lrism-mpi -largtable2 -lfftw3f -lfftw3 -lblas -lm \
         -o $(RISMROOT)/bin/$@

librism.a:
	cd ../rism && make serial

librism-mpi.a:
	cd ../rism && make mpi

.PHONY: clean
clean:
	-rm -r $(OBJMPI) $(OBJSRL) $(RISMROOT)/bin/awrism $(RISMROOT)/bin/awrism.MPI
	-rmdir .o.mpi .o.serial
