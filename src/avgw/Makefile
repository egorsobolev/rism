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
	$(CC) $(CFLAGS) $(MPICFLAGS) -I $(RISMROOT)/include -c -o $@ $<

$(OBJMPI):      | .o.mpi

$(OBJSRL):      | .o.serial

.o.mpi:
	mkdir .o.mpi

.o.serial:
	mkdir .o.serial

mkavgw:  $(OBJSRL) librism.a
	$(CC) $(LFLAGS) -L$(RISMROOT)/lib $(OBJSRL) -lrism -largtable2 \
         -lfftw3f -lfftw3 -lblas -lm -o $(RISMROOT)/bin/$@

mkavgw.MPI:      $(OBJMPI) librism-mpi.a
	$(CC) $(LFLAGS) $(MPILFLAGS) -L$(RISMROOT)/lib $(OBJMPI) \
         -lrism-mpi -largtable2 -lfftw3f -lfftw3 -lblas -lm \
         -o $(RISMROOT)/bin/$@

librism.a:
	cd ../rism && make serial

librism-mpi.a:
	cd ../rism && make mpi

.PHONY: clean
clean:
	-rm -r $(OBJMPI) $(OBJSRL) $(RISMROOT)/bin/mkavgw $(RISMROOT)/bin/mkavgw.MPI
	-rmdir .o.mpi .o.serial
