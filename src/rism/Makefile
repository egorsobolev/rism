include ../config.mk

SRC=$(wildcard *.c)
OBJMPI = $(addprefix .o.mpi/,$(SRC:.c=.o))
OBJSRL = $(addprefix .o.serial/,$(SRC:.c=.o))

all:    serial mpi

serial: librism.a

mpi:    librism-mpi.a

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

librism.a:	$(OBJSRL)
	ar -r ../../lib/$@ $(OBJSRL)

librism-mpi.a:	$(OBJMPI)
	ar -r ../../lib/$@ $(OBJMPI)

.PHONY:	clean
clean:
	-rm -rf .o.mpi .o.serial ../../lib/librism.a ../../lib/librism-mpi.a
