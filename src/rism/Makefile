include ../config.mk

SRC=$(wildcard *.c)
OBJMPI = $(addprefix .o.mpi/,$(SRC:.c=.o))
OBJSRL = $(addprefix .o.serial/,$(SRC:.c=.o))

all:    serial mpi

serial: librism.a

mpi:    librism-mpi.a

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

librism.a:	$(OBJSRL)
	ar -r $(RISMROOT)/lib/$@ $(OBJSRL)

librism-mpi.a:	$(OBJMPI)
	ar -r $(RISMROOT)/lib/$@ $(OBJMPI)

.PHONY:	clean
clean:
	-rm -r $(OBJMPI) $(OBJSRL) $(RISMROOT)/lib/librism.a $(RISMROOT)/lib/librism-mpi.a
	-rmdir .o.mpi .o.serial
