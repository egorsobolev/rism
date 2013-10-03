include ../config.mk

SRC=$(wildcard *.c)
OBJ=$(SRC:.c=.o)

all:	mkavgw

%.o:	%.c
	$(CC) $(CFLAGS) -I $(RISMROOT)/include -c -o $@ $<

mkavgw:	$(OBJ) librism.a
	$(CC) $(LFLAGS) -L$(RISMROOT)/lib $(OBJ) -lrism -largtable2 \
	 -lfftw3f -lfftw3 -lblas -lm -o $(RISMROOT)/bin/$@

librism.a:
	cd ../rism && make

.PHONY:	clean
clean:
	-rm -f $(RISMROOT)/bin/mkavgw
	-rm -f $(OBJ)
