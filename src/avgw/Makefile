RISMROOT=/home/egor/c-proj/rism

CC=gcc
CFLAGS=-O3

#CC=icc
#CFLAGS=-O3 -I/usr/include/x86_64-linux-gnu
#LFLAGS=-L/usr/lib/x86_64-linux-gnu

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

clean:
	rm -f $(RISMROOT)/bin/mkavgw
	rm -f $(OBJ)
