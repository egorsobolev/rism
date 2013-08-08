RISMROOT=/home/egor/c-proj/rism
CC=gcc
CFLAGS=-O3

SRC=$(wildcard *.c)
OBJ=$(SRC:.c=.o)

all:	straj

%.o:	%.c
	$(CC) $(CFLAGS) -I $(RISMROOT)/include -c -o $@ $<

straj:	$(OBJ) librism.a
	$(CC) -L$(RISMROOT)/lib $(OBJ) -lrism -largtable2 -lblas -lm\
	 -o $(RISMROOT)/bin/$@

librism.a:
	cd ../rism && make

clean:
	rm -f $(RISMROOT)/bin/straj
	rm -f $(OBJ)
