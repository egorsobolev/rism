RISMROOT=/home/egor/c-proj/rism
CC=gcc
CFLAGS=-O3

SRC=$(wildcard *.c)
OBJ=$(SRC:.c=.o)

%.o:	%.c
	$(CC) $(CFLAGS) -I $(RISMROOT)/include -c -o $@ $<

librism.a:	$(OBJ)
	ar -r $(RISMROOT)/lib/$@ $(OBJ)

all:    librism.a
