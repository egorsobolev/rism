RISMROOT=/home/egor/c-proj/rism

CC=gcc
CFLAGS=-O3

#CC=icc
#CFLAGS=-O3 -I/usr/include/x86_64-linux-gnu
#LFLAGS=-L/usr/lib/x86_64-linux-gnu

all:	librism.a straj dhist

clean:  clean_librism.a clean_straj clean_dhist

librism.a:
	cd rism && make

straj:  librism.a
	cd straj && make

dhist:  librism.a
	cd dhist && make

clean_librism.a:
	cd rism && make clean

clean_straj:
	cd straj && make clean

clean_dhist:
	cd dhist && make clean
