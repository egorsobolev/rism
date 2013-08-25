RISMROOT=/home/egor/c-proj/rism

CC=gcc
CFLAGS=-O3

#CC=icc
#CFLAGS=-O3 -I/usr/include/x86_64-linux-gnu
#LFLAGS=-L/usr/lib/x86_64-linux-gnu

all:	librism.a straj mkhist mkavgw

clean:  clean_librism.a clean_straj clean_mkhist clean_mkavgw

librism.a:
	cd rism && make

straj:  librism.a
	cd straj && make

mkhist:  librism.a
	cd dhist && make

mkavgw:	librism.a
	cd avgw && make

clean_librism.a:
	cd rism && make clean

clean_straj:
	cd straj && make clean

clean_mkhist:
	cd dhist && make clean

clean_mkavgw:
	cd avgw && make clean
