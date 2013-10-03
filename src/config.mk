RISMROOT=/home/egor/c-proj/rism

CC=gcc
CFLAGS=-O3

#CC=icc
#CFLAGS=-O3 -I/usr/include/x86_64-linux-gnu
#LFLAGS=-L/usr/lib/x86_64-linux-gnu

# IMPB openmpi
MPICFLAGS = -I/usr/lib/openmpi/include -DMPI
MPILFLAGS = -L/usr/lib/openmpi/lib -lmpi -lopen-rte -lopen-pal -ldl \
    -Wl,--export-dynamic -lnsl -lutil -ldl

