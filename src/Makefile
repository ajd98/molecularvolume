# Makefile for volume calculation tools

CC=gcc
IDIR=../include
CFLAGS=-I$(IDIR)

_DEPS=queue.h 
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))

OBJ = floodfill3d.o queue.o

%.o: %.c $(DEPS)
	$(CC) -c -O3 -o $@ $< $(CFLAGS)

volume: $(OBJ)
	ar cr libfloodfill3d.a $(OBJ)
	python setup.py build_ext --inplace
	mv *.so ../

install: $(volume)
