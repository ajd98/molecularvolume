# Makefile for volume calculation tools

CC=gcc
IDIR=../include
CFLAGS=-I$(IDIR) -Wall -L. -g -std=gnu99 -fPIC

_DEPS=queue.h 
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))

OBJ = floodfill3d.o queue.o

volume.so: libfloodfill3d.a 
	python setup.py build_ext --inplace
	mv *.so ../

%.o: %.c $(DEPS)
	$(CC) -c -O3 -o $@ $< $(CFLAGS)

libfloodfill3d.a: $(OBJ)
	ar cr libfloodfill3d.a $(OBJ)

install: volume.so

clean:
	-rm test *.so *.a *.o

test: libfloodfill3d.a
	$(CC) test.c -O3 -o test $(CFLAGS) -lfloodfill3d -lm

.PHONY: install clean 
