# Makefile for volume calculation tools

CC=gcc
IDIR=../include
CFLAGS=-I$(IDIR) -Wall -L. -g

_DEPS=queue.h 
DEPS=$(patsubst %,$(IDIR)/%,$(_DEPS))

OBJ = floodfill3d.o queue.o

%.o: %.c $(DEPS)
	$(CC) -c -O3 -o $@ $< $(CFLAGS)

libfloodfill3d.a: $(OBJ)
	ar cr libfloodfill3d.a $(OBJ)

volume: libfloodfill3d.a 
	python setup.py build_ext --inplace
	mv *.so ../

test: libfloodfill3d.a
	$(CC) test.c -O3 -o test $(CFLAGS) -lfloodfill3d -lm

install: volume

uninstall:
	if [ -a test ]; then\
      rm test ;\
	fi;
	rm *.so *.a *.o
