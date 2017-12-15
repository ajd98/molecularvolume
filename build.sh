#!/bin/sh
rm *.o *.a *.so
export LIBRARY_PATH=$(pwd):${LIBRARY_PATH}

gcc -O3 -c floodfill3d.c -lm -o floodfill3d.o
gcc -O3 -c queue.c -o queue.o
ar cr libfloodfill3d.a queue.o floodfill3d.o

python setup.py build_ext --inplace
python test.py
