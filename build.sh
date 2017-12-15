#!/bin/sh
rm *.so
export LIBRARY_PATH=$(pwd):${LIBRARY_PATH}
gcc -O3 -c -fPIC floodfill3d.c -lm -o libfloodfill3d.so
python setup.py build_ext --inplace
python test.py
