#!/bin/sh
export LIBRARY_PATH=$(pwd):${LIBRARY_PATH}
gcc -O3 -c -fPIC recursivefill.c -lm -o librecursivefill.so
python setup.py build_ext --inplace
python test.py
