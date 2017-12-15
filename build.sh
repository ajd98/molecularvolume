#!/bin/sh
export LIBRARY_PATH=$(pwd):${LIBRARY_PATH}
gcc -c -fPIC recursivefill.c -lm -o librecursivefill.so
python setup.py build_ext --inplace
python test.py
