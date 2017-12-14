#!/bin/sh
python setup.py build_ext --inplace
python test.py
