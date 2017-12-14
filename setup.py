import numpy
from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(
    ext_modules=cythonize("volume.pyx"),
    include_dirs=[numpy.get_include()]
)
