from distutils.core import setup
from Cython.Build import cythonize

setup(
    name = "volume",
    ext_modules = cythonize("volume.pyx")
)
