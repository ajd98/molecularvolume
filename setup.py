import numpy
from distutils.core import setup, Extension
from Cython.Build import cythonize
import os

recursivefilldir = os.path.abspath(os.getcwd())
print(recursivefilldir)

ext_modules=[
    Extension("volume",
        sources=["volume.pyx"],
        libraries=["m", "recursivefill"],
        library_dirs=[recursivefilldir],
        extra_compile_args=["-O3"]
    )
]


setup(
    name="volume",
    ext_modules=cythonize(ext_modules),
    include_dirs=[numpy.get_include()],
    library_dirs=[recursivefilldir]
)
