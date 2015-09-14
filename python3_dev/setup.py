from distutils.core import setup
from Cython.Build import cythonize
import numpy as np

#setup(
#    ext_modules=[
#        Extension("my_module", ["my_module.c"],
#                  include_dirs=[numpy.get_include()]),
#    ],
#)

setup(
    include_dirs = [np.get_include()],   
    ext_modules = cythonize("scale/*.pyx", annotate=True)
)
