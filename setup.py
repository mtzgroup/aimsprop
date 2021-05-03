from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    name="aimsprop",
    ext_modules = cythonize("aimsprop/*.pyx", annotate=True),
    include_path = [numpy.get_include()]
)
