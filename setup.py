from setuptools import setup
from Cython.Build import cythonize
import numpy

setup(
    name="aimsprop",
    ext_modules = cythonize("aimsprop/*.pyx", annotate=True),
    include_path = [numpy.get_include()],
    install_requires = ["numpy", "cython", "matplotlib"],
    scripts = ["scripts/need_for_speed_3.py"],
)
