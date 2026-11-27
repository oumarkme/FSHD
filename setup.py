"""Builds the Cython extensions; all other metadata is in pyproject.toml."""

from Cython.Build import cythonize
import numpy
from setuptools import Extension, setup

EXTENSIONS = [
    Extension(
        f"fshdetector.{name}",
        [f"src/fshdetector/{name}.pyx"],
        include_dirs=[numpy.get_include()],
        define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
        extra_compile_args=["-O3"],
    )
    for name in ("_screen", "_score")
]

setup(
    ext_modules=cythonize(
        EXTENSIONS,
        compiler_directives={
            "language_level": 3,
            "boundscheck": False,
            "wraparound": False,
            "cdivision": True,
        },
    )
)
