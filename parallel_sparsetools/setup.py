from setuptools import setup, Extension
import numpy as np
import os

setup(
    name="Parallel Sparse matrix calculations",
    version="0.0.1",
    author="Martin Claus",
    packages=['parallel_sparsetools'],
    install_requires=['numpy'],
    ext_modules=[
        Extension('_parallel_sparsetools',
                  sources=[os.path.join('parallel_sparsetools', p) for p in ("parallel_csr.cxx", "parallel_sparsetools.cxx")],
                  include_dirs=['parallel_sparsetools', np.get_include()],
                  extra_compile_args=['-fopenmp', '-g', '-O3'],
                  extra_link_args=['-fopenmp'],
                 )
    ],
)
