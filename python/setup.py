#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension
import numpy

pyfet_sources = [
    'pyfet.i',
    'pyfet.c',
]

ext_pyfet = Extension(
    '_pyfet',
    sources=pyfet_sources,
    libraries=['gsl', 'gslcblas', 'fetmodel'],
    library_dirs=['/opt/local/lib', '/Users/motohisa/local/lib',
                  '../src/.libs',
                  ],
    include_dirs=['/opt/local/include',
                  '/usr/local/include'],
)

if __name__ == '__main__':
    setup(
        name="pyfet",
        version="0.0.0",
        description="FET models",
        ext_modules=[ext_pyfet],
        include_dirs=[numpy.get_include()],
        py_modules=['pyfet'],
    )
