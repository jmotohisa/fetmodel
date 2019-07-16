#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension

pycfet_sources = [
    'pycfet.i',
    'pycfet.c',
]

ext_pycfet = Extension(
    '_pycfet',
    sources=pycfet_sources,
    libraries=['gsl', 'gslcblas', 'cfet'],
    library_dirs=['/opt/local/lib', '/Users/motohisa/local/lib',
                  ],
    include_dirs=['/opt/local/include',
                  '/usr/local/include'],
)

if __name__ == '__main__':
    setup(
        name="pycfet",
        version="0.0.0",
        description="cyrindcial MESFET and MOSFET",
        ext_modules=[ext_pycfet],
        py_modules=['pycfet'],
    )
