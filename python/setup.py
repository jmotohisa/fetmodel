#!/usr/bin/env python
# -*- coding: utf-8 -*-

from distutils.core import setup, Extension

cfet_sources = [
    'cfet.i',
    'cfef.c',
    '../ccm.c',
    '../ctl-io.h',
]

ext_cfet = Extension(
    '_cfet',
    sources=cfet_sources,
    libraries=['gsl', 'gslcblas'],
    library_dirs=['/opt/local/lib', ],
    include_dirs=['/opt/local/include', '/usr/local/include'],
)

if __name__ == '__main__':
    setup(
        name="cfet",
        version="0.0.0",
        description="cyrindcial MESFET and MOSFET",
        ext_modules=[ext_cfet],
        py_modules=['cfet'],
    )
