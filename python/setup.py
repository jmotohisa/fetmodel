#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from distutils.core import setup, Extension
import numpy

SRC_PATH = os.path.relpath(os.path.join(os.path.dirname(__file__), "src"))


fetmodel_sources = [
    'fetmodel.i',
    'fetmodel.c',
]

ext_fetmodel = Extension(
    '_fetmodel',
    sources=fetmodel_sources,
    libraries=['gsl', 'gslcblas', 'fetmodel'],
    library_dirs=['/opt/local/lib', '/Users/motohisa/local/lib',
                  '../src/.libs',
                  ],
    include_dirs=['/opt/local/include',
                  '/usr/local/include'],
)


# class EggInfoCommand(egg_info):

#     def run(self):
#         if "build" in self.distribution.command_obj:
#             build_command = self.distribution.command_obj["build"]

#             self.egg_base = build_command.build_base

#             self.egg_info = os.path.join(
#                 self.egg_base, os.path.basename(self.egg_info))
#         egg_info.run(self)


if __name__ == '__main__':
    setup(
        name="fetmodel",
        version="0.0.0",
        description="FET models",
        ext_modules=[ext_fetmodel],
        include_dirs=[numpy.get_include()],
        py_modules=['fetmodel'],
        # package_dir={
        #     "": SRC_PATH, },
        # cmdclass={
        #     "egg_info": EggInfoCommand,
    )
