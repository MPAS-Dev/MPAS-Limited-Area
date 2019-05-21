#!/usr/bin/env python3

from distutils.core import setup

setup(name='MPAS-Limited-Area',
       version = '0.1',
       description = 'Python application for creating limited area MPAS meshes',
       packages = ['LimitedArea'],
       scripts = ['LimitedArea/bin/limited-area'],
      )
