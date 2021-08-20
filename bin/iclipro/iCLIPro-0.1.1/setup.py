#!/usr/bin/env python

from distutils.core import setup

try:
    import pysam
except ImportError:
    sys.stderr.write("Setup script for iCLIPro failed to import 'pysam'.\n")
    sys.stderr.write("Please install pysam and then try again.\n")
    sys.exit(1)

try:
    import matplotlib
except ImportError:
    sys.stderr.write("Setup script for iCLIPro failed to import 'matplotlib'.\n")
    sys.stderr.write("Please install matplotlib and then try again.\n")
    sys.exit(1)

setup(
    name='iCLIPro',
    version = file("iCLIPro/_version.py").readline().split('=')[1].strip().strip('"'),
    description = 'iCLIP read overlap test - ' 
        'indentification and correction of '
        'potential systematic errors in iCLIP data.',
    author = 'Tomaz Curk',
    author_email = 'tomaz.curk@fri.uni-lj.si',
    url = 'http://www.biolab.si/iCLIPro',

    classifiers = [
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)',
        'Operating System :: POSIX',
        'Programming Language :: Python'
    ],

    requires = ['pysam (>= 0.6)', 'matplotlib'],

    py_modules = [
        'iCLIPro._version',
        'iCLIPro.bed',
        'iCLIPro.cmd_utils',
    ],

    scripts = [
        'scripts/iCLIPro',
        'scripts/iCLIPro_bam_splitter',
    ],
)

