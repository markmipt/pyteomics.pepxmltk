#!/usr/bin/env python

'''
setup.py file for pyteomics.pyteomics
'''

from distutils.core import setup

version = open('VERSION').readline().strip()

setup(
    name = 'pyteomics.pepxmltk',
    version = version,
    description      = '''A utility for creation of pepXML files from Python objects and TandemXML files.''',
    long_description = (''.join(open('README').readlines())),
    author           = 'Mark Ivanov & Lev Levitsky',
    author_email     = 'pyteomics@googlegroups.com',
    url              = 'http://hg.theorchromo.ru/pyteomics',
    requires         = [line.strip() for line in open('requirements.txt')],
    classifiers      = ['Intended Audience :: Science/Research',
                        'Programming Language :: Python :: 2.7',
                        'Programming Language :: Python :: 3',
                        'Topic :: Education',
                        'Topic :: Scientific/Engineering :: Bio-Informatics',
                        'Topic :: Scientific/Engineering :: Chemistry',
                        'Topic :: Scientific/Engineering :: Physics',
                        'Topic :: Software Development :: Libraries'],
    license          = 'License :: OSI Approved :: Apache Software License',
    packages = ['pyteomics'],
    package_dir = {'pyteomics': 'pyteomics'},
    package_data={'pyteomics': ['templates/template.jinja']},
    scripts=['pyteomics/scripts/pepxmltk.py', 'pyteomics/scripts/runtandem']
    )
