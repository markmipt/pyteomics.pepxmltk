#!/usr/bin/env python

'''
setup.py file for pyteomics.pepxmltk
'''

from setuptools import setup, find_packages

version = open('VERSION').readline().strip()

setup(
    name                 = 'pyteomics.pepxmltk',
    version              = version,
    description          = '''A utility for creation of pepXML files from Python objects and TandemXML files.''',
    long_description     = (''.join(open('README.md').readlines())),
    author               = 'Mark Ivanov & Lev Levitsky',
    author_email         = 'pyteomics@googlegroups.com',
    url                  = 'http://hg.theorchromo.ru/pyteomics',
    install_requires     = ['pyteomics[XML]', 'jinja2'],
    classifiers          = ['Intended Audience :: Science/Research',
                            'Programming Language :: Python :: 2.7',
                            'Programming Language :: Python :: 3',
                            'Topic :: Education',
                            'Topic :: Scientific/Engineering :: Bio-Informatics',
                            'Topic :: Scientific/Engineering :: Chemistry',
                            'Topic :: Scientific/Engineering :: Physics',
                            'Topic :: Software Development :: Libraries'],
    license              = 'License :: OSI Approved :: Apache Software License',
    packages             = find_packages(),
    namespace_packages   = ['pyteomics'],
    package_data         = {'pyteomics': ['templates/*']},
    entry_points         = {'console_scripts': ['pepxmltk.py = pyteomics.pepxmltk:main',
                                                'runtandem = pyteomics.runtandem:main']}
    )
