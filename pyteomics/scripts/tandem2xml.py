#! python
from pyteomics.tandem2xml import convert
from os import path
from sys import argv

convert(path.abspath(argv[1]), path.abspath(argv[2]))