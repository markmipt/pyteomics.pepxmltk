#! python
from pyteomics.pepxmltk import convert
from os import path
from sys import argv

try:
    fdr = float(argv[-1]) / 100
    out_file = argv[-2]
    files = [path.abspath(f) for f in argv[1:-2]]
except:
    fdr = None
    out_file = argv[-1]
    files = [path.abspath(f) for f in argv[1:-1]]
convert(files, path.abspath(out_file), fdr)