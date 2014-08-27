#! python
from pyteomics.pepxmltk import convert
from os import path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='*', help='list of pepXML or tXML files')
parser.add_argument('-o', required=True, help='path to output file')
parser.add_argument('-f', default=None, help='fdr, %%')
args =  parser.parse_args()

fdr = args.f
if fdr:
    fdr /= 100
convert(args.files, path.abspath(args.o), fdr=fdr)
