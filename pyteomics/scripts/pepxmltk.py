#! python
from pyteomics.pepxmltk import convert
from os import path
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='*', help='list of pepXML or tXML files')
parser.add_argument('-o', help='path to output file')
parser.add_argument('-f', default=None, help='fdr, %%')
args = parser.parse_args()

fdr = args.f
if fdr:
    fdr = float(fdr) / 100
if args.o:
    convert(args.files, path.abspath(args.o), fdr=fdr)
elif not any(path.splitext(path.splitext(filename)[0])[-1] == '.pep'
        for filename in args.files):
    for filename in args.files:
        convert((path.abspath(filename), ),
                path.abspath(filename).split('.t.xml')[0] + '.pep.xml', fdr=fdr)
else:
    convert(args.files[:-1], path.abspath(args.files[-1]), fdr=fdr)
