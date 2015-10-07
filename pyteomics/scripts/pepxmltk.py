#! python
from pyteomics.pepxmltk import convert
from os import path
import argparse

parser = argparse.ArgumentParser(
    description='Convert X!Tandem files to pepXML, or filter pepXML files.',
    epilog='''

Example usage
-------------

* Convert single file to pepXML:

  $ pepxmltk.py input.t.xml output.pep.xml

* Convert multiple files, derive names from input names:

  $ pepxmltk.py input1.t.xml input2.t.xml ... inputN.t.xml 

* Convert multiple files, join results into a single pepXML file:

  $ pepxmltk.py input1.t.xml input2.t.xml ... inputN.t.xml --out output.pep.xml
    
* Join multiple pepXML files and filter them to 1% FDR:

  $ pepxmltk.py input1.pep.xml input2.pep.xml ... inputN.pep.xml output.pep.xml --fdr 1

''',
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('files', nargs='+', help='list of pepXML or tXML files')
parser.add_argument('-o', '--out', help='path to output file')
parser.add_argument('-f', '--fdr', default=None, help='FDR, %% (only works with pepXML input files)')
args = parser.parse_args()

fdr = args.fdr
if fdr:
    fdr = float(fdr) / 100
if args.out:
    convert(args.files, path.abspath(args.out), fdr=fdr)
elif not any(path.splitext(path.splitext(filename)[0])[-1] == '.pep'
        for filename in args.files):
    for filename in args.files:
        convert((path.abspath(filename), ),
                path.abspath(filename).split('.t.xml')[0] + '.pep.xml', fdr=fdr)
else:
    convert(args.files[:-1], path.abspath(args.files[-1]), fdr=fdr)
