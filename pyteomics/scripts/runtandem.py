#! python
from __future__ import print_function
from lxml import etree
import sys
import os
from glob import glob
from shutil import move
import subprocess
import argparse

tandem = r'tandem.exe'
#tandem2xml = r'Tandem2XML'

def inputxml(path, params):
    with open(path, 'w') as inputxml:
        inputxml.write("<?xml version=\"1.0\"?>\n<bioml>\n")
        parameters = sorted(params)
        for par in parameters:
            inputxml.write("<note type=\"input\" label=\"%s\">%s</note>\n"
                  % (par, params[par]))
        inputxml.write("</bioml>")
    return path

def taxonomy_xml(database, path, taxon):
    with open(path, 'w') as xml:
        xml.write("<?xml version=\"1.0\"?>\n")
        xml.write("<bioml label=\"x! taxon-to-file matching list\">\n")
        xml.write("\t<taxon label=\"%s\">\n" % taxon)
        xml.write("\t\t<file format=\"peptide\" URL=\"%s\" />\n" % database)
        xml.write("\t</taxon>\n")
        xml.write("</bioml>")
    return path

def runtandem(folder, params, db, spectra=None, convert=True, overwrite=False):
    if isinstance(spectra, list):
        return [
            runtandem(folder, params, db, s, convert, overwrite)
            for s in spectra
            ]
    elif isinstance(spectra, str):
        print("runtandem: searching files by mask %s..." % spectra)
        splist = glob(spectra)
        print(len(splist), "file(s) found.")
        if len(splist) > 1:
            return runtandem(folder, params, db, splist, convert, overwrite)
        elif splist:
            params["spectrum, path"] = splist[0]

    if os.path.exists(params):
        params = build_dict(params)
    else:
        print("Could not find the input file %s" % params)
        return

    if not os.path.isdir(folder):
        print("Creating %s..." % folder)
        os.makedirs(folder)
    params["list path, taxonomy information"] = taxonomy_xml(
            db, os.path.join(folder, "taxonomy.xml"), "python")
    params["protein, taxon"] = "python"
    params["output, path"] = os.path.join(folder, "output.xml")
    inp = inputxml(os.path.join(folder, "input.xml"), params)
    print("Using the following parameters:")
    print(open(inp).read())
    print("Running X!Tandem...")
    subprocess.call([tandem, inp], stdout=sys.stderr)
    print("Collecting results...")
    outs = glob(os.path.join(folder, "output*t.xml"))
    out = [os.path.split(f)[1] for f in outs]
    out.sort()
    namestub = os.path.split(spectra)[1].rsplit('.', 1)[0]
    txml_name =  namestub + '.t.xml'
    txml = os.path.join(folder, txml_name)
    if not overwrite:
        i = 1
        while os.path.isfile(txml):
            txml = os.path.join(folder, '{}_{}.t.xml'.format(namestub, i))
            i += 1
    move(os.path.join(folder, out[-1]), txml)
    if convert:
        print("Converting results...")
        pepxml_file = txml.rsplit('.t.xml', 1)[0] + '.pep.xml'
        if not subprocess.call([tandem2xml, txml, pepxml_file]):
            print("Task finished. The results are at {}.".format(pepxml_file))
            return pepxml_file
        else:
            print("Conversion failed.")
            return txml
    else:
        print("Task finished. The results are at {}.".format(txml))
        return txml

def build_dict(xml):
    params = {}
    tree = etree.parse(xml)
    for note in tree.iter():
        if "label" in note.attrib:
            if note.text is None:
                params[note.attrib["label"]] = ''
            else:
                params[note.attrib["label"]] = note.text
    if "list path, default parameters" in params:
        del params["list path, default parameters"]
    return params

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='XML file with search parameters')
    parser.add_argument('dir', help='Directory to store the results')
    parser.add_argument('db', help='FASTA database for search')
    parser.add_argument('spectra', nargs='*',
            help='Any number of data files to search (globs supported)')
    parser.add_argument('--noconvert', dest='convert', action='store_false',
            help='Do not convert results to pepXML')
    parser.add_argument('--overwrite', action='store_true',
            help='Overwrite previous files instead of renaming the new ones')
    parser.add_argument('--tandem.exe', dest='tandem', metavar='FILE',
            help='X!Tandem executable')
    parser.add_argument('--tandem2xml', metavar='FILE',
            help='Tandem2XML converter executable')

    args = parser.parse_args()
    print("Starting runtandem...")

    if args.db is not None and not os.path.exists(args.db):
        print("Could not find the database %s" % args.db)
        sys.exit(1)

    spectra = args.spectra or None
    tandem = args.tandem or tandem
    tandem2xml = args.tandem2xml or tandem2xml
    runtandem(args.dir, args.input, args.db, spectra,
                args.convert, args.overwrite)
