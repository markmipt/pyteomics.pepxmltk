from lxml import etree
import os
import sys
from glob import glob
from shutil import move
import subprocess
import logging
import argparse

_tandem = os.environ.get('TANDEMEXE')
_tandem2xml = os.environ.get('TANDEM2XML')

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

def runtandem(folder, params, db, spectra=None, convert=True, overwrite=False,
        tandem=_tandem, tandem2xml=_tandem2xml):
    if tandem is None:
        logging.error('TANDEM executable not set')
        return
    if convert and tandem2xml is None:
        logging.error('TANDEM2XML executable not set')
        return

    if isinstance(params, str):
        if os.path.exists(params):
            params = build_dict(params)
        else:
            logging.error('Input file %s does not exist.' % params)
            return

    if isinstance(spectra, list):
        return [
            runtandem(folder, params, db, s, convert, overwrite, tandem, tandem2xml)
            for s in spectra
            ]
    else:
        if spectra is None:
            spectra = params.get("spectrum, path")
            if spectra is None:
                logging.error('No spectra to process.')
                return
            else:
                logging.info('No "spectra" argument specified, using input file.')
        logging.info("runtandem: processing file %s...", spectra)
        params["spectrum, path"] = spectra
        
    if not os.path.isdir(folder):
        logging.info("Creating %s..." % folder)
        os.makedirs(folder)
    params["list path, taxonomy information"] = taxonomy_xml(
            db, os.path.join(folder, "taxonomy.xml"), "python")
    params["protein, taxon"] = "python"
    params["output, path"] = os.path.join(folder, "output.xml")
    inp = inputxml(os.path.join(folder, "input.xml"), params)
    logging.debug("Using the following parameters:")
    with open(inp) as fi:
        for line in fi:
            logging.debug(line.rstrip())
    logging.info("Running X!Tandem...")
    with open(os.devnull, 'w') as devnull:
        tandemout = (sys.stderr if
            logging.getLogger().getEffectiveLevel() < logging.INFO
            else devnull)
        code = subprocess.call([tandem, inp], stdout=tandemout)
        if code:
            logging.error('Search failed for file %s (exit code %s)', spectra, code)
            return code
        logging.info("Collecting results...")
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
            logging.info("Converting results...")
            pepxml_file = txml.rsplit('.t.xml', 1)[0] + '.pep.xml'
            if not subprocess.call([tandem2xml, txml, pepxml_file],
                    stderr=tandemout, stdout=tandemout):
                logging.info("Task finished. The results are at {}.".format(pepxml_file))
                return pepxml_file
            else:
                logging.error("Conversion failed for file %s", spectra)
                return txml
        else:
            logging.info("Task finished. The results are at {}.".format(txml))
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

def main():
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
    parser.add_argument('-v', '--verbosity', action='count', default=0,
            help='Increase output verbosity')

    levels = [logging.ERROR, logging.INFO, logging.DEBUG]
    args = parser.parse_args()
    level = 2 if args.verbosity > 2 else args.verbosity
    logging.basicConfig(format='%(levelname)5s: %(asctime)s %(message)s',
            datefmt='[%H:%M:%S]', level=levels[level])

    logging.info("Starting runtandem...")
    spectra = args.spectra or None
    tandem = args.tandem or _tandem
    tandem2xml = args.tandem2xml or _tandem2xml
    if tandem is None:
        logging.error("X!Tandem executable not specified. "
                "Use --tandem.exe or set the TANDEMEXE variable")
        sys.exit(2)
    if args.convert and tandem2xml is None:
        logging.error("Tandem2XML command or executable not specified. "
                "Use --tandem2xml or set the TANDEM2XML variable")
        sys.exit(2)
    if args.db is not None and not os.path.exists(args.db):
        logging.error("Could not find the database: %s" % args.db)
        sys.exit(1)

    runtandem(args.dir, args.input, args.db, spectra,
                args.convert, args.overwrite, tandem, tandem2xml)
