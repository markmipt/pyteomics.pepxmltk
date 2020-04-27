from lxml import etree
import os
import sys
import subprocess
import logging
import argparse
from . import pepxmltk

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


def get_free_name(txml):
    i = 1
    folder, initial = os.path.split(txml)
    base, ext = os.path.splitext(initial)
    if base.endswith('.t'):
        base, ext2 = os.path.splitext(base)
        ext = ext2 + ext
    while os.path.isfile(txml):
        txml = os.path.join(folder, '{}_{}{}'.format(base, i, ext))
        i += 1
    if i > 1:
        logging.debug("Changed file name to %s", txml)
    else:
        logging.debug("File name left unchanged: %s", txml)
    return txml


def runtandem(folder, params, db, spectra=None, convert=True, overwrite=False,
        tandem=_tandem, tandem2xml=_tandem2xml):
    if tandem is None:
        logging.error('TANDEM executable not set')
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
    if db:
        params["list path, taxonomy information"] = taxonomy_xml(db, os.path.join(folder, "taxonomy.xml"), "python")
        params["protein, taxon"] = "python"
    else:
        logging.info('Database not provided, using input parameters.')
    if convert:
        params["output, parameters"] = "yes"
        params["output, proteins"] = "yes"
        params["output, performance"] = "yes"
        params.setdefault("protein, cleavage N-terminal mass change", "1.00794")
        params.setdefault("protein, cleavage C-terminal mass change", "17.00305")

    namestub = os.path.split(spectra)[1].rsplit('.', 1)[0]
    txml_basename = namestub + '.t.xml'
    txml = os.path.join(folder, txml_basename)
    params["output, path hashing"] = "no"
    if not overwrite:
        txml = get_free_name(txml)
    params["output, path"] = txml
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

        if convert:
            logging.info("Converting results...")
            base, ext = os.path.splitext(txml)
            pepxml_file = base.rsplit('.t', 1)[0] + '.pep.xml'
            if tandem2xml:
                if not subprocess.call([tandem2xml, txml, pepxml_file], stderr=tandemout, stdout=tandemout):
                    logging.info("Task finished. The results are at {}.".format(pepxml_file))
                    return pepxml_file
                else:
                    logging.error("Conversion failed for file %s", spectra)
                    return txml
            else:
                logging.info('Using pepxmltk for result conversion.')
                pepxmltk.convert([txml], pepxml_file)
                logging.info("Task finished. The results are at {}.".format(pepxml_file))
                return pepxml_file
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
    params.pop("list path, default parameters", None)
    return params


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='XML file with search parameters', required=True)
    parser.add_argument('-o', '--dir', help='Directory to store the results. Defaults to current directory.')
    parser.add_argument('-db', '--fasta', help='FASTA database for search. If given, '
        'taxonomy file will be generated and added to parameters. If omitted, taxonomy must be correctly set '
        'in the input parameters.')
    parser.add_argument('spectra', nargs='*', help='Any number of data files to search')
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
        logging.error("X!Tandem executable not specified. Use --tandem.exe or set the TANDEMEXE variable")
        sys.exit(2)

    if args.fasta is not None and not os.path.exists(args.fasta):
        logging.error("Could not find the database: %s" % args.fasta)
        sys.exit(1)

    directory = args.dir or os.getcwd()

    runtandem(directory, args.input, args.fasta, spectra, args.convert, args.overwrite, tandem, tandem2xml)
