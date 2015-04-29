from lxml import etree
import os
import sys
from glob import glob
from shutil import move
import subprocess
import logging

tandem = os.environ.get('TANDEMEXE')
tandem2xml = os.environ.get('TANDEM2XML')

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
        tandem=tandem, tandem2xml=tandem2xml):
    if tandem is None:
        logging.error('TANDEM executable not set')
        return
    if convert and tandem2xml is None:
        logging.error('TANDEM2XML executable not set')
        return
    if isinstance(spectra, list):
        return [
            runtandem(folder, params, db, s, convert, overwrite)
            for s in spectra
            ]
    elif isinstance(spectra, str):
        logging.info("runtandem: searching files by mask %s...", spectra)
        splist = glob(spectra)
        logging.info("{} file(s) found.".format(len(splist)))
        if len(splist) > 1:
            return runtandem(folder, params, db, splist, convert, overwrite)
        elif splist:
            if os.path.exists(params):
                params = build_dict(params)
                params["spectrum, path"] = splist[0]
            else:
                logging.error('Could not load input file from %s.' % params)
                return
        else:
            logging.error('Nothing to process.')
            return
    else:
        logging.error("Could not find the input file %s", params)
        return

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
