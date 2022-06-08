from __future__ import absolute_import
from pyteomics import tandem, parser, mass, pepxml
from copy import copy
from collections import OrderedDict
import jinja2
import os
from lxml import etree
import argparse
from xml.sax import saxutils
import logging
import shutil
import sys


class Modifications:
    def __init__(self, input_parameters):
        self.modifications = []
        self.std_aa_mass = copy(mass.std_aa_mass)
        self.variable_mods = []
        self.std_aa_mass['N_term'] = float(input_parameters['protein, cleavage N-terminal mass change'])
        self.std_aa_mass['C_term'] = float(input_parameters['protein, cleavage C-terminal mass change'])
        self.get_modifications_from_params(input_parameters)
        self.info_about_xtandem_term_modifications()
        self.sort_modifications()
        self.change_std_aa_mass()
        self.calculate_modification_masses()
        self.add_lowercase_for_term_modifications()

    @staticmethod
    def get_modification_dict(mod, variable):
        modification = dict()
        modification['massdiff'], modification['aminoacid'] = mod.split('@')
        modification['massdiff'] = float(modification['massdiff'])
        modification['variable'] = variable
        if modification['aminoacid'] == '[':
            modification['terminus'] = 'N'
            modification['aminoacid'] = None
        elif modification['aminoacid'] == ']':
            modification['terminus'] = 'C'
            modification['aminoacid'] = None
        else:
            modification['terminus'] = None
        return modification

    def get_modifications_from_params(self, input_parameters):
        for mod in set(
                input_parameters.get('residue, potential modification mass', ',').split(',')
                    + (input_parameters.get('refine, potential modification mass', '').split(',')
                        if input_parameters.get('refine', 'no') == 'yes' else [])
                    + (['42.010565@['] if input_parameters.get('protein, quick acetyl', 'yes') == 'yes' else [])):
            if mod:
                self.modifications.append(self.get_modification_dict(mod, True))
        for mod in set(
                input_parameters.get('residue, modification mass', ',').split(',')
                + (input_parameters.get('refine, modification mass', '').split(',')
                   if input_parameters.get('refine', 'no') == 'yes' else [])):
            if mod:
                self.modifications.append(self.get_modification_dict(mod, False))

    @staticmethod
    def is_equal_modification(modification1, modification2):
        return (modification1['aminoacid'] == modification2.split('@')[1] and
                modification1['massdiff'] - float(modification2.split('@')[0]) < 1e-5)

    def sort_modifications(self):
        self.modifications.sort(key=lambda item: item['variable'])

    def change_std_aa_mass(self):
        for modification in self.modifications:
            if not modification['variable']:
                if modification['aminoacid']:
                    self.std_aa_mass[modification['aminoacid']] += modification['massdiff']
            else:
                break

    def calculate_modification_masses(self):
        for modification in self.modifications:
            if modification['aminoacid']:
                modification['mass'] = (self.std_aa_mass[modification['aminoacid']]
                        + (modification['massdiff'] if modification['variable'] else 0))
            else:
                modification['mass'] = (self.std_aa_mass['%s_term' % modification['terminus']]
                        + (modification['massdiff'] if modification['variable'] else 0))

    def info_about_xtandem_term_modifications(self):
        xtandem_nterm_default = ['-17.0265@C', '-18.0106@E', '-17.0265@Q']
        for mod in xtandem_nterm_default:
            if not any(self.is_equal_modification(modification, mod)
                    for modification in self.modifications):
                temp_mod = self.get_modification_dict(mod, True)
                self.modifications.append(temp_mod)

    def add_lowercase_for_term_modifications(self):
        for modification in self.modifications:
            if modification['terminus']:
                modification['terminus_lower'] = modification['terminus'].lower()


class Psm:
    def __init__(self, psm_tandem, proteases, mods):
        fims = psm_tandem['support']['fragment ion mass spectrum']
        try:
            self.spectrum = fims['note'].replace('\n', '')
        except KeyError:
            self.spectrum = fims['id']
        try:
            self.rt = float(psm_tandem['rt'])
        except (KeyError, TypeError, ValueError):
            self.rt = None
        self.hit_rank = 1
        self.start_scan = psm_tandem['support']['fragment ion mass spectrum']['id']
        self.end_scan = psm_tandem['support']['fragment ion mass spectrum']['id']
        self.precursor_neutral_mass = psm_tandem['mh'] - mass.nist_mass['H+'][0][0]
        self.assumed_charge = psm_tandem['z']
        self.sequence = psm_tandem['protein'][0]['peptide']['seq']
        self.peptide_prev_aa = psm_tandem['protein'][0]['peptide']['pre'][-1].replace('[', '-')
        self.peptide_next_aa = psm_tandem['protein'][0]['peptide']['post'][0].replace(']', '-')
        self.protein, self.protein_descr = self.get_protein_info(psm_tandem['protein'][0]['note'])
        self.num_tot_proteins = len(psm_tandem['protein'])
        self.num_matched_ions = sum(v for k, v in psm_tandem['protein'][0]['peptide'].items() if '_ions' in k)
        self.tot_num_ions = (len(self.sequence) - 1) * sum(
                1 for k in psm_tandem['protein'][0]['peptide'] if '_ions' in k
                ) * max(self.assumed_charge - 1, 1)
        self.calc_neutral_mass = (psm_tandem['protein'][0]['peptide']['mh'] - mass.nist_mass['H+'][0][0])
        self.massdiff = self.precursor_neutral_mass - self.calc_neutral_mass
        self.num_tol_term = self.calc_num_tol_term(proteases)
        self.num_missed_cleavages = psm_tandem['protein'][0][
                'peptide']['missed_cleavages']
        self.start = psm_tandem['protein'][0]['peptide']['start']
        self.alternative_proteins = []
        for prot in psm_tandem['protein'][1:]:
            alt_protein = {}
            alt_protein['dbname'], alt_protein['descr'] = self.get_protein_info(prot['note'])
            alt_protein['num_tol_term'] = self.calc_num_tol_term(proteases)
            alt_protein['peptide_prev_aa'] = prot['peptide']['pre'][-1].replace('[', '-')
            alt_protein['peptide_next_aa'] = prot['peptide']['post'][0].replace(']', '-')
            self.alternative_proteins.append(alt_protein)
        self.modifications = []
        self.mod_label_n = ''
        self.mod_label_c = ''
        for mod in psm_tandem['protein'][0]['peptide'].get('aa', []):
            temp_info = self.get_modification_info(mod, mods)
            if temp_info:
                self.modifications.append(temp_info)
        self.mod_label = ((' ' + self.mod_label_n) if self.mod_label_n else '') + (
                (' ' + self.mod_label_c) if self.mod_label_c else '')
        score_list = ['hyperscore',
                      'nextscore',
                      'b_score',
                      'y_score',
                      'c_score',
                      'z_score',
                      'a_score',
                      'x_score',
                      'expect',
                      'sumI']
        self.scores = OrderedDict.fromkeys(score_list, 0)
        for k in psm_tandem['protein'][0]['peptide']:
            if k in self.scores:
                self.scores[k] = psm_tandem['protein'][0]['peptide'][k]
        self.scores['sumI'] = psm_tandem['sumI']
        for k, v in self.scores.copy().items():
            self.scores[k.replace('_', '')] = self.scores.pop(k)

    def get_modification_info(self, modification, mods):
        position = modification['at'] - self.start + 1
        aa = self.sequence[position - 1]
        flag = 1
        if modification['at'] == 1 and abs(
                modification['modified'] - 42.0106) <= 0.001:
            self.mod_label_n = 'mod_nterm_mass="43.0184"'
        for m in mods.modifications:
            if abs(modification['modified'] - m['massdiff']) <= 0.001:
                if position == 1 and m['terminus'] == 'N':
                    self.mod_label_n = 'mod_nterm_mass="%s"' % (m['mass'])
                    flag = 0
                elif position == len(self.sequence) and m['terminus'] == 'C':
                    self.mod_label_c = ' mod_cterm_mass="%s"' % (m['mass'])
                    flag = 0
        if flag:
            return {'position': position, 'mass': mass.std_aa_mass[aa] + modification['modified']}
        else:
            return None

    def calc_num_tol_term(self, proteases):
        num_tol_term = 0
        if any(self.sequence[-1] in protease.cut and self.peptide_next_aa not in protease.no_cut
               for protease in proteases):
            num_tol_term += 1
        if any(self.sequence[0] not in protease.no_cut and self.peptide_prev_aa in protease.cut
               for protease in proteases):
            num_tol_term += 1
        return num_tol_term

    @staticmethod
    def get_protein_info(protein_label):
        if len(protein_label.split(' ', 1)) == 2:
            return protein_label.split(' ', 1)
        else:
            return protein_label, ''


class Protease:
    def __init__(self, cleavage_rule):
        self.cleavage_rule = cleavage_rule
        self.name = None
        self.cut = None
        self.no_cut = None
        self.sense = None
        self.get_info()

    def get_info(self):
        self.name = self.cleavage_rule  # TODO
        for protease in self.cleavage_rule.split(','):
            protease = protease.replace('X', ''.join(parser.std_amino_acids))
            c_term_rule, n_term_rule = protease.split('|')
            self.sense = self.get_sense(c_term_rule, n_term_rule)
            if self.sense == 'C':
                self.cut, self.no_cut = self.get_cut(c_term_rule, n_term_rule)
            else:
                self.cut, self.no_cut = self.get_cut(n_term_rule, c_term_rule)

    @staticmethod
    def get_sense(c_term_rule, n_term_rule):
        if '{' in c_term_rule:
            return 'N'
        elif '{' in n_term_rule:
            return 'C'
        else:
            if len(c_term_rule) <= len(n_term_rule):
                return 'C'
            else:
                return 'N'

    @staticmethod
    def get_cut(cut, no_cut):
        aminoacids = set(parser.std_amino_acids)
        cut = ''.join(aminoacids & set(cut))
        if '{' in no_cut:
            no_cut = ''.join(aminoacids & set(no_cut))
            return cut, no_cut
        else:
            no_cut = ''.join(set(parser.std_amino_acids) - set(no_cut))
            return cut, no_cut


def easy_write_pepxml(input_files, path_to_output, valid_psms=None):
    n = len(input_files)
    if path_to_output in input_files:
        with open(path_to_output) as f:
            tmp_lines = f.readlines()
    with open(path_to_output, 'w') as output_file:
        for i, infile in enumerate(input_files):
            unlocked = (i == 0)
            if infile != path_to_output:
                with open(infile) as input_file:
                    lines = input_file.readlines()
            else:
                lines = tmp_lines
            for line in lines:
                if '<spectrum_query' in line:
                    check_line = saxutils.unescape(line.split('spectrum="')[1].split('" ')[0], {'&quot;': '"'})
                    if valid_psms and check_line not in valid_psms:
                        unlocked = False
                    else:
                        unlocked = True
                if unlocked:
                    output_file.write(line)
                if '</spectrum_query>' in line:
                    unlocked = (i == n - 1)


def convert(files, path_to_output, fdr=None):
    path_to_file = os.path.abspath(files[0])
    path_to_output = os.path.abspath(path_to_output)
    parameters = dict()
    for _, elem in etree.iterparse(path_to_file, tag='group'):
        if elem.attrib['type'] == 'parameters':
            parameters[elem.attrib['label']] = OrderedDict(
                sorted((sub.attrib['label'], getattr(sub.text, 'strip', lambda: '')())
                    for sub in elem.iterchildren()))
        elem.clear()
    proteases = [Protease(rule) for rule in
            parameters['input parameters']['protein, cleavage site'].split(',')]
    modifications = Modifications(parameters['input parameters'])

    templatevars = {
        'parameters': parameters,
        'proteases': proteases,
        'path_to_file': path_to_file,
        'path_to_output': path_to_output,
        'modifications': [m for m in modifications.modifications if m['aminoacid']],
        'term_modifications': [m for m in modifications.modifications if not m['aminoacid']],
        'psms': (Psm(psm_tandem, proteases, modifications)
            for path_to_file in files for psm_tandem in tandem.read(path_to_file))
    }
    write(**templatevars)
    if fdr:
        merge_pepxml([path_to_output], path_to_output, fdr)


def merge_pepxml(input_files, path_to_output, fdr=None):
    if fdr:
        logging.info('Applying FDR %.1%.', fdr)
        psms = set()
        for infile in input_files:
            f = pepxml.filter(infile, fdr=fdr)
            psms.update(psm['spectrum'] for psm in f)
        easy_write_pepxml(input_files, path_to_output, psms)
    else:
        logging.info('Merging files into %s.', path_to_output)
        easy_write_pepxml(input_files, path_to_output, None)


def write(**template_vars):
    templateloader = jinja2.FileSystemLoader(os.path.join(os.path.dirname(__file__), 'templates'))
    templateenv = jinja2.Environment(loader=templateloader, autoescape=True)
    template_file = "template.jinja"
    template = templateenv.get_template(template_file)

    with open(template_vars['path_to_output'], 'w') as output:
        output.write(template.render(template_vars))


def main():
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
    files = parser.add_mutually_exclusive_group(required=True)
    files.add_argument('files', nargs='*', help='list of X!Tandem XML files', metavar='FILE', default=[])
    files.add_argument('--pepxml', nargs='+', help='list of pepXML files', metavar='FILE')
    parser.add_argument('-o', '--out', help='path to output file')
    parser.add_argument('-f', '--fdr', default=None, help='FDR, %% (only works with pepXML input files)')
    parser.add_argument('-v', '--verbosity', action='count', default=0, help='Increase output verbosity')
    args = parser.parse_args()
    levels = [logging.ERROR, logging.INFO, logging.DEBUG]
    level = 2 if args.verbosity > 2 else args.verbosity
    logging.basicConfig(format='%(levelname)5s: %(asctime)s %(message)s',
            datefmt='[%H:%M:%S]', level=levels[level])

    fdr = args.fdr
    if fdr:
        fdr = float(fdr) / 100

    if args.pepxml:
        logging.debug('Processing pepXML: %s', args.pepxml)
        if args.out:
            out = args.out
            infiles = args.pepxml
        else:
            out = args.pepxml[-1]
            infiles = args.pepxml[:-1]
        if not infiles:
            logging.info('Nothing to do.')
            sys.exit(0)
        elif len(infiles) == 1 and not fdr:
            logging.info('Copying %s to %s', infiles[0], out)
            shutil.copy(infiles[0], out)
        else:
            merge_pepxml(infiles, out, fdr)

    elif args.files:
        logging.debug('Processing files: %s', args.files)
        if args.out:
            convert(args.files, args.out)
        elif len(args.files) == 2:
            logging.info('Converting %s to %s.', *args.files)
            convert([args.files[0]], args.files[1], fdr)
        else:
            for f in args.files:
                out = os.path.splitext(f)[0] + ".pep.xml"
                logging.info('Converting %s to %s.', f, out)
                convert(f, out, fdr)
    logging.info('Done.')
