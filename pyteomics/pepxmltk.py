from pyteomics import tandem, parser, mass, pepxml
from copy import copy
from collections import OrderedDict
import jinja2
from os import path
from lxml import etree
import argparse


class Modifications:
    def __init__(self, input_parameters):
        self.modifications = []
        self.std_aa_mass = copy(mass.std_aa_mass)
        self.variable_mods = []
        self.std_aa_mass['N_term'] = float(
                input_parameters['protein, cleavage N-terminal mass change'])
        self.std_aa_mass['C_term'] = float(
                input_parameters['protein, cleavage C-terminal mass change'])
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
        for mod in set(input_parameters.get(
            'residue, potential modification mass', ',').split(',')
                + (input_parameters.get(
                    'refine, potential modification mass', '').split(',')
                   if input_parameters.get('refine', 'no') == 'yes' else [])
                + (['42.010565@['] if input_parameters.get(
                    'protein, quick acetyl', 'yes') == 'yes' else [])):
            if mod:
                self.modifications.append(self.get_modification_dict(mod, True))
        for mod in set(input_parameters.get(
            'residue, modification mass', ',').split(',')
                + (input_parameters.get(
                    'refine, modification mass', '').split(',')
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
                    self.std_aa_mass[modification['aminoacid']
                            ] += modification['massdiff']
            else:
                break

    def calculate_modification_masses(self):
        for modification in self.modifications:
            if modification['aminoacid']:
                modification['mass'] = (
                        self.std_aa_mass[modification['aminoacid']]
                        + (modification['massdiff']
                            if modification['variable'] else 0))
            else:
                modification['mass'] = (
                        self.std_aa_mass['%s_term' % modification['terminus']]
                        + (modification['massdiff']
                            if modification['variable'] else 0))

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
        self.spectrum = psm_tandem['support']['fragment ion mass spectrum'
                ]['note'].replace('\n', '')
        try:
            self.rt = float(psm_tandem['rt'])
        except (KeyError, TypeError):
            self.rt = None
        self.hit_rank = 1
        self.start_scan = psm_tandem['support']['fragment ion mass spectrum']['id']
        self.end_scan = psm_tandem['support']['fragment ion mass spectrum']['id']
        self.precursor_neutral_mass = psm_tandem['mh'] - mass.calculate_mass('H+')
        self.assumed_charge = psm_tandem['z']
        self.sequence = psm_tandem['protein'][0]['peptide']['seq']
        self.peptide_prev_aa = psm_tandem['protein'][0]['peptide']['pre'][-1]
        self.peptide_prev_aa.replace('[', '-')
        self.peptide_next_aa = psm_tandem['protein'][0]['peptide']['post'][0]
        self.peptide_next_aa.replace(']', '-')
        self.protein, self.protein_descr = self.get_protein_info(
                psm_tandem['protein'][0]['note'])
        self.num_tot_proteins = len(psm_tandem['protein'])
        self.num_matched_ions = sum(
                v for k, v in psm_tandem['protein'][0]['peptide'].items()
                if '_ions' in k)
        self.tot_num_ions = (len(self.sequence) - 1) * sum(
                1 for k in psm_tandem['protein'][0]['peptide'] if '_ions' in k
                ) * max(self.assumed_charge - 1, 1)
        self.calc_neutral_mass = (psm_tandem['protein'][0]['peptide']['mh'] -
                mass.calculate_mass('H+'))
        self.massdiff = self.precursor_neutral_mass - self.calc_neutral_mass
        self.num_tol_term = self.calc_num_tol_term(proteases)
        self.num_missed_cleavages = psm_tandem['protein'][0][
                'peptide']['missed_cleavages']
        self.start = psm_tandem['protein'][0]['peptide']['start']
        self.alternative_proteins = []
        for prot in psm_tandem['protein'][1:]:
            alt_protein = {}
            alt_protein['dbname'], alt_protein['descr'] = self.get_protein_info(
                    prot['note'])
            alt_protein['num_tol_term'] = self.calc_num_tol_term(proteases)
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
            return {'position': position,
                    'mass': mass.std_aa_mass[aa] + modification['modified']}
        else:
            return None

    def calc_num_tol_term(self, proteases):
        num_tol_term = 0
        if any(self.sequence[-1] in protease.cut and
                self.peptide_next_aa not in protease.no_cut
               for protease in proteases):
            num_tol_term += 1
        if any(self.sequence[0] not in protease.no_cut and
                self.peptide_prev_aa in protease.cut
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
    unlocked = True
    if path_to_output in input_files:
        tmp_lines = open(path_to_output, 'r').readlines()
    output_file = open(path_to_output, 'w')
    for infile in input_files:
        if infile != path_to_output:
            input_file = open(infile, 'r')
            lines = input_file.readlines()
            input_file.close()
        else:
            lines = tmp_lines
        for line in lines:
            if '<spectrum_query' in line:
                if valid_psms and line.split(
                        'spectrum="')[1].split('" ')[0] not in valid_psms:
                    unlocked = False
                else:
                    unlocked = True
            if unlocked:
                output_file.write(line)
            if '</spectrum_query>' in line:
                unlocked = True
        unlocked = False
    output_file.close()


def convert(files, path_to_output, fdr=None):
    if path.splitext(path.splitext(files[0])[0])[-1] == '.t':
        for idx, path_to_file in enumerate(files):
            if idx == 0:
                path_to_file = path.abspath(path_to_file)
                path_to_output = path.abspath(path_to_output)
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

        templatevars = {'parameters': parameters,
                        'proteases': proteases,
                        'path_to_file': path_to_file,
                        'path_to_output': path_to_output,
                        'modifications': [m for m in modifications.modifications
                            if m['aminoacid']],
                        'term_modifications': [m for m in
                            modifications.modifications if not m['aminoacid']],
                        'psms': (Psm(psm_tandem, proteases, modifications)
                            for path_to_file in files for psm_tandem in tandem.read(path_to_file))
        }
        write(**templatevars)
        input_files = (path_to_output, )
    else:
        input_files = files
    if fdr:
        psms = set()
        for infile in input_files:
            f = pepxml.filter(infile, fdr=fdr)
            psms.update(psm['spectrum'] for psm in f)
        easy_write_pepxml(input_files, path_to_output, psms)
    elif len(input_files) > 1:
        easy_write_pepxml(input_files, path_to_output, None)

def write(**template_vars):
    templateloader = jinja2.PackageLoader('pyteomics.pepxmltk')
    templateenv = jinja2.Environment(loader=templateloader, autoescape=True,
            extensions=['jinja2.ext.autoescape'])
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
