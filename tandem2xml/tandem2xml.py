from pyteomics import tandem, parser, mass
from copy import copy
from collections import OrderedDict
import jinja2


class Modifications():
    def __init__(self, path_to_file, input_parameters):
        self.modifications = []
        self.std_aa_mass = copy(mass.std_aa_mass)
        self.fixed_mods = []
        self.variable_mods = []
        self.std_aa_mass['N_term'] = float(input_parameters['protein, cleavage N-terminal mass change'])
        self.std_aa_mass['C_term'] = float(input_parameters['protein, cleavage C-terminal mass change'])
        for psm_tandem in tandem.read(path_to_file):
            for mod in psm_tandem['protein'][0]['peptide'].get('aa', []):
                modification = dict()
                modification['aminoacid'] = mod['type']
                modification['massdiff'] = mod['modified']

                if mod['at'] == psm_tandem['protein'][0]['peptide']['start']:
                    modification['terminus'] = 'N'
                elif mod['at'] == psm_tandem['protein'][0]['peptide']['end']:
                    modification['terminus'] = 'C'
                else:
                    modification['terminus'] = None

                modification['variable'] = 'Y'
                if modification not in self.modifications:
                    self.modifications.append(modification)
        self.get_modifications_from_params(input_parameters)
        self.variable_info()
        self.group_terminus()
        self.remove_dupcated_modifications()
        self.drop_wrong_terminus()
        self.sort_modifications()
        self.change_std_aa_mass()
        self.calculate_modification_masses()
        self.info_about_xtandem_term_modifications()
        self.add_lowercase_for_term_modifications()

    def get_modifications_from_params(self, input_parameters):
        for mod in input_parameters['residue, potential modification mass'].split(','):
            self.variable_mods.append((mod.split('@')))
        for mod in input_parameters['residue, modification mass'].split(','):
            self.fixed_mods.append((mod.split('@')))
        if input_parameters.get('protein, quick acetyl', 'yes') == 'yes':
            self.variable_mods.append(['42.010565', '['])
        if input_parameters['refine'] == 'yes':
            for mod in input_parameters['refine, potential modification mass'].split(',')\
                    + input_parameters['refine, potential N-terminus modifications'].split(',')\
                    + input_parameters['refine, potential C-terminus modifications'].split(','):
                if (mod.split('@')) not in self.variable_mods:
                    self.variable_mods.append((mod.split('@')))
            for mod in input_parameters['refine, modification mass'].split(','):
                if (mod.split('@')) not in self.fixed_mods:
                    self.fixed_mods.append((mod.split('@')))

    @staticmethod
    def is_equal_modification(modification1, modification2):
        if isinstance(modification2, list):
            if round(modification1['massdiff'] - float(modification2[0]), 4) == 0:
                if (not modification1['terminus'] and modification1['aminoacid'] == modification2[1]) \
                    or (modification1['terminus'] == 'N' and modification2[1] == '[') \
                        or (modification1['terminus'] == 'C' and modification2[1] == ']'):
                            return True
        elif isinstance(modification2, dict):
            if modification1['aminoacid'] == modification2['aminoacid'] \
                and modification1['massdiff'] == modification2['massdiff']:
                return True

    def sort_modifications(self):
        self.modifications = sorted(self.modifications, key=lambda item: item['variable'])

    def change_std_aa_mass(self):
        for modification in self.modifications:
            if modification['variable'] == 'N':
                if modification['aminoacid']:
                    self.std_aa_mass[modification['aminoacid']] += modification['massdiff']
            else:
                break

    def calculate_modification_masses(self):
        for modification in self.modifications:
            if modification['aminoacid']:
                modification['mass'] = round(self.std_aa_mass[modification['aminoacid']]\
                                       + (modification['massdiff'] if modification['variable'] == 'Y' else 0), 5)
            else:
                modification['mass'] = round(self.std_aa_mass['%s_term' % (modification['terminus'], )]\
                                       + (modification['massdiff'] if modification['variable'] == 'Y' else 0), 5)

    def info_about_xtandem_term_modifications(self):
        for modification in self.modifications:
            if modification['terminus'] and modification['aminoacid']:
                modification['close_label'] = ' symbol="^" /><!--X! Tandem %s-terminal AA variable modification-->'\
                                        % (modification['terminus'].lower())
            else:
                modification['close_label'] = ' />'

    def remove_dupcated_modifications(self):
        self.modifications = [dict(t) for t in set([tuple(d.items()) for d in self.modifications])]

    def drop_wrong_terminus(self):
        for modification_copy in list(self.modifications):
            if not modification_copy['terminus']:
                for modification in list(self.modifications):
                    if modification['terminus'] and self.is_equal_modification(modification_copy, modification):
                        self.modifications.remove(modification)

    def variable_info(self):
        for modification in self.modifications:
            for mod in self.fixed_mods:
                if self.is_equal_modification(modification, mod):
                    modification['variable'] = 'N'
                    break

    def group_terminus(self):
        for modification in self.modifications:
            if modification['terminus'] and any(self.is_equal_modification(modification, mod)
                                                for mod in self.fixed_mods + self.variable_mods):
                modification['aminoacid'] = None

    def add_lowercase_for_term_modifications(self):
        for modification in self.modifications:
            if modification['terminus']:
                modification['terminus_lower'] = modification['terminus'].lower()


class Psm:
    def __init__(self, psm_tandem, proteases):
        self.spectrum = psm_tandem['support']['fragment ion mass spectrum']['note']
        self.start_scan = psm_tandem['support']['fragment ion mass spectrum']['id']
        self.end_scan = psm_tandem['support']['fragment ion mass spectrum']['id']
        self.precursor_neutral_mass = psm_tandem['protein'][0]['peptide']['mh'] - mass.calculate_mass('H+')
        self.assumed_charge = psm_tandem['z']
        self.sequence = psm_tandem['protein'][0]['peptide']['seq']
        self.peptide_prev_aa = psm_tandem['protein'][0]['peptide']['pre'][-1]
        self.peptide_next_aa = psm_tandem['protein'][0]['peptide']['post'][0]
        self.protein, self.protein_descr = self.get_protein_info(psm_tandem['protein'][0]['label'])
        self.num_tot_proteins = len(psm_tandem['protein'])
        self.num_matched_ions = sum(v for k, v in psm_tandem['protein'][0]['peptide'].iteritems() if '_ions' in k)
        self.tot_num_ions = (len(self.sequence) - 1) * \
                            sum(1 for k in psm_tandem['protein'][0]['peptide'] if '_ions' in k)
        self.calc_neutral_mass = mass.calculate_mass(self.sequence) \
                                 + sum([mod['modified'] for mod in psm_tandem['protein'][0]['peptide'].get('aa', [])])
        self.massdiff = round(self.precursor_neutral_mass - self.calc_neutral_mass, 6)
        self.num_tol_term = self.calc_num_tol_term(proteases)
        self.num_missed_cleavages = psm_tandem['protein'][0]['peptide']['missed_cleavages']
        self.start = psm_tandem['protein'][0]['peptide']['start']
        self.alternative_proteins = []
        for prot in psm_tandem['protein'][1:]:
            alt_protein = dict()
            alt_protein['dbname'], alt_protein['descr'] = self.get_protein_info(prot['label'])
            alt_protein['num_tol_term'] = self.calc_num_tol_term(proteases)
            self.alternative_proteins.append(alt_protein)
        self.modifications = []
        for mod in psm_tandem['protein'][0]['peptide'].get('aa', []):
            self.modifications.append(self.get_modification_info(mod))

        score_list = ['hyperscore',
                      'nextscore',
                      'bscore',
                      'yscore',
                      'cscore',
                      'zscore',
                      'ascore',
                      'xscore',
                      'expect']
        self.scores = dict.fromkeys(score_list, 0)
        for k in psm_tandem['protein'][0]['peptide']:
            if k in self.scores:
                self.scores[k] = psm_tandem['protein'][0]['peptide'][k]

    def get_modification_info(self, modification):
        position = modification['at'] - self.start + 1
        aa = self.sequence[position - 1]
        return {'position': position, 'mass': mass.std_aa_mass[aa] + modification['modified']}

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
        return protein_label.split(' ', 1)


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


def convert(path_to_file, path_to_output):
    parameters = dict()
    params = tandem.iterfind(path_to_file, 'group[type="parameters"]', recursive=True)
    for param in params:
        parameters[param['label']] = OrderedDict(sorted({v['label']: (v['note'] if 'note' in v else "")
                                                         for v in param['note']}.items(), key=lambda (k, v): k))
    proteases = [Protease(rule) for rule in parameters['input parameters']['protein, cleavage site'].split(',')]
    modifications = Modifications(path_to_file, parameters['input parameters'])
    psms = (Psm(psm_tandem, proteases) for psm_tandem in tandem.read(path_to_file))
    templateloader = jinja2.FileSystemLoader(searchpath="../templates/")
    templateenv = jinja2.Environment(loader=templateloader)
    template_file = "template.jinja"
    template = templateenv.get_template(template_file)

    templatevars = {'parameters': parameters,
                    'proteases': proteases,
                    'path_to_file': path_to_file,
                    'path_to_output': path_to_output,
                    'modifications': [m for m in modifications.modifications if m['aminoacid']],
                    'term_modifications': [m for m in modifications.modifications if not m['aminoacid']],
                    'psms': psms
    }
    output = open(path_to_output, 'w')
    output.write(template.render(templatevars))


if __name__ == '__main__':
    from sys import argv
    from os import path

    convert(path.abspath(argv[1]), path.abspath(argv[2]))