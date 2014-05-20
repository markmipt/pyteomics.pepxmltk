from pyteomics import tandem, parser, mass
import jinja2


class Psm:
    def __init__(self, psm_tandem):
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
        self.num_tol_term = self.calc_num_tol_term()
        self.num_missed_cleavages = psm_tandem['protein'][0]['peptide']['missed_cleavages']

        self.alternative_proteins = []
        for prot in psm_tandem['protein'][1:]:
            alt_protein = dict()
            alt_protein['dbname'], alt_protein['descr'] = self.get_protein_info(prot['label'])
            alt_protein['num_tol_term'] = self.calc_num_tol_term()
            self.alternative_proteins.append(alt_protein)

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

    @staticmethod
    def calc_num_tol_term():  # TODO
        return 2

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
        parameters[param['label']] = {v['label']: (v['note'] if 'note' in v else "") for v in param['note']}
    proteases = (Protease(rule) for rule in parameters['input parameters']['protein, cleavage site'].split(','))
    psms = (Psm(psm_tandem) for psm_tandem in tandem.read(path_to_file))
    templateloader = jinja2.FileSystemLoader(searchpath="../templates/")
    templateenv = jinja2.Environment(loader=templateloader)
    template_file = "template.jinja"
    template = templateenv.get_template(template_file)

    templatevars = {'parameters': parameters,
                    'proteases': proteases,
                    'path_to_file': path_to_file,
                    'modifications': [],
                    'psms': psms
    }
    output = open(path_to_output, 'w')
    output.write(template.render(templatevars))


if __name__ == '__main__':
    from sys import argv
    from os import path

    convert(path.abspath(argv[1]), path.abspath(argv[2]))