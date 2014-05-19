from pyteomics import tandem
from pyteomics import parser
import jinja2


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

        def get_cut(cut, no_cut):
            aminoacids = set(parser.std_amino_acids)
            cut = ''.join(aminoacids & set(cut))
            if '{' in no_cut:
                no_cut = ''.join(aminoacids & set(no_cut))
                return cut, no_cut
            else:
                no_cut = ''.join(set(parser.std_amino_acids) - set(no_cut))
                return cut, no_cut

        for protease in self.cleavage_rule.split(','):
            protease = protease.replace('X', ''.join(parser.std_amino_acids))
            c_term_rule, n_term_rule = protease.split('|')
            self.sense = get_sense(c_term_rule, n_term_rule)
            if self.sense == 'C':
                self.cut, self.no_cut = get_cut(c_term_rule, n_term_rule)
            else:
                self.cut, self.no_cut = get_cut(n_term_rule, c_term_rule)


def convert(path_to_file, path_to_output):
    parameters = dict()
    params = tandem.iterfind(argv[1], 'group[type=parameters]', recursive=True)
    for param in params:
        parameters[param['label']] = {v['label']: (v['note'] if 'note' in v else "") for v in param['note']}
    proteases = (Protease(rule) for rule in parameters['input parameters']['protein, cleavage site'].split(','))
    templateloader = jinja2.FileSystemLoader(searchpath="../templates/")
    templateenv = jinja2.Environment(loader=templateloader)
    template_file = "template.jinja"
    template = templateenv.get_template(template_file)

    templatevars = {'parameters': parameters,
                    'proteases': proteases,
                    'path_to_file': path_to_file,
                    'modifications': [],
                    'psms': []
                    }
    output = open(path_to_output, 'w')
    output.write(template.render(templatevars))

if __name__ == '__main__':
    from sys import argv
    from os import path
    convert(path.abspath(argv[1], argv[2]))