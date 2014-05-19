from pyteomics import tandem
import jinja2


def convert(path_to_file):
    parameters = dict()
    params = tandem.iterfind(argv[1], 'group[type=parameters]', recursive=True)
    for param in params:
        parameters[param['label']] = {v['label']: (v['note'] if 'note' in v else "") for v in param['note']}

    templateloader = jinja2.FileSystemLoader(searchpath="../templates/")
    templateenv = jinja2.Environment(loader=templateloader)
    template_file = "template.jinja"
    template = templateenv.get_template(template_file)

if __name__ == '__main__':
    from sys import argv
    convert(argv[1])