from pyteomics import tandem
from sys import argv
import jinja2

parameters = dict()
params = tandem.iterfind(argv[1], 'group[type=parameters]', recursive=True)
for param in params:
    parameters[param['label']] = {v['label']: (v['note'] if 'note' in v else "") for v in param['note']}

templateLoader = jinja2.FileSystemLoader(searchpath="../templates/")
templateEnv = jinja2.Environment(loader=templateLoader)
TEMPLATE_FILE = "template.jinja"
template = templateEnv.get_template(TEMPLATE_FILE)
