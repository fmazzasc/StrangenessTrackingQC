from distutils.log import debug
import os
import json
from joblib import Parallel, delayed

os.environ['O2_INFOLOGGER_OPTIONS'] = 'floodProtection=0'
debug_level = "info"
base_path = "/data/fmazzasc/its_data/sim/MB_PbPB"
str_tracking_params = ''

def get_hbf_utils_from_json(json_file):
    out_list = []
    json_dict = json.load(open(json_file))
    json_items = json_dict['stages']
    for item in json_items:
        if item['name'].startswith('svfinder'):
            cmd_list = item['cmd'].split(' ')
            for cmd_ind,cmd in enumerate(cmd_list):
                if cmd.startswith('--configKeyValues'):

                    cfg = cmd_list[cmd_ind+1]
                    cfg_list = cfg.split(';')

                    for item_cfg in cfg_list:
                        if(item_cfg.startswith('"')):
                            item_cfg = item_cfg[1:]
                        if item_cfg.startswith('HBFUtils.'):
                            out_list.append(item_cfg)
                    return ";".join(out_list)


tf_paths = []
dirs = os.listdir(base_path)
config_key_values = str_tracking_params
config_key_values += ";" + get_hbf_utils_from_json(f"{base_path}/workflow.json")
for dire in dirs:
    if not dire.startswith('tf'):
        continue
    path = base_path + '/' + dire
    tf_paths.append(path)
    os.chdir(path)

def run_strangeness_tracking(dire, debug_level='info', config_key_values=""):
    os.chdir(dire)
    os.system(f'o2-strangeness-tracking-workflow -b --infologger-severity {debug_level}  --configKeyValues "{config_key_values}"')

results = Parallel(n_jobs=len(tf_paths))(delayed(run_strangeness_tracking)(dire, debug_level, config_key_values) for dire in tf_paths)