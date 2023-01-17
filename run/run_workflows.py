from distutils.log import debug
import os
import json
from joblib import Parallel, delayed
os.environ['O2_INFOLOGGER_OPTIONS'] = 'floodProtection=0'


## configurable parameters
debug_level = "debug"
base_path = "/data/fmazzasc/its_data/sim/xi_new"
str_tracking_params = 'strtracker.mVertexMatching=true'
aod_producer_params = '--reco-mctracks-only 1 --aod-writer-keep dangling --aod-writer-resfile AO2D --aod-writer-resmode "UPDATE" --run-number 300000 \
                       --aod-timeframe-id ${ALIEN_PROC_ID}001 -b --run --condition-not-after 3385078236000 --shm-segment-size ${SHMSIZE:-50000000000} \
                       --info-sources ITS,MFT,MCH,TPC,ITS-TPC,MFT-MCH,ITS-TPC-TOF,TPC-TOF,FT0,FDD,TPC-TRD,ITS-TPC-TRD,ITS-TPC-TRD-TOF,CTP,FV0,EMC,MID \
                       --lpmp-prod-tag ${ALIEN_JDL_LPMPRODUCTIONTAG:-unknown} \
                       --anchor-pass ${ALIEN_JDL_LPMANCHORPASSNAME:-unknown} --anchor-prod ${ALIEN_JDL_MCANCHOR:-unknown} \
                       --combine-source-devices \
                       --enable-strangeness-tracking'



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
config_key_values_str = str_tracking_params
config_key_values_str += ";" + get_hbf_utils_from_json(f"{base_path}/workflow.json")
for dire in dirs:
    if not dire.startswith('tf2'):
        continue
    # if int(dire.split('tf')[1]) > 1:
    #     continue
    path = base_path + '/' + dire
    tf_paths.append(path)
    os.chdir(path)

def run_strangeness_tracking(dire, debug_level='info', config_key_values=""):
    os.chdir(dire)
    os.system(f'o2-strangeness-tracking-workflow -b --infologger-severity {debug_level}  --configKeyValues "{config_key_values_str}"')

def run_aod_producer(dire, params , debug_level='info'):
    os.chdir(dire)
    os.system(f'o2-aod-producer-workflow -b {params}  --infologger-severity {debug_level}')
 


results = Parallel(n_jobs=len(tf_paths))(delayed(run_strangeness_tracking)(dire, debug_level, config_key_values_str) for dire in tf_paths)
# results = Parallel(n_jobs=len(tf_paths))(delayed(run_aod_producer)(dire, aod_producer_params, debug_level) for dire in tf_paths)