from distutils.log import debug
import os
import sys
from joblib import Parallel, delayed


debug_level = "debug"
base_path = "/data/fmazzasc/its_data/sim/xi"

def run_strangeness_tracking(dire, debug_level='info'):
    os.chdir(dire)
    os.system(f'o2-strangeness-tracking-workflow -b --severity {debug_level}')

tf_paths = []
dirs = os.listdir(base_path)
for dire in dirs:
    if not dire.startswith('tf'):
        continue
    path = base_path + '/' + dire
    tf_paths.append(path)
    os.chdir(path)
    files_list = os.listdir(path)
    if os.path.islink('o2sim_geometry.root'):
        os.unlink('o2sim_geometry.root')
    if os.path.islink('o2sim_grp.root'):
        os.unlink('o2sim_grp.root')
    if os.path.islink('matbud.root'):
        os.unlink('matbud.root')
    geom_file = [f for f in files_list if (f.endswith('_geometry.root') and f.startswith('sgn'))][0]
    grp_file = [f for f in files_list if (f.endswith('_grp.root') and f.startswith('sgn'))][0]
    os.symlink(geom_file, 'o2sim_geometry.root')
    os.symlink(grp_file, 'o2sim_grp.root')
    os.symlink(base_path + "/matbud.root", 'matbud.root')

results = Parallel(n_jobs=len(tf_paths))(delayed(run_strangeness_tracking)(dire, debug_level) for dire in tf_paths)
