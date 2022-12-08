#!/bin/bash

#
# A example workflow MC->RECO->AOD for a simple pp min bias
# production, targetting test beam conditions.

# make sure O2DPG + O2 is loaded
export O2DPG_ROOT=/home/fmazzasc/alice/O2DPG

[ ! "${O2DPG_ROOT}" ] && echo "Error: This needs O2DPG loaded" && exit 1
[ ! "${O2_ROOT}" ] && echo "Error: This needs O2 loaded" && exit 1

# ----------- LOAD UTILITY FUNCTIONS --------------------------
. ${O2_ROOT}/share/scripts/jobutils.sh

# ----------- START ACTUAL JOB  -----------------------------

NWORKERS=${NWORKERS:-64}
MODULES="--skipModules ZDC"
SIMENGINE=${SIMENGINE:-TGeant4}
NSIGEVENTS=${NSIGEVENTS:-600}
NTIMEFRAMES=${NTIMEFRAMES:-20}
INTRATE=${INTRATE:-50000}
SYSTEM=${SYSTEM:-pp}
ENERGY=${ENERGY:-13500}
[[ ${SPLITID} != "" ]] && SEED="-seed ${SPLITID}" || SEED=""

# create workflow
${O2DPG_ROOT}/MC/bin/o2dpg_sim_workflow.py -eCM ${ENERGY} -col ${SYSTEM} -gen external -j ${NWORKERS} -ns ${NSIGEVENTS} -tf ${NTIMEFRAMES} -confKey "Diamond.width[2]=6." -e ${SIMENGINE} ${SEED} -mod "--skipModules ZDC" \
        -ini ${O2DPG_ROOT}/MC/config/PWGLF/ini/GeneratorLFXi${SYSTEM}.ini -field -5

# run workflow
# python3 apply_cuts_to_json.py
${O2DPG_ROOT}/MC/bin/o2_dpg_workflow_runner.py -f workflow_mod.json -tt aod --cpu-limit 64
