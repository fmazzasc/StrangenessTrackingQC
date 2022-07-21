#!/bin/sh

import json
import os

# Reconstruction selections
sv_cuts = 'svertexer.minCosPAXYMeanVertex=-1;svertexer.minDCAToPV=0.001;svertexer.minCosPA=-1;svertexer.maxChi2=5;svertexer.maxDCAXYToMeanVertexV0Casc=10.;svertexer.maxDCAXYCasc=0.9;svertexer.maxDCAZCasc=0.9;svertexer.minRDiffV0Casc=0.1;svertexer.checkV0Hypothesis=false;svertexer.checkCascadeHypothesis=false;svertexer.maxPVContributors=3;svertexer.minCosPAXYMeanVertexCascV0=-1'
# ab_cuts = 'tpcitsMatch.requireToReachLayerAB=6;tpcitsMatch.minContributingLayersAB = 1; tpcitsMatch.maxABLinksOnLayer=100;tpcitsMatch.maxABFinalHyp=100;tpcitsMatch.cutABTrack2ClChi2=100;tpcitsMatch.nABSigmaY=10;tpcitsMatch.nABSigmaZ=10'
ab_cuts = 'tpcitsMatch.cutABTrack2ClChi2=50'

#read json

#read json

json_dict = json.load(open('workflow.json'))
json_items = json_dict['stages']
for item in json_items:
    if item['name'].startswith('svfinder'):
        cmd_list = item['cmd'].split(' ')
        for cmd_ind,cmd in enumerate(cmd_list):
            if cmd.startswith('--configKeyValues'):
                cmd_list[cmd_ind+1] = cmd_list[cmd_ind+1][:-1] + ';' + sv_cuts + cmd_list[cmd_ind+1][-1:]
                new_cmd = ' '.join(cmd_list)
                item['cmd'] = new_cmd
    
    if item['name'].startswith('itstpcMatch'):
        cmd_list = item['cmd'].split(' ')
        for cmd_ind,cmd in enumerate(cmd_list):
            if cmd.startswith('--configKeyValues'):
                cmd_list[cmd_ind+1] = cmd_list[cmd_ind+1][:-2] + ';' + ab_cuts + cmd_list[cmd_ind+1][-2:]
                new_cmd = ' '.join(cmd_list)
                item['cmd'] = new_cmd

json.dump(json_dict, open('workflow_mod.json','w'), indent=4)
                            