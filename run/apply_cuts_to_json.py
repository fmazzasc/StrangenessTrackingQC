#!/bin/sh

import json
import os

# Reconstruction selections
sv_cuts = 'svertexer.minCosPAXYMeanVertex=-1;svertexer.minDCAToPV=0.001;svertexer.minCosPA=-1;svertexer.maxChi2=5;svertexer.checkV0Hypothesis=false'

#read json

json_dict = json.load(open('workflow.json'))
json_items = json_dict['stages']
for item in json_items:
    if item['name'].startswith('svfinder'):
        cmd_list = item['cmd'].split(' ')
        for cmd_ind,cmd in enumerate(cmd_list):
            if cmd.startswith('--configKeyValues'):
                cmd_list[cmd_ind+1] = cmd_list[cmd_ind+1][:-1] + ';' + sv_cuts + '"'
                new_cmd = ' '.join(cmd_list)
                item['cmd'] = new_cmd

json.dump(json_dict, open('workflow_mod.json','w'), indent=4)                