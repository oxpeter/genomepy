#!/usr/bin/python

""" A module to build the file paths of all database files needed for the other modules.
It first looks for the config file pathways.config, and if it cannot find that, then it
builds one itself.

Later improvements should include storing the object as a pickled object for faster
access.
"""

import os

def read_pathways(config_file):
    pathway_dict = {}
    config_h = open(config_file, 'rb')
    config_g = ( (line.split()[0], line.split()[1]) for line in config_h )
    for fileid, pathway in config_g:
        pathway_dict[fileid] = pathway
    return pathway_dict

def construct_pathways(config_file):
    pathway_dict = {}
    # find the directory containing the OGS files:
    tmpfile = os.path.dirname(gffparser.__file__) + '/dbsearch.tmp'

    cmd1 = 'find /home/antqueen/genomics/genomes -type d -iname *C.bir* -print > ' + tmpfile
    cmd2 = 'find /Volumes/antqueen/genomics/genomes -type d -iname *C.bir* -print >> ' + tmpfile
    os.system()
    return pathway_dict

def pick_latest(ogs_dir):
    return latest_file

if __name__ == "__main__":
    config_file = os.path.dirname(gffparser.__file__) + "/pathways.config"
    if os.path.exists(config_file):
        pathway_dict = read_pathways(config_file)
    else:
        pathway_dict = construct_pathways(config_file)