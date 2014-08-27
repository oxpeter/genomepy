#!/usr/bin/python

""" A module to build the file paths of all database files needed for the other modules.
It first looks for the config file pathways.config, and if it cannot find that, then it
builds one itself.

Later improvements should include storing the object as a pickled object for faster
access.
"""

import os
import re
from subprocess import Popen, PIPE
from operator import itemgetter


def read_pathways(config_file):
    pathway_dict = {}
    config_h = open(config_file, 'rb')
    config_g = ( (line.split()[0], line.split()[1]) for line in config_h )
    for fileid, pathway in config_g:
        pathway_dict[fileid] = pathway
    return pathway_dict

def write_config(pathway_dict, config_file):
    config_h = open(config_file, 'w')
    for fileid in pathway_dict:
        config_h.write(" ".join([fileid,pathway_dict[fileid]]) + '\n')
    config_h.close()

def find_files():
    "creates a list of all potentially needed database file paths"
    pathway_dict = {}
    # find the directories containing the OGS and ortholog files:
    cmd1 = [ 'find', os.sep.join(['','home','antqueen','genomics','genomes']), '-type', 'd', '-iname', '*C.bir*', '-print' ]
    cmd2 = [ 'find', os.sep.join(['','Volumes','antqueen','genomics','genomes']), '-type', 'd', '-iname', '*C.bir*', '-print' ]
    cmd3 = [ 'find', os.sep.join(['','home','antqueen','genomics','experiments','analyses','BGI20120208_Genome']), '-type', 'd', '-iname', '*_o*', '-print' ]
    cmd4 = [ 'find', os.sep.join(['','Volumes','antqueen','genomics','experiments','analyses','BGI20120208_Genome']), '-type', 'd', '-iname', '*_o*', '-print' ]
    # collect output of commands:
    out1_h = Popen(cmd1, stdout=PIPE)
    out2_h = Popen(cmd2, stdout=PIPE)
    out3_h = Popen(cmd3, stdout=PIPE)
    out4_h = Popen(cmd4, stdout=PIPE)
    # parse output into a list:
    out1 = re.findall("([^\\n]*)\\n", out1_h.communicate()[0])
    out2 = re.findall("([^\\n]*)\\n", out2_h.communicate()[0])
    out3 = re.findall("([^\\n]*)\\n", out3_h.communicate()[0])
    out4 = re.findall("([^\\n]*)\\n", out4_h.communicate()[0])
    filingcabinet = []
    for path in out1 + out2 + out3 + out4: # check each result
        if os.path.exists(path):    # removes any blank find results
            ls_h = Popen([ 'ls', path ], stdout=PIPE)
            # create a list of all visible filenames in the directory
            ls_out = re.findall("([^\\n]*)\\n", ls_h.communicate()[0])
            for filename in ls_out:
                filingcabinet.append(os.sep.join([path, filename]))

    return filingcabinet

def find_latest(filelist, seqtype, fextension):
    "finds the most recent version of a given file in filelist"
    shortlist = []
    for filename in filelist:
        pattern = '(' + seqtype + ')\.[Vv]?([0-9]*\.[0-9]*\.?[0-9]*)[\._]?(' + fextension + ')'
        handle = re.search( pattern, filename )
        if handle is not None:
            version = handle.group(2)
            shortlist.append((version,filename))
    if len(shortlist) == 0:
        shortlist = [("","Not_found")]
    return sorted(shortlist, key=itemgetter(0))[-1][1]

def construct_pathways(config_file):
    filelist = find_files()
    pathway_dict = {}
    # get latest file for the OGS and assembly:
    pathway_dict['gff'] = find_latest(filelist, 'OGS', "gff") #
    pathway_dict['cds'] = find_latest(filelist, 'OGS', 'cds')
    pathway_dict['pep'] = find_latest(filelist, 'OGS', 'pep')
    pathway_dict['gtf'] = find_latest(filelist, 'OGS', 'gtf')    
    pathway_dict['gff'] = find_latest(filelist, 'OGS', 'lcl.gff')
    pathway_dict['lclcds'] = find_latest(filelist, 'OGS', 'lcl.cds')
    pathway_dict['lclpep'] = find_latest(filelist, 'OGS', 'lcl.pep')
    pathway_dict['ass'] = find_latest(filelist, 'assembly', 'fa')
    pathway_dict['assgi'] = find_latest(filelist, 'assembly', 'gi.fa')
    pathway_dict['goterms'] = find_latest(filelist, 'OGS', 'GOterms.list')
    # get files for other purposes:
    kegg_patt = '1\.keg'
    orthologs_patt = 'BGI\.orthologs\.list'
    for file in filelist:   
        kegg_h = re.search(kegg_patt, file)
        ortho_h = re.search(orthologs_patt, file)
        if kegg_h is not None:  
            pathway_dict['kegg'] = file
        elif ortho_h is not None:
            pathway_dict['ortho'] = file
    
    write_config(pathway_dict, config_file)     
    return pathway_dict

def import_paths():
    print __file__
    config_file = os.path.dirname(__file__) + "/pathways.config"
    if os.path.exists(config_file):
        pathway_dict = read_pathways(config_file)
    else:
        pathway_dict = construct_pathways(config_file)
    return pathway_dict

if __name__ == "__main__":
    import_paths()