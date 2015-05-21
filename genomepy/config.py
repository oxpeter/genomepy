#!/usr/bin/python

""" A module to build the file paths of all database files needed for the other modules.
It first looks for the config file pathways.config, and if it cannot find that, then it
builds one itself.

Later improvements should include storing the object as a pickled object for faster
access.
"""

import os, sys
import re
import time
import cPickle

from subprocess import Popen, PIPE
from operator import itemgetter
from pkg_resources import resource_filename, resource_exists
#import pkgutil

def verbalise(arg1, *args):
    # define escape code: '\x1b[31m  %s  \x1b[0m'
    colordict = {'R':'\x1b[31m', 'G':'\x1b[32m',
         'Y':'\x1b[33m' ,'B':'\x1b[34m', 'M':'\x1b[35m' , 'C':'\x1b[36m' }
    if arg1 in colordict:
        argstring = " ".join([str(arg) for arg in args])
        if sys.stdout.isatty():
            color_code = colordict[arg1]
            end_color = '\x1b[0m'
        else:
            color_code = ""
            end_color = ""
    else:
        argstring = " ".join([arg1] + [arg for arg in args])
        color_code = ""
        end_color = ""

    print "%s%s%s" % (color_code, argstring, end_color)

def check_verbose(v=True):
    "allow optional printing with color conversion capability!"
    global verbalise
    if v:
        verbalise = verbalise
    else:
        verbalise = lambda *a: None

    return verbalise

def read_pathways(config_file):
    """
    Extract the pathways of the filetypes specified in the config file.
    """
    pathway_dict = {}
    config_h = open( config_file, 'rb')
    # split each line into two (variable and path), if no comment symbol # is present:
    config_g = ( (line.split()[0], line.split()[1]) for line in config_h if re.search('#',line) == None)
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
        pattern = '(' + seqtype + ')\.[Vv]?([0-9]*\.[0-9]*\.?[0-9]*)[\._]?(' + fextension + ')$'
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
    pathway_dict['gff'] = find_latest(filelist, 'OGS', "gff")
    pathway_dict['cds'] = find_latest(filelist, 'OGS', 'cds')
    pathway_dict['pep'] = find_latest(filelist, 'OGS', 'pep')
    pathway_dict['gtf'] = find_latest(filelist, 'OGS', 'gtf')
    pathway_dict['lclgff'] = find_latest(filelist, 'OGS', 'lcl.gff')
    pathway_dict['lclcds'] = find_latest(filelist, 'OGS', 'lcl.cds')
    pathway_dict['lclpep'] = find_latest(filelist, 'OGS', 'lcl.pep')
    pathway_dict['ass'] = find_latest(filelist, 'assembly', 'fa')
    pathway_dict['assgi'] = find_latest(filelist, 'assembly', 'gi.fa')
    pathway_dict['assone'] = find_latest(filelist, 'assembly', 'singleline.fa')
    pathway_dict['goterms'] = find_latest(filelist, 'OGS', 'GOterms.list')
    pathway_dict['ncbipep'] = find_latest(filelist, 'OGS', 'ncbi.annotated.pep')
    pathway_dict['iprlist'] = find_latest(filelist, 'iprscan', 'gene.list')
    # get files for other purposes:
    kegg_patt = '1\.keg'
    kegg_convert_patt = 'KEGG_orthologs.list'
    orthologs_patt = 'BGI\.orthologs\.list'
    kegg_ps_patt = '2\.sub\.info'
    for file in filelist:
        kegg_h = re.search(kegg_patt, file)
        ortho_h = re.search(orthologs_patt, file)
        keggc_h = re.search(kegg_convert_patt, file)
        keggp_h = re.search(kegg_ps_patt, file)
        if kegg_h is not None:
            pathway_dict['kegg'] = file
        elif ortho_h is not None:
            pathway_dict['ortho'] = file
        elif keggc_h is not None:
            pathway_dict['keggortho'] = file
        elif keggp_h is not None:
            pathway_dict['keggpathways'] = file
    write_config(pathway_dict, config_file)
    return pathway_dict

def rebuild_config():
    config_path = resource_filename('genomepy', 'data/pathways.cfg')
    print "Rebuilding", config_path
    pathway_dict = construct_pathways(config_path)
    return pathway_dict

def import_paths():
    #config_file = pkgutil.get_data('genomepy', 'data/pathways.cfg')
    #print config_file
    #config_str = resource_string('genomepy', 'data/pathways.cfg')
    #nothing_str = resource_string('genomepy', 'data/nothing.cfg')
    config_exists = resource_exists('genomepy', 'data/pathways.cfg')
    config_path = resource_filename('genomepy', 'data/pathways.cfg')

    print "Using config file from %s" % (config_path)
    if config_exists:
        pathway_dict = read_pathways(config_path)
    else:
        pathway_dict = construct_pathways(config_path)
    return pathway_dict


##### FILE MANIPULATION #####
def pickle_jar(obj, fname):
    pklfile = ".".join([fname, 'pkl'])
    apicklefile = open(pklfile, 'wb')
    cPickle.dump(obj, apicklefile, -1)
    apicklefile.close()

def open_pickle_jar(fname):
    pklfile = ".".join([fname, 'pkl'])
    apicklefile = open(pklfile, 'rb')
    loaded_obj = cPickle.load(apicklefile)
    apicklefile.close()
    return loaded_obj

def pickle_gtypes(lineage_gtypes, indiv_gtypes, fname='pickle'):
    # try to pickle the gtypes!
    pklfile = ".".join([fname, 'aGtypes.pkl'])
    apicklefile = open(pklfile, 'wb')
    cPickle.dump(lineage_gtypes, apicklefile, -1)
    apicklefile.close()

    pklfile = ".".join([fname, 'iGtypes.pkl'])
    ipicklefile = open(pklfile, 'wb')
    cPickle.dump(indiv_gtypes, ipicklefile, -1)
    ipicklefile.close()

def unpickle_gtypes(fname='pickle'):
    # and now to unpickle!
    pklfile = ".".join([fname, 'aGtypes.pkl'])
    apicklefile = open(pklfile, 'rb')
    lineage_gtypes = cPickle.load(apicklefile)
    apicklefile.close()

    pklfile = ".".join([fname, 'iGtypes.pkl'])
    ipicklefile = open(pklfile, 'rb')
    indiv_gtypes = cPickle.load(ipicklefile)
    ipicklefile.close()

    return lineage_gtypes, indiv_gtypes

def file_block(filehandle,  block, number_of_blocks=1000):
    """
    This code adapted from:
    http://xor0110.wordpress.com/2013/04/13/how-to-read-a-chunk-of-lines-from-a-file-in-python/

    Written by Nic Werneck

    A generator that splits a file into blocks and iterates
    over the lines of one of the blocks.

    usage:
    filehandle = open(filename)
    number_of_chunks = 100
    for chunk_number in range(number_of_chunks):
        for line in file_block(filehandle, number_of_chunks, chunk_number):
            process(line)
    """


    assert 0 <= block and block < number_of_blocks
    assert 0 < number_of_blocks

    filehandle.seek(0,2)
    file_size = filehandle.tell()

    ini = file_size * block / number_of_blocks
    end = file_size * (1 + block) / number_of_blocks

    if ini <= 0:
        filehandle.seek(0)
    else:
        filehandle.seek(ini-1)
        filehandle.readline()

    while filehandle.tell() < end:
        yield filehandle.readline()

def create_log(args, outdir=None, outname='results'):
    ## create output folder and log file of arguments:
    timestamp = time.strftime("%b%d_%H.%M")
    if outdir:
        newfolder = os.path.realpath(outdir)
        if os.path.exists(newfolder) is False:  # check to see if folder already exists...
            os.mkdir(newfolder)
        filename = newfolder + '/' + outname + '.' + timestamp + ".log"
    else:
        root_dir = os.getcwd()
        newfolder = root_dir + "/" + outname + "." + timestamp
        os.mkdir(newfolder)  # don't have to check if it exists, as timestamp is unique
        filename = newfolder + "/" + outname + "." + timestamp + ".log"

    log_h = open(filename, 'w')
    log_h.write( "File created on %s\n" % (timestamp) )
    log_h.write( "Program called from %s\n" % (os.getcwd()) )
    log_h.write( "%s\n\n" % (' '.join(sys.argv)) )
    for arg in str(args)[10:-1].split():
        log_h.write( "%s\n" % (arg) )
    log_h.close()
    return filename

def make_a_list(geneobj, col_num=0, readfile=True):
    """given a path, list, dictionary or string, convert into a list of genes.
    col_num specifies the column from which to extract the gene list from."""

    # genefile can be given as a list of genes (say, from find_degs()... ), or as
    # a path to a file containing a list of genes.
    # The following builds a dictionary of genes from either input:
    if type(geneobj) is list:   # allows filtering from hicluster generated list of results.
        genelist = {}.fromkeys(geneobj,1)
    elif type(geneobj) is dict:
        genelist = geneobj # this of course will only work if the keys are the genes!
    elif type(geneobj) is str:   # assuming a filepath...
        if re.search("/",geneobj) is None:
            genelist = {}.fromkeys(geneobj.split(','),1)
        else:   # is a file path
            if re.search(",",geneobj) or not readfile: # can turn of so that a single file path is returned as the string
                genelist = {}.fromkeys(geneobj.split(','),1)
            else:
                genefile_h = open(geneobj, 'rb')
                genelist = {}   # will be a dict of names { 'Cbir01255':1, 'CbirVgq':1, ... }
                                # used a dictionary to automatically remove any name duplications
                filegen = [ line.split() for line in genefile_h if len(line) > 0]

                genefile_h.close()

                for colset in filegen:
                    try:
                        genelist[colset[col_num]]=1
                    except IndexError:
                        verbalise("R", "Column %d not found in %s" % (col_num, str(colset)))
    else:
        genelist = {}

    return genelist


if __name__ == "__main__":
    rebuild_config()