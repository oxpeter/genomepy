#!/usr/bin/env python

import re
import os
import sys
import math

import argparse

def count_deflines(fastafile):
    "counts number of sequences are in a fasta file"
    fasta_h = open(fastafile, 'rb')

    counter = 0
    for line in fasta_h:
        if re.search('^>', line) is not None:
            counter += 1

    fasta_h.close()
    return counter

def split_fasta(fastafile, numfiles):
    "splits fastafile into numfiles even sized fastafiles"

    numseqs = count_deflines(fastafile)
    seqlimit = math.ceil( 1. * numseqs / numfiles ) # num seqs per split file

    fasta_h = open(fastafile, 'rb')

    line = ''
    for f in range(numfiles):
        filepref = os.path.splitext(fastafile)[0]
        fasta_f = open('.'.join([filepref,str(f),'fasta']), 'w')
        counter = 0
        fasta_f.write(line)
        for line in fasta_h:
            if re.search('^>', line) is not None:
                counter += 1
                if counter == seqlimit:
                    break
            fasta_f.write(line)
        fasta_f.close()
    fasta_h.close

def blastall(fastafile, numfiles, database, blastype='blastp'):
    "does blast of split fastafiles against database"
    for f in range(numfiles):
        filepref = os.path.splitext(fastafile)[0]
        fasta_f = '.'.join([filepref,str(f),'fasta'])
        cmd =   blastype + ' -db ' + database + \
                ' -query ' + fasta_f + \
                ' -outfmt 6 -out ' + filepref + str(f) + '.blastp.tsv &'
        os.system(cmd)
        # better tracking of this could be achieved using os.fork() and
        # dropping the '&', but this is beyond my current abilities





if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Speeds up all v all blastp search")

    # input options
    parser.add_argument("-I", "--input_file", type=str, help="The peptide fasta file (query file)")
    parser.add_argument("-D", "--database", type=str, default='/home/antqueen/genomics/indices/orthomcl/CbirAmelDmel', help="The blast database to use (target db)")
    parser.add_argument("-b", "--blast_type", type=str, default='blastp', help="The blast algorithm to use. (default = blastp)")
    parser.add_argument("-p", "--num_threads", type=int, default=1, help="number of threads to distribute blast over")

    args = parser.parse_args()


    ## parse files to set the working directory for saving files
    #  parse input file:
    fullname = os.path.realpath(args.input_file)
    filename = os.path.basename(args.input_file)
    filepath = os.path.dirname(os.path.realpath(args.input_file))
    # parse database path:
    dbfull   = os.path.realpath(args.database)
    # parse blast output name and dir:
    filepref = os.path.splitext(fullname)[0]

    print "splitting %s into %d files..." % (filename, args.num_threads)
    split_fasta(fullname, args.num_threads)
    print "split fasta files saved in dir: %s" % (filepath)
    print "running blastp for all files"
    print "results saved as %s.##.blastp.tsv" % (filepref)
    blastall(fullname, args.num_threads, dbfull, blastype=args.blast_type)

