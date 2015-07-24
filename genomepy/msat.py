#!/usr/bin/env python
""" a series of classes and functions for dealing with msats and RadSeq data"""

import re
import sys

from tempfile import mkstemp
from shutil import move
from os import remove, close
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna

########################

class Msat(object):
    """ A class for creating basic microsatellite objects"""

    def __init__(self, name, scaffold, startpos, endpos):
        self.name = name

        self.scaffold = scaffold
        self.startpos = startpos
        self.endpos = endpos

        if endpos > startpos:
            self.size = endpos - startpos
        else:
            self.length = startpos - endpos

    def __str__(self):
        return "%s:%s(%d-%d)" % (self.name,self.scaffold, self.startpos, self.endpos)

    __repr__ = __str__

    def __len__(self):
        return (self.size)

class Radsat(object):
    """ an object representing an individual's genotype at a small region of the genome,
    based on RadSeq data

    loci must be a list of positions (as integers)
    gtypes must be a list of strings "AG" or "CC" etc
    """
    def __init__(self, indiv, scaffold, loci, gtypes, equiv="unassigned"):

        # check that input is correct:
        if len(loci) != len(gtypes):
            print "Number of loci does not equal the number of genotypes for Radsat", indiv
        else:
            self.gtypes = dict(zip(loci, gtypes))

        # calculate % loci heterozygous:
        hetloc = 0
        homloc = 0
        for locus in self.gtypes:
            if self.gtypes[locus][0] == self.gtypes[locus][1]:
                homloc += 1
            else:
                hetloc += 1
        self.hetpc = 100.0 * hetloc / (hetloc + homloc)

        # define other attributes:
        self.id = equiv
        self.size = len(self.gtypes)
        self.name = indiv
        self.scaffold = scaffold
        self.startpos =  min(loci)
        self.endpos = max(loci)
        self.length = self.endpos - self.startpos

    def __str__(self):
        return "%s\t%s:%s(%d-%d) het:%.2f%%" % (self.name, self.id, self.scaffold, self.startpos, self.endpos, self.hetpc)

    __repr__ = __str__

########################

def window_wiping(file_path):
    "For file in file_path, this fn replaces windows line returns with unix returns"
    #Create temp file
    fh, abs_path = mkstemp()
    new_file = open(abs_path,'w')
    old_file = open(file_path)

    for line in old_file:
        new_file.write(line.replace("\r", "\n"))

    #close temp file
    new_file.close()
    close(fh)
    old_file.close()
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

def upload_msats(filename):
    """ creates a list of msat objects"""
    mdict = {}

    robj = open(filename, 'rb')
    rgen = [ (line.split()[0], line.split()[1], int(line.split()[2]), int(line.split()[3])) for line in robj]
    robj.close()

    for mname, mscaf, mstart, mstop in rgen:
        mdict[mname] = Msat(mname, mscaf, mstart, mstop)

    return mdict

def create_pcr_product(filename):
    "Given an input file with left and right primers, creates a false PCR product"
    fasta_handle = open(filename, 'rb')
    newfile_path = filename + "_PCRproduct.fa"
    newfile_handle = open(newfile_path, 'w')
    for line in fasta_handle:
        def_match = re.search("(>)", line)
        left_match = re.search("LEFT", line)
        right_match = re.search("RIGHT", line)
        if left_match is not None:
            switch = False
        if right_match is not None:
            switch = True

        if left_match is not None:
            defline = line
        elif switch is False:
            left_seq = Seq(line.strip(), generic_dna)
        else:
            right_seq = Seq(line.strip(), generic_dna).reverse_complement()
            newfile_handle.write( "%s%sNNNNNNNNNNNNNNNNNNNN%s\n"  %  (defline, left_seq, right_seq) )
    fasta_handle.close()
    newfile_handle.close()

def run_radsat(gtypefile="No_solitary_loh.PILSIP.iGTYPES.gtype"):
    cmmd, filename = sys.argv
    msats = upload_msats(filename)


    radsat_dict = {}

    for msat_key in msats:
        radsat_locusdict = {}
        radsat_gtypedict = {}

        msat = msats[msat_key]
        print "creating Radsats for msat", msat

        gtobj = open(gtypefile, 'rb')

        gtgen = [ line.split() for line in gtobj if line.split()[5] >= 10 and line.split()[6] >= 10 ]
        gtobj.close()

        for details in gtgen:
            if details[1] == msat.scaffold and  msat.startpos - 10000 < int(details[2]) < msat.endpos + 10000:
                try:
                    radsat_locusdict[details[0]].append(int(details[2]))
                    radsat_gtypedict[details[0]].append(details[3]+details[4])
                except:
                    radsat_locusdict[details[0]] = [int(details[2])]
                    radsat_gtypedict[details[0]] = [details[3]+details[4]]

        # create Radsats:
        for indiv in radsat_locusdict:
            radsat_dict[(indiv,msat.name)] = Radsat(indiv, msat.scaffold, radsat_locusdict[indiv], radsat_gtypedict[indiv], equiv=msat.name)

    for radsat in radsat_dict:
        print radsat_dict[radsat], radsat_dict[radsat].size

########################

if __name__ == '__main__':

    cmmd, filename = sys.argv
    create_pcr_product(filename)


