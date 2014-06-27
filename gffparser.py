#!/usr/bin/env python
""" A module for parsing GFF files """

import sys
import cPickle
import re
import os

import argparse
from BCBio import GFF
from Bio import SeqIO
from Bio.Alphabet import generic_dna, IUPAC
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
#from reportlab.lib import colors
#from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
import progressbar              # from Nilton Volpato

import gff2bed2
import genematch

####### CLASSES ############################################################

class Splicer(object):
    """a class that allows analysis of splice junctions from a gff file and tophat
    junction files"""
    def __init__(self, gff_file="/Volumes/Genome/armyant.OGS.V1.8.6_lcl.gff", chosenscaffold='All'):
        ## construct canonical and alternate splice junctions from gff file:
        ## saved in 2 dictionaries: canonical = { 'scaffold2':{(5292,5344):'cbir_02775',():,():,():}
        ##                          alternate = { 'scaffold2':{(5292,5671):'cbir_02775',():,():,():}

        # extract information from gff file:
        gff_h = open(gff_file, 'rb')
        print 'reading gff file...'
        if chosenscaffold == 'All':
            #         scaf             exon_start            exon_end              gene_name                      strand
            exons = ((line.split()[0], int(line.split()[3]), int(line.split()[4]), line.split()[8].split('=')[1], line.split()[6]) for line in gff_h if line.split()[2] == 'CDS')
        else:   # can specify a single scaffold to analyse (much faster, hopefully!)
            exons = ((line.split()[0], int(line.split()[3]), int(line.split()[4]), line.split()[8].split('=')[1], line.split()[6]) for line in gff_h if line.split()[2] == 'CDS' and line.split()[0] == chosenscaffold)


        print "gff file read."

        self.canonical = {}
        self.alternative = {}
        self.boundaries = {}

        geneid = 'initiating_string'
        junctionlist = [1,2]  # this is to prevent an error when first starting the for loop

        # setup progress bar for this rather long process that is to follow:
        bar = progressbar.ProgressBar(maxval=100262, \
        widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.ETA()]) # can also use progressbar.Percentage()
        count=0
        bar.update(count)

        for exon in exons:
            # add junction start & stop to list for creation of alternate splice junctions
            if exon[3] == geneid: # ie, still analysing the same gene
                junctionlist += [min(exon[2], exon[1]), max(exon[2], exon[1])]

            else:  # ie, a new gene has started
                ## determine junctions, add to dictionaries and reset list:
                junctionlist.sort() # sort by first element in tuple
                # extract canonical splice sites from list of junctions:
                for posn in range(len(junctionlist)):
                    if posn == 0:
                        continue
                    elif posn % 2 == 0:
                        try:
                            self.canonical[scaffold][(junctionlist[posn-1],junctionlist[posn])] = geneid
                        except KeyError:
                            self.canonical[scaffold] = {(junctionlist[posn-1],junctionlist[posn]):geneid}
                        try:
                            self.boundaries[scaffold][junctionlist[posn-1]] = geneid
                            self.boundaries[scaffold][junctionlist[posn]]   = geneid
                        except KeyError:
                            self.boundaries[scaffold] = {junctionlist[posn-1] : geneid}
                            self.boundaries[scaffold] = {junctionlist[posn] : geneid}

                # extract alternative splice sites from junction list:
                for posn in range(len(junctionlist)):
                    if posn % 2 != 0:
                        for secposn in range(len(junctionlist))[posn:]:
                            if secposn % 2 == 0:
                                try:
                                    self.alternative[scaffold][(junctionlist[posn],junctionlist[secposn])] = geneid
                                except KeyError:
                                    self.alternative[scaffold] = {(junctionlist[posn],junctionlist[secposn]):geneid}

                # reset variables:
                scaffold = exon[0]
                geneid = exon[3]
                junctionlist = [min(exon[2], exon[1]), max(exon[2], exon[1])]
            count+=1
            bar.update(count)
        bar.finish()
        gff_h.close()

    def __str__(self):
        """bigstring = "Canonical and alternate junctions saved:\n"
        bigstring +=  "Scaf\tjunction\tgene\n"
        for scaf in self.canonical.keys()[0:3]:
            for jnc in self.canonical[scaf].keys()[0:6]:
                bigstring += "%-12s %-12s %s\n" % (scaf, jnc, self.canonical[scaf][jnc])"""
        return "%d scaffolds:%d junctions" % (len(self.canonical), len(self.canonical[self.canonical.keys()[0]]))

    __repr__ = __str__

    def map_junctions(self, junc_file, chosenscaffold='All'):
        ## parse junctions from bed file and compare to object's dictionary of junctions

        junc_h = open(junc_file, 'rb')
        #           scaf, frag_start, junc_start, junc_end, frag_end (in gff positioning)
        if chosenscaffold == 'All':
            juncgen = ((line.split()[0],                               # 0) scaf
                        int(line.split()[1]) + 1,                           # 1) frag_start
                        int(line.split()[1]) + int(line.split()[10].split(',')[0]),    # 2) junc_start
                        int(line.split()[1]) + int(line.split()[11].split(',')[1]) + 1,# 3) junc_end
                        int(line.split()[1]) + int(line.split()[11].split(',')[1]) + int(line.split()[10].split(',')[1])
                        ) for line in junc_h if line[0] != 't')
        else:
            juncgen = ((line.split()[0],                               # scaf
                        int(line.split()[1]) + 1,                           # frag_start
                        int(line.split()[1]) + int(line.split()[10].split(',')[0]),    # junc_start
                        int(line.split()[1]) + int(line.split()[11].split(',')[1]) + 1,# junc_end
                        int(line.split()[1]) + int(line.split()[11].split(',')[1]) + int(line.split()[10].split(',')[1])
                        ) for line in junc_h if line.split()[0] == chosenscaffold and line[0] != 't')


        # create dictionaries for storing the results:
        self.canonicals = {}
        self.alternatives = {}
        self.novels = {}
        self.readjusts = {}

        self.canon_genes = {}
        self.alt_genes = {}

        ccount = 0
        acount = 0
        ncount = 0
        rcount = 0
        count = 0

        for junction in juncgen:
            count += 1
            if junction[0] not in self.canonical: # mainly if Splicer constructed from single scaffold
                continue
            elif (junction[2], junction[3]) in self.canonical[junction[0]]:
                ccount += 1
                try:
                    self.canonicals[junction[0]][(junction[2], junction[3])] = self.canonical[junction[0]][(junction[2], junction[3])]
                except KeyError:
                    self.canonicals[junction[0]] = {(junction[2], junction[3]):self.canonical[junction[0]][(junction[2], junction[3])]}
            elif (junction[2], junction[3]) in self.alternative[junction[0]]:
                acount += 1
                geneid = self.alternative[junction[0]][(junction[2], junction[3])]
                try:
                    self.alternatives[junction[0]][(junction[2], junction[3])] = geneid
                except KeyError:
                    self.alternatives[junction[0]] = {(junction[2], junction[3]):geneid}
                try:
                    self.alt_genes[geneid] += [(junction[2], junction[3])]
                except:
                    self.alt_genes[geneid] = [(junction[2], junction[3])]
            elif junction[2] in self.boundaries[junction[0]] or junction[3] in self.boundaries[junction[0]]:
                # check to see if one of the junction boundaries is still canonical:
                rcount += 1
                try:
                    geneid = self.boundaries[junction[0]][junction[2]]
                except KeyError:
                    try:
                        geneid = self.boundaries[junction[0]][junction[3]]
                    except KeyError:
                        geneid = 'Unknown'
                try:
                    self.readjusts[junction[0]][(junction[2], junction[3])] = geneid
                except KeyError:
                    self.readjusts[junction[0]] = {(junction[2], junction[3]):geneid}
                try:
                    self.alt_genes[geneid] += [(junction[2], junction[3])]
                except:
                    self.alt_genes[geneid] = [(junction[2], junction[3])]
            else: # neither boundary matches a canonical boundary
                ncount += 1
                try:
                    self.novels[junction[0]][(junction[2], junction[3])] = True
                except KeyError:
                    self.novels[junction[0]] = {(junction[2], junction[3]):True}
        #print "%d junctions analysed" % (count)
        junc_h.close()
        return ccount, acount, rcount, ncount


class My_gff(object):
    "an object for quick assessment of where a GENE/SNP lies"
    def __init__(self, gff_file="/Volumes/Genome/armyant.OGS.V1.8.6.gff"):

        self.genedict = {}
        self.genecount = 0
        self.geneids = {}
        self.genescaf = {}

        gff_h = open(gff_file,'rb')
        for line in gff_h:
            if line.split()[2] == 'mRNA':
                self.genecount += 1
                scaf = line.split()[0]
                start = min(int(line.split()[3]), int(line.split()[4]))
                stop = max(int(line.split()[3]), int(line.split()[4]))
                geneid = re.search('ID=([^;]*);', line.split('\t')[8]).groups()[0]
                self.genescaf[geneid] = scaf
                try:
                    genename = re.search('Name=([^;]*);', line.split('\t')[8]).groups()[0]
                except:
                    genename = "No_ortholog"
                self.geneids[geneid] = genename
                if scaf in self.genedict:
                    self.genedict[scaf][geneid] = (start, stop)
                else:
                    self.genedict[scaf] = {geneid:(start, stop)}

    def __contains__(self, locus):
        scaf, posn = locus
        ingene = False
        try:
            for gene in self.genedict[scaf]:
                if self.genedict[scaf][gene][0] <= posn <= self.genedict[scaf][gene][1]:
                    ingene = True
        except KeyError:
            ingene = False

        return ingene

    def __str__(self):
        return "%d scaffolds, %d genes" % (len(self.genedict),self.genecount)

    __repr__ = __str__

    def __len__(self):
        return (self.genecount)

    def whichscaf(self, geneid):
        return self.genescaf[geneid]

    def gene(self, locus):
        scaf, posn = locus
        ingene = False
        try:
            for gene in self.genedict[scaf]:
                if self.genedict[scaf][gene][0] <= posn <= self.genedict[scaf][gene][1]:
                    ingene = gene
        except KeyError:
            ingene = None

        return ingene

    def nameit(self, geneid):
        try:
            return self.geneids[geneid]
        except KeyError:
            return "Not_found"

####### FUNCTIONS ############################################################

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


##### GFF FILE MANIPULATION #####

def gff2gtf(gff_file):
    "converts a gff file to gtf"
    gtf_file = os.path.splitext(gff_file)[0] + ".gtf"
    gtf_h = open(gtf_file, 'w')
    gff_h = open(gff_file, 'rb')

    for line in gff_h:
        cds_patt = '(.*)Parent=([A-Za-z0-9_\(\)]*)'
        parentline = re.search(cds_patt, line)
        if parentline is not None:
            newline = parentline.groups()[0] + '\tgene_id "' + parentline.groups()[1] + '"; transcript_id "' + parentline.groups()[1] + ';\n'
            gtf_h.write(newline)
    gtf_h.close()
    gff_h.close()

def gtf2gff(gtf_file):
    "converts a gtf file to gff using perl script"
    gff_file = os.path.splitext(gtf_file)[0] + ".gff"
    cmd = "perl /Users/POxley/scripts/Genome_Analysis/gtf2gff3.pl " + gtf_file + " > " + gff_file
    os.system(cmd)

def extended_gff(gff_file):
    "converts OGS gff file into extended gff with genes and exons"
    gff_h = open(gff_file, 'rb')
    bed_h = open(gff_file[:-4] + ".ex.gff", 'w')

    count = 0
    for line in gff_h:
        count += 1

        fields = line.split()
        scaf = fields[0]
        comment = fields[1]
        type = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        score = fields[5]
        strand = fields[6]
        frame = fields[7]
        definition = " ".join(fields[8:])

        if type == 'mRNA':
            geneid_search = re.search("ID=(Cbir[A-Za-z0-9_\.\(\)\/]*);", line)
            geneid = geneid_search.group(1)

            genefields = [scaf, comment, 'gene', str(start), str(end), score, \
                        strand, frame, definition.replace('ID=','ID=g_'),"\n" ]
            mrnaline = line[:-1] + "Parent=g_" + geneid + ";\n"

            bed_h.write("\t".join(genefields))
            bed_h.write(mrnaline)
            count = 0
        if type == 'CDS':
            geneid = re.search("Parent=(Cbir[A-Za-z0-9_\.\(\)\/]*);", line).group(1)
            bed_h.write(line)
            exonline = line.replace("CDS", "exon").replace("Parent=", "ID=" + geneid + "_" + str(count) + ";Parent=")
            bed_h.write(exonline)
    bed_h.close()
    gff_h.close()

def make_bed(gff_file):
    """ Takes a gff_file and creates a bed format including the scaffolds/contigs that
    do not have any features associated with them"""

    gff_fh = open(gff_file, 'rU')

    # get transcript annotation
    tinfo, einfo = gff2bed2.ParseAnno(gff_fh)
    # write into bed format
    scaf_dict = gff2bed2.WriteBED(tinfo, einfo, return_obj=True, gff_name=gff_file)

    # create dictionary of scaffolds not annotated, and their sizes
    missing_dict = {}
    assembly_handle = open("/Volumes/Genome/Cbir.assembly.v3.0_singleline.fa", 'r')
    for line in assembly_handle:
        def_search = re.search('>(lcl\|[a-zC0-9]*)', line)
        line_len = len(line)
        if def_search is not None:
            scaf_id = def_search.groups()[0]
        else:
            missing_dict[scaf_id] = line_len

    # append these onto the bed file, using the first 3 columns
    bed_name = gff_file + ".bed"
    bed_handle = open(bed_name, 'a')
    for scaf in missing_dict:
        bed_handle.write( "%s\t%d\t%d\tNonGene\t.\t+\t%d\t%d\t0\t1\t%d\t%d\n" % (scaf, 0, missing_dict[scaf], 0, missing_dict[scaf], 0, missing_dict[scaf]) )
    bed_handle.close()

def trim_untranslated(gff_file):
    """removes 5' and 3' untranslated sequence information from gff file for later
    comparisons of different isoforms.
    gff file MUST be sorted for this to work."""

    gff_h = open(gff_file, 'rb')
    newgff = open(gff_file[:-4] + ".trimmed.gff", 'w')
    firstline = True
    minpos = 1000000000
    maxpos = -10
    cds_lines = ""
    exon_lines = ""

    for line in gff_h:
        if len(line) > 1:
            feature = line.split()[2]
        else:
            continue
        if feature == "exon":
            prevexon = line.split()
        elif feature == "CDS":
            thiscds = line.split()
            if min(int(thiscds[3]), int(thiscds[4])) < minpos:
                minpos = min(int(thiscds[3]), int(thiscds[4]))
            if max(int(thiscds[3]), int(thiscds[4])) > maxpos:
                maxpos = max(int(thiscds[3]), int(thiscds[4]))
            cds_lines += line
            if prevexon[3] != thiscds[3]:
                prevexon[3] = thiscds[3]
            if prevexon[4] != thiscds[4]:
                prevexon[4] = thiscds[4]
            exon_lines += "\t".join(prevexon) + "\n"
        elif firstline:
            gene_line = line.split()
            firstline = False
        elif feature == "gene":
            # print features of previous gene:
            gene_line[3] = str(minpos)
            gene_line[4] = str(maxpos)
            mrna_line[3] = str(minpos)
            mrna_line[4] = str(maxpos)
            newgff.write("\t".join(gene_line) + "\n")
            newgff.write("\t".join(mrna_line) + "\n")
            newgff.write(exon_lines)
            newgff.write(cds_lines)

            # reset values for next gene:
            minpos = 1000000000
            maxpos = -10
            gene_line = line.split()
            cds_lines = ""
            exon_lines = ""

        elif feature == "mRNA":
            mrna_line = line.split()

        else:   # ie, if three_prime_UTR or five_prime_UTR...
            continue
    # write final gene to file:
    gene_line[3] = str(minpos)
    gene_line[4] = str(maxpos)
    mrna_line[3] = str(minpos)
    mrna_line[4] = str(maxpos)
    newgff.write("\t".join(gene_line) + "\n")
    newgff.write("\t".join(mrna_line) + "\n")
    newgff.write(exon_lines)
    newgff.write(cds_lines)

    newgff.close()
    gff_h.close()

def fix_gtf(gtf_file):
    "edits cufflinks gtf file to allow Transdecoder to maintain gene identity"
    gtf_h = open(gtf_file, 'rb')
    newgtf = open(gtf_file[:-4] + ".edit.gtf", 'w')
    for line in gtf_h:
        geneid = re.search('gene_id "([^"]*)"', line).group(1)
        isoform = re.search('transcript_id "([^"]*)"', line).group(1)
        if geneid == "":
            newline = line.replace('gene_id "','gene_id "' +  isoform)
        elif re.search("Cbir", isoform) is not None:
            newline = line.replace('transcript_id "', 'transcript_id "' + geneid + ".")
        else:
            newline = line
        newgtf.write(newline)
    gtf_h.close()
    newgtf.close()

def bed2gtf(bedfile):
    """convert Transdecoder bedfiles to gtf and gff
    prefix indicates the geneid prefix used by cufflinks
    """

    bed_h = open(bedfile, 'rb')
    gtf_h = open(bedfile[:-4] + ".gtf", 'w')
    gff_h = open(bedfile[:-4] + ".gff", 'w')

    genenames = {}
    genesizes = {}
    genestrands = {}



    for line in bed_h:
        if len(line.split()) < 6:
            continue
        fields = line.split()
        scaf = fields[0]
        feat = 'transdecoder'
        gene_start  = int(fields[1]) + 1
        gene_end    = int(fields[2])
        strand      = fields[5]
        labels      = fields[3].split(";")

        try:
            gene_iso    = re.search("(\w+\.\w+)\.([A-Za-z0-9_\(\)\.]*)", labels[1])
            unexpressed = re.search("ID=(Cbir[A-Za-z0-9_\(\)\.]*)", labels[0]) # genes that are not expressed won't have the cufflinks id
        except IndexError:
            print fields[3], labels

        if unexpressed is not None:
            geneid      = unexpressed.group(1)
            isoform     = geneid
        else:
            try:
                geneid      = gene_iso.group(1)
                isoform     = gene_iso.group(0)
            except:
                print line
                print labels[0]
                die


        # update gene information:
        try:
            genenames[geneid] += isoform + ">"
        except KeyError:
            genenames[geneid] = isoform + ">"

        if geneid in genesizes:
            if min(gene_start, gene_end) < min(genesizes[geneid]):
                genesizes[geneid][genesizes[geneid].index(min(genesizes[geneid]))] = min(gene_start, gene_end)
            if max(gene_start, gene_end) > max(genesizes[geneid]):
                genesizes[geneid][genesizes[geneid].index(max(genesizes[geneid]))] = max(gene_start, gene_end)
        else:
            genesizes[geneid] = [gene_start, gene_end]

        genestrands[geneid] = strand

        # get exon info
        exon_starts = [int(x) for x in fields[11].split(",")]
        exon_lengths= [int(x) for x in fields[10].split(",")]
        exon_number = 1

        # create gtf and gff files:
        transcript_line = "\t".join([scaf, "gffparser", "transcript", str(gene_start), str(gene_end), '1', strand, ".", 'gene_id "' + geneid + '"; transcript_id "' + isoform + '";\n'])
        mRNA_line = "\t".join([scaf, "gffparser", "mRNA", str(gene_start), str(gene_end), ".", strand, ".", "ID=" + isoform + ";Parent=" + geneid + "\n"])
        gtf_h.write(transcript_line)
        gff_h.write(mRNA_line)

        for exon_start, exon_length in zip(exon_starts, exon_lengths):
            gtf_line = '\t'.join([scaf, "gffparser", "exon", str(exon_start + gene_start), str(exon_start + gene_start + exon_length - 1), '1', strand, '.', 'gene_id "' + geneid + '"; transcript_id "' + isoform + '"; exon_number "' + str(exon_number) + '";\n'])
            gff_CDS  = '\t'.join([scaf, "gffparser", "CDS",  str(exon_start + gene_start), str(exon_start + gene_start + exon_length - 1), '.', strand, '.', 'ID=cds.' + str(exon_number) + isoform + ';Parent=' + isoform + '\n'])
            gtf_h.write(gtf_line)
            gff_h.write(gff_CDS)
            exon_number += 1


    #convert gene name to Cbir version if available:
    for geneid in genenames:
        try:
            genename = re.search("(Cbir[A-Za-z0-9_\(\)]+)", genenames[geneid]).group(1)
        except:
            genename = geneid
        gene_line = '\t'.join([scaf, "gffparser", "gene", str(min(genesizes[geneid])), str(max(genesizes[geneid])), '.', genestrands[geneid], '.', 'ID=' + geneid + ';Name=' + genename + '\n'])
        gff_h.write(gene_line)

    gff_h.close()
    gtf_h.close()

def strip_duplicates(bedfile):
    "removes duplicate gene entries, and replacing ID with the cufflinks appended Cbir ID"
    bed_h = open(bedfile, 'rb')
    newbed_h = open(bedfile[:-3] + "uniq.bed", 'w')
    for line in bed_h:
        fields = line.split()
        if fields[12] == ".":
            newbed_h.write("\t".join(fields[:12]) + "\n")
            continue
        cbirID = re.search("Cbir", fields[15])
        if cbirID is not None:
            newbed_h.write("\t".join(fields[12:]) + "\n")
    bed_h.close()
    newbed_h.close()

def parse_names(genelist, gffobj):
    "takes a gene ID or list of gene IDs and returns the gene name from the gff file"
    ## Create dictionary of gene ids and their corresponding names:
    namedict = {}

    for scaf in gffobj:
        for feature in gffobj[scaf].features:
            try:
                namedict[feature.qualifiers['ID'][0]] = feature.qualifiers['Name'][0]
            except:
                #print feature.qualifiers['ID']
                namedict[feature.qualifiers['ID'][0]] = feature.qualifiers['ID'][0]

    ## Pull out names from supplied gene list:
    output_dict = {}
    if type(genelist) == list:
        for g in genelist:
            try:
                output_dict[g] = namedict[g]
            except KeyError:
                print g, "was not found."
    else:
        try:
            output_dict[genelist] = namedict[genelist]
        except KeyError:
            print gene, "was not found."

    return output_dict

def assemble_dict(in_file="/Volumes/Genome/armyant.OGS.V1.8.6.gff", in_seq_file="/Volumes/Genome/Cbir.assembly.v3.0.fa", features_only=False):
    """ Takes a gff file and the genome assembly and creates a SeqRecord dictionary.

    INPUT:
    in_file = path/to/file.gff
    in_seq_file = path/to/assembly.file.fa

    OUTPUT:
    gff_dict['scaffold'] = SeqRecord_obj
    """
    in_seq_handle = open(in_seq_file)
    seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
    in_seq_handle.close()

    gff_dict = {}

    in_handle = open(in_file)

    if features_only:
        for rec in GFF.parse(in_handle, base_dict=seq_dict):
            if len(rec.features)>=1:
                gff_dict[rec.id] = rec
    else:
        for rec in GFF.parse(in_handle, base_dict=seq_dict):
            gff_dict[rec.id] = rec

    return gff_dict

def create_gff(gff_dict, out_file):
    """ Given the gff dictionary, this function creates a new GFF file"""
    print "Writing to", out_file

    out_handle = open(out_file, 'w')
    rec_list = []

    for rec in gff_dict:
        rec_list.append(gff_dict[rec])

    GFF.write(rec_list, out_handle)

    out_handle.close()

def isolate_mRNA(gff_p):
    """takes a gff file and removes the gene header, allowing conversion and viewing of
    all isoforms on JBrowse"""

    newfile =  os.path.splitext(gff_p)[0] + ".mRNA.gff"

    gff_h = open(gff_p, 'rb')
    newfile_h = open(newfile, 'w')

    for line in gff_h:
        if line.split()[2] == "mRNA":
            defsearch = re.search("(.*);Parent=.+", line)
            newfile_h.write(defsearch.groups()[0] + "\n")
        elif line.split()[2] == "exon":
            newfile_h.write(line)

    gff_h.close()
    newfile_h.close()
    return newfile

def reinstate_cbir(cuffcompare_gtf):
    "for valid transcripts, replace cuffcompare id with original Cbir id."
    gtf_h = open(cuffcompare_gtf, 'rb')
    new_gtf = open(cuffcompare_gtf[:-4] + ".reins.gtf", 'w')
    details = {}
    for line in gtf_h:
        fields = line.split('\t')
        detail_list = fields[8].strip().split(';')
        try:
            details = dict([(x.strip().split(' ')[0], x.strip().split(' ')[1]) for x in detail_list if x != ''])
        except IndexError:
            print detail_list
            die

        geneid = re.search("(Cbir[A-Za-z0-9_\(\)\.\/]*)", details['oId'])
        if geneid is not None:
            details['transcript_id'] = '"%s"' % (geneid.group(1))

        new_attributes = 'gene_id %(gene_id)s; transcript_id %(transcript_id)s;' % details
        for attrib in details:
            if attrib != 'gene_id' and attrib != 'transcript_id':
                new_attributes += ' %s %s;' % (attrib, details[attrib])
        new_fields = "\t".join(fields[:8]) + '\t' + new_attributes + '\n'
        new_gtf.write(new_fields)
    gtf_h.close()
    new_gtf.close()

##### SEQRECORD FUNCTIONS #####

def count_genes(gff_dict):
    "Returns the number of features (ie genes) in a SeqRecord obj"
    count = 0
    for rec in gff_dict:
        count += len(gff_dict[rec].features)
    return count

def show_features(gff_dict, scaf):
    "Reports the features found in a given SeqRecord for the given scaffold"
    for feature in gff_dict[scaf].features:
        print "Gene:", feature.qualifiers['ID'][0],
        try:
            print "%s\t%s" % (feature.qualifiers['Name'][0], feature.location)
        except KeyError:
            print feature.location

def gene_info(gff_dict, geneid):
    'returns salient features of a given gene'


def snp_in_gene(scaf, pos, gff_dict):
    """returns salient features of a given SNP.

    INPUT:
    scaf = scaffold SNP is located on
    pos = position in scaffold of SNP (NB uses python base 0!!!)
    gff_dict = dictionary of SeqRecords (must contain the given scaffold as a key)

    OUTPUT:
    cds_pos = the distance into the coding sequence of the SNP
    exon =  the exon the SNP is located on
    codon = the nt sequence of the codon the SNP is on
    frame = the position the SNP is in within the codon (position 1,2 or 3)
    ref_nt = the nucleotide in the reference assembly at the SNP position
    """



    for feat in gff_dict[scaf].features:
        if pos in feat:

            # print gene details for match:
            #print feat.qualifiers['ID'][0],
            try:
                gene_name = feat.qualifiers['Name'][0]
            except KeyError:
                gene_name = None

            try:
                gene_id   = feat.qualifiers['ID'][0]
            except KeyError:
                gene_id   = None

            # collect info about frame posn, distance into gene etc:
            before_len = 0
            before_cnt = 0
            for subfeat in feat.sub_features:
                if pos in subfeat.location:     # ie, pos is in this exon:
                    if feat.strand == 1:
                        exon_pos = pos - subfeat.location.start + 1
                    else:
                        exon_pos = subfeat.location.end - pos
                elif feat.strand == 1 and subfeat.location.start < pos: # ie, exon precedes the pos' exon:
                    before_len += len(subfeat.location)
                    before_cnt += 1
                elif feat.strand == -1 and subfeat.location.start > pos: # ie, exon precedes the pos' exon:
                    before_len += len(subfeat.location)
                    before_cnt += 1

            cds_pos = exon_pos + before_len
            frame =  {1:1,2:2,0:3}[cds_pos % 3]
            seq_start = cds_pos - frame
            seq_end = seq_start + 3
            exon = before_cnt +  1
            #print "SNP position: %d\texon: %d\t" % (cds_pos, exon)


            rec_seq = gff_dict[scaf].seq
            if feat.strand == 1:
                codon = feat.extract(rec_seq)[seq_start:seq_end]
                ref_nt = rec_seq[pos:pos+1]
                #print "Codon: %s\t Codon Position: %d\t Ref. nt: %s\t" % (codon, frame, ref_nt)
            else:
                feat_seq = get_sequence(scaf, gff_dict, feat)
                codon = feat_seq[seq_start:seq_end]
                ref_nt = rec_seq[pos:pos+1].reverse_complement()
                #print "Codon: %s\t Codon Position: %d\t Ref. nt: %s\t" % (codon, frame, ref_nt)
            #print "SNP lies", feat.location.end - pos + 1, "from the start of the gene (including introns)."
            SNP_dict = {"gene_id":gene_id, "gene_name":gene_name, "cds_pos":cds_pos, "exon":exon, "codon":codon, 'codon_str':str(codon), "frame":frame, "ref_nt":ref_nt, 'ref_nt_str':str(ref_nt)}
            break
    else:
        SNP_dict = {"gene_id":None, "gene_name":None, "cds_pos":None, "exon":None, "codon":None, 'codon_str':None, "frame":None, "ref_nt":None, 'ref_nt_str':None}
        # It may be necessary to report -ve results in the future. Perhaps something like:
        #else:
        #    SNP_dict[feat] = (None, None, None, None, None)

    return SNP_dict

def get_sequence(scaf, gff_dict, feature):
    """ There appears to be a bug in the .extract function of features, such that genes
    in the reverse strand have their sequences reverse complemented correctly, but the
    exons can appear in the wrong order. This is my attempt to fix this."""

    sequence = ""
    if feature.strand == 1:
        sequence = feature.extract(gff_dict[scaf].seq)
    else:
        subseq_dict = {}
        for subfeature in feature.sub_features:
            subseq_dict[int(subfeature.location.start)] = subfeature.extract(gff_dict[scaf].seq)
        for exon in sorted(subseq_dict.keys(), reverse=True):
            sequence += subseq_dict[exon]

    sequence.alphabet = IUPAC.ambiguous_dna
    return sequence

def parse_go(gene, gofile='/Volumes/Genome/Genome_analysis/Gene_Ontology/armyant.OGS.V1.5.GOterms.list'):
    "for a given gene or genelist, returns a dictionary of all GO terms associated with the gene"

    # make dictionary of go terms:
    go_dict = {}
    go_handle = open(gofile, 'rb')
    columns = [ line.split('\t') for line in go_handle ]
    for columnset in columns:
        go_dict[columnset[0]] = {}
        for element in columnset[2:]:
            gopattern = "GO:([0-9]*)"
            defpattern = "[a-z].*"
            rego = re.search(gopattern, element)
            defgo = re.search(defpattern, element)
            if rego is not None and defgo is not None:
                defline = defgo.group().split(';')
                #print defline
                gotype = defline[1].split(' ')[1][0] + defline[1].split(' ')[1][1]
                godef =  defline[0]
                go_dict[columnset[0]][rego.group()] =  (godef, gotype)

    #for element in go_dict:
    #    print element, go_dict[element]

    #determine if one gene or many:
    output_dict = {}
    if type(gene) == list:
        for g in gene:
            try:
                output_dict[g] = go_dict[g]
            except KeyError:
                #print g, "was not found."
                output_dict[g] = {"GO:######":("None listed","NA")}
    else:
        try:
            output_dict[gene] = go_dict[gene]
        except KeyError:
            #print gene, "was not found."
            output_dict[gene] = {"GO:######":("None listed","NA")}
    return output_dict



##### INITIATION & MISC #####

def prime_variables():
    return "/Volumes/Genome/armyant.OGS.V1.8.5.gff", "/Volumes/Genome/Cbir.assembly.v3.0.fa"

def mutate_codon(old_seq, frame, new_nt):
    """Given a Bio.seq.seq sequence, a new nucleotide, and its position in the sequence,
    mutate_codon returns a new Bio.seq.seq sequence with the new nt substituted at the
    given position.
    """
    seq_str = str(old_seq)[0:frame-1] + new_nt + str(old_seq)[frame:]
    new_seq = Seq(seq_str, generic_dna)
    return new_seq

def draw_gene(feature):
    "Produces a visual representation of a gene! (Hopefully!)"
    gd_diagram = GenomeDiagram.Diagram(feature.qualifiers['ID'][0])
    gd_track_for_features = gd_diagram.new_track(1, name="Exons")
    gd_feature_set = gd_track_for_features.new_set()

    #for subfeature in feature.sub_features:
    color = colors.blue
    gd_feature_set.add_feature(feature, color=color, label=True)

    gd_diagram.draw(format='linear', orientation='landscape', pagesize='A4', fragments=1, start=0, end=len(feature))
    gd_diagram.write("test_gene.img.pdf", "PDF")


def test():
    """ The basic commands I have been running of '__main__'
    """

    in_file = "/Volumes/Genome/armyant.OGS.V1.8.5.gff"
    out_file = "/Volumes/Genome/armyant.OGS.V1.8.5.gff"
    in_seq_file = "/Volumes/Genome/Cbir.assembly.v3.0.fa"


    # set scaffold and SNP to check:
    cmdargs = sys.argv
    if len(cmdargs) == 3:
        cmd, scaf, cpos = sys.argv
        pos = int(cpos) - 1
    else:
        scaf = "scaffold61"
        pos = 40905
    print "Position has been converted to python base 0!!!"

    # create SeqRecord file from GFF and Assembly:
    gff_dict = assemble_dict(in_file, in_seq_file)

    # count the number of genes in the SeqRecord dictionary:
    gene_num = count_genes(gff_dict)
    print "there are", gene_num, "gene features in this dictionary."

    #for feature in gff_dict[scaf].features:
    #    print "%s: %s" % (feature.qualifiers['ID'][0], feature.location)

    print "\nResults for SNP %s:%d\n" % (scaf, pos)
    result_dict = snp_in_gene(scaf, pos, gff_dict)

    for feature in result_dict:
        cds_pos, exon, cseq, frame, ref_nt = result_dict[feature]
        prot = cseq.translate()

        print feature.qualifiers['ID'][0],
        try:
            print feature.qualifiers['Name'][0], feature.location
        except KeyError:
            print feature.location
        print "SNP position: %d\texon: %d\t" % (cds_pos, exon)
        print "Codon: %s (%s) Codon Position: %d\t Ref. nt: %s\t" % (cseq, prot, frame, ref_nt)

        new_seq = mutate_codon(cseq, frame, "A")
        print new_seq, new_seq.translate()

def create_polymorphism_files():
    print "Assembling SeqRecord from gff file..."
    gff_dict = assemble_dict()
    snp_dict = {}

    print "Unpickling NSL_PILSIP genotypes..."
    aGtypes, iGtypes = unpickle_gtypes('/Volumes/Genome/RAD-Tags/NSL_PILSIP.RD_15')

    print "Scanning for polymorphisms..."
    for scaf, pos in aGtypes['LINEC']:
        SNP1 = aGtypes['LINEC'][(scaf, pos)][0]
        SNP2 = aGtypes['LINEC'][(scaf, pos)][1]

        # have to remove indels at this point :(
        if SNP1 == "I" or SNP2 == "I" or SNP1 == "J" or SNP2 == "J":
            continue

        results_dict = snp_in_gene(scaf, pos, gff_dict)
        for feature in results_dict:
            cds_pos, exon, ref_cdn, frame, ref_nt = results_dict[feature]
            ref_ptn = ref_cdn.translate()
            snp1_cdn = mutate_codon(ref_cdn, frame, SNP1)
            snp2_cdn = mutate_codon(ref_cdn, frame, SNP2)
            snp1_ptn = snp1_cdn.translate()
            snp2_ptn = snp2_cdn.translate()

            snp_dict[(scaf, pos, feature.qualifiers['ID'][0], len(feature), cds_pos)] = ( (ref_nt, ref_cdn, ref_ptn), (SNP1, snp1_cdn, snp1_ptn), (SNP2, snp2_cdn, snp2_ptn) )

    print len(snp_dict), "results found"
    print "Pickling results..."
    pickle_jar(snp_dict, "/Volumes/Genome/RAD-Tags/LineC_polymorphisms")

    print "writing results to tab delimited file..."
    wobj = open( "/Volumes/Genome/RAD-Tags/LineC_polymorphisms.list", 'w' )
    wobj.write("Scaffold\tPos\tAccession\tCDS_pos (%)\tPolymorphism1\tPolymorphism2\n")
    for scaf, pos, id, cds_len, cds_pos in snp_dict:
        if str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1]) == str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][1][1]) and str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1]) == str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][2][1]):
            pass
        elif str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1]) == str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][1][1]):
            wobj.write( "%s\t%d\t%s\t%d (%.2f%%)\t%s>%s (%s>%s)\n" % (scaf, pos, id, cds_pos, 100.0 * cds_pos / cds_len, snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][2][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][2], snp_dict[(scaf, pos, id, cds_len, cds_pos)][2][2])  )
        elif str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1]) == str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][2][1]):
            wobj.write( "%s\t%d\t%s\t%d (%.2f%%)\t%s>%s (%s>%s)\n" % (scaf, pos, id, cds_pos, 100.0 * cds_pos / cds_len, snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][1][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][2], snp_dict[(scaf, pos, id, cds_len, cds_pos)][1][2])  )
        else:
            wobj.write( "%s\t%d\t%s\t%d (%.2f%%)\t%s>%s (%s>%s)\t%s>%s (%s>%s)\n" % (scaf, pos, id, cds_pos, 100.0 * cds_pos / cds_len, snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][1][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][2], snp_dict[(scaf, pos, id, cds_len, cds_pos)][1][2], snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][2][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][2], snp_dict[(scaf, pos, id, cds_len, cds_pos)][2][2])  )
    wobj.close()

    print "All done!"

def find_next_isoform(gene_id, gff_file="/Volumes/Genome/armyant.OGS.V1.8.5.gff"):
    "returns the next number in the list of isoforms for a given gene"
    isoforms = []

    gff_handle = open(gff_file, 'rb')
    for line in gff_handle:
        gid = re.search(gene_id, line)
        giso = re.search('-R([A-Z]+)', line)
        if gid is not None and giso is not None:
            isoforms.append(giso.groups()[0])

    last_iso = isoforms.sorted()[-1]
    # research how to determine next letter in Alphabet
    return last_iso

##### FASTA FILE FUNCTIONS #####

def count_chromosomes(assembly="/Volumes/Genome/Cbir.assembly.v3.0_gi.fa"):
    "uses fasta file of genome to produce tab delimited list of chromosomes and their size"
    newfile = os.path.splitext(assembly)[0] + ".chrom.list"
    newfile_h = open(newfile, 'w')
    assembly_handle = open(assembly, 'rb')
    csize = 0
    chrom = None
    for line in assembly_handle:
        if re.search('>', line) is not None:
            if chrom is not None:
                newfile_h.write( "%s\t%d\n" % (chrom, csize) )
            chrom = line.split()[0][1:]
            csize = 0
        else:
            csize += len(line.strip())
    newfile_h.write( "%s\t%d\n" % (chrom, csize) )

    newfile_h.close()
    assembly_handle.close()

##### CUFFLINKS ASSOCIATED FUNCTIONS #####

def find_cuffjoined(gtf_file):
    "looks for genes joined together by cuffmerge"
    allid = {}
    file_h = open(gtf_file)
    for line in file_h:
        gid = re.search("gene_id \"([A-Za-z0-9_]*)\";", line)
        nearest = re.search("nearest_ref \"([A-Za-z0-9_\(\)]*)\";", line)
        if gid is not None and nearest is not None:
            if gid.groups()[0] not in allid:
                allid[gid.groups()[0]] = [nearest.groups()[0]]
            elif nearest.groups()[0] not in allid[gid.groups()[0]]:
                allid[gid.groups()[0]].append(nearest.groups()[0])
    file_h.close()

    dblid = {}
    for gid in allid:
        if len(allid[gid]) > 1:
            dblid[gid] = allid[gid]

    return dblid

def separate_cuffjoined(gtf_file, dblid):
    """finds all loci in the gtf_file that are joined according to dblid, and alters the
    gene_id to create new gene_ids for all joined genes (effectively separating them again)
    """
    # set variables:
    newfile =  os.path.splitext(gtf_file)[0] + ".separated.gtf"
    count = 0

    # open files:
    newfile_h = open(newfile, 'w')
    file_h = open(gtf_file)

    for line in file_h:
        gid = re.search("gene_id \"([A-Za-z0-9_]*)\";", line)
        nearest = re.search("nearest_ref \"([A-Za-z0-9_\(\)]*)\";", line)
        if gid.groups()[0] in dblid:
            if nearest is None:
                newid = gid.groups()[0] + "novel"
                count += 1
            else:
                newid = gid.groups()[0] + nearest.groups()[0]
        else:
            newid = gid.groups()[0]
        newfile_h.write( line.replace(gid.groups()[0], newid) )

    file_h.close()
    newfile_h.close()
    print count
    return newfile


##### INTEGRATED PROGRAMS ##################

def investigate(args, doblast=True):
    #gffobj = assemble_dict(in_file=args.gff_file, in_seq_file=args.genome_file)
    print "assembling monster"
    go_monster = genematch.GO_maker(args.GO_file)
    print "assembing my gff"
    mygff = My_gff(args.gff_file)

    genefile_h = open(args.input_file, 'rb')
    genelist = [ line.split()[0] for line in genefile_h]
    genefile_h.close()

    reportfile_h = open(args.output_file, 'w')


    # setup progress bar for this rather long process that is to follow:
    bar = progressbar.ProgressBar(maxval=len(genelist), \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.ETA()]) # can also use progressbar.Percentage()
    count=0

    print "\nAnalysing %d genes:" % (len(genelist))
    bar.update(count)
    for gene in genelist:
        #print "*" * 7, "Processing", gene, "*" * 7
        genegos = go_monster.findem(gene)
        genename = mygff.nameit(gene)
        reportfile_h.write("%s\n%-12s'%s'\n%s\n" % ("*" * 70, gene, genename, "*" * 70 ) )
        reportfile_h.write("".join([ "\t%s %s %s\n" % (g, genegos[g][1], genegos[g][0]) for g in genegos ]))
        reportfile_h.write("-" * 70 + "\n")
        if doblast:
            geneseq = genematch.extractseq(gene)
            geneblast = genematch.blast_ncbi(geneseq, queryterms='(("cerapachys biroi"[Organism]) OR "drosophila melangaster"[Organism]) OR "caenorhabditis elegans"[Organism]')
            results = genematch.blast_results(geneblast,3)
            for alignment,hsp in results:
                try:
                    title = re.search('^[^>]+', alignment.title).group(0)
                except:
                    title = alignment.title
                reportfile_h.write( title + "\n" )
                reportfile_h.write( "Score: %d\tBits: %d\tE-value: %d\n" %
                    (hsp.score, hsp.bits, hsp.expect) )
                reportfile_h.write( "id: %d(%.2f%%)\t+ve: %d(%.2f%%)\n" %
                    (hsp.identities, 100.0 * hsp.identities / alignment.length, hsp.positives,
                    100.0 * hsp.positives / alignment.length) )
        reportfile_h.write("\n\n")
        count+=1
        bar.update(count)
    bar.finish()
    reportfile_h.close()



def methylation_analysis(args):
    print "assembling gffobj..."
    gffobj = assemble_dict(in_file=args.gff_file, in_seq_file=args.genome_file)

    print "assembling intronchecker..."
    intronchecker = My_gff(args.gff_file)
    print "Intron checker gff created:", intronchecker

    print "assembling go monster..."
    go_monster = genematch.GO_maker(args.GO_file)

    out_h = open(args.output_file, 'w')
    out_h.write( "%-20s%-10s%-7s%-10s%-6s%-6s%-12s%-5s%-8s%-5s\n" % ("scaf", "posn","strand","coverage","freqC", "freqT", "gene_id", "exon", "cds_pos", "codon_str") )
    out_h.close()

    # collect details of immensity of task before you:
    cmd = "wc " + args.input_file + " > results.tmp"
    os.system(cmd)
    size = open("result.tmp", 'rb').readline().split()[0]
    print "There are %d lines to be processed (approx %d hours)." % (size, size/1000000)
    os.system("rm result.tmp")
    # setup progress bar
    bar = progressbar.ProgressBar(maxval=size, \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.ETA()])
    count=0


    print "extracting SNP information"
    deg_h = open(args.input_file, 'rb')

    deg_h.next()
    for line in deg_h:
        count += 1
        bar.update(count)
        scaf = line.split()[1]
        posn = int(line.split()[2])
        strand = line.split()[3]
        coverage = line.split()[4]
        freqC = line.split()[5]
        freqT = line.split()[6]


        SNPdict = snp_in_gene(scaf, posn, gffobj)
        # SNP_dict = {"gene_id":None, "gene_name":None, "cds_pos":None, "exon":None, "codon":None, 'codon_str':None, "frame":None, "ref_nt":None, 'ref_nt_str':None}
        if SNPdict["gene_id"] is not None:
            scaf = line.split()[1]
            posn = int(line.split()[2])
            genegos = go_monster.findem(SNPdict["gene_id"])  # was... parse_go(SNPdict["gene_id"])
            # genegos = {"GO:######":("GO function","GO definition")}

            out_h = open(args.output_file, 'a')
            out_h.write( "%-20s %-10d %-4s %-6s %-8s %-8s" % (scaf,posn,strand,coverage,freqC,freqT) \
                + " %(gene_id)-15s %(exon)-4d %(cds_pos)-7d %(codon_str)-5s " % (SNPdict) \
                + "   ".join([ g + " " + genegos[g][1] + " " + genegos[g][0] for g in genegos ]) + "\n" )
            out_h.close()
        else:
            # find out if SNP is in an intron:
            ingene = intronchecker.gene((scaf, posn))
            #print "%r" % (ingene)
            if ingene: # ie, SNP lies on an intron of a gene
                genegos = go_monster.findem(SNPdict["gene_id"])
                out_h = open(args.output_file, 'a')
                out_h.write( "%-20s %-10d %-4s %-6s %-8s %-8s" % (scaf,posn,strand,coverage,freqC,freqT) \
                    + " %-15s %-4d %-7s %-5s " % (ingene, 0, 'n/a', 'n/a' ) \
                    + "   ".join([ g + " " + genegos[g][1] + " " + genegos[g][0] for g in genegos ]) + "\n" )
                out_h.close()
            else:
                out_h = open(args.output_file, 'a')
                out_h.write( "%-20s %-10d %-4s %-6s %-8s %-8s" % (scaf,posn,strand,coverage,freqC,freqT) \
                    + " %-17s %-4d %-7s %-5s \n" % ("Intergenic", -1, 'n/a', 'n/a' ) )
                out_h.close()
    else:
        sys.stdout.write("%d lines processed.\n" % (count) )
        sys.stdout.flush()
    bar.finish()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="Various GFF and gene-file manipulations")
    parser.add_argument("-o", "--output_file", type=str, default="output.list", help="File to save results to")
    parser.add_argument("-i", "--input_file", type=str,  help="File to analyse")
    parser.add_argument("-g", "--gff_file", type=str, default="/Volumes/Genome/armyant.OGS.V1.8.6_lcl.gff", help="GFF file for analyses")
    parser.add_argument("-f", "--genome_file", type=str, default="/Volumes/Genome/Cbir.assembly.v3.0_gi.fa", help="Genome fasta file for analyses")
    parser.add_argument("-G", "--GO_file", type=str, default='/Volumes/Genome/Genome_analysis/Gene_Ontology/armyant.OGS.V1.5.GOterms.list', help="GO file for analyses")
    parser.add_argument("-I", "--investigate", action='store_true',  help="analyse a list of genes")
    parser.add_argument("-b", "--blastoff", action='store_true',  help="turns off blast search for investigate option")
    parser.add_argument("-m", "--methylation", action='store_true',  help="perform methylation analysis")
    parser.add_argument("-S", "--splicing", action='store_true',  help="perform alternate splicing analysis")
    args = parser.parse_args()

    if args.investigate:
        investigate(args, doblast=not(args.blastoff))

    if args.methylation:
        methylation_analysis(args)

    if args.splicing:
        trial = Splicer()
        print trial
        bedfile = '/Volumes/Genome/transcriptomes/BroodSwap/controls/C16/Tophat/tophat_F2_20131210/junctions.bed'
        print "scaffold                 canonical   alt  readj  novel"
        sumalt = 0
        sumreadj = 0
        summatch = 0
        sumnovel = 0
        for scaf in trial.canonical:
            match, alt, readj, novel = trial.map_junctions(bedfile, chosenscaffold=scaf)
            summatch += match
            sumalt += alt
            sumreadj += readj
            sumnovel += novel
            #if alt > 0:
                #print "%-25s:%-4d %-4d %-4d %-4d %s" % (scaf, match, alt, readj, novel, sorted(trial.alternatives[scaf].values()))
            #else:
            print "%-25s:%-4d %-4d %-4d %-4d" % (scaf, match, alt, readj, novel)
        grandsum = float(summatch + sumalt + sumreadj + sumnovel)/100
        print ("Total canonical:  %d (%.2f%%)\n\
                Total alternate:  %d (%.2f%%)\n\
                Total readjusted: %d (%.2f%%)\n\
                Total novel:      %d (%.2f%%)" % (summatch, summatch/grandsum, sumalt, sumalt/grandsum, sumreadj, sumreadj/grandsum, sumnovel, sumnovel/grandsum))



