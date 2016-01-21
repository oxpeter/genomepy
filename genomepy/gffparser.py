#!/usr/bin/env python
""" A module for parsing GFF files """

import sys
import cPickle
import re
import os
import datetime
import random

import argparse
from BCBio import GFF
from Bio import SeqIO
from Bio.Alphabet import generic_dna, IUPAC
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
#from reportlab.lib import colors
#from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram

import gff2bed2
from genomepy import genematch, config
from genomepy.config import pickle_jar, open_pickle_jar, pickle_gtypes, unpickle_gtypes





####### CLASSES ############################################################

class Splicer(object):
    """a class that allows analysis of splice junctions from a gff file and tophat
    junction files"""
    def __init__(self, gff_file="", chosenscaffold='All'):
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

class FastaLibrary(object):
    def __init__(self, fastafile, buildtime=False):
        # initialise variables
        t0 = datetime.datetime.now()
        qc = 0
        self.seqlib ={}
        self.fastafile = fastafile
        geneseq = ""

        # read genome fasta file
        handle = open(fastafile, 'rb')
        for line in handle:
            if line[0] == '>':
                try:
                    self.seqlib[query.group(2)] = geneseq
                except UnboundLocalError:
                    pass
                query = re.search( '>(gi\|[0-9]+\|ref\|)?([\w\.]+)\|?', line)
                geneseq = ""
            else:
                geneseq += line.strip()
        else:
            self.seqlib[query.group(2)] = geneseq
        handle.close()

        # report time (if requested)
        if buildtime:
            t1 = datetime.datetime.now()
            verbalise("B", "Time to build library: %s" % datetime.timedelta.total_seconds(t1-t0))

    def __getitem__(self, key):
        return self.seqlib[key]

    def __contains__(self, key):
        return key in self.seqlib

    def __repr__(self):
        return "[FastaLibrary object] %r" % (self.fastafile)

    def __str__(self):
        return "[FASTA LIBRARY]\nFrom file %s\n%d sequences" % (self.fastafile, len(self.seqlib))

class GffFeature(object):
    "a sub-element of a gff library"
    def __init__(self, line="" ):
        self.grandparent = None
        self.parent = None
        self.children = {}
        self.grandchildren = {}

        self.master = None  # currently, this is only for gene-level features

        self.atts = parse_atts(line)
        self.flds = parse_cols(line)
        # flds:  "scaf","source","type","start","end","score","strand", "phase"

        if self.flds:       # mark whether details were actually present to form feature
            self.empty = False
        else:
            self.empty = True

        # assign a unique id:
        # NB: cds features do not have unique ids in the NCBI attribute field. Needs to be
        # combined with start position to become unique.
        if not self.empty:
            if self.flds['type'] == 'CDS':
                suffix = "_%s" % str(self.flds['start'])
            else:
                suffix = ""

            try:
                id = self.atts['ID']
            except KeyError:
                id = random.randint(1000000,9999999)
            self.id = "%s%s" % (id, suffix)

    def __repr__(self):
        return "%r %r %r" % (self.id, self.atts, self.flds)

    def __str__(self):
        try:
            name = self.atts['Name']
        except KeyError:
            name = "---"
        try:
            gene = self.atts['gene']
        except KeyError:
            gene = "---"

        if self.parent:
            parent = self.parent
        else:
            parent = "---"
        try:
            ftype = self.flds['type']
        except KeyError:
            ftype = "---"

        return "id:%s name:%s gene:%s parent:%s type:%s " % (self.id, name, gene, parent, ftype )

class GffLibrary(object):
    def __init__(self, gff_file, fasta_file=None):
        # information on data source
        self.gff_file = gff_file
        self.load_assembly(fasta_file)

        # dictionaries to access subfeatures
        self.featlib = {}   # lists partitioned according to feature type
        self.scaflib = {}   # lists partitioned according to scaffold then feature type
                            # (includes an 'all' category that contains everything)

        self.namelib  = {}  # access subfeatures by subfeature name (LOC123, XP_1234, XM_1234)
        self.idlib    = {}  # access subfeatures by subfeature id
        self.locuslib = {}  # FUTURE: for possible faster searching by position?

        handle = open(gff_file, 'rb')
        for line in handle:
            # create subfeature instance
            subfeature = GffFeature(line)
            if subfeature.empty:
                continue

            scaf = subfeature.flds['scaf']
            ftype = subfeature.flds['type']

            # assign subfeature to all Gff libraries
            # featlib
            if ftype in self.featlib:
                self.featlib[ftype].append(subfeature)
            else:
                self.featlib[ftype] = [subfeature]

            # scaflib
            if scaf in self.scaflib:
                if ftype in self.scaflib[scaf]:
                    self.scaflib[scaf][ftype].append(subfeature)
                    self.scaflib[scaf]['all'].append(subfeature)
                else:
                    self.scaflib[scaf][ftype] = [subfeature]
                    self.scaflib[scaf]['all'].append(subfeature) # 'all' category already created
            else:
                self.scaflib[scaf] = {'all':[subfeature], ftype:[subfeature]}


            # idlib
            self.idlib[subfeature.id] = subfeature

            # namelib
            if 'Name' in subfeature.atts:
                if  subfeature.atts['Name'] in self.namelib:
                    self.namelib[subfeature.atts['Name']].append(subfeature)
                else:
                    self.namelib[subfeature.atts['Name']] = [subfeature]

            # find parent and children subfeatures and connect them
            if 'Parent' in subfeature.atts:
                subfeature.parent = self.idlib[subfeature.atts['Parent']]
                self.idlib[subfeature.atts['Parent']].children[subfeature.id] = subfeature

                if subfeature.parent.parent:
                    subfeature.grandparent = subfeature.parent.parent
                    subfeature.grandparent.grandchildren[subfeature.id] = subfeature

    def __repr__(self):
        return "%d features from library %r" % (len(self.idlib), self.gff_file)

    def __str__(self):
        has_ass = str(self.fastalib) if self.fastalib else "Fasta Library: None"
        ret_str = "[GFF LIBRARY]\nFrom file %s\n%s\n%s\nTotal features: %d"
        return  ret_str % (os.path.basename(self.gff_file),
                            has_ass,
                            "\n".join([ f + ": " + str(len(self.featlib[f])) for f in self.featlib ]),
                            len(self.idlib))

    def __contains__(self, item):
        if item in self.namelib or item in self.idlib:
            return True
        else:
            return False

    def __len__(self):
        return len(self.idlib)

    def build_master_gene(self):
        """
        For each gene, extract all exons, merge overlapping exons, construct framework of
        maximum bounds for each exon, assign exon number (5' to 3'). Add to master_gene
        dic. This method does not alter any existing attributes of the gff instance
        (unless the master_gene attribute had already been created).
        """
        mastercount = 0
        for gene in self.featlib['gene']:

            exons = []
            for gc in gene.grandchildren.values():
                if gc.flds['type'] == 'exon':
                    exons.append((gc.flds['start'], gc.flds['end'], gc.flds['strand']))
            exonlist = list(set(exons))

            if len(exonlist) == 0:
                continue

            # create master overlapping exons:
            while len(exonlist[0]) == 3:
                newlists = [],[] # 1st list = non-overlapping exons, 2nd list = overlapping
                for x in exonlist:
                    newlists[overlaps(x, exonlist[0])].append(x)
                all_positions = [xpos for tup in newlists[1] for xpos in tup[:2] ]
                exonlist = newlists[0] + [(min(*all_positions), max(*all_positions))]

            # FUTURE: rather than store as a list, create GffFeature instances for each
            # master exon and store them instead.
            gene.master = sorted([ x[:2] for x in exonlist ],
                                        key=lambda i: i[0],
                                        reverse=(True if gene.flds['strand']=='-' else False))
            mastercount += 1
        return mastercount

    def findfeatures(self, locus, strand=['+','-'], ftype='all'):
        """looks to see if locus is within features. Returns dictionary of all types
        that contain the locus with lists of those feature instances.

        this method replaces the old ingene and inexon methods, providing greater
        specificity and flexibility.

        OUTPUT from  ftype='all' (NB: actual feature instance is appended, not just its id):
        found['gene'] = [gene01]
        found['mRNA'] = [rna01, rna02]
        found['exon'] = [exon01, exon07]
        found['CDS']  = [cds07]
        """

        scaf, posn = locus
        found = {}
        try:
            for feat in self.scaflib[scaf][ftype]:
                start = min(feat.flds['start'], feat.flds['end'])
                end = max(feat.flds['start'], feat.flds['end'])

                if feat.flds['strand'] in strand:
                    if start <= posn <= end:
                        if feat.flds['type'] in found:
                            found[feat.flds['type']].append(feat)
                        else:
                            found[feat.flds['type']] = [feat]
        except KeyError:
            pass

        return found

    def extractseq(self, featurename, buffer=0, cds=False, translate=False,
                    boundaries=False, trim_from=0, trim_to=999999999):
        """
        For the given feature, return the genomic or cds DNA sequence. Buffer will return
        x bp up- and down-stream of the gene, but only in genomic DNA (not for cds option).
        If boundaries is True, it will return a list of cumulative relative positions of
        exon boundaries for the cds sequence instead (and will force cds to be True)
        """
        # check that the sequence library was supplied:
        if not self.fastalib:
            return {">No genome assembly was supplied to the gff library" % featurename:"NNNN" }

        if boundaries:
            cds = True
            translate = False

        features = self._get_cdfeatures(featurename, cds)
        seq_dic = {}

        if not features:
            return {">%s not found in gff library" % featurename:"NNNN" }

        elif not cds:
            for i,f in enumerate(features):
                seq, defline = self._get_seq(f, buffer)
                if seq and translate:
                    seq_dic[defline] = seq[trim_from:trim_to].translate()
                elif seq:
                    seq_dic[defline] = seq[trim_from:trim_to]
                else:
                    seq_dic[">%s seq%d %s not found" % (featurename, i, f)] = None

        else:
            for feat in features:
                for isoform in feat:
                    cdsseq = []
                    # pull out some useful identifiers from first grandchild to create defline:
                    defline = ">%s" % self.idlib[isoform[0].atts['Parent']]
                    for att in ['gene', 'protein_id', 'transcript_id', 'product']:
                        if att in self.idlib[isoform[0].atts['Parent']].atts:
                            " ".join([defline, self.idlib[isoform[0].atts['Parent']].atts[att]])

                    # collect each sequence:
                    for child in isoform:
                        if child.flds['type'] == 'CDS':
                            # save start location for sorting exons
                            cdsseq.append( (child.flds['start'],self._get_seq(child)[0]) )

                    reverse = True if child.flds['strand'] == '-' else False

                    if boundaries:
                        seq = [ len(p[1]) for p in sorted(cdsseq,
                                                            key=lambda x: x[0],
                                                            reverse=reverse) ]
                    else:
                        seq = reduce(lambda x,y: x+y,
                                        [ p[1] for p in sorted(cdsseq,
                                        key=lambda x: x[0],
                                        reverse=reverse)])[trim_from:trim_to]
                    if translate:
                        seq = seq.translate()
                    seq_dic[defline] = seq

        return seq_dic


    def _get_cdfeatures(self, featurename, cds=False):
        # get feature from requested name (returns a list):
        if featurename in self.idlib:
            feature = [self.idlib[featurename]]
        elif featurename in self.namelib:
            feature = self.namelib[featurename]
        else:
            return None

        f_list = []
        if cds:
            for i,f in enumerate(feature):
                # if cds requested, iterate through the children of the gene to paste
                # their sequences together.
                if cds and f.flds['type'] == 'mRNA':
                    cdfeatures = [ [ elem for elem in f.children.values() if elem.flds['type'] == 'CDS'],]
                elif cds and f.flds['type'] == 'gene':
                    cdfeatures = [ [elem for elem in c.children.values() if elem.flds['type'] == 'CDS'] for c in f.children.values() ]
                else:
                    cdfeatures = None
                f_list.append(cdfeatures)
        else:
            f_list = feature

        return f_list

    def _get_seq(self, f, buffer=0):
        if f.flds['scaf'] in self.fastalib:
            start = min(f.flds['start'], f.flds['end']) - buffer - 1 # to make python count
            end   = max(f.flds['start'], f.flds['end']) + buffer
            if start < 0:
                start = 0
            seq   = Seq(self.fastalib[f.flds['scaf']][start:end ])

            # check if reverse complement needs to be determined:
            if f.flds['strand']=='-':
                seq_rc = Seq(seq, generic_dna).reverse_complement()
                seq    = seq_rc

            defline = ">%s [%s, %d:%d]" % (f, f.flds['scaf'], start, end)
        else:
            seq = Seq("")
            defline = None
        return seq, defline

    def fetch_promoter(self, featurename, upstream=450, downstream=50):
        """
        For the given feature, return the sequence corresonding to the promoter region
        specified
        """
        # get feature from requested name (returns a list):
        if featurename in self.idlib:
            feature = [self.idlib[featurename]]
        elif featurename in self.namelib:
            feature = self.namelib[featurename]
        else:
            return {">%s" % featurename:"%s not found in gff library" % featurename}

        seq_list = {}
        for i,f in enumerate(feature):
            # the strand determines which position will be the start, and which direction
            # is upstream and downstream, and whether the reverse complement needs to
            # be extracted:
            if f.flds['scaf'] in self.fastalib:
                if f.flds['strand']=='+':
                    start = min(f.flds['start'], f.flds['end']) - upstream - 1
                    end   = min(f.flds['start'], f.flds['end']) + downstream
                    if start < 0:
                        start = 0
                    seq   = self.fastalib[f.flds['scaf']][start:end]
                else:
                    start  = max(f.flds['start'], f.flds['end']) - downstream - 1
                    end    = max(f.flds['start'], f.flds['end']) + upstream
                    if start < 0:
                        start = 0
                    seq_rc = Seq(self.fastalib[f.flds['scaf']][start:end], generic_dna)
                    seq    = seq_rc.reverse_complement()

                seq_list[">%s [%s, %d:%d](searched:%s result:%d)" % (f, f.flds['scaf'], start, end, featurename, i)] = seq
            else:
                seq_list[">%s seq%d %s" % (featurename, i, f)] = "%s not found in fasta library" % (f.flds['scaf'])
        return seq_list

    def findnearest(self, locus, ftype='all'):
        """
        Finds the closest features of type ftype, both up- and downstream of the given
        locus.
        """
        # first check if locus is within a feature:
        is_within = self.findfeatures(locus, ftype=ftype)

        if len(is_within) > 0:
            feat_list = [ f for eachlist in is_within.values() for f in eachlist ]
            return (0, feat_list), (0, feat_list)

        else:
            closestup   = 999999999
            closestdown = 999999999
            upgene   = []
            downgene = []

            # next check if scaffold contains any features:
            if locus[0] not in self.scaflib:
                return ((closestup, []), (closestdown, []))

            # now find closest of all features on this scaffold:
            for feature in self.scaflib[locus[0]][ftype]:
                if 0 < locus[1] - max(feature.flds['start'],feature.flds['end']) < closestup:
                    closestup = locus[1] - max(feature.flds['start'],feature.flds['end'])
                    upgene = [feature]
                elif locus[1] - max(feature.flds['start'],feature.flds['end']) == closestup:
                    upgene.append(feature)

                if 0 < min(feature.flds['start'],feature.flds['end']) - locus[1] < closestdown:
                    closestdown = min(feature.flds['start'],feature.flds['end']) - locus[1]
                    downgene = [feature]
                elif min(feature.flds['start'],feature.flds['end']) - locus[1] == closestdown:
                    downgene.append(feature)
            upstream = (closestup, upgene)
            downstream = (closestdown, downgene)

            return upstream, downstream

    def load_assembly(self, fasta_file):
        if fasta_file:
            self.fastalib = FastaLibrary(fasta_file, False) # genome seq lib for sequence extraction
        else:
            self.fastalib = None

####### DEPRECATED CLASS (MAINTAINED CURRENTLY FOR BACKWARDS COMPATABILITY) ############

class My_gff(object):
    "an object for quick assessment of where a GENE/SNP lies"
    def __init__(self, gff_file="", primary_key='ID', highest='gene'):

        self.gff_origin = gff_file
        self.primary_key = primary_key
        self.highest_level = highest

        self.genedict = {}  # a dictionary of scaffolds, containing a dictionary of genes
        self.featurecount = 0  # number of features stored in class instance
        self.genenames = {} # dictionary of useful feature names
        self.genescaf = {}  # dictionary to find which scaf each feature is on
        self.exondict = {}  # dictionary to list all exons for a feature
        self.strandinfo = {} # dict to find which strand a feature is on
        self.toplevel = {}  # top level dictionary (typically mRNA is the feature, so
        #                     this dic allows clustering of mRNAs into their corr. genes
        self.geneid = {}    # in case ID is not the primary identifier, to find parents
        self.mrna2gene = {} # allow identification of gene parent of an mRNA feature
        self.id2gene = {}
        self.id2rna = {}
        self.master = {}    # for storing master gene exon positions (ordered)

        gff_h = open(gff_file,'rb')
        error_log = {}
        for line in gff_h:
            attr = parse_atts(line)
            fields = parse_cols(line)
            if not fields:
                continue
            start = fields['start']
            stop = fields['end']
            try:
                pkey = attr[primary_key]
            except KeyError:
                if fields['type'] in error_log:
                    error_log[fields['type']] += 1
                else:
                    error_log[fields['type']] = 1
                continue

            if fields['type'] == 'gene':
                self.toplevel[attr['Name']] = {}
                for title, attribute in [("mRNA", []),
                                        ('Longname', attr['gene']),
                                        ('locus', attr['Name']),
                                        ('ID', attr['ID']),
                                        ('scaffold', fields['scaf']),
                                        ]:

                    try:
                        self.toplevel[attr['Name']][title] = attribute

                    except KeyError:
                        continue

            if fields['type'] == 'mRNA':
                self.featurecount += 1
                self.geneid[attr['ID']] = pkey

                # add mRNA to list of mRNAs for given gene:
                if 'gene' in attr:
                    if attr['gene'] in self.toplevel:
                        if self.toplevel[attr['gene']]['ID'] == attr['Parent']:
                            self.toplevel[attr['gene']]['mRNA'].append(pkey) ## attr[pkey]
                            self.mrna2gene[pkey] = attr['gene']
                if fields['strand'] == "+":
                    self.strandinfo[pkey] = 1
                else:
                    self.strandinfo[pkey] = 0
                self.genescaf[pkey] = fields['scaf']

                # find most appropriate field for gene name:
                for f in ["Product", "Name", "Dbxref", "ID"]:
                    if f in attr:
                        self.genenames[pkey] = attr[f]
                        break
                else:
                    self.genenames[pkey] = "--------"

                scaf = fields['scaf']
                if scaf in self.genedict:
                    self.genedict[scaf][pkey] = (start, stop)
                else:
                    self.genedict[scaf] = {pkey:(start, stop)}

            if fields['type'] == 'CDS':
                parentid = self.geneid[attr['Parent']]
                if parentid in self.exondict:
                    self.exondict[parentid].append((start, stop, fields['phase']))
                else:
                    self.exondict[parentid] = [(start, stop, fields['phase'])]


            ############## /\ /\ OLD VERSION /\ /\ #############################
            ############## \/ \/ NEW VERSION \/ \/ #############################


        if len(error_log) > 0:
            print "%d feature%s could not be incorporated into gff:\n%s" % (len(error_log),
                    's' if len(error_log)>1 else '',
                    "\n".join([ i + ": " + str(error_log[i]) + " errors" for i in error_log ]))

    def __contains__(self, locus):
        "looks to see if locus is within gene. Does not consider introns/exons"
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
        gffname = os.path.basename(self.gff_origin).split('.')[0][:10]
        return "%-10s: %d scaffolds, %d genes" % (gffname,
                            len(self.genedict),self.featurecount)

    def __repr__(self):
        return "My_gff: %r %r %r" % (self.gff_origin, self.primary_key, self.highest_level)

    def __len__(self):
        return (self.featurecount)

    def build_master_gene(self):
        """
        For each gene, extract all exons, merge overlapping exons, construct framework of
        maximum bounds for each exon, assign exon number (5' to 3'). Add to master_gene
        dic. This method does not alter any existing attributes of the gff instance
        (unless the master_gene attribute had already been created).
        """
        # get all exons:
        for gene in self.toplevel:

            exons = []
            for mrna in self.toplevel[gene]['mRNA']:
                exons += self.exondict[mrna]
            exonlist = list(set(exons))

            if len(exonlist) == 0:
                continue

            # create master overlapping exons:
            while len(exonlist[0]) == 3:
                newlists = [],[] # 1st list = non-overlapping exons, 2nd list = overlapping
                for x in exonlist:
                    newlists[overlaps(x, exonlist[0])].append(x)
                all_positions = [xpos for tup in newlists[1] for xpos in tup[:2] ]
                exonlist = newlists[0] + [(min(*all_positions), max(*all_positions))]

            self.master[gene] = sorted([ x[:2] for x in exonlist ],
                                        key=lambda i: i[0],
                                        reverse=(not self.strandinfo[mrna]))

    def whichscaf(self, pkey, exact=True):
        try:
            scaf = self.genescaf[pkey]
        except KeyError:
            if exact:
                scaf = None
            else:
                allgenes = ( gene for scaf in self.genedict for gene in self.genedict[scaf] )
                matches = []
                for gene in allgenes:
                    if re.search(pkey, gene):
                        matches.append(gene)
                if len(matches) == 1:
                    scaf = self.genescaf[matches[0]]
                else:
                    strand = None

        return scaf

    def whichstrand(self, pkey, exact=True):
        try:
            strand = self.strandinfo[pkey]
        except KeyError:
            if exact:
                strand = None
            else:
                allgenes = ( gene for scaf in self.genedict for gene in self.genedict[scaf] )
                matches = []
                for gene in allgenes:
                    if re.search(pkey, gene):
                        matches.append(gene)
                if len(matches) == 1:
                    strand = self.strandinfo[matches[0]]
                else:
                    strand = None
        return strand

    def gene(self, locus):
        scaf, posn = locus
        ingene = False
        try:
            for pkey in self.genedict[scaf]:
                if self.genedict[scaf][pkey][0] <= posn <= self.genedict[scaf][pkey][1]:
                    ingene = pkey
        except KeyError:
            ingene = None

        return ingene

    def ingene(self, locus):
        """looks to see if locus is within gene. Does not consider introns/exons.
        If true, returns the gene name"""
        scaf, posn = locus
        ingene = []
        try:
            for gene in self.genedict[scaf]:
                if self.genedict[scaf][gene][0] <= posn <= self.genedict[scaf][gene][1]:
                    ingene.append(gene)
        except KeyError:
            ingene = []

        return ingene

    def inexon(self, locus, posn):
        """checks to see if position is in an exonic region.
        Can either give a scaf + posn (gene=False), and it will find if it's in a gene,
        and if it is, if it's in an exon. Or you can give a posn + gene (set gene=True),
        and it will quickly check if it's in the exon of that gene"""
        if locus in self.geneid:
            pkey = [locus]
        else:
            pkey = self.ingene((locus,posn))
            if len(pkey) == 0:
                return False

        # now we have the pkey...
        inexon = { i:None for i in pkey }
        for g in pkey:
            for exon in self.exondict[g]:
                if min(exon) <= posn <= max(exon):
                    inexon[g] = True
            else:
                inexon[g] = False
        return inexon

    def closest_splice(self, locus, posn):
        closest_distance = None
        pkey = self.ingene((locus, posn))
        if self.inexon(posn, locus):
            closest_distance = 999999999999
            for start, end in self.exondict[pkey]:
                if abs(posn - start) < closest_distance:
                    closest_distance =  abs(posn - start)
                if abs(end - posn) < closest_distance:
                    closest_distance =  abs(end - posn)
        return closest_distance

    def whichtranscript(self, locus):
        pass

    def findnearest(self, scaffold, hitpos):
        """ looks for the nearest gene to a given locus """

        pkey = self.ingene((scaffold, hitpos))
        if pkey:
            upstream = (0, pkey)
            downstream = (0, pkey)

        else:
            closestup = 999999999
            closestdown = 999999999
            upgene = 'No gene upstream on scaffold'
            downgene = 'No gene downstream on scaffold'
            if scaffold not in self.genedict:
                return ((closestup, upgene), (closestdown, downgene))
            for gene in self.genedict[scaffold]: # check distances between locus and genes
                if 0 < hitpos - self.genedict[scaffold][gene][0] < closestup:
                    closestup = hitpos - self.genedict[scaffold][gene][0]
                    upgene = gene
                if 0 < self.genedict[scaffold][gene][0] - hitpos < closestdown:
                    closestdown = self.genedict[scaffold][gene][0] - hitpos
                    downgene = gene
            upstream = (closestup, upgene)
            downstream = (closestdown, downgene)

        return upstream, downstream

    def nameit(self, pkey):
        try:
            return self.genenames[pkey]
        except KeyError:
            return None

####### FUNCTIONS ######################################################################

def define_parameters():
    parser = argparse.ArgumentParser(description="Various GFF and gene-file manipulations")
    ## output options
    parser.add_argument("-d", "--directory", type=str,
                        help="Specify the directory to save results to")
    parser.add_argument("-o", "--output_file", type=str, default="output.list",
                        help="Name of file to save results to")
    parser.add_argument("-q", "--quiet", action='store_true', default=False,
                        help="Print fewer messages and output details")

    ## input options
    parser.add_argument("-i", "--input_file", type=str,
                        help="File to analyse")
    parser.add_argument("-g", "--gff", type=str, default=dbpaths['gff'],
                        help="GFF file for analyses")
    parser.add_argument("-f", "--genome_file", type=str, default=dbpaths['ass'],
                        help="Genome fasta file")
    parser.add_argument("-G", "--GO_file", type=str, default=dbpaths['goterms'],
                        help=".list file containing GO terms for each gene")

    ## analysis options
    parser.add_argument("--show_go", action='store_true', default=False,
                        help="run GO analyses")
    parser.add_argument("-I", "--investigate", action='store_true',
                        help="Analyse a list of genes")
    parser.add_argument("-b", "--blastoff", action='store_true',
                        help="Turns off blast search for investigate option")
    parser.add_argument("-M", "--methylation", action='store_true',
                        help="Perform methylation analysis")
    parser.add_argument("-S", "--splicing", type=str,
                        help="Perform alternate splicing analysis on given file")
    parser.add_argument("--bed2gtf", type=str,
                        help="Converts from bed to both ex.gff and gtf formats")
    parser.add_argument("--gff2bed", action='store_true',
                        help="Converts from gff to bed format")
    parser.add_argument("-t", "--trim", type=str,
                        help="Trims UTRs from gff file")
    parser.add_argument("-s", "--strip", type=str,
                        help="Removes duplicate transcripts from bed file")
    parser.add_argument("--upstream", type=int, default=450,
                        help="Set a value for upstream boundary")
    parser.add_argument("--downstream", type=int, default=50,
                        help="Set a value for downstream boundary")
    parser.add_argument("--promoters", type=str, default="",
                        help="""Find promoter sequence for specified gene. If "all"
                        is specified then all genes in gff will be searched.
                        """)
    parser.add_argument('-c', "--column", type=int, default=0,
                        help="Set column number for extracting list values from. (Default=0)")

    return parser

##### GFF FILE PARSING  ############

def parse_atts(line):
    """
    Parses a single line from a gff file.
    The function takes the 9th column (containing all the feature attributes) and returns
    a dictionary of { attribute_name : attribute_value }
    """
    if len(line) == 0:
        return None
    if line[0] == '#':
        return None
    cols = line.split()
    attr = { pair.split('=')[0]:pair.split('=')[1] for pair in " ".join(cols[8:]).split(';')  }
    return attr

def parse_cols(line):
    """
    Takes a gff line and extracts the information from the first 8 columns, returning
    them in a dictionary.
    """
    if len(line) == 0:
        return None
    if line[0] == '#':
        return None
    cols = line.split()
    try:
        attrs = { "scaf":cols[0],       "source":cols[1],   "type":cols[2],
                  "start":int(cols[3]),  "end":int(cols[4]),
                  "score":cols[5],  "strand":cols[6],   "phase":cols[7]
                }
    except ValueError:
        cols = line.split('\t')
        try:
            attrs = { "scaf":cols[0],       "source":cols[1],   "type":cols[2],
                  "start":int(cols[3]),  "end":int(cols[4]),
                  "score":cols[5],  "strand":cols[6],   "phase":cols[7]
                }
        except ValueError:
            verbalise("R", "Error parsing line:\n%s" % line)
    return attrs


##### GFF FILE MANIPULATION ########

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

def make_bed(gff_file, assembly_file):
    """ Takes a gff_file and creates a bed format including the scaffolds/contigs that
    do not have any features associated with them"""

    gff_fh = open(gff_file, 'rU')

    # get transcript annotation
    tinfo, einfo = gff2bed2.ParseAnno(gff_fh)
    # write into bed format
    scaf_dict = gff2bed2.WriteBED(tinfo, einfo, return_obj=True, gff_name=gff_file)

    # create dictionary of scaffolds not annotated, and their sizes
    missing_dict = {}
    assembly_handle = open(assembly_file, 'r')
    for line in assembly_handle:
        def_search = re.search('>(lcl\|)?([a-zC0-9]*)', line)
        line_len = len(line)
        if def_search:
            scaf_id = def_search.group(1)
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
    genelist = {}
    scaffoldlist = {}

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

        #try:
        #    gene_iso    = re.search("(\w+)\.([A-Za-z0-9_\(\)\.]*)", labels[0])
        #    unexpressed = re.search("ID=(Cbir[A-Za-z0-9_\(\)\.]*)", labels[0]) # genes that are not expressed won't have the cufflinks id
        #except IndexError:
        #    print fields[3], labels

        #if unexpressed is not None:
        #    geneid      = unexpressed.group(1)
        #    isoform     = geneid
        #else:
        #    try:
        #        geneid      = gene_iso.group(1)
        #        isoform     = gene_iso.group(2)
        #    except:
        #        print line
        #        print labels[0]
        #        die
        geneid = labels[1]
        isoform = labels[0][3:]
        scaffoldlist[geneid] = scaf
        # update gene information:
        #try:
        #    genenames[geneid] += isoform + ">"
        #except KeyError:
        #    genenames[geneid] = isoform + ">"

        if geneid in genesizes:
            genesizes[geneid] += [gene_start, gene_end]

            #if min(gene_start, gene_end) < min(genesizes[geneid]):
            #    genesizes[geneid][genesizes[geneid].index(min(genesizes[geneid]))] = min(gene_start, gene_end)
            #if max(gene_start, gene_end) > max(genesizes[geneid]):
            #    genesizes[geneid][genesizes[geneid].index(max(genesizes[geneid]))] = max(gene_start, gene_end)
        else:
            genesizes[geneid] = [gene_start, gene_end]

        genestrands[geneid] = strand

        # get exon info
        exon_starts = [int(x) for x in fields[11].split(",")]
        exon_lengths= [int(x) for x in fields[10].split(",")]
        exon_number = 1

        # create gtf and gff files:
        transcript_line = "\t".join([scaf, "gffparser", "transcript", str(gene_start), str(gene_end), '1', strand, ".", 'gene_id "' + geneid + '"; transcript_id "' + isoform + '";\n'])
        gene_line = '\t'.join([scaf, "gffparser", "gene", str(min(genesizes[geneid])), str(max(genesizes[geneid])), '.', genestrands[geneid], '.', 'ID=' + geneid + ';Name=' + geneid + '\n'])
        mRNA_line = "\t".join([scaf, "gffparser", "mRNA", str(gene_start), str(gene_end), ".", strand, ".", "ID=m_" + isoform + ";Parent=" + geneid + "\n"])

        gtf_h.write(transcript_line)
        gff_h.write(mRNA_line)

        for exon_start, exon_length in zip(exon_starts, exon_lengths):
            gtf_line = '\t'.join([scaf, "gffparser", "exon", str(exon_start + gene_start), str(exon_start + gene_start + exon_length - 1), '1', strand, '.', 'gene_id "' + geneid + '"; transcript_id "' + isoform + '"; exon_number "' + str(exon_number) + '";\n'])
            gff_CDS  = '\t'.join([scaf, "gffparser", "CDS",  str(exon_start + gene_start), str(exon_start + gene_start + exon_length - 1), '.', strand, '.', 'ID=cds.' + str(exon_number) + isoform + ';Parent=m_' + isoform + '\n'])
            gtf_h.write(gtf_line)
            gff_h.write(gff_CDS)
            exon_number += 1


    #convert gene name to Cbir version if available:
    for geneid in genesizes:
        gene_line = '\t'.join([scaffoldlist[geneid], "gffparser", "gene", str(min(genesizes[geneid])), str(max(genesizes[geneid])), '.', genestrands[geneid], '.', 'ID=' + geneid + ';Name=' + geneid + '\n'])
        gff_h.write(gene_line)

    gff_h.close()
    gtf_h.close()

def highest_cbir(file):
    "Finds the current highest number assigned to Cbir gene ids"
    pep_h = open(file, 'rb')

    highest = 1
    for line in pep_h:
        num = re.search("(Cbir[A-Za-z0-9_\(\)\.\/]*)", line)
        if num is not None:
            if int(num.group(1)) > highest:
                highest = int(num.group(1))

    return highest

def assign_cbir(gtffile):
    """Takes cuffcompare gtf file and extracts XLOC id along with Cbir id.

    """

    gtf_h = open(gtffile, 'rb')
    locus = {}
    details = {}
    new_cbir = highest_cbir() + 1
    for line in gtf_h:
        fields = line.split('\t')
        detail_list = fields[8].strip().split(';')
        try:
            details = dict([(x.strip().split(' ')[0], x.strip().split(' ')[1]) for x in detail_list if x != ''])
        except IndexError:
            print detail_list
            continue
        try:
            geneid = re.search("(Cbir[A-Za-z0-9_\(\)\.\/]*)", details['nearest_ref'])
        except KeyError:
            geneid = re.search("(Cbir[A-Za-z0-9_\(\)\.\/]*)", details['oId'])

        if geneid is not None and details['class_code'][1:-1] in ("=","c","j"):
            locus[details['gene_id'][1:-1]] = geneid.group(1)
        elif details['class_code'][1:-1] in ("=","c","j", "o"):
            locus[details['gene_id'][1:-1]] = details['gene_id']
        elif details['class_code'][1:-1] in ("x","i","u"):
            locus[details['gene_id'][1:-1]] = 'Cbir_' + str(new_cbir)
            new_cbir += 1
        else:
            print details['class_code'],details['gene_id'],geneid.group(1),'###',
    return locus

def rename_loci(gtffile, bedfile):
    "transforms bedfile XLOC gene ids, replacing them with existing or new Cbir ids."
    locus = assign_cbir(gtffile)
    oldbed = open(bedfile, 'rb')
    newbed = open(bedfile[:-3] + 'cbir.bed', 'w')

    for line in oldbed:
        oldid = re.search(';(XLOC_[0-9]*)[^;]*;', line)
        if oldid is not None:
            try:
                newid = locus[oldid.group(1)]
            except KeyError:
                newid = oldid.group(0)[1:-1]
            oldname = oldid.group(0)
        else:
            newid = ""
            oldname = "notfoundatall"
        newline = line.replace(oldname, ';' + newid + ';')
        newbed.write(newline)

    oldbed.close()
    newbed.close()

def strip_duplicates(bedfile):
    "removes duplicate gene entries, and replacing ID with the cufflinks appended Cbir ID"
    bed_h = open(bedfile, 'rb')
    newbed_h = open(bedfile[:-3] + "uniq.bed", 'w')
    reject_h = open(bedfile[:-3] + "reject.bed", 'w')
    genedic = {}

    for line in bed_h:
        fields = line.split()
        # check that gene is not already annotated:
        if (tuple(fields[0:3]),tuple(fields[5:])) in genedic:
            reject_h.write(line)
        else:
            newbed_h.write(line)
            genedic[(tuple(fields[0:3]),tuple(fields[5:]))] = True

    bed_h.close()
    newbed_h.close()
    reject_h.close()

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

def assemble_dict(in_file, in_seq_file, features_only=False):
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
            continue
        try:
            geneid = re.search("(Cbir[A-Za-z0-9_\(\)\.\/]*)", details['nearest_ref'])
        except KeyError:
            geneid = re.search("(Cbir[A-Za-z0-9_\(\)\.\/]*)", details['oId'])

        if geneid is not None and details['class_code'] == "=":
            details['transcript_id'] = '"%s.%s"' % (details['gene_id'].strip('"'),geneid.group(1))
        else:
            details['transcript_id'] = '"%s.%s"' % (details['gene_id'].strip('"'),details['transcript_id'].strip('"'))

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

    NB: there can be multiple genes sharing the same locus, (usually by being on
    different strands), so each gene needs to be tested (which will slow performance,
    sadly).
    """

    SNP_dict = {"gene_id":None, "gene_name":None, "cds_pos":None,
                    "codon":None, 'codon_str':None, "frame":None,
                    "ref_nt":None, 'ref_nt_str':None, "exon":None}

    for feat in gff_dict[scaf].features:
        if pos in feat:

            # print gene details for match:
            #print feat.qualifiers['ID'][0],

            if 'pseudo' in feat.qualifiers:
                continue

            if 'partial' in feat.qualifiers:
                continue
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

            try:
                cds_pos = exon_pos + before_len
            except UnboundLocalError:
                """
                verbalise("R", "ERROR for ", gene_name, gene_id)
                verbalise("M", before_len, before_cnt)
                verbalise("Y", pos)
                verbalise("G", feat.location)
                verbalise("R", dir(feat))

                verbalise("G", [ str(sf.location) for sf in feat.sub_features ])
                verbalise("R", feat.sub_features)
                """
                return SNP_dict
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
            SNP_dict['gene_id']=(gene_id)
            SNP_dict['gene_name']=(gene_name)
            SNP_dict['cds_pos']=(cds_pos)
            SNP_dict['exon']=(exon)
            SNP_dict['codon']=(codon)
            SNP_dict['codon_str']=(str(codon))
            SNP_dict['frame']=(frame)
            SNP_dict['ref_nt']=(ref_nt)
            SNP_dict['ref_nt_str']=(str(ref_nt))
            break
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

def parse_go(gene, gofile):
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

def overlaps(x1, x2):
    "returns true if the exons overlap in position. Does not check if on same scaffold!"
    x1_1 = min(x1[:2])
    x1_2 = max(x1[:2])
    x2_1 = min(x2[:2])
    x2_2 = max(x2[:2])

    return  x2_1 <= x1_1 <= x2_2 or x1_1 <= x2_1 <= x1_2

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

def test(dbpaths):
    """ The basic commands I have been running of '__main__'
    """

    in_file =  dbpaths['gff']
    out_file = dbpaths['gff']
    in_seq_file = dbpaths['ass']


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

def create_polymorphism_files(gtypefile=os.path.realpath('./NSL_PILSIP.RD_15')):
    print "Assembling SeqRecord from gff file..."
    gff_dict = assemble_dict()
    snp_dict = {}

    print "Unpickling NSL_PILSIP genotypes..."
    aGtypes, iGtypes = unpickle_gtypes(gtypefile)

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
    pickle_jar(snp_dict, os.path.realpath("./LineC_polymorphisms"))

    print "writing results to tab delimited file..."
    wobj = open( os.path.realpath("./LineC_polymorphisms.list"), 'w' )
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

def find_next_isoform(gene_id, gff_file):
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

def count_chromosomes(assembly):
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

    genelist = config.make_a_list(args.input_file)

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
    verbalise("Y", "Intron checker gff created:", intronchecker)

    print "assembling go monster..."
    go_monster = genematch.GO_maker(args.GO_file)

    out_h = open(args.output_file, 'w')
    out_h.write( "%-20s%-10s%-7s%-10s%-6s%-6s%-12s%-5s%-8s%-5s\n" % ("scaf", "posn","strand","coverage","freqC", "freqT", "gene_id", "exon", "cds_pos", "codon_str") )
    out_h.close()

    # collect details of immensity of task before you:
    cmd = "wc " + args.input_file

    wc = os.popen(cmd)
    data = wc.readline()
    wc.close()
    size = int( data.split()[0] )
    verbalise("Y",
        "There are %d lines to be processed (approx %d hours)." % (size, size/1000000)
            )


    count=0
    t0 = datetime.datetime.now()

    print "extracting SNP information"
    deg_h = open(args.input_file, 'rb')

    deg_h.next()
    for line in deg_h:
        count += 1
        if count % 10000 == 0:
            t1 = datetime.datetime.now()
            tdiff = datetime.timedelta.total_seconds(t1-t0)
            speed = 1.0 * count / tdiff  # (counts/s)
            remaining = size - count
            tremaining = datetime.timedelta(seconds=remaining / speed)
            print "\r%d/%d (%.2f%%) complete. Time remaining ~ %s            " % (count,size, 100.0*count/size, tremaining),
            sys.stdout.flush()
        scaf = line.split()[1]
        posn = int(line.split()[2])
        #strand = line.split()[3]
        #coverage = line.split()[4]
        #freqC = line.split()[5]
        #freqT = line.split()[6]
        strand=""
        coverage=""
        freqC=""
        freqT=""

        SNPdict = snp_in_gene(scaf, posn, gffobj)
        # SNP_dict = {"gene_id":None, "gene_name":None, "cds_pos":None,
        #             "exon":None,    "codon":None,     "codon_str":None,
        #             "frame":None,   "ref_nt":None,    "ref_nt_str":None}

        if SNPdict["gene_name"]:
            scaf = line.split()[1]
            posn = int(line.split()[2])
            if args.show_go:
                genegos = go_monster.findem(SNPdict["gene_name"])
                # This was... parse_go(SNPdict["gene_id"])
                # genegos = {"GO:######":("GO function","GO definition")}

            out_h = open(args.output_file, 'a')
            snp_details ="%-20s %-10d %-4s %-6s %-8s %-8s" % (scaf,posn,strand,coverage,freqC,freqT)
            gene_details= " %(gene_name)-15s %(exon)-4d %(cds_pos)-7d %(codon_str)-5s " % (SNPdict)
            if args.show_go:
                go_details="   ".join([g + " " + genegos[g][1] + " " + genegos[g][0] for g in genegos])
            else:
                go_details=""

            out_h.write(snp_details +  gene_details + go_details + "\n")
            out_h.close()

        else:
            # find out if SNP is in an intron:
            ingene = intronchecker.gene((scaf, posn))
            if ingene: # ie, SNP lies on an intron of a gene
                if args.show_go:
                    genegos = go_monster.findem(SNPdict["gene_name"])
                genename = intronchecker.nameit(ingene)
                snp_details="%-20s %-10d %-4s %-6s %-8s %-8s" % (scaf,posn,strand,coverage,freqC,freqT)
                gene_details=" %-15s %-4d %-7s %-5s " % (genename, 0, 'n/a', 'n/a' )
                if args.show_go:
                    go_details="   ".join([g+" "+genegos[g][1]+" "+genegos[g][0] for g in genegos])
                else:
                    go_details=""

                out_h = open(args.output_file, 'a')
                out_h.write(snp_details +  gene_details + go_details + "\n")
                out_h.close()

            else:
                snp_details="%-20s %-10d %-4s %-6s %-8s %-8s" % (scaf,posn,strand,coverage,freqC,freqT)
                gene_details=" %-17s %-4d %-7s %-5s" % ("Intergenic", -1, 'n/a', 'n/a' )

                out_h = open(args.output_file, 'a')
                out_h.write( snp_details + gene_details + "\n" )
                out_h.close()

    else:
        sys.stdout.write("%d lines processed.\n" % (count) )
        sys.stdout.flush()

if __name__ == '__main__':

    dbpaths = config.import_paths()
    parser = define_parameters()

    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output_file)

    # allow use of config file presets to find reference files:
    if args.gff in dbpaths.keys():
        args.gff = dbpaths[args.gff]

    if args.genome_file in dbpaths.keys():
        args.genome_file = dbpaths[args.genome_file]

    if args.investigate:
        investigate(args, doblast=not(args.blastoff))

    if args.methylation:
        methylation_analysis(args)

    if args.splicing:
        trial = Splicer()
        print trial
        bedfile = args.splicing
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

    if args.gff2bed:
        make_bed(args.input, args.genome_file)

    if args.trim:
        trim_untranslated(args.trim)

    if args.bed2gtf:
        bed2gtf(args.bed2gtf)

    if args.strip:
        strip_duplicates(args.strip)

    if args.promoters:
        # build gff library
        verbalise("B", "Building gene library...")
        gfflib = GffLibrary(args.gff, args.genome_file)
        #verbalise("G", gfflib.fastalib)
        verbalise("G", gfflib)
        verbalise("B", "Extracting promoters...")
        # extract promoters
        all_promoters = {}
        if args.promoters == 'all':
            for gene in ( g.id for g in gfflib.featlib['gene']):
                all_promoters[gene] = gfflib.fetch_promoter(gene,
                                                            upstream=args.upstream,
                                                            downstream=args.downstream)
        else:
            for gene in config.make_a_list(args.promoters, args.column):
                all_promoters[gene] = gfflib.fetch_promoter(gene,
                                                            upstream=args.upstream,
                                                            downstream=args.downstream)

        # print to file:
        outname = logfile[:-3] + 'promoters.fasta'
        errorlog = logfile[:-3] + 'errors.log'
        handle = open(outname, 'w')
        errhandle = open(errorlog, 'a')
        errcount = 0
        for p in all_promoters:
            for defline in all_promoters[p]:
                if str(all_promoters[p][defline])[-7:] == 'library':
                    errhandle.write('%s\n%s\n' % (defline, str(all_promoters[p][defline])))
                    errcount += 1
                else:
                    handle.write('%s\n%s\n' % (defline, str(all_promoters[p][defline])))
        handle.close()
        errhandle.close()
        if errcount > 0:
            verbalise("R", "%d genes could not have their sequence extracted" % errcount)




