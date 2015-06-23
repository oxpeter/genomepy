#!/usr/bin/python

""" A series of functions allowing interaction between local gff/cds/pep/fa files and
Tuxedo suite result files and ncbi's BLAST engines
"""
import os
import sys
import re
import argparse
from pkg_resources import resource_filename, resource_exists


import numpy as np
import Bio.Blast.NCBIWWW as ncbi
import Bio.Blast.NCBIXML as xml
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, IUPAC
from scipy.stats import fisher_exact
from scipy import interpolate
import matplotlib.pyplot as plt

from genomepy import config

###### INITIALISE THE FILE PATHS NEEDED FOR ANALYSIS #######################

dbpaths = config.import_paths()

############################################
class GO_item(object):
    "holds the relevant data about a GO term"
    def __init__(self, id, name, namespace, alt_id=[], definition="", subsets=[], synonym=""):

        self.id = self.process_goid(id)

        self.name = name

        # check namespace is appropriate to one of the three categories:
        if namespace.upper() in ["MOLECULAR_FUNCTION","MOLECULAR FUNCTION", "MF"]:
            self.namespace = "MF"
        elif namespace.upper() in ["BIOLOGICAL_PROCESS","BIOLOGICAL PROCESS", "BP"]:
            self.namespace = "BP"
        elif namespace.upper() in ["CELLULAR_COMPONENT","CELLULAR COMPONENT", "CL"]:
            self.namespace = "CL"

        self.altids = [self.process_goid(a) for a in alt_id]

        self.definition = definition
        self.subsets = subsets

        if isinstance(synonym, list):
            self.synonym = synonym
        else:
            self.synonym = [synonym]

    def __str__(self):
        return "GO:%s [%s] '%s'" % (self.id, self.namespace, self.name)

    __repr__ = __str__

    def process_goid(self, id):
        if id[:3] == "GO:":
            try:
                b = int(id[3:])
            except ValueError:
                verbalise("R", "GO term must be numeric! ('GO:' prefix allowed)")
            return id[3:]

        else:
            try:
                b = int(id)
            except ValueError:
                verbalise("R", "GO term must be numeric! ('GO:' prefix allowed)")
            return id

class GO_maker(object):
    "an object for quick access to all GO terms for a gene set"
    def __init__(self, gofile=dbpaths['goterms']):

        self.go_repo = {}
        self.go_dict = {}
        self.go_defs = {}

        go_handle = open(gofile, 'rb')
        for line in go_handle:
            gene = line.split()[0]
            self.go_dict[gene] = {}
            patt = "GO:([0-9]+);?(\s(?!GO:)[\w; ]+)?"
            go_pairs = re.findall(patt, line)
            for pair in go_pairs:
                goterm = pair[0] # ie the GO number GO:001432
                defpattern = " ([-A-Za-z ]*)\; ([A-Za-z ]*)"
                defgo = re.search(defpattern, pair[1])  # the type and definition of the go term
                if len(goterm)>0 and defgo:
                    # parse GO term into type (eg Molecular Process)
                    # and definition (eg Redox Reaction)
                    gotype = defgo.group(2).split()[0][0] + defgo.group(2).split()[1][0]
                    godef =  defgo.group(1)
                    self.go_dict[gene][goterm] =  (godef, gotype)
                    if goterm in self.go_repo:
                        self.go_repo[goterm].append(gene)
                    else:
                        self.go_repo[goterm] = [gene]
                    self.go_defs[goterm] = (godef, gotype)
                elif len(goterm)>0:
                    self.go_dict[gene][goterm] =  ("-unavailable-", "-unavailable-")
                    if goterm in self.go_repo:
                        self.go_repo[goterm].append(gene)
                    else:
                        self.go_repo[goterm] = [gene]
                    self.go_defs[goterm] = ("-unavailable-", "-unavailable-")

    def count_goterms(self):
        "returns the number of unique go terms stored in object"
        return len(self.go_repo)

    def count_genes(self):
        "returns the number of unique gene ids stored"
        return len(self.go_dict)

    def fetch_genes(self, goterm):
        "given a goterm, returns the genes that have that go term"
        try:
            return self.go_repo[goterm]
        except KeyError:
            #print "No genes with %s found!" % (goterm)
            return ['NoGeneFound']

    def fetch_gos(self, genelist):
        "yields subsequent dictionaries of go terms associated with each geneid in genelist"
        for geneid in genelist:
            try:
                yield self.go_dict[geneid]
            except KeyError:
                yield {"GO:######":("None listed","None listed")}

    def findem(self, geneid):
        "given a gene, return the dictionary of GO terms associated with it"
        try:
            output_dict = self.go_dict[geneid]
        except KeyError:
            #print gene, "was not found."
            output_dict = {"GO:######":("None listed","None listed")}

        return output_dict

    def define_go(self, goterm):
        'given a GO number, return the information about it'
        try:
            return self.go_defs[goterm]
        except KeyError:
            return ("GO not found","GO not found")

class Fisher_square(object):
    def __init__(self, TwG, TwoG, NwG, NwoG, id, id_info=None):
        self.TwG  = int(TwG)  # Number of Target Genes with GO terms
        self.TwoG = int(TwoG) # Number of Target Genes without GO terms
        self.NwG  = int(NwG)  # Number of Non-Target genes with GO terms
        self.NwoG = int(NwoG) # Number of Non-Target genes without GO terms

        self.id = id
        self.id_info = id_info  # with idea that this contains GO_item instance

        self.p = None
        self.odds = None

    def __str__(self):
        # create label:
        if self.id_info:
            header = "%s\n" % (self.id_info)
        else:
            header = "%s\n" % (self.id)

        # create results lines:
        if self.p and self.odds:
            results = "%20s %.3f\n%20s %.3f\n" %   ("odds-ratio:", self.odds,
                                "P-value:", self.p)
        else:
            results = "No tests performed\n"

        # create table of values:
        body = "%-22s %s\n%20s %-7d %d\n%20s %-7d %d\n" % ("Has GO", "Doesn't Have GO",
                                "Target Set:", self.TwG, self.TwoG,
                                "Non-Target Set:", self.NwG, self.NwoG,
                                )

        fullstring = header + body + results
        return fullstring

    __repr__ = __str__

    def fisher_test(self, alt='greater'):
        try:
            self.odds, self.p = fisher_exact([
                [self.TwG, self.TwoG],
                [self.NwG, self.NwoG]
            ], alternative=alt)
        except ValueError:
            verbalise("R", "Error performing Fishers Exact:", self)

############################################
def define_arguments():
    parser = argparse.ArgumentParser(description=
            "A module to perform a variety of gene term related analyses")

    ### input options ###
    # logging options:
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="print fewer messages and output details")
    parser.add_argument("-o", "--output", type=str, default='genematch.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")

    # data file options:
    parser.add_argument("-i", "--input", type=str,
                        help="file to analyse")
    parser.add_argument("-c", "--column", type=int, default=0,
                        help="column to extract data from")
    parser.add_argument("-f", "--datafile", type=str,
                        help="file containing gene terms")
    parser.add_argument("-g", "--gofile", type=str,
                        help=".obo file containing go terms,definitions etc")

    # analysis options:
    parser.add_argument("-E", "--enrichment", action='store_true',default=False,
                        help="perform GO term enrichment on file")
    parser.add_argument("-T", "--min_terms", type=int, default=2,
                        help="the minimum number of terms to report in enrichment results")
    parser.add_argument("-Q", "--max_q", type=int, default=0.05,
                        help="the maximum q_value to report in enrichment results")


    return parser

def set_go_item(attributes):
    if set(['id', 'name', 'namespace', 'def']) <= set(attributes.keys()):
        go_item = GO_item(  attributes['id'],
                            attributes['name'],
                            attributes['namespace'],
                            attributes['alt_id'],
                            attributes['def'],
                            attributes['subset'],
                            attributes['synonym'],
                        )
    else:
        verbalise("R", "CANNOT BUILD GO ITEM. INSUFFICIENT ATTRIBUTES:\n", " ".join([str(p) for p in attributes.keys()]))
        return None
    return go_item

def reset_attributes():
    return {'alt_id':[], 'subset':[],
            'synonym':[], 'subsetdef':[],
            'is_a':[], 'consider':[],
            'xref':[], 'intersection_of':[],
            'relationship':[], 'replaced_by':[],
            'disjoint_from':[], 'subsetdef':[],
            }

def build_go_terms(file):
    """
    This function takes the given OBO file from the Gene Ontology Consortium and
    parses it into a series of dictionaries containing GO_item instances.
    """
    id_dic = {}

    handle = open(file, 'rb')
    attributes = reset_attributes()
    for line in handle:
        posn = line.find(":")
        if posn > 0:
            if line[:posn] in attributes:
                try:
                    attributes[line[:posn]].append(line.strip()[posn+2:])
                except AttributeError:
                    verbalise("R", "CANNOT APPEND TERM FROM LINE:", line)
                    verbalise("R", "\n".join([str(p) for p in attributes.items()]))
                    exit()
            else:
                attributes[line[:posn]] = line.strip()[posn+2:]
        elif line[:2] == "[T":
            go_item = set_go_item(attributes)
            if go_item:
                id_dic[go_item.id] = go_item
            attributes = reset_attributes()
    else:
        go_item = set_go_item(attributes)
        if go_item:
            id_dic[go_item.id] = go_item

    handle.close()
    return id_dic

def extract_locus(fname, col_num, datacol=False):
    """ given a file name and a column number, extract_locus() will create three lists:
    scaffold,
    locus start pos
    locus end pos

    if datacol is provided, the contents of that column will be added to the extra list:
    gdata

    designed for use with cufflinks output files.
    """
    scaf = []
    gbeg = []
    gend = []
    gdata = []

    fobj = open(fname, 'rb')
    for line in fobj:
        col = line.split()
        try:
            scaf.append( re.search('[sCcafold]*[0-9]+', col[col_num]).group() )
            gbeg.append( int(re.search(':(.*)-', col[col_num]).groups()[0]) )
            gend.append( int(re.search('-(.*)', col[col_num]).groups()[0]) )
            if datacol:
                gdata.append( col[datacol] )
        except:
            try:
                scaf.append( re.search('[sCcafold]*[0-9]+', col[col_num - 2]).group() )
                gbeg.append( int(re.search(':(.*)-', col[col_num - 2]).groups()[0]) )
                gend.append( int(re.search('-(.*)', col[col_num - 2]).groups()[0]) )
                if datacol:
                    gdata.append( col[datacol] )
            except AttributeError:
                print "### failed to parse:"
                print line

    fobj.close()
    return scaf, gbeg, gend, gdata

def snp2gene(scaffold, pos, gff=dbpaths['gff']):
    """ given a scaffold and nt position, snp2gene will determine if the nt is within a
    gene, and if so, whether it is within an exon or intron"""

    geneid = 'intergenic'
    geneloc = 'non-coding'

    fobj = open(gff, 'rb')
    for line in fobj:
        col = line.split()
        if col[0] == scaffold:
            if col[2] == "mRNA":
                if int(col[3])<=int(col[4]):
                    if float(col[3]) <= pos <= float(col[4]):
                        geneid = re.search('ID=([^;]*);', col[8]).groups()[0]
                else:
                    if float(col[4]) <= pos <= float(col[3]):
                        geneid = re.search('ID=([^;]*);', col[8]).groups()[0]

            if col[2] == "CDS":
                if float(col[3]) <= pos <= float(col[4]):
                    geneloc = 'coding (exonic)'

    fobj.close()

    return (geneid, geneloc)

def locus2gene(scaflist, gbeglist, gendlist, gdatalist=False, gff=dbpaths['gff'], comprehensive=True ):
    """ For each locus as defined by the three input lists, locus2gene() will pull out the
    overlapping genes from the OGS gff file. If an additional list of corresponding data
    is provided, it will be appended to the results dictionary.
    """
    cuffgenes = {}

    for result in range(len(scaflist)):
        if result % 1000 == 0:
            print "%d genes matched of %d" % (result, len(scaflist))
        cur_scaf = scaflist[result]
        cur_gbeg = gbeglist[result]
        cur_gend = gendlist[result]
        if gdatalist:
            cur_gdata = gdatalist[result]
        else:
            cur_gdata = 0
        fobj = open(gff, 'rb')
        for line in fobj:
            col = line.split()
            if col[2] == "mRNA":
                if col[0] == cur_scaf:
                    if float(col[3]) <= cur_gend and float(col[4]) >= cur_gbeg:
                        try:
                            if (cur_scaf, cur_gbeg) in cuffgenes:
                                cuffgenes[(cur_scaf, cur_gbeg, 2)] = (re.search('ID=([^;]*);', col[8]).groups()[0], cur_scaf, cur_gbeg, cur_gend, cur_gdata)
                            else:
                                cuffgenes[(cur_scaf, cur_gbeg)] = (re.search('ID=([^;]*);', col[8]).groups()[0], cur_scaf, cur_gbeg, cur_gend, cur_gdata)
                                if not comprehensive:
                                    break
                        except AttributeError:
                            print col[8]
        fobj.close()

    return cuffgenes

def findgene(fname, dbpaths=dbpaths):
    """ Takes the cuffdiff results file (from fname), and compares it with the current
    C.biroi gff file, returning as a dictionary any genes that overlap the cuffdiff result.

    OUTPUT:
        dict[scaffold, isoform_start_pos] = GeneID, scaffold, start_pos, end_pos, first_state_FPKM, second_state_FPKM, log2_diff
    """
    scaf = []
    gbeg = []
    gend = []
    gfor = []
    gsta = []
    gdif = []
    cuffgenes = {}

    fobj = open(fname)
    for line in fobj:
        col = line.split()
        scaf.append( re.search('[sCcafold]*[0-9]+', col[3]).group() )
        gbeg.append( int(re.search(':(.*)-', col[3]).groups()[0]) )
        gend.append( int(re.search('-(.*)', col[3]).groups()[0]) )
        gfor.append(float(col[7]))
        gsta.append(float(col[8]))
        gdif.append(float(col[9]))

    fobj.close()
    print "Significant transcripts read"


    for result in range(len(scaf)):
        cur_scaf = scaf[result]
        cur_gbeg = gbeg[result]
        cur_gend = gend[result]
        cur_gfor = gfor[result]
        cur_gsta = gsta[result]
        cur_gdif = gdif[result]
        fobj = open(dbpaths['gff'])
        for line in fobj:
            col = line.split()
            if col[2] == "mRNA":
                if col[0] == cur_scaf:
                    if float(col[3]) <= cur_gend and float(col[4]) >= cur_gbeg:
                        try:
                            cuffgenes[(cur_scaf, cur_gbeg)] = (re.search('ID=([^;]*);', col[8]).groups()[0], cur_scaf, cur_gbeg, cur_gend, cur_gfor, cur_gsta, cur_gdif)
                        except AttributeError:
                            print col[8]
        fobj.close()

    return cuffgenes

def blast_ncbi(geneseq, blasttype='blastp', db='nr', queryterms='("formicidae"[Organism]) OR ("drosophila"[Organism]) OR ("caenorhabditis elegans"[Organism])'):
    """ performs a blast search on the NCBI online server for the sequence geneseq.
    Also (currently) parses result and gives top ten non-hypothetical results """

    return ncbi.qblast(blasttype, db, geneseq, expect=2, hitlist_size=10, entrez_query=queryterms)

def blast_results(blast_results, num_results=10):
    parsed_results = xml.parse(blast_results)
    p_result = parsed_results.next()
    counter = num_results
    #print "Number of alignments = ", len(p_result.alignments)
    for alignment in p_result.alignments:
        for hsp in alignment.hsps:
            yield (alignment,hsp)
            """
            if re.search('hypothetical protein', alignment.title) == None:
                print alignment.title
                print "Score: %d\tBits: %d\tE-value: %d" % (hsp.score, hsp.bits, hsp.expect)
                print "id: %d(%.2f%%)\t+ve: %d(%.2f%%)" % (hsp.identities, 100.0 * hsp.identities / alignment.length, hsp.positives, 100.0 * hsp.positives / alignment.length)
            """
            counter -= 1
        if counter == 0:
            break

def extractseq(geneID, type='pep', OGS=dbpaths['cds'][:-4], startpos=0, endpos=-1):
    """ extracts sequence of geneID from the current annotations. type is cds,  pep or fasta.
    """
    #fname = dbpaths[type]
    geneseq = ""
    fname = ".".join([OGS, type])
    fobj = open(fname)
    for line in fobj:
        if line[0] == '>':
            query = re.search( geneID + '[\s\n]', line)
        if query:
            thisline = fobj.next()

            while thisline[0] != '>':
                geneseq += thisline.strip()
                try:
                    thisline = fobj.next()
                except StopIteration:
                    break
            else:
                break
    fobj.close()
    return geneseq[startpos:endpos]

def fly_orthologs(genelist, dbpaths=dbpaths):
    fobj = open(dbpaths['ortho'])
    orthodict = {}
    for line in fobj:
        if re.search("(Cbir[A-Za-z0-9_\(\)\.\/]*)", line) != None and re.search('DROME', line) != None:
            cbgene = re.search("(Cbir[A-Za-z0-9_\(\)\.\/]*)", line).group()
            if  cbgene in genelist:
                orthodict[cbgene] = re.search('DME\|(.*)_DROME', line).groups()[0]
    fobj.close()
    return orthodict

def gff_init():
    """ To speed up gene matching, this will create a gff dictionary """
    pass

def extract_promoter():
    """
    # To extract Vg promoters from all ants, I used the following under __main__:
    """
    filename = "/Users/POxley/Documents/Analysis/Data/Genes/Vg/Vg_promoters.info"
    file_h = open(filename, 'rb')
    results_dict = {}
    file_h.next()
    for line in file_h:
        # get all the parameters from the file:
        gene_id  = line.split()[0]
        spp      = line.split()[1]
        scaf     = line.split()[2]
        polarity = line.split()[3]
        startpos = int(line.split()[4])
        endpos   = int(line.split()[5])

        # to search the concatenated genome file, need to add spp to scaffold name:
        defcheck = spp + "_" + scaf
        sequence = extractseq(defcheck, type='fa', OGS='/Volumes/Genome/Ant_Genomes/EightAntGenomes', startpos=startpos, endpos=endpos)
        #defcheck = spp + "_" + scaf
        #sequence = extractseq(defcheck, type='fa', OGS='/Volumes/Genome/Non-ant_genomes/Hymen.N.vit.genome', startpos=startpos, endpos=endpos)

        # if polarity is -ve, create reverse complement:
        if polarity == '-':
            seq_obj = Seq(sequence, generic_dna)
            sequence = seq_obj.reverse_complement()
            polarity = "REVERSE COMPLEMENT"

        # display results
        results_dict[gene_id] =  ">%s %s %s %s %d-%d\n%s" % (gene_id, spp, scaf, polarity, startpos, endpos, sequence)
        print results_dict[gene_id]

    #for gene in results_dict:
    #    print results_dict[gene]
    return results_dict

######## Gene Ontology and KEGG pathway analyses ######################
def go_finder(genelist):
    pass

def go_enrichment(genelist, background_go=None, go_info=None, return_odds=False, gofile=dbpaths['goterms']):
    """
    returns f_squares = [fisher_square1, fisher_square2...]
    """
    # if background list of GO terms is not supplied, make one:
    if background_go:
        go_obj = background_go
    else:
        go_obj = GO_maker(gofile)

    f_squares = []  # where the Fisher's Exact Tables will be held (and returned)

    god = {}
    goodds = {}
    gopval = {}

    for godict in go_obj.fetch_gos(genelist):
        for goterm in godict:
            try:
                god[goterm] += 1
            except KeyError:
                god[goterm] = 1

    for goterm in god:
        if goterm == "GO:######":
            continue
        dip = god[goterm]

        go_w_deg = dip
        go_wo_deg = len(go_obj.fetch_genes(goterm)) - dip
        nogo_w_deg = len(genelist) - dip
        nogo_wo_deg = go_obj.count_genes() - len(go_obj.fetch_genes(goterm)) - len(genelist) + dip

        if isinstance(go_info,dict):
            f_squares.append(Fisher_square(go_w_deg, nogo_w_deg, go_wo_deg, nogo_wo_deg, goterm, go_info[goterm]))
        else:
            f_squares.append(Fisher_square(go_w_deg, nogo_w_deg, go_wo_deg, nogo_wo_deg, goterm))

    verbalise("B", "Performing Fisher's Exact Test on %d squares" % (len(f_squares)))
    for square in f_squares:
        square.fisher_test()

    return f_squares

def cbir_to_kegg(genelist, dbpaths=dbpaths, reversedic=False):
    "Converts Cbir gene IDs to Kegg orthologs (KOs)"
    kegg_h = open(dbpaths['keggortho'], 'rb')

    keggd = {}
    keggcount = {}

    for line in kegg_h:
        if len(line.split()) == 2: # ie, the line has a cbir entry and KO
            keggd[line.split()[0]] = line.split()[1]
            keggcount[line.split()[1]] = True
        else:
            #keggd[line.split()[0]] = 'na'
            pass

    keggdic = {}

    if reversedic:  # create {geneid:ko}
        for gene in genelist:
            try:
                keggdic[gene] = keggd[gene]
            except KeyError:
                keggdic[gene] = None

    else:       # create {ko:geneid}
        for gene in genelist:
            try:
                keggdic[keggd[gene]] = gene
            except KeyError:
                pass
    return len(keggd), keggdic

def ko_to_pathway(ko, dbpaths=dbpaths):
    "given a kegg ortholog (ko), returns the kegg pathway it's part of"
    kmod_h = open(dbpaths['kegg'], 'rb')

    koterm = ''
    kdescription = 'no description found'
    for line in kmod_h:
        if line[0] == 'C':     # Find details of Kegg Pathway
            try:
                ksearch = re.search("C +([0-9]*) *(.*)\[PATH", line)
                pathwayid = ksearch.group(1)
                pathwaydef = ksearch.group(2)
            except:
                pathwayid = 'cannot parse'
                pathwaydef = 'cannot parse'

        elif line[0] == 'D':    # Find details of Kegg term
            ksearch = re.search("(K[0-9]*) *(.*) \[?", line)
            try:
                koterm = ksearch.group(1)
                kodesc = ksearch.group(2)
            except:
                koterm = 'cannot parse'
                kodesc = 'cannot parse'

        if ko == koterm:
            kdescription = '%s \n(in %s pathway)' % (kodesc, pathwaydef)

    return kdescription

def cbir_to_pathway(geneobj):
    "given a gene (or list of genes), returns the kegg pathways associated"

    gene_ko = {}
    if type(geneobj) is list or type(geneobj) is dict:
        dictlen, kodic = cbir_to_kegg(geneobj, reversedic=True)
        for gene in geneobj:
            gene_ko[gene] = ko_to_pathway(kodic[gene])
    else:   # assume it is a single gene:
        dictlen, kodic = cbir_to_kegg([geneobj], reversedic=True)
        gene_ko[geneobj] = ko_to_pathway(kodic[geneobj])
    return gene_ko

def kegg_module_enrichment(genelist, dbpaths=dbpaths):
    """
    For  a given genelist, determines which KEGG modules are significantly
    enriched. Returns both module enrichment P-values, and upper-level classification
    P-values (as dictionaries).
    NB - will not likely show many, as modules are very small. Better to use BRITE pathways
    instead (via kegg_pathway_enrichment() ).
    """
    size, kegglist = cbir_to_kegg(genelist)

    #print "Creating Kegg module library"
    ## create Kegg module dictionary:
    kmod_h = open(dbpaths['kegg'], 'rb')

    kmodlist = []
    kmodd = {}
    kmodcount = 0
    kmcount = {}
    kmod = 'none'

    kgroupb = {}
    kgroupbcount = {}
    kgroupblist = []
    kbcount = {}
    b_group = 'none'

    kgroupc = {}
    kgroupccount = {}
    kccount = {}
    kgroupclist = []
    c_group = 'none'

    kod = {}    # this is to count how many KOs there are in the list

    for line in kmod_h:
        """if line[0] == 'B':  #   higher level functional description (eg 'Energy Metabolism')
            kbcount[b_group] = len(kgroupbcount)# add old b_group count to dictionary
            kgroupbcount = {}                   # reset counter
            try:
                b_group = re.search("<b>(.*)</b>", line).group(1)
                kgroupblist.append(b_group)
            except:
                b_group = 'none'
        """
        if line[0] == 'C':    # descriptive function (eg Carbon fixation)
            kccount[c_group] = len(kgroupccount)# add old c_group count to dictionary
            kgroupccount = {}                   # reset counter
            try:
                c_group = re.search("C *(.*)", line.trim()).group(1)
                kgroupclist.append(c_group)
            except:
                c_group = 'none'

        elif line[0] == 'D':     # Kegg module (eg M00165 Reductive pentose phosphate cycle)
            kmcount[kmod] = kmodcount
            kmodcount = 0
            try:
                ksearch = re.search("(M[0-9]*) *(.*)\[PATH", line)
                kmod = ksearch.group(1)
                kmoddef = ksearch.group(2)
                kmodlist.append(kmod)
            except:
                kmod = 'none'
                kmoddef = 'none'

        elif line[0] == 'E':    # Kegg term
            try:
                ko = re.search("(K[0-9]*)", line).group(1)

                kmodcount += 1
                kgroupbcount[ko] = True
                kgroupccount[ko] = True
                try:
                    kmodd[ko].append(kmod)
                    kgroupb[ko].append(b_group)
                    kgroupc[ko].append(c_group)
                except:
                    kmodd[ko] = [kmod]
                    kgroupb[ko] = [b_group]
                    kgroupc[ko] = [c_group]
            except:
                pass
    kccount[c_group] = len(kgroupccount)# add final old c_group count to dictionary
    kmcount[kmod] = kmodcount           # add final kmmod group count to dictionary


    kmodtestsize = len(kmodlist)    # for multiple testing correction
    kgrpbtestsize = len(kgroupblist)
    kgrpctestsize = len(kgroupclist)

    #print "calculating Fisher's exact test"
    # count number of kegglist KOs are in each kegg module, perform Fisher's exact test
    dip = 0
    dipcount = {}               # possibly unnecessary
    kmodenrich = {}
    kmododds   = {}

    for mod in kmodlist:
        for ko in kegglist:
            if ko not in kmodd: # some KOs do not exist in a module.
                pass
            elif mod in kmodd[ko]:
                dip += 1
        dipcount[mod] = dip     # possibly unnecessary
        oddrat, pval = fisher_exact([
            [dip, len(kegglist) - dip],
            [kmcount[mod]-dip, len(keggcount) - len(kegglist) - kmcount[mod] + dip]
        ])
        if pval < 0.05:
            print "%s\n         In Path  Not in Path\nDEG    :  %-7d %d\nnon-DEG:  %-7d %d\n%.4f\n" % (mod, dip, len(kegglist) - dip,kmcount[mod]-dip, len(keggcount) - len(kegglist) - kmcount[mod] + dip, pval )
        kmododds[mod]   = oddrat
        kmodenrich[mod] = pval
        dip = 0     # reset for next module
    #sys.stdout.write("%s\n         In Path  Not in Path\nDEG    :  %-7d %d\nnon-DEG:  %-7d %d\n%.4f\n" % (mod, dip, len(kegglist) - dip,kmcount[mod]-dip, len(keggcount) - len(kegglist) - kmcount[mod] + dip, pval ) )
    #print type(pval), pval
    #print type(kmodtestsize), kmodtestsize

    dip = 0
    kcenrich = {}
    kcodds = {}

    for fn in kgroupclist:
        for ko in kegglist:
            if fn in kgroupc[ko]:
                dip += 1
        #sys.stdout.write("         In Path  Not in Path\nDEG    :%-7d %d\nnonDEG:%-7d %d" % (dip, len(kegglist) - dip, kccount[fn]-dip, len(keggcount) - len(kegglist) - kccount[fn] + dip ))
        oddrat, pval = fisher_exact([
            [dip, len(kegglist) - dip],
            [kccount[fn]-dip, size - len(kegglist) - kccount[fn] + dip]
        ])
        #sys.stdout.write(pval)
        #sys.stdout.flush()
        kcodds[fn]   = oddrat
        kcenrich[fn] = pval / kgrpctestsize
        dip = 0     # reset for next module

    #print kmodenrich
    return kmodenrich, kcenrich
    ## Fisher's Exact Test:
    #               In Pathway:         Not in Pathway:                                         SUM:
    #   DEG     :   dip                 len(kegglist) - dip                                     len(kegglist)
    #   non-DEG :   kmcount[mod]-dip    len(keggcount) - len(kegglist) - kmcount[mod] + dip     len(keggcount) - len(kegglist)
    #   SUM     :   kmcount[mod]        len(keggcount) - kmcount[mod]                           len(keggcount)
    #

    pass

def kegg_pathway_enrichment(degs, negs, dbpaths=dbpaths, show_all=True, pthresh=0.01):
    """
    For  a given list of differentially expressed genes (and the background of non
    differentially expressed genes), determines which KEGG pathways are significantly
    enriched. Returns both module enrichment P-values, and upper-level classification
    P-values (as dictionaries).
    """

    deg_num_ko, deg_keggs = cbir_to_kegg(degs)
    neg_num_ko, neg_keggs = cbir_to_kegg(negs)

    print "%-4d kegg pathways from %d DEGs" % (len(deg_keggs), len(degs) )
    print "%-4d kegg pathways from %d nonDEGs" % (len(neg_keggs), len(negs) )

    # create dictionary of kegg pathways {pathwaytype:{pathway:[ko1,ko2,ko3]}}
    pathwaytype_dict = {}
    pathway_dict = {}
    pathway_lookup = {}

    print "extracting pathways..."
    ko1_h = open(dbpaths['kegg'], 'rb')
    for line in ko1_h:
        if line[0] == 'B':   # Kegg path type eg: B  <b>Replication and repair</b>
            pathtype_f = re.search('B.*<b>(.*)<', line)
            if pathtype_f is not None:
                pathtype = pathtype_f.group(1)
            else:
                pathtype = 'unknown'
            pathwaytype_dict[pathtype] = {}
        elif line[0] == 'C':     # Kegg Pathway eg: 01200 Carbon metabolism [PATH:ko01200]
            pathway_f = re.search("C +([0-9]*) *(.*)\[PATH", line)
            if pathway_f is not None:
                pathway_id = pathway_f.group(1)
                pathway_name = pathway_f.group(2)
            else:
                pathway_id = 'unknown'
                pathway_name = 'unknown'
            pathway_dict[pathway_id] = {}
            pathway_lookup[pathway_id] = pathway_name
        elif line[0] == 'D':   # Kegg term eg: K00844  HK; hexokinase [EC:2.7.1.1]
            koterm_f = re.search("(K[0-9]*)", line)
            if koterm_f is not None:
                koterm = koterm_f.group(1)
            else:
                koterm = 'unknown'
            pathwaytype_dict[pathtype][koterm] = 1
            pathway_dict[pathway_id][koterm] =  1


    print "calculating enrichment..."
    pathwaytype_ps = {}
    pathway_ps = {}
    # count number of degs and negs in each pathway:
    for pathwaytype in pathwaytype_dict:
        pwtsize = len(pathwaytype_dict)
        degs_in_path = sum([1 for ko in pathwaytype_dict[pathwaytype] if ko in deg_keggs])
        negs_in_path = sum([1 for ko in pathwaytype_dict[pathwaytype] if ko in neg_keggs])
        degs_not_in = len(deg_keggs) - degs_in_path
        negs_not_in = len(neg_keggs) - negs_in_path

        oddrat, pval = fisher_exact([ [degs_in_path, degs_not_in],
                                    [negs_in_path, negs_not_in] ],
                                      alternative='greater')
        pathwaytype_ps[pathwaytype] = pval

        if pval < pthresh:
            print "%s\n         \
            In Path  Not in Path\n\
            DEG    :  %-7d %d\n\
            non-DEG:  %-7d %d\n\
            Odds Ratio:%.3f\n\
            P-value:%.4f\n" % (pathwaytype,degs_in_path,degs_not_in,negs_in_path,negs_not_in,
            oddrat, pval)


    for pathway in pathway_dict:
        pwtsize = len(pathway_dict)
        degs_in_path = sum([1 for ko in pathway_dict[pathway] if ko in deg_keggs])
        negs_in_path = sum([1 for ko in pathway_dict[pathway] if ko in neg_keggs])
        degs_not_in = len(deg_keggs) - degs_in_path
        negs_not_in = len(neg_keggs) - negs_in_path

        oddrat, pval = fisher_exact([ [degs_in_path, degs_not_in],
                                    [negs_in_path, negs_not_in] ],
                                      alternative='greater')
        pathway_ps[pathway + ' ' + pathway_lookup[pathway]] = pval

    ## Fisher's Exact Test:
    #               In Pathway:         Not in Pathway:
    #   DEG     :   degs_in_path        degs_not_in
    #   non-DEG :   negs_in_path        negs_not_in
    #

    return pathwaytype_ps, pathway_ps

def collect_kegg_pathways(minsize=0, filename=dbpaths['keggpathways']):
    """
    output: smallrefined_dict = { pathwayname:[list of Cbir_genes] }
    """
    kp_h = open(filename, 'rb')
    bigpathway_dict = {}
    smallpathway_dict = {}
    for line in kp_h:
        if line.split()[0] == 'C': # eg C  Glycan Metabolism
            bigpath_c = " ".join(line.split()[1:])
            bigpathway_dict[bigpath_c] = []
        elif line.split()[0] == 'D': # eg D   M00072 Oligosaccharyltransferase [PATH:map00510]
            smallpath_c = " ".join(line.split()[1:])
            smallpathway_dict[smallpath_c] = []
        elif line.split()[0] == 'E': # eg E   Cbir_09279 POLD1; DNA poliymerase delata subunit 1 [EC:2.7.7.7]
            bigpathway_dict[bigpath_c].append(line.split()[1])
            smallpathway_dict[smallpath_c].append(line.split()[1])

    bigrefined_dict = {}
    smallrefined_dict = {}
    for path in bigpathway_dict:
        if len(bigpathway_dict[path]) >= minsize:
            bigrefined_dict[("KEGG pathway", path)] = bigpathway_dict[path]
    for path in smallpathway_dict:
        if len(smallpathway_dict[path]) >= minsize:
            smallrefined_dict[("KEGG module", path)] = smallpathway_dict[path]

    return bigrefined_dict, smallrefined_dict

def collect_ipr_pathways(minsize=0, filename=dbpaths['iprlist']):
    """
    output: iprrefined_dict = { pathwayname:[list of Cbir_genes] }
    """
    ipr_h = open(filename, 'rb')
    ipr_dict = {}
    ipr_defs_dict = {}
    for line in ipr_h:
        geneid = line.split()[0]
        iprids = re.findall('IPR[0-9]+', line)
        iprdefs = re.findall('(IPR[0-9]*);\s*([^\t\n]*)', line)
        for ipr in iprids:
            if ipr in ipr_dict:
                ipr_dict[ipr].append(geneid)
            else:
                ipr_dict[ipr] = [geneid]
        for ipr in iprdefs:
            ipr_defs_dict[ipr[0]] = ipr
    # remove any IPRs with less than the minsize number of Cbir genes specified:
    iprrefined_dict = {}
    for ipr in ipr_dict:
        if len(ipr_dict[ipr]) >= minsize:
            iprrefined_dict[ipr_defs_dict[ipr]] = ipr_dict[ipr]

    return iprrefined_dict

def collect_go_pathways(minsize=0, filename=dbpaths['goterms']):
    """
    output: iprrefined_dict = { pathwayname:[list of Cbir_genes] }
    """
    go_h = open(filename, 'rb')
    go_dict = {}
    go_defs_dict = {}
    for line in go_h:
        geneid = line.split()[0]
        goids = re.findall('GO:[0-9]+', line)
        godefs = re.findall("(GO:[0-9]*);\s*([^;]*);\s([\w ]*)", line)
        for go in goids:
            if go in go_dict:
                go_dict[go].append(geneid)
            else:
                go_dict[go] = [geneid]

        for go in godefs:
            go_defs_dict[go[0]] = (go[0],go[1],go[2])

    # remove any GOs with less than the minsize number of Cbir genes specified:
    gorefined_dict = {}
    for go in go_dict:
        if len(go_dict[go]) >= minsize:
            gorefined_dict[go_defs_dict[go]] = go_dict[go]

    return gorefined_dict

def cbir_ncbi(geneobj, dbpaths=dbpaths):
    "extract NCBI gene name for given Cbir reference ids"
    ncbi_h = open(dbpaths['ncbipep'], 'rb')

    ncbi_d = {}
    for line in ncbi_h:
        defsearch = re.search("\[(Cbir[A-Za-z0-9_\(\)\.\/]*)]", line)
        if defsearch is not None:
            ncbi_d[defsearch.group(1)] = line.strip()[1:]

    gene_gi = {}
    if type(geneobj) is list or type(geneobj) is dict:
        for gene in geneobj:
            try:
                gene_gi[gene] = ncbi_d[gene]
            except:
                gene_gi[gene] = 'no NCBI item'
    else:   # assume it is a single gene:
        try:
            gene_gi[geneobj] = ncbi_d[gene]
        except:
            gene_gi[geneobj] = 'no NCBI item'
    return gene_gi

######## Statistical Functions #########################################
def p_to_q(pvalues, display_on=False, cut1s=False, set_pi_hat=False):
    """
    Given the list of pvalues, convert to pFDR q-values.
    According to Storey and Tibshirani (2003) PNAS 100(16) : 9440

    returns: { p-value:q-value }

    """
    # because fisher's exact test gives highly skewed pvalue dists (with P of 1)
    # it may be necessary to remove the 1s before analysing
    if cut1s:
        pvalues = [ps for ps in pvalues if ps < 1]

    # order p-values:
    pvalues.sort()

    # estimate pi0:
    # evaluate pi0 across the range of lambda:
    lamrange = np.arange(0,0.95,0.01)
    #pbeaters = [ sum( p > lam for p in pvalues) for lam in lamrange ]
    #denominator = [ (len(pvalues) * (1 - lam)) for lam in lamrange ]
    pi0_lam = [ (sum( p > lam for p in pvalues) / (len(pvalues) * (1 - lam))) for lam in lamrange ]
    #pi0_hardway = []

    #for i in range(len(pbeaters)):
    #    pi0_hardway += [ pbeaters[i] / denominator[i] ]
    #if pi0_lam != pi0_hardway:
    #    print "\n\n\npi0_lam is not the same as pi0_hardway!\n\n"
    #print "pi0_hardway length:", len(pi0_hardway)
    #print "p_values size:", len(pvalues)
    # fit cubic spline to data, then calculate value of pi0 for lambda = 1:
    tck = interpolate.splrep(lamrange, pi0_lam, s=3)
    splinecurve = interpolate.splev(np.arange(0,1.0,0.01), tck, der=0)
    pi0_hat = interpolate.splev(1, tck, der=0)
    tck_half = 0
    if pi0_hat > 1:
        tck_half = interpolate.splrep(lamrange[:85], pi0_lam[:85], s=3)
        spline_half = interpolate.splev(np.arange(0,1.0,0.01), tck_half, der=0)
        pi0_hat_half = interpolate.splev(1, tck_half, der=0)
        pi0_hat = pi0_hat_half
        verbalise("R", "pi0_hat > 1! Likely skewed P-value distribution. Converting to ", pi0_hat_half)
    if set_pi_hat:
        pi0_hat = set_pi_hat
    if display_on:
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        try:
            n, bins, patches = ax1.hist(pvalues, bins=20, facecolor='green', label="P-values")
        except IndexError:
            ax1.plot(pvalues)
        plt.title('distribution of P-values')
        ax1.set_xlabel('lambda / P-value')
        ax1.set_ylabel('distribution #')
        plt.legend(loc=4)
        ax2 = ax1.twinx()
        ax2.plot(lamrange, pi0_lam, 'ro', np.arange(0,1.0,0.01), splinecurve, 'r', label='pi0_hat, s=3' )
        if tck_half != 0:
            ax2.plot(lamrange[:95], spline_half[:95], 'b', label='lambda < 0.85')
        ax2.set_ylabel('pi0_hat(lambda)')
        #ax1.plot(t, s1, 'b-')
        plt.legend(loc=1)
        plt.show()


    q_pm =  pi0_hat * pvalues[-1]    #  q(pm)
    # creates an ordered list of q(p(i)) values.
    q_pi_list = [q_pm] + [ (pi0_hat * len(pvalues)*pvalues[i])/i for i in range(len(pvalues)-1,1,-1)]
    # "The estimated q value for the ith most significant feature is q(p(i))"
    q_val = {}
    for i in range(len(pvalues)):
        q_val[pvalues[-1 * (i+1)]] = min(q_pi_list[:i+1])

    return q_val


########################################################################
def mainprog():
    cmmd, fname  = sys.argv

    cuffgenes = findgene(fname)
    print "Number of results:", len(cuffgenes)

    cbgenelist = []
    for result in cuffgenes:
        #print "Gene: %s Expr: %.3f\t%.3f\t%.3f\t(%s:%d-%d)" % (cuffgenes[result][0], cuffgenes[result][4], cuffgenes[result][5], cuffgenes[result][6], cuffgenes[result][1], cuffgenes[result][2], cuffgenes[result][3])
        #print cuffgenes[result][0]
        cbgenelist.append(cuffgenes[result][0])
    flylogs = fly_orthologs(cbgenelist)
    for gene in flylogs:
        print flylogs[gene]

    Cbi2 = extractseq('scaffold176',type='fa', OGS=dbpaths['ass'][:-3], startpos=92000, endpos=96000)
    Cbi3 = extractseq('scaffold197',type='fa', OGS=dbpaths['ass'][:-3], startpos=1552000, endpos=1557000)
    Cbi6 = extractseq('scaffold314',type='fa', OGS=dbpaths['ass'][:-3], startpos=347000, endpos=352000)
    crispr_f = '/Volumes/Genome/CRISPR/mini_genome.fa'
    crispr_h = open(crispr_f, 'w')
    crispr_h.write( ">scaf176 Cbi2\n%s\n>scaf197 Cbi3\n%s\n>scaf314 Cbi6\n%s\n" % (Cbi2, Cbi3, Cbi6) )
    crispr_h.close()

if __name__ == '__main__':
    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    if args.enrichment:
        verbalise("B", "Parsing GO terms")
        id_dic = build_go_terms(args.gofile)
        verbalise("G", [str(p) for p in id_dic.items()][:5])
        verbalise("B",
             "performing GO enrichment on file %s using GO file %s" % (os.path.basename(args.input), os.path.basename(args.datafile)))

        # perform Fisher's Exact Test for enrichment of GO terms:
        geneset = config.make_a_list(args.input, col_num=args.column)
        f_squares = go_enrichment(geneset, go_info=id_dic, return_odds=False, gofile=args.datafile)
        p_list = [p.p for p in f_squares]
        q_dic = p_to_q(p_list, display_on=True, cut1s=False, set_pi_hat=False)

        # prepare result strings:
        header = "Significantly enriched GO terms (Q-value <= 0.05)\n"
        resultstring = "\n".join([str(s) for s in f_squares if s.TwG >= args.min_terms and q_dic[s.p] <= args.max_q])
        resultlist = "\n".join(sorted([str(s.id_info) for s in f_squares if s.TwG >= args.min_terms and q_dic[s.p] <= args.max_q], key=lambda x: x[13]))
        #resultstring = "\n".join(sorted([str(id_dic[p[0]]) for p in pvalues.items() if q_dic[p[1][0]] <= 0.05 ], key=lambda x: x[1]))

        # output to screen:
        verbalise("M", header)
        verbalise(resultstring)
        verbalise("G", resultlist)

        # output to file:
        outfile = open(logfile[:-3] + "out", 'w')
        outfile.write(header)
        outfile.write(resultlist)
        outfile.write("\n\nFisher's Exact Test results for enriched GO terms:\n")
        outfile.write(resultstring)
        outfile.close()