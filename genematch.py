#!/usr/bin/python

""" A series of functions allowing interaction between local gff/cds/pep/fa files and
Tuxedo suite result files and ncbi's BLAST engines
"""

import sys
import re

import Bio.Blast.NCBIWWW as ncbi
import Bio.Blast.NCBIXML as xml
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, IUPAC
from scipy.stats import fisher_exact



############################################

def extract_locus(fname, col_num, datacol=False):
    """ given a file name and a column number, extract_locus() will create three lists:
    scaffold,
    locus start pos
    locus end pos

    if datacol is provided, the contents of that column will be added to the extra list:
    gdata
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

def snp2gene(scaffold, pos, gff="/Volumes/Genome/armyant.OGS.V1.8.4.gff"):
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

def locus2gene(scaflist, gbeglist, gendlist, gdatalist=False, gff="/Volumes/Genome/armyant.OGS.V1.8.4.gff", comprehensive=True ):
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

def findgene(fname):
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
        fobj = open("/Volumes/Genome/armyant.OGS.V1.8.3.gff")
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

def extractseq(geneID, type='pep', OGS='/Volumes/Genome/armyant.OGS.V1.8.6', startpos=0, endpos=-1):
    """ extracts sequence of geneID from the current annotations. type is cds,  pep or fasta.
    """
    geneseq = ""
    fname = ".".join([OGS, type])
    fobj = open(fname)
    for line in fobj:
        query = re.search('>(lcl.)?([A-Za-z_\.0-9]*)', line)
        if query != None:
            if query.groups()[1] == geneID:
                thisline = fobj.next()

                while re.search('>', thisline) == None:
                    geneseq += thisline.strip()
                    try:
                        thisline = fobj.next()
                    except StopIteration:
                        break
                else:
                    break
    fobj.close()
    return geneseq[startpos:endpos]

def fly_orthologs(genelist):
    fobj = open('/Volumes/Genome/Genome_analysis/Orthologs/BGI.orthologs.list')
    orthodict = {}
    for line in fobj:
        if re.search('Cbir_[0-9]*', line) != None and re.search('DROME', line) != None:
            cbgene = re.search('Cbir_[0-9]*', line).group()
            if  cbgene in genelist:
                orthodict[cbgene] = re.search('DME\|(.*)_DROME', line).groups()[0]
    fobj.close()
    return orthodict

def gff_init():
    """ To speed up gene matching, this will create a gff dictionary """
    pass


######## Gene Ontology and KEGG pathway analyses ######################
def go_finder(genelist):
    pass

def kegg_enrichment(genelist):
    """
    For  a given genelist, determines which KEGG modules ( ~ pathways ) are significantly
    enriched. Returns both module enrichment P-values, and upper-level classification
    P-values (as dictionaries).
    """
    #print "Converting Genes to KEGGs"
    # convert genelist to KEGG list:
    kegg_h = open('/Volumes/Genome/Genome_analysis/KEGG_pathways/KEGG_orthologs.list', 'rb')

    keggd = {}
    keggcount = {}

    for line in kegg_h:
        if len(line.split()) == 2: # ie, the line has a cbir entry and KO
            keggd[line.split()[0]] = line.split()[1]
            keggcount[line.split()[1]] = True
        else:
            #keggd[line.split()[0]] = 'na'
            pass

    kegglist = []
    for gene in genelist:
        try:
            kegglist.append(keggd[gene])
        except KeyError:
            pass

    #print "Creating Kegg module library"
    ## create Kegg module dictionary:
    kmod_h = open('/Volumes/Genome/Genome_analysis/KEGG_pathways/ko00002.keg', 'rb')

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
        if dip > 0:
            print "%s\n         In Path  Not in Path\nDEG    :  %-7d %d\nnon-DEG:  %-7d %d\n%.4f\n" % (mod, dip, len(kegglist) - dip,kmcount[mod]-dip, len(keggcount) - len(kegglist) - kmcount[mod] + dip, pval )
        kmododds[mod]   = oddrat
        kmodenrich[mod] = pval * kmodtestsize
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
            [kccount[fn]-dip, len(keggcount) - len(kegglist) - kccount[fn] + dip]
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

########################################################################
def mainprog():
    print "welcome back to python"
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

    Cbi2 = extractseq('scaffold176',type='fa', OGS='/Volumes/Genome/Cbir.assembly.v3.0', startpos=92000, endpos=96000)
    Cbi3 = extractseq('scaffold197',type='fa', OGS='/Volumes/Genome/Cbir.assembly.v3.0', startpos=1552000, endpos=1557000)
    Cbi6 = extractseq('scaffold314',type='fa', OGS='/Volumes/Genome/Cbir.assembly.v3.0', startpos=347000, endpos=352000)
    crispr_f = '/Volumes/Genome/CRISPR/mini_genome.fa'
    crispr_h = open(crispr_f, 'w')
    crispr_h.write( ">scaf176 Cbi2\n%s\n>scaf197 Cbi3\n%s\n>scaf314 Cbi6\n%s\n" % (Cbi2, Cbi3, Cbi6) )
    crispr_h.close()

        #geneID = cuffgenes[result][0]
        #geneseq = extractseq(geneID, 'pep')
        #blastresult = blastgene(geneseq)
        #print "#" * 45

    """
    # To extract Vg promoters from all ants, I used the following under __main__:

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
    """


if __name__ == '__main__':

    #cmmd, fname, scaf, startpos, endpos  = sys.argv

    orco  = ["scaffold520", 276743, 276870]
    itr   = ["scaffold50", 1377329, 1377469]
    site1 = ["scaffold197", 1420598, 1420746]
    genelist = [orco, itr, site1]

    for site in genelist:
        sequence = extractseq(site[0], type='fa', OGS='/Volumes/Genome/Cbir.assembly.v3.0', startpos=site[1]-1000, endpos=site[2]+1000)
        print ">%s\n%s\n" % (site[0], sequence)
