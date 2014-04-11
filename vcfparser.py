#!/usr/bin/env python
""" A script to parse vcf allele information
"""

import argparse
import math
import os
import re
import sys

from Bio.Alphabet import generic_dna, IUPAC
import matplotlib.pyplot as plt
import numpy
from scipy import stats

import GFFparser as gp
import qvalue

########################################################################################

def parse_lines(category,filename):
    vcf_h = open(filename, 'rb')
    result_h = open( filename + "_results.list", 'w')

    for line in vcf_h:
        if re.search("#", line):
            result_h.write(line)
        else:
            scaffold = line.split()[0]
            pos      = int(line.split()[1])
            all1     = line.split()[3]
            all2     = line.split()[4]
            if all2 == ".":             # in case only reference alleles found at site
                all2 = all1
            info     = [ cat.split('=') for cat in line.split()[7].split(';')  ]

            #print info
            info_dict = {}
            for item_pair in info:
                info_dict[item_pair[0]]=item_pair[1].split(',')

            if category == 'DP4':
                all1dp = int(info_dict[category][0]) + int(info_dict[category][1])
                all2dp = int(info_dict[category][2]) + int(info_dict[category][3])
                requested = "%s\t%s" % (all1dp,all2dp)
            else:
                requested = '\t'.join(info_dict[category])
            result_h.write( "%s\t%d\t%s\t%s\t%s\n" % (scaffold, pos, all1, all2, requested) )

    vcf_h.close()
    result_h.close()

def fishers_exact(Forg_all1, Stat_all1, Forg_all2, Stat_all2):
    oddsratio, pvalue = stats.fisher_exact([[Forg_all1, Stat_all1], [Forg_all2, Stat_all2]])
    return pvalue

def binomial_prob(all1, all2, cumulative=True):
    n = all1 + all2
    r = min([all1, all2])
    if cumulative:
        # cumulative probability:
        cumsum = 0
        for i in range(int(r)):
            nCi = math.factorial(n) / ( math.factorial(i) * math.factorial(n-i) )
            exp = 0.5 ** n
            try:
                cumsum += nCi * exp
            except:
                cumsum = None # ???? Not sure how to handle this - may still just skip these entirely
                break
    else:
        # discrete probability
        nCr = math.factorial(n) / ( math.factorial(r) * math.factorial(n-r) )
        exp = 0.5 ** n
        try:
            Prob = nCr * exp
        except:
            if n >= 1012:
                Prob = None
            else:
                print "Could not calculate. n = %d\tnCr = %d\texp = %e" % (n, nCr, exp)
                Prob = None

    return cumsum

def find_ratios(file1):
    "finds the allele depths from the pre-processed files"
    snps = {}
    file1_h = open(file1, 'rb')
    for line in file1_h:
        if not re.search("#", line):
            columns = [ col for col in line.split() ]
            snps[columns[0],columns[1]] = ((columns[2],int(columns[4])), (columns[3],int(columns[5])))
    file1_h.close()
    return snps

def fdr_q(plist):
    "calculates the false discovery rate threshold for a list of p-values"
    m = len(plist)
    plist.sort()

    print "Min P-value: %e\nMax P-value: %e (%d tests)" % ( plist[0], plist[-1], m )

    # calculate FDR:
    maxp = -1
    for j in range(m):
        if plist[j] > j * 0.05 / m:
            maxp = plist[j]
            break
    print "Q threshold = %e" % ( maxp )

    return maxp

def find_q(plist):
    "given (unsorted) list of p-values, returns dictionary of { p-value:q-value }"

    pv_array = numpy.array(sorted(plist))
    qv_array = qvalue.estimate(pv_array)
    print qv_array
    qlist = list(qv_array)
    qdict = dict(zip(plist, qlist))

    # Fdict calculates expected number of false positives 'F' if the given q value is set
    # as the threshold.  Fdict[q-value] = Expected False Positive Number
    Fdict = {}
    for q in qlist:
        Fdict[q] = qlist.index(q) * q
    return qdict, Fdict

def critical_q(F_dict, crit=1):
    below = -10
    above = 1000000
    aboveq= -999
    belowq= -999

    if len(F_dict) < 2:
        print "all q-values are the same"
        best_value = -999
    else:

        for q in F_dict:
            if F_dict[q] == crit:
                print "q-value at which expected number of false positives is 1 is", q
                best_value = q
                break
            elif below < F_dict[q] < crit:
                below = F_dict[q]
                belowq = q
            elif above > F_dict[q] > crit:
                above = F_dict[q]
                aboveq = q
        else:
            print "q-value at which expected number of false positives is 1 is between %.4f (%.2f) and %.4f (%.2f)" % (belowq, below, aboveq, above)
            best_value = (aboveq + belowq) / 2
    return best_value


def compare_allele_ratios_binomial(file1, file2, gff_obj, makeplot=False, min_depth=10 ):
    """compares the ratios of allele depths in two files. Alleles must be identical in
    both files."""
    # collect allele depths:
    snps_1 = find_ratios(file1)
    print "Allele depths calculated for", file1
    snps_2 = find_ratios(file2)
    print "Allele depths calculated for", file2
    print ""

    comp_loci = {}
    comp_depths = {}
    for locus in snps_1:
        if locus in snps_2:                                 # ie, if the locus is present in both files
            if snps_1[locus][1][0] == snps_2[locus][1][0]:  # checking that alt allele is the same.
                comp_loci[locus]   = (snps_1[locus][0][0], snps_2[locus][1][0])  # stores the ref and alt alleles
                comp_depths[locus] = (snps_1[locus][0][1], snps_2[locus][0][1], snps_1[locus][1][1], snps_2[locus][1][1])
                # comp_depths[locus] in order: Forg_all1, Stat_all1, Forg_all2,  Stat_all2
    pvalues = {}
    combined = {}
    pFvalues = {}
    pSvalues = {}

    for locus in comp_depths:
        F1 = float(comp_depths[locus][0])
        S1 = float(comp_depths[locus][1])
        F2 = float(comp_depths[locus][2])
        S2 = float(comp_depths[locus][3])

        # to account for rare reads of ref or alternate alleles in polymorphisms.
        proceed = True
        for value in [F1,F2,S1,S2]:
            if  value <= min_depth: # also can use 0 < value <= min_depth
                proceed = False
                break

        if proceed:
            pforg = binomial_prob(F1, F2)
            pstat = binomial_prob(S1, S2)
            if pforg is not None and pstat is not None: # large values give overflow errors and return None
                pvalues[locus] = ( pforg, pstat )
                pFvalues[locus] = pforg
                pSvalues[locus] = pstat
                pstar = float(pforg) * float(pstat)
                newp = pstar * ( 1 - float(numpy.log(pstar + 0.000000000001)) )  # from James Theiler
                combined[locus] = newp


    # perform FDR calculation on combined tests:
    #plist = sorted(combined.values())
    q_value = fdr_q(combined.values())
    combined_qvalues, combined_F = find_q(combined.values())
    print "Combined pvalues: len - %d  min - %.4f  max - %.4f " % ( len(combined), min(combined.values()), max(combined.values()) )

    forg_qvalues, forg_F = find_q(pFvalues.values())
    print "Forg pvalues: len - %d  min - %.4f  max - %.4f " % ( len(pFvalues), min(pFvalues.values()), max(pFvalues.values()) )
    stat_qvalues, stat_F = find_q(pSvalues.values())
    print "Forg pvalues: len - %d  min - %.4f  max - %.4f " % ( len(pSvalues), min(pSvalues.values()), max(pSvalues.values()) )

    # report q value at which the expecte false positive rate is 1:
    print "\nFORAGING:"
    q_critF = critical_q(forg_F)
    print "STATARY:"
    q_critS = critical_q(stat_F)
    print "COMBINED:"
    q_critB = critical_q(combined_F)


    ## add all information to a dictionary for reporting:
    all_loci = {}
    for locus in combined:
        genefeat = gp.snp_in_gene(locus[0], int(locus[1]), gff_obj)
        all_loci[locus] = {'EFP_stat':stat_F[stat_qvalues[pvalues[locus][1]]], 'EFP_forg':forg_F[forg_qvalues[pvalues[locus][0]]],'pC':combined[locus] , 'pF':pvalues[locus][0] , 'pS':pvalues[locus][1] , 'qC':combined_qvalues[combined[locus]] , 'qF':forg_qvalues[pvalues[locus][0]] , 'qS':stat_qvalues[pvalues[locus][1]], 'F1':comp_depths[locus][0], 'S1':comp_depths[locus][1], 'F2':comp_depths[locus][2], 'S2':comp_depths[locus][3], 'REF':comp_loci[locus][0], 'ALT':comp_loci[locus][1], 'gene_id':genefeat["gene_id"], 'gene_name':genefeat["gene_name"] }

    if makeplot:
        plot_pvalues(combined.values())
        fvalues = []
        svalues = []
        for key in pvalues:
            fvalues.append(pvalues[key][0])
            svalues.append(pvalues[key][1])
        plot_pvalues(pvalues.values())
        plot_pvalues(fvalues)
        plot_pvalues(svalues)

    return  all_loci

def find_editing(all_loci):
    "takes the pvalues and qvalues and log-expression values and finds something useful"
    pass

def compare_allele_ratios(file1, file2, gff_obj):
    """compares the ratios of allele depths in two files. Alleles must be identical in
    both files."""
    snps_1 = find_ratios(file1)
    print "allele depth dictionary calculated for", file1
    snps_2 = find_ratios(file2)
    print "allele depth dictionary calculated for", file2

    comp_loci = {}
    comp_depths = {}
    for locus in snps_1:
        if locus in snps_2:
            if snps_1[locus][1][0] == snps_2[locus][1][0]:
                comp_loci[locus]   = (snps_1[locus][0][0], snps_2[locus][1][0])
                comp_depths[locus] = (snps_1[locus][0][1], snps_2[locus][0][1], snps_1[locus][1][1], snps_2[locus][1][1])

    print "Performing Fisher's Exact Test for all loci found in both files. Please wait a moment."
    pvalues = {}
    for locus in comp_depths:
        # in order: Forg_all1, Stat_all1, Forg_all2, Stat_all2
        pvalues[locus] = fishers_exact(comp_depths[locus][0],comp_depths[locus][1],comp_depths[locus][2],comp_depths[locus][3])

    # perform FDR calculation on combined tests:
    q_value = fdr_q(pvalues.values())
    fishers_qvalues, fishers_F = find_q(pvalues.values())
    print "Summary of Fisher's Exact p-values: len - %d  min - %.4f  max - %.4f " % ( len(pvalues), min(pvalues.values()), max(pvalues.values()) )

    # report q value at which the expected false positive rate is 1:
    print "\nFisher's Exact Test:"
    q_critF = critical_q(fishers_F)

    ## Add all information to dictionary all_loci for reporting:
    all_loci = {}
    for locus in pvalues:
        genefeat = gp.snp_in_gene(locus[0], int(locus[1]), gff_obj)
        all_loci[locus] = { 'pv':pvalues[locus], 'qv':fishers_qvalues[pvalues[locus]], 'expF':fishers_F[fishers_qvalues[pvalues[locus]]] ,  'F1':comp_depths[locus][0], 'S1':comp_depths[locus][1], 'F2':comp_depths[locus][2], 'S2':comp_depths[locus][3], 'REF':comp_loci[locus][0], 'ALT':comp_loci[locus][1], 'gene_id':genefeat["gene_id"], 'gene_name':genefeat["gene_name"] }


    """sig_loci = {}
    for locus in pvalues:
        if pvalues[locus] < maxp:
            sig_loci[locus] = (pvalues[locus], comp_loci[locus], comp_depths[locus])

    plot_pvalues(plist)"""
    return all_loci


def binary_out(all_loci, filename):
    """ Takes the significant loci from binary analysis and produces useful output files.
    """

    print "\n%d loci printed to file %s." % (len(all_loci), filename)
    all_loci_h = open(filename, 'w')

    all_loci_h.write("%-12s%-9s%-7s%-7s%-8s%-10s%-7s%-7s%-7s%-5s%-5s%-10s%-7s%-7s%-6s%-7s\n" % ("scaffold","pos","gene_id","REF","ALT","F_1:F_2","log2","p(F)\tq(F)","EFP(F)","S_1:S_2","log2","p(S)","q(S)","EFP(S)","p(C)","q(C)"))  # create header in file.

    for locus in all_loci:
        logF = numpy.log2( (all_loci[locus]['F1'] / (all_loci[locus]['F2'] + 0.00001) ) )
        logS = numpy.log2( (all_loci[locus]['S1'] / (all_loci[locus]['S2'] + 0.00001) ) )

        try:
            all_loci_h.write( "%-18s%-9s" % (locus[0], locus[1]) \
            + "%(gene_id)-12s%(REF)-5s    %(ALT)-5s   %(F1)3d:%(F2)-3d\t" % all_loci[locus] \
            + "%-9.5f" % (logF) \
            + "%(pF).5f %(qF).5f %(EFP_forg)%d %(S1)3d:%(S2)-3d\t" % all_loci[locus] \
            + "%-9.5f" % (logS) \
            + "%(pS).5f %(qS).5f %(EFP_stat)%d  %(pC).5f  %(qC).5f\n" % all_loci[locus] )
        except:
            print all_loci[locus]
            print asfdsaf

    all_loci_h.close()
    cmd_line = "sort -k 3,3 -k2,2n  " + filename +  "  > Results.binary.sorted.list"
    os.system(cmd_line)

def fishers_output(all_loci, filename):
    """ Takes significant loci from fisher's exact test and outputs some useful info."""
    print "\n%d loci printed to file %s." % (len(all_loci), filename)
    all_loci_h = open(filename, 'w')

    all_loci_h.write("scaffold          pos      gene_id     REF      ALT   F_1:F_2\tlog2\t\tS_1:S_2\tlog2\tp(FET)\tq(FET)\n")  # create header in file.

    for locus in all_loci:
        logF = numpy.log2( (all_loci[locus]['F1'] / (all_loci[locus]['F2'] + 0.00001) ) )
        logS = numpy.log2( (all_loci[locus]['S1'] / (all_loci[locus]['S2'] + 0.00001) ) )

        try:
            all_loci_h.write( "%-18s%-9s" % (locus[0], locus[1]) \
            + "%(gene_id)-12s%(REF)-5s    %(ALT)-5s   %(F1)3d:%(F2)-3d   " % all_loci[locus] \
            + "%-9.5f" % (logF) \
            + "%(S1)3d:%(S2)-3d\t" % all_loci[locus] \
            + "%-9.5f" % (logS) \
            + "%(pv).5f  %(qv).5f %(expF)d\n" % all_loci[locus] )
        except:
            print all_loci[locus]
            print asfdsaf

    all_loci_h.close()
    cmd_line = "sort -k 12,12n -k3,3  " + filename +  "  > Results.fishers.sorted.list"
    os.system(cmd_line)

def allele_specific_gene_expression(all_loci, filename, gff_obj):
    all_loci_h = open(filename, 'w')
    all_loci_h.write("scaffold\tpos\tP(F)\tP(S)\tREF\tALT\tF_1 : F_2\tS_1 : S_2\trat_F\trat_S\n")  # create header in file.

    # find all alleles, their expression bias, and assign them to their gene:
    genelist = {}   # keys are the gene ids.
    for locus in all_loci:
        F1 = float(all_loci[locus][2][0])
        S1 = float(all_loci[locus][2][1])
        F2 = float(all_loci[locus][2][2])
        S2 = float(all_loci[locus][2][3])
        genefeat = gp.snp_in_gene(locus[0], int(locus[1]), gff_obj)
        logF = numpy.log2( (F1 / (F2 + 0.00001) ) )
        logS = numpy.log2( (S1 / (S2 + 0.00001) ) )
        if genefeat["gene_id"] in genelist:
            genelist[genefeat["gene_id"]]['S'].append(logF)
            genelist[genefeat["gene_id"]]['S'].append(logS)
            genelist[genefeat["gene_id"]]['FS'].append(logF * logS)
            genelist[genefeat["gene_id"]]['pS'].append(all_loci[locus][1][0])
            genelist[genefeat["gene_id"]]['pS'].append(all_loci[locus][1][1])
            genelist[genefeat["gene_id"]]['pFS'].append(all_loci[locus][0])

        else:
            genelist[genefeat["gene_id"]] = {'S':[logF, logS], 'FS':[logF * logS], 'pS':[all_loci[locus][1][1], all_loci[locus][1][0]], 'pFS':[all_loci[locus][0]]}

        all_loci_h.write( "%-18s%-10s%-11s\t%.5f\t%.5f\t%.5f\t%.5f\t%-3s%-3s%d:%d\t%d:%d\t%.3f\t%.3f\t%.3f\t%s\n" % (locus[0], locus[1], genefeat["gene_id"], all_loci[locus][1][0], all_loci[locus][1][1], all_loci[locus][0], all_loci[locus][3], all_loci[locus][4][0], all_loci[locus][4][1], F1, F2, S1, S2, logF, logS, logF * logS, genefeat["gene_name"]) ) # for binomial probability

    all_loci_h.close()
    cmd_line = "sort -k 1,1 -k2,2n " + filename + " > Results.all_loci.binary.sorted.list"
    os.system(cmd_line)

    # output those genes with all alleles showing biased expression:
    O_biased = {}
    FS_biased = {}
    for geneid in genelist:
        if min(genelist[geneid]['S']) > 0.5 and ( max(genelist[geneid]['S']) - min(genelist[geneid]['S']) < 1 ) and len(genelist[geneid]['pS']) > 2 and ( all( x > 0 for x in genelist[geneid]['S'] ) or all( x < 0 for x in genelist[geneid]['S'] ) ):
            O_biased[geneid] = genelist[geneid]['S']
        if len(genelist[geneid]['FS']) > 2 and ( all(x < 0 for x in genelist[geneid]['FS']) ):
            FS_biased[geneid] = genelist[geneid]['FS']

    for geneid in O_biased:
        print "Biased(one-way):  %-10s  %d\t%r" % (geneid, len(O_biased[geneid]), O_biased[geneid])
    for geneid in FS_biased:
        print "Biased(two-ways): %-10s  %d\t%r" % (geneid, len(FS_biased[geneid]), FS_biased[geneid])


def plot_pvalues(plist):
    x = range(len(plist))
    #plt.ylim(0,0.05)
    #plt.xlim(0,1000)
    plt.plot(x, sorted(plist))
    plt.show()


########################################################################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parses VCF files. To perform allele ratio test, first run on each VCF file seperately.")
    parser.add_argument("vcf_file", type=str, help="The VCF file for parsing")
    parser.add_argument("-C", "--cat", type=str, dest="category", default='DP4', help="Specify the category to parse")
    parser.add_argument("-r", "--ratios", type=str, dest="ratio_file", default=False, help="second vcf file for comparing allele ratios using Fisher's Exact test\n(Requires vcf_results.list files)")
    args = parser.parse_args()




    if args.ratio_file:

        ## create gff object for gene identification:
        print "\nCreating gff object. Please wait a moment."
        gff_obj = gp.assemble_dict(in_file="/Volumes/Genome/armyant.OGS.V1.8.6_lcl.gff", in_seq_file="/Volumes/Genome/Cbir.assembly.v3.0_gi.fa", features_only=False)
        print "gff object created.\n"


        ## FISHERS EXACT TEST:
        #all_loci = compare_allele_ratios(args.vcf_file, args.ratio_file, gff_obj)
        #print "%d loci found." % (len(all_loci))
        #newfile = args.vcf_file + "_sig_loci.fishers.list"
        #fishers_output(all_loci, newfile)
        #print "### Fisher's Exact Test analysis complete ###\n"

        ## BINARY ANALYSIS:
        all_loci = compare_allele_ratios_binomial(args.vcf_file, args.ratio_file, gff_obj, makeplot=True)
        print "%d loci found." % (len(all_loci))
        newfile = args.vcf_file + "_sig_loci.binary.list"
        binary_out(all_loci, newfile)
        print "### Binary Test analysis complete ###\n"

        ## look for genes with all loci with biased expression:
        #newfile = args.vcf_file + "_all_loci.binary.list"
        #allele_specific_gene_expression(all_loci, newfile, gff_obj)



    else:
        parse_lines(args.category, args.vcf_file)



