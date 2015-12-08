#!/usr/bin/env python
""" A module for developing appropriate qPCR primers for a given gene"""

import os
import argparse
import re
import datetime

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, IUPAC
import progressbar

from genomepy import genematch, gffparser, cris, config

################################################################
def define_arguments():
    parser = argparse.ArgumentParser(description="Creates qPCR appropriate primers from a list of genes")
    ## admin functions:
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="print fewer messages and output details")
    parser.add_argument("-o", "--output", type=str, default='qpcr.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")
    parser.add_argument("-D", "--display_on", action='store_true',default=False,
                        help="display graph results (eg for p value calculation)")

    ## input options:
    parser.add_argument("-i", "--gene_list", type=str,
                        help="""List of genes to design primers for. Can be a
                        comma-delimited list, or the path to a file containing gene names
                        in the first column.""")
    parser.add_argument("-g", "--gff_file", type=str, default=dbpaths['gff'],
                        help="GFF file for analyses")
    parser.add_argument("-a", "--assembly_file", type=str, default=dbpaths['ass'],
                        help="genome fasta file")
    parser.add_argument("-c", "--cds_file", type=str, default=dbpaths['cds'],
                        help="Fasta file of coding sequences")
    parser.add_argument("-P", "--PCR", type=str,
                        help="""file containing primer pairs for analysis (not yet
                        implemented""")

    return parser

def check_PCR(primer1, primer2, gffobj, pcr_name="C_biroi", outputfile="~/tempfile", id_thresh=16, gap_thresh=1):
    """creates pseudo PCR product from primers, and blasts C.biroi to look for multiple
    PCR products."""
    primer1_seq = Seq(primer1, generic_dna)
    primer2_seq = Seq(primer2, generic_dna).reverse_complement()
    PCR_product = primer1_seq + Seq("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", generic_dna) + primer2_seq

    # blast C. biroi genome with PCR_product:
    outfile = outputfile + ".blastout.tmp"
    cris.blastseq(PCR_product, seqnames=pcr_name, outpath=outfile, outfmt='crisprtab')
    blast_dict = cris.parse_crispr_blast(
                        outfile,
                        gffobj,
                        id_thresh=id_thresh,
                        gap_thresh=gap_thresh,
                        shownearest=True
                        )
    cris.report_results(outfile + 'final_results.info', blast_dict)

def parse_p3(input_file):
    p3_h = open(input_file, 'rb')
    for line in p3_h:
        matchname = re.search("PRIMER PICKING RESULTS FOR (.*)", line)
        matchleft = re.search("LEFT PRIMER +([0-9]+) +([0-9]+)[ 0-9\.]*([AGCTN]+)", line)
        matchright =re.search("RIGHT PRIMER +([0-9]+) +([0-9]+)[ 0-9\.]*([AGCTN]+)", line)
        if matchname:
            name = matchname.group(1)
        if matchleft:
            primer1 = matchleft.group(3)
            primer1_start = int(matchleft.group(1))
        if matchright:
            primer2 = matchright.group(3)
            primer2_end = int(matchright.group(1)) + int(matchright.group(2))
            yield name, primer1, primer2, primer1_start, primer2_end

def parse_fold(fold_output):
    "pulls out the sequences with acceptable level of secondary structure (delta)G > -9 kJ/mol"

    fold_h = open(fold_output, 'rb')
    out_h  = open(fold_output + ".acceptable.info", 'w')

    for line in fold_h:
        defmatch = re.search(">(lcl\|)?([^ ]+) ([ACTGN]+) ([ACTGN]+)",line)
        deltaGmatch  = re.search("\(( ?.[0-9]+.*)\)", line)
        if defmatch:
            geneproduct = defmatch.group(2)
            primer1     = defmatch.group(3)
            primer2     = defmatch.group(4)

        if deltaGmatch:
            if float(deltaGmatch.group(1)) > -8:
                out_h.write("%-20s %-22s %-22s\n%s" % (geneproduct, primer1, primer2, line) )

    fold_h.close()
    out_h.close()


def main(args, logfile):
    genelist = config.make_a_list(args.gene_list)
    print "\nCreating genome parsing objects..."
    gffobj = gffparser.assemble_dict()
    verbalise("B", "Assembling gff file...")
    quickinfo = gffparser.My_gff()

    gfflib = gffparser.GffLibrary(args.gff_file, args.assembly_file)

    # create output file (input file of primer3), in case of previous existence
    # NB: will overwrite existing file!!!
    p3_input = logfile[:-3] + "p3_input.txt"
    p3_h = open(p3_input, 'w')
    p3_h.close()

    cumcount = 0
    for geneid in genelist:
        ## get sequence:
        #seq = genematch.extractseq(geneid)
        seq = gfflib.extractseq(geneid, cds=True)
        boundaries = gfflib.cds_boundaries(geneid) # this is a dictionary
        for isoform in boundaries:
            cumlen = 0
            exon_boundaries = []
            for xlen in boundaries[isoform]:
                cumlen += xlen
                exon_boundaries.append(cumlen)
        # current setting therefore only designs primers for the last isoform!

        """
        ## find exon boundaries:
        scaf = quickinfo.whichscaf(geneid)
        exon_details = [ subfeat.location  for feature in gffobj[scaf].features for subfeat in feature.sub_features if feature.qualifiers['ID'][0] == geneid ]
        exon_order = []

        # sort the exon lengths according to position and strand orientation:
        for exon in exon_details:
            exon_order.append((len(exon), exon.start))
            if exon.strand == -1:
                qrev = True
            else:
                qrev = False
        exon_order.sort(key=lambda posn: posn[1], reverse=qrev)
        exon_boundaries = []

        # create exon boundary relative positions by cumulative addition:
        cumlen = 0
        for length,posn in exon_order:
            cumlen += length
            exon_boundaries.append(cumlen)
        """


        print "preparing %s for primer design" % (geneid)
        count = 1

        for boundary in exon_boundaries[:-1]:
            # format for primer3 input file:
            input_txt = "SEQUENCE_ID=%s_exon%d\n" % (geneid, count) + \
            "SEQUENCE_TEMPLATE=%s\n" % (seq) + \
            "SEQUENCE_TARGET=%d,2\n" % (boundary) +\
            "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/home/peter/lib/primer3_config/\n" + \
            "PRIMER_TASK=generic\n" + \
            "PRIMER_PICK_LEFT_PRIMER=1\n" + \
            "PRIMER_PICK_INTERNAL_OLIGO=0\n" + \
            "PRIMER_PICK_RIGHT_PRIMER=1\n" + \
            "PRIMER_OPT_SIZE=20\n" + \
            "PRIMER_MIN_SIZE=18\n" + \
            "PRIMER_MAX_SIZE=22\n" + \
            "PRIMER_MAX_NS_ACCEPTED=0\n" + \
            "PRIMER_PRODUCT_SIZE_RANGE=50-150\n" + \
            "P3_FILE_FLAG=0\n" + \
            "PRIMER_PAIR_MAX_DIFF_TM=1.0\n" + \
            "PRIMER_EXPLAIN_FLAG=1\n" + \
            "PRIMER_OPT_TM=60\n" + \
            "PRIMER_MAX_TM=60\n" + \
            "PRIMER_MIN_TM=58\n" + \
            "PRIMER_MAX_END_GC=2\n" + \
            "=\n"


            # append to primer3 input file:
            p3_h = open(p3_input, 'a')
            p3_h.write(input_txt)
            p3_h.close()

            count += 1
            cumcount += 1

    ## run primer3:
    #print "Running primer3 on %d exon boundaries" % (cumcount)
    p3_out = logfile[:-3] + "p3_out.info"
    cmd = "primer3_core -format_output -output=%s %s" % (p3_out, p3_input)
    os.system(cmd)

    ## blast primer3 primer pairs to check for specificity, and check for tertiary structure:
    # parse primer3 output to extract primer pairs:
    print "\nAnalysing ~%d primer pairs:" % (cumcount * 5)

    p3_gen = parse_p3(p3_out)


    #bar = progressbar.ProgressBar(maxval=cumcount*5, \
    #widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.ETA()])
    count=0
    t0 = datetime.datetime.now()
    #bar.update(count)
    for name, primer1, primer2, start, finish in p3_gen:
        count += 1
        #print "Processing primer pair for %s" % (name)
        # blast primer pair:
        blastout = logfile[:-3] + ".".join([name, primer1, 'blastn'])
        check_PCR(primer1, primer2, pcr_name=name, outputfile=blastout, gffobj=quickinfo)

        ## check blastn output for viable PCR products:
        #hitnumber, hitseq = cris.showblast(filename=blastout + ".blastout.tmp", \
        #    threshold=80, gaps_allowed=3, shownearest=False, \
        #    seq="_unknown_", out_path=blastout + ".parsed.info")

        # prepare to check for tertiary structure:
        # extract PCR product and save to fasta file:
        geneid = re.search("(.*)_exon", name).group(1)
        pcr_seq = genematch.extractseq(geneid, type='cds', startpos=start , endpos=finish )
        fastafile = logfile[:-3] + 'fa'
        defline = " ".join([name, primer1, primer2])
        cris.make_fasta(pcr_seq, seqnames=defline, outpath=fastafile, append='a')

        # calculate time left:
        t1 = datetime.datetime.now()
        tdiff = (t1-t0)
        trate = tdiff / count
        tleft = trate * cumcount * 5 - tdiff

        os.system.write("\rAnalysed %d/%d (%.2f%%) pairs. Approx %s remaining         \r" % (count,
                                                        cumcount * 5,
                                                        100.0 * count / (cumcount * 5),
                                                        tleft ))
        os.system.flush()
        #bar.update(count)
    #bar.finish()

    ## run VIENNA RNAfold to check tertiary structure
    print "Running RNAfold..."
    viennaout = args.output + ".RNAfold.info"
    cmd = "RNAfold -p -d2 --noLP -T 60 < " + fastafile + " > " + viennaout
    os.system(cmd)

    ## cull VIENNA output to acceptable structure list, output primers:
    parse_fold(viennaout)


################################################################

if __name__ == "__main__":

    dbpaths = config.import_paths()

    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    if args.PCR:
        pass
    else:
        main(args, logfile)