#!/usr/bin/env python
""" A module for developing appropriate qPCR primers for a given gene"""

import os
import sys
import argparse
import re
import datetime
import tempfile
import itertools

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, IUPAC

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
    parser.add_argument("-m", "--max_distance", type=int, default=500,
                        help="""when checking primer mismatches, only report primer
                        pairs that are less than max_distance from one another
                        [ default = 500 ]""")
    parser.add_argument("-p", "--pcr", type=str,
                        help="""file containing primer pairs for analysis (not yet
                        implemented""")



    return parser

def check_pcr(primer1, primer2, gffobj, pcr_name="C_biroi", outputfile="~/tempfile", id_thresh=16, gap_thresh=1):
    """
    uses bwa to align both primers to the genome. The only filtering at this step is
    done at the level of the bwa parameters.
    """

    # set file paths:
    temp_dir = tempfile.mkdtemp()
    temp_primer1 = os.path.join(temp_dir, "primer1.fq")
    temp_primer2 = os.path.join(temp_dir, "primer2.fq")
    temp_sai1 = os.path.join(temp_dir, "primer1.sai")
    temp_sai2 = os.path.join(temp_dir, "primer2.sai")
    bwa_idx = '/Volumes/antqueen/genomics/indices/bwa/Cbir.assembly.v3.0.fa'

    # create temp fq file for primer1 and primer2:
    handle1 = open(temp_primer1, 'w')
    handle1.write("@primer1\n%s\n+\n%s\n" % (primer1, '~' * len(primer1)))
    handle1.close()

    handle2 = open(temp_primer2, 'w')
    handle2.write("@primer1\n%s\n+\n%s\n" % (primer2, '~' * len(primer2)))
    handle2.close()

    # set bwa alignment parameters:
    bwa_aln   = "bwa aln -n 5 -k 4 -N -t 30 %s %s > %s"
    bwa_samse = "bwa samse -n 10000000000  %s %s %s  > %s"

    status_return1 = os.system(bwa_aln % (bwa_idx, temp_primer1,temp_sai1))
    status_return2 = os.system(bwa_samse % (bwa_idx,
                                            temp_sai1,
                                            temp_primer1,
                                            temp_sai1[:-3]+'sam'))
    status_return3 = os.system(bwa_aln % (bwa_idx, temp_primer2,temp_sai2))
    status_return4 = os.system(bwa_samse % (bwa_idx,
                                            temp_sai2,
                                            temp_primer2,
                                            temp_sai2[:-3]+'sam'))

    if status_return1 > 0 or status_return2 > 0 or status_return3 > 0 or status_return4 > 0:
        verbalise("R", "Error performing bwa alignment for %s and %s" % (temp_primer1,
                                                                            temp_primer2))

    return (temp_sai1[:-3]+'sam', temp_sai2[:-3]+'sam')

def check_pcr_blast(primer1, primer2, gffobj, pcr_name="C_biroi", outputfile="~/tempfile", id_thresh=16, gap_thresh=1):
    """creates pseudo PCR product from primers, and blasts C.biroi to look for multiple
    PCR products. This function is now deprecated in favor of the much faster check_pcr,
    which uses bwa."""
    primer1_seq = Seq(primer1, generic_dna)
    primer2_seq = Seq(primer2, generic_dna).reverse_complement()
    PCR_product = primer1_seq + Seq("NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", generic_dna) + primer2_seq
    pcr_input_file = outputfile + ".pcr_product_input.tmp"

    # blast C. biroi genome with PCR_product:
    outfile = outputfile + ".blastout.tmp"
    cris.blastseq(PCR_product, seqnames=pcr_name, inpath=pcr_input_file,
                    outpath=outfile, outfmt='crisprtab')
    blast_dict = cris.parse_crispr_blast(
                        outfile,
                        gffobj,
                        id_thresh=id_thresh,
                        gap_thresh=gap_thresh,
                        shownearest=True
                        )
    cris.report_results(outfile + 'final_results.info', blast_dict)

def get_matches(sam):
    """
    extracts all aligned matches from SAM file
    returns dictionary:

    { scaffold:[
            (position, orientation, cigar, mismatches), ...
                ]
    }

    """
    matches_dic = {}
    handle = open(sam, 'rb')

    for line in handle:
        if line[0] == '@':
            continue
        elif len(line.split()) < 12:
            continue
        else:
            # calculate values for best match:
            elements = line.split()
            scaf = elements[2]
            cigar = elements[5]
            if elements[1] == '16':
                orientation = '-'
            else:
                orientation = '+'
            position = int(elements[3])
            mismatches = int(elements[4]) # not entirely true, as it represents map score.

            # add to dictionary:
            matches_dic[scaf] = [(position,orientation,cigar,mismatches),]

            # evaluate all secondary alignments:
            for tag in elements[11:]:
                if tag[:2] == 'XA':
                    step1 = [ aln.split(',') for aln in tag[5:].split(';')]
                    matches = [ (t[0],t[1][0],int(t[1][1:]),t[2],int(t[3])) for t in step1[:-1] ]
                    for scaf,orientation,position,cigar,mismatches in matches:
                        if scaf in matches_dic:
                            matches_dic[scaf].append((position,orientation,cigar,mismatches))
                        else:
                            matches_dic[scaf] = [(position,orientation,cigar,mismatches),]
    handle.close()
    return matches_dic

def pretty_pairs(p1, p2):
    "create a visualisation of the pcr product to evaluate its quality"
    if p1[0] < p2[0]:
        pleft = p1
        pright = p2
    else:
        pleft = p2
        pright = p1

    # set direction of template:
    if pleft[1] == '-':
        cleft = "%s:%d>>" % (pleft[2],pleft[3])
    else:
        cleft = "<<%s:%d" % (pleft[2],pleft[3])

    if pright[1] == '+':
        cright = "%s:%d>>" % (pright[2],pright[3])
    else:
        cright = "<<%s:%d" % (pright[2],pright[3])

    return "%s---%d---%s" % (cleft,pright[0]-pleft[0],cright)

def find_pairs(match1, match2, dist=500):
    """
    looks at two position dictionaries to see if there are any pairs within dist of each
    other. Also considers pairs arising from within each match dictionary separately.
    """
    # create a merged dictionary with all positions labeled by scaffold:
    combined = match1
    for k in match2.keys():
        if k in combined:
            combined[k] += match2[k]

    pairs = []
    for s in combined:
        if len(combined[s]) < 2:
            continue    # if only one match on a scaffold, then it cannot form a PCR product!
        else:
            for p1,p2 in itertools.combinations(combined[s], 2):
                if 0 < abs(p1[0] - p2[0]) <= dist:
                    pairs.append((s, min(p1[0], p2[0]), p1[3] + p2[3], pretty_pairs(p1, p2)))

    return pairs

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
        else:
            geneproduct = "None found"
            primer1 = "Primer1"
            primer2 = "Primer2"

        if deltaGmatch:
            if float(deltaGmatch.group(1)) > -8:
                out_h.write("%-20s %-22s %-22s\n%s" % (geneproduct, primer1, primer2, line) )

    fold_h.close()
    out_h.close()

def display_pcr_products(parsed_data):
    """
    for colorful display of potential pcr products!
    """
    for p in sorted(parsed_data, key=lambda x: x[2]):
        if p[2] == 0:
            verbalise("G", p)
        elif 0 < p[2] <= 2:
            verbalise("C", p)
        elif 2 < p[2] <= 4:
            verbalise("Y", p)
        elif 4 < p[2] <= 6:
            verbalise("M", p)
        elif 6 < p[2]:
            verbalise("R", p)


def main(args, logfile):
    genelist = config.make_a_list(args.gene_list)

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
        isoseqs = gfflib.extractseq(geneid, cds=True)
        boundaries = gfflib.extractseq(geneid, boundaries=True) # this is a dictionary
        exon_boundaries = {}
        for isoform in boundaries:
            cumlen = 0
            exon_boundaries[isoform] = []
            for xlen in boundaries[isoform]:
                cumlen += xlen
                exon_boundaries[isoform].append(cumlen)

        for isoform in isoseqs:

            verbalise("B", "preparing %s for primer design (%s)" % (isoform, geneid))
            count = 1

            for boundary in exon_boundaries[isoform][:-1]:
                # format for primer3 input file:
                input_txt = "SEQUENCE_ID=%s_exon%d\n" % (isoform, count) + \
                "SEQUENCE_TEMPLATE=%s\n" % (isoseqs[isoform]) + \
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




    print "\nAnalysing ~%d primer pairs:" % (cumcount * 5)

    # parse the primer3 output file
    p3_results = parse_p3(p3_out)

    # set up a dictionary to collate all the analyses performed on the primer pairs
    """
    all_results[(genename, primer1, primer2)] = {  'primers':(),
                                                   'pcr_products':(),
                                                   'rna_fold':(),

                                                }
    """
    all_results = {}


    count=0
    t0 = datetime.datetime.now()

    for name, primer1, primer2, start, finish in p3_results:
        count += 1
        all_results[(name, primer1, primer2)] = {'primers':(name,primer1,primer2,start,finish)}

        # align primers to genome:
        blastout = logfile[:-3] + ".".join([name.split(' ')[0][4:10], primer1, 'blastn'])
        sam1, sam2 = check_pcr(primer1, primer2, pcr_name=name, outputfile=blastout, gffobj=gfflib)

        # check for non-specific primer pairs:
        match1 = get_matches(sam1)
        match2 = get_matches(sam2)
        pairs = find_pairs(match1, match2, dist=args.max_distance)

        # clean up temporary files
        temp_dir = os.path.dirname(sam1)
        for f in os.listdir(temp_dir):
            os.remove(os.path.join(temp_dir,f))
        os.rmdir(temp_dir)

        # save results  - display using display_pcr_products()
        all_results[(name, primer1, primer2)]['pcr_products'] = sorted(pairs, key=lambda x: x[2])

        # extract PCR product and save to fasta file for tertiary structure analysis:
        isoform = re.search("(.*)_exon", name).group(1)
        pcr_seq_dic = gfflib.extractseq(isoform.split(" ")[2][5:], cds=True, trim_from=start , trim_to=finish )
        try:
            pcr_seq = pcr_seq_dic[isoform]
        except KeyError:
            verbalise("R", isoform, isoform.split(" ")[2][5:], "\n\n", pcr_seq_dic)
            exit()
        fastafile = logfile[:-3] + 'fa'
        defline = " ".join([name, primer1, primer2])
        cris.make_fasta(pcr_seq, seqnames=defline, outpath=fastafile, append='a')

        # calculate time left:
        t1 = datetime.datetime.now()
        tdiff = (t1-t0)
        trate = tdiff / count
        tleft = trate * cumcount * 5 - tdiff
        sys.stdout.write("\nAnalysed %d/%d (%.2f%%) pairs. Approx %s remaining\n\n" % (count,
                                                        cumcount * 5,
                                                        100.0 * count / (cumcount * 5),
                                                        tleft ))
        sys.stdout.flush()

    ## run VIENNA RNAfold to check tertiary structure
    print "Running RNAfold..."
    viennaout = logfile[:-3] + "RNAfold.info"
    cmd = "RNAfold -p -d2 --noLP -T 60 < " + fastafile + " > " + viennaout
    if os.system(cmd) > 0:
        verbalise("R", "Error running RNAfold")

    ## cull VIENNA output to acceptable structure list, output primers:
    parse_fold(viennaout)


################################################################

if __name__ == "__main__":

    dbpaths = config.import_paths()

    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    if args.pcr:
        pass
    else:
        main(args, logfile)