#!/usr/bin/env python
""" A module for developing appropriate qPCR primers for a given gene"""

import os
import argparse

import genematch
import gffparser

################################################################

def get_genelist(input_file):
    genelist = [line.split()[0] for line in open(input_file)]
    return genelist

def main(args):
    genelist = get_genelist(args.input_file)
    print "\nCreating genome parsing objects..."
    gffobj = gffparser.assemble_dict()
    quickinfo = gffparser.My_gff()

    # create output file (input file of primer3), in case of previous existence
    # NB: will overwrite existing file!!!
    p3_input = args.output_file + ".p3_input.txt"
    p3_h = open(p3_input, 'w')
    p3_h.close()


    for geneid in genelist:
        ## get sequence:
        seq = genematch.extractseq(geneid,type='cds')
        ## find exon boundaries:
        scaf = quickinfo.whichscaf(geneid)
        exon_details = [ subfeat.location  for feature in gffobj[scaf].features for subfeat in feature.sub_features if feature.qualifiers['ID'][0] == geneid ]
        exon_order = []
        # sort the exon lengths according to position and strand orientation:
        for exon in exon_details:
            exon_order.append((len(exon), exon.start))
            if exon.strand == -1:
                qrev = True
        exon_order.sort(key=lambda posn: posn[1], reverse=qrev)
        exon_boundaries = []
        # create exon boundary relative positions by cumulative addition:
        cumlen = 0
        for length,posn in exon_order:
            cumlen += length
            exon_boundaries.append(cumlen)
        print "###", geneid, exon_order, "\n", exon_boundaries
        count = 1
        for boundary in exon_boundaries:
            # format for primer3 input file:
            input_txt = "SEQUENCE_ID=%s_exon%d\n" % (geneid, count) + \
            "SEQUENCE_TEMPLATE=%s\n" % (seq) + \
            "SEQUENCE_TARGET=%d,2\n" % (boundary) +\
            "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/Users/POxley/apps/primer3-2.3.6/primer3_config/\n" + \
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
    ## run primer3:
    cmd = "primer3_core -format_output -output=%s %s" % (args.output_file, p3_input)
    os.system(cmd)



################################################################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Creates qPCR appropriate primers from a list of genes")
    parser.add_argument("-o", "--output_file", type=str, default="p3_output.info", help="File to save results to")
    parser.add_argument("-i", "--input_file", type=str,  help="File to analyse")
    parser.add_argument("-g", "--gff_file", type=str, default="/Volumes/Genome/armyant.OGS.V1.8.6.gff", help="GFF file for analyses")
    parser.add_argument("-c", "--cds_file", type=str, default="/Volumes/Genome/armyant.OGS.V1.8.6.cds", help="Genome fasta file for analyses")
    args = parser.parse_args()




    main(args)