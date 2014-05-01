#!/usr/bin/env python

""" Imports a whitespace-delimited expression matrix and produces a hierarchically
clustered heatmap. Can also build the expression matrix from cufflinks data file.

Multiple filtering options and matrix normalisation methods are now implemented.

T-test and ANOVA testing are also available to test individual rows/columns.

PCA and read distribution plots can be displayed to check for normalcy, clustering etc.



"""

import re
import argparse

def get_tfs(infile, psize=600, limit=10):
    # store TF positions:
    weed_h = open(infile, 'rb')
    tf_tracker = {}
    tf = 'eg'
    tf_seq = {}
    homolog_count = 0
    for line in weed_h:
        if re.search('([0-9]+)\)', line) is not None:
            tf = int(re.search('([0-9]+)\)', line).groups()[0])
            homolog_count = 0
        if len(line.split()) == 2 :
            try:
                pos = int(line.split()[1])
            except ValueError:
                continue
            if tf <= limit:
                try:
                    tf_tracker[homolog_count].append( (tf, pos + psize) )
                except KeyError:
                    tf_tracker[homolog_count] = [ (tf, pos + psize) ]
                homolog_count += 1
        elif len(line.split()) == 5:
            try:
                pos = int(line.split()[3])
            except ValueError:
                continue
            if tf <= limit:
                tf_seq[tf] = line.split()[0]
                try:
                    tf_tracker[homolog_count].append( (tf, pos + psize) )
                except KeyError:
                    tf_tracker[homolog_count] = [ (tf, pos + psize) ]
                homolog_count += 1

    weed_h.close()

    for sample in tf_tracker:
        tf_tracker[sample].sort(key=lambda posn: posn[1])

    return tf_tracker, tf_seq

def write_tfs(outfile, tf_posns, tf_seqs):

    colorme = {1:'\033[1;31m|1|\033[0m',2:'\033[1;32m|2|\033[0m',3:'\033[1;33m|3|\033[0m',4:'\033[1;34m|4|\033[0m',5:'\033[1;35m|5|\033[0m',6:'\033[1;36m|6|\033[0m',7:'\033[1;37m|7|\033[0m',8:'\033[1;38m|8|\033[0m',9:'\033[1;39m|9|\033[0m',10:'\033[1;30m|10|\033[0m',11:'\033[1;34m|11|\033[0m'}

    out_h = open(outfile, 'w')
    for tf in tf_seqs:
        out_h.write("%s = %s\n" % (colorme[tf], tf_seqs[tf]))
    #out_h.write(str(tf_posns))
    out_h.write("     " + "-" * 200 + '\n')
    for sample in tf_posns:
        #out_h.write(str(tf_posns[sample]) + '\n')
        out_h.write("%2s:: " % (sample) )
        prev_val = 0
        sumdif = 0
        for tup in tf_posns[sample]:
            diff = (tup[1] - prev_val) / 3
            prev_val = tup[1]
            sumdif += diff
            out_h.write("-" * diff)
            out_h.write("%s" % (colorme[tup[0]]))
            #out_h.write("%-2s" % ("ABCDEFGHIJKLMNOPQRSTUVWXYZzbcdefghijklmnopqrstuvwxyz"[tup[0]]))
        out_h.write("-" * ((600 - prev_val)/3) + "###" + str(sumdif) + '\n')
    out_h.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Converts output from weederH into a viewable format")
    parser.add_argument("-w", "--weeder_file", type=str, help="The weeder output file for analysing")
    args = parser.parse_args()

    tf_posns, tf_seqs = get_tfs(args.weeder_file)
    outfile = args.weeder_file + ".info"
    write_tfs(outfile, tf_posns, tf_seqs)
