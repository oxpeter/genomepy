#!/usr/bin/env python

"""
Takes a set of alternately spliced exons (for example), and randomly selects another exon
that is at the same position and within 1.2 fold expression value of the target exon.
"""

import os
import sys
import re
import argparse
import itertools
import datetime

import pandas as pd
import numpy as np
import random
import matplotlib.pyplot as plt

from genomepy import config
from genomepy import gffparser as gp

###### INITIALISE THE FILE PATHS NEEDED FOR ANALYSIS #######################

dbpaths = config.import_paths()

############################################################################

def define_arguments():
    parser = argparse.ArgumentParser(description=
            """
            Takes a set of alternately spliced exons (for example), and randomly
            selects another exon that is at the same position and within 1.2 fold
            expression value of the target exon. """)

    ### input options ###
    # logging options:
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="print fewer messages and output details")
    parser.add_argument("-o", "--output", type=str, default='genematch.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")

    # data file options:
    parser.add_argument("input", type=str, nargs=1,
                        help="file containing exon splicing feature to compare")
    parser.add_argument("--expression", type=str, nargs='+',
                        help="gene expression file for gene comparisons")
    parser.add_argument("-c", "--column", type=int, default=0,
                        help="column in which gene names are found (default = 0)")
    parser.add_argument("-g", "--gff", type=str, default=dbpaths['gff'],
                        help="gff file for extracting gene names (default available)")

    # analysis options:
    parser.add_argument("-D", "--display_off", action='store_true', default=False,
                        help="do not display graphs")
    parser.add_argument("-I", "--iterations", type=int, default=1,
                        help="number of random lists to generate (all appended to same file)")
    parser.add_argument("-t", "--alt_type", type=str, default='se',
                        help="""type of alternate splicing event you wish to analyse.
                        Options include se (skipped exons), mxe (mutually exclusive exons),
                        ri (retained introns), a5s (alternate 5' splicing) a3s (alternate
                        3' splicing)""")

    return parser

def parse_miso(line, gff, type):
    """
    parses a gff line from the miso output, allowing parsing of the ID attribute to
    identify which exons have undergone an alternate splicing event.
    """
    # define the relevant regex pattern that successfully identifies a given AS event:
    pattern = {'se':'\.se', 'mxe':'\.mxe\d',
                'a5s':'\.coreAndExt', 'a3s':'\.coreAndExt',
                'ri':'\.withRI', 'any':'\.'}

    fields = gp.parse_cols(line)
    if fields:
        if fields['type'] == 'exon':
            atts = gp.parse_atts(line)
            # see if gff line is an AS exon
            miso_s = re.search(pattern[type], atts['ID'])
            if miso_s:
                # find genes sharing location
                gene_s = gff.ingene((fields['scaf'], fields['start']))
                if len(gene_s) == 0:
                    gene_s = gff.ingene((fields['scaf'], fields['end']))
                    if len(gene_s) == 0:
                        return None

                return gene_s, fields['scaf'], fields['start'], fields['end']
            else:
                return None


def parse_miso_new(line, gff, type):
    """
    parses a gff line from the miso output, allowing parsing of the ID attribute to
    identify which exons have undergone an alternate splicing event.
    """
    # define the relevant regex pattern that successfully identifies a given AS event:
    pattern = {'se':'\.se', 'mxe':'\.mxe\d',
                'a5s':'\.coreAndExt', 'a3s':'\.coreAndExt',
                'ri':'\.withRI'}

    fields = gp.parse_cols(line)
    if fields:
        if fields['type'] == 'exon':
            atts = gp.parse_atts(line)

            # see if gff line is an AS exon
            miso_s = re.search(pattern[type], atts['ID'])
            if miso_s:
                # find genes sharing location
                gene_s = gff.findfeatures((fields['scaf'], fields['start']),
                                            strand=[fields['strand']],
                                            ftype='gene')
                if len(gene_s) == 0:
                    gene_s = gff.findfeatures((fields['scaf'], fields['start']),
                                                strand=[fields['strand']],
                                                ftype='gene')
                    if len(gene_s) == 0:
                        return None

                return ([ g.atts['Name'] for g in gene_s['gene']],
                        fields['scaf'],
                        fields['start'],
                        fields['end'])
            else:
                return None

def parse_expression(file):
    df = pd.read_csv(file, header=0, sep=" ", index_col=0)
    df['mean'] = df.mean(axis=1)
    return df

if __name__ == '__main__':
    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)
    outfile = logfile[:-3] + "out"

    verbalise("B", "Building master gene exon positions...")
    t0 = datetime.datetime.now()
    cbir = gp.My_gff(args.gff, primary_key='gene')
    cbir.build_master_gene()
    t2 = datetime.datetime.now()
    verbalise("M", "Time to complete old gff: %s" % datetime.timedelta.total_seconds(t2-t0))

    verbalise("M", "Building new GFF with masters...") ###################################
    t0 = datetime.datetime.now() ###################################
    cbirlib = gp.GffLibrary(args.gff) ###################################
    cbirlib.build_master_gene()       ###################################
    t2 = datetime.datetime.now() ###################################
    verbalise("M", "Time to complete new gff: %s" % datetime.timedelta.total_seconds(t2-t0))

    verbalise("B", "Extracting alternate exon positions...")
    # extract all alternately spliced exons:
    alt_exons = {}
    handle = open(args.input[0], 'rb')
    altered = []
    altered_new = [] ###################################
    for line in handle:
        data = parse_miso(line, cbir, args.alt_type)
        datanew = parse_miso_new(line, cbirlib, args.alt_type) ###################################
        if data:
            altered.append(data)
        if datanew: ###################################
            altered_new.append(datanew) ###################################
    handle.close()

    verbalise("M", "length of old altered: %d" % len(altered)) ##########################
    verbalise("M", "length of new altered: %d" % len(altered_new)) #############
    verbalise("M", altered[:5])
    verbalise("M", altered_new[:5])
    verbalise("M", sum( 1 for i in altered if i not in altered_new), "missed in new")#############

    no_master = []
    no_match = []
    gene_tracker = {}
    redundant_tracker = {}
    redundant_exons = []
    alt_transcripts = {}
    for gene, scaf, start, end in altered:
        for g in gene:
            # keep track of genes that appear repeatedly:
            if g in gene_tracker:
                gene_tracker[g] += 1
            else:
                gene_tracker[g] = 0

            if g in cbir.master:
                found = False
                for mex in cbir.master[g]:
                    if mex[0] <= start <= mex[1] or mex[0] <= end <= mex[1]:
                        found = True
                        if (g, mex) in redundant_tracker:
                            redundant_exons.append((g,mex[0], mex[1]))
                            break
                        else:
                            # Alternately-spliced exon found!!
                            exon_pos = cbir.master[g].index((mex[0], mex[1]))
                            alt_exons[g,gene_tracker[g]] = exon_pos

                            # determine which transcript it came from:
                            alt_transcripts[g, gene_tracker[g]] = []
                            for mrna in cbir.toplevel[g]['mRNA']:
                                #if (start, end) in cbir.exondict[mrna]:
                                alt_transcripts[g, gene_tracker[g]].append((mrna, exon_pos))

                            # mark gene and exon as used:
                            redundant_tracker[g, mex] = True
                            break
                else:
                    if not found:
                        no_match.append((g,scaf,start,end))
            else:
                no_master.append(g)

    exons_found = sum( 1 for t in altered for g in t[0] )
    if exons_found != len(alt_exons) + len(no_master) + len(no_match) + len(redundant_exons):
        verbalise("R", "Some exons are not accounted for!!")
        verbalise("R", "%d total exons found from parse_miso step" % exons_found)
        verbalise("R", "%d redundant exons found" % (len(redundant_exons)))
    elif exons_found == 0:
        verbalise("R", "No exons found! Exiting.")
        exit()
    if len(no_master) > 0:
        verbalise("R", "%d exons did not have a gene with a master" % (len(no_master)))
        verbalise("R", " ".join([str(g) for g in no_master[:10]]), "...")
    if len(no_match) > 0:
        verbalise("R", "%d exons did not match an exon in the master gene" % (len(no_match)))
        verbalise("R", "\n".join([str(g) for g in no_match[:6]]), "...")
    verbalise("G", "%d %s exons found and sent for partnering." % (len(alt_exons), args.alt_type))
    verbalise("B", "Picking a match for each alternately spliced exon...")

    count=0
    t0 = datetime.datetime.now()
    for iter in range(args.iterations):
        # report progress:
        count += 1
        t1 = datetime.datetime.now()
        tdiff = datetime.timedelta.total_seconds(t1-t0)
        speed = 1.0 * count / tdiff  # (counts/s)
        remaining = args.iterations - count
        tremaining = datetime.timedelta(seconds=remaining / speed)
        sys.stdout.write( "\r%d/%d (%.2f%%) complete. Time remaining ~ %s            " % (count,args.iterations, 100.0*count/args.iterations, tremaining))
        sys.stdout.flush()


        verbalise("Y",
            "\r%d/%d (%.2f%%) complete. Time remaining ~ %s            " % (count,
                    args.iterations,
                    100.0*count/args.iterations,
                    tremaining))


        # initialise result variables:
        matched_exons = {}
        expression = parse_expression(args.expression[0])
        option_size = []
        target_expression = []

        for gene,event in alt_exons:
            # get expression of target gene:
            target = expression.loc[gene]['mean']
            target_expression.append(target)

            # find all genes with expression within 1.2 fold of target:
            filtered = expression.loc[(expression["mean"] <= 1.2 * target) & (expression["mean"] >= 0.83 * target)].index

            # remove genes with no master gene positions (primarily non-coding RNA):
            have_master = [ g for g in filtered if g in cbir.master ]

            # remove genes that also have alternate exons:
            gene_options = [ g for g in have_master if (g,0) not in alt_exons and len(cbir.master[g]) >= alt_exons[gene,event] ]
            option_size.append(len(gene_options))

            # select a matched gene randomly, then append the corresponding exon to match_dic:
            retry = 1
            while retry:
                partner = random.choice(gene_options)
                try:
                    matched_exons[gene,event] = (partner,
                                            cbir.toplevel[partner]['scaffold'],
                                            cbir.master[partner][alt_exons[gene,event]])
                except IndexError:
                    retry += 1
                else:
                    retry = 0
                # set break in case we get stuck:
                if retry > 50 * len(gene_options):
                    print "No exon match found for gene ", gene
                    break

        # report results:
        verbalise("G", "Average number of matched genes to choose from: %.2f" % np.mean(option_size))
        verbalise("G", "Average expression of genes with alternately spliced exons (%s): %.2f" % ( args.alt_type, np.mean(target_expression)))

        unmatched_exons = [ g for g,i in alt_exons if (g,i) not in matched_exons ]
        if len(unmatched_exons) > 0:
            verbalise("R", "Could not find a match for %d exons:")
            verbalise("R", unmatched_exons)

        # output list of pairs:
        handle = open(outfile, 'a')
        handle.write("###\n")
        for gene_id, exon in alt_exons.items():
            t_scaf = cbir.toplevel[gene_id[0]]['scaffold']
            t_positions = cbir.master[gene_id[0]][exon]
            partner = matched_exons[gene_id][0]
            p_scaf = matched_exons[gene_id][1]
            p_positions = matched_exons[gene_id][2]
            handle.write("%-15s %-3d %-12s %-10d %-10d %-15s %-12s %-10d %d\n" % (gene_id[0],
                exon, t_scaf,
                t_positions[0], t_positions[1],
                partner, p_scaf, p_positions[0], p_positions[1] ))
        handle.close()

        # output list of transcripts containing alternate splicing:
        #alt_transcripts[g, gene_tracker[g]] = [(mrna, exon_position)]
        handle = open(outfile[:-3] + 'transcripts.out', 'w')
        for gene, idx in alt_transcripts:
            handle.write("%-16s %s\n" % (gene,
                                    " ".join([ str(x) for x in alt_transcripts[gene,idx]])))
        handle.close()

        if not args.display_off:
            plt.hist(target_expression, bins=50, alpha=0.5)
            plt.title("Expression (VST) of %d genes with alternately-spliced exons (%s)\n(mean = %.2f)" % (len(alt_exons), args.alt_type, np.mean(target_expression) ))
            plt.show()

            plt.hist(option_size, bins=50, alpha=0.5)
            plt.title("Number of comparable genes from which to choose for %s alternately-spliced exons\n(mean = %.2f)" % (args.alt_type, np.mean(option_size)))
            plt.show()



