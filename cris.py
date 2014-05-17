#!/usr/bin/python

""" A program to find CRISPR sites and analyse their suitability for use.

"""

import sys
import re
import os
import types
import argparse

import Bio.Blast.NCBIXML as xml
from Bio.Seq import Seq


import genematch

def findCRISPRsites(sequence_file='/Volumes/Genome/Cbir.assembly.v3.0_singleline.fa', pattern='([Gg][Gg].[Gg][Gg].{16}.[Gg][Gg])|([Cc][Cc].{19}[Cc][Cc].[Cc][Cc])'):
    " extract all CRISPR site matches from the genome as tuples containing each CRISPR sequence: "
    allmatches = []
    scaffmatches = []


    gfile = open(sequence_file, 'rb')
    #pattern = '(.{20}.[Gg][Gg]).{1,99}([Cc][Cc].{20})'
    #pattern = '([Gg][Gg].{18}.[Gg][Gg]).{1,30}([Cc][Cc].{18}[Cc][Cc])'
    #first_round pattern = '([Cc][Cc].{19}[Cc][Cc]).{1,20}([Gg][Gg].{19}[Gg][Gg])'
    #pattern = '([Gg][Gg].[Gg][Gg].{16}.[Gg][Gg])|([Cc][Cc].{19}[Cc][Cc].[Cc][Cc])'
    defpat = '([sCcafold]*[0-9]+)'
    #pattern = '(CGTAACAGTGTGTCGTGGCTCCCAGGGAGAGAAAGAGAGAGAGA)(GA)' # test pattern from inotocin
    linematches = [ (re.search(defpat, line), re.findall(pattern, line) ) for line in gfile ]
    for defmatch, sitematchlist in linematches:
        if not isinstance(defmatch, types.NoneType):
            currentscaf = defmatch.groups()[0]
        if len(sitematchlist) != 0:
            scaffmatches.append((currentscaf, sitematchlist))
            allmatches.extend(sitematchlist)
    gfile.close()
    #print scaffmatches[0:1]
    return allmatches, scaffmatches

def make_fasta(seq, seqnames=None, outpath='/Users/POxley/blastinput.tmp', append='w'):
    "Takes a list of sequences and creates a fasta formatted file"

    # create an input file for blasting, depending on if seq is a list or one sequence.
    if isinstance(seq, list):
        if seqnames:
            namedic = dict(zip(seq,seqnames))
        wobj = open(outpath, append)
        for item in seq:
            if seqnames:
                defline = ">lcl|sequence" + namedic[item] + "\n"
            else:
                position = str(seq.index(item))
                defline = ">lcl|sequence" + position + "\n"
            seqline = str(item) + "\n"
            wobj.write(defline)
            wobj.write(seqline)
        wobj.close()

    else:
        wobj = open(outpath, append)
        if seqnames:
            defline = ">lcl|" + seqnames + '\n'
        else:
            defline = ">lcl|sequence1\n"
        seqline = str(seq) + "\n"
        wobj.write(defline)
        wobj.write(seqline)
        wobj.close()

def blastseq(seq, seqnames=None, inpath='/Users/POxley/blastinput.tmp' , outpath="/Users/POxley/blastoutput.tmp", outfmt='tab'):
    "blast sequence given in input file against C.biroi genome"

    make_fasta(seq, seqnames, inpath)

    if outfmt == 'tab':
        # for tabular output:
        blastcmd_tab = 'blastn -query ' + inpath + ' -db /Volumes/Genome/BLAST_databases/nucleotide/Cbir.v3.0 -outfmt 6 -word_size 7 -evalue 1000000 -out ' + outpath
        os.system(blastcmd_tab)

    elif outfmt == 'xml':
        # for xml output:
        blastcmd_xml = 'blastn -query ' + inpath + ' -db /Volumes/Genome/BLAST_databases/nucleotide/Cbir.v3.0 -outfmt 5 -word_size 7 -evalue 1000000 -out ' + outpath
        os.system(blastcmd_xml)

    else:
        # if you want blast results that are viewer friendly:
        blastcmd_vis = 'blastn -query ' + inpath + \
            ' -db /Volumes/Genome/BLAST_databases/nucleotide/Cbir.v3.0 -out ' + \
            ' -out ' + outpath + \
            ' -outfmt 3 -word_size 7 -evalue 1000000'
        os.system(blastcmd_vis)

def findnearest(scaffold, hitpos, gff_p="/Volumes/Genome/armyant.OGS.V1.8.6.gff"):
    gffdict = {}
    gffobj = open(gff_p, 'rb')
    for line in gffobj:
        col = line.split()
        if col[2] == "mRNA":
            if col[0] not in gffdict:
                gffdict[col[0]] = []
            gffdict[col[0]].append((int(col[3]), int(col[4]), col[6], " ".join(col[8:])))
    gffobj.close()

    if scaffold not in gffdict:
        downstream = (999999, '@')
        upstream = (999999, '@')
    else:
        downstream = (999999, '@')
        upstream = (999999, '@')
        for startpos, stoppos, direction, gene in gffdict[scaffold]:
            if startpos <= hitpos <= stoppos:
                upstream = (0, gene)
                downstream = (0, gene)
                break
            if hitpos < startpos and startpos - hitpos < upstream[0]:
                upstream =  (startpos - hitpos, direction)
            if stoppos < hitpos and hitpos - stoppos < downstream[0]:
                downstream = (hitpos - stoppos, direction)
    return downstream, upstream

def scaflength(scaffold):
    scaffseq = genematch.extractseq(scaffold, type='fa', OGS='/Volumes/Genome/Cbir.assembly.v3.0')
    return len(scaffseq)

def showblast(filename='/Users/POxley/blastoutput.tmp', threshold=80, gaps_allowed=3, shownearest=True, seq="_unknown_", out_path="/Users/POxley/Documents/Analysis/Data/People/Buck/Second_Round/results.info"):
    "parses blast results in xml format (such as those created by blastseq()"

    blastresults = open(filename, 'rb')
    parsed_results = xml.parse(blastresults)

    #print "BLAST results parsed. Extracting hits with identities >= %d%%..." % (threshold)
    hitnumber = {}
    hitseq = {}
    hitcount = 0
    results_h = open(out_path, "a")
    for hit in parsed_results:
        #print "*** should be recording now ***"
        sys.stdout.flush()
        results_h.write( "\nHits for query %s (%s) above threshold of %d%%:\n" % (hit.query, seq, threshold) )
        for alignment in hit.alignments:
            #print "in alignment"
            #sys.stdout.flush()
            for hsp in alignment.hsps:
                #print "in hsp"
                #sys.stdout.flush()
                if seq is not "_unknown_":
                    seq_len = len(seq)
                    pcid = 100.0 * hsp.identities / seq_len
                else:
                    seq_len = hit.query_length
                    pcid = 100.0 * hsp.identities / seq_len
                if pcid >= threshold and hsp.gaps <= gaps_allowed:
                    hitcount += 1
                    #results_h.write( "%s\n" % (alignment.title) )
                    a_scaffold = re.search('[scCafold]*[0-9]*', alignment.title).group()
                    results_h.write( "%s\t\tsubject start: %d\n" % (a_scaffold, hsp.sbjct_start) )
                    if hit.query in hitnumber:
                        hitnumber[hit.query] += 1
                        hitseq[hit.query] = hsp.query
                    else:
                        hitnumber[hit.query] = 1
                        hitseq[hit.query] = hsp.query
                    results_h.write( "id: %d (%.2f%%)\tgaps: %d\tmismatches: %d\n" % (hsp.identities, pcid, hsp.gaps, seq_len - hsp.gaps - hsp.identities) )
                    if shownearest:
                        upstream, downstream = findnearest(a_scaffold, hsp.sbjct_start)
                        results_h.write( "distance to end of scaffold is %d\n" % ( scaflength(a_scaffold) - hsp.sbjct_start) )
                        results_h.write( "Closest gene upstream: %d (%s)\nClosest gene downstream: %d (%s)\n" % (upstream[0], upstream[1], downstream[0], downstream[1]) )
                    results_h.write( "query:  %s\n" % (hsp.query) )
                    results_h.write( "        %s\n" % (hsp.match) )
                    results_h.write( "genome: %s\n\n" % (hsp.sbjct) )
    results_h.write( "### %d matches for sequence %s\n\n" % (hitcount, seq) )
    results_h.close()
    blastresults.close()
    return hitnumber, hitseq

def pcr_creator():
	"""Take two primers and create a pseudo-PCR product by stitching the two together
	(with the second reverse complemented)"""
	return None

if __name__ == '__main__':

    ##### CLI Argument Parser #####
    parser = argparse.ArgumentParser(description="Finds and assesses CRISPR sites")
    parser.add_argument("-s", "--sequence_file", type=str, dest="seq_path", default='/Volumes/Genome/Cbir.assembly.v3.0_singleline.fa', help="Specify the .fa file containing the sequence to be searched")
    parser.add_argument("-p", "--pattern", type=str, dest="pattern", default='([Gg][Gg].[Gg][Gg].{16}.[Gg][Gg])|([Cc][Cc].{19}[Cc][Cc].[Cc][Cc])', help="Specify the regular expression string to search with.")
    parser.add_argument("-o", "--output_file", type=str, dest="out_path", default="/Users/POxley/Documents/Analysis/Data/People/Buck/crispr.info", help="Specify a file in which to save results.")
    parser.add_argument("-g", "--gaps_allowed", type=int, dest="gaps_allowed", default=0, help="Specify the number of gaps allowed in reported sequence matches.\nDefault = 0")
    parser.add_argument("-t", "--threshold", type=int, dest="threshold", default=80, help="Specify the minimum %%id to report in sequence matches.\nDefault = 80")
    #parser.add_argument("fastq_dir", type=str, help="The directory containing all the fastq files for analysis.")

    args = parser.parse_args()
    ###############################

    listcount = 0
    uniquedict = {}
    dupdict = {}
    CRISPRmatches, scafreport = findCRISPRsites(args.seq_path, args.pattern)
    print "\nNumber of scaffolds with matches found is %d" % (len(scafreport))
    for scaf, hitlist in scafreport:
        for seqresults in hitlist:
            seqresult = seqresults[0] + seqresults[1]
            listcount += 1
            if seqresult in uniquedict:
                del uniquedict[seqresult]
                dupdict[seqresult] = True
            elif seqresult in dupdict:
                pass
            else:
                uniquedict[seqresult] = ""

    resultslist = uniquedict

    """ the original code was looking for two matches and comparing both to all others.
    This led to having to duplicate the dictionary, with each sequence as key and its
    complement as value, and vice versa. After checking for the uniqueness of both seqs,
    the dictionary uniquedict was condensed into resultslist to eliminate the duplication.

    resultslist = {}

    for scaf, hitlist in scafreport:
        for seqresult in hitlist:
            if seqresult[0] in uniquedict:
                resultslist[seqresult[0]] = uniquedict[seqresult[0]]"""


    #for sequence in resultslist:
        #print sequence, resultslist[sequence]

    print listcount, "CRISPR sites found in sequence file"
    print len(resultslist) , "(uniqe CRISPR sites found)\n"

    # blast each pair of CRISPR sequences against the C. biroi genome:
    print "Blasting queries to C. biroi genome."
    print "Number of hits above threshold for each query will be saved to %s" % (args.out_path)
    count = 1
    for seq1 in resultslist:
        # seq1 = "GGN" + seq1  #   To search for extended sequence match
        print "%d) %s " % (count, seq1)
        count += 1
        blastseq(seq1)
        #print "sequence blasted"
        blasthits, blastseqs = showblast(seq=seq1, threshold=args.threshold, gaps_allowed=args.gaps_allowed, out_path=args.out_path)
        #print "parsed results saved to file."
        results_h = open(args.out_path, "a")
        for seqid in blasthits:
            if blasthits[seqid] < 2:
                results_h.write( "%s: %s\n" % (seqid, blastseqs[seqid]) )
        results_h.close()
