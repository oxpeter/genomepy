#!/usr/bin/python

""" A program to find CRISPR sites and analyse their suitability for use.

"""

import sys
import re
import os
import types
import argparse
import time
import datetime

import Bio.Blast.NCBIXML as xml
from Bio.Seq import Seq

from genomepy import config, gffparser, genematch

###### INITIALISE THE FILE PATHS NEEDED FOR ANALYSIS #######################

dbpaths = config.import_paths()
verbalise = config.check_verbose(True) # mainly to allow running in Django

############################################

def define_arguments():
    parser = argparse.ArgumentParser(description="Finds and assesses CRISPR sites")

    ### input options ###
    # logging options:
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="print fewer messages and output details")
    parser.add_argument("-o", "--output", type=str, default='crispr_design.out',
                        help="specify the filename to save results to")
    parser.add_argument("-d", "--directory", type=str,
                        help="specify the directory to save results to")
    parser.add_argument("-D", "--display_on", action='store_true',default=False,
                        help="display graph results (eg for p value calculation)")

    ##### CLI Argument Parser #####
    parser.add_argument("-s", "--sequence_file", type=str, dest="seq_path",
                        default=dbpaths['assone'],
                        help="Specify a fasta file containing the sequence to be searched")
    parser.add_argument("--gene", type=str, default="",
                        help="""Specify a gene or gene feature to be searched. The genomic
                        DNA from the feature will be used to search for gRNA sites. Eg:
                        LOC105281428_Q""")
    parser.add_argument("-p", "--pattern", type=str, dest="pattern",
                        default='(?=([Gg][Gg].{18}.[Gg][Gg]))|(?=([Cc][Cc]..{18}[Cc][Cc]))',
                        help="Specify the regular expression string to search with.")
    parser.add_argument("-r", "--guide_rna", type=str,
                        help="""Use a pre-defined guide RNA sequence:\nSX - short
                        exogenous\nLX - long exogenous\nSN - short endogenous\nLN -
                        long endogenous""")
    parser.add_argument("-g", "--gaps_allowed", type=int, dest="gaps_allowed", default=0,
                        help="""Specify the number of gaps allowed in reported sequence
                        matches.\nDefault = 0""")
    parser.add_argument("-i", "--id_thresh", type=int,  default=18,
                        help="Specify the minimum id to report in sequence matches.\nDefault = 18")
    parser.add_argument("-t", "--threshold", type=int, dest="threshold", default=76,
                        help="""Specify the minimum %%id to report in sequence matches.\n
                        Default = 76""")
    parser.add_argument("-x", "--exogenous", action='store_true',
                        help="""Use if sequence is exogenous, to do second blast with GG
                        added to seq.""")
    parser.add_argument('--threads', type=int, default=1,
                        help='Select the number of threads to perform calculation on (default = 1)')

    #parser.add_argument("fastq_dir", type=str, help="The directory containing all the fastq files for analysis.")

    return parser

def findCRISPRsites(sequence_file=dbpaths['ass'], pattern='(?=(.{20}.[Gg][Gg]))|(?=([Cc][Cc]..{20}))'):
    " extract all CRISPR site matches from the genome as tuples containing each CRISPR sequence: "
    defpat = '([sCcafold]*[0-9]+)'

    allmatches = []
    scaffmatches = []
    defline = 'init'
    seqline = ''
    seqdict = {}

    # create dictionary of sequences keyed from their deflines
    gfile = open(sequence_file, 'rb')
    for line in gfile:
        if line[0] == '>':
            if defline != 'init': # ie, skip the parsing on the first occurrence of '>'
                # parse defline:
                defparsed = re.search(defpat, defline)
                if defparsed is not None:
                    definition = defparsed.group(1)
                else:
                    definition = defline.split()[0]

                # add sequence to dictionary
                seqdict[definition] = seqline
                # search line for patterns:
                all = [m[0] + m[1] for m in re.findall(pattern, seqdict[definition])]

                if len(all) > 0:
                    scaffmatches.append((definition, all))
                    allmatches.extend(all)

            # reset for next sequence:
            defline = line.strip()
            seqline = ''
        else:
            seqline += line.strip()
    else:
        if defline != 'init': # ie, skip the parsing on the first occurrence of '>'
            # parse defline:
            defparsed = re.search(defpat, defline)
            if defparsed is not None:
                definition = defparsed.group(1)
            else:
                definition = defline.split()[0]

            # add sequence to dictionary
            seqdict[definition] = seqline
            # search line for patterns:
            all = [m[0] + m[1] for m in re.findall(pattern, seqdict[definition])]

            if len(all) > 0:
                scaffmatches.append((definition, all))
                allmatches.extend(all)

    gfile.close()

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
                defline = ">sequence" + namedic[item] + "\n"
            else:
                position = str(seq.index(item))
                defline = ">sequence" + position + "\n"
            seqline = str(item) + "\n"
            wobj.write(defline)
            wobj.write(seqline)
        wobj.close()

    else:
        wobj = open(outpath, append)
        if seqnames:
            defline = ">" + seqnames + '\n'
        else:
            defline = ">%s\n" % (seq[:25])
        seqline = str(seq) + "\n"
        wobj.write(defline)
        wobj.write(seqline)
        wobj.close()

def blastseq(seq, seqnames=None, inpath='/Users/POxley/blastinput.tmp' , outpath="/Users/POxley/blastoutput.tmp", outfmt='tab', threads=1):
    "blast sequence given in input file against C.biroi genome"
    t0 = datetime.datetime.now()
    make_fasta(seq, seqnames, inpath)

    if outfmt == 'tab':
        # for tabular output:
        blastcmd = 'blastn -query ' + inpath + ' -db ' + dbpaths['blastnuc'] + '/Cbir.v3.0 -outfmt 6 -word_size 4 -evalue 1000000 -out ' + outpath + ' -num_threads ' + str(threads)

    elif outfmt == 'xml':
        # for xml output:
        blastcmd = 'blastn -query ' + inpath + ' -db ' + dbpaths['blastnuc'] + '/Cbir.v3.0 -outfmt 5 -word_size 4 -evalue 1000000 -out ' + outpath + ' -num_threads ' + str(threads)

    elif outfmt == 'crisprtab':
        # for tab delimited values used in crispr search parsing:
        blastcmd = 'blastn -query ' + inpath + ' -db ' + dbpaths['blastnuc'] + "/Cbir.v3.0 -outfmt '6 qseqid qlen nident mismatch gaps sseq sseqid sstart send qstart' -word_size 4 -evalue 1000000 -out " + outpath + ' -num_threads ' + str(threads)

    else:
        # if you want blast results that are viewer friendly:
        blastcmd_vis = 'blastn -query ' + inpath + \
            ' -db ' + dbpaths['blastnuc'] + '/Cbir.v3.0 -out ' + \
            ' -out ' + outpath + \
            ' -outfmt 3 -word_size 4 -evalue 1000000' + \
            ' -num_threads ' + str(threads)
    os.system(blastcmd)
    t1 = datetime.datetime.now()
    #verbalise("B", "Time to blast sequence: %s" % datetime.timedelta.total_seconds(t1-t0))

def findnearest(scaffold, hitpos, gff_p=dbpaths['gff']):
    "deprecated. Now use gffparser gff_obj.findnearest(scaffold, posn)"
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
    scaffseq = genematch.extractseq(scaffold, type='fa', OGS=dbpaths['ass'])
    return len(scaffseq)

def parse_crispr_blast(blast_results, gffobj, id_thresh=18, gap_thresh=0, shownearest=True):
    "parses the tab delimited blastn results to filter out sequences that are no good"
    results_dict = {}
    blastresults_h = open(blast_results, 'rb')
    blastparser = (line.split() for line in blastresults_h)
    # for each sequence, create a list of blast results (each result contained as a dict):
    for qseqid, qlen, nident, mismatch, gaps, sseq, sseqid, sstart, send, qstart in blastparser:
        # check to see if result has PAM sequence or not:
        if sseq[-2:].upper() == 'GG':
            haspam = True
            cutsite = int(send) - 3
        elif sseq[:2].upper() == 'CC':
            haspam = True
            cutsite = int(sstart) + 3
        else:
            haspam = False
            cutsite = int(sstart)

        ## exonic = gffobj.inexon(cutsite, sseqid)

        if int(nident) >= id_thresh and int(gaps) <= gap_thresh:
            if shownearest:
                upstream, downstream = gffobj.findnearest((sseqid, cutsite), ftype='gene')
            else:
                upstream = (1,'not searched')
                downstream = (1,'not searched')

            # check to see if cutsite is in exon:
            exonic = False
            if upstream[0] == 0 or downstream[0] == 0:
                in_exon = gffobj.findfeatures((sseqid,cutsite), ftype='exon')
                if len(in_exon) > 0:
                    exonic = True


            curr_result = { 'downdist':downstream[0],
                            'downgene':" ".join([str(f) for f in downstream[1]]),
                            'updist':upstream[0],
                            'upgene':" ".join([str(f) for f in upstream[1]]),
                            'qlen':int(qlen),
                            'nident':int(nident),
                            'mismatch':len(qseqid)-int(nident),
                            'gaps':int(gaps),
                            'qseq':"%s%s" % (' ' * int(qstart), sseq),
                            'sseqid':sseqid,
                            'sstart':int(sstart),
                            'send':int(send),
                            'exonic':exonic,
                            'cutsite':cutsite,
                            'haspam':haspam}

            if qseqid in results_dict:
                results_dict[qseqid].append(curr_result)
            else:
                results_dict[qseqid] = [curr_result]

    blastresults_h.close()
    return results_dict

def showblast(blast_results='/Users/POxley/blastoutput.tmp', threshold=80, gaps_allowed=3, shownearest=True, seq="_unknown_", out_path="/Users/POxley/tmp.info"):
    "parses blast results in xml format (such as those created by blastseq()"

    blastresults_h = open(blast_results, 'rb')
    parsed_results = xml.parse(blastresults_h)

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
    blastresults_h.close()
    return hitnumber, hitseq

def pcr_creator():
	"""Take two primers and create a pseudo-PCR product by stitching the two together
	(with the second reverse complemented)"""
	return None

def create_log(args):
    """
    DEPRECATED: USE CONFIG.CREATE_LOGS INSTEAD
    """

    ## create output folder and log file of arguments:
    timestamp = time.strftime("%b%d_%H.%M")
    if not args.output:
        root_dir = os.getcwd()
        newfolder = root_dir + "/crispr." + timestamp + "." + args.guide_rna
        os.mkdir(newfolder)  # don't have to check if it exists, as timestamp is unique
        filename = newfolder + "/crispr." + timestamp + ".log"
    else:
        newfolder = os.path.realpath(args.filename)
        if os.path.exists(newfolder) is False:  # check to see if folder already exists...
            os.mkdir(newfolder)
        filename = newfolder + '/' + os.path.basename(args.filename) + '.' + timestamp + ".log"

    log_h = open(filename, 'w')
    log_h.write( "File created on %s\n" % (timestamp) )
    log_h.write( "Program called from %s\n\n" % (os.getcwd()) )
    for arg in str(args)[10:-1].split():
        log_h.write( "%s\n" % (arg) )
    log_h.close()
    return filename

def report_results(filename, blast_dict, header=True):
    results_h = open(filename, "a")
    if header:
        results_h.write( "   %-26s %-20s %-10s || %-4s %-4s %-4s %-5s || %-5s %s %s %s %s\n" % ('qseq','sseqid','cutsite','nident','mismatch','gaps','PAM','exonic?','downgene','downdist','upgene','updist') )
    for seqid in blast_dict:
        results_h.write( "\n## %s - %d bp (%d matches) %s\n" % (seqid, len(seqid), len(blast_dict[seqid]), '#' * 24) )
        for match in blast_dict[seqid]: # each match is a dictionary
            results_h.write( "  %(qseq)-26s %(sseqid)-20s %(cutsite)-10d || %(nident)-4d %(mismatch)-4d %(gaps)-4d PAM:%(haspam)-5s || %(exonic)-5s %(downgene)s %(downdist)-10d %(upgene)s %(updist)d\n" % (match) )
    results_h.close()

def main(args, logfile, verbalise):
    if args.guide_rna.lower() == 'sx':  # short exogenous
        args.pattern = '(?=(.{18}.[Gg][Gg]))|(?=([Cc][Cc]..{18}))'
        args.exogenous = True
    if args.guide_rna.lower() == 'lx':  # long exogenous
        args.pattern = '(?=(.{20}.[Gg][Gg]))|(?=([Cc][Cc]..{20}))'
        args.exogenous = True
    if args.guide_rna.lower() == 'sn':  # short endogenous
        args.pattern = '(?=([Gg][Gg].{16}.[Gg][Gg]))|(?=([Cc][Cc].{16}[Cc][Cc].[Cc][Cc]))'
    if args.guide_rna.lower() == 'ln':  # long endogenous
        args.pattern = '(?=([Gg][Gg].{18}.[Gg][Gg]))|(?=([Cc][Cc]..{18}[Cc][Cc]))'
    else:
        args.guide_rna = 'manual'

    ## gffobj = gffparser.My_gff() ##
    gffobj = gffparser.GffLibrary(dbpaths['gff'], dbpaths['ass'])
    verbalise("Y", gffobj)

    error_log = []
    if args.gene:
        sequences = gffobj.extractseq(args.gene, buffer=0)
        sequencefile = logfile[:-3] + 'input_file.fasta'
        handle = open(sequencefile, 'w')
        for defline in sequences:
            handle.write("%s\n%s\n" % (defline, sequences[defline]))
            if sequences[defline] == "":
                error_log.append( sequences[defline] )
            elif re.search("library$", str(sequences[defline])):
                error_log.append( sequences[defline] )
        handle.close()
        # change seq_path so that newly generated fasta files will be used:
        args.seq_path = sequencefile

        if len(error_log) > 0:
            verbalise("R", "\n".join(error_log))

        handle = open(args.seq_path, 'rb')
        for line in handle:
            verbalise("C", line.strip())
        handle.close()

    listcount = 0
    uniquedict = {}
    dupdict = {}
    CRISPRmatches, scafreport = findCRISPRsites(args.seq_path, args.pattern)
    verbalise("G", "\n%d scaffold%s with matches found." % (len(scafreport), 's' if len(scafreport) > 1 else '' ))
    for scaf, hitlist in scafreport:
        for seqresult in hitlist:
            listcount += 1
            if seqresult in uniquedict:
                del uniquedict[seqresult]
                dupdict[seqresult] = True
            elif seqresult in dupdict:
                pass
            else:
                uniquedict[seqresult] = ""

    resultslist = uniquedict




    verbalise("G", "%d CRISPR sites found in sequence file (%d unique)" % (listcount, len(resultslist)))



    blast_file = logfile[:-3] + 'blast.info'
    report_results(blast_file, {}, header=True)
    verbalise("Y", "Number of hits above threshold for each query will be saved to %s" % (blast_file))
    count = 1
    tstart = time.time()
    pbcount = 0
    for seq1 in resultslist:
        # progress tracker:
        if pbcount > 0:
            time_elapsed = time.time() - tstart
            time_per_cycle = 1. *  time_elapsed / pbcount
            time_remaining = (len(resultslist) - pbcount) * time_per_cycle
            m, s = divmod(time_remaining, 60)
            h, m = divmod(m, 60)
            if pbcount > 0:
                print '\r>> %d of %d genes complete. ETA to completion: %d:%02d:%02d                      ' % (pbcount, len(resultslist), h, m, s),
            sys.stdout.flush()
        else:
            time_remaining = len(resultslist) * 360
            m, s = divmod(time_remaining, 60)
            h, m = divmod(m, 60)
            print '>> 0 of %d genes complete. Completion in approx: %d:%02d:%02d                      ' % (len(resultslist), h, m, s),
            sys.stdout.flush()

        # perform blast, parse results, and write to file:
        pbcount += 1
        count += 1

        blastseq(seq1,
                inpath=logfile[:-3]+'blast_in.tmp',
                outpath=logfile[:-3]+'blast_out.tmp',
                outfmt='crisprtab',
                threads=args.threads)

        blast_dict = parse_crispr_blast(logfile[:-3]+'blast_out.tmp',
                                        gffobj,
                                        id_thresh=args.id_thresh,
                                        gap_thresh=args.gaps_allowed,
                                        shownearest=True)

        if count == 2:
            report_results(logfile[:-3]+'final_results.info', blast_dict)
        else:
            report_results(logfile[:-3]+'final_results.info', blast_dict, False)


    if args.exogenous: # have to blast sequence with the extra T7 GG appended to the search
        tstart = time.time()
        pbcount = 0
        for seq1 in resultslist:
            # progress tracker:
            if pbcount > 0:
                time_elapsed = time.time() - tstart
                time_per_cycle = 1. *  time_elapsed / pbcount
                time_remaining = (len(resultslist) - pbcount) * time_per_cycle
                m, s = divmod(time_remaining, 60)
                h, m = divmod(m, 60)
                if pbcount > 0:
                    print '\r>> %d of %d genes complete. ETA to completion: %d:%02d:%02d                      ' % (pbcount, len(resultslist), h, m, s),
                sys.stdout.flush()
            else:
                time_remaining = len(resultslist) * 90
                m, s = divmod(time_remaining, 60)
                h, m = divmod(m, 60)
                print '\n>> 0 of %d genes complete. Completion in approx: %d:%02d:%02d                      ' % (len(resultslist), h, m, s),
                sys.stdout.flush()

            # determine whether sequence is rev comp or straight, and append GG/CC:

            if seq1[-2:].upper() == 'GG' and seq1[:2].upper() == 'CC':
                print "sequence %s guides in both directions!" % seq1
                seq1 = 'GG' + seq1
            elif seq1[-2:].upper() == 'GG':
                seq1 = 'GG' + seq1
            elif seq1[:2].upper() == 'CC':
                seq1 += 'CC'

            pbcount += 1
            count += 1
            blastseq(seq1, inpath=logfile[:-3]+'blast_in.tmp',
                        outpath=logfile[:-3]+'blast_out.tmp', outfmt='crisprtab',
                        threads=args.threads)
            blast_dict = parse_crispr_blast(logfile[:-3]+'blast_out.tmp',
                                            gffobj,
                                            id_thresh=args.id_thresh,
                                            gap_thresh=args.gaps_allowed,
                                            shownearest=True)

            if count == 2:
                report_results(logfile[:-3]+'final_results.info', blast_dict)
            else:
                report_results(logfile[:-3]+'final_results.info', blast_dict, False)




    os.system('rm  ' + logfile[:-3]+'blast_in.tmp' )
    os.system('rm  ' + logfile[:-3]+'blast_out.tmp' )

if __name__ == '__main__':

    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    main(args, logfile, verbalise)

    ###############################

