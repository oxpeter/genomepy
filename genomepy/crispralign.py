#!/usr/bin/env python
" script to align miseq reads to detect indels caused by ZFN/TALEN/CRISPR "

# created by ben matthews 12/2013
# bmatthews@rockefeller.edu
#
# modified by Peter Oxley 12/24/2013
# poxley@rockefeller.edu

# pre-requisites: python packages
# pysam, genematch, argparse
#
# pre-requisites: other tools
# samtools, gsnap, pysamstats


# usage: crispralign.py [-h] [-D IDX_PATH] [-d IDX_NAME] [-r GEN_REF]
#                      [-o OUT_DIR] [-B] [-R] [-s]
#                      fastq_dir
#
# Process fastq files for CRISPR analysis
#
# positional arguments:
#   fastq_dir             The directory containing all the fastq files for
#                         analysis.
#
# optional arguments:
#   -h, --help            show this help message and exit
#   -D IDX_PATH, --index_path IDX_PATH
#                         Specify the reference genome index path. (Default =
#                         /Volumes/Genome/CRISPR/genome_idx)
#   -d IDX_NAME, --index_file IDX_NAME
#                         Specify the reference genome index name. (Default =
#                         'mini_genome')
#   -r GEN_REF, --genome_ref GEN_REF
#                         Specify the fasta file of the genome to use for stats.
#                         (Default = mini_genome.fa)
#   -o OUT_DIR, --directory OUT_DIR
#                         Specify an output directory. (Default = current
#                         working directory)
#   -B, --build_idx       Build the genome index file from specified fasta file
#   -R, --build_ref       Build the reference genome fasta file (will overwrite
#                         existing genome_ref)
#   -s, --align_stats     perform additional statistical analyses


# output:
# sorted .bam files for each sample

# input:
# directory of .fastq files with numbered samples denoted by *SX*.fastq (the default output from an Illumina MiSeq instrument)

import sys, re, csv, glob, os
import time, datetime # this is for determining how long the run will take
import argparse

import pysam

from genomepy import genematch, config


def define_arguments():
    # CLI parser to get variables that are more likely to change run to run:
    parser = argparse.ArgumentParser(description="Process fastq files for CRISPR analysis.\nCurrently works very well using defaults, requiring only the input file and path to fastq files. Ie:\ncrispralign.py -BRs -i <input_file> <fastq_dir>")

    ### input options ###
    # logging options:
    parser.add_argument("-q", "--quiet", action='store_true',default=False,
                        help="print fewer messages and output details")
    parser.add_argument("-f", "--output", type=str, default='lorf2.out',
                        help="specify the filename to save results to")
    parser.add_argument("--directory", type=str, default=os.getcwd(),
                        help="specify the directory to save results to")


    parser.add_argument("-D", "--index_path", type=str, dest="idx_path",
                        default="./genome_idx",
                        help="Specify the reference genome index path.\n(Default = /Volumes/Genome/CRISPR/genome_idx)")
    parser.add_argument("-d", "--index_file", type=str, dest="idx_name",
                        default="mini_genome",
                        help="Specify the reference genome index name.\n(Default = 'mini_genome')")
    parser.add_argument("-r", "--genome_ref", type=str, dest="gen_ref",
                        default="./mini_genome.fa",
                        help="Specify the fasta file of the genome to use for stats.\n(Default = mini_genome.fa)")

    parser.add_argument("fastq_dir", type=str, nargs=1,
                        help="The directory containing all the fastq files for analysis.")
    parser.add_argument("-i", "--input_file",  type=str, dest="in_file", default=False,
                        help="Specify an input file with site positions.\n(Format: scaffold    start    stop)")
    parser.add_argument("-b", "--buffer",  type=int, default=1000,
                        help="Specify the buffer (in bp) when building genome_ref")
    parser.add_argument("-B", "--build_idx", action="store_true",
                        help="Build the genome index file from specified fasta file")
    parser.add_argument("-R", "--build_ref", action="store_true",
                        help="Build the reference genome fasta file (will overwrite existing genome_ref)")
    parser.add_argument("-s", "--align_stats", action="store_true",
                        help="perform additional statistical analyses")
    parser.add_argument("-a", "--skip_alignment", action="store_true",
                        help="turns off read alignment")

    return parser


def verbal_system(arg_str):
    "Runs os.system, checks return value, and only prints something if value is not 0"
    r_val = os.system(arg_str)
    if r_val != 0:
        verbalise("R", arg_str, "did not work! Error:", r_val)
    return r_val

def build_index(gen_ref, idx_path, idx_name):
    if not os.path.isdir(idx_path):
        os.mkdir(idx_path)
    cmd_line = "gmap_build -d " + idx_name + " -D " + idx_path + " " + gen_ref
    print "Running: ", cmd_line
    print os.system(cmd_line)

def build_reference(gen_ref, pos_stats, buffer=1000):

    scaf_lengths = {}
    ref_h = open(gen_ref, 'w')
    # extract sequences from reference genome (identified under OGS):
    # adds a buffer of 1000bp to either end
    for scaffold in pos_stats:
        sequence = genematch.extractseq(scaffold, type='fa', OGS='/Volumes/antqueen/genomics/genomes/C.biroi/Cbir.assembly.v3.0', startpos=pos_stats[scaffold][0]-buffer, endpos=pos_stats[scaffold][1]+buffer)
        ref_h.write( ">%s\n%s\n" % (scaffold, sequence) )
        scaf_lengths[scaffold] = len(sequence)
    ref_h.close()
    return scaf_lengths

def read_input(in_file):
    """Reads scaffolds and positions of each site to be analysed, as well as the files from the
    miSeq to look at. Note the hash to indicate which lanes you are using.

    Format:
    scaffold1   startpos1   endpos1
    scaffold2   startpos2   endpos2

    #S1 S2 S3 S4

    """
    gene_pos = {}
    lanes = []

    in_file_h = open(in_file, 'rb')
    for line in in_file_h:
        if re.search("#", line):
            lane_files = line.replace("#","").split()
            for id in lane_files:
                lanes.append(id)
        elif re.search("\w",line) is None:  # ignores blank lines
            continue
        else:
            print line
            gene_pos[line.split()[0]] = (min([ int(x) for x in line.split()[1:] ]),
                                        max([ int(x) for x in line.split()[1:] ])
                                         )
            # WARNING: will not be able to have multiple sites in the same scaffold this way!

    in_file_h.close()

    return gene_pos, lanes

def align_reads(lanes, readfiles, idx_path, idx_name, out_dir, fq_dir):
    "Use gsnap to align fastq reads to indexed reference genome"
    # cycle through each of the samples
    cycles = len(lanes)
    t0 = time.time() # returns current time in seconds.

    for i, lane in enumerate(lanes):
        verbalise("M", "\n##### Running lane %s #####" % ( lane ))
        fileprefix = '/'.join([out_dir , lane])

        # this will be the file name of the bam file created
        bamFile= '.'.join([fileprefix,'sorted.bam']);
        verbalise("B", "BamFile = ", bamFile)

        # this finds all of the .fastq files associated with the sample number
        fq_path = fq_dir + '/*' + lane + '*.fastq*'
        reads = glob.glob(fq_path)

        # gsnap only works on unzipped files. Therefore need to check files are .fastq:
        for read in reads[:]:
            if re.search("gz$", read) is not None:
                verbal_system('gunzip -v ' + read)

        # reset the files (which should now all be .fastq):
        fq_path = fq_dir + '/*' + lane + '*.fastq'
        reads = glob.glob(fq_path)

        if len(reads) == 0:
            verbalise("R", "No fastq files found for this id. Skipping.")
            continue

        # first, let's align the reads with gsnap - set to allow indels of +/- 250bp
        # http://research-pub.gene.com/gmap/
        # you will have to copmile gmap/gsnap to allow for long reads like so:
        # ./configure MAX_READLENGTH=750; make; make install


        # gsnap alignment (specs for -d and -D can be supplied at CLI)
        print os.system('cat ' + " ".join(reads) + ' | gsnap -t 4 -y 250 -z 250 -A sam -D ' + idx_path + ' -d ' + idx_name + ' --input-buffer-size 10000 --sam-multiple-primaries > ' + fileprefix + '.sam')

        # convert from .sam file to .bam file
        verbal_system('samtools view -b -S ' + fileprefix + '.sam > ' + fileprefix + '.bam')

        # sort .bam file
        verbal_system('samtools sort ' + fileprefix + '.bam ' + fileprefix + '.sorted')

        # index bam file
        verbal_system('samtools index ' + fileprefix + '.sorted.bam')

        # calculate remaining time:
        j = i + 1
        t1 = time.time()
        t_taken = t1 - t0
        t_rate = t_taken / j
        t_remaining = datetime.timedelta(seconds=( (cycles - j) * t_rate ))
        verbalise("Y", "Sample %d of %d completed. Time remaining: %s" % (j, cycles, t_remaining))

def cigar_stats(bamFile, scaffold, startpos, endpos):
        # use the pysam package to open the .bam.indel file for writing data on the proportion
        out_file = open('.'.join([bamFile,'indel']), "a")

        # open bam file for reading
        bamFP = pysam.Samfile(bamFile, 'rb');

        # get all reads that align to our chromosome of interest
        site1 = bamFP.fetch(scaffold,startpos,endpos)

        # cycle through all reads that mapped to our chromosome
        c0count = 0
        for read in site1:
            if( not( read.is_unmapped ) ):   #if it's mapped
                cigarLine=read.cigar;


                # check the cigar string for insertions or deletions and write the length to our .indel file
                for (cigarType,cigarLength) in cigarLine:
                    try:
                        if(  cigarType == 0):
                            c0count += 1;
                        elif(cigarType == 1): #insertions
                            out_file.write(",".join(['d', str(cigarLength), scaffold + '\n']));
                        elif(cigarType == 2): #deletion
                            out_file.write(",".join(['d', str(-1*cigarLength), scaffold + '\n']));

                    except:
                        verbalise("R", "Problem!");

        verbalise("Y", "%-15s: %d" %  (scaffold, c0count))
        out_file.close()

def alignment_stats(lanes, readfiles, out_dir, gen_ref, startpos, endpos):
    """
    EVERYTHING HERE IS OPTIONAL - BASIC QUANTIFICATION OF INDEL SIZE/FREQUENCY
    """

    for lane in lanes:

        fileprefix = '/'.join([out_dir , lane])

        bamFile= '.'.join([fileprefix,'sorted.bam']);
        # use pysamstats to count the number of reads at each bp that have aligned with an insertion or deletion
        verbal_system('pysamstats --type variation ' + fileprefix + '.sorted.bam --fasta ' + gen_ref + '>' + fileprefix + '.variant.stats')
        verbalise("M", "Read depths for %s:"  % (lane))
        for scaffold in pos_stats:
            cigar_stats(bamFile, scaffold, startpos, endpos)

if __name__ == '__main__':

    parser = define_arguments()
    args = parser.parse_args()

    verbalise = config.check_verbose(not(args.quiet))
    logfile = config.create_log(args, outdir=args.directory, outname=args.output)

    gene_stats, lanes = read_input(args.in_file)


    if args.build_ref:
        verbalise("B",  "constructing reference file")
        scaf_lengths = build_reference(args.gen_ref, gene_stats, args.buffer)
    else:
        scaf_lengths = {}
        # calculate size of mini-genome:
        for scaf in gene_stats:
            scaf_lengths[scaf] = 0
            handle = open(args.gen_ref, 'rb')
            for line in handle:
                if line[0] == ">":
                    continue
                else:
                    scaf_lengths[scaf] += len(line.strip())

    if args.build_idx:
        verbalise("B",  "building index")
        build_index(args.gen_ref, args.idx_path, args.idx_name)

    # collect all .fq files:
    readfiles = os.listdir(args.fastq_dir[0]);

    # perform alignments!
    if not args.skip_alignment:
		verbalise("B", "Aligning reads\n", "#"  * 40)
		align_reads(lanes, readfiles, args.idx_path, args.idx_name,
		            args.directory, args.fastq_dir[0])
    if args.align_stats:
        pos_stats = {}
        for scaf in gene_stats:
            pos_stats[scaf] = (0, scaf_lengths[scaf])

        alignment_stats(lanes, readfiles, args.directory, args.gen_ref, pos_stats[scaffold][0], pos_stats[scaffold][1])

