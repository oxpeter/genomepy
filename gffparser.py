#!/usr/bin/env python
""" A module for parsing GFF files """

import sys
import cPickle
import re
import os

from BCBio import GFF
from Bio import SeqIO
from Bio.Alphabet import generic_dna, IUPAC
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram

import gff2bed2


##### FILE MANIPULATION #####

def pickle_jar(obj, fname):
    pklfile = ".".join([fname, 'pkl'])
    apicklefile = open(pklfile, 'wb')
    cPickle.dump(obj, apicklefile, -1)
    apicklefile.close()

def open_pickle_jar(fname):
    pklfile = ".".join([fname, 'pkl'])
    apicklefile = open(pklfile, 'rb')
    loaded_obj = cPickle.load(apicklefile)
    apicklefile.close()
    return loaded_obj

def pickle_gtypes(lineage_gtypes, indiv_gtypes, fname='pickle'):
    # try to pickle the gtypes!
    pklfile = ".".join([fname, 'aGtypes.pkl'])
    apicklefile = open(pklfile, 'wb')
    cPickle.dump(lineage_gtypes, apicklefile, -1)
    apicklefile.close()

    pklfile = ".".join([fname, 'iGtypes.pkl'])
    ipicklefile = open(pklfile, 'wb')
    cPickle.dump(indiv_gtypes, ipicklefile, -1)
    ipicklefile.close()

def unpickle_gtypes(fname='pickle'):
    # and now to unpickle!
    pklfile = ".".join([fname, 'aGtypes.pkl'])
    apicklefile = open(pklfile, 'rb')
    lineage_gtypes = cPickle.load(apicklefile)
    apicklefile.close()

    pklfile = ".".join([fname, 'iGtypes.pkl'])
    ipicklefile = open(pklfile, 'rb')
    indiv_gtypes = cPickle.load(ipicklefile)
    ipicklefile.close()

    return lineage_gtypes, indiv_gtypes


##### FILE CONVERSIONS #####

def gff2gtf(gff_file):
    "converts a gff file to gtf"
    gtf_file = os.path.splitext(gff_file)[0] + ".gtf"
    gtf_h = open(gtf_file, 'w')
    gff_h = open(gff_file, 'rb')

    for line in gff_h:
        cds_patt = '(.*)Parent=([A-Za-z0-9_\(\)]*)'
        parentline = re.search(cds_patt, line)
        if parentline is not None:
            newline = parentline.groups()[0] + '\tgene_id "' + parentline.groups()[1] + '"; transcript_id "' + parentline.groups()[1] + ';\n'
            gtf_h.write(newline)
    gtf_h.close()
    gff_h.close()

def gtf2gff(gtf_file):
    "converts a gtf file to gff using perl script"
    gff_file = os.path.splitext(gtf_file)[0] + ".gff"
    cmd = "perl /Users/POxley/scripts/Genome_Analysis/gtf2gff3.pl " + gtf_file + " > " + gff_file
    os.system(cmd)

def make_bed(gff_file):
    """ Takes a gff_file and creates a bed format including the scaffolds/contigs that
    do not have any features associated with them"""

    gff_fh = open(gff_file, 'rU')

    # get transcript annotation
    tinfo, einfo = gff2bed2.ParseAnno(gff_fh)
    # write into bed format
    scaf_dict = gff2bed2.WriteBED(tinfo, einfo, return_obj=True, gff_name=gff_file)

    # create dictionary of scaffolds not annotated, and their sizes
    missing_dict = {}
    assembly_handle = open("/Volumes/Genome/Cbir.assembly.v3.0_singleline.fa", 'r')
    for line in assembly_handle:
        def_search = re.search('>(lcl\|[a-zC0-9]*)', line)
        line_len = len(line)
        if def_search is not None:
            scaf_id = def_search.groups()[0]
        else:
            missing_dict[scaf_id] = line_len

    # append these onto the bed file, using the first 3 columns
    bed_name = gff_file + ".bed"
    bed_handle = open(bed_name, 'a')
    for scaf in missing_dict:
        bed_handle.write( "%s\t%d\t%d\tNonGene\t.\t+\t%d\t%d\t0\t1\t%d\t%d\n" % (scaf, 0, missing_dict[scaf], 0, missing_dict[scaf], 0, missing_dict[scaf]) )
    bed_handle.close()


##### GFF FILE MANIPULATION #####

def assemble_dict(in_file="/Volumes/Genome/armyant.OGS.V1.8.6.gff", in_seq_file="/Volumes/Genome/Cbir.assembly.v3.0.fa", features_only=False):
    """ Takes a gff file and the genome assembly and creates a SeqRecord dictionary.

    INPUT:
    in_file = path/to/file.gff
    in_seq_file = path/to/assembly.file.fa

    OUTPUT:
    gff_dict['scaffold'] = SeqRecord_obj
    """
    in_seq_handle = open(in_seq_file)
    seq_dict = SeqIO.to_dict(SeqIO.parse(in_seq_handle, "fasta"))
    in_seq_handle.close()

    gff_dict = {}

    in_handle = open(in_file)

    if features_only:
        for rec in GFF.parse(in_handle, base_dict=seq_dict):
            if len(rec.features)>=1:
                gff_dict[rec.id] = rec
    else:
        for rec in GFF.parse(in_handle, base_dict=seq_dict):
            gff_dict[rec.id] = rec

    return gff_dict

def create_gff(gff_dict, out_file):
    """ Given the gff dictionary, this function creates a new GFF file"""
    print "Writing to", out_file

    out_handle = open(out_file, 'w')
    rec_list = []

    for rec in gff_dict:
        rec_list.append(gff_dict[rec])

    GFF.write(rec_list, out_handle)

    out_handle.close()

def isolate_mRNA(gff_p):
    """takes a gff file and removes the gene header, allowing conversion and viewing of
    all isoforms on JBrowse"""

    newfile =  os.path.splitext(gff_p)[0] + ".mRNA.gff"

    gff_h = open(gff_p, 'rb')
    newfile_h = open(newfile, 'w')

    for line in gff_h:
        if line.split()[2] == "mRNA":
            defsearch = re.search("(.*);Parent=.+", line)
            newfile_h.write(defsearch.groups()[0] + "\n")
        elif line.split()[2] == "exon":
            newfile_h.write(line)

    gff_h.close()
    newfile_h.close()
    return newfile


##### SEQRECORD FUNCTIONS #####

def count_genes(gff_dict):
    "Returns the number of features (ie genes) in a SeqRecord obj"
    count = 0
    for rec in gff_dict:
        count += len(gff_dict[rec].features)
    return count

def show_features(gff_dict, scaf):
    "Reports the features found in a given SeqRecord for the given scaffold"
    for feature in gff_dict[scaf].features:
        print "Gene:", feature.qualifiers['ID'][0],
        try:
            print "%s\t%s" % (feature.qualifiers['Name'][0], feature.location)
        except KeyError:
            print feature.location

def snp_in_gene(scaf, pos, gff_dict):
    """returns salient features of a given SNP.

    INPUT:
    scaf = scaffold SNP is located on
    pos = position in scaffold of SNP (NB uses python base 0!!!)
    gff_dict = dictionary of SeqRecords (must contain the given scaffold as a key)

    OUTPUT:
    cds_pos = the distance into the coding sequence of the SNP
    exon =  the exon the SNP is located on
    codon = the nt sequence of the codon the SNP is on
    frame = the position the SNP is in within the codon (position 1,2 or 3)
    ref_nt = the nucleotide in the reference assembly at the SNP position
    """



    for feat in gff_dict[scaf].features:
        if pos in feat:

            # print gene details for match:
            #print feat.qualifiers['ID'][0],
            try:
                gene_name = feat.qualifiers['Name'][0]
                gene_id   = feat.qualifiers['ID'][0]
            except KeyError:
                gene_name = None
                gene_id   = None
            # collect info about frame posn, distance into gene etc:
            before_len = 0
            before_cnt = 0
            for subfeat in feat.sub_features:
                if pos in subfeat.location:
                    if feat.strand == 1:
                        exon_pos = pos - subfeat.location.start + 1
                        continue
                    else:
                        exon_pos = subfeat.location.end - pos
                        continue
                if feat.strand == 1 and subfeat.location.start < pos:
                    before_len += len(subfeat.location)
                    before_cnt += 1
                elif feat.strand == -1 and subfeat.location.start > pos:
                    before_len += len(subfeat.location)
                    before_cnt += 1

            cds_pos = exon_pos + before_len
            frame =  {1:1,2:2,0:3}[cds_pos % 3]
            seq_start = cds_pos - frame
            seq_end = seq_start + 3
            exon = before_cnt +  1
            #print "SNP position: %d\texon: %d\t" % (cds_pos, exon)


            rec_seq = gff_dict[scaf].seq
            if feat.strand == 1:
                codon = feat.extract(rec_seq)[seq_start:seq_end]
                ref_nt = rec_seq[pos:pos+1]
                #print "Codon: %s\t Codon Position: %d\t Ref. nt: %s\t" % (codon, frame, ref_nt)
            else:
                feat_seq = get_sequence(scaf, gff_dict, feat)
                codon = feat_seq[seq_start:seq_end]
                ref_nt = rec_seq[pos:pos+1].reverse_complement()
                #print "Codon: %s\t Codon Position: %d\t Ref. nt: %s\t" % (codon, frame, ref_nt)
            #print "SNP lies", feat.location.end - pos + 1, "from the start of the gene (including introns)."
            SNP_dict = {"gene_id":gene_id, "gene_name":gene_name, "cds_pos":cds_pos, "exon":exon, "codon":codon, 'codon_str':str(codon), "frame":frame, "ref_nt":ref_nt, 'ref_nt_str':str(ref_nt)}
            break
    else:
        SNP_dict = {"gene_id":None, "gene_name":None, "cds_pos":None, "exon":None, "codon":None, 'codon_str':None, "frame":None, "ref_nt":None, 'ref_nt_str':None}
        # It may be necessary to report -ve results in the future. Perhaps something like:
        #else:
        #    SNP_dict[feat] = (None, None, None, None, None)

    return SNP_dict

def get_sequence(scaf, gff_dict, feature):
    """ There appears to be a bug in the .extract function of features, such that genes
    in the reverse strand have their sequences reverse complemented correctly, but the
    exons can appear in the wrong order. This is my attempt to fix this."""

    sequence = ""
    if feature.strand == 1:
        sequence = feature.extract(gff_dict[scaf].seq)
    else:
        subseq_dict = {}
        for subfeature in feature.sub_features:
            subseq_dict[int(subfeature.location.start)] = subfeature.extract(gff_dict[scaf].seq)
        for exon in sorted(subseq_dict.keys(), reverse=True):
            sequence += subseq_dict[exon]

    sequence.alphabet = IUPAC.ambiguous_dna
    return sequence

##### INITIATION & MISC #####

def prime_variables():
    return "/Volumes/Genome/armyant.OGS.V1.8.5.gff", "/Volumes/Genome/Cbir.assembly.v3.0.fa"

def mutate_codon(old_seq, frame, new_nt):
    """Given a Bio.seq.seq sequence, a new nucleotide, and its position in the sequence,
    mutate_codon returns a new Bio.seq.seq sequence with the new nt substituted at the
    given position.
    """
    seq_str = str(old_seq)[0:frame-1] + new_nt + str(old_seq)[frame:]
    new_seq = Seq(seq_str, generic_dna)
    return new_seq

def draw_gene(feature):
    "Produces a visual representation of a gene! (Hopefully!)"
    gd_diagram = GenomeDiagram.Diagram(feature.qualifiers['ID'][0])
    gd_track_for_features = gd_diagram.new_track(1, name="Exons")
    gd_feature_set = gd_track_for_features.new_set()

    #for subfeature in feature.sub_features:
    color = colors.blue
    gd_feature_set.add_feature(feature, color=color, label=True)

    gd_diagram.draw(format='linear', orientation='landscape', pagesize='A4', fragments=1, start=0, end=len(feature))
    gd_diagram.write("test_gene.img.pdf", "PDF")

def test():
    """ The basic commands I have been running of '__main__'
    """

    in_file = "/Volumes/Genome/armyant.OGS.V1.8.5.gff"
    out_file = "/Volumes/Genome/armyant.OGS.V1.8.5.gff"
    in_seq_file = "/Volumes/Genome/Cbir.assembly.v3.0.fa"


    # set scaffold and SNP to check:
    cmdargs = sys.argv
    if len(cmdargs) == 3:
        cmd, scaf, cpos = sys.argv
        pos = int(cpos) - 1
    else:
        scaf = "scaffold61"
        pos = 40905
    print "Position has been converted to python base 0!!!"

    # create SeqRecord file from GFF and Assembly:
    gff_dict = assemble_dict(in_file, in_seq_file)

    # count the number of genes in the SeqRecord dictionary:
    gene_num = count_genes(gff_dict)
    print "there are", gene_num, "gene features in this dictionary."

    #for feature in gff_dict[scaf].features:
    #    print "%s: %s" % (feature.qualifiers['ID'][0], feature.location)

    print "\nResults for SNP %s:%d\n" % (scaf, pos)
    result_dict = snp_in_gene(scaf, pos, gff_dict)

    for feature in result_dict:
        cds_pos, exon, cseq, frame, ref_nt = result_dict[feature]
        prot = cseq.translate()

        print feature.qualifiers['ID'][0],
        try:
            print feature.qualifiers['Name'][0], feature.location
        except KeyError:
            print feature.location
        print "SNP position: %d\texon: %d\t" % (cds_pos, exon)
        print "Codon: %s (%s) Codon Position: %d\t Ref. nt: %s\t" % (cseq, prot, frame, ref_nt)

        new_seq = mutate_codon(cseq, frame, "A")
        print new_seq, new_seq.translate()

def create_polymorphism_files():
    print "Assembling SeqRecord from gff file..."
    gff_dict = assemble_dict()
    snp_dict = {}

    print "Unpickling NSL_PILSIP genotypes..."
    aGtypes, iGtypes = unpickle_gtypes('/Volumes/Genome/RAD-Tags/NSL_PILSIP.RD_15')

    print "Scanning for polymorphisms..."
    for scaf, pos in aGtypes['LINEC']:
        SNP1 = aGtypes['LINEC'][(scaf, pos)][0]
        SNP2 = aGtypes['LINEC'][(scaf, pos)][1]

        # have to remove indels at this point :(
        if SNP1 == "I" or SNP2 == "I" or SNP1 == "J" or SNP2 == "J":
            continue

        results_dict = snp_in_gene(scaf, pos, gff_dict)
        for feature in results_dict:
            cds_pos, exon, ref_cdn, frame, ref_nt = results_dict[feature]
            ref_ptn = ref_cdn.translate()
            snp1_cdn = mutate_codon(ref_cdn, frame, SNP1)
            snp2_cdn = mutate_codon(ref_cdn, frame, SNP2)
            snp1_ptn = snp1_cdn.translate()
            snp2_ptn = snp2_cdn.translate()

            snp_dict[(scaf, pos, feature.qualifiers['ID'][0], len(feature), cds_pos)] = ( (ref_nt, ref_cdn, ref_ptn), (SNP1, snp1_cdn, snp1_ptn), (SNP2, snp2_cdn, snp2_ptn) )

    print len(snp_dict), "results found"
    print "Pickling results..."
    pickle_jar(snp_dict, "/Volumes/Genome/RAD-Tags/LineC_polymorphisms")

    print "writing results to tab delimited file..."
    wobj = open( "/Volumes/Genome/RAD-Tags/LineC_polymorphisms.list", 'w' )
    wobj.write("Scaffold\tPos\tAccession\tCDS_pos (%)\tPolymorphism1\tPolymorphism2\n")
    for scaf, pos, id, cds_len, cds_pos in snp_dict:
        if str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1]) == str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][1][1]) and str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1]) == str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][2][1]):
            pass
        elif str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1]) == str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][1][1]):
            wobj.write( "%s\t%d\t%s\t%d (%.2f%%)\t%s>%s (%s>%s)\n" % (scaf, pos, id, cds_pos, 100.0 * cds_pos / cds_len, snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][2][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][2], snp_dict[(scaf, pos, id, cds_len, cds_pos)][2][2])  )
        elif str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1]) == str(snp_dict[(scaf, pos, id, cds_len, cds_pos)][2][1]):
            wobj.write( "%s\t%d\t%s\t%d (%.2f%%)\t%s>%s (%s>%s)\n" % (scaf, pos, id, cds_pos, 100.0 * cds_pos / cds_len, snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][1][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][2], snp_dict[(scaf, pos, id, cds_len, cds_pos)][1][2])  )
        else:
            wobj.write( "%s\t%d\t%s\t%d (%.2f%%)\t%s>%s (%s>%s)\t%s>%s (%s>%s)\n" % (scaf, pos, id, cds_pos, 100.0 * cds_pos / cds_len, snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][1][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][2], snp_dict[(scaf, pos, id, cds_len, cds_pos)][1][2], snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][2][1], snp_dict[(scaf, pos, id, cds_len, cds_pos)][0][2], snp_dict[(scaf, pos, id, cds_len, cds_pos)][2][2])  )
    wobj.close()

    print "All done!"

def find_next_isoform(gene_id, gff_file="/Volumes/Genome/armyant.OGS.V1.8.5.gff"):
    "returns the next number in the list of isoforms for a given gene"
    isoforms = []

    gff_handle = open(gff_file, 'rb')
    for line in gff_handle:
        gid = re.search(gene_id, line)
        giso = re.search('-R([A-Z]+)', line)
        if gid is not None and giso is not None:
            isoforms.append(giso.groups()[0])

    last_iso = isoforms.sorted()[-1]
    # research how to determine next letter in Alphabet
    return last_iso

##### FASTA FILE FUNCTIONS #####

def count_chromosomes(assembly="/Volumes/Genome/Cbir.assembly.v3.0_gi.fa"):
    "uses fasta file of genome to produce tab delimited list of chromosomes and their size"
    newfile = os.path.splitext(assembly)[0] + ".chrom.list"
    newfile_h = open(newfile, 'w')
    assembly_handle = open(assembly, 'rb')
    csize = 0
    chrom = None
    for line in assembly_handle:
        if re.search('>', line) is not None:
            if chrom is not None:
                newfile_h.write( "%s\t%d\n" % (chrom, csize) )
            chrom = line.split()[0][1:]
            csize = 0
        else:
            csize += len(line.strip())
    newfile_h.write( "%s\t%d\n" % (chrom, csize) )

    newfile_h.close()
    assembly_handle.close()

##### CUFFLINKS ASSOCIATED FUNCTIONS #####

def find_cuffjoined(gtf_file):
    "looks for genes joined together by cuffmerge"
    allid = {}
    file_h = open(gtf_file)
    for line in file_h:
        gid = re.search("gene_id \"([A-Za-z0-9_]*)\";", line)
        nearest = re.search("nearest_ref \"([A-Za-z0-9_\(\)]*)\";", line)
        if gid is not None and nearest is not None:
            if gid.groups()[0] not in allid:
                allid[gid.groups()[0]] = [nearest.groups()[0]]
            elif nearest.groups()[0] not in allid[gid.groups()[0]]:
                allid[gid.groups()[0]].append(nearest.groups()[0])
    file_h.close()

    dblid = {}
    for gid in allid:
        if len(allid[gid]) > 1:
            dblid[gid] = allid[gid]

    return dblid

def separate_cuffjoined(gtf_file, dblid):
    """finds all loci in the gtf_file that are joined according to dblid, and alters the
    gene_id to create new gene_ids for all joined genes (effectively separating them again)
    """
    # set variables:
    newfile =  os.path.splitext(gtf_file)[0] + ".separated.gtf"
    count = 0

    # open files:
    newfile_h = open(newfile, 'w')
    file_h = open(gtf_file)

    for line in file_h:
        gid = re.search("gene_id \"([A-Za-z0-9_]*)\";", line)
        nearest = re.search("nearest_ref \"([A-Za-z0-9_\(\)]*)\";", line)
        if gid.groups()[0] in dblid:
            if nearest is None:
                newid = gid.groups()[0] + "novel"
                count += 1
            else:
                newid = gid.groups()[0] + nearest.groups()[0]
        else:
            newid = gid.groups()[0]
        newfile_h.write( line.replace(gid.groups()[0], newid) )

    file_h.close()
    newfile_h.close()
    print count
    return newfile


if __name__ == '__main__':

    isodict = assemble_dict(in_file="/Volumes/Genome/transcriptomes/BroodSwap/controls/C16/Cuffmerge/merged.filtered.separated.gff", in_seq_file="/Volumes/Genome/Cbir.assembly.v3.0_gi.fa")
    for feature in isodict['lcl|scaffold581'].features:
        print "ID: %s\tName: %s\tType: %s\tStrand: %r" % (feature.qualifiers['ID'][0], feature.qualifiers['Name'][0], feature.type, feature.strand)
        for subfeat in feature.sub_features:
            if subfeat.qualifiers['Name'][0] == 'Cbir_12163':
                print "ID: %s\tName: %s\tType: %s\tParent: %s\tStrand: %r" % (subfeat.qualifiers['ID'][0], subfeat.qualifiers['Name'][0], subfeat.type, subfeat.qualifiers['Parent'], subfeat.strand)
                mRNA_seq = get_sequence('lcl|scaffold581', isodict, subfeat)
                print "Seq: %s\tLength: %d\tProt: %s\tLength: %d" % (mRNA_seq[0:10], len(mRNA_seq), mRNA_seq.translate()[0:10], len(mRNA_seq.translate()))
                longestORF = 0
                for frame in range(3):
                    for prot in mRNA_seq[frame:].translate().split('*'):
                        if len(prot) > longestORF:
                            longestORF = len(prot)
                            longestSeq = prot
                print "Longest Prot: %s\tLength: %d" % (longestSeq[0:10], len(longestSeq))

    #print subfeat.qualifiers


    """
    gtf_file = "/Volumes/Genome/transcriptomes/BroodSwap/controls/C16/Cuffmerge/merged.filtered.gtf"
    dblid = find_cuffjoined(gtf_file)
    for gene in dblid:
        print "%s\t%s" % (gene, dblid[gene])
    print len(dblid)

    sep_gtf = separate_cuffjoined(gtf_file, dblid)
    gtf2gff(sep_gtf)
    """

    """
    gff_dict = assemble_dict("/Volumes/Genome/armyant.OGS.V1.8.5_lcl.gff")
    print "GFF dictionary parsed. Writing new GFF file..."
    create_gff(gff_dict, "/Volumes/Genome/armyant.OGS.V1.8.5_recreated.gff")
    """

    """
    in_file = "/Users/POxley/Documents/Analysis/Data/Transcriptome/alternate_splicing/iReckon/test_files/COI.gb"
    out_file = "/Users/POxley/Documents/Analysis/Data/Transcriptome/alternate_splicing/iReckon/test_files/COI.gff"
    in_handle = open(in_file)
    out_handle = open(out_file, "w")

    seq_obj = SeqIO.parse(in_handle, "genbank")
    print type(seq_obj)
    for thing in seq_obj:
        print thing
        GFF.write(thing, out_handle)

    in_handle.close()
    out_handle.close()

            """
