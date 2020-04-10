from __future__ import print_function
from collections import defaultdict
import intervaltree
import sys
import argparse
import pyfaidx
from math import floor

"""
Reads in a FASTA reference,
a VCF containing isolate variants,
and a tidied-GFF file
and returns a new file with four columns:
SEQ_ID  POS  BASE  MOLECULE FEATURE

where SEQ_ID is the FASTA identifier of the sequence the variant is in,
POS is the 1-based position, BASE is the DNA base or amino acid at that location,
MOLECULE is one of "DNA" or "PROTEIN",
and FEATURE is the GFF feature that the variant falls within (if any)
"""



base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

codon_to_AA = defaultdict(str, { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
    })

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reference", dest="fasta",
            help="A FASTA reference genome.", type=str, required=True, default=None)
    parser.add_argument("-i", "--variants", dest="variants",
            help="A VCF (or the first 7 tab-separated columns) of variants to translate.", default=None)
    parser.add_argument("-g", "--tidy-gff", dest="gff",
            help="A tidified GFF file to find genomic start locations of peptides.", required=True, type=str, default=None)
    parser.add_argument("-s", "--sample", default=None, dest="sample",
            help="A sample name to use in Tidy output.", required=False)

    return parser.parse_args()


"""
Takes the reference sequence, the position within it, the ref / alt alleles,
the peptide name, the codon number within the protein, and the frame (a number between 0-2 indicating
where within the codon we are) and returns the modified amino acid.
"""
def translate_variant_site(seq, position, alt, strand, phase, feature_name):
    if strand == "-":
        print("Reverse-complementing seq to match - strand", file=sys.stderr)
        seq = "".join([base_complement[i] for i in seq[::-1]])

    return

def tidy_dna_site(seq_id, position, ref, alt, feature_name):
    MOL = "DNA"
    if len(ref) > len(alt):
        alt = "DEL"
    elif len(alt) > len(ref):
        alt = "INS"
    
    return "\t".join([seq_id, str(position), alt, MOL, feature_name])


def translate_seq(seq):
    pep = []
    for i in range(0, len(seq)+1, 3):
        codon = seq[i:i+3]
        aa = codon_to_AA[codon]
        pep.append(aa)
        #print(codon, aa)
    return "".join(pep).strip()

def nuc_to_aa_index(pos, feature_start):
    return floor((pos - feature_start) / 3)

def nuc_to_aa_phase(pos, feature_start):
    return (pos -feature_start) % 3

if __name__ == "__main__":
    
    #pep = ""
    #for i in range(0, len(sys.argv[1])+1, 3):
    #    pep += codon_to_AA[ sys.argv[1][i:i+3] ]
    #print(pep)

    args = parse_args()

    fasta = pyfaidx.Fasta(args.fasta)

    #ref	name	gene	start	end	strand	phase	feature_type
    gff_features = intervaltree.IntervalTree()
    feature_peptide_cache = defaultdict(str)
    with open(args.gff, "r") as gfi:
        for line in gfi:
            if not line.startswith("#"):
                tokens = line.strip().split("\t")
                if tokens[7] == "CDS":
                    tokens[3] = int(tokens[3])
                    tokens[4] = int(tokens[4])
                    tokens[6] = int(tokens[6])
                    gff_features.addi(tokens[3], tokens[4], tuple(tokens))
    header = "SEQID\tpos\tbase\tmolecule\tfeature"
    print(header)
    if args.variants is not None:
        variant_list = []
        with open(args.variants, "r") as ifi:
            for line in ifi:
                if not line.startswith("#"):
                    line = line.strip()
                    tokens = line.split("\t")
                    """
                    A 4-tuple, containing the SEQID, POS, REF and ALT fields.
                    """
                    v = (tokens[0], int(tokens[1]), tokens[3], tokens[4])
                    if gff_features.overlaps(v[1]):
                        feats = gff_features.at(v[1])
                        for feat in feats:
                            feat_data = feat.data
                            feat_id = "_".join([str(s) for s in feat_data])
                            feat_ref = feat_data[0]
                            gene_name = feat_data[2]
                            feat_start = feat_data[3]
                            feat_end = feat_data[4]
                            feat_strand = feat_data[5]
                            feat_phase = feat_data[6]
                            seq = fasta[feat_ref][feat_start-1-feat_phase:feat_end-1].seq
                            codon_number = nuc_to_aa_index(v[1], feat_start)
                            codon_base = nuc_to_aa_phase(v[1], feat_start)
                            if feat_id not in feature_peptide_cache:
                                feature_peptide_cache[feat_id] = translate_seq(seq)
                            aa_seq = feature_peptide_cache[feat_id]
                            print("Located CDS variant. Translating.", feat_id, "codon:", codon_number, "base", codon_base, v, file=sys.stderr)
                            print(tidy_dna_site(v[0], v[1], v[2], v[3], gene_name))
                            #translate_variant_site(seq, v[1], v[3], feat_strand, feat_phase)

                    else:
                        print(tidy_dna_site(v[0], v[1], v[2], v[3], "Noncoding"))
                        print("Non-coding var", v, file=sys.stderr)

