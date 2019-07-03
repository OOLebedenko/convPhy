import os
import pandas as pd
from collections import defaultdict
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
import Bio.Data.CodonTable
from Bio import SeqIO

path_to_genbank = "/home/olebedenko/mtb_ref/sequence.gb"
gb_record = SeqIO.read(open(path_to_genbank, "r"), "genbank")


def annotate_snp(snps_path, path_to_snps_out_csv):
    snps = open(snps_path, "r")
    snps_dir = snps_path.rsplit("/", 1)[0]
    annotated_snps = pd.DataFrame()

    for ind, line in enumerate(snps):
        position = line.split("\t")[1]
        alt = line.split("\t")[2]
        ref = line.split("\t")[3]
        variant_info = annotate(gb_record, ref, position, alt)
        temp = pd.DataFrame([variant_info])
        annotated_snps = pd.concat([annotated_snps, temp])
    annotated_snps.to_csv(path_to_snps_out_csv, index=False)


def effect(codon1, alt_codon):
    standard_table = Bio.Data.CodonTable.standard_dna_table
    all = standard_table.forward_table
    all['TAA'] = 'Stop'
    all['TGA'] = 'Stop'
    all['TAG'] = 'Stop'
    all['ATG'] = 'M'
    try:
        acid = all[str(codon1)]
        if all[str(codon1)] == all[str(alt_codon)]:
            effect = "Syn"
        else:
            effect = "Nonsyn"
        acid += "/" + all[str(alt_codon)]
    except:
        effect, acid = 'unk', 'unk'
    return effect, acid


def alt_complement(alt):
    s = ""
    if alt == 'C':
        s = 'G'
    elif alt == 'G':
        s = 'C'
    elif alt == 'T':
        s = 'A'
    elif alt == 'A':
        s = 'T'
    return s


def codons_def(seq, pos, start, end, alt, strand):

    if strand == 1:
        if (pos - start) % 3 == 0:
            codon1 = str(seq[pos - start:pos - start + 3])
            codon2 = alt + codon1[1] + codon1[2]
        elif (pos - start) % 3 == 1:
            codon1 = str(seq[pos - start - 1:pos - start + 2])

            try:
                codon2 = codon1[0] + alt + codon1[2]
            except:
                codon2 = codon1[0] + alt
                print(pos - start) % 3, codon1, seq, pos, pos - start, start - end, '[', start, end, ']', alt, strand
        elif (pos - start) % 3 == 2:
            codon1 = str(seq[pos - start - 2:pos - start + 1])
            codon2 = codon1[0] + codon1[1] + alt
    elif strand == -1:
        var1 = alt_complement(alt)
        rev_seq = seq.reverse_complement()
        acid = rev_seq.translate(table=11, to_stop=True)
        codon_no = int(round((end - pos - 1) / 3 + 0.5))
        if (end - pos - 1) % 3 == 0:
            codon1 = str(rev_seq[end - pos - 1:end - pos + 2])
            codon2 = var1 + str(rev_seq[end - pos]) + str(rev_seq[end - pos + 1])
        elif (end - pos - 1) % 3 == 1:
            codon1 = str(rev_seq[end - pos - 2:end - pos + 1])
            codon2 = str(rev_seq[end - pos - 2]) + var1 + str(rev_seq[end - pos])
        elif (end - pos - 1) % 3 == 2:
            codon1 = str(rev_seq[end - pos - 3:end - pos])
            codon2 = str(rev_seq[end - pos - 3]) + str(rev_seq[end - pos - 2]) + var1
    return codon1, codon2


def annotate(gb_record, ref, position, alt):
    variant_info = {}
    flag = 0
    pos = int(position) - 1
    variant_info['position'] = position
    variant_info['ref'] = ref
    variant_info['alt'] = alt
    for feature in gb_record.features:
        start = feature.location.nofuzzy_start
        end = feature.location.nofuzzy_end
        if start <= pos and pos <= end:

            feature_keys = feature.qualifiers.keys()
            strand = feature.location.strand

            if feature.type != 'source' and feature.type != 'gene': # source ase used only in first feature line in genbank file for H37Rv
                if feature.type == 'CDS':
                    flag = 1 # flag = 0 in two cases First, if not gene and second not in pos (one flag for two similar cases)
                    feature_ref_sequence = gb_record.seq[start:end]
                    check_constant_reference = feature_ref_sequence[pos - start:pos - start + len(ref)] == ref
                    codon_number = int(round((pos - start) / 3 + 0.5))
                    pos_in_gene = pos - start

                    if check_constant_reference:
                        for key in feature_keys:
                            if key == 'locus_tag':
                                locus_tag = feature.qualifiers[key][0]
                                variant_info['locus_tag'] = locus_tag
                            elif key == 'gene':
                                gene = feature.qualifiers[key][0].replace("'", "")
                                variant_info['gene'] = gene

                        if len(alt) == 1 and len(ref) == 1:
                            codon1, codon2 = codons_def(feature_ref_sequence, pos, start, end, alt, strand)
                            eff, acid = effect(codon1, codon2)
                        else:
                            codon1, codon2 = ('.', '.')
                            eff, acid = ('.', '.')
                        if strand == -1:
                            codon_number = int(round((end - pos - 1) / 3 + 0.5))
                            pos_in_gene = end - pos - 1
                        variant_info['codon_number'] = codon_number
                        variant_info['pos_in_gene'] = pos_in_gene
                        variant_info['codon'] = str(codon1) + "/" + str(codon2)
                        variant_info['effect'] = eff
                        variant_info['amin_acid'] = acid
                        variant_info.update({'ref': ref, 'alt': alt, 'position': position})
                        return variant_info
                    else:
                        replaced_ref = 'Nucleotide position {pos} in the genome sequence has been corrected {ref}:{new_ref}'.\
                                        format(pos=pos, ref=ref, new_ref=feature_ref_sequence[pos - start:pos - start + len(ref)])
                        return replaced_ref, ref, position, alt, feature_ref_sequence[pos - start:pos - start + len(ref)]
                else:
                    return "{feature_type}: RNA or exactly unknown about transation this sequence type".format(feature_type=feature.type), ref, position, alt
        if start >= pos:
            break
    if flag == 0:
        print(feature)
        return {'ref': ref, 'alt': alt, 'position': position, 'effect': 'intergenic'}

if __name__ == '__main__':
    # path_to_genbank="/home/olebedenko/mtb_ref/sequence.gb"
    # gb_record = SeqIO.read(open(path_to_genbank,"r"), "genbank")\

    print annotate(gb_record, "A",467,"G")
    # standard_table = Bio.Data.CodonTable.standard_dna_table
    # all = standard_table.forward_table
    # test = []
    # for ind, feature in enumerate(gb_record.features):
    #         if 0 < ind:
    #             start = feature.location.nofuzzy_start
    #             end = feature.location.nofuzzy_end
    #             strand = feature.location.strand
    #             feature_keys = feature.qualifiers.keys()
    #             test.append(feature.type)

    # print(set(test))
