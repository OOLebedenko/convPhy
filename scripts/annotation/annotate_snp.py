import os
import pandas as pd
from collections import defaultdict
from Bio.SeqFeature import FeatureLocation
from Bio import SeqIO
import Bio.Data.CodonTable
from Bio import SeqIO

#path_to_genbank = "/home/olebedenko/mtb_ref/sequence.gb"


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
    standard_table = Bio.Data.CodonTable.standard_dna_table
    all = standard_table.forward_table
    all['TAA'] = 'Stop'
    all['TGA'] = 'Stop'
    all['TAG'] = 'Stop'
    all['ATG'] = 'M'
    # try:
    if strand == 1:
        if (pos - start) % 3 == 0:
            codon1 = str(seq[pos - start:pos - start + 3])
            # print (pos-start)%3, codon1,seq, pos-start, start-end, pos,'[', start,end,']', alt, strand
            codon2 = alt + codon1[1] + codon1[2]
        # print gen[i][0], gen[i][3], codon1+"/"+codon2+"["+codon21+"]",quality,dp
        elif (pos - start) % 3 == 1:
            codon1 = str(seq[pos - start - 1:pos - start + 2])

            try:
                codon2 = codon1[0] + alt + codon1[2]
            except:
                codon2 = codon1[0] + alt
                print(pos - start) % 3, codon1, seq, pos, pos - start, start - end, '[', start, end, ']', alt, strand
            # print gen[i][0], gen[i][3], codon1+"/"+codon2+"["+codon21+"]",quality,dp
        elif (pos - start) % 3 == 2:
            codon1 = str(seq[pos - start - 2:pos - start + 1])
            # print (pos-start)%3,codon1,seq, pos, pos-start, start-end,'[',start,end,']', alt, strand
            codon2 = codon1[0] + codon1[1] + alt
    elif strand == -1:
        var1 = alt_complement(alt)
        rev_seq = seq.reverse_complement()
        acid = rev_seq.translate(table=11, to_stop=True)
        codon_no = int(round((end - pos - 1) / 3 + 0.5))
        if (end - pos - 1) % 3 == 0:
            codon1 = str(rev_seq[end - pos - 1:end - pos + 2])
            # print codon1,seq, pos, end-pos, start-end,'[',start,end,']', alt, strand
            codon2 = var1 + str(rev_seq[end - pos]) + str(rev_seq[end - pos + 1])
        elif (end - pos - 1) % 3 == 1:
            codon1 = str(rev_seq[end - pos - 2:end - pos + 1])
            # print codon1,seq, pos, end-pos, start-end,'[',start,end,']', alt, strand
            codon2 = str(rev_seq[end - pos - 2]) + var1 + str(rev_seq[end - pos])
        elif (end - pos - 1) % 3 == 2:
            codon1 = str(rev_seq[end - pos - 3:end - pos])
            # print codon1,seq, pos, end-pos, start-end,'[',start,end,']',alt, strand
            codon2 = str(rev_seq[end - pos - 3]) + str(rev_seq[end - pos - 2]) + var1
    # except:
    #	codon1,codon2='unk','unk'
    # print end-pos, codon_no, rev_seq[end-pos-5:end-pos-1],rev_seq[end-pos-1],rev_seq[end-pos:end-pos+5], codon_no,codon1, all[codon1], acid[codon_no-1], all[codon1] == acid[codon_no-1]
    # print 'NEW AMIN = ' + all[codon2], codon2
    return codon1, codon2


def annotate(gb_record, ref, position, alt):
    # position=position-1
    variant_info = {}
    #	columns = ['snp_id', 'sample_number', 'position', 'ref', 'alt', 'codon', 'codon_number', 'amin_acid', 'effect', 'quality', 'DP', 'intergenic', 'VDB', 'AF1', 'AC1', 'DP4', 'MQ', 'FQ', 'PV4', 'GT', 'PL', 'GQ']
    #	for j in columns:
    #		variant_info[j]='NULL'
    second_method = 0
    flag = 0
    pos = int(position) - 1
    variant_info['position'] = position
    variant_info['ref'] = ref
    variant_info['alt'] = alt
    for ind, feature in enumerate(gb_record.features):
        start = feature.location.nofuzzy_start
        end = feature.location.nofuzzy_end
        strand = feature.location.strand
        keys = feature.qualifiers.keys()
        if start <= pos and pos <= end:
            if feature.type != 'source' and feature.type != 'gene':
                # print gene_info
                flag = 1
                sliced_sequense = gb_record.seq[start:end]
                second_check = sliced_sequense[pos - start:pos - start + len(ref)] == ref
                codon_number = int(round((pos - start) / 3 + 0.5))
                pos_in_gene = pos - start
                for key in keys:
                    if key == 'locus_tag':
                        locus_tag = feature.qualifiers[key][0]
                        variant_info['locus_tag'] = locus_tag
                    elif key == 'gene':
                        gene = feature.qualifiers[key][0].replace("'", "")
                        variant_info['gene'] = gene
                    elif key == 'db_xref':
                        for j in feature.qualifiers[key]:
                            if j.split(':')[0] == 'GI':
                                GI = j.split(':')[1]
                            elif j.split(':')[0] == 'GeneID':
                                GeneID = j.split(':')[1]
                # else:
                #	variant_info['gene']= None
                if second_check:
                    second_method = 1
                    if len(alt) == 1 and len(ref) == 1:
                        codon1, codon2 = codons_def(sliced_sequense, pos, start, end, alt, strand)
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
                    return 'Something goes wrong', ref, position, alt, sliced_sequense[
                                                                       pos - start:pos - start + len(ref)]
        if start >= pos:
            break
    if flag == 0:
        return {'ref': ref, 'alt': alt, 'position': position, 'effect': 'intergenic'}
#
#if __name__ == '__main__':
#
#
#	path_to_genbank="/home/olebedenko/mtb_ref/sequence.gb"
# 	gb_record = SeqIO.read(open(path_to_genbank,"r"), "genbank")
# 	print annotate(gb_record, "C",2155168,"G")

