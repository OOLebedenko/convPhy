import pandas as pd
import Bio.Data.CodonTable
from Bio import SeqIO
from tqdm import tqdm


def annotate_snp(snps_path, path_to_snps_out_csv, path_to_genbank):
    gb_record = SeqIO.read(open(path_to_genbank, "r"), "genbank")
    snps = pd.read_csv(snps_path)
    positions = snps['pos']
    refs = snps['ref']    
    alts = snps['alt']
    annotated_snps = pd.DataFrame()

    for ref, pos, alt in tqdm(zip(refs, positions, alts)):
        variant_info = annotate(gb_record, ref, pos, alt)
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

        elif (pos - start) % 3 == 2:
            codon1 = str(seq[pos - start - 2:pos - start + 1])
            codon2 = codon1[0] + codon1[1] + alt
        codon_number = int(round((pos - start) / 3 + 0.5))
        pos_in_gene = pos - start

    elif strand == -1:
        var1 = alt_complement(alt)
        rev_seq = seq.reverse_complement()
        if (end - pos - 1) % 3 == 0:
            codon1 = str(rev_seq[end - pos - 1:end - pos + 2])
            codon2 = var1 + str(rev_seq[end - pos]) + str(rev_seq[end - pos + 1])
        elif (end - pos - 1) % 3 == 1:
            codon1 = str(rev_seq[end - pos - 2:end - pos + 1])
            codon2 = str(rev_seq[end - pos - 2]) + var1 + str(rev_seq[end - pos])
        elif (end - pos - 1) % 3 == 2:
            codon1 = str(rev_seq[end - pos - 3:end - pos])
            codon2 = str(rev_seq[end - pos - 3]) + str(rev_seq[end - pos - 2]) + var1
        codon_number = int(round((end - pos - 1) / 3 + 0.5))
        pos_in_gene = end - pos - 1
    return codon1, codon2, codon_number, pos_in_gene

def get_cds_variant_info(sequence, pos, start, end, alt, strand):
    codon1, codon2, codon_number, pos_in_gene = codons_def(sequence, pos, start, end, alt, strand)
    eff, acid = effect(codon1, codon2)
    codon = str(codon1) + "/" + str(codon2)
    codon_info = {'codon_number' : int(codon_number), \
                         'pos_in_gene': pos_in_gene, 'codon': codon, 'effect': eff, \
                         'amino_acid': acid, 'alt': alt}
    return codon_info


def annotate(gb_record, ref, position, alt):
    pos = int(position) - 1
    variant_info = { 'position': position, 'ref': ref, 'alt': alt}
    for feature in gb_record.features:
        start = feature.location.nofuzzy_start
        end = feature.location.nofuzzy_end
        if pos > end:
            continue
        feature_gb_record_sequence = gb_record.seq[start:end]
        ref_gb = feature_gb_record_sequence[pos - start:pos - start + len(ref)]
        strand = feature.location.strand
        locus_tag = feature.qualifiers.get('locus_tag')
        locus_tag = locus_tag[0] if locus_tag else None
        gene = feature.qualifiers.get('gene') 
        gene = gene[0] if gene else None
        feature_type = feature.type
        if start <= pos:
            if feature_type  == 'source':
                continue
            elif feature_type  == 'gene':
                continue
            elif ref_gb != ref:
                replaced_ref = "It's not SNP. Nucleotide position {pos} in the genome sequence has been corrected {ref}:{new_ref}".\
                               format(pos=pos, ref=ref, new_ref=ref_gb)
                variant_info.update({"type" : replaced_ref})
                return variant_info
            elif feature_type  == 'CDS':
                codon_info = get_cds_variant_info(feature_gb_record_sequence, pos, start, end, alt, strand)
                variant_info.update({'type':feature_type, 'locus_tag': locus_tag, 'gene': gene})
                variant_info.update(codon_info)
                return variant_info
            else: 
                variant_info.update({"type" : "{feature_type}: RNA or exactly unknown about translation this sequence type".\
                                    format(feature_type=feature_type ), 'locus_tag': locus_tag, 'gene': gene})                
                return variant_info
        else:
            variant_info.update({'effect': 'intergenic'})
            return variant_info
