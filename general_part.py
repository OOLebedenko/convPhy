#!/usr/bin/env python3
import argparse
import os
from scripts.general import *
from scripts.core.phenotype_prediction import get_genotype_dict
from scripts.annotation.annotate_snp import annotate_snp

parser = argparse.ArgumentParser(description='search SNP ')
parser.add_argument("-i", "--input", required=True, help="path to dir with all input files")
parser.add_argument("-o", "--out_dir", default="./output_phyc/", help="path to dir_output")
parser.add_argument("-geno", "--genotype_prediction", default=False, type=bool)
parser.add_argument("-pheno", "--phenotype_prediction", default=False, type=bool)
parser.add_argument("-convphy", default=False, type=bool)
parser.add_argument("-p_value", default=False, type=bool)
parser.add_argument("-neg", "--negative_phenotype", type=str)
parser.add_argument("-pos", "--positive_phenotype", type=str)
#parser.add_argument("-ref", default=False, type=str, help="path to reference genbank file to make an annotation of variants" )

args = parser.parse_args()
#print (args)
# input files
run_dir = os.path.join(os.path.abspath(args.input))
info_pos = os.path.join(run_dir, "info_pos.txt")
SNPs_in = os.path.join(run_dir, 'SNPs.txt')
phylip_in = os.path.join(run_dir, 'farhat.phy')
raxml_in = os.path.join(run_dir, "raxml_tree.nh")

if args.negative_phenotype:
	S_in = os.path.join(run_dir, args.negative_phenotype)
else:
	S_in = os.path.join(run_dir, "S_states.txt")
if args.positive_phenotype:
	R_in = os.path.join(run_dir, args.positive_phenotype)
else:
	R_in = os.path.join(run_dir, "R_states.txt")

# check input files existance
for filepath in [info_pos,SNPs_in,phylip_in, raxml_in, R_in,S_in]: 
	if not os.path.isfile(filepath):
		print ("Error! File '%s', which is necessary for executing, doesn't exist!" % filepath)
		raise SystemExit(0)


# output files
out_dir = os.path.join(os.path.abspath(args.out_dir))
os.makedirs(out_dir, exist_ok=True)
os.chdir(run_dir)

# run genotype prediction
if not args.genotype_prediction:
    run_raxml(out_dir, raxml_in, phylip_in)
path_to_ancestor_phylip = os.path.join(out_dir, 'raxml', 'RAxML_marginalAncestralStates.nh')

# run phenotype prediction
if not args.phenotype_prediction:
    run_phenotype_prediction(out_dir, phylip_in, path_to_ancestor_phylip, R_in, S_in)

path_to_ancestor_S_phenotype = os.path.join(out_dir, "phenotype_prediction", 'negative_phenotype.txt')
path_to_ancestor_R_phenotype = os.path.join(out_dir, "phenotype_prediction", 'positive_phenotype.txt')
name_of_S = read_file_by_line(S_in)
name_of_R = read_file_by_line(R_in)
ancestor_S_phenotype = read_file_by_line(path_to_ancestor_S_phenotype)
ancestor_R_phenotype = read_file_by_line(path_to_ancestor_R_phenotype)
genotype_dict = get_genotype_dict(phylip_in, path_to_ancestor_phylip)

# run convPhy
if not args.convphy:
    run_phyc(out_dir, name_of_R, name_of_S, ancestor_S_phenotype,
             ancestor_R_phenotype, info_pos, genotype_dict)

path_to_R_S = os.path.join(out_dir, "phyc", 'pos.txt')
R_S = read_file_by_line(path_to_R_S, split_by_any_space_separater=True)

# run permutation test
if not args.p_value:
    run_p_value(out_dir, R_S, info_pos)

# annotation SNPs
snps_path = os.path.join(out_dir, "p_value", "p_value.txt")
path_to_snps_out_csv = os.path.join(out_dir, "p_value", "annotated_snps.csv")
#annotate_snp(snps_path, path_to_snps_out_csv, args.ref)
