import argparse
from scripts.general import *
from scripts.core.phyc import *
from scripts.annotation.annotate_snp import annotate_snp

parser = argparse.ArgumentParser(description='search SNP ')
parser.add_argument("-i", "--input", required=True, help="path to dir with of all input files")
parser.add_argument("-o", "--out_dir", default="./output_phyc/", help="path to dir_output")
parser.add_argument("-geno", "--genotype_prediction", default=False, type=bool)
parser.add_argument("-pheno", "--phenotype_prediction", default=False, type=bool)
parser.add_argument("-conphy", default=False, type=bool)
parser.add_argument("-p_value", default=False, type=bool)


args = parser.parse_args()
# # входные файлы
run_dir = os.path.join(os.path.abspath(args.input))
info_pos = os.path.join(run_dir, "info_pos.txt")
SNPs_in = os.path.join(run_dir, 'SNPs.txt')
phylip_in = os.path.join(run_dir, 'farhat.phy')
raxml_in = os.path.join(run_dir, "raxml_tree.nh")
R_in = os.path.join(run_dir, 'R_states')
S_in = os.path.join(run_dir, 'S_states')

# выходные файлы
out_dir = os.path.join(os.path.abspath(args.out_dir))
os.makedirs(out_dir, exist_ok=True)

os.chdir(run_dir)

# Подготовка файлов и запуск raxml
if not args.geno:
    path_to_anc_phy = run_raxml(out_dir, raxml_in, phylip_in)
path_to_anc_phy = os.path.join(out_dir, 'raxml', 'RAxML_marginalAncestralStates.nh')

# Подготовка файлов и запуск phenotype_prediction.py
if not args.pheno:
    name_of_R, name_of_S, names_of_ancestral_S, names_of_ancestral_R, genotype = \
        run_phenotype_prediction(out_dir, phylip_in, path_to_anc_phy, R_in, S_in)

# Подготовка файлов и запуск phyC
if not args.phyc:
    R_S = run_phyc(out_dir, name_of_R, name_of_S, names_of_ancestral_S,
         names_of_ancestral_R, info_pos, genotype) #TODO: change on ancestral labeled tree

# Подготовка файлов и запуск p_value
out_dir = os.path.join(args.out_dir, "p_value")
if not args.p_vaue:
    with open(os.path.join(args.out_dir, "phyc", "pos.txt")) as f_in:
        R_S = [line.strip().split() for line in f_in]
    run_p_value(out_dir, R_S, info_pos)

# Аннотация ОНП

snps_path = os.path.join(args.out_dir, "p_value", "p_value.txt")
path_to_snps_out_csv = os.path.join(args.out_dir, "p_value", "annotated_snps.csv")
annotate_snp(snps_path, path_to_snps_out_csv)
