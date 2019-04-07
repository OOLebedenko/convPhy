import argparse
from scripts.general import *
from scripts.core.phyc import *

parser = argparse.ArgumentParser(description='search SNP ')
parser.add_argument("-i", "--input", required=True, help="path to dir with of all input files")
parser.add_argument("-o", "--out_dir", default="./output_phyc/", help="path to dir_output")
# parser.add_argument("-r","--reference", default="/home/olebedenko/mnt_olebedenko/ref.fasta", help="path to reference")

args = parser.parse_args()
# входные файлы
run_dir = os.path.join(os.path.abspath(args.input))
# выходные файлы
out_dir = os.path.join(os.path.abspath(args.out_dir))
os.makedirs(out_dir, exist_ok=True)
os.chdir(run_dir)
info_pos = os.path.join(run_dir, "info_pos.txt")
SNPs_in = os.path.join(run_dir, 'SNPs.txt')
phylip_in = os.path.join(run_dir, 'farhat.phy')
raxml_in = os.path.join(run_dir, "raxml_tree.nh")
R_in = os.path.join(run_dir, 'R_states')
S_in = os.path.join(run_dir, 'S_states')

# тестируемые позиции в геноме
# pos = map(lambda line: int(line.strip().split("\t")[0]), open(info_pos))

# Подготовка файлов и запуск raxml
# path_to_anc_phy = run_raxml(out_dir, raxml_in, phylip_in)
path_to_anc_phy = os.path.join(out_dir, 'raxml', 'RAxML_marginalAncestralStates.nh')

# Подготовка файлов и запуск phenotype_prediction.py
name_of_R, name_of_S, names_of_ancestral_S, names_of_ancestral_R, genotype = \
    run_phenotype_prediction(out_dir, raxml_in, phylip_in, path_to_anc_phy, R_in, S_in)

# Подготовка файлов и запуск phyC
run_phyc(out_dir, name_of_R, name_of_S, names_of_ancestral_S,
         names_of_ancestral_R, info_pos, raxml_in, genotype)
