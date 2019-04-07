import argparse
from scripts.general import *

parser = argparse.ArgumentParser(description='search SNP ')
parser.add_argument("-i", "--input", required=True, help="path to dir with of all input files")
parser.add_argument("-o", "--out_dir", default="./output_phyc/", help="path to dir_output")
# parser.add_argument("-r","--reference", default="/home/olebedenko/mnt_olebedenko/ref.fasta", help="path to reference")

args = parser.parse_args()
# входные файлы
run_dir = os.path.join(os.path.abspath(args.input))
os.chdir(run_dir)

info_pos = run_dir + "info_pos.txt"
SNPs_in = run_dir + 'SNPs.txt'
#phylip_in = run_dir + "farhat.phy"
phylip_in = './phyc_test/input/farhat.phy'
#raxml_in = run_dir + "raxml_tree.nh"
raxml_in = './phyc_test/input/raxml_tree.nh'
# S_in = run_dir + "S_states"
# R_in = run_dir + "R_states"
R_in = './phyc_test/input/R_states'
S_in = './phyc_test/input/S_states'

# тестируемые позиции в геноме
#pos = map(lambda line: int(line.strip().split("\t")[0]), open(info_pos))

# выходные файлы
out_dir =os.path.join(os.path.abspath(args.out_dir))
os.makedirs(out_dir, exist_ok=True)

# Подготовка файлов и запуск raxml
#path_to_anc_phy = run_raxml(out_dir, raxml_in, phylip_in)

# Подготовка файлов и запуск phenotype_prediction.py
path_to_anc_phy = './phyc_test/input/raxml/RAxML_marginalAncestralStates.nh'
run_phenotype_prediction(os.getcwd(), raxml_in, phylip_in, path_to_anc_phy, R_in, S_in)

