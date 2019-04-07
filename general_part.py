import argparse, subprocess, os, shutil, unicodedata
import annotation
import core

parser = argparse.ArgumentParser(description='search SNP ')
parser.add_argument("input", help="path to dir with of all input files")
parser.add_argument("-vcf", default=False, help="path to dir with vcf files")
parser.add_argument("-o", "--out_dir", default="./output_phyc/", help="path to dir_output")
# parser.add_argument("-r","--reference", default="/home/olebedenko/mnt_olebedenko/ref.fasta", help="path to reference")

args = parser.parse_args()


# входные файлы
run_dir = os.path.abspath(args.input) + "/"

os.chdir(run_dir)

if args.vcf:
    core.create_phylip.create_phy_pos_SNPs_files(args.vcf)
else:
    info_pos = run_dir + "info_pos.txt"
    SNPs_in = run_dir + 'SNPs.txt'

    phylip_in = run_dir + "farhat.phy"

raxml_in = run_dir + "raxml_tree.nh"
S_in = run_dir + "S_states"
R_in = run_dir + "R_states"

# тестируемые позиции в геноме
pos = map(lambda line: int(line.strip().split("\t")[0]), open(info_pos))

# выходные файлы
out_dir = os.path.abspath(args.out_dir) + "/"


# Подготовка файлов и запуск raxml
out_dir_raxml = os.path.abspath(args.out_dir) + "/raxml/"
if not os.path.exists(out_dir_raxml):
    os.makedirs(out_dir_raxml)

os.chdir(out_dir_raxml)

subprocess.call(["raxmlHPC", "-f", "A", "-t", raxml_in, "-s", phylip_in, "-m", "GTRGAMMA", "-n", "nh"],
                stdout=os.chdir(out_dir_raxml))

raxml_node_labeled_in = out_dir_raxml + "RAxML_nodeLabelledRootedTree.nh"

# Подготовка файлов и запуск phenotype_prediction.py


out_dir_pheno = out_dir + "/phenotype_prediction/"
if not os.path.exists(out_dir_pheno):
    os.makedirs(out_dir_pheno)

os.chdir(out_dir_pheno)
anc_phy = open(out_dir_raxml + "RAxML_marginalAncestralStates.nh")
leaves_phy = open(phylip_in)

all_nodes_phy = (anc_phy.read() + leaves_phy.read()).split("\n")

core.pheno.pheno(raxml_node_labeled_in, fasta_list, R_in, S_in)

# Подготовка файлов и запуск phyC


out_dir_phyc = out_dir + "phyc/"
if not os.path.exists(out_dir_phyc):
    os.makedirs(out_dir_phyc)

os.chdir(out_dir_phyc)

if not os.path.exists("./input"):
    os.makedirs("./input")
if not os.path.exists("./output"):
    os.makedirs("./output")
os.chdir("./input")

shutil.copy('../../phenotype_prediction/positive_phenotype.txt', './positive_phenotype.txt')
shutil.copy('../../phenotype_prediction/negative_phenotype.txt', './negative_phenotype.txt')
shutil.copy('../../raxml/RAxML_marginalAncestralStates.nh', './ancestral.txt')

shutil.copy(phylip_in, './farhat.phy')
shutil.copy(raxml_node_labeled_in, 'RAxML_nodeLabelledRootedTree.nh')
os.rename("RAxML_nodeLabelledRootedTree.nh", "raxml_tree.nh")

os.chdir(out_dir_phyc)

in_R_phenotype = './input/positive_phenotype.txt'
in_S_phenotype = './input/negative_phenotype.txt'

core.phyc.phyc(pos, SNPs_in, in_R_phenotype, in_S_phenotype, R_in, S_in, phylip_in)

# Подготовка файлов и запуск p_value

logger.info("Preparing of input files and output directory for pvalue run")
out_dir_p_value = out_dir + "/p_value/"
if not os.path.exists(out_dir_p_value):
    os.makedirs(out_dir_p_value)

os.chdir(out_dir_p_value)
in_pos = out_dir_phyc + "output/pos.txt"
info_pos = run_dir + "info_pos.txt"
core.pval.pval(in_pos, info_pos)

# Аннотация ОНП


snps_path = out_dir_p_value + "p_value.txt"
annotation.snps_annotation.snps_ann(snps_path)

# '/home/mrotkevich/Data/Spondilit/filter_final'