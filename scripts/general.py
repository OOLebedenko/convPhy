# -*- coding: utf-8 -*-
import subprocess
import os
from scripts.core.phenotype_prediction import get_phenotype_all_nodes
from scripts.core.phyc import phyc

# Подготовка файлов и запуск raxml
def run_raxml(out_dir, raxml_in, phylip_in):
    out_dir_raxml = os.path.join(os.path.abspath(out_dir), "raxml")
    os.makedirs(out_dir_raxml, exist_ok=True)
    os.chdir(out_dir_raxml)
    subprocess.call(["raxmlHPC", "-f", "A", "-t", raxml_in, "-s", phylip_in, "-m", "GTRGAMMA", "-n", "nh"],
                    stdout=os.chdir(out_dir_raxml))
    path_to_anc_phy = out_dir_raxml + "RAxML_marginalAncestralStates.nh"
    return path_to_anc_phy

# Подготовка файлов и запуск phenotype_prediction.py
def run_phenotype_prediction(out_dir, raxml_in, phylip_in, path_to_anc_phy, R_in, S_in):
    out_dir_pheno = os.path.join(out_dir, "phenotype_prediction")
    in_R_phenotype = open(os.path.join(out_dir_pheno, 'positive_phenotype.txt'), "w")
    in_S_phenotype = open(os.path.join(out_dir_pheno, 'negative_phenotype.txt'), "w")
    os.makedirs(out_dir_pheno, exist_ok=True)
    os.chdir(out_dir_pheno)
    name_of_R, name_of_S, names_of_ancestral_S, names_of_ancestral_R, genotype = \
        get_phenotype_all_nodes(raxml_in, phylip_in, path_to_anc_phy, R_in, S_in)
    in_S_phenotype.write('\n'.join(name_of_S + names_of_ancestral_S))
    in_R_phenotype.writelines('\n'.join(name_of_R + names_of_ancestral_R))
    return  name_of_R, name_of_S, names_of_ancestral_S, names_of_ancestral_R, genotype


# Подготовка файлов и запуск phyC
def run_phyc(out_dir,name_of_R, name_of_S, names_of_ancestral_S,
             names_of_ancestral_R, info_pos, raxml_in, genotype):
    out_dir_phyc = os.path.join(out_dir, "phyc")
    os.makedirs(out_dir_phyc, exist_ok=True)
    os.chdir(out_dir_phyc)
    phyc(name_of_R, name_of_S, names_of_ancestral_S,
         names_of_ancestral_R, info_pos, raxml_in, genotype)
#
# out_dir_phyc = out_dir + "phyc/"
#     os.makedirs(out_dir_phyc)
#
# os.chdir(out_dir_phyc)
#
# if not os.path.exists("./input"):
#     os.makedirs("./input")
# if not os.path.exists("./output"):
#     os.makedirs("./output")
# os.chdir("./input")
#
# shutil.copy('../../phenotype_prediction/positive_phenotype.txt', './positive_phenotype.txt')
# shutil.copy('../../phenotype_prediction/negative_phenotype.txt', './negative_phenotype.txt')
# shutil.copy('../../raxml/RAxML_marginalAncestralStates.nh', './ancestral.txt')
#
# shutil.copy(phylip_in, './farhat.phy')
# shutil.copy(raxml_node_labeled_in, 'RAxML_nodeLabelledRootedTree.nh')
# os.rename("RAxML_nodeLabelledRootedTree.nh", "raxml_tree.nh")
#
# os.chdir(out_dir_phyc)
#
# in_R_phenotype = './input/positive_phenotype.txt'
# in_S_phenotype = './input/negative_phenotype.txt'
#
# core.phyc.phyc(pos, SNPs_in, in_R_phenotype, in_S_phenotype, R_in, S_in, phylip_in)
#
# # Подготовка файлов и запуск p_value
#
# logger.info("Preparing of input files and output directory for pvalue run")
# out_dir_p_value = out_dir + "/p_value/"
# if not os.path.exists(out_dir_p_value):
#     os.makedirs(out_dir_p_value)
#
# os.chdir(out_dir_p_value)
# in_pos = out_dir_phyc + "output/pos.txt"
# info_pos = run_dir + "info_pos.txt"
# core.pval.pval(in_pos, info_pos)
#
# # Аннотация ОНП
#
#
# snps_path = out_dir_p_value + "p_value.txt"
# annotation.snps_annotation.snps_ann(snps_path)
#
# # '/home/mrotkevich/Data/Spondilit/filter_final'