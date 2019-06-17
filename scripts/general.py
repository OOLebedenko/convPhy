# -*- coding: utf-8 -*-
import subprocess
import os
from scripts.core.phenotype_prediction import get_phenotype_all_nodes
from scripts.core.phyc import phyc
from scripts.core.p_value import get_p_value


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
def run_phenotype_prediction(out_dir, phylip_in, path_to_anc_phy, R_in, S_in):
    out_dir_pheno = os.path.join(out_dir, "phenotype_prediction")
    raxml_in = os.path.join(out_dir, "raxml", "RAxML_nodeLabelledRootedTree.nh")
    os.makedirs(out_dir_pheno, exist_ok=True)
    in_R_phenotype = open(os.path.join(out_dir_pheno, 'positive_phenotype.txt'), "w")
    in_S_phenotype = open(os.path.join(out_dir_pheno, 'negative_phenotype.txt'), "w")
    os.chdir(out_dir_pheno)
    name_of_R, name_of_S, names_of_ancestral_S, names_of_ancestral_R, genotype = \
        get_phenotype_all_nodes(raxml_in, phylip_in, path_to_anc_phy, R_in, S_in)
    in_S_phenotype.write('\n'.join(name_of_S + names_of_ancestral_S))
    in_R_phenotype.writelines('\n'.join(name_of_R + names_of_ancestral_R))
    return  name_of_R, name_of_S, names_of_ancestral_S, names_of_ancestral_R, genotype


# Подготовка файлов и запуск phyC
def run_phyc(out_dir,name_of_R, name_of_S, names_of_ancestral_S,
             names_of_ancestral_R, info_pos, genotype):
    raxml_in = os.path.join(out_dir, "raxml", "RAxML_nodeLabelledRootedTree.nh")
    out_dir_phyc = os.path.join(out_dir, "phyc")
    os.makedirs(out_dir_phyc, exist_ok=True)
    os.chdir(out_dir_phyc)
    R_S = phyc(name_of_R, name_of_S, names_of_ancestral_S,
         names_of_ancestral_R, info_pos, raxml_in, genotype)
    return R_S

#Подготовка файлов и запуск p_value
def run_p_value(out_dir, R_S, info_pos):
    out_dir_p_value = os.path.join(out_dir,"p_value")
    os.makedirs(out_dir_p_value, exist_ok=True)
    os.chdir(out_dir_p_value)
    p_value = './p_value.txt'
    p_val = get_p_value(R_S, info_pos)
    with open(p_value, "w") as f_out:
        f_out.writelines(p_val)


# # Аннотация ОНП
#
#
# snps_path = out_dir_p_value + "p_value.txt"
# annotation.snps_annotation.snps_ann(snps_path)
#