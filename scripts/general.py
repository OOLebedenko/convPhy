import subprocess
import os
import pandas as pd
from scripts.core.phenotype_prediction import get_phenotype_all_nodes
from scripts.core.phyc import phyc
from scripts.core.p_value import get_p_value


def read_file_by_line(path_to_file, split_by_any_space_separater=False):
    with open(path_to_file) as f_in:
        if split_by_any_space_separater:
            out_list = [line.strip().split() for line in f_in]
        else:
            out_list = [line.strip() for line in f_in]
    return out_list

def run_raxml(out_dir, raxml_in, phylip_in):
    out_dir_raxml = os.path.join(os.path.abspath(out_dir), "raxml")
    os.makedirs(out_dir_raxml, exist_ok=True)
    os.chdir(out_dir_raxml)
    subprocess.call(["raxmlHPC", "-f", "A", "-t", raxml_in, "-s", phylip_in, "-m", "GTRGAMMA", "-n", "nh"],
                    stdout=os.chdir(out_dir_raxml))


def run_phenotype_prediction(out_dir, phylip_in, path_to_anc_phy, R_in, S_in):
    out_dir_pheno = os.path.join(out_dir, "phenotype_prediction")
    raxml_in = os.path.join(out_dir, "raxml", "RAxML_nodeLabelledRootedTree.nh")
    os.makedirs(out_dir_pheno, exist_ok=True)
    ancestor_R_phenotype = open(os.path.join(out_dir_pheno, 'positive_phenotype.txt'),
                                "w")  # TODO  make with context manager
    ancestor_S_phenotype = open(os.path.join(out_dir_pheno, 'negative_phenotype.txt'), "w")
    os.chdir(out_dir_pheno)
    name_of_R, name_of_S, names_of_ancestral_S, names_of_ancestral_R = \
        get_phenotype_all_nodes(raxml_in, phylip_in, path_to_anc_phy, R_in, S_in)
    ancestor_S_phenotype.write('\n'.join(names_of_ancestral_S))
    ancestor_R_phenotype.writelines('\n'.join(names_of_ancestral_R))


def run_phyc(out_dir, name_of_R, name_of_S, names_of_ancestral_S,
             names_of_ancestral_R, info_pos, genotype_dict):
    raxml_in = os.path.join(out_dir, "raxml", "RAxML_nodeLabelledRootedTree.nh")
    os.makedirs(os.path.join(out_dir, "phyc"), exist_ok=True)
    R_S = phyc(name_of_R, name_of_S, names_of_ancestral_S,
               names_of_ancestral_R, info_pos, raxml_in, genotype_dict)
    R_S.to_csv( os.path.join(out_dir, "phyc", 'pos.csv'), index=False)

def run_p_value(out_dir, path_to_R_S):

    R_S = pd.read_csv(os.path.join(out_dir, "phyc", 'pos.csv'))
    out_dir_p_value = os.path.join(out_dir, "p_value")
    os.makedirs(out_dir_p_value, exist_ok=True)
    p_value_out_csv = os.path.join(out_dir_p_value, 'p_value.csv')
    p_value = get_p_value(R_S)
    p_value.to_csv(p_value_out_csv, index=False)
