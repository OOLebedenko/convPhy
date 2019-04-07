import os
from ete3 import Tree
from collections import defaultdict
from scripts.core.phenotype_prediction import *
from scripts.general import *


def _prepare_pos(info_pos):
    pos = [int(line.strip().split("\t")[0]) for line in open(info_pos)]
    SNPs = [line.strip().split("\t")[1] for line in open(info_pos)]
    return pos, SNPs

def _prepare_structure(names_of_nodes, genotype, SNPs, ancestral=False):
    INFO_POS = defaultdict(lambda: defaultdict(list))
    for name in names_of_nodes:
        sequence = genotype[name]
        for ind in range(len(sequence)):
            if sequence[ind] == SNPs[ind]:
                if ancestral:
                    INFO_POS[ind][name].append(sequence[ind])
                else:
                    INFO_POS[ind][sequence[ind]].append(name)
    return INFO_POS

def traverse(leaves_INFO_POS,ancestral_INFO_POS, names_of_ancestral_R, tree, point, snp):
    R_list = []
    S_list = []
    for lleaf in leaves_INFO_POS[point][snp]:
        node = tree & lleaf
        ans_node = node.up
        while ancestral_INFO_POS[point][ans_node.name] == [snp]:
            node = ans_node
            ans_node = node.up
            if node.is_root():
                break
        if node.is_leaf():
            continue

        elif node.name in names_of_ancestral_R:
            R_list.append(node.name)
        else:
            S_list.append(node.name)
    l1 = R_list
    l2 = S_list
    return l1, l2


def phyc(name_of_R, name_of_S, names_of_ancestral_S, names_of_ancestral_R, info_pos, raxml_in):
    tree = Tree(raxml_in, format=1)
    names_of_ancestral = names_of_ancestral_R + names_of_ancestral_S
    resistant_branch = open('./pos.txt', 'w') #TODO: remove hardcode
    pos, SNPs = _prepare_pos(info_pos)
    R_INFO_POS = _prepare_structure(name_of_R, genotype, SNPs)
    S_INFO_POS = _prepare_structure(name_of_S, genotype, SNPs)
    ancestral_INFO_POS = _prepare_structure(names_of_ancestral, genotype, SNPs, True)

    for point in range(len(pos)):
        print(point)
        R1, S1 = traverse(R_INFO_POS, ancestral_INFO_POS, names_of_ancestral_R, tree, point, SNPs[point])
        R2, S2 = traverse(S_INFO_POS, ancestral_INFO_POS, names_of_ancestral_R, tree, point, SNPs[point])
        l1 = len(set(R1 + R2))
        l2 = len(set(S1 + S2))
        resistant_branch.write(str(point) + "\t" + str(l1) + "\t" + str(l2) + "\n")


if __name__ == "__main__":
    in_S_phenotype = './phyc_test/input/negative_phenotype.txt'
    in_R_phenotype = './phyc_test/input/positive_phenotype.txt'

    R_in = './phyc_test/input/R_states'
    S_in = './phyc_test/input/S_states'


    info_pos = './phyc_test/input/info_pos.txt'

    phylip_in = './phyc_test/input/farhat.phy'
    nodes_seq = './phyc_test/input/ancestral.txt'
    raxml_in = './phyc_test/input/raxml_tree.nh'
    path_to_anc_phy = './phyc_test/input/raxml/RAxML_marginalAncestralStates.nh'



    genotype = get_genotype_dict(phylip_in, path_to_anc_phy)
    name_of_R, name_of_S, names_of_ancestral_S, names_of_ancestral_R = \
        run_phenotype_prediction(os.getcwd(), raxml_in, phylip_in, path_to_anc_phy, R_in, S_in)
    pos, SNPs = _prepare_pos(info_pos)

    phyc(name_of_R, name_of_S, names_of_ancestral_S, names_of_ancestral_R, info_pos, raxml_in)

