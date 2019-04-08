from ete3 import Tree
from tqdm import tqdm


def get_ancestor_phenotype_state_dict(samples_phenotype_file, state):
    phenotype_state_dict = {}
    group_state = []
    with open(samples_phenotype_file) as f_in_pheno:
        for line in f_in_pheno:
            name = line.strip()
            group_state.append(name)
            phenotype_state_dict[name] = state
    return group_state, phenotype_state_dict


def get_genotype_dict(samples_genotype_file, ancestor_phenotype_file):
    genotype_dict = {}
    with open(samples_genotype_file) as f_in_geno, open(ancestor_phenotype_file) as f_anc_geno:
        union_genotype = f_in_geno.readlines() + f_anc_geno.readlines()
        lines_without_spaces = filter(None, (line.strip() for line in union_genotype))
        for line in lines_without_spaces:
            row = line.split()
            id = row[0]
            sequence = row[1]
            genotype_dict[id] = sequence
    return genotype_dict


def count_sequence_distance(sequence_1, sequence_2):
    distance = 0
    for (x, y) in zip(sequence_1, sequence_2):
        if x != y:
            distance += 1
    return distance


def get_ancestor_phenotype(tree, genotype, phenotype):
    ancestor_phenotype = {}
    ancestor_phenotype['R'] = []
    ancestor_phenotype['S'] = []
    pair_nodes = tree.get_tree_root().get_descendants(strategy="levelorder")
    for i in tqdm(range(-1, -len(pair_nodes) + 2, -2), desc="run phenotype prediction"):
        nodei = pair_nodes[i]
        nodej = pair_nodes[i - 1]
        anc = nodei.up.name
        if count_sequence_distance(genotype[nodei.name], genotype[anc]) < \
                count_sequence_distance(genotype[nodej.name], genotype[anc]):
            ancestor_phenotype[phenotype[nodei.name]].append(anc)
            phenotype[anc] = phenotype[nodei.name]
        else:
            ancestor_phenotype[phenotype[nodej.name]].append(anc)
            phenotype[anc] = phenotype[nodej.name]
    return ancestor_phenotype


def get_phenotype_all_nodes(tree_nh, phylip_in, anc_phy, in_R_states, in_S_states):
    tree = Tree(tree_nh, format=1)
    genotype = get_genotype_dict(phylip_in, anc_phy)
    name_of_R, phenotype = get_ancestor_phenotype_state_dict(in_R_states, "R")
    name_of_S, phenotype_S = get_ancestor_phenotype_state_dict(in_S_states, "S")
    phenotype.update(phenotype_S)
    ancestor_phenotype = get_ancestor_phenotype(tree, genotype, phenotype)
    name_of_ancestral_S = name_of_S + ancestor_phenotype['S']
    name_of_ancestral_R = name_of_R + ancestor_phenotype['R']
    return name_of_R, name_of_S, name_of_ancestral_S, name_of_ancestral_R, genotype
