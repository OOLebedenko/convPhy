from ete3 import Tree


def get_ancestor_phenotype_state_dict(samples_phenotype_file, state):
    phenotype_state_dict = {}
    group_state = []
    with open(samples_phenotype_file) as f_in_pheno:
        for line in f_in_pheno:
            name = line.strip()
            group_state.append(name)
            phenotype_state_dict[name] = state
    return group_state, phenotype_state_dict


def get_genotype_list(samples_genotype_file, ancestor_phenotype_file):
    genotype_list = {}
    with open(samples_genotype_file) as f_in_geno, open(ancestor_phenotype_file) as f_anc_geno:
        union_genotype = f_in_geno.readlines() + f_anc_geno.readlines()
        lines_without_spaces = filter(None, (line.strip() for line in union_genotype))
        for line in lines_without_spaces:
            row = line.split()
            id = row[0]
            sequence = row[1]
            genotype_list[id] = sequence
    return genotype_list


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
    for i in range(-1, -len(pair_nodes) + 2, -2):
        nodei = pair_nodes[i]
        nodej = pair_nodes[i - 1]
        anc = nodei.up.name
        print(anc)
        if count_sequence_distance(genotype[nodei.name], genotype[anc]) < \
                count_sequence_distance(genotype[nodej.name], genotype[anc]):
            ancestor_phenotype[phenotype[nodei.name]].append(anc)
            phenotype[anc] = phenotype[nodei.name]
        else:
            ancestor_phenotype[phenotype[nodej.name]].append(anc)
            phenotype[anc] = phenotype[nodej.name]
    return ancestor_phenotype


def get_phenotype_all_nodes(tree_nh, phylo_phy, anc_phy, in_R_states, in_S_states):
    in_R_phenotype = open('positive_phenotype.txt', "w")
    in_S_phenotype = open('negative_phenotype.txt', "w")
    tree = Tree(tree_nh, format=1)
    genotype = get_genotype_list(phylo_phy, anc_phy)
    name_of_R, phenotype = get_ancestor_phenotype_state_dict(in_R_states, "R")
    name_of_S, phenotype_S = get_ancestor_phenotype_state_dict(in_S_states, "S")
    phenotype.update(phenotype_S)
    ancestor_phenotype = get_ancestor_phenotype(tree, genotype, phenotype)
    negative_phenotype_all_nodes = name_of_S + ancestor_phenotype['S']
    positive_phenotype_all_nodes = name_of_R + ancestor_phenotype['R']
    in_S_phenotype.write('\n'.join(negative_phenotype_all_nodes))
    in_R_phenotype.writelines('\n'.join(positive_phenotype_all_nodes))
    return negative_phenotype_all_nodes, positive_phenotype_all_nodes


if __name__ == '__main__':
    in_R_states = './phyc_test/input/R_states'
    in_S_states = './phyc_test/input/S_states'
    phylo_phy = './phyc_test/input/farhat.phy'
    tree_nh = './phyc_test/input/raxml_tree.nh'
    anc_phy = './phyc_test/input/raxml/RAxML_marginalAncestralStates.nh'
    get_phenotype_all_nodes(tree_nh, phylo_phy, anc_phy, in_R_states, in_S_states)



