# convPhy 

## Project description

In this repository we present tool for discovery convergent mutations associated 
with drug resistant phenotype of *M. tberculosis*. 
convPhy consider previously reconstruction phylogenetic tree. 
For each testing SNPs convPhy predict genotype and phenotype all ancestral nodes.
And than significance SNPs check by permutation test:

![](https://drive.google.com/uc?id=13venwzghWsDKznleVb5IaJsLxJnvvwEZ)  


## Dependencies:

          ![RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)
          python3 pip install tqdm
          python3 pip install ete3
          python3 pip install six
          python3 pip install numpy
          

## Run convPhy

### Preparing input files

    The next files with **fixed names** are necessary to run convPhy:
        (if your files have another name you should rename it to run convPhy)
        - info_pos.txt 
            tab separated files **without header** with three columns: 
                snp position, snp, reference
                like this:
                    1	G	T
                    4	G	A
                    682	A	T
        - convPhy.phy
            concatenated SNPs in phylib format 
            (format for many phylogenetic tree reconstruction program)
        - raxml_tree.nh
            previously reconstructed phylogenetic tree with concatenated SNPs as input
        - R_states
            list of all positive (drug resistant) phenotype samples
            like this:
                TB0001
                TB0003
                TB0005
        - S_states
            list of all negative (drug sensitive) phenotype samples
            like this:
                TB0002
                TB0004
                
    All input files placed in one directory.

You you may prepare all input files by own or use **create_input.py** script for creation info_pos.txt and convPhy.phy:
    python3 create_input.py -vcf path_to_directory_with_vcf_files -outgroup outgroup_vcf_files -o output_directory
        path_to_directory_with_vcf_files - directory with vcf files in your project
        outgroup vcf - vcf file for outgroup which you use for phylogenetic tree reconstruction
        output directory - directory with input files vy default named "input" 
        
### Get started
    python3 convphy.py  -i path_to_input_directory -o path_to_output_directory -path_to_genbank path_to_genbankfile
    

## Output files are placed into for directory
    raxml - files with result genotype prediction by RAxML
    phenotype_prediction:
        positive_phenotype.txt - all positive (drug resistant) phenotype nodes on
                                 phylogenetic tree including ancestral nodes
        negative_phenotype.txt - all negative (drug sensitive) phenotype nodes on
                                 phylogenetic tree including ancestral nodes
    phyc:
        pos.txt - file with number of positive (resistant) and negative (sensitive)
                  branches for all nucleotide positions
    p_value:
        p_value.txt - significant p_value positions (threshold 0.05)
        annotated_snps.csv - annotated significant snps 
