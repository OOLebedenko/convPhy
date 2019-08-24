import os
from collections import defaultdict
import operator
import argparse
from tqdm import tqdm
from glob import glob


def add_to_pull_for_one_vcf_file(path_to_vcf):
    ref_pull = {}
    alt_pull = {}
    with open(path_to_vcf) as f:
        for line in f:
            if line[0] != '#':
                row = line.strip().split()
                pos = row[1]
                ref = row[3]
                alt = row[4]
                if len(alt) == 1 and len(ref) == 1:
                    ref_pull[int(pos)] = ref
                    alt_pull[int(pos)] = alt
    return ref_pull, alt_pull


def add_to_pull_for_multiple_vcf_file(path_to_directory_with_vcf_files, path_to_outgroup=None):
    vcf_files = glob(os.path.join(path_to_directory_with_vcf_files, "*.vcf"))
    ref_pull_all = {}
    alt_pull_all = {}
    for path_to_vcf in vcf_files:
        temp_ref, temp_alt = add_to_pull_for_one_vcf_file(path_to_vcf)
        ref_pull_all.update(temp_ref)
        alt_pull_all.update(temp_alt)
    if path_to_outgroup:
        temp_ref, temp_alt = add_to_pull_for_one_vcf_file(path_to_outgroup)
        ref_pull_all.update(temp_ref)
        alt_pull_all.update(temp_alt)
    return ref_pull_all, alt_pull_all


def create_sample_dict(path_to_vcf):
    sample_dict = {}
    positions = []
    alts = []
    with open(path_to_vcf) as f:
        for line in f:
            if line[0] != '#':
                row = line.strip().split()
                pos = row[1]
                ref = row[3]
                alt = row[4]
                if len(alt) == 1 and len(ref) == 1:
                    positions.append(int(pos))
                    alts.append(alt)
        name = os.path.basename(path_to_vcf).split(".")[0]
        sample_dict[name] = positions, alts
    return sample_dict


def write_phylip(path_to_directory_with_vcf_files, path_to_out_dir, f_out_name, path_to_outgroup=None):
    ref_pull_all, alt_pull_all = add_to_pull_for_multiple_vcf_file(path_to_directory_with_vcf_files, path_to_outgroup)
    vcf_files = glob(os.path.join(path_to_directory_with_vcf_files, "*.vcf"))
    if path_to_outgroup:
        vcf_files.append(path_to_outgroup)
    sorted_position = sorted(alt_pull_all.keys())
    with open(os.path.join(path_to_out_dir, f_out_name), "w") as phylip:
        phylip.write('{number_of_samples} {number_of_snps}\n'.format(number_of_samples=len(vcf_files),
                                                                     number_of_snps=len(alt_pull_all)))
        for path_to_vcf in tqdm(vcf_files):
            sample_dict = create_sample_dict(path_to_vcf)
            name = os.path.basename(path_to_vcf).split(".")[0]
            phylip.write(name + " ")
            for position in sorted_position:
                try:
                    k = sample_dict[file][0].index(position)
                    phylip.write(sample_dict[file][1][k])
                except:
                    phylip.write(ref_pull_all[position])
            phylip.write("\n")


def create_snps_file(path_to_directory_with_vcf_files, path_to_out_directory):
    ref_pull_all, alt_pull_all = add_to_pull_for_multiple_vcf_file(path_to_directory_with_vcf_files)
    snps = [alt_pull_all[i] for i in sorted(alt_pull_all)]
    path_to_out_snps = os.path.join(path_to_out_directory, "SNPs.txt")
    with open(path_to_out_snps, "w") as out_snps:
        out_snps.write("".join(snps))


def create_info_pos(path_to_directory_with_vcf_files, path_to_out_directory):
    ref_pull_all, alt_pull_all = add_to_pull_for_multiple_vcf_file(path_to_directory_with_vcf_files)
    sorted_position = [str(key) for key in sorted(ref_pull_all.keys())]
    snps = [alt_pull_all[i] for i in sorted(alt_pull_all)]
    ref = [ref_pull_all[i] + "\n" for i in sorted(ref_pull_all)]
    path_to_out_info_pos = os.path.join(path_to_out_directory, "info_pos.txt")
    with open(path_to_out_info_pos, "w") as info_pos_out:
        info_pos_out.writelines(map(lambda line: "\t".join(line), zip(sorted_position, snps, ref)))



if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='make input files for convPhy')
	parser.add_argument("-vcf", "--path-to-directory-with-vcf-files", required=True, help="path to dir with vcf files")
	parser.add_argument("-outgroup", "--path-to-outgroup", required=True, help="path to outgroup")
	parser.add_argument("-out-phy", "--name-phylip-file", default="convphy.phy", help="path to outg phylip")
	parser.add_argument("-o", "--path_to_out_directory", default="input", help="path to out directory")

	args = parser.parse_args()

	os.makedirs(args.path_to_out_directory, exist_ok=True)

	write_phylip(args.path_to_directory_with_vcf_files, args.path_to_out_directory, args.name_phylip_file, args.path_to_outgroup)
	create_snps_file(args.path_to_directory_with_vcf_files, args.path_to_out_directory)
	create_info_pos(args.path_to_directory_with_vcf_files, args.path_to_out_directory)
