from scripts.creation_input_files.create import *


parser = argparse.ArgumentParser(description='make input files for convPhy')
parser.add_argument("-vcf", "--path-to-directory-with-vcf-files", required=True, help="path to dir with vcf files")
parser.add_argument("-outgroup", "--path-to-outgroup", required=True, help="path to outgroup")
parser.add_argument("-out-phy", "--name-phylip-file", default="convphy.phy", help="path to outg phylip")
parser.add_argument("-o", "--path_to_out_directory", default="input", help="path to out directory")

args = parser.parse_args()

os.makedirs(args.path_to_out_directory, exist_ok=True)

write_phylip(args.path_to_directory_with_vcf_files, args.path_to_out_directory, args.name_phylip_file, args.path_to_outgroup)
create_snps_file(args.path_to_directory_with_vcf_files, args.path_to_out_directory, args.path_to_outgroup)
create_info_pos(args.path_to_directory_with_vcf_files, args.path_to_out_directory, args.path_to_outgroup)
