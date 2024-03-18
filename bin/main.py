import os
import argparse
import logging

import utils
import parse_vcf_data as pvd
import vcf_analysis as va
import mutational_signature as mut
import vcf_analysis as va


def get_args():
    parser = argparse.ArgumentParser(prog="vcf_analysis")
    parser.add_argument("vcf", help = "Input vcf file")
    parser.add_argument("-n", dest="name", default=None, help="Output file name without the format")
    parser.add_argument("-o", dest="output_dir", default="../outputs", help="Output directory path")
    parser.add_argument("-q", dest="quality", default=None, help="Lower quality threshold")
    parser.add_argument("--get-mutational-signatures", dest="mut", default=None, help="Run SigProfilerExtractor on vcf input")
    parser.add_argument("--bin", dest="bins", default=20, type=int, help="Number of bins for each chromosome for visualisation purposes")

    args = parser.parse_args()
    vcf = args.vcf
    output = args.output_dir
    name = args.name
    quality = args.quality
    mut = args.mut
    bin = args.bins

    return (vcf, output, name, quality, mut, bin)


if __name__ == "__main__":
    VCF, OUTPUT, NAME, QUALITY, MUT, BIN = get_args()

    if NAME is None:
        NAME = utils.get_basename(VCF)

    # checks if output exists or create one
    output_dir = utils.check_output_dir(OUTPUT)
    # creates path for vcf
    output_vcf = output_dir + NAME + "_tr.csv"
    # sets location for plot output, checks if exists or creates one
    output_plot_dir = utils.check_output_dir(output_dir + '/Plot_outputs')
    # creates directory for bin plots by chromosome
    if not (os.path.exists(output_plot_dir + '/Bin_by_Chr_Plots')):
        os.mkdir(output_plot_dir + '/Bin_by_Chr_Plots')
        logging.info('Creating output directory')

    ########## BASIC VCF INFO ##########
    pvd.read_basics_vcf(VCF)
    vcf_temp = pvd.prepare_vcf_df(VCF, output_vcf)

    ########## LABEL SUBSTITUTIONS ##########
    output_mut = output_dir + NAME + "mut_sign_SBS.csv"
    vcf_mut = pvd.get_mutational_signature(vcf_temp, output_mut, QUALITY)

    ########## RUN SIGPROEX FOR MUT SIGNATURE ##########
    if MUT is not None:
        sig_out = "../SigProEx_outputs"
        sig_in = "../inputs"
        if not (os.path.exists(sig_out)):
            os.mkdir(sig_out)
        mut.get_signatures(sig_in, sig_out)
    
    ########## BIN GENOME ##########
    output_bin = output_dir + NAME + "_bins.csv"
    vcf_bin = pvd.bin_genome(vcf_temp, output_bin, BIN)

    ########## GENERATE PLOTS ##########
    va.count_mutations(vcf_temp, output_plot_dir)
    va.substitution_type(vcf_mut, output_plot_dir)
    va.variant_distribution(vcf_bin, output_plot_dir)
    