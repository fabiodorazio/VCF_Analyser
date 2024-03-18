import pandas as pd; import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns; sns.set_theme(color_codes=True)

sns.set_palette('Set2')

def count_mutations(vcf, output_plot_dir):
    '''
    Generates barplots for counts of genotypes per Mb
    Generates barplots for count of SNP, DEL and INS total and by patient
    Generates barplots of mutation counts per Mb by patient
    '''
    # FIGURE 1
    # slices the df and sums values in columns. Divides by 1m to obtain the count per Mb
    gen_count = vcf[['HOM_REF_COUNT', 'HOM_ALT_COUNT', 'HET_COUNT']].sum() / 1000000
    # creates a df and reset the index
    gen_count = pd.DataFrame(gen_count, columns=['Count / Mb']).reset_index()
    # renames columns
    gen_count.columns = ['GENOTYPE', 'Count / Mb']
    # plot
    sns.barplot(gen_count, x='GENOTYPE', y='Count / Mb')
    plt.savefig(f'{output_plot_dir}/Figure1.png')
    plt.close()

    # FIGURE 2
    # Groups by variant type and sums values. Divides by 1m to obtain the count per Mb
    var_type = vcf.groupby('TYPE').size() / 1000000
    # creates a df and reset the index
    var_type = pd.DataFrame(var_type, columns=['Var Type / Mb']).reset_index()
    # plot
    sns.barplot(var_type, x='TYPE', y='Var Type / Mb')
    plt.savefig(f'{output_plot_dir}/Figure2.png')
    plt.close()

    # FIGURE 3
    # Length of INDELS
    vcf_len = vcf[['TYPE', 'LENGTH']]
    # discard SNPs
    vcf_len = vcf_len[vcf_len['TYPE'].isin(['DEL','INS'])]
    # plot
    sns.histplot(vcf_len, x='LENGTH', hue='TYPE')
    plt.savefig(f'{output_plot_dir}/Figure3.png')
    plt.close()

    # FIGURE 4
    # slice the df
    mut_count = vcf[vcf.columns[-3:]].sum() / 1000000
    mut_count = pd.DataFrame(mut_count, columns=['MUTATION COUNT / Mb']).reset_index()
    mut_count.columns = ['Sample', 'MUTATION COUNT / Mb']
    sns.barplot(mut_count, x='Sample', y='MUTATION COUNT / Mb')
    plt.savefig(f'{output_plot_dir}/Figure4.png')
    plt.close()


def substitution_type(vcf_mut, output_plot_dir):
    '''
    This function generates a barplot of % of SBS subsitution in each patient for the 6 canonical categories
    '''
    # FIGURE 5
    # slices the df
    subset_vcf = vcf_mut[vcf_mut.columns[-4:]]
    # # groups by mutational signature and sum the values in each category
    s = subset_vcf.groupby('MUTATION SIGN').sum()
    s.reset_index(inplace=True)
    # melts the df for plotting
    s = pd.melt(s, id_vars='MUTATION SIGN')
    # calculates percentage of base substitution
    s['Percentage of SBS substitutions'] = s['value'].apply(lambda x:(x/s.sum()['value'])*100)
    # generates plot
    #sns.set_palette('Set2')
    sns.barplot(s, x='MUTATION SIGN', y='Percentage of SBS substitutions',
            hue='variable')
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    plt.savefig(f'{output_plot_dir}/Figure5.png', bbox_inches='tight')
    plt.close()

    
def variant_distribution(vcf_bin, output_plot_dir):
    '''
    This function generates a barplot of mutation per Mb for individual bins
    Iterates over chromosomes and calculates the sum of mutations by patient in each bin
    '''
    # FIGURE 6
        # group by bin and sum count of variants for each patient
    # iterate over chromosomes
    #vcf_bin = pd.read_csv(vcf_bin)
    for chr in vcf_bin['CHROM'].unique():
        # slices the df on unique chromosomes and selects only columns for count
        vcf_bin_chr = vcf_bin[vcf_bin['CHROM'] == chr]
        vcf_bin_chr = vcf_bin_chr[vcf_bin_chr.columns[-4:]]
        # groups by bin and calculates the sum
        vcf_bin_chr = vcf_bin_chr.groupby('BIN', sort=False).sum()
        vcf_bin_chr.reset_index(inplace=True)
        # melts the df for plotting
        vcf_bin_chr = pd.melt(vcf_bin_chr, id_vars='BIN')
        # rename columns
        vcf_bin_chr.columns = ['BIN', 'sample', 'Count']
        # Plot
        sns.barplot(vcf_bin_chr, x='BIN', y='Count', hue='sample')
        plt.xticks(rotation=90)
        plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
        plt.savefig(f'{output_plot_dir}/Bin_by_Chr_Plots/chr{str(chr)}_Figure6.png', bbox_inches='tight')
