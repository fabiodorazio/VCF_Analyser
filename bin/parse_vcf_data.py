
import sys
import pysam
from pysam import VariantFile
import pandas as pd
import numpy as np
import logging


def read_basics_vcf(vcf_file):
    '''
    Uses pysam to parse the vcf and print some basic info of the vcf
    '''
    with VariantFile(vcf_file, 'r') as vcf:
        # reports basic info from the vcf
        n_chromosomes = len(list(vcf.header.contigs))
        filters = list(vcf.header.filters)
        info = list(vcf.header.info)
        samples = list(vcf.header.samples)
        print(f"There is {n_chromosomes} in the vcf file, filter classes are: {filters}, Info classes are: {info}, and there are: {samples} samples")
        

def retrieve_sample_columns(vcf_file):
    '''
    Works out the number of samples in the vcf and creates a list of names that will be passed as column names
    '''
    with VariantFile(vcf_file, 'r') as vcf:
        n_samples = len(list(vcf.header.samples))
        # retrieves sample columns
        sample_columns = []
        for sample in range(1,n_samples+1):
            sample_columns.append('sample'+str(sample))
    
        return(sample_columns)


def determine_genotype(genotype, type='het'):
    '''
    Determines the number of hom or hets genotypes for each variant but does not distinguish by sample
    Hom can be ref or alt
    '''
    het = []
    hom_ref = []
    hom_alt = []
    # splits each genotype in list on and checks for equivalence between the elements 
    for element in genotype:
        s_genotype = element.split('|')
        if all(x == '0' for x in s_genotype):
            hom_ref.append(element)
        elif all(x == '1' for x in s_genotype):
            hom_alt.append(element)
        else:
            het.append(element)
    # returns length of list for each variant
    if type == 'het':
        return len(het)
    elif type == 'hom_ref':
        return len(hom_ref)
    elif type == 'hom_alt':
        return len(hom_alt)


def var_length(var):
    '''
    Returns length of deletions or insertions by counting the number of characters in ref or alt respectively
    '''
    if var['TYPE'] == 'DEL':
        return -len(var['REF'])
    elif var['TYPE'] == 'INS':
         return len(var['ALT'])
    else:
        return None

def prepare_vcf_df(vcf_file, outfolder):
    '''
    Performs initial data transformation
    Imports the vcf file as data table 
    Works out the type of variant, their length and genotypes by patient
    '''
    s = retrieve_sample_columns(vcf_file)
    header = "CHROM POS ID REF ALT QUAL FILTER INFO FORMAT".split() + s
    vcf = pd.read_csv(vcf_file, delimiter='\t', comment='#', names=header)
    # creates a dictionary for the info column
    if vcf['INFO'] is not None:
        vcf['INFO'] = vcf['INFO'].str.split(";").apply(lambda x: dict([y.split("=") for y in x]))
        # add TYPE column for types of mutations
        #vcf['TYPE'] = vcf['INFO'].apply(lambda x:x['VT']) # this will only add either INDEL or SNP
        for pos in range(0,vcf.shape[0]):
            ref = vcf.iloc[pos,:].loc['REF']
            alt = vcf.iloc[pos,:].loc['ALT']
            if len(ref) > len(alt):
                vcf.at[pos,'TYPE'] = 'DEL'
            elif len(ref) < len(alt):
                vcf.at[pos,'TYPE'] = 'INS'
            elif ',' in alt:
                vcf.at[pos,'TYPE'] = 'MUL ALT'
            elif alt == '<DEL>':
                vcf.at[pos,'TYPE'] = 'SV'
            elif alt == '<INV>':
                vcf.at[pos,'TYPE'] = 'INV'
            else:
                vcf.at[pos,'TYPE'] = 'SNP'
    # adds LENGTH COLUMN for INS, DEL and SNP
    vcf['LENGTH'] = vcf.apply(lambda x:var_length(x), axis=1)
    # adds GENOTYPE and COUNT columns 
    vcf['GENOTYPE'] = vcf.apply(lambda row: [row['sample1'], row['sample2'], row['sample3']], axis=1)
    vcf['HOM_REF_COUNT'] = vcf['GENOTYPE'].apply(lambda x: determine_genotype(x, type='hom_ref'))
    vcf['HOM_ALT_COUNT'] = vcf['GENOTYPE'].apply(lambda x: determine_genotype(x, type='hom_alt'))
    vcf['HET_COUNT'] = vcf['GENOTYPE'].apply(lambda x: determine_genotype(x, type='het'))

    # adds genetype count by patient
    sample_count = vcf[['sample1', 'sample2', 'sample3']]
    # renames columns
    names = ['sample1_count', 'sample2_count', 'sample3_count']
    sample_count.columns = names
    # replaces genotypes with count: 0 for hom reference and 1 for het or hom alternative 
    sample_count = sample_count.replace({'0|0':0, '0|1':1, '1|1':1, '1|0':1})
    # concatenates data frames
    vcf = pd.concat([vcf, sample_count],axis=1)
    vcf.to_csv(outfolder, index=False)
    #vcf.to_csv('../outputs/vcf.csv', index=False)

    return(vcf)



def get_mutational_signature(vcf, outpath, quality=None):
    '''
    Slices the previous data frame on variant quality/filter
    Subsets only SBS
    Mutational substitutions are generated from ref and alt columns and combined in thr 6 canonical classes where the pyrimidine is mutated
    '''
    # filters dataset based on variant quality and select only SNP
    if quality is not None:
        vcf_pass = vcf[(vcf['QUAL'] > quality) & (vcf['TYPE'] == 'SNP')]
    else:
        vcf_pass = vcf[(vcf['FILTER'] == 'PASS') & (vcf['TYPE'] == 'SNP')]
    # adds mutational substitution
    vcf_pass['MUTATION SIGN'] = vcf_pass['REF'] + '>' + vcf_pass['ALT']
    # combines the substitutions from the alternative strand
    repl_dic = {'G>T': 'C>A', 'G>C':'C>G', 'G>A':'C>T', 'A>T':'T>A', 'A>G':'T>C', 'A>C':'T>G'}
    vcf_pass['MUTATION SIGN'] = vcf_pass['MUTATION SIGN'].replace(repl_dic)
    vcf_pass.to_csv(outpath, index=False)
    #vcf_pass.to_csv('../outputs/vcf_pass.csv', index=False)

    return(vcf_pass)


def bin_genome(vcf, outpath, n_bins):
    '''
    Binning for easier visualisation
    for each chromosome determines the number of bases 
    the dataframe is binned based on POS
    '''
    # creates list of bin labels
    bin_labels = []
    for bin_n in range(1,n_bins + 1):
        bin_lab = 'Bin' + str(bin_n)
        bin_labels.append(bin_lab)
    # iterates over each chromosome
    each_chr = []
    for chr in vcf['CHROM'].unique():
        # slices df for each chromosome
        chr_df = vcf[vcf['CHROM'] == chr]    
        bin = []
        start = 0
        # max size of each chromosome
        tot = vcf.iloc[-1]['POS']
        bins_size = tot // n_bins
        # bins the dataframe based on coordinates in POS
        chr_df['BIN'] = pd.cut(chr_df['POS'], bins=n_bins, labels=bin_labels)
        # appends each chromosome df to list
        each_chr.append(chr_df)

    # concatenate dataframes into a single df
    vcf_binned = pd.concat(each_chr)
    # save
    vcf_binned.to_csv(outpath, index=False)

    return(vcf_binned)



