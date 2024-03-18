'''
Determine the mutational signature of the vcf file using SigProfilerExtractor
'''

from SigProfilerExtractor import sigpro as sig

# Run the command above if executing for the first time
'''
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install('GRCh38')
'''

def get_signatures(input_dir, output_dir):
    # to get input from vcf files
    path_to_vcf_file = input_dir
    data = path_to_vcf_file
    sig.sigProfilerExtractor("vcf", output_dir, data, reference_genome="GRCh38", minimum_signatures=1, maximum_signatures=3)


