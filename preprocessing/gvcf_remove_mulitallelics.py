


import pysam
import gzip

def extract_genotypes(vcf_file):
    with pysam.VariantFile(vcf_file) as vcf:
        for record in vcf:
            genotypes = []
            for sample in record.samples.values():
                genotype = sample["GT"]
                if genotype is not None:
                    genotypes.extend(genotype.split('/'))
            yield record, genotypes

def filter_and_save_vcf(input_vcf, output_vcf):
    with pysam.VariantFile(input_vcf) as vcf_in, gzip.open(output_vcf, 'wt') as vcf_out:
        header = str(vcf_in.header)
        vcf_out.write(header)
        
        for record, genotypes in extract_genotypes(input_vcf):
            unique_genotypes = set(genotypes)
            if len(unique_genotypes) <= 2:
                vcf_out.write(str(record))

input_vcf_file = 'input.vcf.gz'
output_vcf_file = 'filtered.vcf.gz'

filter_and_save_vcf(input_vcf_file, output_vcf_file)
