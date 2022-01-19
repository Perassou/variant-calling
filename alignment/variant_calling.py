#! /usr/bin/env python3

'''
Author: Marie Peras
Date: 01/19/2022

Variant calling pipeline: the workflow takes short read data (fastq files) as input and outputs variant calls (a vcf file)
'''

import argparse
import logging
from subprocess import call
import sys

logging.basicConfig(level=logging.INFO,
                    format='[%(filename)s:%(lineno)s - %(funcName)s()] [%(asctime)s %(levelname)s] %(message)s')
baselog = logging.getLogger(__name__)
baselog.setLevel(logging.INFO)

def run_bwa(ref_fa: str, read1_fq: str, read2_fq: str, output_name: str) -> str: 
    '''
    Run BWA mem program 
    :param ref_fa: str
        Path to the indexed reference to use for the alignment 
    :param read1_fq: str
        Path to the forward fastq file 
    :param read2_fq: str
        Path to the reverse fastq file 
    :param output_name: str
        Name of the output sam file
    '''
    logger = baselog.getChild('run BWA mem')
    logger.info(f'Launching BWA mem on: {read1_fq} and {read2_fq} using the {ref_fa} as reference')
    call('bwa mem {} {} {} > {}'.format(ref_fa, read1_fq, read2_fq, output_name), shell=True)
    return output_name

def samtools_view(sam_file: str, output_name: str) -> str:
    '''
    Run Samtools view program to convert a bam file to a sam file 
    :param sam_file: str
        Path to the sam file 
    :param output_name: str
        Name of the output sam file
    '''
    logger = baselog.getChild('run samtools view')
    logger.info(f'Convert the sam file: {sam_file} into a bam file: {output_name}')
    call('samtools view -S -b {} > {}'.format(sam_file, output_name), shell=True)
    return output_name

def samtools_sort(bam_file: str, output_name: str) -> str:
    '''
    Run Samtools sort program to sort the alignments from the bam file
    :param bam_file: str
        Path to the bam file to sort
    :param output_name: str 
        Name of the sorted bam file 
    '''
    logger = baselog.getChild('run samtools sort')
    logger.info(f'Sort the alignments from {bam_file} and generate {output_name}')
    call('samtools sort {} -o {}'.format(bam_file, output_name), shell=True)
    return output_name

def run_freebayes(reference: str, aligned_bam: str, output_name: str) -> str:
    '''
    Run Freebayes program and generate the vcf file 
    :param reference: str
        Path to the reference fasta file
    :param aligned_bam: str
        Path to the aligned bam file 
    :param output_name: str 
        Name of the vcf file 
    '''
    logger = baselog.getChild('run freebayes')
    logger.info(f'RUn freebayes on {aligned_bam} using the reference {reference }and generate {output_name}')
    call('freebayes -f {} {} > {}'.format(reference, aligned_bam, output_name), shell = True)
    return output_name


def run_pipeline(args):
    '''
    Run the different programs to generate a vcf file 
    '''
    reference = args.index
    ref = args.reference 
    read1_fq = args.read1_fq
    read2_fq = args.read2_fq
    sam_file = run_bwa(reference, read1_fq, read2_fq, 'aln-pe.sam')
    bam_file = samtools_view(sam_file, 'aln-pe.bam')
    bam_file_sorted = samtools_sort(bam_file, 'aln-pe_sorted.bam')
    vcf = run_freebayes(ref, bam_file_sorted, 'output.vcf')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run variant calling pipeline')
    parser.add_argument('-i', '--index', dest='index', help='Index from the bwa index command', required=True)
    parser.add_argument('-1', '--read1_fq', dest='read1_fq', help='Forward read in a fastq format', required=True)
    parser.add_argument('-2', '--read2_fq', dest='read2_fq', help='Reverse read in a fastq format', required=True)
    parser.add_argument('-r', '--reference', dest='reference', help='Reference genome in a fasta format', required=True)

    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout, level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')

    run_pipeline(args)