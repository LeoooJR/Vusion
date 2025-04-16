#!/usr/bin/python3

"""
HOW TO RUN :

python3 CombineVCF2Leaves-4.py \
-r ~/Bureau/8072/Data/REF/human_g1k_v37_decoy.fasta.fai \
-V ST,45-CEREBMD-12_S12_ST_final.vcf \
-V VD,45-CEREBMD-12_S12_VD_final.vcf \
-V PL,45-CEREBMD-12_S12_PL_final.vcf \
-V HC,45-CEREBMD-12_S12_HC_final.vcf \
-V VS,45-CEREBMD-12_S12_VS_final.vcf \
-V FL,45-CEREBMD-12_S12_VS_FILT3R_filtered.vcf \
--pileup 45-CEREBMD-12_S12_CT.tsv \
-o test.vcf \
--sample_id test \
-t 10,30,40,60,70,80,20,30,50,1 \


DICTIONNARY :
'chr:pos:REF:ALT':{'VC':{
                         'VAF': {'PL': 8.968609865470851},
                         'GT': {'PL': '0/0'},
                         'FILTER': {'PL': 'PASS'},
                         'INFO': {'PL': ...;SVLEN=34;SVTYPE=INS'},
                         'FORMAT': {'PL': 'GT:AD'},
                         'SAMPLE': {'PL': '0/0:203,20'},
                         'TRC': {'PL': 223.0},
                         'ARC': {'PL': 20.0},
                         'RRC': {'PL': 203.0}},
                    'VT': 'INS',
                    'final_metrics': {'TRC': 223, 'TRC-': '-1', 'TRC+': '-1', \
                    'RRC-': '-1', 'RRC+': '-1', 'VCN': 1, 'VCI': 'PL', \
                    'ARR': '8.96861', 'LOW': 0, 'VAR': 'LSC', 'ARC': '-1,-1,20.0',\
                    'RRC': '-1,-1,203.0', 'GT': '0/0', 'BRC': '-1', 'BRR': '-1', \
                    'BRE': '-1', 'BKG': '-1', 'ARC+': '-1', 'ARC-': '-1', 'SBM': '-1',\
                    'SBP': '-1', 'FILTER': 'PASS', 'PIL': 'N'},
                    'vcf_fields': ['VAR', 'BKG', 'TRC', 'RRC', 'ARC', 'BRC', 'ARR', \
                    'BRR', 'BRE', 'SBP', 'SBM', 'LOW', 'VCI', 'VCN','PIL']}}
"""

from callers import VariantCallerRepository
import errors
import files as io
from loguru import logger
import os
import sys
from variants import VariantsRepository

import pprint

import utils as functions

def combine(params):

    # ===========================================================================================
    # Initiate constant variables
    # ===========================================================================================
    SBM = 0.95
    MAX_THRESHOLD = 100.0
    MIN_THRESHOLD = 0.0

    # Create a variant caller repository
    # This repository will be used to check if the variant callers are supported
    # and to get the variant caller object for each variant caller
    caller_repository = VariantCallerRepository()

    # Create a variants repository
    # This repository will be used to store the variants and their information
    variants_repository = VariantsRepository(rescue=params.rescue)

    # ===========================================================================================
    # Check mandatory options
    # ===========================================================================================

    # Check reference genome index
    try:
        fasta_index = io.FastaIndex(path=params.reference, lazy=False)
    except errors.FastaIndexError:
        logger.error(f"{params.reference} is not a valid FASTA index.")
        raise SystemExit(f"{params.reference} is not a valid FASTA index.")
    
    logger.success(f"Fasta index {params.reference} has been successfully checked.")
    
    # Check pileup
    try:
        pileup = io.Pileup(path=params.pileup, lazy=True)
        variants_repository.set_pileup(pileup)
    except errors.PileupError:
        logger.error(f"{params.pileup} is not a valid PILEUP.")
        raise SystemExit(f"{params.pileup} is not a valid PILEUP.")
    
    logger.success(f"Pileup {params.pileup} has been successfully checked.")

    # Check VCFs
    vcfs: dict = {}

    for vcf in params.vcfs:

        input = vcf.split(',')

        if len(input) != 2:
            logger.error('Wrong number of argument in --vcf option.')
            logger.error(f'Error was raised by: {input}.')
            raise ValueError("Wrong number of argument in --vcf option.")
        
        try:
            int(input[0])
            logger.error("Wrong type of argument in --vcf option.")
            logger.error(f"Error was raised by: {input[1]}.")
            raise SystemExit("Wrong type of argument in --vcf option.")
        except ValueError:
            if not caller_repository.is_supported(input[0]):
                logger.error(f"{input[0]} caller is not supported in --vcf option.")
                raise SystemExit("Caller not supported in --vcf options.")
            
        try:
            vcfs[input[0]] = {"vcf": io.VCF(path=input[1], caller=caller_repository.get_VC(input[0]), lazy=True), "index": None}
        except (errors.VCFError, errors.VariantCallerError) as e:
            if isinstance(e,errors.VCFError):
                logger.error(f"{vcf} is not a valid VCF.")
                logger.error(f"Error: {e}")
                raise SystemExit(f"{vcf} is not a valid VCF.")
            else:
                logger.error(f"{input[0]} is not a supported variant caller.")
                logger.error(f"Error: {e}")
                raise SystemExit(f"{input[0]} is not a supported variant caller.")

        logger.debug(f"Variant Callers inputed: {input[0]}")

    # Check if all mandatory option are given and modify variable depending of given options

    if params.disable_strand_bias:
        SBM = 2

    thresholds: list[str] = params.thresholds.split(',')

    # Check that we have 10 values in thresholds
    if len(thresholds) != 10:
        logger.error(f"Invalid number of values in --threshold option.")
        raise SystemExit(f"Invalid number of values in --threshold option.")
        
    #Check that all values can be converted to floats
    try:
        thresholds: list[float] = list(map(float, thresholds)) 
    except ValueError:
        logger.error("Invalid values in thresholds option.")
        raise SystemExit("Invalid values in thresholds option.")
        
    if not any(list(map(lambda value: value >= MIN_THRESHOLD and value <= MAX_THRESHOLD, thresholds))):
        logger.error("Option --thresholds values cannot be equal or higher than 100, or lower than 0.")
        raise SystemExit("Option --thresholds values cannot be equal or higher than 100, or lower than 0.")

    # Check that 6 first given threshold are unique
    if len(thresholds[0:6]) != len(set(thresholds[0:6])):
        logger.error("Option --thresholds six first values must be unique.")
        raise SystemExit("Option --thresholds six first values must be unique.")

    # Check that second group of threshold values are unique
    if len(thresholds[6:9]) != len(set(thresholds[6:9])):
        logger.error("Option --thresholds values 7, 8 and 9 must be unique.")
        raise SystemExit("Option --thresholds values 7, 8 and 9 must be unique.")
    
    # Sort the first 6 values and the last 3 values
    # This is done to make sure that the values are in the right order
    thresholds[0:6] = sorted(thresholds[0:6])
    thresholds[6:9] = sorted(thresholds[6:9])

    logger.debug(f"Thresholds: {thresholds}")

    # ============================================================================================
    # Parse VCFs
    # ============================================================================================

    # Populate the variants repository from the VCFs
    # This will create a dictionary of variants with the following structure:
    # 'chr':{
    #     '(vcf_pos, pileup_pos)':{
    #         'ref:alt':{
    #             'collection':{},
    #             'type': '',
    #             'display': '',
    #         }
    #     }
    # }
    variants_repository.populate(vcfs=vcfs)

    # ===========================================================================================
    # Process variants with Pileup
    # ===========================================================================================

    # Normalize variants with common metrics
    # Use the pileup to normalize the variants
    # This will create a dictionary of variants with the following structure:
    # 'chr':{
    #     '(vcf_pos, pileup_pos)':{
    #         'ref:alt':{
    #             'collection':{},
    #             'type': '',
    #             'display': '',
    #             'filter' 'REJECTED|FAIL|PASS',
    #             'sample':{},
    #         }
    #     }
    # }
    variants, ITD = variants_repository.normalize(sample=params.sample, pileup=pileup, thresholds=thresholds, length_indels=params.length_indels, sbm=SBM, sbm_homozygous=params.sbm_homozygous)

    # ===========================================================================================
    # Process complex variants without Pileup : INV,MNV and CSV
    # ===========================================================================================
    variants: dict = functions.process_without_pileup(variants=variants, lookups=variants_repository.cache["complex"].union(ITD), thresholds=thresholds, sbm=SBM, sbm_homozygous=params.sbm_homozygous)

    # ===========================================================================================
    # rescuing rejected calls w/ FILTER <PASS> (optional; if any)
    # ===========================================================================================
    if params.rescue:

        # ===========================================================================================
        # Process rejected variant
        # ===========================================================================================
        rejected: dict = functions.process_without_pileup(variants=variants, lookups=variants_repository.cache["rejected"], thresholds=thresholds, sbm=SBM, sbm_homozygous=params.sbm_homozygous)

    # ===========================================================================================
    # Write VCFs
    # ===========================================================================================

    # Create a writer object to write the VCF file
    writter: io.GenomicWritter = io.GenomicWritter(file=params.output)

    # write rejected calls to file
    # with open(os.path.join(os.getcwd(), f"{params.sample}_failed.vcf"), mode='w',encoding='utf-8') as OUT_TRASH_FILE:
    #     OUT_TRASH_FILE.write(VCF_HEADER)
    #     ordered_variant_key = functions.order_var(rejected, contigs["contig"])
        
    #     for variant_key in ordered_variant_key:
    #         OUT_TRASH_FILE.write(functions.print_var(variant_key, rejected, 'final_metrics'))

    # Write the VCF file
    writter.writeVCF(contigs=fasta_index.get_contigs(), variants=variants, samples=[params.sample], thresholds=thresholds)
