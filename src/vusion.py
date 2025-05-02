#!/usr/bin/python3

from callers import VariantCallerRepository
import errors
import files as io
from loguru import logger
from variants import VariantsRepository

def combine(params):

    # ===========================================================================================
    # Initiate constant variables
    # ===========================================================================================
    # Set the strand bias metric
    SBM: float = 2.0 if params.disable_strand_bias else 0.95
    # Maximum and minimum threshold values
    MAX_THRESHOLD: float = 100.0
    MIN_THRESHOLD: float = 0.0

    # Create a variant caller repository
    # This repository will be used to check if the variant callers are supported
    # and to get the variant caller object for each variant caller
    callers = VariantCallerRepository()

    # Create a variants repository
    # This repository will be used to store the variants and their information
    variants = VariantsRepository(sample=params.sample, rescue=params.rescue)

    # ===========================================================================================
    # Check mandatory options
    # ===========================================================================================

    # Check reference genome index
    try:
        fai = io.FastaIndex(path=params.reference, lazy=False)
    except errors.FastaIndexError as e:
        logger.error(f"{params.reference} is not a valid FASTA index.")
        logger.error(f"Error: {e}")
        raise SystemExit(f"{params.reference} is not a valid FASTA index.")
    
    logger.success(f"Fasta index {params.reference} has been successfully checked.")
    
    # Check pileup
    try:
        pileup = io.Pileup(path=params.pileup, sample=params.sample, lazy=True)
        variants.pileup = pileup
    except errors.PileupError as e:
        logger.error(f"{params.pileup} is not a valid PILEUP.")
        logger.error(f"Error: {e}")
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
            if not callers.is_supported(input[0]):
                logger.error(f"{input[0]} caller is not supported in --vcf option.")
                raise SystemExit("Caller not supported in --vcf options.")
            
        try:
            vcfs[input[0]] = {"vcf": io.VCF(path=input[1], caller=callers.get_VC(input[0]), lazy=True), "index": None}
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
    logger.debug("Collecting all variants.")
    
    variants.populate(vcfs=vcfs)

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
    logger.debug("Calculation of final metrics.")

    variants.normalize(pileup=pileup, thresholds=thresholds, length_indels=params.length_indels, sbm=SBM, sbm_homozygous=params.sbm_homozygous)

    # ===========================================================================================
    # Write VCF(s)
    # ===========================================================================================

    # Create a writer object to write the VCF file
    writter: io.GenomicWritter = io.GenomicWritter(process=0)

    if params.intermediate_results and params.rescue:

        logger.debug((f"Writting VCF file of rejected variants in {params.output}."))

        # Write the VCF file
        writter.write(output=params.output, 
                    template="vcf", 
                    collection=variants.repository, 
                    lookups=variants.rejected_variants, 
                    sample=variants.sample, 
                    contigs=fai.contigs, 
                    thresholds=thresholds,
                    suffix="rejected")
        
        logger.success(f"VCF file of rejected variants successfully written to {params.output}")

    logger.debug(f"Writting VCF file in {params.output}.")

    # Write the VCF file
    writter.write(output=params.output, 
                  template="vcf", 
                  collection=variants.repository, 
                  lookups=variants.common_variants | variants.complex_variants, 
                  sample=variants.sample, 
                  contigs=fai.contigs, 
                  thresholds=thresholds)

    logger.success(f"VCF file successfully written to {params.output}")