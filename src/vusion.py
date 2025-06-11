from callers import VariantCallerRepository
import exceptions as exceptions
import files as io
from loguru import logger
import os
from variants import VariantsRepository

try:
    from icecream import ic
    ic.configureOutput(prefix='->', includeContext=False)
except ImportError:  # Graceful fallback if IceCream isn't installed.
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a)  # noqa

def combine(params):

    # ===========================================================================================
    # Initiate constant variables
    # ===========================================================================================
    # Set the strand bias metric
    SBM: float = 2.0 if params.disable_strand_bias else 0.95

    logger.debug(f"Thresholds: {params.thresholds}")

    # Create a variant caller repository
    # This repository will be used to check if the variant callers are supported
    # and to get the variant caller object for each variant caller
    callers = VariantCallerRepository()

    # Create a variants repository
    # This repository will be used to store the variants and their information
    variants = VariantsRepository(sample=params.sample, rescue=params.rescue, intermediate_results=(params.output if params.intermediate_results else ''))

    # ===========================================================================================
    # Check mandatory options
    # ===========================================================================================

    # Check output directory
    if not os.path.isdir(params.output):
        logger.error(f"No such directory: '{params.output}'")
        raise SystemExit(f"No such directory: '{params.output}'")
    else:
        if not os.access(params.output, os.W_OK):
            logger.error(f"Write permissions are not granted for the directory: {params.output}")
            raise SystemExit(f"Write permissions are not granted for the directory: {params.output}")

    # Check reference genome index
    try:
        fai = io.FastaIndex(path=params.reference, lazy=False)
    except exceptions.FastaIndexError as e:
        logger.error(f"{params.reference} is not a valid FASTA index.")
        logger.error(f"{e}")
        raise SystemExit(f"{params.reference} is not a valid FASTA index.")
    
    logger.success(f"Fasta index {params.reference} has been successfully checked.")
    
    # Check pileup
    try:
        pileup = io.Pileup(path=params.pileup, sample=params.sample, lazy=True)
        variants.pileup = pileup
    except exceptions.PileupError as e:
        logger.error(f"{params.pileup} is not a valid PILEUP.")
        logger.error(f"{e}")
        raise SystemExit(f"{params.pileup} is not a valid PILEUP.")
    
    logger.success(f"Pileup {params.pileup} has been successfully checked.")

    # Check VCFs
    vcfs: dict = {}

    for vcf in params.vcfs:
        
        if len(vcf) == 3:
            id, path, yaml = vcf
            logger.debug(f"YAML config file {yaml} provided for the VCF {id}")
            config_file = io.Config(path=yaml, lazy=False)

            try:
                callers.add(id=id, recipe=config_file.params)
            except exceptions.VariantCallerPluginError as e:
                raise SystemExit(e)

        else:
            id, path = vcf

        if not callers.is_supported(id):
            logger.error(f"{id} variant caller is not supported in --vcf option.")
            raise SystemExit(f"{id} variant caller not supported in --vcf options.")
            
        try:
            vcfs[id] = {"vcf": io.VCF(path=path, caller=callers.get_VC(id), lazy=True), "index": None}
        except (exceptions.VCFError, exceptions.VariantCallerError) as e:
            if isinstance(e,exceptions.VCFError):
                logger.error(f"{vcf} is not a valid VCF.")
                logger.error(f"{e}")
                raise SystemExit(f"{vcf} is not a valid VCF.")
            else:
                logger.error(f"{id} is not a supported variant caller.")
                logger.error(f"{e}")
                raise SystemExit(f"{id} is not a supported variant caller.")

        logger.debug(f"Variant Callers inputed: {id}")

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

    variants.normalize(thresholds=params.thresholds, length_indels=params.length_indels, sbm=SBM, sbm_homozygous=params.sbm_homozygous)

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
                    thresholds=params.thresholds,
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
                  thresholds=params.thresholds)

    logger.success(f"VCF file successfully written to {params.output}")