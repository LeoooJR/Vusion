from callers import VariantCallerRepository
from console import stdout_console
from rich.panel import Panel
import exceptions as exceptions
import files as io
from loguru import logger
import os
from variants import VariantsRepository

def supervisor(params: object) -> None:
    """
    Combine variant calls from multiple callers.
    
    Args:
        params: Command line parameters.
        
    Raises:
        SystemExit: If errors occur during the execution.
    """

    # ===========================================================================================
    # Initiate constant variables
    # ===========================================================================================
    # Set the strand bias metric
    SBM: float = 2.0 if params.disable_strand_bias else 0.95

    # Trace the thresholds
    logger.debug(f"Thresholds: {params.thresholds}")

    # Create a variant caller repository
    # This repository will be used to check if the variant callers are supported
    # and to get the variant caller object for each variant caller
    callers = VariantCallerRepository()

    # Create a VCF repository
    # This repository will be used to store the VCF files
    vcfs = io.VCFRepository()

    # Create a variants repository
    # This repository will be used to store the variants and their informations
    variants = VariantsRepository(sample=params.sample, rescue=params.rescue, intermediate_results=(params.output if params.intermediate_results else ''))

    # ===========================================================================================
    # Check mandatory options and arguments
    # ===========================================================================================

    # Check if the output directory exists
    if not os.path.isdir(params.output):
        logger.error(f"No such directory: '{params.output}'")
        raise SystemExit(f"No such directory: '{params.output}'")
    # Check if the output directory is writable
    else:
        if not os.access(params.output, os.W_OK):
            logger.error(f"Write permissions are not granted for the directory: {params.output}")
            raise SystemExit(f"Write permissions are not granted for the directory: {params.output}")

    # Check if the reference genome index is valid
    try:
        # Create a fasta index object
        fai = io.FastaIndex(path=params.reference, lazy=False)
        # Trace the success
        logger.success(f"Fasta index {params.reference} has been successfully checked.")
    except exceptions.FastaIndexError as e:
        logger.error(f"{params.reference} is not a valid FASTA index: {e}")
        raise SystemExit(f"{params.reference} is not a valid FASTA index: {e}")
    
    # Check if the pileup is valid
    try:
        # Create a pileup object
        pileup = io.Pileup(path=params.pileup, sample=params.sample, lazy=True)
        # Trace the success
        logger.success(f"Pileup {params.pileup} has been successfully checked.")
        # Store the pileup object in the variants repository
        variants.pileup = pileup
    # Catch an error if the pileup is not valid
    except exceptions.PileupError as e:
        logger.error(f"{params.pileup} is not a valid PILEUP: {e}")
        raise SystemExit(f"{params.pileup} is not a valid PILEUP: {e}")
    
    # Check if the VCFs are valid
    for vcf in params.vcfs:
        # Check if a YAML config file is provided
        if len(vcf) == 3:
            # Store the id, path and yaml config file
            id, path, yaml = vcf
            # Trace
            logger.debug(f"YAML config file {yaml} provided for the VCF {path}")
            # Create a config file object
            config_file = io.Config(path=yaml, lazy=False)
            try:
                # Create a variant caller plugin object
                plugin = io.VariantCallerPlugin(id=id, config=config_file)
                # Add the plugin to the caller repository
                callers.add(plugin)
            # Raise an error if the variant caller plugin is not valid
            except exceptions.VariantCallerPluginError as e:
                raise SystemExit(e)
        # If no YAML config file is provided
        else:
            # Store the id and path of the VCF
            id, path = vcf
        # Check if the variant caller is supported
        if not callers.is_supported(id):
            logger.error(f"{id} variant caller is not supported in --vcf option.")
            raise SystemExit(f"{id} variant caller not supported in --vcf options.")
        # Try to create a VCF object
        try:
            vcfs.add(item=(id, io.VCF(path=path, caller=callers.get_VC(id), lazy=True)))
            # Trace the success
            logger.debug(f"Variant Callers inputed: {id}")
        except (exceptions.VCFError, exceptions.VariantCallerError) as e:
            # If the error is a VCFError, means that the VCF is not valid
            if isinstance(e,exceptions.VCFError):
                # Trace the error
                logger.error(f"{vcf} is not a valid VCF: {e}")
                raise SystemExit(f"{vcf} is not a valid VCF: {e}")
            # If the error is a VariantCallerError, means that the variant caller is not supported
            else:
                # Trace the error
                logger.error(f"{id} is not a supported variant caller: {e}")
                raise SystemExit(f"{id} is not a supported variant caller: {e}")

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
    # Trace
    logger.debug("Collecting all variants.")
    # Try to populate the variants repository
    try:
        variants.populate(vcfs=vcfs)
    except exceptions.VariantCallerPluginError as e:
        raise SystemExit(e)

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
    # Trace
    logger.debug("Calculation of final metrics.")
    # Normalize the variants
    variants.normalize(thresholds=params.thresholds, length_indels=params.length_indels, sbm=SBM, sbm_homozygous=params.sbm_homozygous)

    # ===========================================================================================
    # Write VCF(s)
    # ===========================================================================================

    # Create a genomic writter object to write the VCF file
    writter: io.GenomicWritter = io.GenomicWritter(process=0)

    # If the intermediate results and the rescue option are enabled
    if params.intermediate_results and params.rescue:

        # Trace
        logger.debug((f"Writting VCF file of rejected variants in {params.output}."))

        # Write the VCF file of rejected variants
        writter.write(output=params.output, 
                    template="vcf", 
                    collection=variants.repository, 
                    lookups=variants.rejected_variants, 
                    sample=variants.sample, 
                    contigs=fai.contigs, 
                    thresholds=params.thresholds,
                    suffix="rejected")
        
        # Trace the success
        logger.success(f"VCF file of rejected variants successfully written to {params.output}")

    # Trace
    logger.debug(f"Writting VCF file in {params.output}.")

    # Write the VCF file of common and complex variants
    writter.write(output=params.output, 
                  template="vcf", 
                  collection=variants.repository, 
                  lookups=variants.common_variants | variants.complex_variants, 
                  sample=variants.sample, 
                  contigs=fai.contigs, 
                  thresholds=params.thresholds)

    # Trace the success
    logger.success(f"VCF file successfully written to {params.output}")

    # Print the success message to standard output stream
    stdout_console.print(Panel.fit(f"VCF successfully generated at '{params.output}'.", title="Success", highlight=True), style="result")