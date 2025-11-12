"""
A module to supervise each step of the program.
"""

import os
from argparse import Namespace as Arguments
from pathlib import Path
from typing import Final

from loguru import logger
from rich.panel import Panel

import exceptions
import files as io
from callers import VariantCallerRepository
from console import print_stdout
from variants import VariantsRepository


def supervisor(context: Arguments) -> int:
    """
    Merge and reconcile variant calls from multiple callers.

    Args:
        context: Command line parameters.

    Raises:
        SystemExit: If errors occur during the execution.
    """

    # ===========================================================================================
    # Initiate constant variables
    # ===========================================================================================
    # Set the strand bias metric
    SBM: Final[float] = 2.0 if context.disable_strand_bias else 0.95

    # Trace the thresholds
    logger.debug(f"Thresholds: {context.thresholds}")

    # Create a variant caller repository
    # This repository will be used to check if the variant callers are supported
    # and to get the variant caller object for each variant caller
    callers: Final[VariantCallerRepository] = VariantCallerRepository()

    # Create a VCF repository
    # This repository will be used to store the VCF files
    vcfs: Final[io.VCFRepository] = io.VCFRepository()

    # Create a variants repository
    # This repository will be used to store the variants and their informations
    variants: Final[VariantsRepository] = VariantsRepository(
        sample=context.sample,
        rescue=context.rescue,
        intermediate_results=(context.output if context.intermediate_results else ""),
    )

    # ===========================================================================================
    # Check mandatory options and arguments
    # ===========================================================================================

    # Check if the output directory exists
    if not Path(context.output).is_dir():
        logger.error(f"No such output directory: '{context.output}'")
        raise SystemExit(f"No such output directory: '{context.output}'")
    # Check if the output directory is writable
    else:
        if not os.access(context.output, os.W_OK):
            logger.error(
                f"Write permissions are not granted for the output directory: {context.output}"
            )
            raise SystemExit(
                f"Write permissions are not granted for the output directory: {context.output}"
            )

    # Check if the reference genome index is valid
    try:
        # Create a fasta index object
        fai: io.FastaIndex = io.FastaIndex(path=context.reference, lazy=False)
        # Trace the success
        logger.success(
            f"Fasta index {context.reference} has been successfully checked."
        )
    except exceptions.FastaIndexError as e:
        logger.error(f"{context.reference} is not a valid FASTA index: {e}")
        raise SystemExit(f"{context.reference} is not a valid FASTA index") from e

    # Check if the pileup is valid
    try:
        # Create a pileup object
        pileup: io.Pileup = io.Pileup(
            path=context.pileup, sample=context.sample, lazy=True
        )
        # Trace the success
        logger.success(f"Pileup {context.pileup} has been successfully checked.")
        # Store the pileup object in the variants repository
        variants.pileup = pileup
    # Catch an error if the pileup is not valid
    except exceptions.PileupError as e:
        logger.error(f"{context.pileup} is not a valid PILEUP file: {e}")
        raise SystemExit(f"{context.pileup} is not a valid PILEUP file") from e

    # Check if the VCFs are valid
    # Iterator is a list of lists with metadatas about the VCFs: [id, path, [yaml]]
    for vcf in context.vcfs:
        # Check if a YAML config file is provided (non-builtin variant callers)
        if len(vcf) == 3:
            # Store the id, path and yaml config file
            id: str  # Variant caller identifier
            path: str  # Path to the VCF file
            yaml: str  # Path to the YAML config file
            id, path, yaml = vcf
            # Trace
            logger.debug(f"YAML config file {yaml} provided for the VCF {path}")
            try:
                # Create a config file object
                config_file: io.Config = io.Config(path=yaml, lazy=False)
            except exceptions.ConfigError as e:
                raise SystemExit(f"{yaml} is not a valid YAML config file: {e}") from e
            try:
                # Create a variant caller plugin object
                plugin: io.VariantCallerPlugin = io.VariantCallerPlugin(
                    id=id, config=config_file
                )
                # Add the plugin to the caller repository
                callers.add(plugin)
            # Raise an error if the variant caller plugin is not valid
            except exceptions.VariantCallerPluginError as e:
                raise SystemExit(
                    f"{id} is not a valid variant caller plugin: {e}"
                ) from e
        # If no YAML config file is provided (builtin variant callers)
        else:
            # Store the id and path of the VCF
            id: str  # Variant caller identifier
            path: str  # Path to the VCF file
            id, path = vcf
        # Check with identifier if the variant caller is supported
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
            if isinstance(e, exceptions.VCFError):
                # Trace the error
                logger.error(f"{vcf} is not a valid VCF: {e}")
                raise SystemExit(f"{vcf} is not a valid VCF") from e
            # If the error is a VariantCallerError, means that the variant caller is not supported
            else:
                # Trace the error
                logger.error(f"{id} is not a supported variant caller: {e}")
                raise SystemExit(f"{id} is not a supported variant caller") from e

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
    logger.debug("Collecting all variants...")
    # Try to populate the variants repository
    try:
        variants.populate(vcfs=vcfs)
    except exceptions.VariantCallerPluginError as e:
        raise SystemExit("Error while collecting all variants") from e

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
    logger.debug("Calculation of final metrics...")
    # Normalize the variants
    variants.normalize(
        thresholds=context.thresholds,
        length_indels=context.length_indels,
        sbm=SBM,
        sbm_homozygous=context.sbm_homozygous,
    )

    # ===========================================================================================
    # Write VCF(s)
    # ===========================================================================================

    # Create a genomic writter object to write the VCF file
    writter: io.GenomicWritter = io.GenomicWritter(process=0)

    # If the intermediate results and the rescue option are enabled
    if context.intermediate_results and context.rescue:

        # Trace
        logger.debug((f"Writting VCF file of rejected variants in {context.output}."))

        # Write the VCF file of rejected variants
        writter.write(
            output=context.output,
            template="vcf",
            collection=variants.repository,
            lookups=variants.rejected_variants,
            sample=variants.sample,
            contigs=fai.contigs,
            thresholds=context.thresholds,
            suffix="rejected",
        )

        # Trace the success
        logger.success(
            f"VCF file of rejected variants successfully written to {context.output}"
        )

    # Trace
    logger.debug(f"Writting VCF file in {context.output}.")

    # Write the VCF file of common and complex variants
    writter.write(
        output=context.output,
        template="vcf",
        collection=variants.repository,
        lookups=variants.common_variants | variants.complex_variants,
        sample=variants.sample,
        contigs=fai.contigs,
        thresholds=context.thresholds,
    )

    # Trace the success
    logger.success(f"VCF file successfully written to {context.output}")

    # Print the success message to standard output stream
    print_stdout(
        Panel.fit(
            f"VCF successfully generated at '{context.output}'.",
            title="Success",
            highlight=True,
        )
    )

    # Return 0 (success) as Unix convention
    return 0
