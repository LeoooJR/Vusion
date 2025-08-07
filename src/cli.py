from __init__ import __version__
from argparse import ArgumentParser
from loguru import logger
from supervisor import supervisor
from sys import argv
import validation

class EntryPoint:
    """
    Entry point for the program.
    
    This class is responsible for parsing the command line arguments and launching the program.
    """

    FUNC = {"call": supervisor}

    def __init__(self):

        # Create the parser for the command line arguments
        self.parser = ArgumentParser(prog="Vusion", 
                                     description="Combine multiple VCF files.")
        
        # Add the version argument
        self.parser.add_argument(
            "-v",
            "--version",
            action="version",
            version=f"Vusion v{__version__}",
        )

        # Add the reference argument
        self.parser.add_argument(
            "-r",
            "--reference",
            type=str,
            dest="reference",
            metavar="FAI",
            required=True,
            help="Path to reference genome fasta index file (.fai)"
        )

        # Add the output argument
        self.parser.add_argument(
            "-o",
            "--output",
            type=str,
            dest="output",
            metavar="PATH",
            required=True,
            help="Path to the output directory",
        )

        # Add the VCF argument
        self.parser.add_argument(
            "-V",
            "--vcf",
            type=str,
            dest="vcfs",
            nargs="+",
            action=validation.ValidateVCFSAction,
            required=True,
            metavar="VCF",
            help="<ID,vcf,[yaml]> path to \
                  VCF file associated with identifier : bcftools (BT), varscan (VS), vardict (VD), pindel (PL), haplotypecaller (HC), FILT3R (FL), deepvariant (DV) \
                  control & hotspot (CS & HS). \
                  Provide a YAML file path for non-integrated variant callers.",
        )

        # Add the pileup argument
        self.parser.add_argument(
            "-p",
            "--pileup",
            dest="pileup",
            type=str,
            required=True,
            metavar="PILEUP",
            help="Path to pileup processed mpileup data file",
        )

        # Add the sample argument
        self.parser.add_argument(
            "-s",
            "--sample",
            type=str,
            dest="sample",
            metavar="STR",
            required=True,
            help="The sample identifier as specified in both bam and vcf files",
        )

        # Add the thresholds argument
        self.parser.add_argument(
            "-t",
            "--thresholds",
            type=str,
            dest="thresholds",
            default=[10, 30, 40, 60, 70, 80, 20, 30, 50, 1],
            metavar="INT [INT ...]",
            required=False,
            action=validation.ValidateThresholdsAction,
            help="Threshold used for variant categorization in VCF callsets",
        )

        # Add the disable strand bias argument
        self.parser.add_argument(
            "--disable-strand-bias",
            action="store_true",
            dest="disable_strand_bias",
            default=False,
            required=False,
            help="Disable strand bias based filtering of variants",
        )

        # Add the rescue argument
        self.parser.add_argument(
            "-R",
            "--rescue",
            action="store_true",
            default=False,
            required=False,
            help="Rescue rejected calls",
        )

        # Add the hotspot argument
        self.parser.add_argument(
            "-H",
            "--hotspot",
            action="store_true",
            default=False,
            required=False,
            help="Combine two combined callsets; One originated from a hotspot feature list",
        )

        # Add the length indels argument
        self.parser.add_argument(
            "-l",
            "--length-indels",
            dest="length_indels",
            required=False,
            default=1,
            action=validation.ValidateIndelsLengthAction,
            metavar="INT",
            help="Set Del/Ins mininmum length for variant not found in pileup. \
                  Variants below this value are considered false positives. \
                  By default, minimum length is set to 1",
        )

        # Add the sbm homozygous argument
        self.parser.add_argument(
            "--sbm-homozygous",
            dest="sbm_homozygous",
            required=False,
            default=(2 / 3),
            action=validation.ValidateSBMAction,
            metavar="FLOAT",
            help="Define sbm limits for homozygous variants. \
                By default, a strand bias is considered present if there is a 2/3 imbalance \
                of reads on one strand and 1/3 on the other.",
        )

        # Add the intermediate results argument
        self.parser.add_argument(
            "--intermediate-results",
            action="store_true",
            dest="intermediate_results",
            default=False,
            required=False,
            help="Should intermediate results be saved."
        )

        # Add the debug argument
        self.parser.add_argument(
            "-d",
            "--debug",
            dest="debug",
            action="store_true",
            help="Should logs be saved?",
            default=False,
        )

        # Set the default function to call
        self.parser.set_defaults(func=self.FUNC["call"])

    def launch(self) -> int:
        """
        Launch the program with command line arguments.
        
        Returns:
            int: Exit code of the program.
        """

        # Parse the command line arguments
        cmd = self.parser.parse_args(argv[1:])
        
        # Remove the default handler
        logger.remove(0)

        # Should the log be saved in a file ?
        if cmd.debug:

            logger.add("vusion.log")

        return cmd.func(params=cmd)
    
    def __str__(self):
        """
        String representation of the program.
        
        Returns:
            str: The name of the program.
        """

        return "Vusion"

    def __repr__(self):
        """
        String representation of the program for debugging.
        
        Returns:
            str: The name of the program.
        """
        
        return "Vusion"
