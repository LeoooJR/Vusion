from __init__ import __version__
from argparse import ArgumentParser, ArgumentTypeError
from loguru import logger
from vusion import combine
from sys import argv
import validation

class Program:

    FUNC = {"combine": combine}

    def __init__(self):

        self.parser = ArgumentParser(prog="Vusion", 
                                     description="Combine multiple VCF files.")
        self.parser.add_argument(
            "-v",
            "--version",
            action="version",
            version=f"Vusion v{__version__}",
        )
        self.parser.add_argument(
            "-r",
            "--reference",
            type=str,
            dest="reference",
            metavar="FAI",
            required=True,
            help="Path to reference genome fasta index file (.fai)"
        )
        self.parser.add_argument(
            "-o",
            "--output",
            type=str,
            dest="output",
            metavar="PATH",
            required=True,
            help="Path to the output directory",
        )
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
        self.parser.add_argument(
            "-p",
            "--pileup",
            dest="pileup",
            type=str,
            required=True,
            metavar="PILEUP",
            help="Path to pileup processed mpileup data file",
        )
        self.parser.add_argument(
            "-s",
            "--sample",
            type=str,
            dest="sample",
            metavar="STR",
            required=True,
            help="The sample identifier as specified in both bam and vcf files",
        )
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
        self.parser.add_argument(
            "-d",
            "--disable-strand-bias",
            action="store_true",
            dest="disable_strand_bias",
            default=False,
            required=False,
            help="Disable strand bias based filtering of variants",
        )
        self.parser.add_argument(
            "-R",
            "--rescue",
            action="store_true",
            default=False,
            required=False,
            help="Rescue rejected calls",
        )
        self.parser.add_argument(
            "-H",
            "--hotspot",
            action="store_true",
            default=False,
            required=False,
            help="Combine two combined callsets; One originated from a hotspot feature list",
        )
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
        self.parser.add_argument(
            "--intermediate-results",
            action="store_true",
            dest="intermediate_results",
            default=False,
            required=False,
            help="Should intermediate results be saved."
        )
        self.parser.add_argument(
            "--verbosity",
            action="store_true",
            default=False,
            required=False,
            help="Should logs be printed to the shell.",
        )

        self.parser.set_defaults(func=self.FUNC["combine"])

    def launch(self) -> int:

        cmd = self.parser.parse_args(argv[1:])

        # Should the log be printed to CLI or saved in a file ?
        if not cmd.verbosity:

            logger.remove(0)

            logger.add("vusion.log")

        return cmd.func(params=cmd)
    
    def __str__(self):

        return "Vusion"

    def __repr__(self):
        
        return "Vusion"
