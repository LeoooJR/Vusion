from __init__ import __version__
from argparse import ArgumentParser
from loguru import logger
from vusion import combine
from sys import argv


class Program:

    FUNC = {"combine": combine}

    def __init__(self):

        self.parser = ArgumentParser(description="Combine multiple VCF files.")
        self.parser.add_argument(
            "-v",
            "--version",
            action="version",
            version=f"AFO-VCF {__version__}",
        )
        self.parser.add_argument(
            "-r",
            "--reference",
            type=str,
            dest="reference",
            metavar="",
            required=True,
            help="<reference.fasta.fai>\
            Path to reference genome fasta index file",
        )
        self.parser.add_argument(
            "-o",
            "--output",
            type=str,
            dest="output",
            metavar="",
            required=True,
            help="Path to the output directory",
        )
        self.parser.add_argument(
            "-V",
            "--vcf",
            type=str,
            dest="vcfs",
            action="append",
            required=True,
            help=" <ID,vcf,[yaml]> path to input \
        VCF file associated with identifier : bcftools (BT), varscan (VS), vardict (VD), pindel (PL), haplotypecaller (HC), FILT3R (FL), deepvariant (DV) \
        control & hotspot (CS & HS). \
        Provide a YAML file path for non-integrated variant callers.",
        )
        self.parser.add_argument(
            "-p",
            "--pileup",
            dest="pileup",
            type=str,
            metavar="",
            required=True,
            help="Path to pileup processed mpileup data file",
        )
        self.parser.add_argument(
            "-s",
            "--sample",
            type=str,
            dest="sample",
            metavar="",
            required=True,
            help="The sample identifier as specified in both bam and vcf files",
        )
        self.parser.add_argument(
            "-t",
            "--thresholds",
            type=str,
            dest="thresholds",
            default="10,30,40,60,70,80,20,30,50,1",
            metavar="",
            required=False,
            help="Threshold used for variant categorization in VCF callsets",
        )
        self.parser.add_argument(
            "-d",
            "--disable_strand_bias",
            action="store_true",
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
        # Si INDEL not in pileup and <l, will be considered as FP
        self.parser.add_argument(
            "-l",
            "--length_indels",
            action="store",
            dest="length_indels",
            required=False,
            default=1,
            type=float,
            help="Set Del/Ins min length (for cases where Del/Ins are not in Pileup)\
                Ex : 3:12580629:C:CT len(Ins) = 1 \
                By default, minimum length is set to 1",
        )
        self.parser.add_argument(
            "--sbm_homozygous",
            action="store",
            required=False,
            default=(2 / 3),
            type=float,
            help="Define sbm limits for homozygous variants. \
                By default, a strand bias is considered present if there is a 2/3 imbalance \
                of reads on one strand and 1/3 on the other strand.",
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

        if not cmd.verbosity:

            logger.remove(0)

            logger.add("vusion.log")

        return cmd.func(params=cmd)
