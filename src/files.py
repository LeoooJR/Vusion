import callers
from collections import defaultdict
import errors
import jinja2
import numpy as np
import os
import pandas as pd
import re


class GenomicWritter:

    def __init__(self, file: str):

        self.file = file

    def writeVCF(
        self,
        contigs: object,
        variants: object,
        samples: list[str],
        thresholds: list[float],
    ):

        def sort_variant(contigs: object, variants: dict) -> object:

            # First sort the positions for each chromosome
            for chromosome in variants:

                variants[chromosome] = {
                    k: v
                    for k, v in sorted(
                        variants[chromosome].items(),
                        key=lambda item: int(item[0].vcf_position),
                    )
                }

            # Then sort the chromosomes
            {
                k: v
                for k, v in sorted(variants.items(), key=lambda item: item[0])
            }

        def format_sample(metrics: dict) -> str:

            return ":".join(
                [
                    (
                        ",".join(
                            [
                                str(metrics[f"{format}+"]),
                                str(metrics[f"{format}-"]),
                                (
                                    str(
                                        metrics[f"{format}+"]
                                        + metrics[f"{format}-"]
                                    )
                                    if metrics[f"{format}+"] != -1
                                    and metrics[f"{format}+"] != -1
                                    else str(metrics[f"{format}"])
                                ),
                            ]
                        )
                        if format in ["RRC", "ARC"]
                        else str(metrics[format])
                    )
                    for format in FORMAT
                ]
            )

        FORMAT: list[str] = [
            "GT",
            "VAR",
            "BKG",
            "TRC",
            "RRC",
            "ARC",
            "BRC",
            "ARR",
            "BRR",
            "BRE",
            "SBP",
            "SBM",
            "LOW",
            "VCI",
            "VCN",
            "PIL",
            "RES",
        ]

        HEADER: list[str] = [
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
        ]

        # Add sample names to the header
        HEADER.extend(samples)

        INFOS: list[str] = ["VAR"]

        # Path to the template directory
        ressources = os.path.join(
            os.path.dirname(os.path.abspath(__file__)), "templates"
        )

        env = jinja2.Environment(loader=jinja2.FileSystemLoader(ressources))

        # Load the header template
        template = env.get_template("header")

        with open(self.file, mode="w") as out:

            # Write the rendered header to the file
            out.writelines(
                template.render(contigs=contigs, thresholds=thresholds)
            )

            header: str = "\t".join(HEADER)

            # Write the column names
            out.write(f"\n#{header}\n")

            # Sort the variants by chromosome and position
            sort_variant(contigs=contigs, variants=variants)

            # Iterate over each level of the variants dictionary
            for chromosome in variants:

                for positions in variants[chromosome]:

                    for mutation in variants[chromosome][positions]:

                        variant: dict = variants[chromosome][positions][
                            mutation
                        ]

                        ref, alt = mutation.split(":")

                        if variant.get("filter", "REJECTED") != "REJECTED":

                            # Write ONLY if normalized metrics are present
                            if "sample" in variant:

                                out.write(
                                    "\t".join(
                                        [
                                            f"chr{chromosome}",  # Chromosome field
                                            str(
                                                positions.vcf_position
                                            ),  # Position field
                                            ".",  # ID field
                                            ref,  # Reference field
                                            alt,  # Alternate field
                                            ".",
                                            variant["filter"],  # Filter field
                                            "=".join(
                                                [INFOS[0], variant["type"]]
                                            ),  # Info field
                                            ":".join(FORMAT),  # Format field
                                            format_sample(variant["sample"]),
                                        ]
                                    )
                                )  # Sample values field

                                out.write("\n")


class GenomicReader:

    pass


class GenomicFile:

    def __init__(self, path: str):

        self.path = path

    def is_empty(self) -> bool:
        """Check if file is empty"""
        return os.path.getsize(self.path) == 0

    def is_file(self) -> bool:
        """Check if path is a file"""
        return os.path.isfile(self.path)

    def get_path(self) -> str:

        return self.path


class VCF(GenomicFile):

    HEADER: dict[str:int] = {
        "CHROM": 0,
        "POS": 1,
        "ID": 2,
        "REF": 3,
        "ALT": 4,
        "QUAL": 5,
        "FILTER": 6,
        "INFO": 7,
        "FORMAT": 8,
        "SAMPLE": 9,
    }

    DEL_FIRST_NC: int = -1

    def __init__(
        self, path: str, caller: callers.VariantCaller, lazy: bool = True
    ):

        super().__init__(path=path)

        self.verify()

        self.caller: callers.VariantCaller = caller

        if not lazy:

            self.parse()

    def get_header(self):

        return self.HEADER

    def parse(self):

        pass

    def verify(self):

        if not self.is_file():

            raise errors.VCFError(f"The file {self.path} does not exist.")

        if self.is_empty():

            raise errors.VCFError(f"The file {self.path} is empty.")

        try:

            with open(self.path, mode="r") as vcf:

                line = vcf.readline()

                if not line:

                    raise errors.VCFError(
                        f"First line of {self.path} is empty."
                    )

                else:

                    # Check if first line start with "#"
                    if line[0] != "#":

                        raise errors.VCFError(
                            f"First line inconsistent with VCF header format"
                        )

        except FileNotFoundError:

            raise errors.VCFError(f"{self.path} is not a valid path")

        except IOError:

            raise errors.VCFError(
                f"An error occurred while reading {self.path}"
            )

    @staticmethod
    def convert(a: object) -> object:
        """Convert variable to appropriate type"""
        try:
            # If the variable contains a '/' or '|' character, it is a genotype information, return the variable as is
            # Else return the variable as an evaluated expression
            return (
                a if sum(list(map(lambda x: x in a, ("/", "|")))) else eval(a)
            )
        except Exception:
            # If the variable cannot be evaluated, return the variable as is
            return a

    def format_to_values(self, values: list[str]) -> dict:
        """map FORMAT string to respective SAMPLE values"""

        # Split the values string into a list of fields
        values: list[str] = values.split(":")

        return {f: VCF.convert(v) for f, v in zip(self.caller.FORMAT, values)}

    def info_to_values(self, values: str) -> dict:

        infos = list(map(lambda item: item.split("="), values.split(";")))

        return {k: v for k, v in infos}

    def VAF(self, variant: str) -> float:

        return self.caller.VAF(variant, self.HEADER)

    def depth(self, variant: str) -> int:

        return self.caller.depth(variant, self.HEADER)

    def arc(self, variant: str) -> tuple[float]:

        return self.caller.arc(variant, self.HEADER)

    def rrc(self, variant: str) -> tuple[float]:

        return self.caller.rrc(variant, self.HEADER)


class Pileup(GenomicFile):

    HEADER: dict[str:int] = {
        "barcode": 0,
        "chromosome": 1,
        "position": 2,
        "reference": 3,
        "depth": 4,
        "A+": 5,
        "A-": 6,
        "T+": 7,
        "T-": 8,
        "C+": 9,
        "C-": 10,
        "G+": 11,
        "G-": 12,
        "N": 13,
        "INS": 14,
        "DEL": 15,
    }

    PLUS_STRAND: list[int] = [0, 2, 4, 6]
    MINUS_STRAND: list[int] = [1, 3, 5, 7]

    DEL_FIRST_NC: int = 1

    def __init__(self, path: str, lazy: bool = True):

        super().__init__(path=path)

        self.verify()

        if not lazy:

            self.parse()

    def get_header(self):

        return self.HEADER

    def parse(self):

        with open(self.path, "r") as pileup:

            indels: dict = {}

            for record in pileup:

                coverage: defaultdict = defaultdict(int)

                chromosome, position, reference, depth, read_base, quality = (
                    record.rstrip("\n").split("\t")
                )

                try:
                    position: int = int(position)
                except ValueError:
                    raise errors.PileupError(f"Incorrect position value {position} in pileup file.")

                indels.setdefault(position, {
                        "insertions": defaultdict(
                            lambda: np.zeros(shape=(1, 2), dtype=np.uint)
                        ),
                        "deletions": defaultdict(
                            lambda: np.zeros(shape=(1, 2), dtype=np.uint)
                        )
                    }
                )

                reference: str = reference.upper()

                read_base: str = re.sub(r"\$", "", re.sub(r"\^.", "", read_base))

                insertions: list[str] = re.findall(
                    r"(\+[0-9]+[ACGTNacgtn]+)", read_base
                )


                for insertion in insertions:

                    size, seq = re.split(r"(?<=\d)(?!.*\d)", insertion)

                    if seq.isupper():

                        indels[position]["insertions"][seq][0][0] += 1

                    else:

                        indels[position]["insertions"][seq.upper()][0][1] += 1

                    read_base: str = re.sub(fr"\{insertion}", "", read_base)

                insstring: str = ";".join(
                    list(
                        map(
                            lambda ins: f"{ins[0]}:{ins[1][0][0]}:{ins[1][0][1]}",
                            indels[position]["insertions"].items(),
                        )
                    )
                ) if len(indels[position]["insertions"].items()) else None

                deletions: list[str] = re.findall(
                    r"(\-[0-9]+[ACGTNacgtn]+)", read_base
                )

                for deletion in deletions:

                    size, seq = re.split(r"(?<=\d)(?!.*\d)", deletion)

                    if seq.isupper():

                        indels[position]["deletions"][seq][0][0] += 1

                    else:

                        indels[position]["deletions"][seq.upper()][0][1] += 1

                    read_base: str = re.sub(fr"\{deletion}", "", read_base)

                for base in read_base:

                    if base == ".":

                        base = reference

                    elif base == ",":

                        base = reference.lower()

                    if base == "A":
                        coverage["A+"] += 1
                    elif base == "a":
                        coverage["A-"] += 1
                    elif base == "T":
                        coverage["T+"] += 1
                    elif base == "t":
                        coverage["T-"] += 1
                    elif base == "C":
                        coverage["C+"] += 1
                    elif base == "c":
                        coverage["C-"] += 1
                    elif base == "G":
                        coverage["G+"] += 1
                    elif base == "g":
                        coverage["G-"] += 1
                    elif base == "N" or base == "n":
                        coverage["N"] += 1
                    elif base == "*":
                        indels[position]["deletions"][base][0][0] += 1
                    else:
                        raise errors.PileupError(
                            f"Unknown base {base} in pileup file"
                        )
                    
                parts = [
                    f"{key}:{values[0][0]}" for key, values in indels[position]["deletions"].items() if key == "*"
                ]
                
                if (position - self.DEL_FIRST_NC) in indels:

                    parts.extend([f"{key}:{values[0][0]}:{values[0][1]}" for key, values in indels[(position - self.DEL_FIRST_NC)]["deletions"].items() if key != "*"])

                    del indels[(position - self.DEL_FIRST_NC)]

                delstring: str = ';'.join(parts) if len(parts) else None

                yield (
                    f"{chromosome}\t{position}\t{reference}\t{depth}\t{coverage['A+']}\t{coverage['A-']}\t{coverage['T+']}\t{coverage['T-']}\t{coverage['C+']}\t{coverage['C-']}\t{coverage['G+']}\t{coverage['G-']}\t{coverage['N']}\t{insstring}\t{delstring}\n"
                )

    def verify(self):

        if not self.is_file():

            raise errors.PileupError(f"The file {self.path} does not exist.")

        if self.is_empty():

            raise errors.PileupError(f"The file {self.path} is empty.")


class VCFIndex(GenomicFile):

    def __init__(self, path: str, lazy: bool = True):
        super().__init__(path)

        self.verify()

    def verify(self):

        pass


class FastaIndex(GenomicFile):

    def __init__(self, path: str, lazy: bool = True):

        super().__init__(path)

        self.verify()

        if not lazy:

            self.parse()

        else:

            self.contigs: pd.DataFrame = pd.DataFrame()

    def get_contigs(self):

        return self.contigs

    def parse(self):

        contigs: list[pd.Series] = []

        with open(self.path, mode="r") as ref:
            for line in ref:
                if line:
                    contigs.append(pd.Series(data=line.strip().split("\t")))

        self.contigs: pd.DataFrame = pd.DataFrame(data=contigs)
        self.contigs.columns = [
            "contig",
            "length",
            "index",
            "pbline",
            "byteline",
        ]
        self.contigs: pd.DataFrame = self.contigs.astype(
            {"contig": "string", "pbline": "uint8", "byteline": "uint8"}
        )

    def verify(self):

        if not self.is_file():

            raise errors.FastaIndexError(
                f"The file {self.path} does not exist."
            )

        if self.is_empty():

            raise errors.FastaIndexError(f"The file {self.path} is empty.")

        try:

            with open(self.path, mode="r") as fasta:

                line = fasta.readline()

                if not line:

                    raise errors.FastaIndexError(
                        f"First line of {self.path} is empty."
                    )

                else:

                    # Check if first line is composed of 5 columns
                    if len(line.split("\t")) != 5:

                        raise errors.FastaIndexError(
                            f"First line inconsistent with Fasta index format"
                        )

        except FileNotFoundError:

            raise errors.FastaIndexError(f"{self.path} is not a valid path")

        except IOError:

            raise errors.FastaIndexError(
                f"An error occurred while reading {self.path}"
            )
