import callers
from collections import defaultdict
import errors
from functools import lru_cache
import jinja2
import numpy as np
import os
import pandas as pd
from typing import Final
import re

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
    
    def __str__(self):
        
        return f"{self.path}"
    
    def __repr__(self):
        
        return f"{self.path}"

class VCF(GenomicFile):

    HEADER: Final[dict[str:int]] = {
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

    DEL_FIRST_NC: Final[int] = -1

    def __init__(
        self, path: str, caller: callers.VariantCaller, lazy: bool = True
    ):

        super().__init__(path=path)

        self.verify()

        self.caller: callers.VariantCaller = caller

        if not lazy:

            self.parse()

    @property
    def header(self):

        return getattr(self, "HEADER", None)
    
    def is_compliant(self, record: list[str]):

        return (len(record) >= 10 and self.caller.is_compliant(record, self.HEADER))

    def parse(self):

        with open(self.path, "r") as vcf:

            for record in vcf:

                yield record

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
    
    def genotype(self, variant: list[str]) -> str:

        try:

            genotype: str = self.caller.genotype(variant, self.HEADER)
        
            if sum([c in genotype for c in ["/", "|"]]):

                return genotype
            
            else:

                raise errors.VCFError("Genotype value does not match the genotype format.")
        
        except (IndexError, NotImplementedError):

            raise errors.VCFError("Genotype cannot be extracted.")

    def VAF(self, variant: list[str]) -> float:

        try: 

            vaf: float = self.caller.VAF(variant, self.HEADER)

            if vaf <= 0.0:

                raise errors.VCFError("Variant allele frequency cannot be zero or negative.")
            
            return vaf
        
        except (IndexError, NotImplementedError):

            raise errors.VCFError("Variant allele frequency cannot be extracted.")

    def depth(self, variant: list[str]) -> int:

        try:

            depth: int = self.caller.depth(variant, self.HEADER)

            if depth <= 0:

                raise errors.VCFError("Coverage cannot be zero or negative.")
            
            return depth
        
        except (IndexError, NotImplementedError):

            raise errors.VCFError("Coverage cannot be extracted.")

    def arc(self, variant: list[str]) -> tuple[float]:

        try:

            arc: tuple[float] = self.caller.arc(variant, self.HEADER)

            if not all([value > 0 for value in arc if value]):

                raise errors.VCFError("Alternate read count cannot be zero or negative.")
            
            return arc

        except (IndexError, NotImplementedError):

            raise errors.VCFError("Alternate read count cannot be extracted.")
        
    def rrc(self, variant: list[str]) -> tuple[float]:

        try:

            rcc: tuple[float] = self.caller.rrc(variant, self.HEADER)
        
            if not all([value > 0 for value in rcc if value]):

                raise errors.VCFError("Reference read count cannot be zero or negative.")
            
            return rcc

        except (IndexError, NotImplementedError):

            raise errors.VCFError("Reference read count cannot be extracted.")

class Pileup(GenomicFile):

    HEADER: Final[dict[str:int]] = {
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

    PLUS_STRAND: Final[list[int]] = [0, 2, 4, 6]
    MINUS_STRAND: Final[list[int]] = [1, 3, 5, 7]

    DEL_FIRST_NC: Final[int] = 1

    def __init__(self, path: str, sample: str, lazy: bool = True):

        super().__init__(path=path)

        #self.verify()

        self.sample: str = sample

        if not lazy:

            self.parse()

    @property
    def header(self):

        return getattr(self, "HEADER", None)

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

                if position < 0:

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

                if not reference.isalpha():

                    raise errors.PileupError(f"Incrorrect reference value {reference} in pileup file.")

                try:
                    depth: int = int(depth)
                except ValueError:
                    raise errors.PileupError(f"Incorrect depth value {depth} in pileup file.")

                if depth < 0:

                    raise errors.PileupError(f"Incorrect depth value {depth} in pileup file.")

                regex: dict[str:str] = {r'\^.': '',
                                        r'\$': '',
                                        '[BD-EH-SU-Z]': 'N',
                                        '[bd-eh-su-z]': 'n'}
                
                for pattern in regex:

                    read_base: str = re.sub(pattern, regex[pattern], read_base)

                insertions: list[str] = re.findall(
                    r"(\+[0-9]+[ATCGNatcgn]+)", read_base
                )

                for insertion in insertions:

                    depth += 1

                    size, seq = re.split(r"(?<=\d)(?!.*\d)", insertion)

                    seq: str = seq[:int(size[1:])]

                    if seq.isupper():

                        indels[position]["insertions"][seq][0][0] += 1

                    else:

                        indels[position]["insertions"][seq.upper()][0][1] += 1

                    read_base: str = re.sub(fr"\{size}{seq}", "", read_base)

                inss: str = ";".join(
                    list(
                        map(
                            lambda ins: f"{ins[0]}:{ins[1][0][0]},{ins[1][0][1]}",
                            indels[position]["insertions"].items(),
                        )
                    )
                ) if len(indels[position]["insertions"].items()) else None

                deletions: list[str] = re.findall(
                    r"(\-[0-9]+[ATCGNatcgn]+)", read_base
                )

                for deletion in deletions:

                    depth += 1

                    size, seq = re.split(r"(?<=\d)(?!.*\d)", deletion)

                    seq: str = seq[:int(size[1:])]

                    if seq.isupper():

                        indels[position]["deletions"][seq][0][0] += 1

                    else:

                        indels[position]["deletions"][seq.upper()][0][1] += 1

                    read_base: str = re.sub(fr"\{size}{seq}", "", read_base)

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

                    parts.extend([f"{key}:{values[0][0]},{values[0][1]}" for key, values in indels[(position - self.DEL_FIRST_NC)]["deletions"].items() if key != "*"])

                delss: str = ';'.join(parts) if len(parts) else None

                if len(indels) == 2:

                    del indels[list(indels.keys())[0]]

                yield (
                    f"{self.sample}\t{chromosome}\t{position}\t{reference}\t{depth}\t{coverage['A+']}\t{coverage['A-']}\t{coverage['T+']}\t{coverage['T-']}\t{coverage['C+']}\t{coverage['C-']}\t{coverage['G+']}\t{coverage['G-']}\t{coverage['N']}\t{inss}\t{delss}\n"
                )

    def verify(self):

        if not self.is_file():

            raise errors.PileupError(f"The file {self.path} does not exist.")

        if self.is_empty():

            raise errors.PileupError(f"The file {self.path} is empty.")
        
        try:

            with open(self.path, mode="r") as pileup:

                line = pileup.readline()

                if not line:

                    raise errors.PileupError(
                        f"First line of {self.path} is empty."
                    )

                else:
                    columns: list[str] = len(line.split("\t"))
                    # Check if first line is composed of 5 or 6 columns
                    # Column 6 is optional in pileup
                    if not len(columns) in [5, 6]:

                        raise errors.PileupError(
                            f"First line inconsistent with Pileup format"
                        )
                    
                    else:

                        try:

                            int(columns[1])
                            int(columns[3])

                        except ValueError:

                            raise errors.PileupError(f"First line inconsistent with Pileup format")
                        
                        if len(columns[2]) != 1 or (not columns[2] in ['A','T','C','G']):
                            
                            raise errors.PileupError(f"First line inconsistent with Pileup format")
                        
                        if len(columns) == 6:
                        
                            if not (columns[5].isascii() and len(columns[5])):

                                raise errors.PileupError(f"First line inconsistent with Pileup format")
                            
        except FileNotFoundError:

            raise errors.PileupError(f"{self.path} is not a valid path")

        except IOError:

            raise errors.PileupError(
                f"An error occurred while reading {self.path}"
            )

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

            self._contigs: pd.DataFrame = pd.DataFrame()

    @property
    def contigs(self):

        return getattr(self, "_contigs", None)
    
    @contigs.setter
    def contigs(self, value: pd.DataFrame):

        if not isinstance(value, pd.DataFrame):

            raise TypeError("value must be a Dataframe.")
        
        self._contigs: pd.DataFrame = value

    @lru_cache(maxsize=23)
    def is_indexed_chromosome(self, chromosome: str) -> bool:

        return chromosome in self._contigs["contig"].values
    
    def is_correct_position(self, chromsome: str, position: int) -> bool:

        if not isinstance(position, int):

            return False
        
        try:

            res: bool = position <= (self._contigs.loc[self._contigs["contig"] == chromsome])["length"]

        except (KeyError, pd.errors.IndexingError):

            res: bool = False
        
        return res

    def parse(self):

        contigs: list[pd.Series] = []

        with open(self.path, mode="r") as ref:

            for line in ref:

                if line:

                    contigs.append(pd.Series(data=line.strip().split("\t")))

        self._contigs: pd.DataFrame = pd.DataFrame(data=contigs)

        self._contigs.columns = [
            "contig",
            "length",
            "index",
            "pbline",
            "byteline",
        ]
        
        self._contigs: pd.DataFrame = self._contigs.astype(
            {"contig": "string", "length": "uint", "pbline": "uint8", "byteline": "uint8"}
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

class GenomicReader:

    def __init__(self, process: int = 0):

        if process < 0:

            raise ValueError("Process parameter cannot be signed integer.")
        
        if not isinstance(process, int):

            raise ValueError("Process parameter must be a unsigned integer.")

        self.process: int = process

    def read(self, file: GenomicFile | list[GenomicFile]):

        if isinstance(file, list):
            
            # Parallel computing
            if self.process:

                pass
            # Sequential computing
            else:

                pass

        else:

            yield from file.parse()


class GenomicWritter:

    def __init__(self, process: int = 0):

        if process < 0:

            raise ValueError("Process parameter cannot be signed integer.")
        
        if not isinstance(process, int):

            raise ValueError("Process parameter must be a unsigned integer.")

        self.process: int = process

    def write(self, output: str, template: str, collection: object | list[object], lookups: set[tuple], sample: str, contigs: object, thresholds: list[float], suffix: str = None):

        def write_pileup():

            pass

        def write_vcf(output: str, contigs: object, variants: object, lookups: set[tuple], samples: list[str], thresholds: list[float]):

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
            ressources: str = os.path.join(
                os.path.dirname(os.path.abspath(__file__)), "templates"
            )

            env = jinja2.Environment(loader=jinja2.FileSystemLoader(ressources))

            # Load the header template
            template = env.get_template("header")

            with open(output, mode="w") as out:

                # Write the rendered header to the file
                out.writelines(
                    template.render(contigs=contigs, thresholds=thresholds)
                )

                header: str = "\t".join(HEADER)

                # Write the column names
                out.write(f"\n#{header}\n")

                # Iterate over each level of the variants dictionary
                for lookup in lookups:

                        variant: dict = variants[lookup[0]][lookup[1]][
                            lookup[2]
                        ]

                        ref, alt = lookup[2].split(":")

                        # Write ONLY if normalized metrics are present
                        if "sample" in variant:

                            out.write(
                                "\t".join(
                                    [
                                        f"chr{lookup[0]}",  # Chromosome field
                                        str(
                                            lookup[1].vcf_position
                                        ),  # Position field
                                        ".",  # ID field
                                        ref,  # Reference field
                                        alt,  # Alternate field
                                        ".", # Qual field
                                        variant["filter"],  # Filter field
                                        "=".join(
                                            [INFOS[0], variant["type"]]
                                        ),  # Info field
                                        ":".join(FORMAT),  # Format field
                                        format_sample(variant["sample"]), # Sample values field
                                    ]
                                )
                            ) 

                            out.write("\n")

        if isinstance(collection, list):

            # Parallel computing
            if self.process:

                pass

            # Sequential computing
            else:

                pass

        else:

            if template == "vcf":

                output: str = os.path.join(output, f"{sample}.{suffix}.vcf" if suffix else f"{sample}.vcf") 

                write_vcf(output=output, 
                         contigs=contigs, 
                         variants=collection,
                         lookups=lookups, 
                         samples=[sample], 
                         thresholds=thresholds)