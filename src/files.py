import ast
import autopep8
from hashlib import sha256
import callers
from config import ConfigParser, Term, Expression
from collections import defaultdict
import enum
import exceptions as exceptions
from functools import lru_cache
import jinja2
from loguru import logger
import numpy as np
import os
import pandas as pd
from pathlib import Path
from typing import Final
import re
import utils

try:
    from icecream import ic
except ImportError:  # Graceful fallback if IceCream isn't installed.
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a)  # noqa

class GenomicFile:

    """Base class for genomic files"""

    def __init__(self, path: str):

        # Path to the file
        self.path = path

    def is_empty(self) -> bool:
        """Check if file is empty"""

        return os.path.getsize(self.path) == 0

    def is_file(self) -> bool:
        """Check if path is a file"""

        return os.path.isfile(self.path)

    def get_path(self) -> str:
        """Get the path of the file.
        
        Returns:
            str: The path of the file.
        """

        return self.path
    
    def basename(self) -> str:
        """Get the basename of the file.
        
        Returns:
            str: The basename of the file.
        """

        return os.path.basename(self.path)
    
    def informations(self) -> dict:
        """Get information about the file.
        
        Returns:
            dict: Dictionary containing file information.
        """

        return utils.file_infos(self.path)
    
    def __str__(self):
        """String representation of the file.
        
        Returns:
            str: The path of the file.
        """
        
        return f"{self.path}"
    
    def __repr__(self):
        """String representation of the file for debugging.
        
        Returns:
            str: The path of the file.
        """
        
        return f"{self.path}"
    
    def __hash__(self):
        """Hash of the file.
        
        Returns:
            int: Hash value of the file path.
        """
        
        return hash(self.path)

class VCF(GenomicFile):

    """Class for VCF files"""

    # VCF header
    # 0-based indexed
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
        self, path: str, caller: callers.VariantCaller, lazy: bool = True, index: str = None
    ):

        super().__init__(path=path)

        # Verify if the file is a valid VCF file before anything
        self.verify()

        # The caller from which the VCF file is generated
        self.caller: callers.VariantCaller = caller

        self.index: VCFIndex = index

        if not lazy:

            self.parse()

    @property
    def header(self):
        """Get the header of the VCF file"""

        return getattr(self, "HEADER", None)
    
    def is_compliant(self, record: list[str]):
        """Check if the record is compliant with the VCF format of the caller"""

        return (len(record) >= 10 and self.caller.is_compliant(record, self.HEADER))

    def parse(self):

        with open(self.path, "r") as vcf:

            for record in vcf:

                yield record

    def verify(self):
        """Check if the VCF file is valid"""

        # Check if the file exists
        if not self.is_file():

            raise exceptions.VCFError(f"The file {self.path} does not exist.")

        # Check if the file is empty
        if self.is_empty():

            raise exceptions.VCFError(f"The file {self.path} is empty.")

        try:

            with open(self.path, mode="r") as vcf:

                # Read the first line of the file
                line = vcf.readline()

                # Check if the first line is empty
                if not line:

                    raise exceptions.VCFError(
                        f"First line of {self.path} is empty."
                    )

                else:

                    # Check if first line start with "#"
                    if line[0] != "#":

                        raise exceptions.VCFError(
                            f"First line inconsistent with VCF header format. Expected '#', got {line[0]}."
                        )

        except FileNotFoundError:

            raise exceptions.VCFError(f"{self.path} is not a valid path")

        except IOError:

            raise exceptions.VCFError(
                f"An error occurred while reading {self.path}"
            )
        
        except Exception as e:

            if isinstance(e, exceptions.VCFError):

                raise

            else:

                raise exceptions.VCFError(f"An unexpected error has occurred when validating VCF file: {ic.format(e)}")

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
        """Convert INFO field string to dictionary of key-value pairs.
        
        Args:
            values (str): INFO field string from VCF record.
            
        Returns:
            dict: Dictionary of INFO field key-value pairs.
        """

        infos = list(map(lambda item: item.split("="), values.split(";")))

        return {k: v for k, v in infos}
    
    def genotype(self, variant: list[str]) -> str:
        """Extract the genotype from the variant record"""

        try:
            # Call the genotype method from the caller
            genotype: str = self.caller.genotype(variant, self.HEADER)

            # Check if the genotype is compliant with the VCF format
            # The genotype should contain either '/' or '|' character
            if sum([c in genotype for c in ["/", "|"]]):

                return genotype
            
            else:

                raise exceptions.VCFError("Genotype value does not match the genotype format.")
        
        except (IndexError, NotImplementedError):

            raise exceptions.VCFError("Genotype cannot be extracted.")
        
        except Exception as e:

            if isinstance(e, exceptions.VCFError):

                raise

            else:

                raise exceptions.VCFError(f"An unexpected error has occurred when extracting genotype value: {ic.format(e)}")

    def VAF(self, variant: list[str]) -> float:
        """Extract the variant allele frequency from the variant record"""

        try: 
            # Call the VAF method from the caller
            vaf: float = self.caller.VAF(variant, self.HEADER)

            # Check if the VAF is a valid value
            if vaf <= 0.0:

                raise exceptions.VCFError("Variant allele frequency cannot be zero or negative.")
            
            return vaf
        
        except (IndexError, NotImplementedError):

            raise exceptions.VCFError("Variant allele frequency cannot be extracted from the variant record.")
        
        except Exception as e:

            if isinstance(e, exceptions.VCFError):

                raise

            else:
                raise exceptions.VCFError(f"An unexpected error has occurred when extracting VAF: {ic.format(e)}")

    def depth(self, variant: list[str]) -> int:
        """Extract the depth from the variant record"""

        try:

            # Call the depth method from the caller
            depth: int = self.caller.depth(variant, self.HEADER)

            # Check if the depth is a valid value
            # The depth should be a positive integer
            if depth <= 0:

                raise exceptions.VCFError("Coverage cannot be zero or negative.")
            
            return depth
        
        except (IndexError, NotImplementedError):

            raise exceptions.VCFError("Coverage cannot be extracted.")
        
        except Exception as e:

            if isinstance(e, exceptions.VCFError):

                raise

            else:

                raise exceptions.VCFError(f"An unexpected error has occurred when extracting depth: {ic.format(e)}")

    def arc(self, variant: list[str]) -> tuple[float]:
        """Extract the alternate read count from the variant record"""

        try:
            # Call the arc method from the caller
            arc: tuple[float] = self.caller.arc(variant, self.HEADER)

            # Check if the alternate read count are valid values
            # Alternate read count should be a positive float
            if not all([value > 0 for value in arc if value]):

                raise exceptions.VCFError("Alternate read count cannot be zero or negative.")
            
            return arc

        except (IndexError, NotImplementedError):

            raise exceptions.VCFError("Alternate read count cannot be extracted.")
        
        except Exception as e:

            if isinstance(e, exceptions.VCFError):

                raise

            else:

                raise exceptions.VCFError(f"An unexpected error has occurred when extracting ARC: {ic.format(e)}")
        
    def rrc(self, variant: list[str]) -> tuple[float]:
        """Extract the reference read count from the variant record"""

        try:

            # Call the rrc method from the caller
            rcc: tuple[float] = self.caller.rrc(variant, self.HEADER)

            # Check if the reference read count are valid values
            # Reference read count should be a positive float
            if not all([value > 0 for value in rcc if value]):

                raise exceptions.VCFError("Reference read count cannot be zero or negative.")
            
            return rcc

        except (IndexError, NotImplementedError):

            raise exceptions.VCFError("Reference read count cannot be extracted.")
        
        except Exception as e:

            if isinstance(e, exceptions.VCFError):

                raise

            else:

                raise exceptions.VCFError(f"An unexpected error has occurred when extracting RRC: {ic.format(e)}")

class Pileup(GenomicFile):

    """Class for pileup files"""

    # Header of the formated pileup
    # 0-based indexed
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

    # Indexs of the bases on plus strand
    PLUS_STRAND: Final[list[int]] = [0, 2, 4, 6]
    # Indexs of the bases on minus strand
    MINUS_STRAND: Final[list[int]] = [1, 3, 5, 7]

    DEL_FIRST_NC: Final[int] = 1

    def __init__(self, path: str, sample: str, lazy: bool = True):

        super().__init__(path=path)

        # Verify if the file is a valid Pileup file before anything
        self.verify()

        # Pileup is matched to a sample
        self.sample: str = sample

        if not lazy:

            self.parse()

    @property
    def header(self):
        """Get the header of the Pileup file"""

        return getattr(self, "HEADER", None)

    def parse(self):
        """Parse the Pileup file"""

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
                    raise exceptions.PileupError(f"Incorrect position value {position} in pileup file.")

                if position < 0:

                    raise exceptions.PileupError(f"Incorrect position value {position} in pileup file.")
                    

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

                    raise exceptions.PileupError(f"Incrorrect reference value {reference} in pileup file.")

                try:
                    depth: int = int(depth)
                except ValueError:
                    raise exceptions.PileupError(f"Incorrect depth value {depth} in pileup file.")

                if depth < 0:

                    raise exceptions.PileupError(f"Incorrect depth value {depth} in pileup file.")

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
                        raise exceptions.PileupError(
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
        """Check if the Pileup file is valid"""

        # Check if the file exists
        if not self.is_file():

            raise exceptions.PileupError(f"The file {self.path} does not exist.")
        
        # Check if the file is empty
        if self.is_empty():

            raise exceptions.PileupError(f"The file {self.path} is empty.")
        
        try:

            with open(self.path, mode="r") as pileup:
                
                # Read the first line of the file
                line = pileup.readline()

                # Check if the first line is empty
                if not line:

                    raise exceptions.PileupError(
                        f"First line of {self.path} is empty."
                    )

                else:

                    columns: list[str] = line.split("\t")

                    # Check if first line is composed of 5 or 6 columns
                    # Column 6 is optional in pileup
                    if not len(columns) in [5, 6]:

                        raise exceptions.PileupError(
                            f"First line inconsistent with Pileup format. Expected 5 or 6 columns, got {len(columns)}."
                        )
                    
                    else:

                        try:
                            # Check if position and depth are integers
                            int(columns[1])
                            int(columns[3])

                        except ValueError:

                            raise exceptions.PileupError(f"First line inconsistent with Pileup format. Position and depth values must be integers.")
                        
                        # Check if reference is a single character and in the list of bases
                        if len(columns[2]) != 1 or (not columns[2] in ['A','T','C','G','N']):
                            
                            raise exceptions.PileupError(f"First line inconsistent with Pileup format. Reference value must be a single character in [A,T,C,G,N].")
                        
                        # If the quality column is present
                        if len(columns) == 6:
                            
                            # Check if it is a string of ASCII characters
                            if not (len(columns[5]) and columns[5].isascii()):

                                raise exceptions.PileupError(f"First line inconsistent with Pileup format. Quality value must be a string of ASCII characters.")
                            
        except FileNotFoundError:

            raise exceptions.PileupError(f"{self.path} is not a valid path")

        except IOError:

            raise exceptions.PileupError(
                f"An error occurred while reading {self.path}"
            )
        
        except Exception as e:

            if isinstance(e, exceptions.PileupError):

                raise

            else:

                raise exceptions.PileupError(f"An unexpected error has occurred when validating Pileup file: {ic.format(e)}")

class VCFIndex(GenomicFile):
    """Class for VCF index files"""

    def __init__(self, path: str, lazy: bool = True):
        
        super().__init__(path)

        self.verify()

    def verify(self):

        pass


class FastaIndex(GenomicFile):

    def __init__(self, path: str, lazy: bool = True):

        super().__init__(path)

        # Verify if the file is a valid Fasta index file before anything
        self.verify()

        if not lazy:

            self.parse()

        else:

            self._contigs: pd.DataFrame = pd.DataFrame()

    @property
    def contigs(self):
        """Get the contigs of the Fasta index file"""

        return getattr(self, "_contigs", None)
    
    @contigs.setter
    def contigs(self, value: pd.DataFrame):
        """Set the contigs of the Fasta index file"""

        if not isinstance(value, pd.DataFrame):

            raise TypeError("value must be a Dataframe.")
        
        self._contigs: pd.DataFrame = value

    @lru_cache(maxsize=23)
    def is_indexed_chromosome(self, chromosome: str) -> bool:
        """Check if the chromosome is indexed in the Fasta index file"""

        return chromosome in self._contigs["contig"].values
    
    def is_correct_position(self, chromsome: str, position: int) -> bool:
        """Check if the position is correct for the given chromosome"""

        if not isinstance(position, int):

            return False
        
        try:

            res: bool = position <= (self._contigs.loc[self._contigs["contig"] == chromsome])["length"]

        except (KeyError, pd.errors.IndexingError):

            res: bool = False
        
        return res

    def parse(self):
        """Parse the Fasta index file"""

        # Each contig is a Series
        contigs: list[pd.Series] = []

        with open(self.path, mode="r") as ref:

            for line in ref:

                if line:

                    contigs.append(pd.Series(data=line.strip().split("\t")))

        # Convert the list of Series to a DataFrame
        self._contigs: pd.DataFrame = pd.DataFrame(data=contigs)

        # Rename the columns
        self._contigs.columns = [
            "contig",
            "length",
            "index",
            "pbline",
            "byteline",
        ]
        
        # Convert the columns to the appropriate types, reducing memory usage
        self._contigs: pd.DataFrame = self._contigs.astype(
            {"contig": "string", "length": "uint", "pbline": "uint8", "byteline": "uint8"}
        )

    def verify(self):
        """Check if the Fasta index file is valid"""

        # Check if the file exists
        if not self.is_file():

            raise exceptions.FastaIndexError(
                f"The file {self.path} does not exist."
            )

        # Check if the file is empty
        if self.is_empty():

            raise exceptions.FastaIndexError(f"The file {self.path} is empty.")

        try:

            with open(self.path, mode="r") as fasta:

                # Read the first line of the file
                line = fasta.readline()

                # Check if the first line is empty
                if not line:

                    raise exceptions.FastaIndexError(
                        f"First line of {self.path} is empty."
                    )

                else:

                    columns: list[str] = line.split("\t")

                    # Check if first line is composed of 5 columns
                    if len(columns) != 5:

                        raise exceptions.FastaIndexError(
                            f"First line inconsistent with Fasta index format. Expected 5 columns, got {len(columns)}."
                        )

        except FileNotFoundError:

            raise exceptions.FastaIndexError(f"{self.path} is not a valid path")

        except IOError:

            raise exceptions.FastaIndexError(
                f"An error occurred while reading {self.path}"
            )
        
        except Exception as e:

            if isinstance(e, exceptions.FastaIndexError):

                raise

            else:

                raise exceptions.FastaIndexError(f"An unexpected error has occurred when validating FASTA index file: {ic.format(e)}")

class Config(GenomicFile):

    def __init__(self, path: str, lazy: bool = True):
        """Initialize a Config object.
        
        Args:
            path (str): Path to the config file.
            lazy (bool, optional): Whether to parse the file immediately. Defaults to True.
        """

        super().__init__(path)

        self.verify()

        self.parser = ConfigParser(path)

        if not lazy:

            self.parse()
    
    def verify(self):
        """Verify that the config file exists and is not empty.
        
        Raises:
            ConfigError: If the file does not exist or is empty.
        """

        # Check if the file exists
        if not self.is_file():

            raise exceptions.ConfigError(
                f"Config file {self.path} does not exist."
            )

        # Check if the file is empty
        if self.is_empty():

            raise exceptions.ConfigError(f"The file {self.path} is empty.")
        
    def parse(self):
        """Parse the config file and load its parameters.
        
        Raises:
            SystemExit: If there is an error loading the config file.
        """

        try:

            self.params = self.parser.load()

        except exceptions.ConfigError as e:

            raise SystemExit(e)

class VariantCallerPlugin(GenomicFile):
    """Class for variant caller plugins"""

    # Maximum size of plugin in MB
    MAX_SIZE = 0.01

    STATES = enum.Enum(value="State",
                      names="unsafe safe")
    
    RESSOURCES = os.path.join(
            utils.get_project_dir(), "templates"
        )
    
    TEMPLATE = "caller"

    def __init__(self, id: str, config: Config, path: str = None):

        self.state = VariantCallerPlugin.STATES.unsafe

        self.id = id

        self.config = config

        self.config_hash: str = sha256(str(config.params).encode()).hexdigest()

        if path:

            super().__init__(path)

        else:

            self.package: Path = utils.get_or_create_config_dir()

            self.source: Path = self.package.joinpath(f"{config.params["caller"]["name"]}.py")

            self.sum: Path = self.package.joinpath(f"{config.params["caller"]["name"]}.sum")

            if self.sum.exists():

                pfile, phash = None, None

                cfile, chash  = None, None

                try:
                    # Read the sums
                    with open(self.sum, mode="r") as sumfile:

                        try:

                            pfile, phash = next(sumfile).split('\t')

                            cfile, chash  = next(sumfile).split('\t')

                        except ValueError:

                            logger.warning(f"Checksum file {self.sum} is corrupted.")

                            raise exceptions.CheckSumFileError(f"Checksum file {self.sum} is corrupted.")

                    # Check if the config file is consistent with the cached config file
                    if (chash.strip('\n') == self.config_hash) and (self.source.exists()):

                        hash = utils.hash_file(self.source)

                        # Check if the plugin file is consistent with the cached plugin file
                        if (f"{config.params["caller"]["name"]}.py" == pfile) and (hash == phash.strip('\n')):

                            logger.debug(f"Using cached file {self.source}")

                            super().__init__(self.source)

                            self.verify()

                        else:

                            logger.warning(f"Cached python file checksum is not consistent with {self.source}.")

                            raise exceptions.CheckSumFileError(f"Cached python file checksum is not consistent with {self.source}.")

                    else:

                        logger.warning(f"Cached config file checksum is not consistent with provided config.")

                        raise exceptions.CheckSumFileError(f"Cached config file checksum is not consistent with provided config.")

                except Exception as e:

                    self.remove()

                    if isinstance(e, (exceptions.CheckSumFileError, exceptions.VariantCallerPluginError)):

                        self.write()

                        super().__init__(self.source)

                    else:

                        raise exceptions.VariantCallerPluginError(f"An unexpected error has occurred when loading variant caller plugin: {e}")
            
            else:

                self.write()

                super().__init__(self.source)

    def remove(self):
        """Remove the plugin file and the checksum file."""

        utils.clean(files=[self.source, self.sum])

    def is_safe(self):
        """Check if the plugin file is safe."""
        return self.state == VariantCallerPlugin.STATES.safe

    def _is_valid_python(self, ast: ast.Module) -> bool:
        """
        Securely check if a file is a valid Python code.

        Args:
            ast (ast.Module): The abstract syntax tree of the plugin file.

        Returns:
            bool: True if the plugin file is valid, False otherwise.
        """

        visitor = utils.PluginPythonChecker()

        visitor.visit(ast)

        if len(visitor.not_safe_calls):

            logger.warning(f"Potentially dangerous call {visitor.not_safe_calls} in {self.path}")

            return False
        
        elif visitor.imports != ["abc", "enum", "typing"]:

            logger.warning(f"Importing forbidden modules {visitor.imports} in {self.path}")

            return False
                
        return True

    def verify(self):
        """Verify the plugin file."""

        infs = self.informations()

        fsize = infs["size"]

        if fsize >= self.MAX_SIZE:

            raise exceptions.VariantCallerPluginError("Abnormally high plugin size detected.")
        
        if not self.path.suffix == ".py":

            raise exceptions.VariantCallerPluginError("File name does not end with the Python extension.")
        
        try:

            with open(self.path, mode='r') as plugin:

                content = plugin.read()

            tree: ast.Module = ast.parse(content)

            if not self._is_valid_python(ast=tree):

                raise exceptions.VariantCallerPluginError(f"The plugin file contains code that is either insecure, erroneous or unmanageable in nature.")

        except FileNotFoundError as e:

            raise exceptions.VariantCallerPluginError(f"Cannot found plugin {self.path} on filesystem.")
        
        except SyntaxError as e:

            raise exceptions.VariantCallerPluginError(f"Syntax error in plugin {self.path}.")
        
        except UnicodeDecodeError:

            raise exceptions.VariantCallerPluginError(f"Invalid UTF-8 encoding for plugin {self.path}.")
        
        except Exception as e:

            if isinstance(e, exceptions.VariantCallerPluginError):

                raise

            raise exceptions.VariantCallerPluginError(f"An unexpected error has occurred when reading plugin {self}: {e}")
        # If it come to this instruction, the plugin is safe
        self.state = VariantCallerPlugin.STATES.safe

    def parse(self):
        """Parse the plugin file."""
        pass

    def write(self) -> str:
        """Write the plugin file."""
        def is_expression(object) -> bool:
            """Check if the object is an expression."""
            return isinstance(object, Expression)
        
        def is_term(object) -> bool:
            """Check if the object is a term."""
            return isinstance(object, Term)
        
        def is_pourcentage(object: Term) -> bool:
            """Check if the object is a pourcentage."""
            if object:
                return object.metadata.unit == '%' if object.metadata else False
            else:
                return False
        
        def is_indexed(object: Term) -> bool:
            """Check if the object is indexed."""
            if object:
                return isinstance(object.metadata.index, int) if object.metadata else False
            else:
                return False
        
        def is_in_format(object: Term, format: str) -> bool:
            """Check if the object is in the format field."""
            if object:
                if object.metadata and object.metadata.header:

                    return object.metadata.header == "format"
                
                else:

                    return object.field in format
            else:
                return False
            
        def is_in_infos(object: Term, infos: str) -> bool:
            """Check if the object is in the info field."""
            if object:
                if object.metadata and object.metadata.header:

                    return object.metadata.header == "info"
                
                else:

                    return object.field in infos
            else:
                 return False
        
        logger.debug(f"Rendering python template for {self.config.params["caller"]["name"]}")

        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(self.RESSOURCES)
        )

        env.tests["expression"] = is_expression

        env.tests["term"] = is_term

        env.tests["pourcentage"] = is_pourcentage

        env.tests["indexed"] = is_indexed

        env.tests["in_format"] = is_in_format

        env.tests["in_infos"] = is_in_infos

        # Load the header template
        template = env.get_template(self.TEMPLATE)

        content: str = template.render(
            name=self.config.params["caller"]["name"],
            info=self.config.params["caller"]["info"],
            format=self.config.params["caller"]["format"],
            genotype=self.config.params["caller"]["genotype"],
            depth=self.config.params["caller"]["depth"],
            vaf=self.config.params["caller"]["vaf"],
            rrc=self.config.params["caller"]["rrc"],
            arc=self.config.params["caller"]["arc"],
        )

        fcontent: str = autopep8.fix_code(content)

        with open(self.source, mode="w") as plugin:

            plugin.writelines(
                fcontent
            )

        plugin_hash = sha256(fcontent.encode()).hexdigest()

        with open(self.sum, mode="w") as sumfile:

            sumfile.write(f"{self.config.params["caller"]["name"]}.py\t{plugin_hash}\n")

            sumfile.write(f"{self.config.params["caller"]["name"]}.yaml\t{self.config_hash}\n")

class CheckSumFile:
    """Class for checksum files"""

    def verify(self):
        """Verify the checksum file."""

        pass

    def parse(self):
        """Parse the checksum file."""

        pass

    def compare(self, file):
        """Compare the checksum file with a file.
        
        Args:
            file: The file to compare with the checksum file.
            
        Raises:
            ValueError: If the file is not a valid file.
        """

        pass

class GenomicReader:

    def __init__(self, process: int = 0):
        """Initialize a GenomicReader object.
        
        Args:
            process (int, optional): Number of processes to use. Defaults to 0.
            
        Raises:
            ValueError: If process is negative or not an integer.
        """

        if process < 0:

            raise ValueError("Process parameter cannot be signed integer.")
        
        if not isinstance(process, int):

            raise ValueError("Process parameter must be a unsigned integer.")

        self.process: int = process

    def read(self, file: GenomicFile | list[GenomicFile]):
        """Read genomic files.
        
        Args:
            file (GenomicFile | list[GenomicFile]): File or list of files to read.
            
        Yields:
            The parsed content of the files.
        """

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
        """Initialize a GenomicWritter object.
        
        Args:
            process (int, optional): Number of processes to use. Defaults to 0.
            
        Raises:
            ValueError: If process is negative or not an integer.
        """

        if process < 0:

            raise ValueError("Process parameter cannot be signed integer.")
        
        if not isinstance(process, int):

            raise ValueError("Process parameter must be a unsigned integer.")

        self.process: int = process

    def write(self, output: str, template: str, collection: object | list[object], lookups: set[tuple], sample: str, contigs: object, thresholds: list[float], suffix: str = None):
        """Write genomic data to output files.
        
        Args:
            output (str): Output directory path.
            template (str): Template to use for writing.
            collection (object | list[object]): Data to write.
            lookups (set[tuple]): Lookup data.
            sample (str): Sample name.
            contigs (object): Contig information.
            thresholds (list[float]): Threshold values.
            suffix (str, optional): Output file suffix. Defaults to None.
        """

        def write_pileup():
            """Write pileup data."""

            pass

        def write_vcf(output: str, contigs: object, variants: object, lookups: set[tuple], samples: list[str], thresholds: list[float]):
            """Write VCF data.
            
            Args:
                output (str): Output directory path.
                contigs (object): Contig information.
                variants (object): Variant data.
                lookups (set[tuple]): Lookup data.
                samples (list[str]): List of sample names.
                thresholds (list[float]): Threshold values.
            """

            def format_sample(metrics: dict) -> str:
                """Format sample metrics for output.
                
                Args:
                    metrics (dict): Sample metrics.
                    
                Returns:
                    str: Formatted metrics string.
                """

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