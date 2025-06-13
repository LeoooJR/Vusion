from abc import ABC, abstractmethod
import autopep8
from config import Term, Expression
import enum
import exceptions as exceptions
from hashlib import sha256
import importlib
import jinja2
from loguru import logger
import os
from pathlib import Path
import sys
from typing import final
import utils

try:
    from icecream import ic
except ImportError:  # Graceful fallback if IceCream isn't installed.
    ic = lambda *a: None if not a else (a[0] if len(a) == 1 else a)  # noqa

class VariantCaller(ABC):

    """ Abstract class for variant callers """
    
    @staticmethod
    @abstractmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:
        """Extract the genotype from the variant."""
        pass

    @staticmethod
    @abstractmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:
        """Calculate the variant allele frequency (VAF)."""
        pass

    @staticmethod
    @abstractmethod
    def depth(variant: list[str], header: dict[int]) -> int:
        """Extract the depth of the variant."""
        pass

    @staticmethod
    @abstractmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """Extract the reference allele counts."""
        pass

    @staticmethod
    @abstractmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """Extract the alternate allele counts."""
        pass

class VariantCallerRepository:

    """ A class to manage the variant callers """

    __slots__ = ["callers"]

    # known variant callers are:
    # deepvariant       (DV)
    # bcftools			(ST)
    # varscan			(VS)
    # vardict			(VD)
    # pindel			(PL)
    # haplotypecaller	(HC)

    def __init__(self):

        # Initialize the supported variant callers
        self.callers: dict[str: VariantCaller] = {
            "BT": BCFTools(),
            "VS": Varscan(),
            "VD": Vardict(),
            "PL": Pindel(),
            "HC": Haplotypecaller(),
            "FL": Filt3r(),
            "DV": DeepVariant(),
        }
    
    # Get the variant caller by its name
    # If the caller is not supported, raise an error
    def get_VC(self, caller: str) -> VariantCaller:
        """Get a variant caller by its name.
        
        Args:
            caller (str): The name of the variant caller.
            
        Returns:
            VariantCaller: The variant caller instance.
            
        Raises:
            VariantCallerError: If the caller is not supported.
        """

        try:

            return self.callers[caller.upper()]

        except KeyError:

            raise exceptions.VariantCallerError(
                f"Variant Caller {caller} not supported."
            )

    def get_BT(self) -> VariantCaller:
        """Get the BCFTools variant caller.
        
        Returns:
            VariantCaller: The BCFTools variant caller instance.
        """

        return self.callers["BT"]

    def get_VS(self) -> VariantCaller:
        """Get the Varscan variant caller.
        
        Returns:
            VariantCaller: The Varscan variant caller instance.
        """

        return self.callers["VS"]

    def get_VD(self) -> VariantCaller:
        """Get the Vardict variant caller.
        
        Returns:
            VariantCaller: The Vardict variant caller instance.
        """

        return self.callers["VD"]

    def get_PL(self) -> VariantCaller:
        """Get the Pindel variant caller.
        
        Returns:
            VariantCaller: The Pindel variant caller instance.
        """

        return self.callers["PL"]

    def get_HS(self) -> VariantCaller:
        """Get the Haplotypecaller variant caller.
        
        Returns:
            VariantCaller: The Haplotypecaller variant caller instance.
        """

        return self.callers["HS"]

    def get_FL(self) -> VariantCaller:
        """Get the Filt3r variant caller.
        
        Returns:
            VariantCaller: The Filt3r variant caller instance.
        """

        return self.callers["FL"]

    def get_DV(self) -> VariantCaller:
        """Get the DeepVariant caller.
        
        Returns:
            VariantCaller: The DeepVariant caller instance.
        """

        return self.callers["DV"]

    def is_supported(self, caller: str) -> bool:
        """Check if a variant caller is supported.
        
        Args:
            caller (str): The name of the variant caller to check.
            
        Returns:
            bool: True if the caller is supported, False otherwise.
        """

        return caller in self.callers
    
    def write(self, recipe, dest) -> str:

        def is_expression(object) -> bool:

            return isinstance(object, Expression)
        
        def is_term(object) -> bool:

            return isinstance(object, Term)
        
        def is_pourcentage(object: Term) -> bool:

            return object.metadata.unit == '%' if object.metadata else False
        
        def is_indexed(object: Term):

            return isinstance(object.metadata.index, int) if object.metadata else False
        
        def is_in_format(object: Term, format: str):

            if object.metadata and object.metadata.header:

                return object.metadata.header == "format"
            
            else:

                return object.field in format
            
        def is_in_infos(object: Term, infos: str):

            if object.metadata and object.metadata.header:

                return object.metadata.header == "info"
            
            else:

                return object.field in infos            
        
        logger.debug(f"Rendering python template for {recipe["caller"]["name"]}")

        # Path to the template directory
        ressources = os.path.join(
            utils.get_project_dir(), "templates"
        )

        env = jinja2.Environment(
            loader=jinja2.FileSystemLoader(ressources)
        )

        env.tests["expression"] = is_expression

        env.tests["term"] = is_term

        env.tests["pourcentage"] = is_pourcentage

        env.tests["indexed"] = is_indexed

        env.tests["in_format"] = is_in_format

        env.tests["in_infos"] = is_in_infos

        # Load the header template
        template = env.get_template("caller")

        content: str = template.render(
            name=recipe["caller"]["name"],
            info=recipe["caller"]["info"],
            format=recipe["caller"]["format"],
            genotype=recipe["caller"]["genotype"],
            depth=recipe["caller"]["depth"],
            vaf=recipe["caller"]["vaf"],
            rrc=recipe["caller"]["rrc"],
            arc=recipe["caller"]["arc"],
        )

        with open(dest, mode="w") as plugin:

            plugin.writelines(
                autopep8.fix_code(content)
            )

        plugin_hash = sha256(content.encode()).hexdigest()

        return plugin_hash
    
    def load(self, id, package, source):

        try:

            if not package in sys.path:

                sys.path.append(str(package))

            module = importlib.import_module(source, package=str(package))

            self.callers[id] = getattr(module, source)()

        except ImportError as e:

            logger.error(e)

            raise exceptions.VariantCallerPluginError(e)

        except AttributeError as e:

            logger.error(e)

            raise exceptions.VariantCallerPluginError(e)
        
        except Exception as e:

            if isinstance(e, ImportError) or isinstance(e, AttributeError):

                raise

            else:

                raise exceptions.VariantCallerPluginError(f"An unexpected error has occurred when loading variant caller plugin: {e}")
            
    def clean(self, files: list[Path]):

        for file in files:

            file.unlink()
    
    def add(self, id: str, recipe: dict):
        
        config_hash = sha256(str(recipe).encode()).hexdigest()

        package: Path = utils.get_or_create_config_dir()

        sum: Path = package.joinpath(f"{recipe["caller"]["name"]}.sum")

        source: Path = package.joinpath(f"{recipe["caller"]["name"]}.py")

        if sum.exists():

            pfile, phash = None, None

            cfile, chash  = None, None

            try:

                with open(sum, mode="r") as sumfile:

                    try:

                        pfile, phash = next(sumfile).split('\t')

                        cfile, chash  = next(sumfile).split('\t')

                    except ValueError:

                        logger.warning(f"Checksum file {sum} is corrupted.")

                        raise exceptions.CheckSumFileError(f"Checksum file {sum} is corrupted.")

                if (chash.strip('\n') == config_hash) and (source.exists()):

                    hash = utils.hash_file(source)

                    if (f"{recipe["caller"]["name"]}.py" == pfile) and (hash == phash.strip('\n')):

                        logger.debug(f"Using cached file {source}")

                    else:

                        logger.warning(f"Cached python file checksum is not consistent with {source}.")

                        raise exceptions.CheckSumFileError(f"Cached python file checksum is not consistent with {source}.")

                else:

                    logger.warning(f"Cached config file checksum is not consistent with provided config.")

                    raise exceptions.CheckSumFileError(f"Cached config file checksum is not consistent with provided config.")

            except Exception as e:

                if isinstance(e, exceptions.CheckSumFileError):

                    plugin_hash: str = self.write(recipe, source)

                    with open(sum, mode="w") as sumfile:

                        sumfile.write(f"{recipe["caller"]["name"]}.py\t{plugin_hash}\n")

                        sumfile.write(f"{recipe["caller"]["name"]}.yaml\t{config_hash}\n")

                else:

                    self.clean(files=[source, sum])

                    raise exceptions.VariantCallerPluginError(f"An unexpected error has occurred when loading variant caller plugin: {ic.format(e)}")
        
        else:

            plugin_hash: str = self.write(recipe, source)

            with open(sum, mode="w") as sumfile:

                sumfile.write(f"{recipe["caller"]["name"]}.py\t{plugin_hash}\n")

                sumfile.write(f"{recipe["caller"]["name"]}.yaml\t{config_hash}\n")

        try:
                            
            self.load(id=id, package=package, source=recipe["caller"]["name"])

        except Exception:
            
            self.clean(files=[source, sum])

            raise
    
    def __len__(self):

        return len(self.callers)

@final
class BCFTools(VariantCaller):

    """ A class to manage the BCFtools variant caller """

    # BCFtools FORMAT
    FORMAT = enum.IntEnum(value="FORMAT",
                          names=','.join(["GT", "PL"]),
                          start=0)

    # Check if the variant is compliant with the BCFtools format.
    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """Check if the variant is compliant with the BCFtools format.
        
        Args:
            variant (list[str]): The variant to check.
            header (dict[str:int]): The header of the variant.
            
        Returns:
            bool: True if the variant is compliant, False otherwise.
        """
        
        return (
            len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == [key.name for key in self.FORMAT]            
            ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
        ) and (
            len(variant[header["INFO"]]) and any(
            list(
                map(
                    lambda x: f"{x}=" in variant[header["INFO"]], ["DP", "DP4"]
                )
            )
        ))
    
    # Extract the genotype from the variant
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:

        return variant[header["SAMPLE"]].split(':')[BCFTools.FORMAT.GT.value]

    # Calculate the variant allele frequency (VAF)
    # The VAF is the depth of the variant allele divided by the total depth
    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:

        total_depth: int = int(
            variant[header["INFO"]].split("DP=")[1].split(";")[0]
        )

        alleles_depth: int = (
            variant[header["INFO"]].split("DP4=")[1].split(";")[0]
        )

        variant_depth: int = float(alleles_depth.split(",")[2]) + float(
            alleles_depth.split(",")[3]
        )

        try:
            vaf: float = variant_depth / total_depth
        except ZeroDivisionError:
            vaf: float = 0.0

        return vaf

    # Extract the depth of the variant
    @staticmethod
    def depth(variant: list[str], header: dict[int]) -> int:

        return int(variant[header["INFO"]].split("DP=")[1].split(";")[0])

    # Extract the reference allele counts
    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:

        alleles_depth: list[str] = (
            variant[header["INFO"]].split("DP4=")[1].split(";")[0]
        )
        depths: list[str] = alleles_depth.split(",")
        rrc_plus: int = int(depths[0])
        rrc_minus: int = int(depths[1])
        rrc: int = rrc_plus + rrc_minus
        return (rrc, rrc_plus, rrc_minus)

    # Extract the alternate allele counts
    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:

        alleles_depth: list[str] = (
            variant[header["INFO"]].split("DP4=")[1].split(";")[0]
        )
        depths: list[str] = alleles_depth.split(",")
        arc_plus: int = int(depths[2])
        arc_minus: int = int(depths[3])
        arc: int = arc_plus + arc_minus
        return (arc, arc_plus, arc_minus)

    def __str__(self):
        return "BCFTools"
    
    def __repr__(self):
        return "BCFTools"

@final
class Varscan(VariantCaller):

    """ A class to manage the Varscan variant caller """

    # Varscan FORMAT
    FORMAT = enum.IntEnum(value="FORMAT",
                          names=','.join(["GT",
                          "GQ",
                          "SDP",
                          "DP",
                          "RD",
                          "AD",
                          "FREQ",
                          "PVAL",
                          "RBQ",
                          "ABQ",
                          "RDF",
                          "RDR",
                          "ADF",
                          "ADR",
                          ]),
                          start=0)

    # Check if the variant is compliant with the Varscan format.
    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """Check if the variant is compliant with the Varscan format.
        
        Args:
            variant (list[str]): The variant to check.
            header (dict[str:int]): The header of the variant.
            
        Returns:
            bool: True if the variant is compliant, False otherwise.
        """

        return len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == [key.name for key in self.FORMAT]            
        ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
    
    # Extract the genotype from the variant
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:

        return variant[header["SAMPLE"]].split(':')[Varscan.FORMAT.GT.value]

    # Calculate the variant allele frequency (VAF)
    # The VAF is the depth of the variant allele divided by the total depth
    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:

        vaf = (
            float((variant[header["SAMPLE"]].split(":")[Varscan.FORMAT.FREQ.value]).rstrip('%')) / 100
        )

        return vaf

    # Extract the depth of the variant
    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:

        return int(variant[header["SAMPLE"]].split(":")[Varscan.FORMAT.SDP.value])

    # Extract the reference allele counts
    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:

        values = variant[header["SAMPLE"]].split(":")

        rrc = int(values[Varscan.FORMAT.RD.value])
        rrc_plus = int(values[Varscan.FORMAT.RDF.value])
        rrc_minus = int(values[Varscan.FORMAT.RDR.value])

        return (rrc, rrc_plus, rrc_minus)

    # Extract the alternate allele counts
    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:

        values = variant[header["SAMPLE"]].split(":")

        arc = int(values[Varscan.FORMAT.AD.value])
        arc_plus = int(values[Varscan.FORMAT.ADF.value])
        arc_minus = int(values[Varscan.FORMAT.ADR.value])

        return (arc, arc_plus, arc_minus)

    def __str__(self):
        return "Varscan"
    
    def __repr__(self):
        return "Varscan"

@final
class Vardict(VariantCaller):

    """ A class to manage the Vardict variant caller """

    # Vardict FORMAT
    FORMAT = enum.IntEnum(value="FORMAT",
                          names=','.join(["GT", "DP", "VD", "AD", "AF", "RD", "ALD"]),
                          start=0)

    # Check if the variant is compliant with the Vardict format.
    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """Check if the variant is compliant with the Vardict format.
        
        Args:
            variant (list[str]): The variant to check.
            header (dict[str:int]): The header of the variant.
            
        Returns:
            bool: True if the variant is compliant, False otherwise.
        """
        
        return (
            len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == [key.name for key in self.FORMAT]           
            ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
        ) and (
            (len(variant[header["INFO"]]))
            and ("AF=" in variant[header["INFO"]])
        )
    
    # Extract the genotype from the variant
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:

        gt: str = variant[header["SAMPLE"]].split(':')[Vardict.FORMAT.GT.value]
        # Manage case when Vardict return 1/0 instead of 0/1
        return "0/1" if gt == "1/0" else gt

    # Calculate the variant allele frequency (VAF)
    # The VAF is the depth of the variant allele divided by the total depth
    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:

        vaf = float(variant[header["INFO"]].split("AF=")[1].split(";")[0])

        return vaf

    # Extract the depth of the variant
    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:

        return int(variant[header["SAMPLE"]].split(":")[Vardict.FORMAT.DP.value])

    # Extract the reference allele counts
    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:

        values: list[str] = variant[header["SAMPLE"]].split(":")

        rrc = int(values[3].split(",")[0])

        rd: list[str] = values[Vardict.FORMAT.RD.value].split(",")
        rrc_plus: int = int(rd[0])
        rrc_minus: int = int(rd[1])

        return (rrc, rrc_plus, rrc_minus)

    # Extract the alternate allele counts
    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:

        values: list[str] = variant[header["SAMPLE"]].split(":")

        arc = int(values[Vardict.FORMAT.AD.value].split(",")[1])

        ald: list[int] = values[Vardict.FORMAT.ALD.value].split(",")
        arc_plus: int = int(ald[0])
        arc_minus: int = int(ald[1])

        return (arc, arc_plus, arc_minus)

    def __str__(self):
        return "Vardict"
    
    def __repr__(self):
        return "Vardict"

@final
class Pindel(VariantCaller):
    
    """ A class to manage the Pindel variant caller """

    # Pindel FORMAT
    FORMAT = enum.IntEnum(value="FORMAT",
                          names=','.join(["GT", "AD"]),
                          start=0)

    # Check if the variant is compliant with the Pindel format.
    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """Check if the variant is compliant with the Pindel format.
        
        Args:
            variant (list[str]): The variant to check.
            header (dict[str:int]): The header of the variant.
            
        Returns:
            bool: True if the variant is compliant, False otherwise.
        """
        
        return len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == [key.name for key in self.FORMAT]    
        ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
    
    # Extract the genotype from the variant
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:

        return variant[header["SAMPLE"]].split(':')[Pindel.FORMAT.GT.value]

    # Calculate the variant allele frequency (VAF)
    # The VAF is the depth of the variant allele divided by the total depth
    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:

        depths = variant[header["SAMPLE"]].split(":")[Pindel.FORMAT.AD.value].split(",")
        try:
            vaf: float = float(depths[1]) / (float(depths[0]) + float(depths[1]))
        except ZeroDivisionError:
            vaf: float = 0.0
            
        return vaf

    # Extract the depth of the variant
    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:

        alleles_depth: list[str] = (
            variant[header["SAMPLE"]].split(":")[Pindel.FORMAT.AD.value].split(",")
        )

        return int(alleles_depth[0]) + int(alleles_depth[1])

    # Extract the reference allele counts
    # Pindel does not provide the reference allele counts
    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:

        return (None, None, None)

    # Extract the alternate allele counts
    # Pindel does not provide the alternate allele counts for each strand
    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:

        arc = int(variant[header["SAMPLE"]].split(":")[Pindel.FORMAT.AD.value].split(",")[1])

        return (arc, None, None)

    def __str__(self):

        return "Pindel"
    
    def __repr__(self):

        return "Pindel"

@final
class Haplotypecaller(VariantCaller):
    
    """ A class to manage the Haplotypecaller variant caller """

    # Haplotypecaller FORMAT
    FORMAT = enum.IntEnum(value="FORMAT",
                          names=','.join(["GT", "AD", "DP", "GQ", "PL"]),
                          start=0)

    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """Check if the variant is compliant with the Haplotypecaller format.
        
        Args:
            variant (list[str]): The variant to check.
            header (dict[str:int]): The header of the variant.
            
        Returns:
            bool: True if the variant is compliant, False otherwise.
        """
        
        return len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == [key.name for key in self.FORMAT]            
        ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
    
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:
        """Extract the genotype from the variant."""

        return variant[header["SAMPLE"]].split(':')[Haplotypecaller.FORMAT.GT.value]

    # The VAF is the depth of the variant allele divided by the total depth
    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:
        """Calculate the variant allele frequency (VAF)."""

        metrics: list[str] = variant[header["SAMPLE"]].split(":")

        total_depth = float(metrics[Haplotypecaller.FORMAT.DP.value])
        alleles_depth = float(metrics[1].split(",")[Haplotypecaller.FORMAT.AD.value])

        try:
            vaf: float = alleles_depth / total_depth
        except ZeroDivisionError:
            vaf: float = 0.0

        return vaf

    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:
        """Extract the depth of the variant."""

        return int(variant[header["SAMPLE"]].split(":")[Haplotypecaller.FORMAT.DP.value])

    # Haplotypecaller does not provide the reference allele counts for each strand
    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """Extract the reference allele counts."""

        return (None, None, None)

    # Haplotypecaller does not provide the alternate allele counts for each strand
    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """Extract the alternate allele counts."""

        arc = int(variant[header["SAMPLE"]].split(":")[Haplotypecaller.FORMAT.AD.value].split(",")[1])

        return (arc, None, None)

    def __str__(self):

        return "Haplotypecaller"
    
    def __repr__(self):

        return "Haplotypecaller"

@final
class Filt3r(VariantCaller):
    
    """ A class to manage the Filt3r variant caller """

    # Filt3r FORMAT
    FORMAT = []

    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:
        """Extract the genotype from the variant."""

        raise NotImplementedError

    # The VAF is the depth of the variant allele divided by the total depth
    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:
        """Calculate the variant allele frequency (VAF)."""

        return float(variant[6].split(";")[4].split("=")[1])

    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:
        """Extract the depth of the variant."""

        infos: list[str] = variant[6].split(";")

        variant_reads: int = int(infos[1].split("=")[1])

        reference_reads: int = int(infos[2].split("=")[1])

        return variant_reads + reference_reads

    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """Extract the reference allele counts."""

        return (None, None, None)

    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """Extract the alternate allele counts."""

        arc = int(variant[6].split(";")[1].split("=")[1])

        return (arc, None, None)

    def __str__(self):

        return "Flit3r"
    
    def __repr__(self):

        return "Flit3r"

@final
class DeepVariant(VariantCaller):
    
    """ A class to manage the DeepVariant variant caller """

    # DeepVariant FORMAT
    FORMAT = enum.IntEnum(value="FORMAT",
                          names=','.join(["GT", "GQ", "DP", "AD", "VAF", "PL"]),
                          start=0)

    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """Check if the variant is compliant with the DeepVariant format.
        
        Args:
            variant (list[str]): The variant to check.
            header (dict[str:int]): The header of the variant.
            
        Returns:
            bool: True if the variant is compliant, False otherwise.
        """
        
        return len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == [key.name for key in self.FORMAT]            
        ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
    
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:
        """Extract the genotype from the variant."""

        return variant[header["SAMPLE"]].split(':')[DeepVariant.FORMAT.GT.value]

    # The VAF is the depth of the variant allele divided by the total depth
    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:
        """Calculate the variant allele frequency (VAF)."""

        return float(variant[header["SAMPLE"]].split(":")[DeepVariant.FORMAT.VAF.value])

    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:
        """Extract the depth of the variant."""

        return float(variant[header["SAMPLE"]].split(":")[DeepVariant.FORMAT.DP.value])

    # DeepVariant does not provide the reference allele counts for each strand
    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """Extract the reference allele counts."""

        return (
            int(variant[header["SAMPLE"]].split(":")[DeepVariant.FORMAT.AD.value].split(",")[0]),
            None,
            None,
        )

    # DeepVariant does not provide the alternate allele counts for each strand
    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """Extract the alternate allele counts."""

        arc = int(variant[header["SAMPLE"]].split(":")[DeepVariant.FORMAT.AD.value].split(",")[1])

        return (arc, None, None)

    def __str__(self):
        return "Deepvariant"

    def __repr__(self):
        return "Deepvariant"
