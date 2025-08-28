from abc import ABC, abstractmethod
import enum
import exceptions as exceptions
import importlib
from loguru import logger
from repository import Repository
import sys
from typing import final

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

class VariantCallerRepository(Repository):

    """ A class to manage the variant callers """

    __slots__ = ["callers", "plugins"]

    # known variant callers are:
    # deepvariant       (DV)
    # bcftools			(ST)
    # varscan			(VS)
    # vardict			(VD)
    # pindel			(PL)
    # haplotypecaller	(HC)

    def __init__(self):

        # Initialize the built-in supported variant callers
        self.callers: dict[str: VariantCaller] = {
            "BT": BCFTools(),
            "VS": Varscan(),
            "VD": Vardict(),
            "PL": Pindel(),
            "HC": Haplotypecaller(),
            "FL": Filt3r(),
            "DV": DeepVariant(),
        }
        # Save the plugins in a dictionary
        self.plugins: dict[str: VariantCaller] = {}
    
    # Get the variant caller by its name
    # If the caller is not supported, raise a VariantCallerError exception
    def get_VC(self, caller: str) -> VariantCaller:
        """
        Get a variant caller by its name.
        
        Args:
            caller (str): The name of the variant caller.
            
        Returns:
            VariantCaller: The variant caller instance.
            
        Raises:
            VariantCallerError: If the caller is not supported.
        """

        # Try to get the variant caller
        try:

            # Return the variant caller
            return self.callers[caller.upper()]

        # Catch a KeyError exception if the caller is not supported
        except KeyError:

            # Raise a VariantCallerError exception if the caller is not supported
            raise exceptions.VariantCallerError(f"Variant Caller not supported: {caller}")

    def get_BT(self) -> VariantCaller:
        """
        Get the BCFTools variant caller.
        
        Returns:
            VariantCaller: The BCFTools variant caller instance.
        """

        return self.get_VC("BT")

    def get_VS(self) -> VariantCaller:
        """
        Get the Varscan variant caller.
        
        Returns:
            VariantCaller: The Varscan variant caller instance.
        """

        return self.get_VC("VS")

    def get_VD(self) -> VariantCaller:
        """
        Get the Vardict variant caller.
        
        Returns:
            VariantCaller: The Vardict variant caller instance.
        """

        return self.get_VC("VD")

    def get_PL(self) -> VariantCaller:
        """
        Get the Pindel variant caller.
        
        Returns:
            VariantCaller: The Pindel variant caller instance.
        """

        return self.get_VC("PL")

    def get_HS(self) -> VariantCaller:
        """
        Get the Haplotypecaller variant caller.
        
        Returns:
            VariantCaller: The Haplotypecaller variant caller instance.
        """

        return self.get_VC("HS")

    def get_FL(self) -> VariantCaller:
        """
        Get the Filt3r variant caller.
        
        Returns:
            VariantCaller: The Filt3r variant caller instance.
        """

        return self.get_VC("FL")

    def get_DV(self) -> VariantCaller:
        """
        Get the DeepVariant caller.
        
        Returns:
            VariantCaller: The DeepVariant caller instance.
        """

        return self.get_VC("DV")

    def is_supported(self, caller: str) -> bool:
        """
        Check if a variant caller is supported.
        
        Args:
            caller (str): The name of the variant caller to check.
            
        Returns:
            bool: True if the caller is supported, False otherwise.
        """

        return caller.upper() in self.callers
    
    def is_plugin(self, caller: str) -> bool:
        """
        Check if a variant caller is a plugin.
        
        Args:
            caller (str): The name of the variant caller to check.
            
        Returns:
            bool: True if the caller is a plugin, False otherwise.
        """        

        return caller.upper() in self.plugins
    
    def load(self, plugin: object):
        """
        Load a variant caller plugin.

        Args:
            plugin (VariantCallerPlugin): The variant caller plugin to load.
        """

        try:

            # If the plugin package is not in the system path
            if not plugin.package in sys.path:

                # Add the plugin package to the system path
                sys.path.append(str(plugin.package))

            # Import the module
            module = importlib.import_module(plugin.config.params["caller"]["name"], package=str(plugin.package))

            # Get the variant caller class
            caller = getattr(module, plugin.config.params["caller"]["name"])()

            # Add the variant caller to the callers
            self.callers[plugin.id] = caller

            # Add the variant caller to the plugins
            self.plugins[plugin.id] = caller

        # Catch an ImportError exception if the module is not found
        except ImportError as e:

            # Trace the error
            logger.error(e)

            # Raise a VariantCallerPluginError exception
            raise exceptions.VariantCallerPluginError(e)

        # Catch an AttributeError exception if the class is not found
        except AttributeError as e:

            logger.error(e)

            raise exceptions.VariantCallerPluginError(e)
        
        # Catch other exceptions
        except Exception as e:

            # If the exception is an ImportError or an AttributeError
            if isinstance(e, ImportError) or isinstance(e, AttributeError):

                # Raise the exception
                raise

            # If the exception is not an ImportError or an AttributeError
            else:

                # Raise a VariantCallerPluginError exception
                raise exceptions.VariantCallerPluginError(f"An unexpected error has occurred when loading variant caller plugin: {e}")
            
    def populate(self, callers: list):
        """
        Populate the repository with variant callers.
        
        Args:
            callers (list): The list of variant callers to add.
        """

        # Loop through the variant callers
        for caller in callers:

            # Add the variant caller to the repository
            self.add(caller)
            
    def add(self, caller: object):
        """
        Add a variant caller to the repository.

        Args:
            caller (object): The variant caller to add.
        """
        
        try:
                            
            self.load(caller)

        except Exception:
            
            caller.remove()

            raise

    def remove(self, caller: object):
        """
        Remove a variant caller from the repository.

        Args:
            caller (object): The variant caller to remove.
        """

        # Remove the variant caller from the plugins
        self.plugins.pop(caller.id)

        # Remove the variant caller from the callers
        self.callers.pop(caller.id)
    
    def __len__(self):
        """
        Return the number of variant callers in the repository.
        """

        return len(self.callers)
    
    def __iter__(self):
        """
        Iterate over the variant callers in the repository.
        """

        yield from self.callers.items()

@final
class BCFTools(VariantCaller):

    """ A class to manage the BCFtools variant caller """

    # Set the FORMAT field keys as an enum
    FORMAT = enum.IntEnum(value="FORMAT",
                          names=','.join(["GT", "PL"]),
                          start=0)

    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """
        Check if the variant is compliant with the BCFtools format.
        
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
    
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:
        """
        Get the genotype of the variant.

        Args:
            variant (list[str]): The variant to get the genotype from.
            header (dict[str:int]): The header of the VCF file.
        """

        return variant[header["SAMPLE"]].split(':')[BCFTools.FORMAT.GT.value]

    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:
        """
        Calculate the variant allele frequency (VAF).

        Args:
            variant (list[str]): The variant to calculate the VAF from.
            header (dict[str:int]): The header of the VCF file.
        """

        total_depth: int = int(
            variant[header["INFO"]].split("DP=")[1].split(";")[0]
        )

        alleles_depth: int = (
            variant[header["INFO"]].split("DP4=")[1].split(";")[0]
        )

        variant_depth: int = int(alleles_depth.split(",")[2]) + int(
            alleles_depth.split(",")[3]
        )

        try:
            vaf: float = variant_depth / total_depth
        except ZeroDivisionError:
            vaf: float = 0.0

        return vaf

    @staticmethod
    def depth(variant: list[str], header: dict[int]) -> int:
        """
        Get the depth of the variant.

        Args:
            variant (list[str]): The variant to get the depth from.
            header (dict[int]): The header of the VCF file.
        """

        return int(variant[header["INFO"]].split("DP=")[1].split(";")[0])

    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the reference allele counts.

        Args:
            variant (list[str]): The variant to get the reference allele counts from.
            header (dict[str:int]): The header of the VCF file.
        """

        alleles_depth: list[str] = (
            variant[header["INFO"]].split("DP4=")[1].split(";")[0]
        )
        depths: list[str] = alleles_depth.split(",")
        rrc_plus: int = int(depths[0])
        rrc_minus: int = int(depths[1])
        rrc: int = rrc_plus + rrc_minus
        return (rrc, rrc_plus, rrc_minus)

    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the alternate allele counts.

        Args:
            variant (list[str]): The variant to get the alternate allele counts from.
            header (dict[str:int]): The header of the VCF file.
        """

        alleles_depth: list[str] = (
            variant[header["INFO"]].split("DP4=")[1].split(";")[0]
        )
        depths: list[str] = alleles_depth.split(",")
        arc_plus: int = int(depths[2])
        arc_minus: int = int(depths[3])
        arc: int = arc_plus + arc_minus
        return (arc, arc_plus, arc_minus)

    def __str__(self):
        """
        Return the string representation of the variant caller.
        """

        return "BCFTools"
    
    def __repr__(self):
        """
        Return the string representation of the variant caller.
        """

        return "BCFTools"

@final
class Varscan(VariantCaller):

    """ A class to manage the Varscan variant caller """

    # Set the FORMAT field keys as an enum
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

    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """
        Check if the variant is compliant with the Varscan format.
        
        Args:
            variant (list[str]): The variant to check.
            header (dict[str:int]): The header of the VCF file.
            
        Returns:
            bool: True if the variant is compliant, False otherwise.
        """

        return len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == [key.name for key in self.FORMAT]            
        ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
    
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:
        """
        Get the genotype of the variant.
        """

        return variant[header["SAMPLE"]].split(':')[Varscan.FORMAT.GT.value]

    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:
        """
        Get the variant allele frequency (VAF).
        """

        vaf = (
            float((variant[header["SAMPLE"]].split(":")[Varscan.FORMAT.FREQ.value]).rstrip('%')) / 100
        )

        return vaf

    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:
        """
        Get the depth of the variant.
        """

        return int(variant[header["SAMPLE"]].split(":")[Varscan.FORMAT.SDP.value])

    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the reference allele counts.
        """

        values = variant[header["SAMPLE"]].split(":")

        rrc = int(values[Varscan.FORMAT.RD.value])
        rrc_plus = int(values[Varscan.FORMAT.RDF.value])
        rrc_minus = int(values[Varscan.FORMAT.RDR.value])

        return (rrc, rrc_plus, rrc_minus)

    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the alternate allele counts.
        """

        values = variant[header["SAMPLE"]].split(":")

        arc = int(values[Varscan.FORMAT.AD.value])
        arc_plus = int(values[Varscan.FORMAT.ADF.value])
        arc_minus = int(values[Varscan.FORMAT.ADR.value])

        return (arc, arc_plus, arc_minus)

    def __str__(self):
        """
        Return the string representation of the variant caller.
        """

        return "Varscan"
    
    def __repr__(self):
        """
        Return the string representation of the variant caller.
        """

        return "Varscan"

@final
class Vardict(VariantCaller):

    """ A class to manage the Vardict variant caller """

    # Set the FORMAT field keys as an enum
    FORMAT = enum.IntEnum(value="FORMAT",
                          names=','.join(["GT", "DP", "VD", "AD", "AF", "RD", "ALD"]),
                          start=0)

    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """
        Check if the variant is compliant with the Vardict format.
        
        Args:
            variant (list[str]): The variant to check.
            header (dict[str:int]): The header of the VCF file.
            
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
    
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:
        """
        Get the genotype of the variant.
        """

        gt: str = variant[header["SAMPLE"]].split(':')[Vardict.FORMAT.GT.value]
        # Manage case when Vardict return 1/0 instead of 0/1
        return "0/1" if gt == "1/0" else gt

    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:
        """
        Get the variant allele frequency (VAF).
        """

        vaf = float(variant[header["INFO"]].split("AF=")[1].split(";")[0])

        return vaf

    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:
        """
        Get the depth of the variant.
        """

        return int(variant[header["SAMPLE"]].split(":")[Vardict.FORMAT.DP.value])

    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the reference allele counts.
        """

        values: list[str] = variant[header["SAMPLE"]].split(":")

        rrc = int(values[3].split(",")[0])

        rd: list[str] = values[Vardict.FORMAT.RD.value].split(",")
        rrc_plus: int = int(rd[0])
        rrc_minus: int = int(rd[1])

        return (rrc, rrc_plus, rrc_minus)

    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the alternate allele counts.
        """

        values: list[str] = variant[header["SAMPLE"]].split(":")

        arc = int(values[Vardict.FORMAT.AD.value].split(",")[1])

        ald: list[int] = values[Vardict.FORMAT.ALD.value].split(",")
        arc_plus: int = int(ald[0])
        arc_minus: int = int(ald[1])

        return (arc, arc_plus, arc_minus)

    def __str__(self):
        """
        Return the string representation of the variant caller.
        """

        return "Vardict"
    
    def __repr__(self):
        """
        Return the string representation of the variant caller.
        """

        return "Vardict"

@final
class Pindel(VariantCaller):
    
    """ A class to manage the Pindel variant caller """

    # Set the FORMAT field keys as an enum
    FORMAT = enum.IntEnum(value="FORMAT",
                          names=','.join(["GT", "AD"]),
                          start=0)

    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """
        Check if the variant is compliant with the Pindel format.
        
        Args:
            variant (list[str]): The variant to check.
            header (dict[str:int]): The header of the VCF file.
            
        Returns:
            bool: True if the variant is compliant, False otherwise.
        """
        
        return len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == [key.name for key in self.FORMAT]    
        ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
    
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:
        """
        Get the genotype of the variant.
        """

        return variant[header["SAMPLE"]].split(':')[Pindel.FORMAT.GT.value]

    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:
        """
        Get the variant allele frequency (VAF).
        """

        depths = variant[header["SAMPLE"]].split(":")[Pindel.FORMAT.AD.value].split(",")
        try:
            vaf: float = float(depths[1]) / (float(depths[0]) + float(depths[1]))
        except ZeroDivisionError:
            vaf: float = 0.0
            
        return vaf

    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:
        """
        Get the depth of the variant.
        """

        alleles_depth: list[str] = (
            variant[header["SAMPLE"]].split(":")[Pindel.FORMAT.AD.value].split(",")
        )

        return int(alleles_depth[0]) + int(alleles_depth[1])

    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the reference allele counts.
        """

        return (None, None, None) # Pindel does not provide the reference allele counts

    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the alternate allele counts.
        """

        arc = int(variant[header["SAMPLE"]].split(":")[Pindel.FORMAT.AD.value].split(",")[1])

        return (arc, None, None) # Pindel does not provide the alternate allele counts for forward and reverse strands

    def __str__(self):
        """
        Return the string representation of the variant caller.
        """

        return "Pindel"
    
    def __repr__(self):
        """
        Return the string representation of the variant caller.
        """

        return "Pindel"

@final
class Haplotypecaller(VariantCaller):
    
    """ A class to manage the Haplotypecaller variant caller """

    # Set the FORMAT field keys as an enum
    FORMAT = enum.IntEnum(value="FORMAT",
                          names=','.join(["GT", "AD", "DP", "GQ", "PL"]),
                          start=0)

    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """
        Check if the variant is compliant with the Haplotypecaller format.
        
        Args:
            variant (list[str]): The variant to check.
            header (dict[str:int]): The header of the VCF file.
            
        Returns:
            bool: True if the variant is compliant, False otherwise.
        """
        
        return len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == [key.name for key in self.FORMAT]            
        ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
    
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:
        """
        Get the genotype of the variant.
        """

        return variant[header["SAMPLE"]].split(':')[Haplotypecaller.FORMAT.GT.value]

    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:
        """
        Get the variant allele frequency (VAF).
        """

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
        """
        Get the depth of the variant.
        """

        return int(variant[header["SAMPLE"]].split(":")[Haplotypecaller.FORMAT.DP.value])

    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the reference allele counts.
        """

        return (None, None, None) # Haplotypecaller does not provide the reference allele counts for each strand

    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the alternate allele counts.
        """

        arc = int(variant[header["SAMPLE"]].split(":")[Haplotypecaller.FORMAT.AD.value].split(",")[1])

        return (arc, None, None) # Haplotypecaller does not provide the alternate allele counts for forward and reverse strands

    def __str__(self):
        """
        Return the string representation of the variant caller.
        """

        return "Haplotypecaller"
    
    def __repr__(self):
        """
        Return the string representation of the variant caller.
        """

        return "Haplotypecaller"

@final
class Filt3r(VariantCaller):
    
    """ A class to manage the Filt3r variant caller """

    # Filt3r FORMAT
    FORMAT = []

    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:
        """
        Get the genotype of the variant.
        """

        raise NotImplementedError

    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:
        """
        Get the variant allele frequency (VAF).
        """

        return float(variant[6].split(";")[4].split("=")[1])

    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:
        """
        Get the depth of the variant.
        """

        infos: list[str] = variant[6].split(";")

        variant_reads: int = int(infos[1].split("=")[1])

        reference_reads: int = int(infos[2].split("=")[1])

        return variant_reads + reference_reads

    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the reference allele counts.
        """

        return (None, None, None) # Filt3r does not provide the reference allele counts for each strand

    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the alternate allele counts.
        """

        arc = int(variant[6].split(";")[1].split("=")[1])

        return (arc, None, None) # Filt3r does not provide the alternate allele counts for forward and reverse strands

    def __str__(self):
        """
        Return the string representation of the variant caller.
        """

        return "Flit3r"
    
    def __repr__(self):
        """
        Return the string representation of the variant caller.
        """

        return "Flit3r"

@final
class DeepVariant(VariantCaller):
    
    """ A class to manage the DeepVariant variant caller """

    # Set the FORMAT field keys as an enum
    FORMAT = enum.IntEnum(value="FORMAT",
                          names=','.join(["GT", "GQ", "DP", "AD", "VAF", "PL"]),
                          start=0)

    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """
        Check if the variant is compliant with the DeepVariant format.
        
        Args:
            variant (list[str]): The variant to check.
            header (dict[str:int]): The header of the VCF file.
            
        Returns:
            bool: True if the variant is compliant, False otherwise.
        """
        
        return len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == [key.name for key in self.FORMAT]            
        ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
    
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:
        """
        Get the genotype of the variant.
        """

        return variant[header["SAMPLE"]].split(':')[DeepVariant.FORMAT.GT.value]

    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:
        """
        Get the variant allele frequency (VAF).
        """

        return float(variant[header["SAMPLE"]].split(":")[DeepVariant.FORMAT.VAF.value])

    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:
        """
        Get the depth of the variant.
        """

        return float(variant[header["SAMPLE"]].split(":")[DeepVariant.FORMAT.DP.value])

    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the reference allele counts.
        """

        return (
            int(variant[header["SAMPLE"]].split(":")[DeepVariant.FORMAT.AD.value].split(",")[0]),
            None,
            None,
        ) # DeepVariant does not provide the reference allele counts for each strand

    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """
        Get the alternate allele counts.
        """

        arc = int(variant[header["SAMPLE"]].split(":")[DeepVariant.FORMAT.AD.value].split(",")[1])

        return (arc, None, None) # DeepVariant does not provide the alternate allele counts for forward and reverse strands

    def __str__(self):
        """
        Return the string representation of the variant caller.
        """

        return "Deepvariant"

    def __repr__(self):
        """
        Return the string representation of the variant caller.
        """

        return "Deepvariant"
