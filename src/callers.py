import errors
from typing import final


class VariantCaller:

    """ Base class for variant callers """

    def __init__(self):

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

        try:

            return self.callers[caller.upper()]

        except KeyError:

            raise errors.VariantCallerError(
                f"Variant Caller {caller} not supported."
            )

    def get_BT(self) -> VariantCaller:

        return self.callers["BT"]

    def get_VS(self) -> VariantCaller:

        return self.callers["VS"]

    def get_VD(self) -> VariantCaller:

        return self.callers["VD"]

    def get_PL(self) -> VariantCaller:

        return self.callers["PL"]

    def get_HS(self) -> VariantCaller:

        return self.callers["HS"]

    def get_FL(self) -> VariantCaller:

        return self.callers["FL"]

    def get_DV(self) -> VariantCaller:

        return self.callers["DV"]

    def is_supported(self, caller: str) -> bool:

        return caller in self.callers
    
    def __len__(self):

        return len(self.callers)

@final
class BCFTools(VariantCaller):

    """ A class to manage the BCFtools variant caller """

    # BCFtools FORMAT
    FORMAT = ["GT", "PL"]

    def __init__(self):

        super().__init__()

    # Check if the variant is compliant with the BCFtools format.
    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):

        return (
            len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == self.FORMAT            
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

        return variant[header["SAMPLE"]].split(':')[0]

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

@final
class Varscan(VariantCaller):

    """ A class to manage the Varscan variant caller """

    # Varscan FORMAT
    FORMAT = [
        "GT",
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
    ]

    def __init__(self):

        super().__init__()

    # Check if the variant is compliant with the Varscan format.
    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):

        return len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == self.FORMAT            
        ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
    
    # Extract the genotype from the variant
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:

        return variant[header["SAMPLE"]].split(':')[0]

    # Calculate the variant allele frequency (VAF)
    # The VAF is the depth of the variant allele divided by the total depth
    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:

        vaf = (
            float(variant[header["SAMPLE"]].split("%")[0].split(":")[-1]) / 100
        )

        return vaf

    # Extract the depth of the variant
    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:

        return int(variant[header["SAMPLE"]].split(":")[2])

    # Extract the reference allele counts
    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:

        values = variant[header["SAMPLE"]].split(":")

        rrc = int(values[4])
        rrc_plus = int(values[10])
        rrc_minus = int(values[11])

        return (rrc, rrc_plus, rrc_minus)

    # Extract the alternate allele counts
    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:

        values = variant[header["SAMPLE"]].split(":")

        arc = int(values[5])
        arc_plus = int(values[12])
        arc_minus = int(values[13])

        return (arc, arc_plus, arc_minus)

    def __str__(self):
        return "Varscan"

@final
class Vardict(VariantCaller):

    """ A class to manage the Vardict variant caller """

    # Vardict FORMAT
    FORMAT = ["GT", "DP", "VD", "AD", "AF", "RD", "ALD"]

    def __init__(self):
        super().__init__()

    # Check if the variant is compliant with the Vardict format.
    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):

        return (
            len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == self.FORMAT            
            ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
        ) and (
            (len(variant[header["INFO"]]))
            and ("AF=" in variant[header["INFO"]])
        )
    
    # Extract the genotype from the variant
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:

        gt: str = variant[header["SAMPLE"]].split(':')[0]
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

        return int(variant[header["SAMPLE"]].split(":")[1])

    # Extract the reference allele counts
    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:

        values: list[str] = variant[header["SAMPLE"]].split(":")

        rrc = int(values[3].split(",")[0])

        rd: list[str] = values[5].split(",")
        rrc_plus: int = int(rd[0])
        rrc_minus: int = int(rd[1])

        return (rrc, rrc_plus, rrc_minus)

    # Extract the alternate allele counts
    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:

        values: list[str] = variant[header["SAMPLE"]].split(":")

        arc = int(values[3].split(",")[1])

        ald: list[int] = values[6].split(",")
        arc_plus: int = int(ald[0])
        arc_minus: int = int(ald[1])

        return (arc, arc_plus, arc_minus)

    def __str__(self):
        return "Vardict"

@final
class Pindel(VariantCaller):
    
    """ A class to manage the Pindel variant caller """

    # Pindel FORMAT
    FORMAT = ["GT", "AD"]

    def __init__(self):
        super().__init__()

    # Check if the variant is compliant with the Pindel format.
    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):

        return len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == self.FORMAT            
        ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
    
    # Extract the genotype from the variant
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:

        return variant[header["SAMPLE"]].split(':')[0]

    # Calculate the variant allele frequency (VAF)
    # The VAF is the depth of the variant allele divided by the total depth
    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:

        depths = variant[header["SAMPLE"]].split(":")[1].split(",")
        try:
            vaf: float = float(depths[1]) / (float(depths[0]) + float(depths[1]))
        except ZeroDivisionError:
            vaf: float = 0.0
            
        return vaf

    # Extract the depth of the variant
    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:

        alleles_depth: list[str] = (
            variant[header["SAMPLE"]].split(":")[1].split(",")
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

        arc = int(variant[header["SAMPLE"]].split(":")[1].split(",")[1])

        return (arc, None, None)

    def __str__(self):

        return "Pindel"

@final
class Haplotypecaller(VariantCaller):
    
    """ A class to manage the Haplotypecaller variant caller """

    # Haplotypecaller FORMAT
    FORMAT = ["GT", "AD", "DP", "GQ", "PL"]

    def __init__(self):
        super().__init__()

    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """Check if the variant is compliant with the Haplotypecaller format."""

        return len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == self.FORMAT            
        ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
    
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:
        """Extract the genotype from the variant."""

        return variant[header["SAMPLE"]].split(':')[0]

    # The VAF is the depth of the variant allele divided by the total depth
    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:
        """Calculate the variant allele frequency (VAF)."""

        metrics: list[str] = variant[header["SAMPLE"]].split(":")

        total_depth = float(metrics[2])
        alleles_depth = float(metrics[1].split(",")[1])

        try:
            vaf: float = alleles_depth / total_depth
        except ZeroDivisionError:
            vaf: float = 0.0

        return vaf

    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:
        """Extract the depth of the variant."""

        return int(variant[header["SAMPLE"]].split(":")[2])

    # Haplotypecaller does not provide the reference allele counts for each strand
    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """Extract the reference allele counts."""

        return (None, None, None)

    # Haplotypecaller does not provide the alternate allele counts for each strand
    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """Extract the alternate allele counts."""

        arc = int(variant[header["SAMPLE"]].split(":")[1].split(",")[1])

        return (arc, None, None)

    def __str__(self):

        return "Haplotypecaller"

@final
class Filt3r(VariantCaller):
    
    """ A class to manage the Filt3r variant caller """

    # Filt3r FORMAT
    FORMAT = []

    def __init__(self):
        super().__init__()

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

@final
class DeepVariant(VariantCaller):
    
    """ A class to manage the DeepVariant variant caller """

    # DeepVariant FORMAT
    FORMAT = ["GT", "GQ", "DP", "AD", "VAF", "PL"]

    def __init__(self):

        super().__init__()

    # Ensure the metrics can be later extracted.
    def is_compliant(self, variant: list[str], header: dict[str:int]):
        """Check if the variant is compliant with the DeepVariant format."""

        return len(variant[header["SAMPLE"]]) and (
            variant[header["FORMAT"]].split(':') == self.FORMAT            
        ) and (len(variant[header["SAMPLE"]].split(':')) == len(self.FORMAT))
    
    @staticmethod
    def genotype(variant: list[str], header: dict[str:int]) -> str:
        """Extract the genotype from the variant."""

        return variant[header["SAMPLE"]].split(':')[0]

    # The VAF is the depth of the variant allele divided by the total depth
    @staticmethod
    def VAF(variant: list[str], header: dict[str:int]) -> float:
        """Calculate the variant allele frequency (VAF)."""

        return float(variant[header["SAMPLE"]].split(":")[-2])

    @staticmethod
    def depth(variant: list[str], header: dict[str:int]) -> int:
        """Extract the depth of the variant."""

        return float(variant[header["SAMPLE"]].split(":")[2])

    # DeepVariant does not provide the reference allele counts for each strand
    @staticmethod
    def rrc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """Extract the reference allele counts."""

        return (
            int(variant[header["SAMPLE"]].split(":")[3].split(",")[0]),
            None,
            None,
        )

    # DeepVariant does not provide the alternate allele counts for each strand
    @staticmethod
    def arc(variant: list[str], header: dict[str:int]) -> tuple[int]:
        """Extract the alternate allele counts."""

        arc = int(variant[header["SAMPLE"]].split(":")[3].split(",")[1])

        return (arc, None, None)

    def __str__(self):
        return "Deepvariant"
