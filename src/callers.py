import errors

class VariantCaller():

    def __init__(self):

        pass

class VariantCallerRepository():

    # known variant callers are:
    # samtools			(ST)
    # varscan			(VS)
    # vardict			(VD)
    # pindel			(PL)
    # haplotypecaller	(HC)
    # smCounter2		(SM; DEPRECATED)
    # control & hotspot (CS & HS) ## when --hotspot option is set;
    # CtlSet & HotSpot should both originate from CombineVCF2final_metrics

    def __init__(self):

        self.callers = {"BT": BCFTools(),
                        "VS": Varscan(),
                        "VD": Vardict(),
                        "PL": Pindel(),
                        "HC": Haplotypecaller(),
                        "FL": Flit3r(),
                        "DV": DeepVariant()}
        
        # self.BT = BCFTools()
        # self.VS = Varscan()
        # self.VD = Vardict()
        # self.PL = Pindel()
        # self.HS = Haplotypecaller()
        # self.FL = Flit3r()
        # self.DV = DeepVariant()

    def get_VC(self, caller: str) -> VariantCaller:

        try:
            return self.callers[caller]
        except KeyError:
            raise errors.VariantCallerError(f"Variant Caller {caller} not supported.")

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

class BCFTools(VariantCaller):

    FORMAT = ["GT", "PL"]

    def __init__(self):
        super().__init__()

    @staticmethod
    def VAF(variant: str) -> float:

        total_depth = float(variant[7].split('DP=')[1].split(';')[0])

        alleles_depth = variant[7].split('DP4=')[1].split(';')[0]

        variant_depth = float(alleles_depth.split(',')[2]) + float(alleles_depth.split(',')[3])

        try:
            vaf = variant_depth / total_depth
        except ZeroDivisionError:
            vaf = 0.0

        return vaf
    
    @staticmethod
    def depth(variant: str) -> int:

        return int(variant[7].split('DP=')[1].split(';')[0])
    
    @staticmethod
    def rrc(variant: str) -> tuple[float]:

        variant_depth = variant[7].split('DP4=')[1].split(';')[0]
        depths = variant_depth.split(',')
        rrc_plus = float(depths[0])
        rrc_minus = float(depths[1])
        rrc = rrc_plus + rrc_minus
        return(rrc, rrc_plus, rrc_minus)
    
    @staticmethod
    def arc(variant: str) -> tuple[float]:

        variant_depth = variant[7].split('DP4=')[1].split(';')[0]
        depths = variant_depth.split(',')
        arc_plus = float(depths[2])
        arc_minus = float(depths[3])
        arc = arc_plus + arc_minus
        return(arc,arc_plus,arc_minus)


class Varscan(VariantCaller):

    FORMAT = ["GT", "GQ", "SDP", "DP", "RD", "AD", "FREQ", "PVAL", "RBQ", "ABQ", "RDF", "RDR", "ADF", "ADR"]

    def __init__(self):
        super().__init__()

    @staticmethod
    def VAF(variant: str) -> float:

        vaf = float(variant[9].split('%')[0].split(':')[-1])/100

        return vaf
    
    @staticmethod
    def depth(variant: str) -> int:

        return int(variant[9].split(':')[2])
    
    @staticmethod
    def rrc(variant: str) -> tuple[float]:

        values = variant[9].split(':')

        rrc = float(values[4])
        rrc_plus = float(values[10])
        rrc_minus = float(values[11])

        return(rrc,rrc_plus,rrc_minus)
    
    @staticmethod
    def arc(variant: str) -> tuple[float]:

        values = variant[9].split(':')

        arc = float(values[5])
        arc_plus = float(values[12])
        arc_minus = float(values[13])

        return(arc,arc_plus,arc_minus)



class Vardict(VariantCaller):

    FORMAT = ["GT", "DP", "VD", "AD", "AF", "RD", "ALD"]

    def __init__(self):
        super().__init__()

    @staticmethod
    def VAF(variant: str) -> float:

        vaf = float(variant[7].split('AF=')[1].split(';')[0])

        return vaf
    
    @staticmethod
    def depth(variant: str) -> int:

        return int(variant[9].split(":")[1])
    
    @staticmethod
    def rrc(variant: str) -> tuple[float]:

        values = variant[9].split(":")

        rrc = float(values[3].split(',')[0])
        rrc_plus = float(values[5].split(',')[0])
        rrc_minus = float(values[5].split(',')[1])

        return(rrc,rrc_plus,rrc_minus)
    
    @staticmethod
    def arc(variant: str) -> tuple[float]:

        values = variant[9].split(":")

        arc = float(values[3].split(',')[1])
        arc_plus = float(values[6].split(',')[0])
        arc_minus = float(values[6].split(',')[1])

        return(arc,arc_plus,arc_minus)

class Pindel(VariantCaller):

    FORMAT = ["GT", "AD"]

    def __init__(self):
        super().__init__()

    @staticmethod
    def VAF(variant: str) -> float:

        depths = variant[9].split(':')[1].split(',')
        vaf = float(depths[1])/(float(depths[0]) + float(depths[1]))
        return vaf
    
    @staticmethod
    def depth(variant: str) -> int:

        return (
        float(variant[9].split(':')[1].split(',')[0]) +
        float(variant[9].split(':')[1].split(',')[1])
        )
    
    @staticmethod
    def rrc(variant: str) -> tuple[float]:

        pass
    
    @staticmethod
    def arc(variant: str) -> tuple[float]:

        return float(variant[9].split(':')[1].split(',')[1])



class Haplotypecaller(VariantCaller):

    FORMAT = ["GT", "AD", "DP", "GQ", "PL"]

    def __init__(self):
        super().__init__()

    @staticmethod
    def VAF(variant: str) -> float:

        total_depth = float(variant[9].split(':')[2])
        variant_depth = float(variant[9].split(':')[1].split(',')[1])
        vaf = variant_depth / total_depth
        return vaf
    
    @staticmethod
    def depth(variant: str) -> int:

        return float(variant[9].split(':')[2])
    
    @staticmethod
    def rrc(variant: str) -> tuple[float]:

        pass
    
    @staticmethod
    def arc(variant: str) -> tuple[float]:

        return float(variant[9].split(':')[1].split(',')[1])



class Flit3r(VariantCaller):

    FORMAT = []

    def __init__(self):
        super().__init__()

    @staticmethod
    def VAF(variant: str) -> float:

        return float(variant[6].split(';')[4].split('=')[1])
    
    @staticmethod
    def depth(variant: str) -> int:

        return (
            float(variant[6].split(';')[1].split('=')[1]) +
            float(variant[6].split(';')[2].split('=')[1])
        )
    
    @staticmethod
    def rrc(variant: str) -> tuple[float]:

        pass
    
    @staticmethod
    def arc(variant: str) -> tuple[float]:

        return float(variant[6].split(';')[1].split('=')[1])



class DeepVariant(VariantCaller):

    FORMAT = ["GT", "GQ", "DP", "AD", "VAF", "PL"]

    def __init__(self):
        super().__init__()

    @staticmethod
    def VAF(variant: str) -> float:

        return variant.split(':')[-2]
    
    @staticmethod
    def depth(variant: str) -> int:

        return variant.split(':')[2]
    
    @staticmethod
    def rrc(variant: str) -> tuple[float]:

        return variant.split(':')[3].split(',')[0]
    
    @staticmethod
    def arc(variant: str) -> tuple[float]:

        return variant.split(':')[3].split(',')[1]

