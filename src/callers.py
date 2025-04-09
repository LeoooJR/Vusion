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
                        "HS": Haplotypecaller(),
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
        except ValueError:
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
            vaf = 0

        return vaf
    
    @staticmethod
    def depth(variant: str) -> int:

        return float(variant[7].split('DP=')[1].split(';')[0])
    
    @staticmethod
    def rcc(variant: str) -> float:

        variant_depth = variant[7].split('DP4=')[1].split(';')[0]
        rrc_plus = float(variant_depth.split(',')[0])
        rrc_minus = float(variant_depth.split(',')[1])
        rrc = rrc_plus + rrc_minus
        return(rrc, rrc_plus, rrc_minus)
    
    @staticmethod
    def arc(variant: str) -> float:

        variant_depth = variant[7].split('DP4=')[1].split(';')[0]
        arc_plus = float(variant_depth.split(',')[2])
        arc_minus = float(variant_depth.split(',')[3])
        arc = arc_plus + arc_minus
        return(arc,arc_plus,arc_minus)


class Varscan(VariantCaller):

    FORMAT = ["GT", "GQ", "SDP", "DP", "RD", "AD", "FREQ", "PVAL", "RBQ", "ABQ", "RDF", "RDR", "ADF", "ADR"]

    def __init__(self):
        super().__init__()

    @staticmethod
    def VAF(variant: str) -> float:

        pass
    
    @staticmethod
    def depth(variant: str) -> int:

        pass
    
    @staticmethod
    def rcc(variant: str) -> float:

        pass
    
    @staticmethod
    def arc(variant: str) -> float:

        pass



class Vardict(VariantCaller):

    FORMAT = ["GT", "DP", "VD", "AD", "AF", "RD", "ALD"]

    def __init__(self):
        super().__init__()

    @staticmethod
    def VAF(variant: str) -> float:

        pass
    
    @staticmethod
    def depth(variant: str) -> int:

        pass
    
    @staticmethod
    def rcc(variant: str) -> float:

        pass
    
    @staticmethod
    def arc(variant: str) -> float:

        pass



class Pindel(VariantCaller):

    FORMAT = ["GT", "AD"]

    def __init__(self):
        super().__init__()

    @staticmethod
    def VAF(variant: str) -> float:

        pass
    
    @staticmethod
    def depth(variant: str) -> int:

        pass
    
    @staticmethod
    def rcc(variant: str) -> float:

        pass
    
    @staticmethod
    def arc(variant: str) -> float:

        pass



class Haplotypecaller(VariantCaller):

    FORMAT = ["GT", "AD", "DP", "GQ", "PL"]

    def __init__(self):
        super().__init__()

    @staticmethod
    def VAF(variant: str) -> float:

        pass
    
    @staticmethod
    def depth(variant: str) -> int:

        pass
    
    @staticmethod
    def rcc(variant: str) -> float:

        pass
    
    @staticmethod
    def arc(variant: str) -> float:

        pass



class Flit3r(VariantCaller):

    FORMAT = []

    def __init__(self):
        super().__init__()

    @staticmethod
    def VAF(variant: str) -> float:

        pass
    
    @staticmethod
    def depth(variant: str) -> int:

        pass
    
    @staticmethod
    def rcc(variant: str) -> float:

        pass
    
    @staticmethod
    def arc(variant: str) -> float:

        pass



class DeepVariant(VariantCaller):

    FORMAT = ["GT", "GQ", "DP", "AD", "VAF", "PL"]

    def __init__(self):
        super().__init__()

    @staticmethod
    def VAF(variant: str) -> float:

        pass
    
    @staticmethod
    def depth(variant: str) -> int:

        pass
    
    @staticmethod
    def rcc(variant: str) -> float:

        pass
    
    @staticmethod
    def arc(variant: str) -> float:

        pass

