class VariantCallerRepository():

    def __init__(self):
        
        self.BT = BCFTools()
        self.VS = Varscan()
        self.VD = Vardict()
        self.PL = Pindel()
        self.HS = Haplotypecaller()
        self.FL = Filt3r()
        self.DV = DeepVariant()

    def get_BT(self):

        return self.BT
    
    def get_VS(self):

        return self.VS
    
    def get_VD(self):

        return self.VD
    
    def get_PL(self):

        return self.PL
    
    def get_HS(self):

        return self.HS
    
    def get_FL(self):

        return self.FL
    
    def get_DV(self):

        return self.DV
    

class VariantCaller():

    def __init__(self):

        pass


class BCFTools(VariantCaller):

    FORMAT = ["GT", "PL"]

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



class Filt3r(VariantCaller):

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

