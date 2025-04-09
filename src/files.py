import callers

class GenomicFile():

    def __init__(self, path: str):

        self.path = path

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
        "SAMPLE": 9
    }

    def __init__(self, path: str, caller: callers.VariantCaller, lazy: bool = True):

        super().__init__(path=path)

        self.verify(self.path)

        self.caller: callers.VariantCaller = caller

        if not lazy:

            self.parse(self.path)

    def get_header(self):

        return self.HEADER
    
    def parse(self):

        pass

    def verify(self):

        pass

    @staticmethod
    def convert(a: object) -> object:
        """ Convert variable to appropriate type """
        try:
            # If the variable contains a '/' or '|' character, it is a genotype information, return the variable as is
            # Else return the variable as an evaluated expression
            return a if sum(list(map(lambda x: x in a,('/','|')))) else eval(a)
        except Exception:
            # If the variable cannot be evaluated, return the variable as is
            return a

    def format_to_values(self, values: str|list[str]) -> dict:
        """ map FORMAT string to respective SAMPLE values """

        # Split the values string into a list of fields
        values: list[str] = values.split(":")
        return {f: VCF.convert(v) for f, v in zip(self.caller.FORMAT, values)}

class Pileup(GenomicFile):

    HEADER: dict[str:int] = {"barcode": 0, 
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
              "Ins": 14, 
              "Del": 15}
    
    PLUS_STRAND: list[int] = [0,2,4,6]
    MINUS_STRAND: list[int] = [1,3,5,7]

    def __init__(self, path: str, lazy: bool = True):

        super().__init__(path=path)

        self.verify(self.path)

        if not lazy:

            self.parse(self.path)

    def get_header(self):

        return self.HEADER
    
    def parse(self):

        pass

    def verify(self):

        pass

class FastaIndex(GenomicFile):

    def __init__(self, path: str, lazy: bool = True):

        super().__init__(path)

        self.verify(self.path)

        if not lazy:

            self.parse(self.path)

    def parse(self):

        pass

    def verify(self, path):

        pass