import callers
import errors
import jinja2
import os

class GenomicWritter():

    def __init__(self, file: str):
        
        self.file = file

    def writeVCF(self, contigs: object, variants: object, samples: list[str], thresholds: list[float]):

        def sort_variant(contigs: object , variants: dict) -> object:

            return {k: v for k, v in sorted(variants.items(), key=lambda item: [item[1]["VC"]["CHROM"], int(item[1]["VC"]["POS"])] )}

        FORMAT = [
        'GT', 'VAR', 'BKG', 'TRC', 'RRC', 'ARC', 'BRC', 'ARR',
        'BRR', 'BRE', 'SBP', 'SBM', 'LOW', 'VCI', 'VCN', 'PIL', 'RES'
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

        HEADER.extend(samples)

        INFOS: list[str] = ["VAR"]

        ressources = os.path.join(os.path.dirname(os.path.abspath(__file__)),'templates')

        env = jinja2.Environment(loader=jinja2.FileSystemLoader(ressources))

        template = env.get_template("header")

        with open(self.file, mode='w') as out:

            out.writelines(template.render(contigs = contigs, thresholds = thresholds))

            out.write(f"\n#{'\t'.join(HEADER)}\n")

            variants = sort_variant(contigs=contigs, variants=variants)

            for variant in variants:

                out.write('\t'.join([variants[variant]["VC"]["CHROM"], # Chromosome field
                        variants[variant]["VC"]["POS"], # Position field
                        '.', # ID field
                        variants[variant]["VC"]["REF"], # Reference field
                        variants[variant]["VC"]["ALT"], # Alternate field
                        '.',
                        variants[variant]["final_metrics"]["FILTER"], # Filter field
                        '='.join([INFOS[0],variants[variant]["VT"]]), # Info field
                        ':'.join(FORMAT), # Format field
                        ':'.join([str(variants[variant]["final_metrics"][f]) for f in FORMAT])])) # Sample values field
                
                out.write('\n')

class GenomicReader():

    pass

class GenomicFile():

    def __init__(self, path: str):

        self.path = path

    def is_empty(self) -> bool:
        """ Check if file is empty """
        return os.path.getsize(self.path) == 0

    def is_file(self) -> bool:
        """ Check if path is a file """
        return os.path.isfile(self.path)
    
    def get_path(self):

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
        "SAMPLE": 9
    }

    def __init__(self, path: str, caller: callers.VariantCaller, lazy: bool = True):

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

            raise errors.VCFError(f"Error: The file {self.path} does not exist.")

        if self.is_empty():

            raise errors.VCFError(f"Error: The file {self.path} is empty.")
        
        try:

            with open(self.path, mode='r') as vcf:

                line = vcf.readline()

                if not line:

                    raise errors.VCFError(f"Error: First line of {self.path} is empty.")
                
                else:

                    # Check if first line start with "#"
                    if line[0] != '#':

                        raise errors.VCFError(f"Error: First line inconsistent with VCF header format")

        except FileNotFoundError:

            raise errors.VCFError(f"{self.path} is not a valid path")
                
        except IOError:

            raise errors.VCFError(f"An error occurred while reading {self.path}")

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
    
    def info_to_values(self, values: str) -> dict:

        infos = list(map(lambda item: item.split('='), values.split(';')))

        return {k: v for k, v in infos}
    
    def VAF(self, variant: str) -> float:

        return self.caller.VAF(variant)

    def depth(self, variant: str) -> int:

        return self.caller.depth(variant)

    def arc(self, variant: str) -> tuple[float]:

        return self.caller.arc(variant)

    def rrc(self, variant: str) -> tuple[float]:

        return self.caller.rrc(variant)

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
              "INS": 14, 
              "DEL": 15}
    
    PLUS_STRAND: list[int] = [0,2,4,6]
    MINUS_STRAND: list[int] = [1,3,5,7]

    def __init__(self, path: str, lazy: bool = True):

        super().__init__(path=path)

        self.verify()

        if not lazy:

            self.parse()

    def get_header(self):

        return self.HEADER
    
    def parse(self):

        pass

    def verify(self):

        if not self.is_file():

            raise errors.PileupError(f"Error: The file {self.path} does not exist.")

        if self.is_empty():

            raise errors.PileupError(f"Error: The file {self.path} is empty.")
        
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

    def parse(self):

        pass

    def verify(self):

        if not self.is_file():

            raise errors.FastaIndexError(f"Error: The file {self.path} does not exist.")

        if self.is_empty():

            raise errors.FastaIndexError(f"Error: The file {self.path} is empty.")
        
        try:

            with open(self.path, mode='r') as fasta:

                line = fasta.readline()

                if not line:
                    
                    raise errors.FastaIndexError(f"Error: First line of {self.path} is empty.")
                
                else:

                    # Check if first line is composed of 5 columns
                    if len(line.split('\t')) != 5:

                        raise errors.FastaIndexError(f"Error: First line inconsistent with Fasta index format")
        
        except FileNotFoundError:

            raise errors.FastaIndexError(f"{self.path} is not a valid path")

        except IOError:

            raise errors.FastaIndexError(f"An error occurred while reading {self.path}")