from collections import deque, Counter, namedtuple
from files import Pileup
from functools import partial, lru_cache
from loguru import logger
import math
import re
import utils as functions

class VariantsRepository():

    def __init__(self, rescue: bool = False):

        # Dictionary to store variants
        self.variants: dict = {}

        self.cache: dict[str:set[str]] = {"common": set(), # Set to store variants possibly found in pileup
                                          "complex": set(), # Set to store complex variants
                                          "rejected": set()} # Set to store rejected variants
        
        self.rescue: bool = rescue

    def set_pileup(self, pileup: Pileup):

        self.pileup = pileup

    @staticmethod
    def get_variant_type(ref: str, alt: str) -> str:

        # Empty string to store variant type
        # Empty string return False
        variant_type: str = ''

        def is_snp(ref: str, alt: str) -> str:

            return "SNV" if len(ref) == 1 and len(alt) == 1 else ''

        def is_ins(ref: str, alt: str) -> str:
            
            return "INS" if len(ref) == 1 and len(alt) > 1 else ''

        def is_del(ref: str, alt: str) -> str:

            return "DEL" if len(ref) > 1 and len(alt) == 1 else ''

        def is_inv(ref: str, alt: str) -> str:

            OLD_CHARS: str = "ACGTacgt"
            REPLACE_CHARS: str = "TGCAtgca"

            rev: str = alt.translate(str.maketrans(OLD_CHARS,REPLACE_CHARS))[::-1]

            return "INV" if len(ref) == len(alt) and ref == rev else ''

        def is_mnv(ref: str, alt: str) -> str:

            OLD_CHARS: str = "ACGTacgt"
            REPLACE_CHARS: str = "TGCAtgca"

            rev: str = alt.translate(str.maketrans(OLD_CHARS,REPLACE_CHARS))[::-1]

            return "MNV" if len(ref) == len(alt) and ref != rev else ''
        
        # List of functions to check variant type
        # Deque is used to pop functions from the left at each iteration
        # until one of them returns a value
        # If no function returns a value, the variant type is set to CSV (complex structural variant)
        funcs: deque = deque([is_snp, is_ins, is_del, is_inv, is_mnv])
        
        try:
            # Loop until one of the functions returns a value
            # Empty string equivalent to False
            while not variant_type:

                variant_type: str = funcs.popleft()(ref, alt)

        # IndexError is raised when all functions have been popped from the deque
        # and none of them returned a value
        # In this case, the variant type is set to CSV (complex structural variant)
        # This is a fallback to avoid infinite loop
        except IndexError:

            logger.warning(f"Could not determine variant type from {ref}:{alt}")

            variant_type: str = "CSV"

        return variant_type
    
    @staticmethod
    @lru_cache(maxsize=7)
    def classification(arr: float, thresholds: tuple[float]) -> str:
        """
        categorize variant based on ALT Read Count Ratio (ARR) thresholds:
        Scale:
        L/P-HO Likely/Probably HOmozygous
        L/P-HE Likely/Probably HEterozygous
        L/P-SC Likely/Probably SubClonal
        [ LSC [ PSC [ PHE [ LHE ] PHE ] PHO ] LHO ]
        |     |     |     |     |     |     |     |
        0    [0]   [1]   [2]   [3]   [4]   [5]  100 (ARR)
        """

        MAX_THRESHOLD = 100.0
        MIN_THRESHOLD = 0.0

        classes: list[str] = ["LSC", "PSC", "PHE", "LHE", "PHE", "PHO", "LHO", "HO"]

        if arr == MIN_THRESHOLD:

            level = classes[0]

        elif arr == MAX_THRESHOLD:

            level = classes[-1]

        else:

            queue: deque = deque(zip(thresholds, classes))

            try:

                while arr >= queue[0][0]:

                    level: str = queue[0][1]

                    queue.popleft()
                
            except IndexError:

                logger.warning(f"Could not determine variant class from ARR with value: {arr}")

                level: str = "UKN"
        
        return level
    
    @staticmethod
    def get_genotype(genotypes: list[str], arr: float, thresholds: list[float]) -> tuple[str]:

        thresholds = thresholds.copy()
        
        thresholds.insert(0, 0.0)

        thresholds.append(100.0)

        level: str = VariantsRepository.classification(arr, tuple(thresholds))

        if level in ["PHE", "LHE"]:

            genotype: str = "0/1"

        elif level in ["PHO", "LHO", "HO"]:

            genotype: str = "1/1"

        else:

            genotype: str = "0/0"

        # Float division
        onset_threshold = math.ceil(len(genotypes) / 2)

        if sum(list(map(lambda x: x >= onset_threshold, dict(Counter(genotypes)).values()))) == 0:

            level: str = "WAR"

        return (genotype, level)
    
    @staticmethod
    def compute_background_metrics(pileup: Pileup, record: list[str], variant: dict, thresholds: list[float]) -> dict[str:float|str]:
        """
        estimate <BRC/R/E> (background read counts/ratio/enrichment)
        BRC : background read counts
        BRR : background read ratio
        BRE : background read enrichment
        Returns : a tuple containing the estimated BRC, BRR, and BRE.
        """

        thresholds.insert(0, 0.0)

        thresholds.insert(-1, 100.0)

        @lru_cache
        def classification(bre: float, arr: float, thresholds: tuple[float]) -> str:

            MAX_THRESHOLD = 100.0
            MIN_THRESHOLD = 0.0

            classes: list[str] = ["LCL", "PCL", "PNO", "LNO"]

            if arr >= 30.0:

                level: str = classes[1]

            else:

                if bre == MIN_THRESHOLD:

                    level = classes[0]

                elif bre == MAX_THRESHOLD:

                    level = classes[-1]

                else:

                    queue: deque = deque(zip(thresholds, classes))

                    try:

                        while bre >= queue[0][0]:

                            level: str = queue[0][1]

                            queue.popleft()
                        
                    except IndexError:

                        logger.warning(f"Could not categorize background signal from BRE with value: {bre}")

                        level: str = "UKN"
            
            return level
        
        pass

    def compute_sample_metrics():

        pass

    def get_filter():

        pass

    def populate(self, vcfs: dict):

        # Loop through each VCF file
        for caller in vcfs:

            # Get the header of the VCF file
            # The header contains the names of the columns in the VCF file
            # It is used to map the values in the VCF file to the correct keys in the dictionary
            header: dict = vcfs[caller]["vcf"].get_header()

            with open(vcfs[caller]["vcf"].get_path(), mode='r') as vcf:

                logger.debug(f"Processing {vcfs[caller]["vcf"].get_path()}")

                for line in vcf:

                    # VCF line structure is like [CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE]
                    record = line.strip().split('\t')

                    # Skip header
                    if line[0] == '#':
                        continue

                    # In haplotype caller skip some weird exceptions
                    # If sample part is empty
                    # If AD:DP is not in FORMAT (variant depth and total depth info)
                    # If total depth is 0
                    if caller == 'HC' and (
                        len(record[9])==0 or
                        not 'AD:DP' in record[8]
                        or record[9].split(':')[2] == '0'
                    ):
                        continue

                    # In pindel skip exceptions
                    # If it is writen INV instead of variant allele
                    # If coverage information is 0
                    if caller == 'PL' and (
                        int(record[9].split(':')[1].split(',')[0])==0 or
                        'INV' in record[4]
                    ):
                        continue

                    # Variant identifier is a string that contains the reference and alternative alleles
                    # It is used to identify the variant in the dictionary
                    ref: str = record[3]
                    alt: str = record[4]

                    if re.match(r"^[ATCG]+$", alt):

                        mutation: str = f"{ref}:{alt}"

                        category: str = VariantsRepository.get_variant_type(ref=ref, alt=alt)
                        
                        # Remove 'chr' from chromosome name if present
                        # This is to avoid issues with different chromosome naming conventions
                        chromosome: str = record[0].removeprefix('chr')
                        
                        if not chromosome in self.variants:

                            self.variants[chromosome] = {}

                        # Set positiion to integer
                        # Integer reduces memory usage in dictionary as key
                        Positions = namedtuple("Positions", ["vcf_position", "pileup_position"])
                        positions = Positions(vcf_position=int(record[1]), pileup_position=((int(record[1]) + self.pileup.DEL_FIRST_NC) if category == "DEL" else int(record[1])))
                            
                        if not positions in self.variants[chromosome]:

                            self.variants[chromosome][positions] = {}

                        # New variant
                        if not mutation in self.variants[chromosome][positions]:

                            # Create a new entry in the dictionary for the variant
                            self.variants[chromosome][positions][mutation] = {"collection": {"REF": ref,
                                                                                             "ALT": alt,
                                                                                             'VAF': {},
                                                                                            'GT': {},
                                                                                            'FILTER': {},
                                                                                            'INFO': {},
                                                                                            'FORMAT': {},
                                                                                            'SAMPLE': {},
                                                                                            'TRC': {},
                                                                                            'TRC-': {},
                                                                                            'TRC+': {},
                                                                                            'ARC':{},
                                                                                            'ARC-': {},
                                                                                            'ARC+': {},
                                                                                            'RRC':{},
                                                                                            'RRC-': {},
                                                                                            'RRC+': {},}}
                            
                            # Create a variable which link to the same memory address as the variant in the global dictionary
                            # It reduce writing complexity and improve readability
                            variant: dict = self.variants[chromosome][positions][mutation]

                            # Define variant type
                            variant['type'] = category

                            display: str = f"{ref}:{alt}"

                            if variant.get("type",'') in ['INV','MNV','CSV']:

                                # Keep trace of the complex variant
                                # Set allow faster lookup in O(1)
                                # instead of O(n)
                                self.cache["complex"].add(f"{chromosome}:{positions}:{mutation}")                                    
                                
                            else:

                                # Update data to be consistent with the pileup
                                # Update <ALT> allele descriptor for <INS> and <DEL>
                                # (Insertion) A/AA = A; (Deletion) AA/A = A
                                if variant['type'] == 'INS':
                                    
                                    # Pileup display for variant
                                    display: str = f"{ref}:{alt[1:]}"

                                # reset <POS> if <DEL>
                                # delA = p+1 in read pile
                                elif variant['type'] == 'DEL':

                                    # Pileup display for variant
                                    display: str = f"{ref[1:]}:{alt}"

                                # Keep trace of the variant
                                # Set allow faster lookup in O(1)
                                # instead of O(n)
                                self.cache["common"].add(f"{chromosome}:{positions}:{mutation}")

                            variant["display"] = display

                        # Variant already present
                        else:

                            # Create a variable which link to the same memory address as the variant in the global dictionary
                            # It reduce writing complexity and improve readability
                            variant: dict = self.variants[chromosome][positions][mutation]

                        variant['collection']['FILTER'][caller] = record[header['FILTER']]
                        variant['collection']['INFO'][caller] = record[header['INFO']]
                        variant['collection']['FORMAT'][caller] = record[header["FORMAT"]]
                        variant['collection']['SAMPLE'][caller] = record[header["SAMPLE"]]

                        gt = record[header["SAMPLE"]].split(':')[0]

                        if caller == "VD":
                            # Save genotype to raise warning if all callers don't report same GT.
                            # First manage case when Vardict return 1/0 instead of 0/1
                            if gt == '1/0':
                                gt = '0/1'

                        variant['collection']['GT'][caller] = gt

                        variant['collection']['VAF'][caller] = 100 * vcfs[caller]["vcf"].VAF(record)

                        # Store TRC, ARC and RRC for each VC
                        variant['collection']['TRC'][caller] = vcfs[caller]["vcf"].depth(record)

                        if caller in ['VS','VD','BT']:

                            variant['collection']['ARC'][caller], variant['collection']['ARC+'][caller], variant['collection']['ARC-'][caller] = vcfs[caller]["vcf"].arc(record)

                            variant['collection']['RRC'][caller], variant['collection']['RRC+'][caller], variant['collection']['RRC-'][caller] = vcfs[caller]["vcf"].rrc(record)

                            variant['collection'][f"TRC+"][caller] = variant['collection']['ARC+'][caller] + variant['collection']['RRC+'][caller]

                            variant['collection'][f"TRC-"][caller] = variant['collection']['ARC-'][caller] + variant['collection']['RRC-'][caller]

                        else:

                            variant['collection']['ARC'][caller] = vcfs[caller]["vcf"].arc(record)[0]

                            variant['collection']['RRC'][caller] = (
                                variant['collection']['TRC'][caller] -
                                variant['collection']['ARC'][caller]
                            )                       

    def normalize(self, sample: str, pileup: Pileup, thresholds: list[float], length_indels: int, sbm: float, sbm_homozygous: float) -> tuple[dict]:

        def get_variants(variants, chromosome):

            return {positions.pileup_position for positions in variants[chromosome]}
        
        cache = functions.Cache(func=get_variants, max_size=1)

        ITD: set[str] = set()

        with open(pileup.get_path(), mode='r') as f:

            for n, line in enumerate(f, start=1):

                datas = line.strip('\n').split('\t')

                if datas[0] == sample:

                    chromosome_pileup: str = datas[pileup.HEADER['chromosome']].removeprefix('chr')
                    position_pileup: int = int(datas[pileup.HEADER['position']])
                    reference_pileup: str = datas[pileup.HEADER['reference']]

                    # Check if variant is reported in one of the VCF files
                    if (reference_pileup != 'N') and (chromosome_pileup in self.variants) and (position_pileup in cache.call(args = [self.variants, chromosome_pileup])):

                        matchs: list[tuple] = [positions for positions in self.variants[chromosome_pileup] if positions.pileup_position == position_pileup]

                        if matchs:

                            for positions in matchs:
                        
                                for mutation in self.variants[chromosome_pileup][positions]:
                                    
                                    # O(1) lookup
                                    if f"{chromosome_pileup}:{positions}:{mutation}" in self.cache["common"]:
                                        
                                        # Create a variable which link to the same memory address as the variant in the global dictionary
                                        # It reduce writing complexity and improve readability
                                        variant: dict = self.variants[chromosome_pileup][positions][mutation]

                                        # depth of coverage at <CHR:POS> = sum(Nt) + #DEL (if any)
                                        coverage = {"strand": {"+": 0,
                                                            "-": 0},
                                                    "total": 0}

                                        for column, value in enumerate(datas[pileup.HEADER["A+"]:pileup.HEADER["N"]], start=0):

                                            try:
                                                coverage['total'] += int(value)
                                                if column in pileup.PLUS_STRAND:
                                                    coverage["strand"]['+'] += int(value)
                                                elif column in pileup.MINUS_STRAND:
                                                    coverage["strand"]['-'] += int(value)
                                            except ValueError:
                                                logger.warning("Uknown coverage value present in pileup file.")
                                                logger.warning(f"Warning was raised by: {value} at line {n} column {column}.")

                                        # manage DEL counts
                                        if datas[pileup.HEADER["DEL"]] != 'None':
                                            if datas[pileup.HEADER["DEL"]][0] == '*':
                                                try:
                                                    coverage["total"] += int(datas[pileup.HEADER["DEL"]].split(':')[1].split(';')[0])
                                                except ValueError:
                                                    logger.warning("Unknow DEL value present in pileup file.")
                                                    logger.warning(f"Warning was raised by: {datas[pileup.HEADER["DEL"]]} at line {n} column {pileup.HEADER["DEL"]}.")
                                            else:
                                                # clintools bug where a DEL does not start w/ *:\d+
                                                # (causing illegal division by zero)
                                                for deletion in datas[pileup.HEADER["DEL"]].split(';'):
                                                    del_cov1, del_cov2 = deletion.split(':')[1].split(',')
                                                    coverage["total"] += (int(del_cov1) + int(del_cov2))

                                        if not 'sample' in variant:
                                            variant['sample'] = {}

                                        variant['sample']['TRC'] = coverage.get("total", 0)
                                        variant['sample']['TRC-'] = coverage["strand"].get('-', 0)
                                        variant['sample']['TRC+'] = coverage["strand"].get('+', 0)

                                        # Add <REF> strand-specific counts for each <ALT> at <CHR:POS>
                                        # ie. the matched key in callset
                                        # RRC : Reference Read Counts
                                        # ARC : Alternative Read Counts
                                        for strand in coverage["strand"]:

                                            variant['sample'][f"RRC{strand}"] = datas[pileup.HEADER[f"{datas[pileup.HEADER['reference']]}{strand}"]]
                                        
                                        variant_count: int = 0

                                        # if <ALT> is an <INDEL>
                                        if variant['type'] in ['DEL','INS']:

                                            pattern: str = r"\b" + variant["display"].split(':')[0 if variant['type'] == 'DEL' else 1] + r"\b:[0-9]+,[0-9]+"

                                            data: str = datas[pileup.HEADER[variant['type']]]

                                            rsearch = re.search(pattern, data)
                                                
                                            if rsearch:
                                            
                                                    # Extract once and split to get counts
                                                    counts = rsearch.group(0).split(':')[-1]

                                                    arc_plus_strand, arc_minus_strand = map(int, counts.split(','))
                                                    
                                                    variant_count: int = arc_plus_strand + arc_minus_strand
                                                    
                                                    #Save indel count per strand in callset ARC+ and ARC-
                                                    for strand in coverage["strand"]:

                                                        arc = f'ARC{strand}'

                                                        variant['sample'].setdefault(arc, '-1')

                                                        variant['sample'][arc] = arc_plus_strand if strand == '+' else arc_minus_strand

                                        # any other <VarType>
                                        else:
                                            # keep 1st character of <ALT> string (approx. counts for non-SNV variants)

                                            for strand in coverage["strand"]:

                                                arc = f'ARC{strand}'

                                                if not arc in variant['sample']:
                                                    variant['sample'][arc] = ''

                                                variant['sample'][arc] = \
                                                    datas[pileup.HEADER[f"{mutation.split(':')[1][0]}{strand}"]]
                                                    
                                                variant_count += int(
                                                    datas[pileup.HEADER[f"{mutation.split(':')[1][0]}{strand}"]]
                                                )
                                        
                                        # --------------------------------------------------------------
                                        # <ALT> not covered in read pile;
                                        # keep trace of rejected variants (for test / rescuing purpose)
                                        # --------------------------------------------------------------
                                        if not variant_count:

                                            if (variant['type'] == 'INS') and (len(variant["collection"]["ALT"])-1 > length_indels):

                                                ITD.add(f"{chromosome_pileup}:{positions}:{mutation}")

                                            else:

                                                variant["filter"] = "REJECTED"
                                                
                                                if self.rescue:
                                                    # Cache variant for saving
                                                    self.cache["rejected"].add(f"{chromosome_pileup}:{positions}:{mutation}")

                                        else:

                                            # --------------------------------------------------------------
                                            # Compute metrics
                                            # --------------------------------------------------------------

                                            # Save number of Variant Callers that found this variant
                                            variant['sample']['VCN'] = len(variant['collection']['VAF'])

                                            # keep trace of used VC identifier(s)
                                            variant['sample']['VCI'] = ','.join(
                                                sorted(variant['collection']['VAF'].keys())
                                            )

                                            #format [R/A]RC for vcf
                                            variant['sample']['ARC'],\
                                            variant['sample']['RRC'] = \
                                                functions.format_rrc_arc(variant)

                                            # Compute VAF from pileup
                                            variant['sample']['ARR'] = format(
                                                (
                                                    (float(variant['sample']['ARC+']) + float(variant['sample']['ARC-'])) / float(variant['sample']['TRC'])
                                                )*100,
                                                '.5f'
                                            )

                                            # Estimate GT
                                            variant['sample']['GT'], variant['sample']['VAR'] = VariantsRepository.get_genotype(genotypes=variant["collection"]["GT"].values(),
                                                                                                                                arr=float(variant["sample"]["ARR"]),
                                                                                                                                thresholds=thresholds[0:6])

                                            variant['sample']['BRC'],\
                                            variant['sample']['BRR'],\
                                            variant['sample']['BRE'] = \
                                                functions.estimate_brc_r_e(variant, datas)

                                            variant['sample']['SBP'],variant['SBM'] = \
                                                functions.estimate_sbm(variant, sbm_homozygous)
                                            
                                            variant['sample']['LOW'] = int(float(variant["sample"]["ARR"]) < thresholds[-1])
                                            variant['sample']['BKG'] = \
                                                functions.categorize_background_signal(variant, thresholds)

                                            variant['sample']['PIL'] = 'Y'
                                            variant['sample']['RES'] = 'N'

                                            variant['filter'] = "FAIL" if ((variant["sample"] in ["PNO", "LNO"]) 
                                                                        or (variant['sample']['LOW'] == 1) 
                                                                        or ((abs(float(variant['sample']['SBP'])) <= 0.05) and (float(variant['sample']['SBM']) >= sbm))) else "PASS"
        
        return (self.variants, ITD)