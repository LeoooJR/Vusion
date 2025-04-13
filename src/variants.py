from collections import deque
from files import Pileup
from loguru import logger
import re
import utils as functions

class VariantsRepository():

    def __init__(self):

        self.variants: dict = {}

        self.trace: set[str] = {}

        self.rev: dict = {}
        
        self.INV_MNV_CSV: set[str] = {}

        self.FLiT3r: set[str] = {}

    @staticmethod
    def get_variant_type(ref: str, alt: str):

        variant_type: str = ''

        funcs: deque = deque([is_snp, is_ins, is_del, is_inv, is_mnv])

        def is_snp(ref: str, alt: str) -> str:

            return "SNV" if len(ref) == 1 and len(alt) == 1 else ''

        def is_ins(ref: str, alt: str) -> str:
            
            return "INS" if len(ref) == 1 and len(alt) > 1 else ''

        def is_del(ref: str, alt: str) -> str:

            return "DEL" if len(ref) > 1 and len(alt) == 1 else ''

        def is_inv(ref: str, alt: str) -> str:

            OLD_CHARS = "ACGTacgt"
            REPLACE_CHARS = "TGCAtgca"

            rev = alt.translate(str.maketrans(OLD_CHARS,REPLACE_CHARS))[::-1]

            return "INV" if len(ref) == len(alt) and ref == rev else ''

        def is_mnv(ref: str, alt: str) -> str:

            OLD_CHARS = "ACGTacgt"
            REPLACE_CHARS = "TGCAtgca"

            rev = alt.translate(str.maketrans(OLD_CHARS,REPLACE_CHARS))[::-1]

            return "INV" if len(ref) == len(alt) and ref != rev else ''
        
        try:

            while not variant_type:

                variant_type: str = funcs.popleft()(ref, alt)

        except IndexError:

            logger.warning(f"Could not determine variant type from {ref}:{alt}")

            variant_type: str = "CSV"

        return variant_type


    def populate(self, vcfs: dict):

        for caller in vcfs:

            with open(vcfs[caller]["vcf"].get_path(), mode='r') as vcf:

                logger.debug(f"Processing {vcfs[caller]["vcf"].get_path()}")

                for line in vcf:

                    # Reminder of what info will containe
                    # vcf line structure is like [CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE]
                    datas = line.strip().split('\t')

                    # Skip header
                    if line[0] == '#':
                        continue

                    # In haplotype caller skip some weird exceptions
                    # If sample part is empty
                    # If AD:DP is not in FORMAT (variant depth and total depth info)
                    # If total depth is 0
                    if caller == 'HC' and (
                        len(datas[9])==0 or
                        not 'AD:DP' in datas[8]
                        or datas[9].split(':')[2] == '0'
                    ):
                        continue

                    # In pindel skip exceptions
                    # If it is writen INV instead of variant allele
                    # If coverage information is 0
                    if caller == 'PL' and (
                        int(datas[9].split(':')[1].split(',')[0])==0 or
                        'INV' in datas[4]
                    ):
                        continue

                    chromosome: str = datas[0].removeprefix('chr')
                    
                    if not chromosome in self.variants:

                        self.variants[chromosome] = {}

                    position = int(datas[1])
                        
                    if not position in self.variants[chromosome]:

                        self.variants[chromosome][position] = {}

                    variant_identifier: str = f"{datas[3]}:{datas[4]}"

                    if not variant_identifier in self.variants[chromosome][position]:

                        self.variants[chromosome][position][variant_identifier] = {{"VC": {'REF': datas[3],
                                                                                        'ALT': datas[4],
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
                                                                                        'RRC+': {},}}}
                    
                        ## 1/ define <VT>
                        self.variants[chromosome][position][variant_identifier]['VT'] = functions.define_variant_type(ref=datas[3], alt=datas[4])

                    variant: dict = self.variants[chromosome][position][variant_identifier]

                    if caller == 'FL' : # Less informations with FLiT3r
                        variant['VC']['FORMAT'][caller] = datas[vcfs[input[0]]["vcf"].get_header()["FORMAT"]]
                        variant['VC']['SAMPLE'][caller] = datas[vcfs[input[0]]["vcf"].get_header()["SAMPLE"]]
                    else:
                        variant['VC']['FILTER'][caller] = datas[vcfs[input[0]]["vcf"].get_header()['FILTER']]
                        variant['VC']['INFO'][caller] = datas[vcfs[input[0]]["vcf"].get_header()['INFO']]
                        variant['VC']['FORMAT'][caller] = datas[vcfs[input[0]]["vcf"].get_header()["FORMAT"]]
                        variant['VC']['SAMPLE'][caller] = datas[vcfs[input[0]]["vcf"].get_header()["SAMPLE"]]

                        GT = datas[9].split(':')[0]

                        if caller == "VD":
                            # Save genotype to raise warning if all callers don't report same GT.
                            # First manage case when Vardict return 1/0 instead of 0/1
                            
                            if GT == '1/0':
                                GT = '0/1'

                        variant['GT'][caller] = GT

                    variant['VC']['VAF'][caller] = 100 * vcfs[caller]["vcf"].VAF(datas)

                    # Store TRC, ARC and RRC for each VC
                    variant['VC']['TRC'][caller] = vcfs[caller]["vcf"].depth(datas)

                    if caller in ['VS','VD','BT']:

                        variant['VC']['ARC'][caller], variant['VC']['ARC+'][caller], variant['VC']['ARC-'][caller] = vcfs[caller]["vcf"].arc(datas)

                        variant['VC']['RRC'][caller], variant['VC']['RRC+'][caller], variant['VC']['RRC-'][caller] = vcfs[caller]["vcf"].rrc(datas)

                        variant['VC'][f"TRC+"][caller] = variant['VC']['ARC+'][caller] + variant['VC']['RRC+'][caller]
                        variant['VC'][f"TRC-"][caller] = variant['VC']['ARC-'][caller] + variant['VC']['RRC-'][caller]

                    else:

                        variant['VC']['ARC'][caller] = vcfs[caller]["vcf"].arc(datas)[0]

                        variant['VC']['RRC'][caller] = (
                            variant['VC']['TRC'][caller] -
                            variant['VC']['ARC'][caller]
                        )

                    # TBI (Must be done only once)
                    # 2/ revert vt norm of variants that should not have been left-aligned (if any)
                    if variant['VT'] == 'DEL':

                        for caller in variant['VC']['VAF']:

                            ## from vt decompose
                            if 'OLD_MULTIALLELIC' in variant['VC']['INFO'][caller]:
                                continue

                            # from vt decompose
                            if 'OLD_VARIANT' in variant['VC']['INFO'][caller]:
                                info_old_variant = variant['VC']['INFO'][caller] \
                                                .split('OLD_VARIANT=')[1] \
                                                .split(',')
                                for old_variant in info_old_variant:
                                    ov_ref, ov_alt = old_variant.split('/')
                                    # DEL is not parsimonious eg. POS:ABB/AB (should be POS+1:BB/B)
                                    if len(ov_alt) != 1:
                                        continue
                                    # Keep trace of old/new DEL descriptor
                                    variant_identifier_updated: str = f"{ov_ref}:{ov_alt}"

                                    if not f"{chromosome}:{position}:{variant_identifier}" in self.rev:
                                        self.rev[f"{chromosome}:{position}:{variant_identifier}"] = []
                                        
                                    self.rev[f"{chromosome}:{position}:{variant_identifier}"].append(f"{chromosome}:{position}:{variant_identifier_updated}")

                                    # DEL is parsimonious;
                                    # clone <hash_key> according to <OLD_VARIANT> description
                                    # call[new_hash] = callsets[hash].copy()
                                    # Update variant
                                    value = self.variants[chromosome][position].pop(variant_identifier)
                                    value["VC"]["REF"] = ov_ref
                                    value["VC"]["ALT"] = ov_alt
                                    self.variants[chromosome][position][variant_identifier_updated] = value
                                    # Sync hash value
                                    variant_identifier = variant_identifier_updated
                                break

                    # ETBI

                    ref, alt = variant_identifier.split(':')

                    if not (variant["VT"] in ['INV','MNV','CSV'] or "FL" in variant["VC"]["VAF"]):

                        # reset <ALT> allele descriptor for <INS> and <DEL>
                        # A/AA = insA; AA/A = delA
                        if variant['VT'] == 'INS':

                            variant_identifier_updated: str = f"{ref}:{alt[1:]}"

                            self.variants[chromosome][position][variant_identifier_updated] = self.variants[chromosome][position].pop(variant_identifier) 

                        elif variant['VT'] == 'DEL':

                            variant_identifier_updated: str = f"{ref[1:]}:{alt}"

                            self.variants[chromosome][position][variant_identifier_updated] = self.variants[chromosome][position].pop(variant_identifier)

                            if len(self.variants[chromosome][position]) == 1:

                                self.variants[chromosome][(position+1)] = self.variants[chromosome].pop(position)

                            else:

                                if not (position+1) in self.variants[chromosome]:
                                    
                                    self.variants[chromosome][(position+1)] = {}

                                self.variants[chromosome][(position+1)][variant_identifier_updated] = self.variants[chromosome][position].pop(variant_identifier)
                            
                            # reset <POS> if <DEL>
                            # delA = p+1 in read pile
                            # position += 1

                            self.trace.add(f"{chromosome}:{position}:{ref}:{alt}")                                    
                        
                    elif variant["VT"] in ['INV','MNV','CSV']:

                        self.INV_MNV_CSV.add(f"{chromosome}:{position}:{ref}:{alt}")
                        
                    elif "FL" in variant["VC"]["VAF"]:

                        self.FLiT3r.add(f"{chromosome}:{position}:{ref}:{alt}")

    def normalize(self, sample: str, pileup: Pileup, thresholds: list[float], length_indels: int, sbm: float, sbm_homozygous: float) -> tuple[dict]:

        ITD: set[str] = {}
        rejected: dict = {}

        with open(pileup.get_path(), mode='r') as f:

            for n, line in enumerate(f, start=1):

                datas = line.strip('\n').split('\t')

                if datas[0] == sample:

                    chromosome: str = datas[1].removeprefix('chr')
                    position: int = int(datas[2])

                    # Check if variant is reported in one of the VCF files
                    if (datas[pileup.HEADER['reference']] != 'N') and (chromosome in self.variants) and (position in self.variants[chromosome]):
                        
                        for variant_identifier in self.variants[chromosome][position]:

                            ref, alt = variant_identifier.split(':')

                            if f"{chromosome}:{position}:{ref}:{alt}" in self.trace:

                                variant: dict = self.variants[chromosome][position][variant_identifier]

                                # depth of coverage at <CHR:POS> = sum(Nt) + #DEL (if any)
                                coverage = {"plus": 0,
                                            "minus": 0,
                                            "total": 0}

                                for column, value in enumerate(datas[5:13], start=0):

                                    try:
                                        coverage['total'] += int(value)
                                        if column in pileup.PLUS_STRAND:
                                            coverage['plus'] += int(value)
                                        elif column in pileup.MINUS_STRAND:
                                            coverage["minus"] += int(value)
                                    except ValueError:
                                        logger.warning("Uknown coverage value present in pileup file.")
                                        logger.warning(f"Warning was raised by: {value} at line {n} column {column}.")

                                if not 'final_metrics' in variant:
                                    variant['final_metrics'] = {}

                                # manage DEL counts
                                if datas[-1] != 'None':
                                    if datas[-1][0] == '*':
                                        try:
                                            coverage["total"] += int(datas[-1].split(':')[1].split(';')[0])
                                        except ValueError:
                                            logger.warning("Unknow DEL value present in pileup file.")
                                            logger.warning(f"Warning was raised by: {datas[-1]} at line {n} column {column}.")
                                    else:
                                        # clintools bug where a DEL does not start w/ *:\d+
                                        # (causing illegal division by zero)
                                        for deletion in datas[-1].split(';'):
                                            del_cov1, del_cov2 = deletion.split(':')[1].split(',')
                                            coverage["total"] += (int(del_cov1) + int(del_cov2))

                                variant['final_metrics']['TRC'] = coverage["total"]
                                variant['final_metrics']['TRC-'] = coverage['minus']
                                variant['final_metrics']['TRC+'] = coverage['plus']

                                # Add <REF> strand-specific counts for each <ALT> at <CHR:POS>
                                # ie. the matched key in callset
                                # RRC : Reference Read Counts
                                # ARC : Alternative Read Counts
                                for strand in ['-','+']:

                                    variant['final_metrics'][f"RRC{strand}"] = datas[pileup.HEADER[f"{datas[pileup.HEADER['reference']]}{strand}"]]

                                # if <ALT> is an <INDEL>
                                variants_count: int = 0

                                if (variant['VT'] == 'DEL') or (variant['VT'] == 'INS'):
                                        
                                    if re.search(r"\b" + variant_identifier.split(':')[1] + r"\b:[0-9]+,[0-9]+" ,\
                                        (datas[pileup.HEADER[variant['VT']]])):
                                    
                                            arc_plus_strand, arc_minus_strand = re.findall(
                                                r"\b" + variant_identifier.split(':')[1] + r"\b:[0-9]+,[0-9]+",
                                                datas[pileup.HEADER[variant['VT']]]
                                                )[0].split(':')[-1].split(',')
                                            
                                            variants_count: int = int(arc_plus_strand) + int(arc_minus_strand)
                                            
                                            #Save indel count per strand in callset ARC+ and ARC-
                                            for strand in ['+','-']:

                                                arc = f'ARC{strand}'

                                                if not arc in variant['final_metrics']:
                                                    variant['final_metrics'][arc] = ''

                                                variant['final_metrics'][arc] = arc_plus_strand if strand == '+' else arc_minus_strand

                                # any other <VarType>
                                else:
                                    # keep 1st character of <ALT> string (approx. counts for non-SNV variants)

                                    if variant_identifier.split(':')[1][0] == 'N':
                                        del variant
                                    else:
                                        for strand in ['+','-']:

                                            arc = f'ARC{strand}'

                                            if not arc in variant['final_metrics']:
                                                variant['final_metrics'][arc] = ''

                                            variant['final_metrics'][arc] = \
                                                datas[pileup.HEADER[f"{variant_identifier.split(':')[1][0]}{strand}"]]
                                            
                                            variants_count += int(
                                                datas[pileup.HEADER[f"{variant_identifier.split(':')[1][0]}{strand}"]]
                                            )
                                
                                # --------------------------------------------------------------
                                # <ALT> not covered in read pile;
                                # keep trace of rejected variants (for test / rescuing purpose)
                                # --------------------------------------------------------------
                                if variants_count == 0:

                                    if (variant['VT'] == 'INS') and (len(variant["VC"]["ALT"])-1 > length_indels):

                                        ITD.add(f"{chromosome}:{position}:{ref}:{alt}")

                                    else:
                                        
                                        rejected[f"{chromosome}:{position}:{ref}:{alt}"] = self.variants[chromosome][position].pop(variant_identifier)

                                    continue

                                # --------------------------------------------------------------
                                # Compute metrics
                                # --------------------------------------------------------------

                                # Save number of Variant Callers that found this variant
                                variant['final_metrics']['VCN'] = len(variant['VC']['VAF'])

                                # keep trace of used VC identifier(s)
                                variant['final_metrics']['VCI'] = ','.join(
                                    sorted(variant['VC']['VAF'].keys(), key=str.lower)
                                )

                                #format [R/A]RC for vcf
                                variant['final_metrics']['ARC'],\
                                variant['final_metrics']['RRC'] = \
                                    functions.format_rrc_arc(variant)

                                # Compute VAF from pileup
                                alt_read_count_plus = variant['final_metrics']['ARC+']
                                alt_read_count_minus = variant['final_metrics']['ARC-']
                                total_read_count = variant['final_metrics']['TRC']
                                variant['final_metrics']['ARR'] = format(
                                    (
                                        (float(alt_read_count_plus) + float(alt_read_count_minus)) /
                                        float(total_read_count)
                                    )*100,
                                    '.5f'
                                )

                                # Estimate GT
                                variant['final_metrics']['VAR'] = \
                                    functions.categorize_variant_type(variant, thresholds)
                                variant['final_metrics']['GT'] = \
                                    functions.estimate_gt(variant)
                                
                                # If variant callers do not agree on genotype change VAR to WAR (warning)
                                variant['final_metrics']['VAR'] = \
                                    functions.compare_gt(variant)

                                variant['final_metrics']['BRC'],\
                                variant['final_metrics']['BRR'],\
                                variant['final_metrics']['BRE'] = \
                                    functions.estimate_brc_r_e(variant, datas)

                                variant['final_metrics']['SBP'],variant['SBM'] = \
                                    functions.estimate_sbm(variant, sbm_homozygous)
                                variant['final_metrics']['LOW'] = \
                                    functions.vaf_user_threshold(variant, thresholds)
                                variant['final_metrics']['BKG'] = \
                                    functions.categorize_background_signal(variant, thresholds)
                                variant['final_metrics']['FILTER'] = \
                                    functions.format_float_descriptors(variant, sbm)

                                variant['final_metrics']['PIL'] = 'Y'
                                variant['final_metrics']['RES'] = 'N'
        
        return (self.variants, ITD, rejected)