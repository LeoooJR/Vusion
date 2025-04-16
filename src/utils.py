#!/bin/python3

from collections import namedtuple
import errors
from loguru import logger
import os
import re
import numpy
from scipy.stats import fisher_exact

class Cache():

    def __init__(self, func, max_size: int = 1):
        
        self.max_size = max_size

        self.cache: dict = {}

        self.func: function = func

    def add():

        pass

    def call(self, args: list[str]):

        if args[1].__hash__:

            if not args[1] in self.cache:

                if self.max_size and len(self.cache) >= self.max_size:

                    self.cache.clear()

                self.cache[args[1]] = self.func(args[0], args[1])
            
            return self.cache[args[1]]

# ===========================================================================================
# Basics functions on dictionary
# ===========================================================================================
def merge_collections(collections: list[tuple|list|dict]) -> tuple|list|dict:
    """
    Merges specified variant dictionaries.
    Parameters:
        dicts: A list of variant dictionaries
    Returns:
        A merged variant dictionary
    """
    if isinstance(collections[0], dict):

        if all(list(map(lambda collection: isinstance(collection, dict), collections))):

            output: dict = {}

            for collection in collections:
                
                output.update(collection)
        
        else:

            raise ValueError(f"Not all of the collections being merged are of the same data type.")

    else:

        raise ValueError(f"{type(collections[0])} cannot be merged.")

    return output


# ===========================================================================================
# Functions on variants
# ===========================================================================================

def estimate_sbm(variant,sbm_homozygous):
    """
    Estimate strand-bias metrics (<SBM>) and
    strand-bias Fisher exact test p-value (<SBP>) from TS(*) used by SLS/LRB.

    Parameters:
        dic (dict): A dictionary containing variant information, including final metrics.
        variant_key (str): The key representing the specific variant in the dictionary.
    Returns:
        tuple: A tuple containing
        the estimated strand-bias Fisher exact test p-value (<SBP>)
        and strand-bias metrics (<SBM>).

    (*) Modified formula is applied to avoid potential illegal division by zero.
    """
    alt_plus_x_ref_minus = int(variant['sample']['ARC+']) * \
                           int(variant['sample']['RRC-'])
    alt_minus_x_ref_plus = int(variant['sample']['ARC-']) * \
                           int(variant['sample']['RRC+'])

    # Check if there is a zero count or if variant is homozygous
    if (
        alt_plus_x_ref_minus==0 or
        alt_minus_x_ref_plus==0 or
        variant['sample']['GT'] == "1/1"
    ):
        strand_max = max(
            int(variant['sample']['ARC+']),
            int(variant['sample']['ARC-'])
        )
        if strand_max !=0:
            strand_ratio = strand_max/(
                int(variant['sample']['ARC+'])+
                int(variant['sample']['ARC-'])
            ) #ratio of reads on one strand to the other
        else:
            strand_ratio = 0 # no ALT allele = no stand bias

        variant['sample']['SBP'] = 0.05
        if strand_ratio >= sbm_homozygous:
            variant['sample']['SBM'] = '1.50000'
        else:
            variant['sample']['SBM'] = '0.50000'

    # Otherwise, perform standard Fisher Test
    else:
        m = max(alt_plus_x_ref_minus,alt_minus_x_ref_plus)
        # no <REF> allele = no strand-bias (case of illegal division by zero)
        # if m == 0:
        #     variant['sample']['SBM'] = '0.50000'
        # else:
            # variant['sample']['SBM'] = format(
            #     (m / (alt_plus_x_ref_minus + alt_minus_x_ref_plus)),
            #     '.5f'
            # )

        # [n11,n12] n1p
        # [n21,n22] n2p
        # np1 np2  npp
        n11 = int(variant['sample']['ARC+'])
        n12 = int(variant['sample']['ARC-'])

        if (variant['type'] == 'INS') or (variant['type'] == 'INV'):
            n21 = (
                int(variant['sample']['TRC+']) -
                int(variant['sample']['ARC+'])
            )
            n22 = (
                int(variant['sample']['TRC-']) -
                int(variant['sample']['ARC-'])
            )
            # 131223 : cas rare avec bug de CT ou TRC+<ARC+...
            if n21<0 :
                n21=0
            if n22<0:
                n22=0
        else:
            n21 = int(variant['sample']['RRC+'])
            n22 = int(variant['sample']['RRC-'])


        variant['sample']['SBM'] = format(
            (m / (alt_plus_x_ref_minus + alt_minus_x_ref_plus)),
            '.5f'
        )
        table_for_fisher_test = numpy.array([[n11, n12], [n21, n22]])
        oddsr, p = fisher_exact(table_for_fisher_test, alternative='two-sided')
        #variant['sample']['SBP'] = format(p, '.5f')
        variant['sample']['SBP'] = round(p,5)

    return(variant['sample']['SBP'],variant['sample']['SBM'])


def categorize_variant_type(variant,thresholds):
    """
    categorize VarType based on ALT Read Count Ratio (ARR) thresholds:

    Parameters:
    - dic (dict): A dictionary containing variant information,
                  including the ALT Read Count Ratio (ARR).
    - variant_key (str): The key representing the specific variant in the dictionary.
    - thresholds (list): A list of thresholds passed in mandatory options

    Returns:
    str: Returns the categorized VarType, which could be one of the following

    Note:

    L/P-HO Likely/Probably HOmozygous
    L/P-HE Likely/Probably HEterozygous
    L/P-SC Likely/Probably SubClonal
    [ LSC [ PSC [ PHE [ LHE ] PHE ] PHO ] LHO ]
    |     |     |     |     |     |     |     |
    0    [0]   [1]   [2]   [3]   [4]   [5]  100 (ARR)
    """
    alt_read_count_ratio = float(variant['sample']['ARR'])

    if alt_read_count_ratio > thresholds[5]:
        variant['sample']['VAR'] = 'LHO'
    elif alt_read_count_ratio > thresholds[4]:
        variant['sample']['VAR'] = 'PHO'
    elif alt_read_count_ratio < thresholds[0]:
        variant['sample']['VAR'] = 'LSC'
    elif alt_read_count_ratio < thresholds[1]:
        variant['sample']['VAR'] = 'PSC'
    elif (alt_read_count_ratio < thresholds[2]) or (alt_read_count_ratio > thresholds[3]):
        variant['sample']['VAR'] = 'PHE'
    else:
        variant['sample']['VAR'] = 'LHE'
    return variant['sample']['VAR']


def estimate_gt(variant):
    """
    Estimate GT from results of categorize VarType function
    GT is determined based on Variant Allele Frequency (VAF) values,
    either from pileup VAF or meanVAF if pileup is not available.

    Parameters:
    - dic (dict): A dictionary containing variant information, including VAF values.
    - variant_key (str): The key representing the specific variant in the dictionary.

    Returns:
    - str: The estimated genotype (GT) for the given variant based on VAF values.
    """
    if len(re.findall('SC',variant['sample']['VAR'])) > 0:
        variant['sample']['GT'] = '0/0'
    elif len(re.findall('HE',variant['sample']['VAR'])) > 0:
        variant['sample']['GT'] = '0/1'
    else:
        variant['sample']['GT'] = '1/1'
    return variant['sample']['GT']


def compare_gt(variant):
    """
    Compare GT identified between VC with pileup or meanVAF (if processed without pileup)
    If 50% or more of VC's GT are discordant to pileup or from meanVAF ,
    change VAR to WAR (warning), GT stay unchanged

    Parameters:
    - dic (dict): A dictionary containing variant information,
                  including GT values from different Variant callers.
    - variant_key (str): The key representing the specific variant in the dictionary.

    Returns:
    The (modified) VAR field in the dictionnary
    """
    callers = list(variant['collection']['VAF'].keys())
    pileup_gt = variant['sample']['GT']

    if 'FL' in callers:
        # If FL in callers, callers cant be agreed
        # Because FL do not compute GT
        variant['sample']['VAR'] = 'WAR'
    else:

        # gt_list = variant['collection']['GT'].copy()
        # gt_list[:] = [x if x != './1' else '0/1' for x in gt_list]
        # => TypeError: unhashable type: 'slice' => normal, jai change la structure du dic pour les GT
        gt_list = [] #List of all GT identified for a variant
        for variant_caller in callers:

            gt_list.append(variant['collection']['GT'][variant_caller])

        # Replace "./1" by "0/1", "1/." by "0/1", and "0/0" by "0/1".
        gt_list[:] = [x if x != './1' else '0/1' for x in gt_list]
        gt_list[:] = [x if x != '1/.' else '0/1' for x in gt_list]
        gt_list[:] = [x if x != '0/0' else '0/1' for x in gt_list]

        for genotype in gt_list:
            ok_gt = ['1/1', '0/1']
            if genotype not in ok_gt:
                #print("!! WARNING !!" + variant_key + " has a GT looking like that " + genotype + "\n")
                pass

        count_genotype = [gt_list.count('1/1'), gt_list.count('0/1')]
        if pileup_gt == '0/0':
            pileup_gt = '0/1'

        if pileup_gt == '1/1':
            gt_pileup_count = count_genotype[0]
        elif pileup_gt == '0/1':
            gt_pileup_count = count_genotype[1]


        genotype_ok = True
        # If 50% of variant caller GT are discordant to pileup,
        # set genotype_ok to False to return a warning
        #if (sum(count_genotype)-gt_pileup_count) / sum(count_genotype) >= 0.5:
        #    genotype_ok = False
        if sum(count_genotype) != gt_pileup_count:
            genotype_ok = False
        if not genotype_ok:
            variant['sample']['VAR'] = 'WAR'


    return variant['sample']['VAR']


def vaf_user_threshold(variant, thresholds):
    """
    Check if the variant allele ratio (VAF) is below a user-defined threshold.

     Parameters:
    - dic (dict): A dictionary containing variant information,
                  including the variant allele ratio (VAF).
    - variant_key (str): The key representing the specific variant in the dictionary.
    - thresholds (list): A list of thresholds passed in mandatory options

    Returns:
    int: Returns 0 if the VAF is equal to or above the user-defined threshold, otherwise returns 1

     Note:
    - The function updates the 'LOW' field in the dictionary to 0 or 1 based on the VAF threshold.
    - The VAF value is expected to be a floating-point number in the 'ARR' field of the dictionary.
    - The user-defined threshold is provided as a list,
    and thresholds[-1] is interpreted as a percentage.
    """

    return int(float(variant['sample']['ARR']) < thresholds[-1])


def format_rrc_arc(variant):
    """
    Format read count information [R/A]RC for VCF (Variant Call Format),
    i.e., generate [R/A]RC=[R/A]RC+,[R/A]RC-,sum().

    Parameters:
    - dic (dict): A dictionary containing variant information,
                  including read count details (RRC, ARC)
    - variant_key (str): The key representing the specific variant in the dictionary.

    Returns:
    Tuple: Returns a tuple containing the formatted [R/A]RC information for
           ARC (ALT Read Count) and RRC (REF Read Count).

    Note:
    - The input dictionary is expected to have 'sample' containing
    'ARC+' and 'ARC-' for ALT Read Count, and 'RRC+' and 'RRC-' for REF Read Count.
    """
    for read_count_type in ['RRC', 'ARC'] :
        # REF + ALT
        read_count_minus = variant['sample'][read_count_type + '-']
        read_count_plus = variant['sample'][read_count_type + '+']

        if read_count_minus!="-1" and read_count_plus != "-1":

            tmp_read_count = int(read_count_minus) + int(read_count_plus)

            variant['sample'][read_count_type] = ','.join(
            [
                str(int(read_count_plus)),
                str(int(read_count_minus)),
                str(tmp_read_count)
            ])
        else:
            
            variant['sample'][read_count_type] = ','.join(
            [
                str(read_count_plus),
                str(read_count_minus),
                str(variant['sample'][read_count_type])
            ])

    return(variant['sample']['ARC'],variant['sample']['RRC'])


def estimate_brc_r_e(variant,pileup_line_info):
    """
    estimate <BRC/R/E> (background read counts/ratio/enrichment)
    BRC : background read counts
    BRR : background read ratio
    BRE : background read enrichment

    Parameters:
    - dic (dict): A dictionary containing variant information,
                  including the ALT Read Count Ratio (ARR).
    - variant_key (str): The key representing the specific variant in the dictionary.
    - pileup_line_info (list): List containing information from pileup

    Returns : a tuple containing the estimated BRC, BRR, and BRE.
    """
    ref_and_alt_read_count = int(variant['sample']['ARC-']) + \
                             int(variant['sample']['ARC+']) + \
                             int(variant['sample']['RRC-']) + \
                             int(variant['sample']['RRC+'])
    ins_read_counts = 0
    total_read_count = variant['sample']['TRC']

    if pileup_line_info[14] != 'None':
        # Get number of read with ins format is A:1,0
        tmp_table_count = re.findall(r'\d+', pileup_line_info[14])
        for tmp_count in tmp_table_count:
            ins_read_counts += int(tmp_count)
    if variant['type'] == 'INS':
        ins_read_counts -= (
            int(variant['sample']['ARC+']) +
            int(variant['sample']['ARC-'])
        )

    # Removing DEL counts if variant is at some position of a DEL.
    # Because we miss valid variants in specific case like that
    # Clintool bug where non-existing deletion is reported and start with A , C, T or G
    if (variant['type'] != 'DEL') and (pileup_line_info[15] != 'None'):
        # Escaping clintool bug where non-existing deletion is reported and start with A , C, T or G
        if pileup_line_info[15][0] == '*':
            del_read_count = int(pileup_line_info[15].strip().split(';')[0].split(':')[1])
            tmp_total_read_count = int(total_read_count) - int(del_read_count)


            variant['sample']['BRC'] = (
                tmp_total_read_count +
                int(ins_read_counts) -
                min([
                    tmp_total_read_count,
                    int(ref_and_alt_read_count)
                    ])
                )

        else:
            variant['sample']['BRC'] = (
                int(total_read_count) +
                int(ins_read_counts) -
                min([
                    int(total_read_count),
                    int(ref_and_alt_read_count)
                ])
            )
    elif variant['type'] != 'DEL':
        variant['sample']['BRC'] = (
            int(total_read_count) +
            int(ins_read_counts) -
            min([
                int(total_read_count),
                int(ref_and_alt_read_count)
            ])
        )
    else:
        # 230413 Not counting snp background for Deletion
        del_read_count = pileup_line_info[15].strip().split(';')
        del_alt_count = 0
        for del_info in del_read_count:
            tmp_del_info = del_info.split(':')
            if (tmp_del_info[0] != '*') and (tmp_del_info[0] != (variant["collection"]["REF"])[1:]):
                del_alt_count += (
                    int(tmp_del_info[1].split(',')[0]) +
                    int(tmp_del_info[1].split(',')[1])
                )
        variant['sample']['BRC'] = del_alt_count


    background_read_counts = int(variant['sample']['BRC'])
    variant['sample']['BRR'] = format(
        background_read_counts / int(total_read_count),
        '.5f'
    )
    background_read_ratio = float(variant['sample']['BRR'])
    alt_read_count_ratio = float(variant['sample']['ARR'])/100
    variant['sample']['BRE'] = background_read_ratio / \
                                                 (alt_read_count_ratio + background_read_ratio)
    # scale ratio from [0-1] to [0-100]
    variant['sample']['BRE'] = format(
        float(variant['sample']['BRE']) * 100,
        '.5f'
    )
    variant['sample']['BRR'] = format(
        float(variant['sample']['BRR']) * 100,
        '.5f'
    )


    return(
        variant['sample']['BRC'],
        variant['sample']['BRR'],
        variant['sample']['BRE']
    )


def categorize_background_signal(variant,thresholds):
    """
    categorize background based on background read enrichment (BRE) thresholds :

    Parameters:
    - dic (dict): A dictionary containing variant information,
                  including background read enrichment (BRE).
    - variant_key (str): The key representing the specific variant in the dictionary.
    - thresholds (list): A list of thresholds passed in mandatory options

    Returns:
    str: Returns the categorized background signal,
         which could be one of the following: 'LNO', 'PNO', 'PCL', or 'LCL'.

    Note:
    - The input dictionary is expected to have 'sample' containing 'ARR' and 'BRE'.
    - update : 230413 If variant ratio is higher that 30% we consider it automatically clean

    L/P-NO Likely/Probably NOisy
    L/P-CL Likely/Probably CLean
    [ LCL [ PCL [ PNO [          LNO          ]
    |     |     |     |                       |
    0    [6]   [7]   [8]                    100 (BRE)
    """
    alt_read_count_ratio = float(variant['sample']['ARR'])
    if alt_read_count_ratio >= 30:
        variant['sample']['BKG'] = 'PCL'
    else:
        background_read_enrichment = float(variant['sample']['BRE'])

        if background_read_enrichment >= thresholds[8]:
            variant['sample']['BKG'] = 'LNO'
        elif background_read_enrichment >= thresholds[7]:
            variant['sample']['BKG'] = 'PNO'
        elif background_read_enrichment >= thresholds[6]:
            variant['sample']['BKG'] = 'PCL'
        else:
            variant['sample']['BKG'] = 'LCL'

    return variant['sample']['BKG']


def format_float_descriptors(variant,sbm):
    """
    format float descriptors and set <CT>-based <FILTER> field

    Parameters:
    - dic (dict): A dictionary containing variant information
    - variant_key (str): The key representing the specific variant in the dictionary.
    - sbm (float): Signal Bias Metric threshold (0.95 or 2 if disable_strand_bias)

    Returns:
    str: Returns the formatted <CT>-based <FILTER> field, which could be 'PASS' or 'FAIL'.
    """
    #for descriptor in ['ARR', 'BRR', 'BRE', 'SBM', 'SBP']:
    #   des = callsets[variant_key]['sample'][descriptor]
    #   callsets[variant_key]['sample'][descriptor] = des
    if len(re.findall('NO', variant['sample']['BKG'])) > 0:
        variant['sample']['FILTER'] = 'FAIL'
    elif variant['sample']['LOW'] == 1:
        variant['sample']['FILTER'] = 'FAIL'
    elif (abs(float(variant['sample']['SBP'])) <= 0.05) and \
        (float(variant['sample']['SBM']) >= sbm):
        variant['sample']['FILTER'] = 'FAIL'
    else:
        variant['sample']['FILTER'] = 'PASS'
    return variant['sample']['FILTER']


def process_without_pileup(variants: dict, lookups: set, thresholds: list[float], sbm, sbm_homozygous):
    """
    Process variant information without pileup data.

    Parameters:
    - dico (dict): A dictionary containing variant information from different Variant callers.
    - thresholds (list): A list of thresholds passed in mandatory options
    - sbm (float): Signal Bias Metric threshold.

    Returns:
    dict: Returns the processed dictionary containing updated final metrics for each variant.

    Note:
    - The function calculates mean values for ARR, TRC, ARC, and RRC from Variant callers
    - It determines VAR and GT based on mean VAF and GT agreement among Variant callers
    - Certain values (BRC, BRR, BRE, BKG) are set to -1 as
    they cannot be determined without pileup data
    - Depending of VC, other values related to them are also set to -1
    (ARC+,ARC-,RRC+,RRC-,SBM,SBP,TRC+,TRC-)
    """
    for lookup in lookups:

        Positions = namedtuple("Position", ["vcf_position", "pileup_position"])

        chromsome, positions, ref, alt = lookup.split(':')

        variant_identifier = f"{ref}:{alt}"

        variant: dict = variants[chromsome][eval(positions)][variant_identifier]

        if not 'sample' in variant:
            variant['sample'] = {}

        # Save number of Variant callers that found this variant
        variant['sample']['VCN'] = len(variant['collection']['VAF'])

        # keep trace of used VC identifier(s)
        variant['sample']['VCI'] = ','.join(
            sorted(
                variant['collection']['VAF'].keys()
            ))


        # Get mean vaf from Variant callers and convert in %
        vaf_sum = sum(variant['collection']['VAF'].values())
        vcn = float(variant['sample']['VCN'])
        mean_arr = vaf_sum / vcn
        variant['sample']['ARR'] = format(mean_arr,'.5f')

        #Filter on vaf
        variant['sample']['LOW'] = vaf_user_threshold(
            variant,
            thresholds
        )
        variant['sample']['VAR'] = categorize_variant_type(
            variant,
            thresholds
        )


        # Get mean TRC from Variant callers
        variant['sample']['TRC'] = int(round(
            sum(variant['collection']['TRC'].values()) /
            variant['sample']['VCN'],
            5))

        # Get mean ARC from Variant callers
        variant['sample']['ARC'] = int(round(
            (sum(variant['collection']['ARC'].values())/
            variant['sample']['VCN']),
            5))

        # Get mean RRC from Variant callers
        variant['sample']['RRC'] = int(round(
            (sum(variant['collection']['RRC'].values())/
            variant['sample']['VCN']),
            5))


        # Determine GT
        # estimate <GT> from CT data
        variant['sample']['GT'] = estimate_gt(variant)
        # 2. Set WAR if VC do not agree
        variant['sample']['VAR'] = compare_gt(variant)

        # Without pileup, impossible to determine certain values - set -1
        variant['sample']['BRC'] = "-1"
        variant['sample']['BRR'] = "-1"
        variant['sample']['BRE'] = "-1"
        variant['sample']['BKG'] = "-1"


        if 'RRC+' not in variant['collection'] or len(variant['collection']['RRC+']) != len(variant['collection']['RRC']):

            undertermined_ct_values = ["ARC+","ARC-","RRC+", "RRC-", "SBM", "SBP", "TRC+", "TRC-"]

            for ct_values in undertermined_ct_values:
                variant['sample'].setdefault(ct_values, "-1")

            variant['sample']['ARC'], variant['sample']['RRC'] = format_rrc_arc(variant)

        else:
            # Get mean TRC+/-,RRC+/-, ARC+/- from Variant callers
            for strand in["+","-"]:
                for count in ['TRC','RRC','ARC']:
                    if count == 'TRC':
                        sum_strand = sum(variant['collection']['RRC'+strand].values()) + sum(variant['collection']['ARC'+strand].values())
                    else:
                        sum_strand = sum(variant['collection'][count+strand].values())
                    variant['sample'][count+strand] = round(
                        sum_strand/
                        variant['sample']['VCN'],
                        5)

            variant['sample']['ARC'],\
            variant['sample']['RRC'] = format_rrc_arc(variant)

            variant['sample']['SBP'],\
            variant['SBM'] = estimate_sbm(variant, sbm_homozygous)

        if variant.get("filter", "") == "REJECTED":
            
            # rescue <PASS> calls only
            # do not rescue SNV (more likely to be VS artefacts)
            # pindel bug where ARC > TRC
            filter: str = "PASS" if ((format_float_descriptors(variant, sbm) == "PASS") and (variant["type"] != "SNV") and float(variant["sample"]["ARR"]) <= 100.0) else "REJECTED"

            if filter == "PASS":

                variant['sample']['RES'] = 'Y'

        else:

            filter: str = format_float_descriptors(variant, sbm)

            variant['sample']['RES'] = 'N'

        variant['filter'] = filter

        variant['sample']['PIL'] = "N"

    return variants
