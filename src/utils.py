#!/bin/python3

import errors
from loguru import logger
import os
import re
import numpy
from scipy.stats import fisher_exact

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
def define_variant_type(ref: str, alt: str):
    """
    define VarType (VT) ie. SNV, MNV, INS/DEL, INV or CSV (Complex Structural Variant)
    Parameters:
    - sequence: String
    - seq_len : a dictionnary containing the length of the reference and the alternative allele
    Return : the variant type
    """
    OLD_CHARS = "ACGTacgt"
    REPLACE_CHARS = "TGCAtgca"
    rev = alt.translate(str.maketrans(OLD_CHARS,REPLACE_CHARS))[::-1]
    if len(ref) == 1 and len(alt) == 1:
        type = 'SNV'
    elif len(ref) == 1 and len(alt) > 1 :
        type = 'INS'
    elif len(ref) > 1 and len(alt) == 1:
        type = 'DEL'
    elif len(ref) == len(alt):
        if ref == rev :
            type = 'INV'
        else:
            # mostly coming from <VD> and both annotated as <Complex>
            type = 'MNV'
    else:
        # mostly coming from <PL> and annotated as either <INV> or <RPL>
        type = 'CSV'
    return type


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
    alt_plus_x_ref_minus = int(variant['final_metrics']['ARC+']) * \
                           int(variant['final_metrics']['RRC-'])
    alt_minus_x_ref_plus = int(variant['final_metrics']['ARC-']) * \
                           int(variant['final_metrics']['RRC+'])

    # Check if there is a zero count or if variant is homozygous
    if (
        alt_plus_x_ref_minus==0 or
        alt_minus_x_ref_plus==0 or
        variant['final_metrics']['GT'] == "1/1"
    ):
        strand_max = max(
            int(variant['final_metrics']['ARC+']),
            int(variant['final_metrics']['ARC-'])
        )
        if strand_max !=0:
            strand_ratio = strand_max/(
                int(variant['final_metrics']['ARC+'])+
                int(variant['final_metrics']['ARC-'])
            ) #ratio of reads on one strand to the other
        else:
            strand_ratio = 0 # no ALT allele = no stand bias

        variant['final_metrics']['SBP'] = 0.05
        if strand_ratio >= sbm_homozygous:
            variant['final_metrics']['SBM'] = '1.50000'
        else:
            variant['final_metrics']['SBM'] = '0.50000'

    # Otherwise, perform standard Fisher Test
    else:
        m = max(alt_plus_x_ref_minus,alt_minus_x_ref_plus)
        # no <REF> allele = no strand-bias (case of illegal division by zero)
        # if m == 0:
        #     variant['final_metrics']['SBM'] = '0.50000'
        # else:
            # variant['final_metrics']['SBM'] = format(
            #     (m / (alt_plus_x_ref_minus + alt_minus_x_ref_plus)),
            #     '.5f'
            # )

        # [n11,n12] n1p
        # [n21,n22] n2p
        # np1 np2  npp
        n11 = int(variant['final_metrics']['ARC+'])
        n12 = int(variant['final_metrics']['ARC-'])

        if (variant['VT'] == 'INS') or (variant['VT'] == 'INV'):
            n21 = (
                int(variant['final_metrics']['TRC+']) -
                int(variant['final_metrics']['ARC+'])
            )
            n22 = (
                int(variant['final_metrics']['TRC-']) -
                int(variant['final_metrics']['ARC-'])
            )
            # 131223 : cas rare avec bug de CT ou TRC+<ARC+...
            if n21<0 :
                n21=0
            if n22<0:
                n22=0
        else:
            n21 = int(variant['final_metrics']['RRC+'])
            n22 = int(variant['final_metrics']['RRC-'])


        variant['final_metrics']['SBM'] = format(
            (m / (alt_plus_x_ref_minus + alt_minus_x_ref_plus)),
            '.5f'
        )
        table_for_fisher_test = numpy.array([[n11, n12], [n21, n22]])
        oddsr, p = fisher_exact(table_for_fisher_test, alternative='two-sided')
        #variant['final_metrics']['SBP'] = format(p, '.5f')
        variant['final_metrics']['SBP'] = round(p,5)

    return(variant['final_metrics']['SBP'],variant['final_metrics']['SBM'])


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
    alt_read_count_ratio = float(variant['final_metrics']['ARR'])

    if alt_read_count_ratio > thresholds[5]:
        variant['final_metrics']['VAR'] = 'LHO'
    elif alt_read_count_ratio > thresholds[4]:
        variant['final_metrics']['VAR'] = 'PHO'
    elif alt_read_count_ratio < thresholds[0]:
        variant['final_metrics']['VAR'] = 'LSC'
    elif alt_read_count_ratio < thresholds[1]:
        variant['final_metrics']['VAR'] = 'PSC'
    elif (alt_read_count_ratio < thresholds[2]) or (alt_read_count_ratio > thresholds[3]):
        variant['final_metrics']['VAR'] = 'PHE'
    else:
        variant['final_metrics']['VAR'] = 'LHE'
    return variant['final_metrics']['VAR']


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
    if len(re.findall('SC',variant['final_metrics']['VAR'])) > 0:
        variant['final_metrics']['GT'] = '0/0'
    elif len(re.findall('HE',variant['final_metrics']['VAR'])) > 0:
        variant['final_metrics']['GT'] = '0/1'
    else:
        variant['final_metrics']['GT'] = '1/1'
    return variant['final_metrics']['GT']


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
    callers = list(variant['VC']['VAF'].keys())
    pileup_gt = variant['final_metrics']['GT']

    if 'FL' in callers:
        # If FL in callers, callers cant be agreed
        # Because FL do not compute GT
        variant['final_metrics']['VAR'] = 'WAR'
    else:

        # gt_list = variant['VC']['GT'].copy()
        # gt_list[:] = [x if x != './1' else '0/1' for x in gt_list]
        # => TypeError: unhashable type: 'slice' => normal, jai change la structure du dic pour les GT
        gt_list = [] #List of all GT identified for a variant
        for variant_caller in callers:

            gt_list.append(variant['VC']['GT'][variant_caller])

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
            variant['final_metrics']['VAR'] = 'WAR'


    return variant['final_metrics']['VAR']


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

    return int(float(variant['final_metrics']['ARR']) < thresholds[-1])


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
    - The input dictionary is expected to have 'final_metrics' containing
    'ARC+' and 'ARC-' for ALT Read Count, and 'RRC+' and 'RRC-' for REF Read Count.
    """
    for read_count_type in ['RRC', 'ARC'] :
        # REF + ALT
        read_count_minus = variant['final_metrics'][read_count_type + '-']
        read_count_plus = variant['final_metrics'][read_count_type + '+']

        if read_count_minus!="-1" and read_count_plus != "-1":
            tmp_read_count = int(read_count_minus) + int(read_count_plus)
            variant['final_metrics'][read_count_type] = ','.join(
            [
                str(int(read_count_plus)),
                str(int(read_count_minus)),
                str(tmp_read_count)
            ])
        else:
            variant['final_metrics'][read_count_type] = ','.join(
            [
                str(read_count_plus),
                str(read_count_minus),
                str(variant['final_metrics'][read_count_type])
            ])

    return(variant['final_metrics']['ARC'],variant['final_metrics']['RRC'])


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
    ref_and_alt_read_count = int(variant['final_metrics']['ARC-']) + \
                             int(variant['final_metrics']['ARC+']) + \
                             int(variant['final_metrics']['RRC-']) + \
                             int(variant['final_metrics']['RRC+'])
    ins_read_counts = 0
    total_read_count = variant['final_metrics']['TRC']

    if pileup_line_info[14] != 'None':
        # Get number of read with ins format is A:1,0
        tmp_table_count = re.findall(r'\d+', pileup_line_info[14])
        for tmp_count in tmp_table_count:
            ins_read_counts += int(tmp_count)
    if variant['VT'] == 'INS':
        ins_read_counts -= (
            int(variant['final_metrics']['ARC+']) +
            int(variant['final_metrics']['ARC-'])
        )

    # Removing DEL counts if variant is at some position of a DEL.
    # Because we miss valid variants in specific case like that
    # Clintool bug where non-existing deletion is reported and start with A , C, T or G
    if (variant['VT'] != 'DEL') and (pileup_line_info[15] != 'None'):
        # Escaping clintool bug where non-existing deletion is reported and start with A , C, T or G
        if pileup_line_info[15][0] == '*':
            del_read_count = int(pileup_line_info[15].strip().split(';')[0].split(':')[1])
            tmp_total_read_count = int(total_read_count) - int(del_read_count)


            variant['final_metrics']['BRC'] = (
                tmp_total_read_count +
                int(ins_read_counts) -
                min([
                    tmp_total_read_count,
                    int(ref_and_alt_read_count)
                    ])
                )

        else:
            variant['final_metrics']['BRC'] = (
                int(total_read_count) +
                int(ins_read_counts) -
                min([
                    int(total_read_count),
                    int(ref_and_alt_read_count)
                ])
            )
    elif variant['VT'] != 'DEL':
        variant['final_metrics']['BRC'] = (
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
            if (tmp_del_info[0] != '*') and (tmp_del_info[0] != (variant["VC"]["REF"])[1:]):
                del_alt_count += (
                    int(tmp_del_info[1].split(',')[0]) +
                    int(tmp_del_info[1].split(',')[1])
                )
        variant['final_metrics']['BRC'] = del_alt_count


    background_read_counts = int(variant['final_metrics']['BRC'])
    variant['final_metrics']['BRR'] = format(
        background_read_counts / int(total_read_count),
        '.5f'
    )
    background_read_ratio = float(variant['final_metrics']['BRR'])
    alt_read_count_ratio = float(variant['final_metrics']['ARR'])/100
    variant['final_metrics']['BRE'] = background_read_ratio / \
                                                 (alt_read_count_ratio + background_read_ratio)
    # scale ratio from [0-1] to [0-100]
    variant['final_metrics']['BRE'] = format(
        float(variant['final_metrics']['BRE']) * 100,
        '.5f'
    )
    variant['final_metrics']['BRR'] = format(
        float(variant['final_metrics']['BRR']) * 100,
        '.5f'
    )


    return(
        variant['final_metrics']['BRC'],
        variant['final_metrics']['BRR'],
        variant['final_metrics']['BRE']
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
    - The input dictionary is expected to have 'final_metrics' containing 'ARR' and 'BRE'.
    - update : 230413 If variant ratio is higher that 30% we consider it automatically clean

    L/P-NO Likely/Probably NOisy
    L/P-CL Likely/Probably CLean
    [ LCL [ PCL [ PNO [          LNO          ]
    |     |     |     |                       |
    0    [6]   [7]   [8]                    100 (BRE)
    """
    alt_read_count_ratio = float(variant['final_metrics']['ARR'])
    if alt_read_count_ratio >= 30:
        variant['final_metrics']['BKG'] = 'PCL'
    else:
        background_read_enrichment = float(variant['final_metrics']['BRE'])

        if background_read_enrichment >= thresholds[8]:
            variant['final_metrics']['BKG'] = 'LNO'
        elif background_read_enrichment >= thresholds[7]:
            variant['final_metrics']['BKG'] = 'PNO'
        elif background_read_enrichment >= thresholds[6]:
            variant['final_metrics']['BKG'] = 'PCL'
        else:
            variant['final_metrics']['BKG'] = 'LCL'

    return variant['final_metrics']['BKG']


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
    #   des = callsets[variant_key]['final_metrics'][descriptor]
    #   callsets[variant_key]['final_metrics'][descriptor] = des
    if len(re.findall('NO', variant['final_metrics']['BKG'])) > 0:
        variant['final_metrics']['FILTER'] = 'FAIL'
    elif variant['final_metrics']['LOW'] == 1:
        variant['final_metrics']['FILTER'] = 'FAIL'
    elif (abs(float(variant['final_metrics']['SBP'])) <= 0.05) and \
        (float(variant['final_metrics']['SBM']) >= sbm):
        variant['final_metrics']['FILTER'] = 'FAIL'
    else:
        variant['final_metrics']['FILTER'] = 'PASS'
    return variant['final_metrics']['FILTER']


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

        chromsome, position, ref, alt = lookup.split(':')

        variant_identifier = f"{ref}:{alt}"

        variant: dict = variants[chromsome][position][variant_identifier]

        if not 'final_metrics' in variant:
            variant['final_metrics'] = {}

        # Save number of Variant callers that found this variant
        variant['final_metrics']['VCN'] = len(variant['VC']['VAF'])

        # keep trace of used VC identifier(s)
        variant['final_metrics']['VCI'] = ','.join(
            sorted(
                variant['VC']['VAF'].keys(),
                key=str.lower
            ))


        # Get mean vaf from Variant callers and convert in %
        vaf_sum = sum(variant['VC']['VAF'].values())
        vcn = float(variant['final_metrics']['VCN'])
        mean_arr = vaf_sum / vcn
        variant['final_metrics']['ARR'] = format(mean_arr,'.5f')

        #Filter on vaf
        variant['final_metrics']['LOW'] = vaf_user_threshold(
            variant,
            thresholds
        )
        variant['final_metrics']['VAR'] = categorize_variant_type(
            variant,
            thresholds
        )


        # Get mean TRC from Variant callers
        variant['final_metrics']['TRC'] = int(round(
            sum(variant['VC']['TRC'].values()) /
            variant['final_metrics']['VCN'],
            5))

        # Get mean ARC from Variant callers
        variant['final_metrics']['ARC'] = int(round(
            (sum(variant['VC']['ARC'].values())/
            variant['final_metrics']['VCN']),
            5))

        # Get mean RRC from Variant callers
        variant['final_metrics']['RRC'] = int(round(
            (sum(variant['VC']['RRC'].values())/
            variant['final_metrics']['VCN']),
            5))


        # Determine GT
        # estimate <GT> from CT data
        variant['final_metrics']['GT'] = estimate_gt(variant)
        # 2. Set WAR if VC do not agree
        variant['final_metrics']['VAR'] = compare_gt(variant)

        # Without pileup, impossible to determine certain values - set -1
        variant['final_metrics']['BRC'] = "-1"
        variant['final_metrics']['BRR'] = "-1"
        variant['final_metrics']['BRE'] = "-1"
        variant['final_metrics']['BKG'] = "-1"


        if 'RRC+' not in variant['VC'] or len(variant['VC']['RRC+']) != len(variant['VC']['RRC']):

            undertermined_ct_values = ["ARC+","ARC-","RRC+","RRC-","SBM","SBP","TRC+","TRC-"]

            for ct_values in undertermined_ct_values:
                variant['final_metrics'][ct_values] = "-1"
                variant['final_metrics']['ARC'], variant['final_metrics']['RRC'] = format_rrc_arc(variant)
        else:
            # Get mean TRC+/-,RRC+/-, ARC+/- from Variant callers
            for strand in["+","-"]:
                for count in ['TRC','RRC','ARC']:
                    if count == 'TRC':
                        sum_strand = sum(variant['VC']['RRC'+strand].values()) + sum(variant['VC']['ARC'+strand].values())
                    else:
                        sum_strand = sum(variant['VC'][count+strand].values())
                        variant['final_metrics'][count+strand] = round(
                            sum_strand/
                            variant['final_metrics']['VCN'],
                            5)

        variant['final_metrics']['ARC'],\
        variant['final_metrics']['RRC'] = format_rrc_arc(variant)

        variant['final_metrics']['SBP'],\
        variant['SBM'] = estimate_sbm(variant,sbm_homozygous)

        variant['final_metrics']['FILTER'] = format_float_descriptors(
            variant,
            sbm
        )
        variant['final_metrics']['PIL'] = "N"
        variant['final_metrics']['RES'] = 'N'

    return variants
