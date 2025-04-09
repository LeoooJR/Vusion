#!/bin/python3

import errors
from loguru import logger
import os
import re
import numpy
from scipy.stats import fisher_exact

# ===========================================================================================
# I/O
# ===========================================================================================

def is_empty(path: str) -> bool:
    """ Check if file is empty """
    return os.path.getsize(path) == 0

def is_file(path: str) -> bool:
      """ Check if path is a file """
      return os.path.isfile(path)

def verify_files(file: str, format: str = "vcf"):

    def verify_fasta(file: str):

        if not is_file(file):

            raise errors.FastaIndexError(f"Error: The file {file} does not exist.")

        if is_empty(file):

            raise errors.FastaIndexError(f"Error: The file {file} is empty.")
        
        try:

            with open(file, mode='r') as fasta:

                line = fasta.readline()

                if not line:
                    
                    raise errors.FastaIndexError(f"Error: First line of {file} is empty.")
                
                else:

                    # Check if first line is composed of 5 columns
                    if len(line.split('\t')) != 5:

                        raise errors.FastaIndexError(f"Error: First line inconsistent with Fasta index format")
        
        except FileNotFoundError:

            raise errors.FastaIndexError(f"{file} is not a valid path")

        except IOError:

            raise errors.FastaIndexError(f"An error occurred while reading {file}")

    def verify_vcf(file: str):

        if not is_file(file):

            raise errors.VCFError(f"Error: The file {file} does not exist.")

        if is_empty(file):

            raise errors.VCFError(f"Error: The file {file} is empty.")
        
        try:

            with open(file, mode='r') as vcf:

                line = vcf.readline()

                if not line:

                    raise errors.VCFError(f"Error: First line of {file} is empty.")
                
                else:

                    # Check if first line start with "#"
                    if line[0] != '#':

                        raise errors.VCFError(f"Error: First line inconsistent with VCF header format")

        except FileNotFoundError:

            raise errors.VCFError(f"{file} is not a valid path")
                
        except IOError:

            raise errors.VCFError(f"An error occurred while reading {file}")
                

    def verify_pileup(file: str):

        if not is_file(file):

            raise errors.PileupError(f"Error: The file {file} does not exist.")

        if is_empty(file):

            raise errors.PileupError(f"Error: The file {file} is empty.")

    if format == "vcf":

        verify_vcf(file)
    
    elif format == "fai":

        verify_fasta(file)

    elif format == "pileup":

        verify_pileup(file)

    else:
        logger.error(f"File format {format} cannot be verified.")
        raise ValueError("File format cannot be verified.")

# ===========================================================================================
# Read VCF
# ===========================================================================================
    
def samtools_vaf(line_dict):
    """
    Functions to get VAF in different file format
    SAMTOOLS
    Get depth value and variant depth on reverse and forward,
    then divide the total of variant per total depth
    """
    depth = float(line_dict[7].split('DP=')[1].split(';')[0])
    depth_variant_info = line_dict[7].split('DP4=')[1].split(';')[0]
    variant_depth = float(depth_variant_info.split(',')[2]) + \
                float(depth_variant_info.split(',')[3])
    vaf = variant_depth/depth
    return vaf


def varscan_vaf(line_dict):
    """
    VARSCAN
    Get VAF in percent in sample info and return it divided by 100
    """
    vaf = float(line_dict[9].split('%')[0].split(':')[-1])/100
    return vaf


def vardict_vaf(line_dict):
    """
    VARDICT
    Get allele frequency in info
    """
    vaf = float(line_dict[7].split('AF=')[1].split(';')[0])
    return vaf


def pindel_vaf(line_dict):
    """
    PINDEL
    Get total depth and variant depth return vaf
    """
    depths = line_dict[9].split(':')[1].split(',')
    vaf = float(depths[1])/(float(depths[0]) + float(depths[1]))
    return vaf


def haplotypecaller_vaf(line_dict):
    """
    HAPLOTYPECALLER
    Get vaf from info
    """
    total_depth = float(line_dict[9].split(':')[2])
    variant_depth = float(line_dict[9].split(':')[1].split(',')[1])
    vaf = variant_depth / total_depth
    # Weirdly Allele Frequency of HC is always 1 or 0.5,
    #so we compute it instead of reporting the given AF
    #vaf = float(line_dict[7].split('AF=')[1].split(';')[0])
    return vaf


def haplotypecaller_coverage(line_dict):
    """
    HAPLOTYPECALLER
    Get TRC from info
    """
    return float(line_dict[9].split(':')[2])


def varscan_coverage(line_dict):
    """
    varscan
    Get TRC from info
    """
    return float(line_dict[9].split(':')[2])


def vardict_coverage(line_dict):
    """
    vardict
    Get TRC from info
    """
    return float(line_dict[9].split(":")[1])


def st_coverage(line_dict):
    """
    samtools
    Get TRC from info
    """
    return float(line_dict[7].split('DP=')[1].split(';')[0])


def pindel_coverage(line_dict):
    """
    pindel
    Get TRC from info
    """
    return (
        float(line_dict[9].split(':')[1].split(',')[0])+
        float(line_dict[9].split(':')[1].split(',')[1])
    )


def varscan_rrc(line_dict):
    """
    varscan
    Get RRC from info
    """
    rrc = float(line_dict[9].split(':')[4])
    rrc_plus = float(line_dict[9].split(':')[10])
    rrc_minus = float(line_dict[9].split(':')[11])
    return(rrc,rrc_plus,rrc_minus)


def vardict_rrc(line_dict):
    """
    vardict
    Get RRC from info
    """
    rrc = float(line_dict[9].split(":")[3].split(',')[0])
    rrc_plus = float(line_dict[9].split(":")[5].split(',')[0])
    rrc_minus = float(line_dict[9].split(":")[5].split(',')[1])
    return(rrc,rrc_plus,rrc_minus)


def st_rrc(line_dict):
    """
    samtools
    Get RRC from info
    """
    depth_variant_info = line_dict[7].split('DP4=')[1].split(';')[0]
    rrc = float(depth_variant_info.split(',')[0]) + float(depth_variant_info.split(',')[1])
    rrc_plus = float(depth_variant_info.split(',')[0])
    rrc_minus = float(depth_variant_info.split(',')[1])
    return(rrc,rrc_plus,rrc_minus)


def varscan_arc(line_dict):
    """
    varscan
    Get ARC from info
    """
    arc = float(line_dict[9].split(':')[5])
    arc_plus = float(line_dict[9].split(':')[12])
    arc_minus = float(line_dict[9].split(':')[13])
    return(arc,arc_plus,arc_minus)


def vardict_arc(line_dict):
    """
    vardict
    Get ARC from info
    """
    arc = float(line_dict[9].split(":")[3].split(',')[1])
    arc_plus = float(line_dict[9].split(":")[6].split(',')[0])
    arc_minus = float(line_dict[9].split(":")[6].split(',')[1])
    return(arc,arc_plus,arc_minus)


def haplotypecaller_arc(line_dict):
    """
    HC
    Get ARC from info
    """
    return float(line_dict[9].split(':')[1].split(',')[1])


def st_arc(line_dict):
    """
    samtools
    Get ARC from info
    """
    depth_variant_info = line_dict[7].split('DP4=')[1].split(';')[0]
    arc = float(depth_variant_info.split(',')[2]) + float(depth_variant_info.split(',')[3])
    arc_plus = float(depth_variant_info.split(',')[2])
    arc_minus = float(depth_variant_info.split(',')[3])
    return(arc,arc_plus,arc_minus)


def pindel_arc(line_dict):
    """
    pindel
    Get ARC from info
    """
    return float(line_dict[9].split(':')[1].split(',')[1])

"""
FILT3R
"""
def FLT3R_vaf(line_dict):
    return float(line_dict[6].split(';')[4].split('=')[1])

def FLT3_arc(line_dict):
    return float(line_dict[6].split(';')[1].split('=')[1])

def FLT3_coverage(line_dict):
    return (
        float(line_dict[6].split(';')[1].split('=')[1])+
        float(line_dict[6].split(';')[2].split('=')[1])
    )
# ===========================================================================================
# Output VCF
# ===========================================================================================
def print_var(variant_id, variant_dict, variant_type):
    """
    Function to print variant in vcf format

    Parameters:
    - variant_id (str): The unique identifier for the variant
    - variant_dict (dict): A dictionary containing variant information
    - type (str): The type of variant information to print (e.g., 'final_metrics')

    Returns:
    str: Returns a formatted string representing the variant in VCF format.

    Note:
    - The function extracts relevant information for CHROM, POS, ID, REF, ALT,
    QUAL, FILTER, VAR, and GT
    - It set ID and QUAL at '.'as placeholders.
    """

    line = '\t'.join([
        variant_dict[variant_id]["VC"]["CHROM"],
        variant_dict[variant_id]["VC"]["POS"],
        '.',
        variant_dict[variant_id]["VC"]["REF"],
        variant_dict[variant_id]["VC"]["ALT"],
        '.'
    ])
    line += '\t' + variant_dict[variant_id][variant_type]['FILTER']
    line += '\tVAR=' + variant_dict[variant_id]['VT']
    line += '\tGT:' + ':'.join(variant_dict[variant_id]['vcf_fields'])
    line += '\t' + variant_dict[variant_id][variant_type]['GT'] + ':'

    fields_info = []
    for field in variant_dict[variant_id]['vcf_fields']:
        fields_info.append(str(variant_dict[variant_id][variant_type][field]))

    line += ':'.join(fields_info)
    line += '\n'
    return line


def print_var4cartagenia(variant_id, variant_dict, variant_type):
    """
    Function to print variant in cartagenia format

    Parameters:
    - variant_id (str): The unique identifier for the variant.
    - variant_dict (dict): A dictionary containing variant information.
    - type (str): The type of variant information to print (e.g., 'final_metrics').

    Returns:
    str: Returns a formatted string representing the variant in Cartagenia format
    """
    fields_to_remove = ['VCI', 'VCN', 'PIL','RES']
    for field in fields_to_remove:
        del variant_dict[variant_id]['final_metrics'][field]
        variant_dict[variant_id]['vcf_fields'].remove(field)


    variant_id_info = variant_id.split(':')
    # Format genotype
    if variant_dict[variant_id][variant_type]['GT'] == '0/0':
        variant_dict[variant_id][variant_type]['GT'] = '0/1'
    # Change ARR, BRR and BRE format
    for variable in ['ARR', 'BRR', 'BRE']:
        if variable in variant_dict[variant_id][variant_type]:
            variant_dict[variant_id][variant_type][variable] = format(
                float(variant_dict[variant_id][variant_type][variable])/100,
                '.5f'
            )

    # Get only total read count for ARC and RRC
    # instead of minus_strand_count,plus_strand_count,total_count
    for variable in ['RRC', 'ARC']:
        if variable in variant_dict[variant_id][variant_type]:
            variant_dict[variant_id][variant_type][variable] = (
                variant_dict[variant_id][variant_type][variable].split(',')[2]
            )
    # Set CHROM POS ID REF ALT QUAL
    # Get CHROM, POS, REF and ALT
    # And set ID and QUAL at '.'
    line = '\t'.join([
        variant_dict[variant_id]["VC"]["CHROM"],
        variant_dict[variant_id]["VC"]["POS"],
        '.',
        variant_dict[variant_id]["VC"]["REF"],
        variant_dict[variant_id]["VC"]["ALT"],
        '.'
    ])
    line += '\t' + variant_dict[variant_id][variant_type]['FILTER']
    line += '\tVAR=' + variant_dict[variant_id]['VT']
    line += '\tGT:' + ':'.join(variant_dict[variant_id]['vcf_fields'])

    line += '\t' + variant_dict[variant_id][variant_type]['GT'] + ':'
    fields_info = []

    for field in variant_dict[variant_id]['vcf_fields']:
        fields_info.append(str(variant_dict[variant_id][variant_type][field]))

    line += ':'.join(fields_info)
    line += '\n'
    return line


def order_var(variants: dict, contigs):
    """
    Order variant for VCF easy reading

    Parameters:
    - variants (dict): variants
    - contigs: contig IDs

    Returns:
    list: Returns a list of ordered variant identifiers based on contig information
    """

    # variant_info_ordered_list = []

    # for contig in sorted(contigs):

    #     tmp_variant_tuple = []

    #     for key in variants:

    #         variant_info = [variants[key]["VC"]["CHROM"], variants[key]["VC"]["POS"], variants[key]["VC"]["REF"], variants[key]["VC"]["ALT"]]
            
    #         if variants["VC"]["CHROM"] == contig:
    #             tmp_variant_tuple.append(variant_info)

    #     if len(tmp_variant_tuple) > 0:
    #         variant_info_ordered_list.append(
    #             sorted(tmp_variant_tuple, key=lambda variant: [int(variant[1]), variant[2]])
    #         )

    # variant_id_ordered_list = []
    # for chr_list_of_variant in variant_info_ordered_list:
    #     for variant in chr_list_of_variant:
    #         variant_id_ordered_list.append(':'.join(variant))
    # return variant_id_ordered_list

    return [k for k, v in sorted(variants.items(), key=lambda item: [item[1]["VC"]["CHROM"], int(item[1]["VC"]["POS"])] )]


# ===========================================================================================
# Basics functions on dictionary
# ===========================================================================================
def merge_variant_dict(dicts):
    """
    Merges specified variant dictionaries.
    Parameters:
        dicts: A list of variant dictionaries
    Returns:
        A merged variant dictionary
    """
    merged_dict = {}
    for d in dicts:
        for key, value in d.items():
            if key in merged_dict:
                print(f"Warning: Key '{key}' already exists. Value '{value}' will be overwritten.")
            merged_dict[key] = value
    return merged_dict


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
        variant_type = 'SNV'
    elif len(ref) == 1 and len(alt) > 1 :
        variant_type = 'INS'
    elif len(ref) > 1 and len(alt) == 1:
        variant_type = 'DEL'
    elif len(ref) == len(alt):
        if ref == rev :
            variant_type = 'INV'
        else:
            # mostly coming from <VD> and both annotated as <Complex>
            variant_type = 'MNV'
    else:
        # mostly coming from <PL> and annotated as either <INV> or <RPL>
        variant_type = 'CSV'
    return variant_type


def estimate_sbm(dic,variant_key,sbm_homozygous):
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
    alt_plus_x_ref_minus = int(dic[variant_key]['final_metrics']['ARC+']) * \
                           int(dic[variant_key]['final_metrics']['RRC-'])
    alt_minus_x_ref_plus = int(dic[variant_key]['final_metrics']['ARC-']) * \
                           int(dic[variant_key]['final_metrics']['RRC+'])

    # Check if there is a zero count or if variant is homozygous
    if (
        alt_plus_x_ref_minus==0 or
        alt_minus_x_ref_plus==0 or
        dic[variant_key]['final_metrics']['GT'] == "1/1"
    ):
        strand_max = max(
            int(dic[variant_key]['final_metrics']['ARC+']),
            int(dic[variant_key]['final_metrics']['ARC-'])
        )
        if strand_max !=0:
            strand_ratio = strand_max/(
                int(dic[variant_key]['final_metrics']['ARC+'])+
                int(dic[variant_key]['final_metrics']['ARC-'])
            ) #ratio of reads on one strand to the other
        else:
            strand_ratio = 0 # no ALT allele = no stand bias

        dic[variant_key]['final_metrics']['SBP'] = 0.05
        if strand_ratio >= sbm_homozygous:
            dic[variant_key]['final_metrics']['SBM'] = '1.50000'
        else:
            dic[variant_key]['final_metrics']['SBM'] = '0.50000'

    # Otherwise, perform standard Fisher Test
    else:
        m = max(alt_plus_x_ref_minus,alt_minus_x_ref_plus)
        # no <REF> allele = no strand-bias (case of illegal division by zero)
        # if m == 0:
        #     dic[variant_key]['final_metrics']['SBM'] = '0.50000'
        # else:
            # dic[variant_key]['final_metrics']['SBM'] = format(
            #     (m / (alt_plus_x_ref_minus + alt_minus_x_ref_plus)),
            #     '.5f'
            # )

        # [n11,n12] n1p
        # [n21,n22] n2p
        # np1 np2  npp
        n11 = int(dic[variant_key]['final_metrics']['ARC+'])
        n12 = int(dic[variant_key]['final_metrics']['ARC-'])

        if (dic[variant_key]['VT'] == 'INS') or (dic[variant_key]['VT'] == 'INV'):
            n21 = (
                int(dic[variant_key]['final_metrics']['TRC+']) -
                int(dic[variant_key]['final_metrics']['ARC+'])
            )
            n22 = (
                int(dic[variant_key]['final_metrics']['TRC-']) -
                int(dic[variant_key]['final_metrics']['ARC-'])
            )
            # 131223 : cas rare avec bug de CT ou TRC+<ARC+...
            if n21<0 :
                n21=0
            if n22<0:
                n22=0
        else:
            n21 = int(dic[variant_key]['final_metrics']['RRC+'])
            n22 = int(dic[variant_key]['final_metrics']['RRC-'])


        dic[variant_key]['final_metrics']['SBM'] = format(
            (m / (alt_plus_x_ref_minus + alt_minus_x_ref_plus)),
            '.5f'
        )
        table_for_fisher_test = numpy.array([[n11, n12], [n21, n22]])
        oddsr, p = fisher_exact(table_for_fisher_test, alternative='two-sided')
        #dic[variant_key]['final_metrics']['SBP'] = format(p, '.5f')
        dic[variant_key]['final_metrics']['SBP'] = round(p,5)

    return(dic[variant_key]['final_metrics']['SBP'],dic[variant_key]['final_metrics']['SBM'])


def categorize_variant_type(dic,variant_key,thresholds):
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
    alt_read_count_ratio = float(dic[variant_key]['final_metrics']['ARR'])

    if alt_read_count_ratio > thresholds[5]:
        dic[variant_key]['final_metrics']['VAR'] = 'LHO'
    elif alt_read_count_ratio > thresholds[4]:
        dic[variant_key]['final_metrics']['VAR'] = 'PHO'
    elif alt_read_count_ratio < thresholds[0]:
        dic[variant_key]['final_metrics']['VAR'] = 'LSC'
    elif alt_read_count_ratio < thresholds[1]:
        dic[variant_key]['final_metrics']['VAR'] = 'PSC'
    elif (alt_read_count_ratio < thresholds[2]) or (alt_read_count_ratio > thresholds[3]):
        dic[variant_key]['final_metrics']['VAR'] = 'PHE'
    else:
        dic[variant_key]['final_metrics']['VAR'] = 'LHE'
    return dic[variant_key]['final_metrics']['VAR']


def estimate_gt(dic,variant_key):
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
    if len(re.findall('SC',dic[variant_key]['final_metrics']['VAR'])) > 0:
        dic[variant_key]['final_metrics']['GT'] = '0/0'
    elif len(re.findall('HE',dic[variant_key]['final_metrics']['VAR'])) > 0:
        dic[variant_key]['final_metrics']['GT'] = '0/1'
    else:
        dic[variant_key]['final_metrics']['GT'] = '1/1'
    return dic[variant_key]['final_metrics']['GT']


def compare_gt(dic,variant_key):
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
    callers = list(dic[variant_key]['VC']['VAF'].keys())
    pileup_gt = dic[variant_key]['final_metrics']['GT']

    if 'FL' in callers:
        # If FL in callers, callers cant be agreed
        # Because FL do not compute GT
        dic[variant_key]['final_metrics']['VAR'] = 'WAR'
    else:

        # gt_list = dic[variant_key]['VC']['GT'].copy()
        # gt_list[:] = [x if x != './1' else '0/1' for x in gt_list]
        # => TypeError: unhashable type: 'slice' => normal, jai change la structure du dic pour les GT
        gt_list = [] #List of all GT identified for a variant
        for variant_caller in callers:

            gt_list.append(dic[variant_key]['VC']['GT'][variant_caller])

        # Replace "./1" by "0/1", "1/." by "0/1", and "0/0" by "0/1".
        gt_list[:] = [x if x != './1' else '0/1' for x in gt_list]
        gt_list[:] = [x if x != '1/.' else '0/1' for x in gt_list]
        gt_list[:] = [x if x != '0/0' else '0/1' for x in gt_list]

        for genotype in gt_list:
            ok_gt = ['1/1', '0/1']
            if genotype not in ok_gt:
                print("!! WARNING !!" + variant_key + " has a GT looking like that " + genotype + "\n")

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
            dic[variant_key]['final_metrics']['VAR'] = 'WAR'


    return dic[variant_key]['final_metrics']['VAR']


def vaf_user_threshold(dic,variant_key,thresholds):
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
    if float(dic[variant_key]['final_metrics']['ARR']) >= thresholds[-1]:
        dic[variant_key]['final_metrics']['LOW'] = 0
    else:
        dic[variant_key]['final_metrics']['LOW'] = 1
    return dic[variant_key]['final_metrics']['LOW']


def format_rrc_arc(dic,variant_key):
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
        read_count_minus = dic[variant_key]['final_metrics'][read_count_type + '-']
        read_count_plus = dic[variant_key]['final_metrics'][read_count_type + '+']

        if read_count_minus!="-1" and read_count_plus != "-1":
            tmp_read_count = int(read_count_minus) + int(read_count_plus)
            dic[variant_key]['final_metrics'][read_count_type] = ','.join(
            [
                str(int(read_count_plus)),
                str(int(read_count_minus)),
                str(tmp_read_count)
            ])
        else:
            dic[variant_key]['final_metrics'][read_count_type] = ','.join(
            [
                str(read_count_plus),
                str(read_count_minus),
                str(dic[variant_key]['final_metrics'][read_count_type])
            ])

    return(dic[variant_key]['final_metrics']['ARC'],dic[variant_key]['final_metrics']['RRC'])


def estimate_brc_r_e(dic,variant_key,pileup_line_info):
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
    ref_and_alt_read_count = int(dic[variant_key]['final_metrics']['ARC-']) + \
                             int(dic[variant_key]['final_metrics']['ARC+']) + \
                             int(dic[variant_key]['final_metrics']['RRC-']) + \
                             int(dic[variant_key]['final_metrics']['RRC+'])
    ins_read_counts = 0
    total_read_count = dic[variant_key]['final_metrics']['TRC']

    if pileup_line_info[14] != 'None':
        # Get number of read with ins format is A:1,0
        tmp_table_count = re.findall(r'\d+', pileup_line_info[14])
        for tmp_count in tmp_table_count:
            ins_read_counts += int(tmp_count)
    if dic[variant_key]['VT'] == 'INS':
        ins_read_counts -= (
            int(dic[variant_key]['final_metrics']['ARC+']) +
            int(dic[variant_key]['final_metrics']['ARC-'])
        )

    # Removing DEL counts if variant is at some position of a DEL.
    # Because we miss valid variants in specific case like that
    # Clintool bug where non-existing deletion is reported and start with A , C, T or G
    if (dic[variant_key]['VT'] != 'DEL') and (pileup_line_info[15] != 'None'):
        # Escaping clintool bug where non-existing deletion is reported and start with A , C, T or G
        if pileup_line_info[15][0] == '*':
            del_read_count = int(pileup_line_info[15].strip().split(';')[0].split(':')[1])
            tmp_total_read_count = int(total_read_count) - int(del_read_count)


            dic[variant_key]['final_metrics']['BRC'] = (
                tmp_total_read_count +
                int(ins_read_counts) -
                min([
                    tmp_total_read_count,
                    int(ref_and_alt_read_count)
                    ])
                )

        else:
            dic[variant_key]['final_metrics']['BRC'] = (
                int(total_read_count) +
                int(ins_read_counts) -
                min([
                    int(total_read_count),
                    int(ref_and_alt_read_count)
                ])
            )
    elif dic[variant_key]['VT'] != 'DEL':
        dic[variant_key]['final_metrics']['BRC'] = (
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
            if (tmp_del_info[0] != '*') and (tmp_del_info[0] != (dic[variant_key]["VC"]["REF"])[1:]):
                del_alt_count += (
                    int(tmp_del_info[1].split(',')[0]) +
                    int(tmp_del_info[1].split(',')[1])
                )
        dic[variant_key]['final_metrics']['BRC'] = del_alt_count


    background_read_counts = int(dic[variant_key]['final_metrics']['BRC'])
    dic[variant_key]['final_metrics']['BRR'] = format(
        background_read_counts / int(total_read_count),
        '.5f'
    )
    background_read_ratio = float(dic[variant_key]['final_metrics']['BRR'])
    alt_read_count_ratio = float(dic[variant_key]['final_metrics']['ARR'])/100
    dic[variant_key]['final_metrics']['BRE'] = background_read_ratio / \
                                                 (alt_read_count_ratio + background_read_ratio)
    # scale ratio from [0-1] to [0-100]
    dic[variant_key]['final_metrics']['BRE'] = format(
        float(dic[variant_key]['final_metrics']['BRE']) * 100,
        '.5f'
    )
    dic[variant_key]['final_metrics']['BRR'] = format(
        float(dic[variant_key]['final_metrics']['BRR']) * 100,
        '.5f'
    )


    return(
        dic[variant_key]['final_metrics']['BRC'],
        dic[variant_key]['final_metrics']['BRR'],
        dic[variant_key]['final_metrics']['BRE']
    )


def categorize_background_signal(dic,variant_key,thresholds):
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
    alt_read_count_ratio = float(dic[variant_key]['final_metrics']['ARR'])
    if alt_read_count_ratio >= 30:
        dic[variant_key]['final_metrics']['BKG'] = 'PCL'
    else:
        background_read_enrichment = float(dic[variant_key]['final_metrics']['BRE'])

        if background_read_enrichment >= thresholds[8]:
            dic[variant_key]['final_metrics']['BKG'] = 'LNO'
        elif background_read_enrichment >= thresholds[7]:
            dic[variant_key]['final_metrics']['BKG'] = 'PNO'
        elif background_read_enrichment >= thresholds[6]:
            dic[variant_key]['final_metrics']['BKG'] = 'PCL'
        else:
            dic[variant_key]['final_metrics']['BKG'] = 'LCL'

    return dic[variant_key]['final_metrics']['BKG']


def format_float_descriptors(dic,variant_key,sbm):
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
    if len(re.findall('NO', dic[variant_key]['final_metrics']['BKG'])) > 0:
        dic[variant_key]['final_metrics']['FILTER'] = 'FAIL'
    elif dic[variant_key]['final_metrics']['LOW'] == 1:
        dic[variant_key]['final_metrics']['FILTER'] = 'FAIL'
    elif (abs(float(dic[variant_key]['final_metrics']['SBP'])) <= 0.05) and \
        (float(dic[variant_key]['final_metrics']['SBM']) >= sbm):
        dic[variant_key]['final_metrics']['FILTER'] = 'FAIL'
    else:
        dic[variant_key]['final_metrics']['FILTER'] = 'PASS'
    return dic[variant_key]['final_metrics']['FILTER']


def process_without_pileup(dico,thresholds,sbm,sbm_homozygous):
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
    for variant_key in list(dico.keys()):
        if 'final_metrics' not in dico[variant_key]:
            dico[variant_key]['final_metrics'] = {}

        # Save number of Variant callers that found this variant
        dico[variant_key]['final_metrics']['VCN'] = len(dico[variant_key]['VC']['VAF'])

        # keep trace of used VC identifier(s)
        dico[variant_key]['final_metrics']['VCI'] = ','.join(
            sorted(
                dico[variant_key]['VC']['VAF'].keys(),
                key=str.lower
            ))


        # Get mean vaf from Variant callers and convert in %
        vaf_sum = sum(dico[variant_key]['VC']['VAF'].values())
        vcn = float(dico[variant_key]['final_metrics']['VCN'])
        mean_arr = vaf_sum / vcn
        dico[variant_key]['final_metrics']['ARR'] = format(mean_arr,'.5f')

        #Filter on vaf
        dico[variant_key]['final_metrics']['LOW'] = vaf_user_threshold(
            dico,
            variant_key,
            thresholds
        )
        dico[variant_key]['final_metrics']['VAR'] = categorize_variant_type(
            dico,
            variant_key,
            thresholds
        )


        # Get mean TRC from Variant callers
        dico[variant_key]['final_metrics']['TRC'] = int(round(
            sum(dico[variant_key]['VC']['TRC'].values()) /
            dico[variant_key]['final_metrics']['VCN'],
            5))

        # Get mean ARC from Variant callers
        dico[variant_key]['final_metrics']['ARC'] = int(round(
            (sum(dico[variant_key]['VC']['ARC'].values())/
            dico[variant_key]['final_metrics']['VCN']),
            5))

        # Get mean RRC from Variant callers
        dico[variant_key]['final_metrics']['RRC'] = int(round(
            (sum(dico[variant_key]['VC']['RRC'].values())/
            dico[variant_key]['final_metrics']['VCN']),
            5))


        # Determine GT
        # 1. estimate <GT> from mean VAF
        # Ancienne version utilisee avant pour le rescue -
        # reviens au meme que le resultat de categorize_variant_type
        #- ne comprends pas pk 2 methodes diff
        # ---
        # if float(dico[variant_key]['VC']['VCR']) > thresholds[4] :
        #     dico[variant_key]['final_metrics']['GT'] = '1/1'
        # elif float(dico[variant_key]['VC']['VCR']) >= thresholds[1] :
        #     dico[variant_key]['final_metrics']['GT'] = '0/1'
        # else:
        #     dico[variant_key]['final_metrics']['GT'] = '0/0'
        # estimate <GT> from CT data
        dico[variant_key]['final_metrics']['GT'] = estimate_gt(dico,variant_key)
        # 2. Set WAR if VC do not agree
        dico[variant_key]['final_metrics']['VAR'] = compare_gt(dico,variant_key)

        # Without pileup, impossible to determine certain values - set -1
        dico[variant_key]['final_metrics']['BRC'] = "-1"
        dico[variant_key]['final_metrics']['BRR'] = "-1"
        dico[variant_key]['final_metrics']['BRE'] = "-1"
        dico[variant_key]['final_metrics']['BKG'] = "-1"


        if 'RRC+' not in dico[variant_key]['VC'] or \
        len(dico[variant_key]['VC']['RRC+']) != len(dico[variant_key]['VC']['RRC']):

            undertermined_ct_values = ["ARC+","ARC-","RRC+","RRC-","SBM","SBP","TRC+","TRC-"]
            for ct_values in undertermined_ct_values:
                dico[variant_key]['final_metrics'][ct_values] = "-1"
            dico[variant_key]['final_metrics']['ARC'],\
            dico[variant_key]['final_metrics']['RRC'] = format_rrc_arc(dico,variant_key)
        else:
            # Get mean TRC+/-,RRC+/-, ARC+/- from Variant callers
            for strand in["+","-"]:
                for count in ['TRC','RRC','ARC']:
                    if count == 'TRC':
                        sum_strand = sum(dico[variant_key]['VC']['RRC'+strand].values()) \
                        + sum(dico[variant_key]['VC']['ARC'+strand].values())
                    else:
                        sum_strand = sum(dico[variant_key]['VC'][count+strand].values())
                    dico[variant_key]['final_metrics'][count+strand] = round(
                        sum_strand/
                        dico[variant_key]['final_metrics']['VCN'],
                        5)

            dico[variant_key]['final_metrics']['ARC'],\
            dico[variant_key]['final_metrics']['RRC'] = format_rrc_arc(dico,variant_key)

            dico[variant_key]['final_metrics']['SBP'],\
            dico[variant_key]['SBM'] = estimate_sbm(dico,variant_key,sbm_homozygous)

        dico[variant_key]['final_metrics']['FILTER'] = format_float_descriptors(
            dico,
            variant_key,
            sbm
        )
        dico[variant_key]['final_metrics']['PIL'] = "N"
        dico[variant_key]['final_metrics']['RES'] = 'N'


    return dico
