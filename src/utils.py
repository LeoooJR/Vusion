#!/bin/python3

import re

class Cache():

    def __init__(self, func, max_size: int = 1):
        
        self.max_size = max_size

        self.cache: dict = {}

        self.func: function = func

    def add(self, args: list[str], key: object):

        if self.max_size and len(self.cache) >= self.max_size:

            self.cache.clear()

        self.cache[key] = self.func(*args)

    def call(self, args: list[str], key: object):

        if key.__hash__:

            if not key in self.cache:

                self.add(args, key)
            
            return self.cache[key]
        
        else:

            raise TypeError("Key is not hashable.")

# ===========================================================================================
# Basics functions on dictionary
# ===========================================================================================

def merge_collections(collections: list[object]) -> object:
    """
    Merge collections.
    Parameters:
        collections: A list of collections to merge.
    Returns:
        A merged collection of the same type as inputed collections.
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
    ref_and_alt_read_count = variant['sample']['ARC-'] + \
                             variant['sample']['ARC+'] + \
                             variant['sample']['RRC-'] + \
                             variant['sample']['RRC+']
    ins_read_counts = 0

    total_read_count = variant['sample']['TRC']

    if pileup_line_info[14] != 'None':
        # Get number of read with ins format is A:1,0
        tmp_table_count = re.findall(r'\d+', pileup_line_info[14])
        for tmp_count in tmp_table_count:
            ins_read_counts += int(tmp_count)

    if variant['type'] == 'INS':
        ins_read_counts -= (
            variant['sample']['ARC+'] +
            variant['sample']['ARC-']
        )

    # Removing DEL counts if variant is at some position of a DEL.
    # Because we miss valid variants in specific case like that
    # Clintool bug where non-existing deletion is reported and start with A, C, T or G
    if (variant['type'] != 'DEL') and (pileup_line_info[15] != 'None'):

        # Escaping clintool bug where non-existing deletion is reported and start with A, C, T or G
        if pileup_line_info[15][0] == '*':

            del_read_count = int(pileup_line_info[15].strip().split(';')[0].split(':')[1])
            tmp_total_read_count = total_read_count - int(del_read_count)


            variant['sample']['BRC'] = (
                tmp_total_read_count +
                ins_read_counts -
                min([
                    tmp_total_read_count,
                    ref_and_alt_read_count
                    ])
                )

        else:

            variant['sample']['BRC'] = (
                total_read_count +
                ins_read_counts -
                min([
                    total_read_count,
                    ref_and_alt_read_count
                ])
            )

    elif variant['type'] != 'DEL':

        variant['sample']['BRC'] = (
            total_read_count +
            ins_read_counts -
            min([
                total_read_count,
                ref_and_alt_read_count
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


    background_read_counts = variant['sample']['BRC']
    
    variant['sample']['BRR'] = round(
        background_read_counts / total_read_count,
        5
    )
    background_read_ratio = variant['sample']['BRR']

    alt_read_count_ratio = variant['sample']['ARR'] / 100
    variant['sample']['BRE'] = background_read_ratio / \
                                                 (alt_read_count_ratio + background_read_ratio)
    # scale ratio from [0-1] to [0-100]
    variant['sample']['BRE'] = round(
        float(variant['sample']['BRE']) * 100,
        5
    )
    
    variant['sample']['BRR'] = round(
        float(variant['sample']['BRR']) * 100,
        5
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
