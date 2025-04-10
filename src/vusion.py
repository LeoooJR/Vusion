#!/usr/bin/python3

"""
HOW TO RUN :

python3 CombineVCF2Leaves-4.py \
-r ~/Bureau/8072/Data/REF/human_g1k_v37_decoy.fasta.fai \
-V ST,45-CEREBMD-12_S12_ST_final.vcf \
-V VD,45-CEREBMD-12_S12_VD_final.vcf \
-V PL,45-CEREBMD-12_S12_PL_final.vcf \
-V HC,45-CEREBMD-12_S12_HC_final.vcf \
-V VS,45-CEREBMD-12_S12_VS_final.vcf \
-V FL,45-CEREBMD-12_S12_VS_FILT3R_filtered.vcf \
--pileup 45-CEREBMD-12_S12_CT.tsv \
-o test.vcf \
--sample_id test \
-t 10,30,40,60,70,80,20,30,50,1 \


DICTIONNARY :
'chr:pos:REF:ALT':{'VC':{
                         'VAF': {'PL': 8.968609865470851},
                         'GT': {'PL': '0/0'},
                         'FILTER': {'PL': 'PASS'},
                         'INFO': {'PL': ...;SVLEN=34;SVTYPE=INS'},
                         'FORMAT': {'PL': 'GT:AD'},
                         'SAMPLE': {'PL': '0/0:203,20'},
                         'TRC': {'PL': 223.0},
                         'ARC': {'PL': 20.0},
                         'RRC': {'PL': 203.0}},
                    'VT': 'INS',
                    'final_metrics': {'TRC': 223, 'TRC-': '-1', 'TRC+': '-1', \
                    'RRC-': '-1', 'RRC+': '-1', 'VCN': 1, 'VCI': 'PL', \
                    'ARR': '8.96861', 'LOW': 0, 'VAR': 'LSC', 'ARC': '-1,-1,20.0',\
                    'RRC': '-1,-1,203.0', 'GT': '0/0', 'BRC': '-1', 'BRR': '-1', \
                    'BRE': '-1', 'BKG': '-1', 'ARC+': '-1', 'ARC-': '-1', 'SBM': '-1',\
                    'SBP': '-1', 'FILTER': 'PASS', 'PIL': 'N'},
                    'vcf_fields': ['VAR', 'BKG', 'TRC', 'RRC', 'ARC', 'BRC', 'ARR', \
                    'BRR', 'BRE', 'SBP', 'SBM', 'LOW', 'VCI', 'VCN','PIL']}}
"""

import callers
import copy
import errors
import files
from hashlib import sha256
from loguru import logger
import os
import pandas as pd
import re
import sys

sys.path.insert(1, os.path.dirname(os.path.abspath(__file__)))
import utils as functions

def combine(params):

    # ===========================================================================================
    # Initiate variables
    # ===========================================================================================
    VCF_HEADER = ''
    SBM = 0.95
    MAX_THRESHOLD = 100.0
    MIN_THRESHOLD = 0.0
    HEADER = {
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

    # rejected variants (if any; will be written to $opt{trash_file})
    rejected = {}
    rev = {}
    INV_MNV_CSV = {}
    FLiT3r = {}
    readpile = {}
    ITD = {}
    #  thresholds = [10, 30, 40, 60, 70, 80, 20, 30, 50, 1]

    repository = callers.VariantCallerRepository()

    # Check reference genome index
    try:
        fasta_index = files.FastaIndex(path=params.reference, lazy=False)
    except errors.FastaIndexError:
        logger.error(f"{params.reference} is not a valid FASTA index.")
        raise SystemExit(f"{params.reference} is not a valid FASTA index.")
    
    logger.success(f"Fasta index {params.reference} has been successfully checked.")
    
    # Check pileup
    try:
        pileup = files.Pileup(path=params.pileup, lazy=True)
    except errors.PileupError:
        logger.error(f"{params.pileup} is not a valid PILEUP.")
        raise SystemExit(f"{params.pileup} is not a valid PILEUP.")
    
    logger.success(f"Pileup {params.pileup} has been successfully checked.")

    # Check VCFs
    vcfs: dict = {}

    for vcf in params.vcfs:

        input = vcf.split(',')

        if len(input) != 2:
            logger.error('Wrong number of argument in --vcf option.')
            logger.error(f'Error was raised by: {input[1]}.')
            raise ValueError("Wrong number of argument in --vcf option.")
        
        try:
            int(input[0])
            logger.error("Wrong type of argument in --vcf option.")
            logger.error(f"Error was raised by: {input[1]}.")
            raise SystemExit("Wrong type of argument in --vcf option.")
        except ValueError:
            if not repository.is_supported(input[0]):
                logger.error(f"{input[0]} caller is not supported in --vcf option.")
                raise SystemExit("Caller not supported in --vcf options.")
            
        try:
            vcfs[input[0]] = {"vcf": files.VCF(path=input[1], caller=repository.get_VC(input[0]), lazy=True), "index": None}
        except (errors.VCFError, errors.VariantCallerError) as e:
            if isinstance(e,errors.VCFError):
                logger.error(f"{vcf} is not a valid VCF.")
                logger.error(f"Error: {e}")
                raise SystemExit(f"{vcf} is not a valid VCF.")
            else:
                logger.error(f"{input[0]} is not a supported variant caller.")
                logger.error(f"Error: {e}")
                raise SystemExit(f"{input[0]} is not a supported variant caller.")

        logger.debug(f"Variant Callers inputed: {input[0]}")

    # Load reference genome dict ie. list of ordered contig names / length
    contigs: list[pd.Series] = []

    with open(fasta_index.get_path(), mode='r') as ref:
        logger.debug(f"Parsing {fasta_index.get_path()}.")
        for line in ref:
            if line:
                contigs.append(pd.Series(data=line.strip().split('\t')))
                
    contigs: pd.DataFrame = pd.DataFrame(data=contigs)
    contigs.columns=["contig", "length", "index", "pbline", "byteline"]
    contigs = contigs.astype({"contig": "string", "pbline": 'uint8', "byteline": "uint8"})

    # Check if all mandatory option are given and modify variable depending of given options

    if params.disable_strand_bias:
        SBM = 2

    if params.thresholds:

        thresholds = params.thresholds.split(',')

        # Check that we have 10 values in thresholds
        if len(thresholds) != 10:
            logger.error(f"Invalid number of values in --threshold option.")
            raise SystemExit(f"Invalid number of values in --threshold option.")
        
        #Check that all values can be converted to floats
        try:
            thresholds = list(map(float, thresholds)) 
        except ValueError:
            logger.error("Invalid values in thresholds option.")
            raise SystemExit("Invalid values in thresholds option.")
        
        if not any(list(map(lambda value: value >= MIN_THRESHOLD and value <= MAX_THRESHOLD, thresholds))):
            logger.error("Option --thresholds values cannot be equal or higher than 100, or lower than 0.")
            raise SystemExit("Option --thresholds values cannot be equal or higher than 100, or lower than 0.")

        # Check that 6 first given threshold are unique
        if len(thresholds[0:6]) != len(set(thresholds[0:6])):
            logger.error("Option --thresholds six first values must be unique.")
            raise SystemExit("Option --thresholds six first values must be unique.")

        # Check that second group of threshold values are unique
        if len(thresholds[6:9]) != len(set(thresholds[6:9])):
            logger.error("Option --thresholds values 7, 8 and 9 must be unique.")
            raise SystemExit("Option --thresholds values 7, 8 and 9 must be unique.")
        
        thresholds[0:6] = sorted(thresholds[0:6])
        thresholds[6:9] = sorted(thresholds[6:9])

        logger.debug(f"Thresholds: {thresholds}")

    for _, caller in enumerate(vcfs):

        if not _:

            calls = {}  # dictionary of unique variant calls

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

                # Create a variant unique id from CHROM, POSITION, REF and ALT
                hash: str = sha256(
                            string=f"{(datas[0]).removeprefix('chr')}:{datas[1]}:{datas[3]}:{'|'.join(datas[4])}".encode()
                        ).hexdigest()
                
                if not hash in calls:

                    calls[hash] = {"VC": {'CHROM': datas[0],
                                        'POS': datas[1],
                                        'REF': datas[3],
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
                                        'RRC+': {},}}
                
                    ## 1/ define <VT>
                    calls[hash]['VT'] = functions.define_variant_type(ref=datas[3], alt=datas[4])

                if caller == 'FL' : # Less informations with FLiT3r
                    calls[hash]['VC']['FORMAT'][caller] = datas[HEADER["FORMAT"]]
                    calls[hash]['VC']['SAMPLE'][caller] = datas[HEADER["SAMPLE"]]
                else:
                    calls[hash]['VC']['FILTER'][caller] = datas[HEADER['FILTER']]
                    calls[hash]['VC']['INFO'][caller] = datas[HEADER['INFO']]
                    calls[hash]['VC']['FORMAT'][caller] = datas[HEADER["FORMAT"]]
                    calls[hash]['VC']['SAMPLE'][caller] = datas[HEADER["SAMPLE"]]

                    GT = datas[9].split(':')[0]

                    if caller == "VD":
                        # Save genotype to raise warning if all callers don't report same GT.
                        # First manage case when Vardict return 1/0 instead of 0/1
                        
                        if GT == '1/0':
                            GT = '0/1'

                    calls[hash]['VC']['GT'][caller] = GT

                calls[hash]['VC']['VAF'][caller] = 100 * vcfs[caller]["vcf"].VAF(datas)

                # Store TRC, ARC and RRC for each VC
                calls[hash]['VC']['TRC'][caller] = vcfs[caller]["vcf"].depth(datas)

                if caller in ['VS','VD','BT']:

                    calls[hash]['VC']['ARC'][caller], calls[hash]['VC']['ARC+'][caller], calls[hash]['VC']['ARC-'][caller] = vcfs[caller]["vcf"].arc(datas)

                    calls[hash]['VC']['RRC'][caller], calls[hash]['VC']['RRC+'][caller], calls[hash]['VC']['RRC-'][caller] = vcfs[caller]["vcf"].rrc(datas)

                    calls[hash]['VC'][f"TRC+"][caller] = calls[hash]['VC']['ARC+'][caller] + calls[hash]['VC']['RRC+'][caller]
                    calls[hash]['VC'][f"TRC-"][caller] = calls[hash]['VC']['ARC-'][caller] + calls[hash]['VC']['RRC-'][caller]

                    # for id, strand in enumerate(['+','-'], start=1):

                    #     calls[hash]['VC'][f'ARC{strand}'][caller] = calls[hash]['VC']['ARC'][caller][id]
                    #     calls[hash]['VC'][f'RRC{strand}'][caller] = calls[hash]['VC']['RRC'][caller][id]

                    #     calls[hash]['VC'][f"TRC{strand}"][caller] = (
                    #         calls[hash]['VC']['RRC'][caller][id] +
                    #         calls[hash]['VC']['ARC'][caller][id]
                    #     )

                else:

                    calls[hash]['VC']['ARC'][caller] = vcfs[caller]["vcf"].arc(datas)

                    calls[hash]['VC']['RRC'][caller] = (
                        calls[hash]['VC']['TRC'][caller] -
                        calls[hash]['VC']['ARC'][caller]
                    )

                # TBI (Must be done only once)
                # 2/ revert vt norm of variants that should not have been left-aligned (if any)
                if calls[hash]['VT'] == 'DEL':

                    for caller in calls[hash]['VC']['VAF']:

                        ## from vt decompose
                        if 'OLD_MULTIALLELIC' in calls[hash]['VC']['INFO'][caller]:
                            continue

                        # from vt decompose
                        if 'OLD_VARIANT' in calls[hash]['VC']['INFO'][caller]:
                            info_old_variant = calls[hash]['VC']['INFO'][caller] \
                                            .split('OLD_VARIANT=')[1] \
                                            .split(',')
                            for old_variant in info_old_variant:
                                ov_ref, ov_alt = old_variant.split('/')
                                # DEL is not parsimonious eg. POS:ABB/AB (should be POS+1:BB/B)
                                if len(ov_alt) != 1:
                                    continue
                                # Keep trace of old/new DEL descriptor
                                new_hash: str = sha256(
                                    string=f"{(datas[0]).removeprefix('chr')}:{datas[1]}:{ov_ref}:{'|'.join(ov_alt)}".encode()
                                ).hexdigest()
                                if not hash in rev:
                                    rev[hash] = []
                                rev[hash].append(new_hash)
                                # DEL is parsimonious;
                                # clone <hash_key> according to <OLD_VARIANT> description
                                # call[new_hash] = callsets[hash].copy()
                                # Update variant
                                variant = calls.pop(hash)
                                variant["VC"]["REF"] = ov_ref
                                variant["VC"]["ALT"] = ov_alt
                                calls[new_hash] = variant
                                # Sync hash value
                                hash = new_hash
                            break

                # ETBI

                if not (calls[hash]["VT"] in ['INV','MNV','CSV'] or "FL" in calls[hash]["VC"]["VAF"]):

                    alt = calls[hash]["VC"]["ALT"]

                    position = datas[1]

                    # reset <ALT> allele descriptor for <INS> and <DEL>
                    # A/AA = insA; AA/A = delA
                    if calls[hash]['VT'] == 'INS':

                        alt = calls[hash]["VC"]["ALT"][1:]

                    elif calls[hash]['VT'] == 'DEL':

                        alt = calls[hash]["VC"]["REF"][1:]
                        # reset <POS> if <DEL>
                        # delA = p+1 in read pile
                        position: int = str(int(datas[1]) + 1)
                                

                    # 1/ link the pileup entry covering dictionary,
                    # and then add the elements from the new dictionary to it. the variant call
                    # Link the variant to a Pileup entry
                    position_index: str = sha256(
                                        string=f"{(datas[0]).removeprefix('chr')}:{position}".encode()
                                    ).hexdigest()
                        
                    if not position_index in readpile:
                            
                        readpile[position_index] = {}

                    readpile[position_index][hash] = alt
                    
                elif calls[hash]["VT"] in ['INV','MNV','CSV']:

                    if not hash in INV_MNV_CSV:

                        INV_MNV_CSV[hash] = calls.pop(hash)
                    
                elif "FL" in calls[hash]["VC"]["VAF"]:

                    if not hash in FLiT3r:

                        FLiT3r[hash] = calls.pop(hash)

    # ===========================================================================================
    # Process exceptions without Pileup : INV,MNV and CSV
    # ===========================================================================================
    INV_MNV_CSV = functions.process_without_pileup(INV_MNV_CSV,thresholds,SBM,params.sbm_homozygous)

    # ===========================================================================================
    # Process FL output without Pileup :
    # ===========================================================================================
    FLiT3r = functions.process_without_pileup(FLiT3r,thresholds,SBM,params.sbm_homozygous)

    # ===========================================================================================
    # Process the rest of variants with Pileup
    # ===========================================================================================

    # --------------------------------
    # process pileup data
    # --------------------------------

    #####################################################################################################################

    with open(pileup.get_path(), mode='r') as pileup:

        # Get header, make it upper case, split it on <tab>,
        # take only 3 first letter (so barcode become BAR)
        # Save the column POSITION for each header part (BAR is first column, CR second...)
        # This will be used for reading lines after
        header = pileup.readline().strip('\n').upper()
        header_dict = {}

        for id, header_part in enumerate(header.split('\t')):
            header_dict[header_part[:3]] = id

        for n, line in enumerate(pileup, start=1):

            datas = line.strip('\n').split('\t')

            # Create a key with format chromosome:POSITION

            position_index: str = sha256(
                                    string=f"{(datas[header_dict['CHR']]).removeprefix('chr')}:{datas[header_dict['POS']]}".encode()
                                ).hexdigest()
            
            # print(position_index in readpile)

            # Check if variant is reported in one of the VCF files
            if (datas[header_dict['REF']] != 'N') and (position_index in readpile):
                
                for key in readpile[position_index]:
                    # depth of coverage at <CHR:POS> = sum(Nt) + #DEL (if any)
                    coverage = {"plus": 0,
                                "minus": 0,
                                "total": 0}
                    plus_strand_POSITIONs = [0,2,4,6]
                    minus_strand_POSITIONs = [1,3,5,7]

                    for column, value in enumerate(datas[5:13], start=0):
                        try:
                            coverage['total'] += int(value)
                            if column in plus_strand_POSITIONs:
                                coverage['plus'] += int(value)
                            elif column in minus_strand_POSITIONs:
                                coverage["minus"] += int(value)
                        except ValueError:
                            logger.warning("Uknown coverage value present in pileup file.")
                            logger.warning(f"Warning was raised by: {value} at line {n} column {column}.")

                    if not 'final_metrics' in calls[key]:
                        calls[key]['final_metrics'] = {}

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
                                coverage += (int(del_cov1) + int(del_cov2))

                    calls[key]['final_metrics']['TRC'] = coverage["total"]
                    calls[key]['final_metrics']['TRC-'] = coverage['minus']
                    calls[key]['final_metrics']['TRC+'] = coverage['plus']

                    # Add <REF> strand-specific counts for each <ALT> at <CHR:POS>
                    # ie. the matched key in callset
                    # RRC : Reference Read Counts
                    # ARC : Alternative Read Counts
                    for strand in ['-','+']:

                        calls[key]['final_metrics'][f"RRC{strand}"] = datas[header_dict[datas[header_dict['REF']] + strand]]

                    # if <ALT> is an <INDEL>
                    VARIANT_COUNTS = 0

                    if (calls[key]['VT'] == 'DEL') or (calls[key]['VT'] == 'INS'):
                    
                        if re.search(r"\b" + readpile[position_index][key] + r"\b:[0-9]+,[0-9]+" ,\
                                    (datas[header_dict[calls[key]['VT']]])):
                            
                            count1, count2 = re.findall(
                                r"\b" + readpile[position_index][key] + r"\b:[0-9]+,[0-9]+",
                                datas[header_dict[calls[key]['VT']]]
                                )[0].split(':')[-1].split(',')
                            
                            VARIANT_COUNTS = int(count1) + int(count2)
                            
                            #Save indel count per strand in callset ARC+ and ARC-
                            for strand in ['+','-']:
                                ARC_STRAND = 'ARC' + strand
                                if not ARC_STRAND in calls[key]['final_metrics']:
                                    calls[key]['final_metrics'][ARC_STRAND] = ''
                                if strand == '+':
                                    calls[key]['final_metrics'][ARC_STRAND] = count1
                                else:
                                    calls[key]['final_metrics'][ARC_STRAND] = count2

                    # any other <VarType>
                    else:
                        # keep 1st character of <ALT> string (approx. counts for non-SNV variants)
                        variant_first_char = readpile[position_index][key][0]
                        if variant_first_char == 'N':
                            del calls[key]
                            continue
                        for strand in ['+','-']:
                            ARC_STRAND = 'ARC' + strand
                            if not ARC_STRAND in calls[key]['final_metrics']:
                                calls[key]['final_metrics'][ARC_STRAND] = ''
                            calls[key]['final_metrics'][ARC_STRAND] = \
                                datas[header_dict[variant_first_char + strand]]
                            VARIANT_COUNTS += int(
                                datas[header_dict[variant_first_char + strand]]
                            )
                    
                    # --------------------------------------------------------------
                    # <ALT> not covered in read pile;
                    # keep trace of rejected calls (for test / rescuing purpose)
                    # --------------------------------------------------------------
                    if VARIANT_COUNTS == 0:

                        if (calls[key]['VT'] == 'INS') and len(calls[key]["VC"]["ALT"])-1 > params.length_indels:

                            ITD[key] = calls[key].copy()

                        else:

                            rejected[key] = calls[key].copy()

                        del calls[key]

                        continue

                    # --------------------------------------------------------------
                    # Compute metrics
                    # --------------------------------------------------------------
                    # Get meanVAF from Variant Callers
                    # 231010 : Mean VAF only used withoutPileup -
                    # not usefull to compute meanVAF in pileup
                    #
                    # calls[key]['final_metrics']['VCR'] = format(
                    #     (sum(calls[key]['VC']['VAF'].values())/
                    #     calls[key]['final_metrics']['VCN']),
                    #     '.5f'
                    # )
                    #round((sum(calls[key]['VC']['VAF'].values())
                    #/calls[key]['final_metrics']['VCN']),5)

                    # Save number of Variant Callers that found this variant
                    calls[key]['final_metrics']['VCN'] = len(calls[key]['VC']['VAF'])
                    # keep trace of used VC identifier(s)
                    calls[key]['final_metrics']['VCI'] = ','.join(
                        sorted(calls[key]['VC']['VAF'].keys(), key=str.lower)
                    )

                    #format [R/A]RC for vcf
                    calls[key]['final_metrics']['ARC'],\
                    calls[key]['final_metrics']['RRC'] = \
                        functions.format_rrc_arc(calls,key)

                    # Compute VAF from pileup
                    alt_read_count_plus = calls[key]['final_metrics']['ARC+']
                    alt_read_count_minus = calls[key]['final_metrics']['ARC-']
                    total_read_count = calls[key]['final_metrics']['TRC']
                    calls[key]['final_metrics']['ARR'] = format(
                        (
                            (float(alt_read_count_plus) + float(alt_read_count_minus)) /
                            float(total_read_count)
                        )*100,
                        '.5f'
                    )

                    # Estimate GT
                    calls[key]['final_metrics']['VAR'] = \
                        functions.categorize_variant_type(calls,key,thresholds)
                    calls[key]['final_metrics']['GT'] = \
                        functions.estimate_gt(calls,key)
                    # If variant callers do not agree on genotype change VAR to WAR (warning)
                    calls[key]['final_metrics']['VAR'] = \
                        functions.compare_gt(calls,key)

                    calls[key]['final_metrics']['BRC'],\
                    calls[key]['final_metrics']['BRR'],\
                    calls[key]['final_metrics']['BRE'] = \
                        functions.estimate_brc_r_e(calls,key,datas)

                    calls[key]['final_metrics']['SBP'],calls[key]['SBM'] = \
                        functions.estimate_sbm(calls,key,params.sbm_homozygous)
                    calls[key]['final_metrics']['LOW'] = \
                        functions.vaf_user_threshold(calls,key,thresholds)
                    calls[key]['final_metrics']['BKG'] = \
                        functions.categorize_background_signal(calls,key,thresholds)
                    calls[key]['final_metrics']['FILTER'] = \
                        functions.format_float_descriptors(calls,key,SBM)

                    calls[key]['final_metrics']['PIL'] = 'Y'
                    calls[key]['final_metrics']['RES'] = 'N'

    # ===========================================================================================
    # Rescue INS not identified by pileup (probably ITD or long insertions), mostly coming from PL
    # And remove rescued INS from Trash
    # Do the same with DEL
    # ===========================================================================================

    ITD = functions.process_without_pileup(ITD,thresholds,SBM,params.sbm_homozygous)

    # ===========================================================================================
    # Merge Pileup and noPileup dict
    # ===========================================================================================
    calls = functions.merge_variant_dict ([calls,INV_MNV_CSV,ITD,FLiT3r])

    # ===========================================================================================
    # Clean calls and rejected variants
    # ===========================================================================================
    # Delete variants called outside of read pile (very unlikely)
    key_to_delete = []
    for variant_key in calls:
        if not 'final_metrics' in calls[variant_key]:
            key_to_delete.append(variant_key)
    for bad_key in key_to_delete:
        del calls[bad_key]

    # Select best descriptor for parcimonious DEL
    for variant_key in rev:
        BRR = 1000
        # Check that $hash_key is not in rejected calls
        if not variant_key in calls:
            if variant_key in rejected:
                del rejected[variant_key]
        else:
            BRR = float(calls[variant_key]['final_metrics']['BRR'])
        for new_key in rev[variant_key]:
            NEW_BRR = 1000
            # Check that $hash_key is not in rejected calls
            if not new_key in calls:
                if new_key in rejected:
                    del rejected[new_key]
            elif BRR == 1000:
                NEW_BRR = float(calls[new_key]['final_metrics']['BRR'])
                BRR = NEW_BRR
            else:
                NEW_BRR = float(calls[new_key]['final_metrics']['BRR'])
                if BRR <= NEW_BRR:
                    del calls[new_key]
                else:
                    del calls[variant_key]
                    BRR = NEW_BRR

    # ===========================================================================================
    # Process Rejected dic
    # ===========================================================================================
    rejected = functions.process_without_pileup(rejected,thresholds,SBM,params.sbm_homozygous)

    # ===========================================================================================
    # rescuing rejected calls w/ FILTER <PASS> (optional; if any)
    # ===========================================================================================
    if params.rescue:
        # rescue <PASS> calls only
        # do not rescue SNV (more likely to be VS artefacts)
        # pindel bug where ARC > TRC
        rescued_keys = [
            key for key in rejected if (
                rejected[key]['final_metrics']['FILTER'] == 'PASS' and \
                rejected[key]['VT'] != 'SNV' and \
                float(rejected[key]['final_metrics']['ARR']) <= 100
            )
        ]

        for key in rescued_keys:
            rejected[key]['final_metrics']['RES'] = 'Y'
            calls[key] = rejected[key].copy()
            del rejected[key]

    # ===========================================================================================
    # set vcf headers/fields for output files
    # ===========================================================================================
    # vcf_fields = [
    #     'VAR', 'BKG', 'TRC', 'RRC', 'ARC', 'BRC', 'ARR',
    #     'BRR', 'BRE', 'SBP', 'SBM', 'LOW','VCI','VCN','PIL','RES'
    # ]

    # for variant_key in calls:
    #     # calls[variant_key]['vcf_fields'] = vcf_fields
    #     # Toutes les cles pointent vers la même liste en mémoire....
    #     # ce qui implique que la modif de lune modifie aussi toutes les autres...
    #     # copy.copy() : chaque cle pointe vers une copie disincte de vcf_fields,
    #     # ce qui permet de les modifier sans affecter les autres cles
    #     calls[variant_key]['vcf_fields'] = copy.copy(vcf_fields)
    # for variant_key in rejected:
    #     rejected[variant_key]['vcf_fields'] = copy.copy(vcf_fields)

    # ===========================================================================================
    # Generate the vcf header
    # Different headers for normal, rejected and cartagenia
    # Separate header in part 1, 2 and 3 only because cartagenia
    # has 2 lines different in from normal header
    # Just merge then part one, two and 3 wih normal header, rejected header and cartagenia header
    # ===========================================================================================

    writter = files.GenomicWritter(file=params.output)
    # write rejected calls to file
    # with open(os.path.join(os.getcwd(), f"{params.sample}_failed.vcf"), mode='w',encoding='utf-8') as OUT_TRASH_FILE:
    #     OUT_TRASH_FILE.write(VCF_HEADER)
    #     ordered_variant_key = functions.order_var(rejected, contigs["contig"])
        
    #     for variant_key in ordered_variant_key:
    #         OUT_TRASH_FILE.write(functions.print_var(variant_key, rejected, 'final_metrics'))

    
    writter.writeVCF(contigs=contigs, variants=calls, samples=[params.sample], thresholds=thresholds)
    # write valid calls to file
    # with open(params.output, 'w',encoding='utf-8') as OUT_FILE:
    #     OUT_FILE.write(VCF_HEADER)
    #     for variant_key in functions.order_var(calls, contigs["contig"]):
    #         OUT_FILE.write(functions.print_var(variant_key, calls, 'final_metrics'))
