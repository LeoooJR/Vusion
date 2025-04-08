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

import copy
import errors
from hashlib import sha256
from loguru import logger
import os
import pandas as pd
import re
import sys

sys.path.insert(1, os.path.dirname(os.path.abspath(__file__)))
import utils as functions

# ===========================================================================================
# Initiate variables
# 231010 - CT replaced by final_metrics
# ===========================================================================================
callsets = {}  # dictionary of unique variant calls
VCF_HEADER = ''
thresholds = [10, 30, 40, 60, 70, 80, 20, 30, 50, 1]

known_callers: frozenset = {'ST', 'VS', 'VD', 'PL', 'HC', 'CS', 'HS', 'FL'}
# known variant callers are:
# samtools			(ST)
# varscan			(VS)
# vardict			(VD)
# pindel			(PL)
# haplotypecaller	(HC)
# smCounter2		(SM; DEPRECATED)
# control & hotspot (CS & HS) ## when --hotspot option is set;
# CtlSet & HotSpot should both originate from CombineVCF2final_metrics

def combine(params):

    # Check reference genome index
    try:
        functions.verify_files(file=params.reference, format="fai")
    except errors.FastaIndexError:
        logger.error(f"{params.reference} is not a valid FASTA index.")
        raise SystemExit
    
    logger.success(f"Fasta index {params.reference} has been successfully checked.")
    
    # Check pileup
    try:
        functions.verify_files(file=params.pileup, format="pileup")
    except errors.PileupError:
        logger.error(f"{params.reference} is not a valid PILEUP.")
        raise SystemExit
    
    logger.success(f"Pileup {params.reference} has been successfully checked.")

    # Check VCFs
    vcfs: dict = {}

    for vcf in params.vcfs:

        input = vcf.split(',')

        if len(input) != 2:
            logger.error('Wrong number of argument in --vcf option.')
            logger.error(f'Error was raised by: {input}.')
            raise ValueError("Wrong number of argument in --vcf option.")
        
        try:
            int(input[0])
            logger.error("Wrong type of argument in --vcf option.")
            logger.error(f"Error was raised by: {input}.")
            raise SystemExit("Wrong type of argument in --vcf option.")
        except ValueError:
            if input[0] not in known_callers:
                logger.error(f"{input[0]} caller is not supported in --vcf option.")
                raise SystemExit("Caller not supported in --vcf options.")
            
        try:
            functions.verify_files(file=vcf, format="vcf")
        except errors.VCFError:
            logger.error(f"{vcf} is not a valid VCF.")
            raise SystemExit
        
        vcfs[input[0]] = {"vcf": input[1], "index": None}

        logger.debug(f"Variant Callers inputed: {input[0]}")

    # Load reference genome dict ie. list of ordered contig names / length
    contigs: list[pd.Series] = []

    with open(params.reference, mode='r') as ref:
        logger.debug(f"Parsing {params.reference}")
        for line in ref:
            contigs.append(pd.Series(data=line.strip().split('\t')))

    contigs: pd.DataFrame = pd.DataFrame(data=contigs, columns=["contig", "length", "index", "pbline", "byteline"])

    contigs = contigs.astype({"contig": "string", "length": "uint", "index": "uint", "pbline": "uint", "byteline": "uint"})

    print(contigs.memory_usage(deep=True))
    print(contigs.dtypes)

    # Check if all mandatory option are given and modify variable depending of given options
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
    }

    if params.disable_strand_bias:
        SBM = 2

    if params.thresholds:

        thresholds = params.thresholds.split(',')

        # Check that we have 10 values in threshold
        if len(thresholds) != 10:
            logger.error(f"Invalid number of values in --threshold option.")
            raise SystemExit
        
        #Check that all values can be converted to floats
        try:
            thresholds = list(map(float, thresholds)) 
        except ValueError:
            logger.error("Invalid values in threshold option.")
            raise SystemExit
        
        if not any(list(map(lambda value: value >= MIN_THRESHOLD and value <= MAX_THRESHOLD, thresholds))):
            logger.error("Option --thresholds values cannot be equal or higher than 100, or lower than 0.")
            raise SystemExit

        # Check that 6 first given threshold are unique
        if len(thresholds[0:6]) != len(set(thresholds[0:6])):
            logger.error("Option --thresholds six first values must be unique.")
            raise SystemExit

        # Check that second group of threshold values are unique
        if len(thresholds[6:9]) != len(set(thresholds[6:9])):
            logger.error("Option --thresholds values 7, 8 and 9 must be unique.")
            raise SystemExit
        
        thresholds[0:6] = sorted(thresholds[0:6])
        thresholds[6:9] = sorted(thresholds[6:9])

        logger.debug(f"Threshold: {thresholds}")

    # rejected variants (if any; will be written to $opt{trash_file})
    rejected = {}
    rev = {}
    INV_MNV_CSV = {}
    FLiT3r = {}
    readpile = {}

    for caller in vcfs:
        # Check caller and set the subparse function to the corresponding function to get VAF
        if caller == 'ST':
            subparse = functions.samtools_vaf
            TRC = functions.st_coverage
            ARC = functions.st_arc
            RRC = functions.st_rrc
        if caller == 'VS':
            subparse = functions.varscan_vaf
            TRC = functions.varscan_coverage
            ARC = functions.varscan_arc
            RRC = functions.varscan_rrc
        if caller == 'VD':
            subparse = functions.vardict_vaf
            TRC = functions.vardict_coverage
            ARC = functions.vardict_arc
            RRC = functions.vardict_rrc
        if caller == 'PL':
            subparse = functions.pindel_vaf
            TRC = functions.pindel_coverage
            ARC = functions.pindel_arc
        if caller == 'HC':
            subparse = functions.haplotypecaller_vaf
            TRC = functions.haplotypecaller_coverage
            ARC = functions.haplotypecaller_arc
        if caller == 'FL':
            subparse = functions.FLT3R_vaf
            TRC = functions.FLT3_coverage
            ARC = functions.FLT3_arc

        with open(vcfs[caller]["vcf"], mode='r') as vcf:

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

                # Create a variant unique id, with chromosome, POSITION,
                # reference and alternative allele
                hash: str = sha256(
                            string=f"{(datas[0]).removeprefix('chr')}:{datas[1]}:{datas[3]}:{'|'.join(datas[4])}".encode()
                        ).hexdigest()

                call: dict = {}

                call['VC'] = {
                    'VAF': {},
                    'GT': {},
                    'FILTER': {},
                    'INFO': {},
                    'FORMAT': {},
                    'SAMPLE': {},
                    'TRC': {},
                    'ARC':{},
                    'RRC':{}
                }

                if caller == 'FL' : # Less informations with FLiT3r
                    call['VC']['FORMAT'][caller] = datas[HEADER["FORMAT"]]
                    call['VC']['SAMPLE'][caller] = datas[HEADER["SAMPLE"]]
                else:
                    call['VC']['FILTER'][caller] = datas[HEADER['FILTER']]
                    call['VC']['INFO'][caller] = datas[HEADER['INFO']]
                    call['VC']['FORMAT'][caller] = datas[HEADER["FORMAT"]]
                    call['VC']['SAMPLE'][caller] = datas[HEADER["SAMPLE"]]

                    # Save genotype to raise warning if all callers don't report same GT.
                    # First manage case when Vardict return 1/0 instead of 0/1
                    # callsets[hash]['VC']['GT'] = []
                    # GT = callsets[hash][caller]['SAMPLE'].split(':')[0]
                    GT = datas[9].split(':')[0]
                    if GT == '1/0':
                        GT = '0/1'
                    #callsets[hash]['VC']['GT'].append(GT)
                    call['VC']['GT'][caller] = GT

                call['VC']['VAF'][caller] = 100 * subparse(datas)

                # Store TRC, ARC and RRC for each VC
                call['VC']['TRC'][caller] = TRC(datas)
                if caller in ['VS','VD','ST']:
                    call['VC']['ARC'][caller] = ARC(datas)[0]
                    call['VC']['RRC'][caller] = RRC(datas)[0]
                    VALUE_STRAND = ['TRC-','TRC+','ARC+','ARC-','RRC+','RRC-']
                    for value in VALUE_STRAND:
                        if value not in call['VC']:
                            call['VC'][value] = {}
                    i = 1
                    for strand in ['+','-']:
                        call['VC']['TRC'+strand][caller] = (
                            RRC(datas)[i] +
                            ARC(datas)[i]
                        )
                        call['VC']['ARC'+strand][caller] = ARC(datas)[i]
                        call['VC']['RRC'+strand][caller] = RRC(datas)[i]
                        i +=1
                else:
                    call['VC']['ARC'][caller] = ARC(datas)
                    call['VC']['RRC'][caller] = (
                        float(TRC(datas))-
                        float(ARC(datas))
                    )

                ## 1/ define <VT>
                call['VT'] = functions.define_variant_type(ref=datas[3], alt=datas[4])

                # 2/ revert vt norm of variants that should not have been left-aligned (if any)
                if call['VT'] == 'DEL':
                    for caller in call['VC']['VAF']:
                        ## from vt decompose
                        if 'OLD_MULTIALLELIC' in call['VC']['INFO'][caller]:
                            continue
                        # from vt decompose
                        if 'OLD_VARIANT' in call['VC']['INFO'][caller]:
                            info_old_variant = call['VC']['INFO'][caller] \
                                            .split('OLD_VARIANT=')[1] \
                                            .split(',')
                            for old_variant in info_old_variant:
                                ov_ref, ov_alt = old_variant.split('/')
                                # DEL is not parsimonious eg. POS:ABB/AB (should be POS+1:BB/B)
                                if len(ov_alt) != 1:
                                    continue
                                # Keep trace of old/new DEL descriptor
                                hash: str = sha256(
                                    string=f"{(datas[0]).removeprefix('chr')}:{datas[1]}:{ov_ref}:{'|'.join(ov_alt)}".encode()
                                ).hexdigest()
                                if not hash in rev:
                                    rev[hash] = []
                                rev[hash].append(hash)
                                # DEL is parsimonious;
                                # clone <hash_key> according to <OLD_VARIANT> description
                                # call[new_hash] = callsets[hash].copy()
                            break

                if not (call["VT"] in ['INV','MNV','CSV'] or "FL" in call["VC"]["VAF"]):

                    if not hash in callsets:

                        callsets[hash] = call

                        variant_info = {}
                        variant_info['CHROM'], variant_info['POS'], \
                        variant_info['REF'], variant_info['ALT']  = key.split(':')

                        # reset <ALT> allele descriptor for <INS> and <DEL>
                        # A/AA = insA; AA/A = delA
                        if callsets[hash]['VT'] == 'INS':
                            var: str = 'ALT'
                        elif callsets[hash]['VT'] == 'DEL':
                            var: str = 'REF'
                        else:
                            var: str = ''

                        if var:
                            variant_info['ALT'] = variant_info[var][1:]

                        # reset <POS> if <DEL>
                        # delA = p+1 in read pile
                        if var == 'REF':
                            position: int = str(int(datas[1]) + 1)

                        # 1/ link the pileup entry covering dictionary,
                        # and then add the elements from the new dictionary to it. the variant call
                        position_index: str = sha256(
                                            string=f"{(datas[0]).removeprefix('chr')}:{position}".encode()
                                        ).hexdigest()
                        
                        if not position_index in readpile:
                            readpile[position_index] = {}

                        readpile[position_index][hash] = datas[4]
                    
                elif call["VT"] in ['INV','MNV','CSV']:

                    if not hash in INV_MNV_CSV:

                        INV_MNV_CSV[hash] = call
                    
                elif "FL" in call["VC"]["VAF"]:

                    if not hash in FLiT3r:

                        FLiT3r[hash] = call

    # ===========================================================================================
    # Process exceptions without Pileup : INV,MNV and CSV
    # ===========================================================================================
    # INV_MNV_CSV = {}
    # for key, value in list(callsets.items()):
    #     if value['VT'] in ['INV','MNV','CSV']:
    #         INV_MNV_CSV[key] = value
    #         del callsets[key]
    # INV_MNV_CSV = {key: callsets[key] for key in callsets if callsets[key]["VT"] in ['INV','MNV','CSV']}
    INV_MNV_CSV = functions.process_without_pileup(INV_MNV_CSV,thresholds,SBM,params.sbm_homozygous)

    # ===========================================================================================
    # Process FL output without Pileup :
    # ===========================================================================================
    # FLiT3r = {}
    # for key, value in list(callsets.items()):
    #     VCI = sorted(value['VC']['VAF'].keys(),key=str.lower)
    #     if 'FL' in VCI:
    #         FLiT3r[key] = value
    #         del callsets[key]

    # FLiT3r = {key: callsets[key] for key in callsets if "FL" in callsets[hash]["VC"]["VAF"]}
    FLiT3r = functions.process_without_pileup(FLiT3r,thresholds,SBM,params.sbm_homozygous)

    # ===========================================================================================
    # Process the rest of variants with Pileup
    # ===========================================================================================
    # readpile = {}

    # ---------------------------------
    # 1/ define pileup data to process
    # ---------------------------------
    # for key in callsets:
    #     variant_info = {}
    #     variant_info['CHROM'], variant_info['POS'], \
    #     variant_info['REF'], variant_info['ALT']  = key.split(':')

    #     # reset <ALT> allele descriptor for <INS> and <DEL>
    #     # A/AA = insA; AA/A = delA
    #     if callsets[key]['VT'] == 'INS':
    #         variant_info['VAR'] = 'ALT'
    #     elif callsets[key]['VT'] == 'DEL':
    #         variant_info['VAR'] = 'REF'
    #     else:
    #         variant_info['VAR'] = ''
    #     if variant_info['VAR'] != '':
    #         variant_info['ALT'] = variant_info[variant_info['VAR']][1:]

    #     # reset <POS> if <DEL>
    #     # delA = p+1 in read pile
    #     if variant_info['VAR'] == 'REF':
    #         variant_info['POS'] = str(int(variant_info['POS']) + 1)

    #     # 1/ link the pileup entry covering dictionary,
    #     # and then add the elements from the new dictionary to it. the variant call
    #     TMP_KEY: str = sha256(
    #                         string=f"{(datas[0]).removeprefix('chr')}:{datas[1]}".encode()
    #                     ).hexdigest()
    #     if not TMP_KEY in readpile:
    #         readpile[TMP_KEY] = {}
    #     readpile[TMP_KEY][key] = variant_info['ALT']

    # --------------------------------
    # process pileup data
    # --------------------------------

    with open(params.pileup, mode='r') as pileup:
        # Get header, make it upper case, split it on <tab>,
        # take only 3 first letter (so barcode become BAR)
        # Save the column POSITION for each header part (BAR is first column, CR second...)
        # This will be used for reading lines after
        header = pileup.readline().strip('\n').upper()
        header_dict = {}
        POSITION = 0

        for header_part in header.split('\t'):
            header_dict[header_part[:3]] = POSITION
            POSITION += 1

        for n, line in enumerate(pileup, start=1):

            datas = line.strip('\n').split('\t')

            # Create a key with format chromosome:POSITION
            # PILEUP_KEY = ':'.join([
            #     datas[header_dict['CHR']],
            #     datas[header_dict['POS']]
            # ])

            position_index: str = sha256(
                                    string=f"{(datas[header_dict['CHR']]).removeprefix('chr')}:{datas[header_dict['POS']]}".encode()
                                ).hexdigest()

            # Check if variant is reported in one of the VCF files
            if (datas[header_dict['REF']] != 'N') and (position_index in readpile):
                
                for key in readpile[position_index]:
                    #a = 0
                    # depth of coverage at <CHR:POS> = sum(Nt) + #DEL (if any)
                    coverage = {"plus": 0,
                                "minus": 0,
                                "total": 0}
                    coverage_plus_strand = 0
                    COVERAGE_MINUS = 0
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

                    if not 'final_metrics' in callsets[key]:
                        callsets[key]['final_metrics'] = {}

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

                    callsets[key]['final_metrics']['TRC'] = coverage["total"]
                    callsets[key]['final_metrics']['TRC-'] = coverage['minus']
                    callsets[key]['final_metrics']['TRC+'] = coverage['plus']

                    # Add <REF> strand-specific counts for each <ALT> at <CHR:POS>
                    # ie. the matched key in callset
                    # RRC : Reference Read Counts
                    # ARC : Alternative Read Counts
                    for strand in ['-','+']:

                        RRC_STRAND = f"RRC{strand}"

                        callsets[key]['final_metrics'][RRC_STRAND] = datas[header_dict[datas[header_dict['REF']] + strand]]

                    # if <ALT> is an <INDEL>
                    VARIANT_COUNTS = 0

                    if (callsets[key]['VT'] == 'DEL') or (callsets[key]['VT'] == 'INS'):
                    
                        if re.search(r"\b" + readpile[position_index][key] + r"\b:[0-9]+,[0-9]+" ,\
                                    (datas[header_dict[callsets[key]['VT']]])):
                            
                            count1, count2 = re.findall(
                                r"\b" + readpile[position_index][key] + r"\b:[0-9]+,[0-9]+",
                                datas[header_dict[callsets[key]['VT']]]
                                )[0].split(':')[-1].split(',')
                            
                            VARIANT_COUNTS = int(count1) + int(count2)
                            
                            #Save indel count per strand in callset ARC+ and ARC-
                            for strand in ['+','-']:
                                ARC_STRAND = 'ARC' + strand
                                if not ARC_STRAND in callsets[key]['final_metrics']:
                                    callsets[key]['final_metrics'][ARC_STRAND] = ''
                                if strand == '+':
                                    callsets[key]['final_metrics'][ARC_STRAND] = count1
                                else:
                                    callsets[key]['final_metrics'][ARC_STRAND] = count2

                    # any other <VarType>
                    else:
                        # keep 1st character of <ALT> string (approx. counts for non-SNV variants)
                        variant_first_char = readpile[position_index][key][0]
                        if variant_first_char == 'N':
                            del callsets[key]
                            continue
                        for strand in ['+','-']:
                            ARC_STRAND = 'ARC' + strand
                            if not ARC_STRAND in callsets[key]['final_metrics']:
                                callsets[key]['final_metrics'][ARC_STRAND] = ''
                            callsets[key]['final_metrics'][ARC_STRAND] = \
                                datas[header_dict[variant_first_char + strand]]
                            VARIANT_COUNTS += int(
                                datas[header_dict[variant_first_char + strand]]
                            )

                    # --------------------------------------------------------------
                    # <ALT> not covered in read pile;
                    # keep trace of rejected calls (for test / rescuing purpose)
                    # --------------------------------------------------------------
                    if VARIANT_COUNTS == 0:
                        rejected[key] = callsets[key].copy()
                        del callsets[key]
                        continue

                    # --------------------------------------------------------------
                    # Compute metrics
                    # --------------------------------------------------------------
                    # Get meanVAF from Variant Callers
                    # 231010 : Mean VAF only used withoutPileup -
                    # not usefull to compute meanVAF in pileup
                    #
                    # callsets[key]['final_metrics']['VCR'] = format(
                    #     (sum(callsets[key]['VC']['VAF'].values())/
                    #     callsets[key]['final_metrics']['VCN']),
                    #     '.5f'
                    # )
                    #round((sum(callsets[key]['VC']['VAF'].values())
                    #/callsets[key]['final_metrics']['VCN']),5)

                    # Save number of Variant Callers that found this variant
                    callsets[key]['final_metrics']['VCN'] = len(callsets[key]['VC']['VAF'])
                    # keep trace of used VC identifier(s)
                    callsets[key]['final_metrics']['VCI'] = ','.join(
                        sorted(callsets[key]['VC']['VAF'].keys(), key=str.lower)
                    )

                    #format [R/A]RC for vcf
                    callsets[key]['final_metrics']['ARC'],\
                    callsets[key]['final_metrics']['RRC'] = \
                        functions.format_rrc_arc(callsets,key)

                    # Compute VAF from pileup
                    alt_read_count_plus = callsets[key]['final_metrics']['ARC+']
                    alt_read_count_minus = callsets[key]['final_metrics']['ARC-']
                    total_read_count = callsets[key]['final_metrics']['TRC']
                    callsets[key]['final_metrics']['ARR'] = format(
                        (
                            (float(alt_read_count_plus) + float(alt_read_count_minus)) /
                            float(total_read_count)
                        )*100,
                        '.5f'
                    )

                    # Estimate GT
                    callsets[key]['final_metrics']['VAR'] = \
                        functions.categorize_variant_type(callsets,key,thresholds)
                    callsets[key]['final_metrics']['GT'] = \
                        functions.estimate_gt(callsets,key)
                    # If variant callers do not agree on genotype change VAR to WAR (warning)
                    callsets[key]['final_metrics']['VAR'] = \
                        functions.compare_gt(callsets,key)

                    callsets[key]['final_metrics']['BRC'],\
                    callsets[key]['final_metrics']['BRR'],\
                    callsets[key]['final_metrics']['BRE'] = \
                        functions.estimate_brc_r_e(callsets,key,datas)

                    callsets[key]['final_metrics']['SBP'],callsets[key]['SBM'] = \
                        functions.estimate_sbm(callsets,key,params.sbm_homozygous)
                    callsets[key]['final_metrics']['LOW'] = \
                        functions.vaf_user_threshold(callsets,key,thresholds)
                    callsets[key]['final_metrics']['BKG'] = \
                        functions.categorize_background_signal(callsets,key,thresholds)
                    callsets[key]['final_metrics']['FILTER'] = \
                        functions.format_float_descriptors(callsets,key,SBM)

                    callsets[key]['final_metrics']['PIL'] = 'Y'
                    callsets[key]['final_metrics']['RES'] = 'N'

#####################################################################################################################

        # ===========================================================================================
        # Rescue INS not identified by pileup (probably ITD or long insertions), mostly coming from PL
        # And remove rescued INS from Trash
        # Do the same with DEL
        # ===========================================================================================
        ITD = {}
        DEL_noPileup = {}
        for key, value in list(rejected.items()):
            if value['VT']=='INS' and len(key.split(':')[3])-1 > params.len_delins:
                # if INS == params.len_delins : likely a PL artefact, do not rescue
                ITD[key] = value
                del rejected[key]

        ITD = functions.process_without_pileup(ITD,thresholds,SBM,params.sbm_homozygous)

        # ===========================================================================================
        # Merge Pileup and noPileup dict
        # ===========================================================================================
        # print("callsets : "+ str(len(callsets)))
        # print("INV_MNV_CSV : "+ str(len(INV_MNV_CSV)))
        # print("ITD : "+ str(len(ITD)))
        callsets = functions.merge_variant_dict ([callsets,INV_MNV_CSV,ITD,FLiT3r])
        # print("callsets final : "+ str(len(callsets)))

        # ===========================================================================================
        # Clean callsets and rejected variants
        # ===========================================================================================
        # Delete variants called outside of read pile (very unlikely)
        key_to_delete = []
        for variant_key in callsets:
            if not 'final_metrics' in callsets[variant_key]:
                key_to_delete.append(variant_key)
        for bad_key in key_to_delete:
            del callsets[bad_key]

        # Select best descriptor for parcimonious DEL
        for variant_key in rev:
            BRR = 1000
            # Check that $hash_key is not in rejected calls
            if not variant_key in callsets:
                if variant_key in rejected:
                    del rejected[variant_key]
            else:
                BRR = float(callsets[variant_key]['final_metrics']['BRR'])
            for new_key in rev[variant_key]:
                NEW_BRR = 1000
                # Check that $hash_key is not in rejected calls
                if not new_key in callsets:
                    if new_key in rejected:
                        del rejected[new_key]
                elif BRR == 1000:
                    NEW_BRR = float(callsets[new_key]['final_metrics']['BRR'])
                    BRR = NEW_BRR
                else:
                    NEW_BRR = float(callsets[new_key]['final_metrics']['BRR'])
                    if BRR <= NEW_BRR:
                        del callsets[new_key]
                    else:
                        del callsets[variant_key]
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
                callsets[key] = rejected[key].copy()
                del rejected[key]


        # ===========================================================================================
        # set vcf headers/fields for output files
        # ===========================================================================================
        vcf_fields = [
            'VAR', 'BKG', 'TRC', 'RRC', 'ARC', 'BRC', 'ARR',
            'BRR', 'BRE', 'SBP', 'SBM', 'LOW','VCI','VCN','PIL','RES'
        ]

        for variant_key in callsets:
            # callsets[variant_key]['vcf_fields'] = vcf_fields
            # Toutes les cles pointent vers la même liste en mémoire....
            # ce qui implique que la modif de lune modifie aussi toutes les autres...
            # copy.copy() : chaque cle pointe vers une copie disincte de vcf_fields,
            # ce qui permet de les modifier sans affecter les autres cles
            callsets[variant_key]['vcf_fields'] = copy.copy(vcf_fields)
        for variant_key in rejected:
            rejected[variant_key]['vcf_fields'] = copy.copy(vcf_fields)



        # ===========================================================================================
        # Generate the vcf header
        # Different headers for normal, rejected and cartagenia
        # Separate header in part 1, 2 and 3 only because cartagenia
        # has 2 lines different in from normal header
        # Just merge then part one, two and 3 wih normal header, rejected header and cartagenia header
        # ===========================================================================================
        VCF_HEADER_PART1 = '##fileformat=VCFv4.0\n'
        for contig in contigs:
            contig_id = contigs[contig][0]
            contig_length = contigs[contig][1]
            VCF_HEADER_PART1 += '##contig=<ID=' + contig_id + ',length=' + contig_length + '>\n'

        VCF_HEADER_PART1 += '##FILTER=<ID=PASS,Description="All filters passed">\n'
        VCF_HEADER_PART1 += '##FILTER=<ID=FAIL,Description="At least one filter failed">\n'
        VCF_HEADER_PART1 += '##INFO=<ID=VAR,Number=1,Type=String,Description="Type of the variant '+ \
        'ie. SNV, MNV, INS, DEL, INV or CSV">\n'
        VCF_HEADER_PART1 += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        VCF_HEADER_PART1 += '##FORMAT=<ID=VAR,Number=1,Type=String,Description="ALT is likely (L) '+\
        'or probably (P) subclonal (SC), heterozygous (HE) or homozygous (HO)">\n'
        VCF_HEADER_PART1 += '##FORMAT=<ID=BKG,Number=1,Type=String,Description="ALT is likely (L) ' +\
        'or probably (P) clean (CL) or noisy (NO)">\n'
        VCF_HEADER_PART1 += '##FORMAT=<ID=TRC,Number=1,Type=Integer,Description="Total #reads covering the variant site">\n'

        # ARC and RRC format is different in cartagenia
        VCF_HEADER_PART2 = '##FORMAT=<ID=RRC,Number=3,Type=Integer,Description="#reads supporting the REF allele ' + \
        'onto the (+) and (-) strand, and total #REF reads">\n'
        VCF_HEADER_PART2 += '##FORMAT=<ID=ARC,Number=3,Type=Integer,Description="#reads supporting the ALT allele ' + \
        'onto the (+) and (-) strand, and total #ALT reads">\n'
        VCF_HEADER_CARTAGENIA_PART2 = '##FORMAT=<ID=RRC,Number=1,Type=Integer,Description="#reads '+ \
        'supporting the REF allele">\n'
        VCF_HEADER_CARTAGENIA_PART2 += '##FORMAT=<ID=ARC,Number=1,Type=Integer,Description="#reads ' +\
        'supporting the ALT allele">\n'

        VCF_HEADER_PART3 = '##FORMAT=<ID=BRC,Number=1,Type=Integer,Description="#reads ' +\
        'supporting neither REF nor ALT alleles">\n'
        VCF_HEADER_PART3 += '##FORMAT=<ID=ARR,Number=1,Type=Float,Description="ALT allele ratio (%)">\n'
        VCF_HEADER_PART3 += '##FORMAT=<ID=BRR,Number=1,Type=Float,Description="Background '+ \
        'read ratio ie. BRC/TRC (%)">\n'
        VCF_HEADER_PART3 += '##FORMAT=<ID=BRE,Number=1,Type=Float,Description="Background ' +\
        'read enrichment ie. BRR/(BRR+ARR) (%)">\n'
        VCF_HEADER_PART3 += '##FORMAT=<ID=SBP,Number=1,Type=Float,Description="Estimated strand bias p-value">\n'
        VCF_HEADER_PART3 += '##FORMAT=<ID=SBM,Number=1,Type=Float,Description="Estimated strand bias TS-metric">\n'
        VCF_HEADER_PART3 += '##FORMAT=<ID=LOW,Number=1,Type=Integer,Description="ARR is GE '+\
        '(resp. LT) ' + str(thresholds[-1]) + '% (0 resp. 1)">\n'

        VCF_HEADER_PART4 = '##FORMAT=<ID=VCI,Number=1,Type=String,Description="VC identifier(s) eg. ST, VS, etc.">\n'
        VCF_HEADER_PART4 += '##FORMAT=<ID=VCN,Number=1,Type=Integer,Description="#VC having called the variant">\n'
        VCF_HEADER_PART4 += '##FORMAT=<ID=PIL,Number=1,Type=String,Description="Describe if the variant '+ \
        'has been found in the pileup or not">\n'
        VCF_HEADER_PART4 += '##FORMAT=<ID=RES,Number=1,Type=String,Description="Describe if the variant has been rescued or not">\n'

        VCF_HEADER_PART5 = '\t'.join([
            '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
            params.sample_id]) + '\n'

        # rejected_header = VCF_HEADER_PART1 + \
        #'##FORMAT=<ID=VCR,Number=1,Type=Float,Description="ALT allele ratio (%)">\n'
        # rejected_header += '##FORMAT=<ID=VCN,Number=1,Type=Integer,Description=" \
        ##VC having called the variant">\n'
        # rejected_header += '##FORMAT=<ID=VCI,Number=.,Type=String,Description=" \
        #VC identifier(s) eg. ST, VS, etc.">\n'
        # rejected_header += '\t'.join([
        #     '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
        #     params.sample_id]) + '\n'
        VCF_HEADER = VCF_HEADER_PART1 + VCF_HEADER_PART2 + VCF_HEADER_PART3 + VCF_HEADER_PART4 + VCF_HEADER_PART5
        VCF_HEADER_cartagenia = VCF_HEADER_PART1 + VCF_HEADER_CARTAGENIA_PART2 + VCF_HEADER_PART3 + VCF_HEADER_PART5



        # write rejected calls to file
        with open(os.path.join(params.output, params.sample, "_failed.vcf"), mode='w',encoding='utf-8') as OUT_TRASH_FILE:
            OUT_TRASH_FILE.write(VCF_HEADER)
            ordered_variant_key = functions.order_var(rejected.keys(), contigs)

            for variant_key in ordered_variant_key:
                OUT_TRASH_FILE.write(functions.print_var(variant_key, rejected, 'final_metrics'))
        OUT_TRASH_FILE.close()

        # write valid calls to file
        with open(params.output_file, 'w',encoding='utf-8') as OUT_FILE:
            OUT_FILE.write(VCF_HEADER)
            for variant_key in functions.order_var(callsets.keys(), contigs):
                OUT_FILE.write(functions.print_var(variant_key, callsets, 'final_metrics'))
        OUT_FILE.close()

        # print variant in Cartagenia compliant format
        # (temporary quicky/dirty addon for test phase purpose)
        cartagenia_file = params.output_file.split('.vcf')[0] + '_cartagenia.vcf'
        with open(cartagenia_file, 'w',encoding='utf-8') as OUT_CARTAGENIA:
            OUT_CARTAGENIA.write(VCF_HEADER_cartagenia)
            for variant_key in functions.order_var(callsets.keys(), contigs):

                OUT_CARTAGENIA.write(
                    functions.print_var4cartagenia(variant_key, callsets, 'final_metrics')
                )
        OUT_CARTAGENIA.close()
