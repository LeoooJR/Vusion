from files import Pileup
from hashlib import sha256
from loguru import logger
import re
import utils as functions

class VariantsRepository():

    def __init__(self):

        pass

    def normalize(self, variants: dict[str:dict], sample: str, pileup: Pileup, thresholds: list[float], length_indels: int, sbm: float) -> tuple[dict]:

        ITD = {}
        rejected = {}

        with open(pileup.get_path(), mode='r') as f:

            for n, line in enumerate(f, start=1):

                datas = line.strip('\n').split('\t')

                if datas[0] == sample:

                    # Create a key with format chromosome:POSITION

                    position_index: str = sha256(
                                            string=f"{(datas[pileup.HEADER['chromosome']]).removeprefix('chr')}:{datas[pileup.HEADER['position']]}".encode()
                                        ).hexdigest()

                    # Check if variant is reported in one of the VCF files
                    if (datas[pileup.HEADER['reference']] != 'N') and (position_index in readpile):
                        
                        for key in readpile[position_index]:

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

                            if not 'final_metrics' in variants[key]:
                                variants[key]['final_metrics'] = {}

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

                            variants[key]['final_metrics']['TRC'] = coverage["total"]
                            variants[key]['final_metrics']['TRC-'] = coverage['minus']
                            variants[key]['final_metrics']['TRC+'] = coverage['plus']

                            # Add <REF> strand-specific counts for each <ALT> at <CHR:POS>
                            # ie. the matched key in callset
                            # RRC : Reference Read Counts
                            # ARC : Alternative Read Counts
                            for strand in ['-','+']:

                                variants[key]['final_metrics'][f"RRC{strand}"] = datas[pileup.HEADER[f"{datas[pileup.HEADER['reference']]}{strand}"]]

                            # if <ALT> is an <INDEL>
                            variants_count: int = 0

                            if (variants[key]['VT'] == 'DEL') or (variants[key]['VT'] == 'INS'):
                                    
                                if re.search(r"\b" + readpile[position_index][key] + r"\b:[0-9]+,[0-9]+" ,\
                                    (datas[pileup.HEADER[variants[key]['VT']]])):
                                
                                        arc_plus_strand, arc_minus_strand = re.findall(
                                            r"\b" + readpile[position_index][key] + r"\b:[0-9]+,[0-9]+",
                                            datas[pileup.HEADER[variants[key]['VT']]]
                                            )[0].split(':')[-1].split(',')
                                        
                                        variants_count: int = int(arc_plus_strand) + int(arc_minus_strand)
                                        
                                        #Save indel count per strand in callset ARC+ and ARC-
                                        for strand in ['+','-']:

                                            arc = f'ARC{strand}'

                                            if not arc in variants[key]['final_metrics']:
                                                variants[key]['final_metrics'][arc] = ''

                                            variants[key]['final_metrics'][arc] = arc_plus_strand if strand == '+' else arc_minus_strand

                            # any other <VarType>
                            else:
                                # keep 1st character of <ALT> string (approx. counts for non-SNV variants)

                                if readpile[position_index][key][0] == 'N':
                                    del variants[key]
                                else:
                                    for strand in ['+','-']:

                                        arc = f'ARC{strand}'

                                        if not arc in variants[key]['final_metrics']:
                                            variants[key]['final_metrics'][arc] = ''

                                        variants[key]['final_metrics'][arc] = \
                                            datas[pileup.HEADER[f"{readpile[position_index][key][0]}{strand}"]]
                                        
                                        variants_count += int(
                                            datas[pileup.HEADER[f"{readpile[position_index][key][0]}{strand}"]]
                                        )
                            
                            # --------------------------------------------------------------
                            # <ALT> not covered in read pile;
                            # keep trace of rejected variants (for test / rescuing purpose)
                            # --------------------------------------------------------------
                            if variants_count == 0:

                                if (variants[key]['VT'] == 'INS') and (len(variants[key]["VC"]["ALT"])-1 > length_indels):

                                    ITD[key] = variants.pop(key)

                                else:

                                    rejected[key] = variants.pop(key)

                                continue

                            # --------------------------------------------------------------
                            # Compute metrics
                            # --------------------------------------------------------------

                            # Save number of Variant Callers that found this variant
                            variants[key]['final_metrics']['VCN'] = len(variants[key]['VC']['VAF'])

                            # keep trace of used VC identifier(s)
                            variants[key]['final_metrics']['VCI'] = ','.join(
                                sorted(variants[key]['VC']['VAF'].keys(), key=str.lower)
                            )

                            #format [R/A]RC for vcf
                            variants[key]['final_metrics']['ARC'],\
                            variants[key]['final_metrics']['RRC'] = \
                                functions.format_rrc_arc(variants,key)

                            # Compute VAF from pileup
                            alt_read_count_plus = variants[key]['final_metrics']['ARC+']
                            alt_read_count_minus = variants[key]['final_metrics']['ARC-']
                            total_read_count = variants[key]['final_metrics']['TRC']
                            variants[key]['final_metrics']['ARR'] = format(
                                (
                                    (float(alt_read_count_plus) + float(alt_read_count_minus)) /
                                    float(total_read_count)
                                )*100,
                                '.5f'
                            )

                            # Estimate GT
                            variants[key]['final_metrics']['VAR'] = \
                                functions.categorize_variant_type(variants,key,thresholds)
                            variants[key]['final_metrics']['GT'] = \
                                functions.estimate_gt(variants,key)
                            
                            # If variant callers do not agree on genotype change VAR to WAR (warning)
                            variants[key]['final_metrics']['VAR'] = \
                                functions.compare_gt(variants,key)

                            variants[key]['final_metrics']['BRC'],\
                            variants[key]['final_metrics']['BRR'],\
                            variants[key]['final_metrics']['BRE'] = \
                                functions.estimate_brc_r_e(variants,key,datas)

                            variants[key]['final_metrics']['SBP'],variants[key]['SBM'] = \
                                functions.estimate_sbm(variants,key,params.sbm_homozygous)
                            variants[key]['final_metrics']['LOW'] = \
                                functions.vaf_user_threshold(variants,key,thresholds)
                            variants[key]['final_metrics']['BKG'] = \
                                functions.categorize_background_signal(variants,key,thresholds)
                            variants[key]['final_metrics']['FILTER'] = \
                                functions.format_float_descriptors(variants, key, sbm)

                            variants[key]['final_metrics']['PIL'] = 'Y'
                            variants[key]['final_metrics']['RES'] = 'N'
        
        return (variants, ITD, rejected)