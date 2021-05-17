#!/usr/bin/env python

'''
This script is to filter the xAtlas calls from a project-level vcf file
The input is an vcf file after joint calling using GLnexus.
Some variant will be called in a sample although the gvcf for that position in that sample is not PASS. 
This script is to convert the genotypes of those positions (indicated by the INFO/FT) into ./.
I'll also check the genotype to remove half-missing calls, e.g. ./1
Also, it allows filtering based on DP on all variants and allelic fraction on het variants
'''

import sys
import gzip
import argparse

parser = argparse.ArgumentParser(description='filter xAtlas joint genotyping calls after GLnexus')

parser.add_argument('-i',metavar='FILENAME', required=True,
                    help='Input file in .vcf or .vcf.gz format')
parser.add_argument('-d',required=False, metavar='MIN_DP', default = 0,
                    help='DP threshold to filter variants')
parser.add_argument('-r',required=False, metavar='MIN_ALLELIC_RATIO', default = 0,
                    help='AF threshold to filter variants')
parser.add_argument('-c',required=False, action='store_true',
                    help='do not generate vcf file, only output the count of variants being assigned as ./. for each variant')
args = parser.parse_args()

inputF = args.i
DP_thres = int(args.d)
ratio_thres = float(args.r)
count_only = args.c

def parseVcf(file):
    for line in file:
        line = line.decode("utf-8").strip()
        if line[0] == '#':
            if not count_only:
                print(line)
        else:
            line = line.split('\t')
            vid = line[2]
            # get variant info
            formats = line[8].split(':')
            FT_idx = formats.index('FT')
            DP_idx = formats.index('DP')
            VR_idx = formats.index('VR')
            samples = line[9:]
            
            # add count for samples that filled with './.'
            filtered_count = 0
            for i in range(len(samples)):
                sample = samples[i].split(':')
                gt = sample[0]
                ft = sample[FT_idx]
                dp = sample[DP_idx]
                vr = sample[VR_idx]
                
                # assign missing calls to failed variants
                if ft != '.':
                    sample[0] = './.'
                    filtered_count += 1
                elif '.' in gt:
                    sample[0] = './.'
                    filtered_count += 1
                    sample[FT_idx] = 'half_missing_genotype'
                elif gt[0] != '0' and gt[-1] != '0' and gt[0] != gt[-1]:
                    sample[0] = './.'
                    filtered_count += 1
                    sample[FT_idx] = 'alt_alt_genotype'
                elif DP_thres > 0 and int(dp) < DP_thres:
                    sample[0] = './.'
                    filtered_count += 1
                    sample[FT_idx] = 'low_DP_than_' + str(DP_thres)
                elif ratio_thres > 0 and gt[0] == '0' and gt[-1] != '0':
                    ratio = int(vr) / int(dp)
                    if ratio < ratio_thres or ratio > (1-ratio_thres):
                        sample[0] = './.'
                        filtered_count += 1
                        sample[FT_idx] = 'low_or_high_alt_fraction_than_' + str(ratio_thres)

                samples[i] = ':'.join(sample)
            if count_only:
                print('\t'.join([vid, str(filtered_count)]))
            else:
                print('\t'.join(line[:9] + samples))
            

def callParseVcf(inputF):                
    if inputF[-4:] == '.vcf':
        with open(inputF) as f:
            parseVcf(f)
    elif inputF[-7:] == '.vcf.gz':  
        with gzip.open(inputF, 'rb') as f:
            parseVcf(f)
    else:
        raise ValueError('vcf file name is not correct')

if __name__ == "__main__":
    callParseVcf(inputF)


