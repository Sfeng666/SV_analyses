from optparse import OptionParser
import math
import sympy
from statistics import median

# 0.1. help and usage 
usage = "usage: %prog [options] args"
description = '''Function: prepare the split input genotype data of CNVs for GEA (using BayeScEnv)'''
version = '%prog 01.02.2023'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--wd",
                    action="store",
                    dest = "wd",
                    help = "working directory",
                    metavar = "PATH")
parser.add_option("--in_snp_dir",
                    action="store",
                    dest = "in_snp_dir",
                    help = "directory of filtered SNPs from all contigs in vcf format (.vcf)",
                    metavar = "PATH")
parser.add_option("--in_del_af",
                    action="store",
                    dest = "in_del_af",
                    help = "file path of allele frequency of filtered deletions (.txt)",
                    metavar = "PATH")
parser.add_option("--in_dup_af",
                    action="store",
                    dest = "in_dup_af",
                    help = "input path of allele frequency of filtered duplications (.txt)",
                    metavar = "PATH")
parser.add_option("--in_sample_size",
                    action="store",
                    dest = "in_sample_size",
                    help = "input path of allelic (auto/X) sample sizes of all samples analyzed (.txt)",
                    metavar = "PATH")
parser.add_option("--in_chr_assign",
                    action="store",
                    dest = "in_chr_assign",
                    help = "input path of assignments of contigs to chromosomes (.txt)",
                    metavar = "PATH")
parser.add_option("--in_sample_order",
                    action="store",
                    dest = "in_sample_order",
                    help = "input path of a file that contains the order of samples as environmental variables (.txt)",
                    metavar = "PATH")
parser.add_option("--out_efs_sp",
                    action="store",
                    dest = "out_efs_sp",
                    help = "output path of intermediate stats of median effective sample size (.txt)",
                    metavar = "PATH") 
parser.add_option("--out_log",
                    action="store",
                    dest = "out_log",
                    help = "output path of the log file (.log)",
                    metavar = "PATH")                                        
(options,args) = parser.parse_args()

# 0.2. introduced variables
wd = options.wd
in_snp_dir = options.in_snp_dir
in_del_af = options.in_del_af
in_dup_af = options.in_dup_af
in_sample_size = options.in_sample_size
in_chr_assign = options.in_chr_assign
in_sample_order = options.in_sample_order
out_efs_sp = options.out_efs_sp
out_log = options.out_log

# 0.3. local variables
subsample_size = 10000

# 0.4. def functions
## function to calculate effective sample size as the expectation of unique chromosomal draws
def expectation_sympy(nr, nc):
    return float(sum(j*math.factorial(nc)*sympy.functions.combinatorial.numbers.stirling(nr, j)/(math.factorial(nc - j)*nc**nr) for j in range(1, nc + 1)))

# 1. input chromosomal assignment information of contigs
chrs = ['X', 'auto']
chr_assign = {x:[] for x in chrs}
contig_assign = {}
with open(in_chr_assign, 'r') as f:
    for line in f:
        if not line.startswith('contig'):
            line = line.strip().split('\t')
            contig = line[0]
            cor_assign = line[2]
            assign = cor_assign.split('-')[0]
            if cor_assign != 'ambiguous':
                chr_assign[assign].append(contig)
                contig_assign[contig] = assign

# 2. extract haplotype size for each sample
sample_size = {}
with open(in_sample_size, 'r') as f:
    i = 0
    for line in f:
        line = line.strip().split('\t')
        if i == 0:
            samples = line
        elif i == 1:
            sample_size['auto'] = {x: int(y) for x,y in zip(samples, line)}
        elif i == 2:
            sample_size['X'] = {x: int(y) for x,y in zip(samples, line)}
        i += 1

# 3. extract the order of samples in the input of environmental variables
sample_order = []
with open(in_sample_order, 'r') as f:
    for line in f:
        if not line.startswith('Population'):
            line = line.strip().split('\t')
            sample_order.append(line[0])

# # 4. calculate effective sample size of each SNP for autosomes and X chromosome of each population
# expectation_dict = {}
# efs_sps = {chr:{sample: [] for sample in sample_order} for chr in chrs}
# efs_sps_median = {sample: {} for sample in sample_order}
# for chr in chrs:
#     for contig in chr_assign[chr]:
#         in_snp = '{0}/{1}_snp_clean.vcf'.format(in_snp_dir, contig)
#         with open(in_snp, 'r') as f:
#             for line in f:
#                 if line.startswith('#CHROM'):
#                     samples = line.strip().split(' ')[-1].split('\t')
#                 elif not line.startswith('#'):
#                     line = line.strip().split('\t')
#                     ref = line[3].upper()
#                     alt = line[4].upper()
#                     alleles = [ref] + alt.split(',')
#                     allele_freq = {x:0 for x in alleles}
#                     ct_sp = {}

#                     ## extract allele frequency of the top 2 alleles across all samples
#                     for i in range(len(line[9:])):
#                         info = line[9:][i]
#                         gt = list(int(x) for x in info.split(':')[0].split('/'))
#                         ct = [int(info.split(':')[1])] + list(int(x) for x in info.split(':')[2].split(','))
#                         dp = int(info.split(':')[3])

#                         if len(set(gt)) == 1:
#                             ct_sp[i] = {alleles[gt[0]]: max(ct)}
#                         else:
#                             ct_sp[i] = {alleles[x]: ct[gt.index(x)] for x in gt}

#                         for ale in alleles:
#                             if ale in ct_sp[i]:
#                                 allele_freq[ale] += ct_sp[i][ale]/dp
#                             else:
#                                 ct_sp[i][ale] = 0

#                     ## sort alleles by the total allele frequency (so that each sample has the same weight of voting despite sequencing depth), and determine the minor allele
#                     bialleles = sorted(allele_freq, key=lambda x: allele_freq[x], reverse=True)[:2]

#                     ## calculate effective sample size of an allele for each sample
#                     for i in range(len(line[9:])):
#                         dp_bi = sum(list(ct_sp[i][ale] for ale in bialleles if ale in ct_sp[i]))
#                         ### only calculate effective sample size for combinations of depth and pool size that haven't been calculated
#                         key_sz = ','.join([str(dp_bi), str(sample_size[chr][samples[i]])])
#                         if key_sz in expectation_dict:
#                             efs = expectation_dict[key_sz]
#                         else:
#                             efs = expectation_sympy(dp_bi, sample_size[chr][samples[i]])
#                             expectation_dict[key_sz] = efs
#                         efs_sps[chr][samples[i]].append(efs)

# # 5. calculate median effective sample size for autosomes and X chromosome of each population
# for chr in efs_sps:
#     for sample in efs_sps[chr]:
#         efs_sps_median[sample][chr] = median(efs_sps[chr][sample])
# with open(out_efs_sp, 'w') as f:
#     header = '\t'.join(['sample', 'auto', 'X']) + '\n'
#     f.write(header)
#     for sample in efs_sps_median:
#         row = '\t'.join([sample, str(efs_sps_median[sample]['auto']), str(efs_sps_median[sample]['X'])]) + '\n'
#         f.write(row)

# for test only
efs_sps_median = {sample: {} for sample in sample_order}
with open(out_efs_sp, 'r') as f:
    i = 0
    for line in f:
        line = line.strip().split('\t')
        if i > 0:
            sample = line[0]
            efs_auto = float(line[1])
            efs_X = float(line[2])
            efs_sps_median[sample]['auto'] = efs_auto
            efs_sps_median[sample]['X'] = efs_X
        i += 1

# 6. generate genotype file of CNVs for BayeScEnv
## define a function to generate genotype file of CNVs for BayeScEnv
def generate_gt_file(in_af, cnv):
    ### transfer allele frequency data into count data for CNVs
    dic_gt = {chr: {sample: {} for sample in sample_order} for chr in chrs}
    with open(in_af, 'r') as f:
        i = 0
        for line in f:
            line = line.strip().split('\t')
            if i == 0:
                samples = list(x.split('_')[1] for x in line[2:] if 'MAF' in x)
            elif i > 0:
                contig = line[0]
                if contig in contig_assign:
                    pos = line[1]
                    afs_cnv = list(float(line[2:][j*2 + 1]) for j in range(len(samples))) # allele frequency of CNVs
                    cts_cnv = list(int(afs_cnv[samples.index(sample)]*efs_sps_median[sample][contig_assign[contig]] + 0.5) for sample in samples) # estimated count of CNV alleles
                    cts_ref = list(int(efs_sps_median[sample][contig_assign[contig]] + 0.5) - cts_cnv[samples.index(sample)] for sample in samples) # estimated count of reference allales
                    for sample in samples:
                        dic_gt[contig_assign[contig]][sample][len(dic_gt[contig_assign[contig]][sample]) + 1] = [':'.join([contig, pos]), \
                            str(cts_cnv[samples.index(sample)] + cts_ref[samples.index(sample)]), \
                            '2', str(cts_cnv[samples.index(sample)]), \
                                str(cts_ref[samples.index(sample)])]  # extract allele count (of top two alleles) of each sample, corrected by median effective sample size
            i += 1

    for chr in chrs:
        out_gt_pos = '{0}/{1}_{2}_pos_cnv_gt.txt'.format(wd, cnv, chr)
        with open(out_gt_pos, 'w') as f_pos, open(out_log, 'a') as f_log:
            # determine the interval length of subsampling by sample size. interval length = numebr of samples
            interval = math.ceil(len(dic_gt[chr][sample_order[0]])/subsample_size)
            subsamplesize_div = math.ceil(len(dic_gt[chr][sample_order[0]])/interval)
            f_log.write('{0}: Number of subsamples for {1} = {2}\n'.format(chr, cnv, interval))

            # output subsampled CNVs
            for sub in range(1, interval + 1):
                out_gt = '{0}/{1}_{2}_{3}_cnv_gt.txt'.format(wd, cnv, chr, sub)
                with open(out_gt, 'w') as f:
                    if sub + (subsamplesize_div - 1)*interval <= len(dic_gt[chr][sample_order[0]]):
                        subsamplesize_actual = subsamplesize_div
                    else:
                        subsamplesize_actual = subsamplesize_div - 1
                    f.write('[loci]={}\n\n'.format(subsamplesize_actual))
                    f.write('[populations]={}\n\n'.format(len(sample_order)))
                    for pop in sample_order:
                        f.write('[pop]={}\n'.format(sample_order.index(pop) + 1))
                        for i in range(subsamplesize_actual):
                            f.write('\t'.join([str(i + 1)] + dic_gt[chr][pop][sub + i*interval][1:]) + '\n')
                            if pop == sample_order[0]: # avoid repeated output to the position file
                                f_pos.write(('\t'.join([str(sub), str(i + 1), dic_gt[chr][pop][sub + i*interval][0]]) + '\n'))
                        f.write('\n')

## 6.1. for duplications
generate_gt_file(in_dup_af, 'dup')
## 6.2. for deletions
generate_gt_file(in_del_af, 'del')

