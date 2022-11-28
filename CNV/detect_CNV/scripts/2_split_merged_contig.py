# -*- coding: utf-8 -*-
################################################################################
##  *
##  *  Function: generate distributions of coverage of CNV-suggestive read pairs of eaach contig]
##  *  Writer:Siyuan Feng
##  *  Mail: siyuanfeng.bioinfo@gmail.com
##  *  Version: 11.11.2022
##  *
################################################################################

import os
from math import ceil

# 0. required directories & files & arguments & functions
## 0.1. directories
dir_out = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/results'
dir_out_merged = '{0}/merged'.format(dir_out)
dir_scripts = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/scripts'
dir_conda = '/home/sfeng77/anaconda3/bin/activate'
dir_in_stage = '/staging/sfeng77/CNV_input'
os.system('mkdir -p {}'.format(dir_out_merged))

## 0.2. arguments
env_conda = 'WGS_analysis'

## 0.3. files
file_sample = '/home/sfeng77/jobs/rename_samples/rename_dict_rm_lowqual.txt' # table of sample with old and new names
file_contig_len = '/home/sfeng77/reference/WT3-2.0/refseq/annotation/GCF_013340165.1_LBDM_Dsuz_2.1.pri_assembly_report.txt' # table including contig length
file_bed_genome = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/scripts/genome.bed' # length of contigs as a input for bedtools complement

# ## 0.4. functions
# ### function to get split sites of a merged contig
# def get_split_sites(bed_gap, contig, num_split, dict_contig_len):
#     split_sites_ideal = list(int(x*dict_contig_len[contig]/num_split) for x in range(1, num_split))
#     gap_sites_all = []
#     with open(bed_gap, 'r') as f:
#         for line in f:
#             line = line.strip().split('\t')
#             gap_sites = list(range(int(line[1]), int(line[2])))
#             gap_sites_all += gap_sites
#     split_sites_actual = list(sorted(gap_sites_all, key = lambda x: abs(x - site))[0] for site in split_sites_ideal) # choose sites that are closest to ideal split sites
#     return(split_sites_actual)

# 1. build a dictionary for renamed samples
dict_rename = {}
with open(file_sample, 'r') as f:
  for line in f:
    line = line.strip().split('\t')
    old_name = line[0]
    new_name = line[1]
    dict_rename[old_name] = new_name
# dict_rename['merged'] = 'merged'

# 2. build a dictionary of contig size
dict_contig_len = {}
with open(file_contig_len, 'r') as f, open(file_bed_genome, 'w') as fw:
    for line in f:
        line = line.strip().split('\t')
        if not "#" in line[0]:
            contig = line[6]
            length = int(line[8])
            dict_contig_len[contig] = length
            fw.write('\t'.join([contig, str(length)]) + '\n')

            file_bed_genome_contig = '{0}/genome_{1}.bed'.format(dir_out_merged, contig) # length of each contig as a input for bedtools genomecov
            with open(file_bed_genome_contig, 'w') as fc:
                fc.write('\t'.join([contig, str(length)]) + '\n')
list_contig_sorted = sorted(dict_contig_len, key = lambda x:dict_contig_len[x], reverse=True) # sort contigs by decreasing size

# 3. split merged input bam files by gaps along each contig
for contig in list_contig_sorted:

    ## 3.1. estimate the largest possible bam file size by scaling contig size
    scale = float(dict_contig_len[contig]/dict_contig_len['NC_050699.1'])
    est_bam_size = 35*scale*29

    ## 3.2. decide how/whether to split merged contigs by estimated bam size
    if est_bam_size >= 500:
        ### determine the number of split by estimated bam size
        num_split = ceil(float(est_bam_size)/500) # use ceil to round up
        split_sites = list(int(x*dict_contig_len[contig]/num_split) for x in range(1, num_split))
        
        # 3.3. generate bed files of distant and everted read pairs of merged sample
        bed_contig_distant = '{0}/{1}_distant.bed'.format(dir_out_merged, contig) # coordinates of inserts - distant
        bed_contig_everted = '{0}/{1}_everted.bed'.format(dir_out_merged, contig) # coordinates of inserts - everted
        bed_contig_distant_reads = '{0}/{1}_distant_reads.bed'.format(dir_out_merged, contig) # coordinates of reads - distant
        bed_contig_everted_reads = '{0}/{1}_everted_reads.bed'.format(dir_out_merged, contig) # coordinates of reads - everted
        bed_contig_distant_cluster = '{0}/{1}_distant_clustered.bed'.format(dir_out_merged, contig) # coordinates of clusters - distant
        bed_contig_everted_cluster = '{0}/{1}_everted_clustered.bed'.format(dir_out_merged, contig) # coordinates of clusters - everted

        # bed_gap_everted = '{0}/{1}_everted_gap.bed'.format(dir_out_merged, contig)
        bed_cov_distant = '{0}/{1}_distant_cov.bed'.format(dir_out_merged, contig)
        bed_cov_everted = '{0}/{1}_everted_cov.bed'.format(dir_out_merged, contig)
        bed_cov_distant_reads = '{0}/{1}_distant_cov_reads.bed'.format(dir_out_merged, contig)
        bed_cov_everted_reads = '{0}/{1}_everted_cov_reads.bed'.format(dir_out_merged, contig)
        bed_cov_distant_cluster = '{0}/{1}_distant_cov_cluster.bed'.format(dir_out_merged, contig)
        bed_cov_everted_cluster = '{0}/{1}_everted_cov_cluster.bed'.format(dir_out_merged, contig)

        with open(bed_contig_distant, 'w') as fd, open(bed_contig_everted, 'w') as fe, \
        open(bed_contig_distant_reads, 'w') as fdr, open(bed_contig_everted_reads, 'w') as fer, \
            open(bed_contig_distant_cluster, 'w') as fdc, open(bed_contig_everted_cluster, 'w') as fec:
            for old_name in dict_rename:
                sample = dict_rename[old_name]
                path_out_sp = '{0}/{1}'.format(dir_out, sample)
                in_distant = '{0}/{1}_{2}_distant_inserts.txt'.format(path_out_sp, sample, contig)
                in_everted = '{0}/{1}_{2}_everted_inserts.txt'.format(path_out_sp, sample, contig)

                with open(in_distant, 'r') as ifd:
                    for line in ifd:
                        line = line.strip().split('\t')
                        contig = line[1]
                        bed_ls = str(int(line[2]) - 1)
                        bed_le = line[3]
                        bed_rs = str(int(line[5]) - 1)
                        bed_re = line[6]
                        fd.write('\t'.join([contig, bed_ls, bed_re]) + '\n')
                        fdr.write('\t'.join([contig, bed_ls, bed_le]) + '\n')
                        fdr.write('\t'.join([contig, bed_rs, bed_re]) + '\n')

                with open(in_everted, 'r') as ife:
                    for line in ife:
                        line = line.strip().split('\t')
                        contig = line[1]
                        bed_ls = str(int(line[2]) - 1)
                        bed_le = line[3]
                        bed_rs = str(int(line[8]) - 1)
                        bed_re = line[9]
                        fe.write('\t'.join([contig, bed_ls, bed_re]) + '\n')
                        fer.write('\t'.join([contig, bed_ls, bed_le]) + '\n')
                        fer.write('\t'.join([contig, bed_rs, bed_re]) + '\n')

                in_distant_cluster = '{0}/{1}_{2}_distant_inserts_clustered.txt'.format(path_out_sp, sample, contig)
                in_everted_cluster = '{0}/{1}_{2}_everted_inserts_clustered.txt'.format(path_out_sp, sample, contig)
                with open(in_distant_cluster, 'r') as ifd:
                    for line in ifd:
                        line = line.strip().split('\t')
                        contig = line[0].split(',')[0]
                        read_ct = int(line[1])
                        if read_ct >= 4:
                            bed_s = str(int(line[0].split(',')[1]) - 1)
                            bed_e = line[0].split(',')[2]
                            fdc.write('\t'.join([contig, bed_s, bed_e]) + '\n')

                with open(in_everted_cluster, 'r') as ife:
                    for line in ife:
                        line = line.strip().split('\t')
                        contig = line[0].split(',')[0]
                        read_ct = int(line[1])
                        if read_ct >= 4:
                            bed_s = str(int(line[0].split(',')[1]) - 1)
                            bed_e = line[0].split(',')[2]
                            fec.write('\t'.join([contig, bed_s, bed_e]) + '\n')
        
        # # 3.4. identify gaps between merged read pairs using bedtools complement
        # os.system('''source {0} {1}
        # bedtools sort -i {2} | bedtools complement -L -i - -g {3} > {4}
        # bedtools sort -i {5} | bedtools complement -L -i - -g {3} > {6}'''.format(dir_conda, env_conda, bed_contig_distant, file_bed_genome, bed_gap_distant, bed_contig_everted, bed_gap_everted))

        # 3.5. generate a readout of coverage of inserts along each contig using bedtools genomecov
        file_bed_genome_contig = '{0}/genome_{1}.bed'.format(dir_out_merged, contig) # length of each contig as a input for bedtools genomecov
        # os.system('''source {0} {1}
        # bedtools sort -i {2} | bedtools genomecov -bga -i - -g {3} > {4}
        # bedtools sort -i {5} | bedtools genomecov -bga -i - -g {3} > {6}'''.format(dir_conda, env_conda, bed_contig_distant, file_bed_genome_contig, bed_cov_distant, bed_contig_everted, bed_cov_everted))

        # # 3.6. generate a readout of coverage of reads along each contig using bedtools genomecov
        # os.system('''source {0} {1}
        # bedtools sort -i {2} | bedtools genomecov -bga -i - -g {3} > {4}
        # bedtools sort -i {5} | bedtools genomecov -bga -i - -g {3} > {6}'''.format(dir_conda, env_conda, bed_contig_distant_reads, file_bed_genome_contig, bed_cov_distant_reads, bed_contig_everted_reads, bed_cov_everted_reads))

        # 3.7. generate a readout of coverage of clusters along each contig using bedtools genomecov
        os.system('''source {0} {1}
        bedtools sort -i {2} | bedtools genomecov -bga -i - -g {3} > {4}
        bedtools sort -i {5} | bedtools genomecov -bga -i - -g {3} > {6}'''.format(dir_conda, env_conda, bed_contig_distant_cluster, file_bed_genome_contig, bed_cov_distant_cluster, bed_contig_everted_cluster, bed_cov_everted_cluster))

        # 4.5. decide split sites based on identified gaps
        # split_sites_distant = get_split_sites(bed_gap_distant, contig, num_split, dict_contig_len)
        # split_sites_everted = get_split_sites(bed_contig_everted, contig, num_split, dict_contig_len)
        # print('split sites distant:', split_sites_distant)
        # print('split sites everted:', split_sites_everted)
