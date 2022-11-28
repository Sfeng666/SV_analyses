# -*- coding: utf-8 -*-
################################################################################
##  *
##  *  Function: generate length distributions of CNV-suggestive inserts of eaach contig]
##  *  Writer:Siyuan Feng
##  *  Mail: siyuanfeng.bioinfo@gmail.com
##  *  Version: 11.18.2022
##  *
################################################################################

import os

# 0. required directories & files & arguments & functions
## 0.1. directories
dir_out = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/results'
dir_out_merged = '{0}/merged'.format(dir_out)
os.system('mkdir -p {}'.format(dir_out_merged))

## 0.2. files
file_sample = '/home/sfeng77/jobs/rename_samples/rename_dict_rm_lowqual.txt' # table of sample with old and new names
file_contig_len = '/home/sfeng77/reference/WT3-2.0/refseq/annotation/GCF_013340165.1_LBDM_Dsuz_2.1.pri_assembly_report.txt' # table including contig length
file_len_distr_distant = '{0}/len_distribution_distant_inserts.txt'.format(dir_out)
file_len_distr_everted = '{0}/len_distribution_everted_inserts.txt'.format(dir_out)
file_len_distr_distant_cluster = '{0}/len_distribution_distant_clusters.txt'.format(dir_out)
file_len_distr_everted_cluster = '{0}/len_distribution_everted_clusters.txt'.format(dir_out)
file_ct_distr_distant_cluster = '{0}/ct_distribution_distant_clusters_merged.txt'.format(dir_out)
file_ct_distr_everted_cluster = '{0}/ct_distribution_everted_clusters_merged.txt'.format(dir_out)

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
with open(file_contig_len, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        if not "#" in line[0]:
            contig = line[6]
            length = int(line[8])
            dict_contig_len[contig] = length
list_contig_sorted = sorted(dict_contig_len, key = lambda x:dict_contig_len[x], reverse=True) # sort contigs by decreasing size

# # 3. calculate whole-genome length distribution of CNV-suggestive inserts and clusters
# with open(file_len_distr_distant, 'w') as fd, open(file_len_distr_everted, 'w') as fe, \
#     open(file_len_distr_distant_cluster, 'w') as fdc, open(file_len_distr_everted_cluster, 'w') as fec:
#     for old_name in dict_rename:
#         sample = dict_rename[old_name]
#         path_out_sp = '{0}/{1}'.format(dir_out, sample)
#         for contig in list_contig_sorted:
#             ## 3.1. add whole-genome estimated duplication or deletion sizes to a file
#             ## output example: CN-Nin    1000
#             in_distant = '{0}/{1}_{2}_distant_inserts.txt'.format(path_out_sp, sample, contig)
#             in_everted = '{0}/{1}_{2}_everted_inserts.txt'.format(path_out_sp, sample, contig)

#             with open(in_distant, 'r') as ifd:
#                 for line in ifd:
#                     line = line.strip().split('\t')
#                     ls = int(line[2])
#                     re = int(line[6])
#                     del_size = re - ls + 1
#                     fd.write('\t'.join([sample, str(del_size)]) + '\n')

#             with open(in_everted, 'r') as ife:
#                 for line in ife:
#                     line = line.strip().split('\t')
#                     ls = int(line[2])
#                     le = int(line[3])
#                     rs = int(line[8])
#                     re = int(line[9])
#                     dup_size = max(ls, le, rs, re) - min(ls, le, rs, re) + 1
#                     fe.write('\t'.join([sample, str(dup_size)]) + '\n')

#             ## 3.2. add whole-genome estimated duplication or deletion cluster sizes (mutiplied by the number of supportive read pairs) to a file
#             ## output example: CN-Nin    1000   5
#             in_distant_cluster = '{0}/{1}_{2}_distant_inserts_clustered.txt'.format(path_out_sp, sample, contig)
#             in_everted_cluster = '{0}/{1}_{2}_everted_inserts_clustered.txt'.format(path_out_sp, sample, contig)
#             with open(in_distant_cluster, 'r') as ifd:
#                 for line in ifd:
#                     line = line.strip().split('\t')
#                     read_ct = line[1]
#                     cs = int(line[0].split(',')[1])
#                     ce = int(line[0].split(',')[2])
#                     csize = ce - cs + 1
#                     fdc.write('\t'.join([sample, str(csize), read_ct]) + '\n')

#             with open(in_everted_cluster, 'r') as ife:
#                 for line in ife:
#                     line = line.strip().split('\t')
#                     read_ct = line[1]
#                     cs = int(line[0].split(',')[1])
#                     ce = int(line[0].split(',')[2])
#                     csize = ce - cs + 1
#                     fec.write('\t'.join([sample, str(csize), read_ct]) + '\n')

# 4. extract the number of CNV-supporting read pairs for each cluster for distribution plotting
os.system("awk '{{print $2}}' {0}/merged_N*_*_distant_inserts_clustered.txt > {1}".format(dir_out_merged, file_ct_distr_distant_cluster))
os.system("awk '{{print $2}}' {0}/merged_N*_*_everted_inserts_clustered.txt > {1}".format(dir_out_merged, file_ct_distr_everted_cluster))

