# -*- coding: utf-8 -*-
################################################################################
##  *
##  *  Function: prepare dag file to run the CNV analyses pipeline]
##  *  Writer:Siyuan Feng
##  *  Mail: siyuanfeng.bioinfo@gmail.com
##  *  Version: 11.04.2022
##  *
################################################################################

import os
from math import ceil

# 0. required directories and files
## 0.1. directories
dir_tar = '/staging/sfeng77/test_pipeline/out'
dir_out = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/results'
dir_scripts = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/scripts'
dir_software_pack = '/home/sfeng77/biosoft/biosoft'
dir_conda = '/home/sfeng77/anaconda3/bin/activate'
dir_in_stage = '/staging/sfeng77/CNV_input'

## 0.2. arguments
env_conda = 'WGS_analysis'

## 0.3. files
file_sample = '/home/sfeng77/jobs/rename_samples/rename_dict_rm_lowqual.txt' # table of sample with old and new names
file_contig = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/scripts/assignment_cor_amb.txt' # list of contigs
file_contig_len = '/home/sfeng77/reference/WT3-2.0/refseq/annotation/GCF_013340165.1_LBDM_Dsuz_2.1.pri_assembly_report.txt' # table including contig length
file_para = '/home/sfeng77/jobs/SV_analyses/CNV/insert_size/plot/parameters_CNV.txt'
file_conda_pack_squid = 'http://proxy.chtc.wisc.edu/SQUID/sfeng77/suzukii_WGS/WGS_analysis_latest.tar.gz'
file_dag = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/results/pipeline.dag'

# 1. build a dictionary for renamed samples
dict_rename = {}
with open(file_sample, 'r') as f:
  for line in f:
    line = line.strip().split('\t')
    old_name = line[0]
    new_name = line[1]
    dict_rename[old_name] = new_name
# dict_rename['merged'] = 'merged'

# 2. build a dictionary for sample-library specific parameters
dict_para = {}
with open(file_para, 'r') as f:
    i = 0
    for line in f:
        line = line.strip().split('\t')
        if i != 0:
            sample = line[0]
            insertSizeCutoff = line[1]
            insertSizeDiffCutoff = line[2]
            dict_para[sample] = [insertSizeCutoff, insertSizeDiffCutoff]
        i += 1

# # 3. build a list of contigs
# list_contig = []
# with open(file_contig, 'r') as f:
#     i = 0
#     for line in f:
#         line = line.strip().split('\t')
#         if i != 0:
#             contig = line[0]
#             list_contig.append(contig)
#         i += 1

# 4. build a dictionary of contig size
dict_contig_len = {}
with open(file_contig_len, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        if not "#" in line[0]:
            contig = line[6]
            length = int(line[8])
            dict_contig_len[contig] = length
list_contig_sorted = sorted(dict_contig_len, key = lambda x:dict_contig_len[x], reverse=True) # sort contigs by decreasing size

# 5. Write the DAG input file (.dag)
with open(file_dag, 'w') as f: 
    # set configuration variables for this dag, to remove default limits on job submissions
    # note!! the assumption to allow unlimited job submissions, is there should not be too many job submissions at the same time.
    f.write('CONFIG {0}/unlimited.config\n'.format(dir_scripts))

    for old_name in dict_rename:
        sample = dict_rename[old_name]
        path_out_sp = '{0}/{1}'.format(dir_out, sample)
        os.system('mkdir -p {0}'.format(path_out_sp))

    for contig in list_contig_sorted:
        for old_name in dict_rename:
            sample = dict_rename[old_name]
            path_out_sp = '{0}/{1}'.format(dir_out, sample)

            basic_disk = 2*1024**2
            basic_memory = 1024
            # ## 5.1. step 1: find and cluster discordant read pairs for each sample
            # ### determine computational resource by scaling contig size
            # request_cpus = 2
            # scale = float(dict_contig_len[contig]/dict_contig_len['NC_050699.1'])
            # request_disk = basic_disk + int(1024**2*scale)
            # request_memory = basic_memory + int(1000*scale)

            # ### estimate the largest possible bam file size by scaling contig size
            # est_bam_size = 35*scale
            # if sample == 'merged':
            #     request_disk = basic_disk + int(1024**2*scale*29)
            #     request_memory = basic_memory + int(1000*scale*29)
            #     est_bam_size = est_bam_size*29

            # ### file transfering method depends on estimated size of input bam file
            # if est_bam_size <= 500: # cutoff of HTC on whether one should input file from staging
            #     job = 'JOB 1_find_cluster_discordant_reads_{0}_{1} {2}/1_find_cluster_discordant_reads.sub\n'.format(sample, contig, dir_scripts)
            #     vars = '''VARS 1_find_cluster_discordant_reads_{0}_{1} \
            #         sample="{0}" contig="{1}" conda_pack_squid="{2}" insert_size_cutoff="{3}" insertSizeDiffCutoff="{4}" job="$(JOB)"\
            #         dir_out="{5}" dir_scripts="{6}" in_bam_compressed="{5}/{0}/{0}_{1}.bam.genozip" software_pack="{7}"\
            #         request_disk="{8}" request_memory="{9}" request_cpus="{10}"\n'''.format(
            #         sample, contig, file_conda_pack_squid, dict_para[sample][0], dict_para[sample][1], 
            #         dir_out, dir_scripts, dir_software_pack, 
            #         request_disk, request_memory, request_cpus)
            #     if sample == 'merged':
            #         pre = 'SCRIPT PRE 1_find_cluster_discordant_reads_{0}_{1} {7}/1_find_cluster_discordant_reads_merge_pre.sh {2} {3} {4} {5} {6} {6}/{0} {1} {8}\n'.format(
            #             sample, contig, dir_conda, env_conda, dir_tar, file_sample, dir_out, dir_scripts, dir_in_stage) 
            #     else:
            #         pre = 'SCRIPT PRE 1_find_cluster_discordant_reads_{0}_{1} {7}/1_find_cluster_discordant_reads_pre.sh {2} {3} {4} {5} {0} {6}/{0} {1} {8}\n'.format(
            #             sample, contig, dir_conda, env_conda, dir_tar, old_name, dir_out, dir_scripts, dir_in_stage)

            # elif est_bam_size > 500:
            #     job = 'JOB 1_find_cluster_discordant_reads_{0}_{1} {2}/1_find_cluster_discordant_reads_large.sub\n'.format(sample, contig, dir_scripts)
            #     vars = '''VARS 1_find_cluster_discordant_reads_{0}_{1} \
            #         sample="{0}" contig="{1}" conda_pack_squid="{2}" insert_size_cutoff="{3}" insertSizeDiffCutoff="{4}" job="$(JOB)"\
            #         dir_out="{5}" dir_scripts="{6}" in_bam_compressed="{5}/{0}/{0}_{1}.bam.genozip" software_pack="{7}"\
            #         request_disk="{8}" request_memory="{9}" request_cpus="{10}" path_in_stage="{11}"\n'''.format(
            #         sample, contig, file_conda_pack_squid, dict_para[sample][0], dict_para[sample][1], 
            #         dir_out, dir_scripts, dir_software_pack, 
            #         request_disk, request_memory, request_cpus,
            #         dir_in_stage) # variables added later
            #     if sample == 'merged':
            #         pre = 'SCRIPT PRE 1_find_cluster_discordant_reads_{0}_{1} {7}/1_find_cluster_discordant_reads_merge_pre_large.sh {2} {3} {4} {5} {6} {6}/{0} {1} {8}\n'.format(
            #             sample, contig, dir_conda, env_conda, dir_tar, file_sample, dir_out, dir_scripts, dir_in_stage) 
            #     else:
            #         pre = 'SCRIPT PRE 1_find_cluster_discordant_reads_{0}_{1} {7}/1_find_cluster_discordant_reads_pre_large.sh {2} {3} {4} {5} {0} {6}/{0} {1} {8}\n'.format(
            #             sample, contig, dir_conda, env_conda, dir_tar, old_name, dir_out, dir_scripts, dir_in_stage)
            # f.writelines([job, vars, pre])

            # ### write a dag to: generate bam for each contig (for next step of merging across samples) and then stored at /staging
            # ### doing this because it is not appropriate to generate bam f            request_disk = basic_disk + int(1024**2*scale)
            # request_disk = basic_disk + int(1024**2*10)
            # request_memory = basic_memory
            # request_cpus = 1

            # job = 'JOB 1.1_generate_bam_contig_{0}_{1} {2}/1.1_generate_bam_contig.sub\n'.format(sample, contig, dir_scripts)
            # vars = '''VARS 1.1_generate_bam_contig_{0}_{1} \
            #     sample="{0}" contig="{1}" conda_pack_squid="{2}" path_in_stage="{3}" job="$(JOB)"\
            #     dir_out="{4}" dir_scripts="{5}" in_bam_rename_sort="{4}/{0}/{0}_sort.bam" in_bam_rename_sort_index="{4}/{0}/{0}_sort.bam.bai"\
            #     request_disk="{6}" request_memory="{7}" request_cpus="{8}"\n'''.format(
            #     sample, contig, file_conda_pack_squid, dir_in_stage, 
            #     dir_out, dir_scripts, 
            #     request_disk, request_memory, request_cpus)
            # f.writelines([job, vars])

        ### write a dag to: 1. merge bam files for each contig across samples; 2. split into smaller portions; 3. find and cluster discordant read pairs for merged library
        scale = float(dict_contig_len[contig]/dict_contig_len['NC_050699.1'])
        request_disk = basic_disk + int(0.1*1024**2*scale*len(dict_rename)*3)
        request_memory = basic_memory + int(100*scale*len(dict_rename))
        request_cpus = 2
        est_bam_size = 90*scale*len(dict_rename)

        ### decide whether to split merged contigs by estimated bam size
        if est_bam_size >= 500:
            large = 'true'
        else:
            large = 'false'

        ### determine the number and positions of split by estimated bam size and local minima of insert coverage
        num_split = ceil(float(est_bam_size)/500) # use ceil to round up
        len_split = int(dict_contig_len[contig]/num_split)
        split_sites = [0] + list(int(x*dict_contig_len[contig]/num_split) for x in range(1, num_split)) + [dict_contig_len[contig]]
        split_sites_range = list([split_site - 0.1*len_split, split_site + 0.1*len_split] for split_site in split_sites[1:-1])
        insert_cov_distant = '{0}/merged/{1}_distant_cov.bed'.format(dir_out, contig)
        insert_cov_everted = '{0}/merged/{1}_everted_cov.bed'.format(dir_out, contig)

        #### define a function to find actual split sites as local minima of insert coverage
        def split_by_cov(insert_cov, split_sites_range, split_sites):
            split_sites_minima = []
            split_sites_range_cov = {}
            with open(insert_cov, 'r') as cd:
                split_site = 0
                split_sites_range_cov.setdefault(split_site, {})
                for line in cd:
                    line  = line.strip().split('\t')
                    start = int(line[1])
                    end = int(line[2])
                    cov = int(line[3])
                    range_start = split_sites_range[split_site][0]
                    range_end = split_sites_range[split_site][1]
                    if start <= range_end:
                        if end >= range_start:
                            split_sites_range_cov[split_site]['{0}_{1}'.format(max(range_start, start), min(range_end, end))] = cov
                    elif start > range_end:
                        if split_site < len(split_sites_range) - 1:
                            split_site += 1
                            split_sites_range_cov.setdefault(split_site, {})
                        elif split_site >= len(split_sites_range) - 1:
                            break

            for split_site in split_sites_range_cov:
                ranges_cov = sorted(split_sites_range_cov[split_site], key=lambda x: split_sites_range_cov[split_site][x]) # sort split site ranges by insert coverage
                cov_minima = split_sites_range_cov[split_site][ranges_cov[0]]
                ranges_cov_minima = list(x for x in ranges_cov if split_sites_range_cov[split_site][x] == cov_minima)
                range_minima_closest = sorted(ranges_cov_minima, key=lambda x: min(abs(int(x.split('_')[0]) - split_sites[split_site + 1]), abs(int(x.split('_')[1]) - split_sites[split_site + 1])))[0]
                split_sites_minima.append(range_minima_closest)
            split_sites_minima = [0] + split_sites_minima + [split_sites[-1]]
            return(split_sites_minima)

        #### find actual split sites for distant inserts
        split_sites_distant = split_by_cov(insert_cov_distant, split_sites_range, split_sites)
        #### find actual split sites for everted inserts
        split_sites_everted = split_by_cov(insert_cov_everted, split_sites_range, split_sites)

        split = 0
        while split < len(split_sites) - 1:
            start_distant = split_sites_distant[split]
            end_distant = split_sites_distant[split + 1]
            start_everted = split_sites_everted[split]
            end_everted = split_sites_everted[split + 1]
            bed_split = '{0}/merged/{1}_{2}.bed'.format(dir_out, contig, split)
            with open(bed_split, 'w') as fb:
                fb.write('\t'.join([contig, str(start_distant), str(end_distant)]))
                fb.write('\t'.join([contig, str(start_everted), str(end_everted)]))

            ### write dag
            job_relation = 'PARENT {0} CHILD 1_find_cluster_discordant_reads_merged_{1}_{2}\n'.format('\t'.join(list('1.1_generate_bam_contig_{0}_{1}'.format(sample, contig) for sample in dict_rename.values())), contig, split)
            job = 'JOB 1_find_cluster_discordant_reads_merged_{0}_{1} {2}/1_find_cluster_discordant_reads_merged.sub\n'.format(contig, split, dir_scripts)
            vars = '''VARS 1_find_cluster_discordant_reads_merged_{0}_{1} \
                contig="{0}" conda_pack_squid="{2}" insert_size_cutoff="{3}" insertSizeDiffCutoff="{4}" path_in_stage="{5}" split="{1}" large="{6}" job="$(JOB)"\
                dir_out="{7}" dir_scripts="{8}" bed_split="{9}" software_pack="{10}"\
                request_disk="{11}" request_memory="{12}" request_cpus="{13}"\n'''.format(
                contig, split, file_conda_pack_squid, dict_para[sample][0], dict_para[sample][1], dir_in_stage, large, 
                dir_out, dir_scripts, bed_split, dir_software_pack,
                request_disk, request_memory, request_cpus)
            f.writelines([job_relation, job, vars])
                
            split += 1






