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
from statistics import mean
from math import ceil

# 0. required directories and files
## 0.1. directories
dir_tar = '/staging/sfeng77/test_pipeline/out'
dir_out = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/results'
dir_out_merged = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/results/merged'
dir_out_global = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/results/global_CNV_AF'
dir_scripts = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/scripts'
dir_software_pack = '/home/sfeng77/biosoft/biosoft'
dir_conda = '/home/sfeng77/anaconda3/bin/activate'
dir_in_stage = '/staging/sfeng77/CNV_input'

## 0.2. arguments
env_conda = 'WGS_analysis'
maf_cutoff = 0.05 # cutoff for average MAF across all samples

## 0.3. files
file_sample = '/home/sfeng77/jobs/rename_samples/rename_dict_rm_lowqual.txt' # table of sample with old and new names
file_contig = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/scripts/assignment_cor_amb.txt' # list of contigs
file_contig_len = '/home/sfeng77/reference/WT3-2.0/refseq/annotation/GCF_013340165.1_LBDM_Dsuz_2.1.pri_assembly_report.txt' # table including contig length
file_para = '/home/sfeng77/jobs/SV_analyses/CNV/insert_size/plot/parameters_CNV.txt'
file_conda_pack_squid = 'http://proxy.chtc.wisc.edu/SQUID/sfeng77/suzukii_WGS/WGS_analysis_latest.tar.gz'
file_dag = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/results/pipeline.dag'
file_global_cluster_del = '{0}/global_clusters_del.txt'.format(dir_out_global) # table of global deletions across all samples
file_global_cluster_dup = '{0}/global_clusters_dup.txt'.format(dir_out_global) # table of global duplications across all samples
file_global_cluster_del_af = '{0}/global_clusters_del_af.txt'.format(dir_out_global) # table of global deletions with average allele frequency above MAF threshold across all samples
file_global_cluster_dup_af = '{0}/global_clusters_dup_af.txt'.format(dir_out_global) # table of global duplications with average allele frequency above MAF threshold across all samples
os.system('mkdir -p {0}'.format(dir_out_global))

# 1. build a dictionary for renamed samples
dict_rename = {}
with open(file_sample, 'r') as f:
  for line in f:
    line = line.strip().split('\t')
    old_name = line[0]
    new_name = line[1]
    dict_rename[old_name] = new_name
# dict_rename['merged'] = 'merged'
sample_order = list(dict_rename.values())

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

        # ### write a dag to: 1. merge bam files for each contig across samples; 2. split into smaller portions; 3. find and cluster discordant read pairs for merged library
        # scale = float(dict_contig_len[contig]/dict_contig_len['NC_050699.1'])
        # request_disk = basic_disk + int(0.4*1024**2*scale*len(dict_rename)*3)
        # request_memory = basic_memory + int(100*scale*len(dict_rename))
        # request_cpus = 2
        # est_bam_size = 35*scale*len(dict_rename) # the estimation has to be the same with that in '2_split_merged_contig.py'!

        # #### define a function to find actual split sites as local minima of insert coverage
        # def split_by_cov(insert_cov, split_sites_range, split_sites):
        #     split_sites_minima = []
        #     split_sites_range_cov = {}
        #     with open(insert_cov, 'r') as cd:
        #         split_site = 0
        #         split_sites_range_cov.setdefault(split_site, {})
        #         for line in cd:
        #             line  = line.strip().split('\t')
        #             start = int(line[1])
        #             end = int(line[2])
        #             cov = int(line[3])
        #             range_start = split_sites_range[split_site][0]
        #             range_end = split_sites_range[split_site][1]
        #             if start <= range_end:
        #                 if end >= range_start:
        #                     split_sites_range_cov[split_site]['{0}_{1}'.format(max(range_start, start), min(range_end, end))] = cov
        #             elif start > range_end:
        #                 if split_site < len(split_sites_range) - 1:
        #                     split_site += 1
        #                     split_sites_range_cov.setdefault(split_site, {})
        #                 elif split_site >= len(split_sites_range) - 1:
        #                     break

        #     for split_site in split_sites_range_cov:
        #         ranges_cov = sorted(split_sites_range_cov[split_site], key=lambda x: split_sites_range_cov[split_site][x]) # sort split site ranges by insert coverage
        #         cov_minima = split_sites_range_cov[split_site][ranges_cov[0]]
        #         ranges_cov_minima = list(x for x in ranges_cov if split_sites_range_cov[split_site][x] == cov_minima)
        #         range_minima_closest = sorted(ranges_cov_minima, key=lambda x: min(abs(int(x.split('_')[0]) - split_sites[split_site + 1]), abs(int(x.split('_')[1]) - split_sites[split_site + 1])))[0]
        #         site_minima_closest = sorted(list(range(int(range_minima_closest.split('_')[0]), int(range_minima_closest.split('_')[1]))), key=lambda x: abs(x - split_sites[split_site + 1]))[0]
        #         split_sites_minima.append(site_minima_closest)
        #     split_sites_minima = [0] + split_sites_minima + [split_sites[-1]]
        #     return(split_sites_minima)

        # ### decide whether to split merged contigs by estimated bam size
        # if est_bam_size >= 500:
        #     large = 'true'

        #     ### determine the number and positions of split by estimated bam size and local minima of insert coverage
        #     num_split = ceil(float(est_bam_size)/500) # use ceil to round up
        #     len_split = int(dict_contig_len[contig]/num_split)
        #     split_sites = [0] + list(int(x*dict_contig_len[contig]/num_split) for x in range(1, num_split)) + [dict_contig_len[contig]]
        #     split_sites_range = list([int(split_site - 0.1*len_split), int(split_site + 0.1*len_split)] for split_site in split_sites[1:-1])
        #     insert_cov_distant = '{0}/{1}_distant_cov.bed'.format(dir_out_merged, contig)
        #     insert_cov_everted = '{0}/{1}_everted_cov.bed'.format(dir_out, contig)

        #     #### find actual split sites for distant inserts
        #     split_sites_distant = split_by_cov(insert_cov_distant, split_sites_range, split_sites)
        #     #### find actual split sites for everted inserts
        #     split_sites_everted = split_by_cov(insert_cov_everted, split_sites_range, split_sites)

        #     split = 0
        #     while split < len(split_sites) - 1:
        #         start_distant = split_sites_distant[split]
        #         end_distant = split_sites_distant[split + 1]
        #         start_everted = split_sites_everted[split]
        #         end_everted = split_sites_everted[split + 1]
        #         bed_split = '{0}/{1}_{2}.bed'.format(dir_out_merged, contig, split)
        #         with open(bed_split, 'w') as fb:
        #             fb.write('\t'.join([contig, str(start_distant), str(end_distant)]) + '\n')
        #             fb.write('\t'.join([contig, str(start_everted), str(end_everted)]) + '\n')

        #         ### write dag
        #         job_relation = 'PARENT {0} CHILD 1_find_cluster_discordant_reads_merged_{1}_{2}\n'.format('\t'.join(list('1.1_generate_bam_contig_{0}_{1}'.format(sample, contig) for sample in dict_rename.values())), contig, split)
        #         job = 'JOB 1_find_cluster_discordant_reads_merged_{0}_{1} {2}/1_find_cluster_discordant_reads_merged.sub\n'.format(contig, split, dir_scripts)
        #         vars = '''VARS 1_find_cluster_discordant_reads_merged_{0}_{1} \
        #             contig="{0}" conda_pack_squid="{2}" insert_size_cutoff="{3}" insertSizeDiffCutoff="{4}" path_in_stage="{5}" split="{1}" large="{6}" job="$(JOB)"\
        #             dir_out="{7}" dir_scripts="{8}" bed_split="{9}" software_pack="{10}"\
        #             request_disk="{11}" request_memory="{12}" request_cpus="{13}"\n'''.format(
        #             contig, split, file_conda_pack_squid, dict_para[sample][0], dict_para[sample][1], dir_in_stage, large, 
        #             dir_out, dir_scripts, bed_split, dir_software_pack,
        #             request_disk, request_memory, request_cpus)
        #         # f.writelines([job_relation, job, vars]) # annotate this line if parental jobs has all completed
        #         f.writelines([job, vars]) # repalce above line with this one if parental jobs has all completed
        #         split += 1        
        # else:
        #     large = 'false'
        #     split = 0
        #     start = 0
        #     end = dict_contig_len[contig]
        #     bed_split = '{0}/merged/{1}_{2}.bed'.format(dir_out, contig, split)
        #     with open(bed_split, 'w') as fb:
        #         fb.write('\t'.join([contig, str(start), str(end)]))

        #     ### write dag
        #     job_relation = 'PARENT {0} CHILD 1_find_cluster_discordant_reads_merged_{1}_{2}\n'.format('\t'.join(list('1.1_generate_bam_contig_{0}_{1}'.format(sample, contig) for sample in dict_rename.values())), contig, split)
        #     job = 'JOB 1_find_cluster_discordant_reads_merged_{0}_{1} {2}/1_find_cluster_discordant_reads_merged.sub\n'.format(contig, split, dir_scripts)
        #     vars = '''VARS 1_find_cluster_discordant_reads_merged_{0}_{1} \
        #         contig="{0}" conda_pack_squid="{2}" insert_size_cutoff="{3}" insertSizeDiffCutoff="{4}" path_in_stage="{5}" split="{1}" large="{6}" job="$(JOB)"\
        #         dir_out="{7}" dir_scripts="{8}" bed_split="{9}" software_pack="{10}"\
        #         request_disk="{11}" request_memory="{12}" request_cpus="{13}"\n'''.format(
        #         contig, split, file_conda_pack_squid, dict_para[sample][0], dict_para[sample][1], dir_in_stage, large, 
        #         dir_out, dir_scripts, bed_split, dir_software_pack,
        #         request_disk, request_memory, request_cpus)
        #     # f.writelines([job_relation, job, vars]) # annotate this line if parental jobs has all completed
        #     f.writelines([job, vars]) # repalce above line with this one if parental jobs has all completed

        # # 5.2. step 2: get the global identity of candidate CNVs within each population sample by comparing them to reference CNVs
        # request_disk = 1024**2
        # request_memory = 1024*2*1.5
        # request_cpus = 2

        # ### combine merged clusters from split files of each contig
        # merged_everted_inserts_clustered_split = "{0}/merged_{1}_*_everted_inserts_clustered.txt".format(dir_out_merged, contig)
        # merged_distant_inserts_clustered_split = "{0}/merged_{1}_*_distant_inserts_clustered.txt".format(dir_out_merged, contig)
        # merged_everted_inserts_clustered = "{0}/merged_{1}_everted_inserts_clustered.txt".format(dir_out_merged, contig)
        # merged_distant_inserts_clustered = "{0}/merged_{1}_distant_inserts_clustered.txt".format(dir_out_merged, contig)
        # os.system('''cat {0} > {1}
        # cat {2} > {3}'''.format(merged_everted_inserts_clustered_split, merged_everted_inserts_clustered, merged_distant_inserts_clustered_split, merged_distant_inserts_clustered))

        # for old_name in dict_rename:
        #     sample = dict_rename[old_name]
        #     path_out_sp = '{0}/{1}'.format(dir_out, sample)
        #     sample_everted_inserts_clustered = "{0}/{1}/{1}_{2}_everted_inserts_clustered.txt".format(dir_out, sample, contig)
        #     sample_distant_inserts_clustered = "{0}/{1}/{1}_{2}_distant_inserts_clustered.txt".format(dir_out, sample, contig)

        #     ### write a dag to compare population CNVs to reference CNVs for each population
        #     job = 'JOB 2_compare_clusters_across_samples_{0}_{1} {2}/2_compare_clusters_across_samples.sub\n'.format(contig, sample, dir_scripts)
        #     vars = '''VARS 2_compare_clusters_across_samples_{0}_{1} \
        #         sample="{1}" contig="{0}" insertSizeDiffCutoff="{2}" job="$(JOB)"\
        #         dir_out="{3}" dir_scripts="{4}" merged_everted_inserts_clustered="{5}" merged_distant_inserts_clustered="{6}" sample_everted_inserts_clustered="{7}" sample_distant_inserts_clustered="{8}" software_pack="{9}"\
        #         request_disk="{10}" request_memory="{11}" request_cpus="{12}"\n'''.format(
        #         contig, sample, dict_para[sample][1], 
        #         dir_out, dir_scripts, merged_everted_inserts_clustered, merged_distant_inserts_clustered, sample_everted_inserts_clustered, sample_distant_inserts_clustered, dir_software_pack,
        #         request_disk, request_memory, request_cpus)
        #     f.writelines([job, vars])

    # ## merge pairwise comparisons between ref and samples to one ref-samples table
    # header = '\t'.join(['contig', 'coordinates_ref', 'count_ref']) \
    #     + '\t' + '\t'.join(list('coordinates_{0}\tcount_{0}'.format(sample) for sample in sample_order)) + '\n'
    # with open(file_global_cluster_del, 'w') as fgcd, open(file_global_cluster_dup, 'w') as fgce:
    #     fgcd.write(header)
    #     fgce.write(header)

    #     for contig in list_contig_sorted:
    #         ### build dict entries for ref CNVs identified from inserts data merged across samples
    #         #### for deletions
    #         dict_cnv_del = {}        
    #         merged_distant_inserts_clustered = "{0}/merged_{1}_distant_inserts_clustered.txt".format(dir_out_merged, contig)
    #         with open(merged_distant_inserts_clustered, 'r') as fmdel:
    #             for line in fmdel:
    #                 line = line.strip().split('\t')
    #                 cnv = line[0]
    #                 ct = line[1]
    #                 dict_cnv_del.setdefault(cnv, [ct, {}])
    #         #### for duplications
    #         dict_cnv_dup = {}
    #         merged_everted_inserts_clustered = "{0}/merged_{1}_everted_inserts_clustered.txt".format(dir_out_merged, contig)
    #         with open(merged_everted_inserts_clustered, 'r') as fmdup:
    #             for line in fmdup:
    #                 line = line.strip().split('\t')
    #                 cnv = line[0]
    #                 ct = line[1]
    #                 dict_cnv_dup.setdefault(cnv, [ct, {}])

    #         ### add each pairwise comparison to global dicts of deletions and duplications
    #         for sample in sample_order:
    #             #### for deletions
    #             pairwise_cluster_del = '{0}/{1}_distant_clusters_{2}.txt'.format(dir_out_merged, contig, sample)
    #             with open(pairwise_cluster_del, 'r') as fpdel:
    #                 for line in fpdel:
    #                     line = line.strip().split('\t')
    #                     # pw_coord = line[0] # this set coordinates could be: 1. coordinates of ref sample (may be same as 2); 2. coordinates of individual sample (may be same as 1); 3. new coordinates inferred from all read pairs from ref and individual samples;
    #                     ref_coord = line[0]
    #                     ref_ct = line[1]
    #                     sp_coord = line[2]
    #                     if sp_coord != 'NA':
    #                         sp_coord = ','.join(sp_coord.split(',')[1:])
    #                     sp_ct = str(int(float(line[3])))
    #                     if ref_coord in dict_cnv_del:
    #                             dict_cnv_del[ref_coord][1].setdefault(sample, [sp_coord, sp_ct])
    #                     # else:
    #                     #         if ref_coord == 'NA' and pw_coord in dict_cnv_del: # there are cases when ref coordinate is 'NA' due to matching criteria not met with individual CNV, but the individual coordinate is the same as that of the merged sample. We would like to recue such matcthes.
    #                     #             dict_cnv_del[pw_coord][1].setdefault(sample, [sp_coord, sp_ct])

    #             #### for duplications
    #             pairwise_cluster_dup = '{0}/{1}_everted_clusters_{2}.txt'.format(dir_out_merged, contig, sample)
    #             with open(pairwise_cluster_dup, 'r') as fpdup:
    #                 for line in fpdup:
    #                     line = line.strip().split('\t')
    #                     # pw_coord = line[0] # this set coordinates could be: 1. coordinates of ref sample (may be same as 2); 2. coordinates of individual sample (may be same as 1); 3. new coordinates inferred from all read pairs from ref and individual samples;
    #                     ref_coord = line[0]
    #                     ref_ct = line[1]
    #                     sp_coord = line[2]
    #                     if sp_coord != 'NA':
    #                         sp_coord = ','.join(sp_coord.split(',')[1:])
    #                     sp_ct = str(int(float(line[3])))
    #                     if ref_coord in dict_cnv_dup:
    #                             dict_cnv_dup[ref_coord][1].setdefault(sample, [sp_coord, sp_ct])
    #                     # else:
    #                     #         if ref_coord == 'NA' and pw_coord in dict_cnv_dup: # there are cases when ref coordinate is 'NA' due to matching criteria not met with individual CNV, but the individual coordinate is the same as that of the merged sample. We would like to recue such macthes.
    #                     #             dict_cnv_dup[pw_coord][1].setdefault(sample, [sp_coord, sp_ct])

    #         ### write dict entries (ref CNVs) that have at least one match across samples to output
    #         #### for deletions
    #         for deletion in dict_cnv_del:
    #             if set(list(dict_cnv_del[deletion][1][sample][0] for sample in dict_cnv_del[deletion][1])) != set(['NA']):
    #                 sample_cols = []
    #                 for sample in sample_order:
    #                     # if sample in dict_cnv_del[deletion][1]:
    #                     try:
    #                         sample_cols += dict_cnv_del[deletion][1][sample]
    #                     except:
    #                         print('del', deletion, contig, sample)
    #                     # else:
    #                     #     sample_cols += ['NA', '0'] # cases of CHTC running error will cause some sample runs being aborted. will be deleted
    #                 cols = [contig, ','.join(deletion.split(',')[1:]), dict_cnv_del[deletion][0]] + sample_cols
    #                 # print(cols)
    #                 fgcd.write('\t'.join(cols) + '\n')

    #         #### for duplications
    #         for duplication in dict_cnv_dup:
    #             if set(list(dict_cnv_dup[duplication][1][sample][0] for sample in dict_cnv_dup[duplication][1])) != set(['NA']):
    #                 sample_cols = []
    #                 for sample in sample_order:
    #                     # if sample in dict_cnv_dup[duplication][1]:
    #                     try:
    #                         sample_cols += dict_cnv_dup[duplication][1][sample]
    #                     except:
    #                         print('dup', duplication, contig, sample)
    #                     # else:
    #                     #     sample_cols += ['NA', '0'] # cases of CHTC running error will cause some sample runs being aborted. will be deleted
    #                 cols = [contig, ','.join(duplication.split(',')[1:]), dict_cnv_dup[duplication][0]] + sample_cols
    #                 # print(cols)
    #                 fgce.write('\t'.join(cols) + '\n')

    #         # 5.3. step 3: obtain allele frequency for reference CNVs in each population (using the Schrider method)
    #         ## generate input CNV coordinates for each contig of each sample, in order to parallelly count read pairs over breakpoints
    #         for sample in sample_order:
    #             found = 0 
    #             ### for deletions
    #             bed_del_sp_contig = '{0}/global_CNV_AF/del_{1}_{2}.bed'.format(dir_out, sample, contig)
    #             with open(bed_del_sp_contig, 'w') as fin_del_sp_contig:
    #                 for deletion in dict_cnv_del:
    #                     if set(list(dict_cnv_del[deletion][1][sample][0] for sample in dict_cnv_del[deletion][1])) != set(['NA']):
    #                         sample_coords = dict_cnv_del[deletion][1][sample][0]
    #                         sample_count = dict_cnv_del[deletion][1][sample][1]
    #                         if sample_coords != 'NA':
    #                             start = str(int(sample_coords.split(',')[0]) - 1) # transfer coordinate from 1-based SAM to 0-based BED 
    #                             end = sample_coords.split(',')[1]
    #                             cols = [contig, start, end, sample_count]
    #                             fin_del_sp_contig.write('\t'.join(cols) + '\n')
    #                             found = 1

    #             ### for duplications
    #             bed_dup_sp_contig = '{0}/global_CNV_AF/dup_{1}_{2}.bed'.format(dir_out, sample, contig)
    #             with open(bed_dup_sp_contig, 'w') as fin_dup_sp_contig:
    #                 for duplication in dict_cnv_dup:
    #                     if set(list(dict_cnv_dup[duplication][1][sample][0] for sample in dict_cnv_dup[duplication][1])) != set(['NA']):
    #                         sample_coords = dict_cnv_dup[duplication][1][sample][0]
    #                         sample_count = dict_cnv_dup[duplication][1][sample][1]
    #                         if sample_coords != 'NA':
    #                             start = str(int(sample_coords.split(',')[0]) - 1) # transfer coordinate from 1-based SAM to 0-based BED
    #                             end = sample_coords.split(',')[1]
    #                             cols = [contig, start, end, sample_count]
    #                             fin_dup_sp_contig.write('\t'.join(cols) + '\n')
    #                             found = 1

    #             if found: # no need to submit a job, if no deletions and duplications were found within the contig of the sample
    #                 ## write a dag to calculate allele frequency for CNVs of each population
    #                 scale = float(dict_contig_len[contig]/dict_contig_len['NC_050699.1'])
    #                 request_disk = basic_disk + int(1024**2*scale*0.8)
    #                 request_memory = basic_memory + int(100*scale*10)
    #                 request_cpus = 2

    #                 job = 'JOB 3_calc_cnv_af_{0}_{1} {2}/3_calc_cnv_af.sub\n'.format(contig, sample, dir_scripts)
    #                 vars = '''VARS 3_calc_cnv_af_{0}_{1} \
    #                     contig="{0}" sample="{1}" conda_pack_squid="{2}" path_in_stage="{3}" job="$(JOB)" \
    #                     dir_out="{4}" dir_scripts="{5}" bed_del_sp_contig="{6}" bed_dup_sp_contig="{7}" \
    #                     request_disk="{8}" request_memory="{9}" request_cpus="{10}"\n'''.format(
    #                     contig, sample, file_conda_pack_squid, dir_in_stage,
    #                     dir_out, dir_scripts, bed_del_sp_contig, bed_dup_sp_contig,
    #                     request_disk, request_memory, request_cpus)
    #                 f.writelines([job, vars])

    # 5.4. step 4: filter CNVs below a minor allele frequency threshold
    header = '\t'.join(['contig', 'coordinates_ref', 'count_ref']) \
        + '\t' + '\t'.join(list('coordinates_{0}\tMAF_{0}'.format(sample) for sample in sample_order)) + '\n'
    with open(file_global_cluster_del_af, 'w') as fdelaf, open(file_global_cluster_dup_af, 'w') as fdupaf:
        fdelaf.write(header)
        fdupaf.write(header)

        ## build a dict of chromosome - reference CNV coordinates - sample CNV coordinates & sample CNV MAF
        ### for deletions
        dict_cnv_del_af = {}
        dict_cnv_sp2ref_del = {} # dict to link sample CNV coordinates to its matched reference CNV coordinates
        with open(file_global_cluster_del, 'r') as fdel:
            i = 0
            for line in fdel:
                if i > 0:
                    line = line.strip().split('\t')
                    contig = line[0]
                    ref_coord = line[1]
                    sample_coords = list(line[3:][x*2] for x in range(len(sample_order)))
                    dict_cnv_del_af.setdefault(contig, {})
                    dict_cnv_del_af[contig][ref_coord] = {sample: [sample_coords[sample_order.index(sample)], 0] for sample in sample_order}

                    dict_cnv_sp2ref_del.setdefault(contig, {})
                    for j in range(len(sample_coords)):
                        if sample_coords[j] != 'NA':
                            dict_cnv_sp2ref_del[contig].setdefault(sample_order[j], {})
                            dict_cnv_sp2ref_del[contig][sample_order[j]][sample_coords[j]] = ref_coord
                i += 1

        for contig in list_contig_sorted:
            for sample in sample_order:
                file_contig_sample_af = '{0}/del_af_{1}_{2}.txt'.format(dir_out_global, sample, contig)
                if os.path.exists(file_contig_sample_af):
                    with open(file_contig_sample_af, 'r') as fspaf:
                        for line in fspaf:
                            line = line.strip().split('\t')
                            sample_coord = ','.join([str(int(line[1]) + 1), line[2]])
                            ref_coord = dict_cnv_sp2ref_del[contig][sample][sample_coord]
                            np = int(line[3])
                            q = int(line[4])
                            del_af = np/(np + q/2) # formula for AF calculation: p = Np/(Np + (q/2))
                            dict_cnv_del_af[contig][ref_coord][sample] = [sample_coord, del_af]

            if contig in dict_cnv_del_af:
                for ref_coord in dict_cnv_del_af[contig]:
                    del_afs = list(dict_cnv_del_af[contig][ref_coord][sample][1] for sample in sample_order)
                    mean_maf = min(mean(del_afs), 1 - mean(del_afs))
                    if mean_maf >= maf_cutoff:
                        sample_cols = []
                        for sample in sample_order:
                            sample_cols += [dict_cnv_del_af[contig][ref_coord][sample][0], str(dict_cnv_del_af[contig][ref_coord][sample][1])]
                        cols = [contig, ref_coord] + sample_cols
                        # print(cols)
                        fdelaf.write('\t'.join(cols) + '\n')


        ### for duplications
        dict_cnv_dup_af = {}
        dict_cnv_sp2ref_dup = {} # dict to link sample CNV coordinates to its matched reference CNV coordinates
        with open(file_global_cluster_dup, 'r') as fdup:
            i = 0
            for line in fdup:
                if i > 0:
                    line = line.strip().split('\t')
                    contig = line[0]
                    ref_coord = line[1]
                    sample_coords = list(line[3:][x*2] for x in range(len(sample_order)))
                    dict_cnv_dup_af.setdefault(contig, {})
                    dict_cnv_dup_af[contig][ref_coord] = {sample: [sample_coords[sample_order.index(sample)], 0] for sample in sample_order}

                    dict_cnv_sp2ref_dup.setdefault(contig, {})
                    for j in range(len(sample_coords)):
                        if sample_coords[j] != 'NA':
                            dict_cnv_sp2ref_dup[contig].setdefault(sample_order[j], {})
                            dict_cnv_sp2ref_dup[contig][sample_order[j]][sample_coords[j]] = ref_coord
                i += 1

        for contig in list_contig_sorted:
            for sample in sample_order:
                file_contig_sample_af = '{0}/dup_af_{1}_{2}.txt'.format(dir_out_global, sample, contig)
                if os.path.exists(file_contig_sample_af):
                    with open(file_contig_sample_af, 'r') as fspaf:
                        for line in fspaf:
                            line = line.strip().split('\t')
                            sample_coord = ','.join([str(int(line[1]) + 1), line[2]])
                            ref_coord = dict_cnv_sp2ref_dup[contig][sample][sample_coord]
                            np = int(line[3])
                            npq = int(line[4])
                            if npq > 0:
                                dup_af = np/(npq/2) # formula for AF calculation: p = Np/(Np+q/2)
                            else:
                                dup_af = 0
                            dict_cnv_dup_af[contig][ref_coord][sample] = [sample_coord, dup_af]

            if contig in dict_cnv_dup_af:
                for ref_coord in dict_cnv_dup_af[contig]:
                    dup_afs = list(dict_cnv_dup_af[contig][ref_coord][sample][1] for sample in sample_order)
                    if max(dup_afs) > 1:
                        dup_afs = list(dup_af/max(dup_afs) for dup_af in dup_afs) # if there is at least one duplication frequency > 1, normalize all population AF by the maximum AF
                    mean_maf = min(mean(dup_afs), 1 - mean(dup_afs))
                    if mean_maf >= maf_cutoff: # filter CNVs by MAF
                        sample_cols = []
                        for sample in sample_order:
                            sample_cols += [dict_cnv_dup_af[contig][ref_coord][sample][0], str(dup_afs[sample_order.index(sample)])]
                        cols = [contig, ref_coord] + sample_cols
                        # print(cols)
                        fdupaf.write('\t'.join(cols) + '\n')
