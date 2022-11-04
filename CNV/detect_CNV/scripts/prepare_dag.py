# -*- coding: utf-8 -*-
################################################################################
##  *
##  *  Function: prepare dag file to run the CNV analyses pipeline]
##  *  Writer:Siyuan Feng
##  *  Mail: siyuanfeng.bioinfo@gmail.com
##  *  Version: 11.02.2022
##  *
################################################################################

import os

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
file_conda_pack_squid = 'http://proxy.chtc.wisc.edu/SQUID/sfeng77/suzukii_WGS/WGS_analysis.tar.gz'
file_dag = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/results/pipeline.dag'

# 1. build a dictionary for renamed samples
dict_rename = {}
with open(file_sample, 'r') as f:
  for line in f:
    line = line.strip().split('\t')
    old_name = line[0]
    new_name = line[1]
    dict_rename[old_name] = new_name
dict_rename['merged'] = 'merged'

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

# 3. build a list of contigs
list_contig = []
with open(file_contig, 'r') as f:
    i = 0
    for line in f:
        line = line.strip().split('\t')
        if i != 0:
            contig = line[0]
            list_contig.append(contig)
        i += 1

# 4. build a dictionary of contig size
dict_contig_len = {}
with open(file_contig_len, 'r') as f:
    for line in f:
        line = line.strip().split('\t')
        if not "#" in line[0]:
            contig = line[6]
            length = int(line[8])
            dict_contig_len[contig] = length 

# 5. Write the DAG input file (.dag)
with open(file_dag, 'w') as f: 
    # set configuration variables for this dag, to remove default limits on job submissions
    # note!! the assumption to allow unlimited job submissions, is there should not be too many job submissions at the same time.
    f.write('CONFIG {0}/unlimited.config\n'.format(dir_scripts))

    for old_name in dict_rename:
        sample = dict_rename[old_name]
        path_out_sp = '{0}/{1}'.format(dir_out, sample)
        os.system('mkdir -p {0}'.format(path_out_sp))

        for contig in list_contig: # test run on a large contig
            ## 5.0. determine computational resource by scaling contig size
            basic_disk = 1.5*1024**2
            basic_memory = 1024
            request_cpus = 2
            scale = float(dict_contig_len[contig]/dict_contig_len['NC_050699.1'])
            request_disk = basic_disk + int(1024**2*scale)
            request_memory = basic_memory + int(1000*scale)

            ## 5.1. estimate the largest possible bam file size by scaling contig size
            est_bam_size = 50*scale

            if sample == 'merged':
                request_disk = basic_disk + int(1024**2*scale*29)
                request_memory = basic_memory + int(1000*scale*29)
                est_bam_size = est_bam_size*29

            ## 5.2. step 1: find and cluster discordant read pairs
            ### file transfering method depends on estimated size of input bam file
            if est_bam_size <= 500: # cutoff of HTC on whether one should input file from staging
                job_1 = 'JOB 1_find_cluster_discordant_reads_{0}_{1} {2}/1_find_cluster_discordant_reads.sub\n'.format(sample, contig, dir_scripts)
                vars_1 = '''VARS 1_find_cluster_discordant_reads_{0}_{1} \
                    sample="{0}" contig="{1}" conda_pack_squid="{2}" insert_size_cutoff="{3}" insertSizeDiffCutoff="{4}" job="$(JOB)"\
                    dir_out="{5}" dir_scripts="{6}" in_bam_compressed="{5}/{0}/{0}_{1}.bam.genozip" software_pack="{7}"\
                    request_disk="{8}" request_memory="{9}" request_cpus="{10}"\n'''.format(
                    sample, contig, file_conda_pack_squid, dict_para[sample][0], dict_para[sample][1], 
                    dir_out, dir_scripts, dir_software_pack, 
                    request_disk, request_memory, request_cpus)
                if sample == 'merged':
                    pre_1 = 'SCRIPT PRE 1_find_cluster_discordant_reads_{0}_{1} {7}/1_find_cluster_discordant_reads_merge_pre.sh {2} {3} {4} {5} {6} {6}/{0} {1}\n'.format(
                        sample, contig, dir_conda, env_conda, dir_tar, file_sample, dir_out, dir_scripts) 
                else:
                    pre_1 = 'SCRIPT PRE 1_find_cluster_discordant_reads_{0}_{1} {7}/1_find_cluster_discordant_reads_pre.sh {2} {3} {4} {5} {0} {6}/{0} {1}\n'.format(
                        sample, contig, dir_conda, env_conda, dir_tar, old_name, dir_out, dir_scripts)

            elif est_bam_size > 500:
                job_1 = 'JOB 1_find_cluster_discordant_reads_{0}_{1} {2}/1_find_cluster_discordant_reads_large.sub\n'.format(sample, contig, dir_scripts)
                vars_1 = '''VARS 1_find_cluster_discordant_reads_{0}_{1} \
                    sample="{0}" contig="{1}" conda_pack_squid="{2}" insert_size_cutoff="{3}" insertSizeDiffCutoff="{4}" job="$(JOB)"\
                    dir_out="{5}" dir_scripts="{6}" in_bam_compressed="{5}/{0}/{0}_{1}.bam.genozip" software_pack="{7}"\
                    request_disk="{8}" request_memory="{9}" request_cpus="{10}" path_in_stage="{11}"\n'''.format(
                    sample, contig, file_conda_pack_squid, dict_para[sample][0], dict_para[sample][1], 
                    dir_out, dir_scripts, dir_software_pack, 
                    request_disk, request_memory, request_cpus,
                    dir_in_stage) # variables added later
                if sample == 'merged':
                    pre_1 = 'SCRIPT PRE 1_find_cluster_discordant_reads_{0}_{1} {7}/1_find_cluster_discordant_reads_merge_pre_large.sh {2} {3} {4} {5} {6} {6}/{0} {1} {8}\n'.format(
                        sample, contig, dir_conda, env_conda, dir_tar, file_sample, dir_out, dir_scripts, dir_in_stage) 
                else:
                    pre_1 = 'SCRIPT PRE 1_find_cluster_discordant_reads_{0}_{1} {7}/1_find_cluster_discordant_reads_pre_large.sh {2} {3} {4} {5} {0} {6}/{0} {1} {8}\n'.format(
                        sample, contig, dir_conda, env_conda, dir_tar, old_name, dir_out, dir_scripts, dir_in_stage)
            f.writelines([job_1, vars_1, pre_1])






