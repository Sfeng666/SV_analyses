# use picard CollectInsertSizeMetrics to calculate metrics and distribution of insert size
import os

# 0. parameters
path_conda = '/home/sfeng77/anaconda3/bin/activate'
env_conda = 'WGS_analysis'
path_tar = '/staging/sfeng77/test_pipeline/out/'
path_out = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/results'
script = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/scripts/1_find_cluster_discordant_reads_pre.sh'
script_merge = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/scripts/1_find_cluster_discordant_reads_merge_pre.sh'
rename_dict = '/home/sfeng77/jobs/rename_samples/rename_dict_rm_lowqual.txt'
contig_assign = '/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/scripts/assignment_cor_amb.txt'

# 1. build a dictionary for renamed samples
dict_rename = {}
with open(rename_dict, 'r') as f:
  for line in f:
    line = line.strip().split('\t')
    old_name = line[0]
    new_name = line[1]
    dict_rename[old_name] = new_name

# 2. build a list of contigs
list_contigs = []
with open(contig_assign, 'r') as f:
  i = 0
  for line in f:
    line = line.strip().split('\t')
    if i != 0:
      contig = line[0]
      list_contigs.append(contig)
    i += 1

# # 3. prepare bam files of each contig of each sample
# for old_name in dict_rename:
#   path_out_sp = '{0}/{1}'.format(path_out, dict_rename[old_name])
#   os.system('mkdir -p {0}'.format(path_out_sp))
#   log = '{0}/log_{1}'.format(path_out_sp, dict_rename[old_name])
#   for contig in list_contigs:
#       os.system("bash {0} {1} {2} {3} {4} {5} {6} {7} >{8} 2>&1".format(script, path_conda, env_conda, path_tar, old_name, dict_rename[old_name], path_out_sp, contig, log))

# 4. prepare bam files of each contig across all samples
path_out_sp = '{0}/merged'.format(path_out)
os.system('mkdir -p {0}'.format(path_out_sp))
log = '{0}/log_merged'.format(path_out_sp)

for contig in list_contigs[:1]:
    os.system("bash {0} {1} {2} {3} {4} {5} {6} {7} >{8} 2>&1".format(script_merge, path_conda, env_conda, path_tar, rename_dict, path_out, path_out_sp, contig, log))
