# use picard CollectInsertSizeMetrics to calculate metrics and distribution of insert size
import os

# 0. parameters
path_conda = '/home/sfeng77/anaconda3/bin/activate'
env_conda = 'WGS_analysis'
path_tar = '/staging/sfeng77/test_pipeline/out/'
path_out = '/home/sfeng77/jobs/SV_analyses/CNV/insert_size/results'
script = '/home/sfeng77/jobs/SV_analyses/CNV/insert_size/scripts/calc_insert_size.sh'
rename_dict = '/home/sfeng77/jobs/rename_samples/rename_dict.txt'

# 1. run calculcation of insert size
with open(rename_dict, 'r') as f:
  for line in f:
    line = line.strip().split('\t')
    old_name = line[0]
    new_name = line[1]
    os.system("bash {0} {1} {2} {3} {4} {5} {6} &".format(script, path_conda, env_conda, path_tar, old_name, new_name, path_out))