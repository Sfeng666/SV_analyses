# function: prepare input table for plotting distribitions of insert size for each sample
import os

# 0. parameters
# 0.1. absolute
rename_dict = '/Users/siyuansmac/bioinfo/project/suzukii_WGS/rename_samples/rename_dict.txt'

# 0.2. relative (to git repo)
path_metrics = 'CNV/insert_size/results'
path_plot = 'CNV/insert_size/plot'

# 0.3. local
out_table = '{0}/insert_size_distribution_all_samples.txt'.format(path_plot)
sample_ola = [] #  sample names of olazcuaga data
sample_pool = [] # sample names of pool lab data

# 1. group samples into two datasets
with open(rename_dict, 'r') as f:
  for line in f:
    line = line.strip().split('\t')
    old_name = line[0]
    new_name = line[1]
    if old_name != 'British-Columbia-BC':
        if old_name != new_name: # all pool lab samples were renamed
            sample_pool.append(new_name)
        else:
            sample_ola.append(new_name)
sample_order = sample_pool + sample_ola

# 2. extract distribution dataframe
with open(out_table, 'w') as f:
    header = '\t'.join(['insert_size', 'FR_read_count', 'sample']) + '\n'
    f.write(header)
    for sample in sample_order:
        sample_table = '{0}/{1}_insert_size_metrics.txt'.format(path_metrics, sample)
        with open(sample_table, 'r') as f_in:
            start_read = 0
            for line in f_in:
                line = line.strip().split('\t')
                if line[0] == 'insert_size':
                    start_read = 1
                    continue
                if start_read == 1 and line[0].isdigit():
                    # for i in range(int(line[1])):
                    #     write_line = '\t'.join([line[0], sample]) + '\n'
                    #     f.write(write_line) 
                    write_line = '\t'.join(line[:2] + [sample]) + '\n'
                    f.write(write_line)

