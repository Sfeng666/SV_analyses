# generate a dag file to resubmit failed jobs nodes (not recognized as failure by htcondor) due to python compatibility
import os

fail_list = "/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/results/rescue_step2_pyfail.txt"
dag = "/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/results/pipeline.dag"
rescue_dag = "/home/sfeng77/jobs/SV_analyses/CNV/detect_CNV/results/rescue_step2_pyfail.dag"
os.system("head -1 {0} > {1}".format(dag, rescue_dag))

distant_comb = []
everted_comb = []
with open(fail_list, 'r') as f:
    for line in f:
        contig = '_'.join(line.strip().split('/')[1].split('_')[:2])
        sample = line.strip().split('/')[1].split('_')[-1].replace('.txt', '')
        if 'distant' in line:
            distant_comb.append('|'.join([contig, sample]))
        elif 'everted' in line:
            everted_comb.append('|'.join([contig, sample]))

diff_set = set(everted_comb) - set(distant_comb)
rescue_set = set(everted_comb) & set(distant_comb) # only jobs that generate null results for both distant and everted clusters could be caused by python issue
print(len(diff_set))
print(len(rescue_set))

for comb in rescue_set:
    os.system('grep -w "2_compare_clusters_across_samples_{0}_{1}" {2} >> {3}'.format(comb.split('|')[0], comb.split('|')[1], dag, rescue_dag))
