from optparse import OptionParser

# 0.1. help and usage 
usage = "usage: %prog [options] args"
description = '''Function: estimate genome-wide Reynolds Fst for CNVs '''
version = '%prog 01.03.2023'
parser = OptionParser(usage=usage,version=version, description = description)
parser.add_option("--wd",
                    action="store",
                    dest = "wd",
                    help = "working directory",
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
parser.add_option("--in_efs",
                    action="store",
                    dest = "in_efs",
                    help = "input path of effective sample sizes ((auto/X)) of all samples (.txt)",
                    metavar = "PATH")
parser.add_option("--in_chr_assign",
                    action="store",
                    dest = "in_chr_assign",
                    help = "input path of assignments of contigs to chromosomes (.txt)",
                    metavar = "PATH")
parser.add_option("--out_log",
                    action="store",
                    dest = "out_log",
                    help = "output path of the log file (.log)",
                    metavar = "PATH")  
(options,args) = parser.parse_args()

# 0.2. introduced variables
wd = options.wd
in_del_af = options.in_del_af
in_dup_af = options.in_dup_af
in_efs = options.in_efs
in_chr_assign = options.in_chr_assign
out_log = options.out_log

# 0.3. def functions
## function to calculate the Reynolds Fst
def Fst_reynolds(p1_afs, p2_afs, size1, size2):
    ### p1_afs and p2_afs are dicitonaries that include allele frequency spectrum in each population
    ### size 1 and size 2 are allelic samples size of each population (the number of auto/X chromosomes, instead of the number of individuals)
    ### simply implementing the Reynolds's formulation (Reynolds et al, 1983), while the variable names is indicating the location of terms in the formulations.

    ### the first term that is shared by both numerator (al) and denominator (al + bl)
    first_term = sum(list((p1_afs[x] - p2_afs[x])**2 for x in p1_afs))/2

    ### the expression at the numerator of the second term that is also shared by both numerator (al) and denominator (al + bl)
    second_term_num2 = size1*(1 - sum(list(x**2 for x in p1_afs.values()))) + size2*(1 - sum(list(x**2 for x in p2_afs.values())))

    ### implement al
    al_frac_num = ((size1 + size2)/2)*second_term_num2
    al_frac_denom = size1*size2*((size1 + size2)/2- 1)
    frac = al_frac_num/al_frac_denom
    al = first_term - frac
    if al < 0:
        al = 0
 
    ### implement al + bl
    albl_frac_num = (size1*size2 - (size1 + size2)/2)*second_term_num2
    albl_frac = albl_frac_num/al_frac_denom
    albl = first_term + albl_frac

    ### seperately return al and al + bl for Fst calculation (weighted average of sites within windows)
    return al, albl
        
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

# 2. read the input of median effective sample sizes of autosomes and X chromosome for each population
efs_sps_median = {chr: {} for chr in chrs}
with open(in_efs, 'r') as f:
    i = 0
    for line in f:
        line = line.strip().split('\t')
        if i > 0:
            sample = line[0]
            efs_auto = float(line[1])
            efs_X = float(line[2])
            efs_sps_median['auto'][sample] = efs_auto
            efs_sps_median['X'][sample] = efs_X
        i += 1

# 3. estimate genome-wide pairwise Fst by weighted averaging point Fst across all CNVs on autosome or X-chromosome
## define a function to estimate genome-wide pairwise Fst
def calc_genomewide_fst(in_af, cnv):
    with open(in_af, 'r') as f:
        i = 0
        for line in f:
            line = line.strip().split('\t')
            if i == 0:
                samples = list(x.split('_')[1] for x in line[2:] if 'MAF' in x)
                matrix_pw_fst = {chr: {'Fst_matrix_num': {sample: [0]*len(samples) for sample in samples},
                'Fst_matrix_denom': {sample: [0]*len(samples) for sample in samples}} for chr in chrs} # matrices of Fst numerators and denominators for autosome and X chromosome
                num_sites = {chr: 0 for chr in chrs}
            elif i > 0:
                contig = line[0]
                if contig in contig_assign:
                    chr = contig_assign[contig]
                    pos = line[1]
                    afs_cnv = list(float(line[2:][j*2 + 1]) for j in range(len(samples))) # allele frequency of CNVs
                    efs_sps = list(efs_sps_median[chr][sample] for sample in samples) # effective sample sizes of all populations
                    num_sites[chr] += 1
                    ### calculate pairwise Fst for each CNV
                    for sample1 in samples:
                        af1 = {'cnv': afs_cnv[samples.index(sample1)], 'ref': 1 - afs_cnv[samples.index(sample1)]}
                        size1 = efs_sps[samples.index(sample1)]
                        Fst_vect_num = []
                        Fst_vect_denom = []
                        for sample2 in samples:
                            af2 = {'cnv': afs_cnv[samples.index(sample2)], 'ref': 1 - afs_cnv[samples.index(sample2)]}
                            size2 = efs_sps[samples.index(sample2)]
                            Fst_frac = Fst_reynolds(af1, af2, size1, size2)
                            al = Fst_frac[0]
                            albl = Fst_frac[1]
                            Fst_vect_num.append(al)
                            Fst_vect_denom.append(albl)
                        matrix_pw_fst[chr]['Fst_matrix_num'][sample1] = list(x + y for x,y in zip(matrix_pw_fst[chr]['Fst_matrix_num'][sample1], Fst_vect_num)) # accumulate numerators of pairwise Fst for further weighted averaging across CNVs
                        matrix_pw_fst[chr]['Fst_matrix_denom'][sample1] = list(x + y for x,y in zip(matrix_pw_fst[chr]['Fst_matrix_denom'][sample1], Fst_vect_denom)) # accumulate denominators of pairwise Fst for further weighted averaging across CNVs
            i += 1

    ### write genome-wide pairwise Fst as a matrix table
    for chr in matrix_pw_fst:
        out_fst_matrix = '{0}/Fst_{1}_{2}.txt'.format(wd, cnv, chr)
        with open(out_fst_matrix, 'w') as f:
            header = '\t'.join([str(num_sites[chr])] + samples) + '\n'
            f.write(header)
            for sample in samples:
                Fst_vect = list(x/y for x,y in zip(matrix_pw_fst[chr]['Fst_matrix_num'][sample], matrix_pw_fst[chr]['Fst_matrix_denom'][sample])) # calculate genome-wide Fst as a weighted average across all CNVs within autusome or X-chromosome
                row = '\t'.join([sample] + list('{:.6f}'.format(x) for x in Fst_vect)) + '\n'
                f.write(row)

## 3.1. for duplications
calc_genomewide_fst(in_dup_af, 'dup')
## 3.2. for deletions
calc_genomewide_fst(in_del_af, 'del')