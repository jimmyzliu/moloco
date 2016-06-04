# multi-coloco

# input files - separate gwas summary files
# can specify region?
# each one with "SNP CHR BP BETA SE P"

# Then do:
# 1. merge summary stats
# 2. calc bfs
# 3. rbfs and ppas for each configuration

# test region
# chr15:78516053-80860978

import gzip
import itertools
import math
import os
import numpy as np
import string
import sys
import time
from scipy import stats

def wakefield_bf(beta,se,p,w = 0.1):
  # calculate Wakefield bayes factor
  z = stats.norm.isf(p/2)
  r = w/(se**2 + w)
  bf = math.sqrt(1-r) * math.exp(z**2/2*r)
  return bf

def make_sigma(cor_mat,v,w):
  sigma = np.diag(v+w)
  for i in range(len(v)):
    for j in range(i+1,len(v)):
      c = cor_mat[i,j]
      sigma[i,j] = c * math.sqrt((v[i] + w[i]) * (v[j] + w[j]))
      sigma[j,i] = c * math.sqrt((v[i] + w[i]) * (v[j] + w[j]))
  return np.matrix(sigma)

def is_pos_def(A):
  #return np.all(np.linalg.eigvals(A) > 0)
  try:
    np.linalg.cholesky(A)
  except:
    return False
  return True

def prod(ls):
  prod = 1.0
  for i in ls:
    prod = prod * i
  return prod

def get_psd(A):
  diag = np.diag(A)
  detA = np.linalg.det(A)
  while not is_pos_def(A):
    #while detA <= 0:
    A = 0.999 * A
    for i in range(np.shape(A)[0]):
      A[i,i] = diag[i]
    #detA = np.linalg.det(A)
  return A

def gather_assocs(in_files,chrom,start,stop):
  # only include snps common to all traits
  # calculate wakefield bf
  # assoc_dict[snp][file1] = [snp,chrom,bp,beta,se,p,bf]
  n_files = len(in_files)
  assoc_dict = {}
  start = int(start)
  stop = int(stop)
  # summary stats and wakefield bfs
  for file in in_files:
    if file.split(".")[-1] == "gz":
      f = gzip.open(file,'r')
      header = gzip.open(file,'r').readline().split()
    else:
      f = open(file,'r')
      header = open(file,'r').readline().split()
    snp_i = header.index("SNP")
    chr_i = header.index("CHR")
    bp_i = header.index("BP")
    p_i = header.index("P")
    se_i = header.index("SE")
    if "OR" in header:
      is_odds = True
      or_i = header.index("OR")
    if "BETA" in header:
      is_odds = False
      beta_i = header.index("BETA")
    
    for i in f:
      line = i.split()
      if line[chr_i] != chrom:
        continue
      if int(line[bp_i]) >= start and int(line[bp_i]) <= stop:
        [snp,chrom,bp,se,p] = [line[j] for j in [snp_i,chr_i,bp_i,se_i,p_i]]
        if is_odds:
          beta = math.log(float(line[or_i]))
        bf = wakefield_bf(float(beta),float(se),float(p))
        if not snp in assoc_dict:
          assoc_dict[snp] = {file:[snp,chrom,bp,float(beta),float(se),float(p),bf]}
        else:
          assoc_dict[snp][file] = [snp,chrom,bp,float(beta),float(se),float(p),bf]
  bad_snps = [snp for snp in assoc_dict if len(assoc_dict[snp]) != n_files]
  for snp in bad_snps:
    del assoc_dict[snp]
  print "Found " + str(len(assoc_dict)) + " SNPs common to all " + str(n_files) + " traits in chr" + chrom + ":" + str(start) + "-" + str(stop)
  return assoc_dict

def adjust_bfs(assoc_dict,configs,n_files,in_files,overlap):
  # adjusted bfs
  # gene correlation matrix
  # if assuming no overlap, cor_mat is identity matrix
  cor_mat = np.identity(n_files)
  if overlap:
    for i in range(n_files):
      for j in range(i+1,n_files):
        beta1 = [assoc_dict[snp][in_files[i]][3] for snp in assoc_dict]
        beta2 = [assoc_dict[snp][in_files[j]][3] for snp in assoc_dict]
        cortest = stats.pearsonr(beta1,beta2)
        if cortest[1] >= 0.01:
          cor_mat[i,j] = 0.0
          cor_mat[j,i] = 0.0
        else:
          cor_mat[i,j] = cortest[0]
          cor_mat[j,i] = cortest[0]
    cor_mat = np.matrix(cor_mat)
  # configs - different adj_bf depending on configs
  # traits abcde...
  # get adj_bf for each config combo
  adj_bf_dict = {config:{snp:0.0 for snp in assoc_dict} for config in configs}
  for snp in assoc_dict:
    v = np.array([assoc_dict[snp][in_file][4]**2 for in_file in in_files])
    betas = np.array([assoc_dict[snp][in_file][3] for in_file in in_files])
    ses = np.array([assoc_dict[snp][in_file][4] for in_file in in_files])
    means = np.array([0.0 for i in range(n_files)])
    w_0 = np.array([0.0 for i in range(n_files)])
    sigma_h0 = make_sigma(cor_mat,v,w_0)
    #
    for config in configs:
      w = [0.0 for i in range(n_files)]
      for i in config:
        w[string.ascii_lowercase.index(i)] = 0.1
      sigma_h1 = make_sigma(cor_mat,v,w)
      if not is_pos_def(sigma_h1):
        sigma_h1 = get_psd(sigma_h1)
      adj_bf = stats.multivariate_normal.pdf(betas, means, sigma_h1) / stats.multivariate_normal.pdf(betas,means,sigma_h0)
      adj_bf_dict[config][snp] = adj_bf
  return adj_bf_dict

def get_final_configs(configs,n_files):
  # list of all possible configurations
  combs = configs[:]
  for iter in range(1,n_files+1):
    #print iter
    new_configs = []
    for config in configs:
      if len(config.split(",")) != iter:
        continue
      for comb in combs:
        overlap = 0
        raw_comb = "".join(comb.split(","))
        raw_config = "".join(config.split(","))
        # 1. check to see no numbers in comb and config overlap
        for i in raw_comb:
          if i in raw_config:
            overlap = 1
            break
        # 2. check that length is less than n_genes
        if len(raw_comb + raw_config) > n_files+1:
          overlap = 1
          continue
        # check that config hasn't been seen
        new_config = sorted([comb,config])
        if not new_config in new_configs and overlap == 0:
          new_configs.append(new_config)
    configs = configs + [",".join(i) for i in new_configs]
  final_configs = []
  for i in configs:
    config = ",".join([k for k in sorted([j for j in i.split(",")])])
    if not config in final_configs:
      final_configs.append(config)
  return final_configs


def do_moloco(adj_bf_dict,final_configs,configs,priors):
  single_bfs = {config:0.0 for config in configs}
  for config in configs:
    n_snps = len(adj_bf_dict[config])
    prior = priors[len(config)-1]
    trait_bfs = []
    for snp in adj_bf_dict[config]:
      trait_bfs.append(adj_bf_dict[config][snp])
    single_bfs[config] = priors[len(config)-1] * sum(trait_bfs)
    
  config_dict = {i:0.0 for i in final_configs}
  for config in final_configs:
    if config in single_bfs:
      config_dict[config] = single_bfs[config]
      continue
    # left side
    left_trait_bfs = 1.0
    for i in config.split(","):
      left_trait_bfs = left_trait_bfs * single_bfs[i]
    # right side
    prior = prod([priors[len(j)-1] for j in config.split(",")]) / priors[len("".join(config).replace(",",""))-1]
    right_trait_bfs = prod([priors[len(j)-1] for j in config.split(",")]) * single_bfs["".join(sorted(config.replace(",","")))]
    config_bf = left_trait_bfs - right_trait_bfs
    config_dict[config] = config_bf
    
  config_ppas = {config:config_dict[config]/sum(config_dict.values()) for config in final_configs}
  moloco_stats = {}
  for config in config_dict:
    moloco_stats[config] = [config_dict[config],config_ppas[config]]
  return moloco_stats

def moloc_iter(in_files,chrom,start,stop,priors,out_file,overlap):
  # calc bf and consolidate sum stats
  n_files = len(in_files)
  assoc_dict = gather_assocs(in_files,chrom,start,stop)
  if len(assoc_dict) == 0:
    "Moving on to the next region"
    return True
  # single configs
  num2alpha = dict(zip(range(0, 26), string.ascii_lowercase))
  configs = [num2alpha[i] for i in range(n_files)]
  for i in range(1,n_files):
    n_causal = i + 1
    sel = "".join([num2alpha[k] for k in [j for j in range(n_files)]])
    configs = configs + ["".join(j) for j in list(itertools.combinations(sel,n_causal))]
  
  adj_bf_dict = adjust_bfs(assoc_dict,configs,n_files,in_files,overlap)
  # all configs
  final_configs = get_final_configs(configs,n_files)
  moloco = do_moloco(adj_bf_dict,final_configs,configs,priors)
  
  out_path = out_file + "." + chrom + "." + str(start) + "." + str(stop) + ".moloco"
  print "Great success! Writing results to " + out_path  
  write_out = open(out_path,'wa')
  print >>write_out, "config logBF PP"
  for config in final_configs:
    print >>write_out, config + " " + str(math.log(moloco[config][0])) + " " + str(moloco[config][1])


def main():
  # chrom = 15
  # start = 78516053
  # stop = 80860978
  #in_files = "lung_cancer.gwax.assoc.gz,bronchitis.gwax.assoc.gz,heart_disease.gwax.assoc.gz"
  start_time = time.time()
  print "\nMOLOCO"
  print "Version 0.1"
  print "Direct complaints to: jliu@nygenome.org\n"
  args = sys.argv[1:]
  do_bed = False
  overlap = True
  for i in range(len(args)):
    if args[i] == "--stats":
      in_files = args[i+1]
      print "Summary statistics: " + in_files
      in_files = in_files.split(",")
      n_files = len(in_files)
      for in_file in in_files:
        if not os.path.exists(in_file):
          print "Error: cannot find " + in_file
          sys.exit()
    if args[i] == "--chr":
      chrom = args[i+1]
      print "Chromosome: " + chrom
    if args[i] == "--from":
      start = args[i+1]
      print "Start position: " + start
    if args[i] == "--to":
      stop = args[i+1]
      print "Stop position: " + stop
    if args[i] == "--priors":
      priors = args[i+1].split(",")
      print "Priors: " + ",".join(priors)
      priors = [float(j) for j in priors]
    if args[i] == "--out":
      out_file = args[i+1]
      print "Output file: " + out_file
    if args[i] == "--no-overlap":
      overlap = False
    if args[i] == "--bed":
      bed_file = args[i+1]
      do_bed = True
      print "Bed file: " + bed_file
      if not os.path.exists(bed_file):
        print "Error: cannot find " + bed_file
        sys.exit()
  print ""
  if not 'in_files' in locals():
    print "Error: cannot find --stats\nBummer"
    sys.exit()
  if not 'chrom' in locals() and not 'bed_file' in locals():
    print "Error: cannot find --chr\nBummer"
    sys.exit()
  if not 'start' in locals() and not 'bed_file' in locals():
    print "Error: cannot find --from\nBummer"
    sys.exit()
  if not 'stop' in locals() and not 'bed_file' in locals():
    print "Error: cannot find --to\nBummer"
    sys.exit()
  if not overlap:
    print "Assuming no studies do not contain overlapping samples"
  if not 'priors' in locals():
    # use default priors
    mu = [10**i for i in range(1,n_files)]
    priors = [1e-4] + [1e-4 / i for i in mu]
    print "No priors specified. Using default priors: " + ",".join([str(i) for i in priors])
  if not 'out_file' in locals():
    out_file = "moloc"
  if len(priors) != n_files:
    print "Error: Number of priors must equal total number of phenotypes"
    sys.exit()
  
  if n_files == 1:
    print "Error: only one set of association statistics found. Need more"
    sys.exit()
  if n_files > 6:
    print "Warning: MOLOCO works best for 6 or fewer traits. I mean, it will still work for more traits, but now might be a good time to go for lunch/take a long walk outside/evaluate life choices"
  
  for i in range(n_files):
    print string.ascii_lowercase[i] + " = " + in_files[i] 
  
  if do_bed:
    if 'chrom' in locals() or 'start' in locals() or 'stop' in locals():
      print "Found bed file. Ignoring --chr --from --to"
    regions = [i.split() for i in open(bed_file,'r')]
    for region in regions:
      [chrom,start,stop] = region
      moloc_iter(in_files,chrom,int(start),int(stop),priors,out_file,overlap)
  else:
    moloc_iter(in_files,chrom,start,stop,priors,out_file,overlap)
  
  end_time = time.time()
  elapse = str(end_time - start_time)
  print "\nJob done in " + elapse + " seconds"


if __name__ == "__main__":
  main()
