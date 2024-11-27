#output matrix
#!/usr/bin/python3

import os
import sys
import numpy as np
import pandas as pd
from collections import defaultdict
import subprocess

#for file in /data01/home/zhaogefei/eqtl/pbmc_1k/snp_out filtered snps
snp_folder = "~/snp_out/"
output_folder = "~/snp_out_filtered/"
os.makedirs(output_folder, exist_ok=True)
onlyfiles = [f for f in os.listdir(snp_folder) if f.endswith(".vcf")]
for f in onlyfiles:
    path_to_file = os.path.join(snp_folder, f)
    #out to _filtered_pass.vcf
    cmd=f"gatk VariantFiltration -R hg38.fa -V {path_to_file} --filter-name FS -filter 'FS > 30.0' --filter-name QD -filter 'QD < 2.0' -O {output_folder}{f.replace('.vcf', '_filtered_pass.vcf')}"
    subprocess.run(cmd, shell=True, check=True)



input_folder = "~/snp_out_filtered/"
output_folder = "~/matrix/"


suffix = "_filtered_pass.vcf"
onlyfiles = [f for f in os.listdir(input_folder) if f.endswith(suffix)]


min_shared_snv = 10
snv_dict = defaultdict(list)
sample_list = []


for f in onlyfiles:
    path_to_file = os.path.join(input_folder, f)
    sample_name = f.split("_")[0]
    sample_list.append(sample_name)
    with open(path_to_file, 'r') as current_file:
        for line in current_file:
            if not line.startswith("#"):
                chrom, position = line.split()[:2]
                snv_index = f"{chrom}__{position}"
                snv_dict[snv_index].append(sample_name)


snv_filtered = {snv: snv_dict[snv] for snv in snv_dict if len(snv_dict[snv]) >= min_shared_snv}
snv_filtered_name_list = list(snv_filtered.keys())


snv_matrix = np.zeros((len(snv_filtered_name_list), len(sample_list)), dtype=int)
snv_index_map = {snv: idx for idx, snv in enumerate(snv_filtered_name_list)}
sample_index_map = {sample: idx for idx, sample in enumerate(sample_list)}

for snv, samples in snv_filtered.items():
    for sample in samples:
        snv_matrix[snv_index_map[snv], sample_index_map[sample]] = 1

snv_df = pd.DataFrame(snv_matrix, index=snv_filtered_name_list, columns=sample_list)
snv_df.index.name = "SNVid"

os.makedirs(output_folder, exist_ok=True)
path_to_saved_snv_df = os.path.join(output_folder, "SNV_matrix.csv")
snv_df.to_csv(path_to_saved_snv_df)
