#!/usr/bin/env python3

import os
import sys
import getopt
import glob
import re
import pandas as pd
import torch
import tensorqtl
from tensorqtl import genotypeio, cis, trans
print(f"PyTorch {torch.__version__}")
print(f"Pandas {pd.__version__}")

def main():

	region   = "Caudate";
	ancestry = "EA";
	group    = "NCS_SCZ";
	age = "13"

	sample = region + "." + ancestry + "." + group + ".age" + age

	d_proj  = "/dcs04/lieber/ds2a/shanzida/Projects/RNA-Seq/TopMed_2022/mRNA/TensorQTL/Caudate/EA/NCS_SCZ/"

	#features = ["Gene", "Exon", "Jxn", "Tx"]
	features = ["Gene"]

	for feature in features:

		lfeature = feature.lower()

		d_snp  = d_proj + "/SNP/"
		d_res  = d_proj + feature + "/eQTLs/"
		if not os.path.exists(d_res):
			os.makedirs(d_res)

		f_expr = d_proj + "Exp/tran_expr." + feature + "." + sample + ".bed.gz"
		f_pos  = d_proj + "Exp/tran_expr." + feature + "." + sample + ".pos"
		f_covs = d_proj + "Data/covs_" + feature + "." + sample + ".txt"

		# load phenotypes and covariates
		phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(f_expr)
		covariates_df  = pd.read_csv(f_covs, sep="\t",  index_col=0)

		chrs = list(range(1, 23))
		chrs.append("X")

		for chr in chrs:

			f_raw  = d_snp  + "snp_gen." + sample + ".chr" + str(chr) + ".traw"
			f_res  = d_res  + "eqtls_" + lfeature + "." + sample
			f_map  = d_res  + "eqtls_" + lfeature + "." + sample + ".cis_qtl_pairs.chr" + str(chr) + ".map"

			raw = pd.read_csv(f_raw, sep="\t")
			raw = raw.set_index("SNP")

			variant_df = raw[["CHR", "POS"]]
			variant_df.columns = ["chrom", "pos"]

			## Write out SNP map
			variant_df.to_csv(f_map, sep="\t", index=False)

			variant_df = variant_df.assign(chrom = "chr" + variant_df["chrom"].astype(str))

			genotype_df = raw.filter(regex="Br")

			genotype_df.to_csv("/dcs04/lieber/ds2a/shanzida/Projects/RNA-Seq/TopMed_2022/mRNA/TensorQTL/Caudate/EA/NCS_SCZ/genotype.txt")
			variant_df.to_csv("/dcs04/lieber/ds2a/shanzida/Projects/RNA-Seq/TopMed_2022/mRNA/TensorQTL/Caudate/EA/NCS_SCZ/variant.txt")

			colnames = []
			for cn in genotype_df.columns:
				colnames.append(re.findall("Br[0-9]+", cn)[0])

			genotype_df.columns = colnames


			cis.map_nominal(genotype_df,
				variant_df,
				phenotype_df.loc[phenotype_pos_df['chr']=='chr' + str(chr)],
				phenotype_pos_df.loc[phenotype_pos_df['chr']=='chr' + str(chr)],	
				f_res,
				covariates_df=covariates_df,
				maf_threshold=0.05,
				window=500000,
				output_dir=d_res,
				run_eigenmt=True,
				write_top=True,
				write_stats=True
			)

						

if __name__ == "__main__":
	main()
