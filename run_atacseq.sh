#!/bin/bash  

# Must bashrc to load nextflow tower token
source ~/.bashrc  
module load nextflow  

export repo_dir=$HOME/neurogenomics/GitRepos/CUT_n_TAG
export project_id=HK5M2BBXY
mkdir -p $repo_dir/processed_data/$project_id


nextflow run nf-core/atacseq \
	--input $repo_dir/raw_data/$project_id/design_assays.csv \
	--genome GRCh37 \
	--outdir $repo_dir/processed_data/$project_id \
	-with-singularity $HOME/atacseq_latest.sif \
	-with-tower \
	-r 1.2.1 \
	-profile imperial \
	-c $repo_dir/hpc_config2 
#   --narrow_peak
