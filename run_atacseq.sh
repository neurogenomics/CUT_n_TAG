#!/bin/bash  

source ~/.bashrc  
module load nextflow  

export repo_dir=$HOME/neurogenomics/GitRepos/CUT_n_TAG
export project_id=HK5M2BBXY
mkdir -p $repo_dir/processed_data/$project_id


nextflow run nf-core/atacseq \
    --input $repo_dir/raw_data/$project_id/design_noindex2.csv \
    --genome GRCh37 \
    --narrow_peak \
    --outdir $repo_dir/processed_data/$project_id \
    -with-singularity $HOME/atacseq_latest.sif \
    -c $repo_dir/hpc_config \
    -r 1.2.1
