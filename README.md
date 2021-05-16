# CUT_n_TAG

Pre-/post-processing pipelines and results for CUT&TAG data generated by the [Neurogenomics Lab @ Imperial College London](https://www.neurogenomics.co.uk/).   


## Results  

[Links to all results](https://neurogenomics.github.io/CUT_n_TAG/scripts/nf-core_atacseq_results.html).
Each subheader is the unique ID of a given sequencing batch, assigned by the [Imperial BRC Genomics Facility](https://www.imperial.ac.uk/medicine/research-and-impact/facilities/genomics-facility/). 

### HK5M2BBXY  

**Description**: Initial test run of four samples (two H3k27ac + two H3k27ame3). Accidentally merged libraries across assay types (H3k27ac/H3k27ame3) during nf-core/atacseq run (will fix).  

#### [MultiQC report](https://neurogenomics.github.io/CUT_n_TAG/web/processed_data/HK5M2BBXY/multiqc/narrowPeak/multiqc_report.html)  

#### [ataqv report](https://neurogenomics.github.io/CUT_n_TAG/web/processed_data/HK5M2BBXY/bwa/mergedLibrary/ataqv/narrowPeak/html/html/index.html)  

#### [NF execution report](https://neurogenomics.github.io/CUT_n_TAG/web/processed_data/HK5M2BBXY/pipeline_info/execution_report.html)  

#### [NF execution timeline](https://neurogenomics.github.io/CUT_n_TAG/web/processed_data/HK5M2BBXY/pipeline_info/execution_timeline.html)  

#### [CUT&TAG vs. ENCODE](https://neurogenomics.github.io/CUT_n_TAG/scripts/visualize_peaks.html)  

Analyses comparing CUT&TAG H3k27ac data (Imperial) vs. ChIP-seq H3k27ac data (ENCODE), both in K562 cell-lines.  


### All samples

[Analysis](https://neurogenomics.github.io/CUT_n_TAG/processed_data/analysis/CnT_analysis.html)

[Additional QC](https://neurogenomics.github.io/CUT_n_TAG/processed_data/analysis/CnT_extra_QC.html)


## Pipelines

### [nf-core/atacseq](https://nf-co.re/atacseq)  

- Platform: nf-core (nextflow + singularity/docker)   

- Discussion on adapting this pipeline for CUT&RUN data.  

 
#### 1. Setup  containers on HPC   

- Docker isn't allowed on HPC by itself because it presents some security risk. Instead, follow [these instructions](https://osf.io/6w7f9/wiki/Computing/) to create a R-based Docker container (Rocker) inside a singularity container.  

- By default singularity bind mounts](https://singularity.lbl.gov/quickstart) `/home/$USER`, `/tmp`, and `$PWD` into your container at runtime.  

```
mkdir -p /rds/general/user/$USER/ephemeral/tmp/  
mkdir -p /rds/general/user/bms20/ephemeral/rtmp/ 
```

- On HPC, Rocker containers can be run through Singularity with a single command much like the native Docker commands, e.g. "singularity exec docker://rocker/tidyverse:latest R"

- By default singularity bind mounts](https://singularity.lbl.gov/quickstart) `/home/$USER`, `/tmp`, and `$PWD` into your container at runtime. 

- !IMPORTANT! You may need to change the path of "/rds/general//user/$USER/home/R/x86_64-redhat-linux-gnu-library/3.6/" to the actualy location of your R library.  

- Run Rocker within singularity  
```
singularity exec -B /rds/general/user/$USER/ephemeral/tmp/:/tmp,/rds/general/user/$USER/ephemeral/tmp/:/var/tmp,/rds/general/user/$USER/ephemeral/rtmp/:/rds/general/user/$USER/home/R/x86_64-redhat-linux-gnu-library/3.6/ --writable-tmpfs docker://rocker/tidyverse:latest R
``` 

#### 2. Download nf-core/atacseq container  

Now you can download the *nfcore/atacseq* singularity container via [DockerHub](https://hub.docker.com/r/nfcore/atacseq/)   

- This will download "atacseq_latest.sif" to your home directory.
`singularity pull docker://nfcore/atacseq:latest` 

- Copy this .sif file to the cacheDir specified in your nextflow config file.  
`scp ~/atacseq_latest.sif /rds/general/user/$USER/projects/neurogenomics-lab/live/.singularity-cache/` 

- Once you have the container downloaded, you can now specify it in the [`-profile`](download the singularity image outside of the pipeline and save in the same dir as the cacheDir path for the singularity option in the custom config file) flag in the main pipeline (see below).  
- More info on this process is on the [lab Wiki](https://osf.io/6w7f9/wiki/Computing/).

#### 3. Prepare nextflow config file  

The config file tells nextflow how to run on Imperial's HPC.    
- `module load nextflow`  
- Copy the config file to the expected location so HPC knows how to run nextflow properly:  
 `scp hpc_config $HOME/.nextflow/config`

#### 4. *Optional*: Register Nexflow Tower  

- Register with [nextflow-tower according to Combiz's instructions](https://combiz.github.io/scflow-manual/example-run.html#enable-nextflow-tower)  
to get real-time reports as the pipeline runs. Once registered, add the token to your config file.  

-  Run the nextflow pipeline. See [here](https://nf-co.re/atacseq/1.2.1/parameters) for all parameter options.  

#### 5. Download the singularity container  

- In theory, *nf-core/atacseq* should download the singularity automatically when it runs.  
However in practice, downloading it this way either takes waaayyy too long, and/or fails entirely. 

- Therefore, per Narun Fancy's recommendation *"download the singularity image outside of the pipeline and save in the same dir as the cacheDir path for the singularity option in the custom config file"*. 
`/rds/general/user/$USER/projects/neurogenomics-lab/live/.singularity-cache`  

- For more info on the `-profile` flag, see [here](https://nf-co.re/atacseq/1.2.1/usage#profile).  

#### 6. Finally run the pipeline!  
- `--input`: Path to [design file](https://nf-co.re/atacseq/1.2.1/usage#input).  
- `--genome`: Genome build your fasrq files are in.  
- `-profile`: Path to [container profile](https://nf-co.re/atacseq/1.2.1/usage#profile).   

`nextflow run nf-core/atacseq --input raw_data/HK5M2BBXY/design.csv --genome GRCh37 -r 1.2.1 -profile /rds/general/user/$USER/projects/neurogenomics-lab/live/.singularity-cache/atacseq_latest.sif` 




### [henipipe](https://pypi.org/project/henipipe/)  
- Platform: python  

### [CUTTAG_tutorial](https://yezhengstat.github.io/CUTTag_tutorial/)  
- Platform: workflowr (R + CLI)  

### [CUT&RUNTools](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1802-4)  
- Platform: CLI- 



## Documentation

Exercepts from the full [BRC Genome help page](https://imperial-genomics-facility.github.io/igf-pipeline-help/index.html)

### File name
Illumina uses the following file name convention for the output fastq files

For example: **samplename_S1_L001_R1_001.fastq.gz**   

- **samplename** : Name of the sample provided in the samplesheet  
- **S1** : Number of sample based on the sample order on the samplesheet  
- **L001** : Lane number of the flowcell  
- **R1** : The read. For e.g. R1 indicates Read 1 and R2 indicates Read 2 of a paired-end run  
- **001** : Its always 001  
- *.fastq.gz* : File extension. Its a gzipped fastq file  

Please check the Illumina BCL2Fastq documentation for more information.
