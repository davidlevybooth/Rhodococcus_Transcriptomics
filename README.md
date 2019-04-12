# Rhodococcus_Transcriptomics
Rhodococcus Transcriptomic Pipeline

This snakemake pipeline aids in downloading and processing of genomic and transcriptomic sequences related to the growth of Rhodococcus rhodochrous strain EP4 on 4-ethylphenol, 4-propylguaiacol and succinate. These data are used in support of Levy-Booth et al., and Fetherolf et al., publications and are not intended to be used as an independent Prokaryotic transcriptomics pipeline. 
To use:

1.	Download or clone repository:

git clone https://github.com/levybooth/Rhodococcus_Transcriptomics.git
 
2.	Create your conda environment (might take a while): 

conda env create --name RREP4 --file environment.yaml

3.	Run pipeline: 

snakemake --snakefile align_genome3.snakefile --configfile config2.yaml

Tools for trimming, alignment and quantification can be selected by editing the config file (config2.yaml)
