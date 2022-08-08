# Analysis code for Hustmyer et al. 2022 #
This repository contains original code used to analyze data for the Hustmyer et.  al. 2022 paper (see also extended methods of Hustmyer et. al 2022) including a fully portable and reproducible ChIP-seq analysis pipeline. Updated versions of the pipeline can be found at https://github.com/mikewolfe/ChIPseq_pipeline. Details about the design of pipeline and instructions on downloading can be found beginning with header 'What the pipeline does', which includes downloadable test files.  

All instructions to re-run ChIP-seq analysis can be found in the config file located at `config/config_hustmyer.yaml` The sample sheet, which indicates how each sample was post-processed (columns: 'postprocessing') is located at `pep/sample_sheet_hustmyer.csv`. If the sample was a ChIP sample, the sample sheet also indicates which named input sample was used to normalize the ChIP sample (column: 'input_sample') and the relevant genome for alignment (column: 'genome'). Single and paired end reads were distinguished in the 'bamCoverage' column (see alignment instructions).  

## Alignment

The genomes that each sample was aligned to were specified in the sample sheet (column: 'genome'). For deletion strains, custom bed files were made with the coordinates of the CDS for genes to be masked from each alignment. These .bed files can be found in the `bed_files` directory and were used together with the parent genome and the `workflow/scripts/combine_fasta.py` script to mask the specified regions in-place by replacing the sequence at those regions with Ns.

- For WT CFT073 strains, specified by the genome and NCBI deposit number NC_004431.1, the bed file containing annotations for each gene was parsed by the `workflow/scripts/parse_genbank.py` script in the pipeline with genome name NC_004431.1. The genbank file for CFT073 was parsed using the values `"gene old_locus_tag locus_tag"`. No genes were masked for this genome.
- For CFT073 *rfaH* deletions, specified by genome name NC_004431.1_dr, genes were masked from NC000431.1 using the .bed file: `bed_files/cft_rfaH.bed`.
- For CFT073 *hns* and *stpA* double deletions, specified by genome name NC_004431.1_2x, genes were masked from NC000431.1 using the .bed file: `bed_files/cft_hns_stpA.bed`
- For CFT073 *hns*, *stpA*, and *rfaH* triple deletions, specified by genome name NC_004431.1_3x, genes were masked from NC000431.1 using the .bed file: `bed_files/cft_hns_stpA_rfaH.bed`
- For CFT073 *stpA* deletions, specified by genome name NC_004431.1_ds, genes were masked from NC000431.1 using the .bed file: `bed_files/cft_stpA.bed`
- For CFT073 *stpA* and *rfaH* double deletions, specified by genome name NC_004431.1_ds_dr, genes were masked from NC000431.1 using the .bed file: `bed_files/cft_stpA_rfaH.bed`
- For CFT073 *hns* deletions, specified by genome name NC_004431.1_dh, genes were masked from NC000431.1 using the .bed file: `bed_files/cft_hns.bed`
- For *E. coli* RL3000 and MG1655, the NCBI accession no. U00096.3 with genome name U00096.3 was used. 

Paired-end reads were denoted in the sample sheet by the number of raw input fastq files for a given sample (i.e. containing entries in both filenameR1 and filenameR2) and the parameters used to specify how coverage calculations are performed for each sample in the bamCoverage_params column. For paired end reads these parameters are `--samFlagInclude 66 --extendReads` whereas For single-end ChIP-exo reads, the column was left blank. 

## Coverage and norm

After alignment, ChIP and input samples were processed using parameters specified in the `coverage_and_norm` module found in the supplied config file. 
* A 5 bp coverage calculation was specified with the tag `resolution: 5`
* All raw data was scaled by the median signal within a given sample (see extended methods) as specified with the tag `within: median`
* Not a number and infinity were dropped from ratio coverage so that bigWigs could be visualized in IGV using the specific tag: `dropNaNsandInfs: true`.
* RobustZ scaling was not used as specified by the tag: `RobustZ: false`.
* Paired end vs single end reads were distinguished by the column: "bamCoverage_params" in sample sheet (see details above)
* A pseudo count was not added to the coverage before normalizing by the median to prevent distortion of the data, as specified by `pseudocount: value: 0`
* For occupancy scaling RNAP ChIP data (see `postprocessing` and extended methods), the maxima for scaled data was determined by the average of 20 fixed tRNA genes. 
	* The coordinates for these genes can be found in `bed_files/top20fix.bed`. The maxima was scaled using both the postprocessing module and the tag `bwtools_fixed_scale` 
* For all occupancy scaled data, the minima was arbitrarily determined by taking the average of the 20 genes with the lowest coverage using the tag `bwtools_query_subtract`. 
	* To look for the 20 lowest genes, the .bed file for the entire genome was used with the tag `results/alignment/process_genbank/NC_004431.1/NC_004431.1.bed`. To specify 20 genes, the tag `number_of_regions` with `value: 20` was used. 
* For some occupancy scaled data (H-NS and StpA, Fig. 6), the maxima was also arbitrarily determined by taking the average of the 20 genes with the highest coverage using the tag `bwtools_query_scale` with the same parameters as `bwtools_query_subtract`.

All manipulations of bigwig files were performed using the script `workflow/scripts/bwtools.py` with parameters described above and specified by the workflow rules in `workflow/rules/coverage_and_norm.smk` and the constraints implied by the config file `config/config_hustmyer.yaml` and sample metadata `pep/sample_sheet_hustmyer.csv`.

## Peak calling
RfaH bound regions were called using macs2 with the `peak_calling` module of the pipeline with the parameters `--broad`
* Although all data was peak caller analyzed for general understanding of the results, only RfaH regions were called using macs2 as detailed in extended methods, Fig. 2A, Fig. S2E, Supplemental Dataset 2   

## Variant calling

All CFT073 input samples were analyzed for variants using `breseq` with the `variant_calling` module against the NC_004431.1 reference genome. 
* These data are reported in Supplemental Dataset 1G and 1H. 

## Bwtools multicompare

Used in Fig. 1F and G to visualize WT H-NS and WT RNAP ChIP/Input in IGV. Generates average coverage bigWig file of WT and H-NS from 3 biological replicates, which can be opened in IGV to create figure.
* The model `mean_hns__wt` filters for the desired WT H-NS IP/Input replicates specifying the sample filter  `filter: compare_mod1 == "1" and not input_sample.isnull()` and the information in column `compare_mod1` in the sample sheet. 
	* Median ChIP/input coverage was used by specifying the appropriate input files with `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_ratio.bw"`
	* Average of 3 replicates was specified by the tag `operation_within: "mean"`
* The model `mean_rnap_wt` instructs similar processing, but for 3 biological replicates of RNAP ChIP in WT as specified with the filter and in sample sheet. 

## Postprocessing 

Various models were used to quantify signal per gene, transcription unit, RfaH regulated transcription unit, log2 ratio ChIP/Input, etc. Filters were used to specify which samples were processed with which models. 
Below is a description of each model and the relevant figure upon which each analysis was visualized. 
* Note, for occupancy scaling, many parameters were set in the `coverage_and_norm` of the config file.
* For signals per gene, the gene annotations from WT CFT073 were processed from the NC_004431.1 NCBI genbank as described above in the alignment section and specified by the output file `results/alignment/process_genbank/NC_004431.1/NC_004431.1.bed`.
		* All other relevant bed file locations are described below. 
		* Model descriptions are also included in the config file. See sample sheet columns for reference as to which samples underwent data processing from which models.
* Model 1 - Provide the coverage every 5 bp for the median scaled log2(ChIP/Input) at every gene in the genome, specified by the tag `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_log2ratio.bw"`. 
	* Tag `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_log2ratio.bw"` controls coverage calculations and scaling.
	* Genes are specified by the tag `regions: "results/alignment/process_genbank/NC_004431.1/NC_004431.1.bed"`. 
	* The tag `summarize: "identity" specifies full resolution (5 bp) read coverage.
		* Filter here applies to all samples (see sample sheet). 
		* Relevant figures and datasets: used for general understanding of data. 
* Model 2 - Provide the average signal per gene coordinate of the median scaled log2(ChIP/Input) at every gene in the genome, specified by the tag `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_log2ratio.bw"`. 
	* Tag `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_log2ratio.bw"` controls coverage calculations and scaling.
	* Genes are specified by the tag`regions: "results/alignment/process_genbank/NC_004431.1/NC_004431.1.bed"`. 
	* The tags `summarize: "single"` and `summary_func: "mean"` specify to calculate the average for each region.. The tag `frac_na: 0.20` specifies that if a gene region contains more than 20% Not a numbers, then the entire region will be reported as NaN.
		* Filter here applies to all samples (see sample sheet).
		* Relevant figures and datasets: Used for general understanding of the data.
* Model 3 - Same as model 1, but median scaled ChIP/input coverage per gene boundary without log2 transformation.
	* Tag `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_ratio.bw"` controls coverage calculations and scaling. 
	* Filter here applies to all samples (see sample sheet). 
	* Relevant figures and datasets: Used for general understanding of the data.
* Model 4 - Same as model 2, but median scaled ChIP/input average signal per gene without log2 transformation. 
	* Tag `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_ratio.bw"` controls coverage calculations and scaling. 
		* Filter here applies to all NC_004431.1 samples (see sample sheet). 
		* Relevant figures and datasets: Fig. 6B, Fig. S1A, Fig. S1D, Fig. S5A, B, G, H, Fig. S6B-C.
* Model 6 - Provides average coverage at gene boundaries that were arbitrarily occupancy scaled using the top 20 genes for each dataset as the maxima scaling factor and subtracting the average 20 lowest genes (see `coverage_and_norm` and `workflow/scripts/bwtools.py` for details).
	* The tag `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_ratio_querysub_queryscale.bw"` provides instructions for the data, where the data was first median scaled, then ratioed, then occupancy scaled using the querysub and queryscale paramters in `coverage_and_norm`
		* Filter here was applied to all NC_004431.1 samples, but only H-NS and StpA CFT073 ChIP-seq IP/Input ratios were analyzed
		* Relevant figures and datasets: Used for general understanding of the data and for dataset 1, Fig. 6C and D, Fig. S6E-H Supplemental Dataset 1E, 1F.
* Model 7 - look at the full coverage for all parsed transcription units in a log2 ratio, median normalize reads, no occupancy scale  
    * Tag `filesignature:"results/coverage_and_norm/bwtools_compare/%s_median_log2ratio.bw"` calculates coverage at the transciption unit boundaries. 
	* Tag `regions: "bed_files/parse_TU.bed"` directs to the boundaries of all transcription units. 
	* Tag `upstream: 300` indicates to calculate coverage 300 bp upstream of the 5' end specified in `bed_files/parse_TU.bed`
		* Filter here applies to all NC_004431.1 samples (see sample sheet), but only WT RNAP and WT H-NS coverage was used in the figure. 
		* Relevant figures: Fig. 1F, Fig. 1G
* Model 15 - Occupancy scale RNAP data per gene 
	* look at the mean coverage for all genes in a ratio, median normalize, then occupancy scale data based on top 20 fixed tRNA genes (for RNAP),20 lowest genes as minima
	* Uses the `bwtools_fixed_scale` tRNA genes for maxima scaling.
	* Tag `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_ratio_querysub_fixedscale.bw"` used to control coverage calculations.
		* Filter here applies to all NC_004431.1 samples (see sample sheet), but data only analyzed for WT RNAP and *-hns-stpA* RNAP occupancy data for Fig. S2B. 
		* Relevent figures: Fig. S2B
* Model 24 - look at the mean coverage for all genes with median/input samples from the hns chip data set in k12
	* Analyzes the IP/Input data per gene for RL3000 samples. 
	* Calculation of coverage specified by the tag `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_ratio.bw"`.
	* Regions for all genes in MG1655 specified by the tag `regions: "results/alignment/process_genbank/U00096.3/U00096.3.bed"`
	* Tag `summarize_func: "mean"` reports the average of the 5 bp coverage per region
		* Filter `filter: 'postprocess_mod24_1 == "1" and not input_sample.isnull()'` specifies H-NS ChIP data collected in RL3000 as specified on sample sheet. 
		* Relevant figures: Fig. S1B, Supplemental Dataset 1B, 1C
* Model 27 - look at the raw coverage IP for RfaH regulated TUs in WT, +/- 1800 nt
	* Calculation of coverage specified by the tag `filesignature: "results/coverage_and_norm/deeptools_coverage/%s_median.bw"`.    
	* Regions for RfaH regulated transcription units specified by the tag `regions: "bed_files/rfaH_regCtrlTUs.bed"`
	* Tag `upstream: 1800` and `downstream: 1800` specifies return of coverage 1800 nt downstream and upstream of specified transcription unit boundary. 
		* Filter `'postprocess_mod25 == "1" and not input_sample.isnull()'` refers to samples in sample sheet for production of relevant figures.
		* Relevant figures and datasets: Fig. 2B-E and Fig. S4, Fig. S6A, S6D. 
			* Each IP was then individually normalized within the TU window,
			using the average of the five lowest signals as background (set to 0) and the average of the five
			highest signals set to 1. Using Microsoft Excel, sequential rolling averages over 25 bp windows (at least 5 times) were
			then used to smooth signals, which also had the consequence of causing the highest peaks to be
			somewhat below 1.
* Model 28 - look at the IP/Input signal per gene for WT H-NS and RNAP at CDS
	* Goal was to calculate median scaled IP/Input coverage at coding regions in CFT073. 
	* Calculation of coverage specified by the tag `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_ratio.bw"`
	* Regions of coding genes specified by tag `regions: "bed_files/NC_004431.1_CDS.bed"`
		* Filter `'postprocess_mod28== "1" and not input_sample.isnull()'` specifies WT RNAP and H-NS in CFT073 ChIP/Input as specified in sample sheet. 
		* Relevant Figs.: Fig. 1B and C, Fig. 1E, Supplemental Dataset 1A, 1B, 1D
* Model 30 - look at the raw coverage ratio for RNAP occupancy at RfaH regulated TUs in all conditions, -300 nt upstream from the predicted transcription start site. 
    * Coverage calculation controlled by the tag `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_ratio_querysub_fixedscale.bw"`
	* Relevant RfaH regulated and control transcription unit boundaries for traveling ratios specified by the tag `regions: "bed_files/rfaH_regCtrlTUs.bed"`
	* Tag `upstream: 300` indicates to calculate coverage 300 bp upstream of the 5' end specified in `bed_files/rfaH_regCtrlTUs.bed`
		* Filter `filter: 'postprocess_mod30 == "1" and not input_sample.isnull()'` specifies all RNAP IPs and inputs in sample sheet.
		* Relevant figures and datasets: Figs. 3-5 Extended Dataset 3. 5 bp coverage plots for relevant TUs were averaged in 100 bp windows in excel. Traveling ratios were calculated in excel as detailed in extended methods and Fig. 3F.
* Model 31 - look at sigma coverage at RfaH regulated transcription units
	* Goal was to calculate Sigma70 IP/Input coverage in 5 bp windows across RfaH regulated transcription units in WT and *-rfaH* strains
	* Coverage calculation controlled by the tag `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_ratio.bw"`
	* Regions of RfaH regulated transcription units specified by the tag `regions: "bed_files/rfaH_regCtrlTUs.bed"`
	* Tag `upstream: 300` indicates to calculate coverage 300 bp upstream of the 5' end specified in `bed_files/rfaH_regCtrlTUs.bed`
		* Filter: `filter: 'postprocess_mod31 == "1" and not input_sample.isnull()'` specifies the relevant Sigma70 IPs in sample sheet. 
		* Relevant Figures: Figs. S5C-F
* Model 34 - look at raw ChIP-exo coverage of StpA and H-NS per gene
	* Goal was to calculate H-NS and StpA ChIP-exo signal per gene in MG1655
	* Coverage calculation specified by `filesignature: "results/coverage_and_norm/deeptools_coverage/%s_median.bw"`, as no input was deposited the raw ChIP was median normalized
		* `filter: 'postprocess_mod33 == "1"'` specifies the relevant ChIP-exo samples in the sample sheet. 
		*  Relevant figures and datasets: Fig. S1E
* Model 35 - 1.4 kb windows ChIP/Input for heat maps
	* Goal was to calculate the coverage of IP/Input in 1.4 kb windows across the entire CFT073 genome to generate heat maps. 
	* Coverage calculation was specified by the tag `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_ratio.bw"`
	* Regions to calculate average coverage specified by `regions: "bed_files/1400bp_windows.bed"` which specify the 1.4 kb tiling windows across the genome. 
		* Filter `filter: 'postprocess_mod35 == "1" and not input_sample.isnull()'` specifies the relevant samples in the sample sheet. 
		* Relevant figures and datasets: Fig. 1D, Fig. 2A (top heat map), Fig. S1C
* Model 36 - look at WT StpA, H-NS, and RNAP coverage in 1.4 kb windows
	* Goal was to calculate the occupancy scaled H-NS and StpA ChIP-seq across 1.4 kb windows.
	* Coverage calculation specified by `filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_ratio_querysub_queryscale.bw"`, which arbitrarily scales to maxima and minima of each replicate as described in `coverage_and_norm`	
	* Regions to calculate average coverage specified by `regions: "bed_files/1400bp_windows.bed"` which specify the 1.4 kb tiling windows across the genome. 
        * Filter `filter: 'postprocess_mod32 == "1" and not input_sample.isnull()'`
		* Relevant figures: Fig. S6A, D heat maps
* Model 37 - Median ratio IP/Input for RfaH rregulated TUs, 1000 bp US and DS
	* Goal was to calculate the Sigma70 and H-NS coverage at RfaH regulated transcription units in the *-hns-stpA* strains
	* Coverage calculation specified by the tag `'filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_ratio.bw"`
	* Regions to calculate average coverage specified by `regions:"bed_files/RfaH_2Xreg.bed"`
	* Tag `upstream: 1000` indicates to calculate coverage 1000 bp upstream of the 5' end specified in `"bed_files/RfaH_2Xreg.bed"`
		* Filter: `filter: 'postprocess_mod37 == "1" and not input_sample.isnull()'` specifies relevant Sigma70 ChIP samples in sample sheet
		* Relevant Figs: S3
* Model 38 - Median ratio IP/Input for RfaH regulated TUs, 8000
	* Goal was to calculate the Sigma70 and H-NS coverage at RfaH regulated transcription units in the *-hns-stpA* strains
	* Coverage calculation specified by the tag `'filesignature: "results/coverage_and_norm/bwtools_compare/%s_median_ratio.bw"`
	* Regions to calculate average coverage specified by `regions:"bed_files/RfaH_2Xreg.bed"`
	* Tag `upstream: 8000` indicates to calculate coverage 8000 bp upstream of the 5' end specified in `"bed_files/RfaH_2Xreg.bed"`
		* Filter: `filter: 'postprocess_mod37 == "1" and not input_sample.isnull()'` specifies relevant Sigma70 ChIP samples in sample sheet
		* Relevant Figs: S3

# What the pipeline does
Starting from raw fastq files this pipeline does the following:

- preprocessing to remove adapters and low quality sequences using [`cutadapt`](https://cutadapt.readthedocs.io/en/stable/)
  and [`trimmomatic`](http://www.usadellab.org/cms/?page=trimmomatic)
- alignment to (a) reference genome(s) using [`bowtie2`](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- read coverage calculation with a variety of normalization options within or between samples 
  with [`deeptools`](https://deeptools.readthedocs.io/en/latest/) and custom python scripts 
- peak calling using either [`macs2`](https://github.com/macs3-project/MACS) or a custom python 
  implementation of [`CMARRT`](https://github.com/mikewolfe/CMARRT_python)
- variant calling against a reference genome on ChIP input samples using
  [`breseq`](https://barricklab.org/twiki/bin/view/Lab/ToolsBacterialGenomeResequencing)


In addition this pipeline uses [`multiqc`](https://multiqc.info/) to compile the 
following quality control metrics into an interactive html report:
- Read QC using [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) both before and after preprocessing
- A large number of ChIP quality control metrics calculated by [`deeptools`](https://deeptools.readthedocs.io/en/latest/)

# Installing the pipeline

This pipeline has three dependencies:
- The package manager [`miniconda`](https://docs.conda.io/en/latest/miniconda.html)
- The workflow management system [`snakemake`](https://snakemake.readthedocs.io/en/stable/index.html)
- The API for working with Portable Encaspulated Projects (PEP) [`peppy`](http://peppy.databio.org/en/latest/)

`miniconda` can be installed following the installation instructions for your system [here](https://docs.conda.io/en/latest/miniconda.html).

Once `miniconda` is installed, both `snakemake` and `peppy` can be installed in their own environment easily using:
```
conda create -n ChIPseq_pipeline snakemake=5.24.2 peppy
conda activate ChIPseq_pipeline
```

**Note** If you are using a computational cluster that requires job management
software, you may want to install that with your environment as well.
For example, if you are using an `htcondor`-managed server you would
instead create your environment like so:
```
conda create -n ChIPseq_pipeline snakemake=5.24.2 peppy htcondor=8.9.5
conda activate ChIPseq_pipeline
```

Now you can pull the pipeline from github using:
```
git clone --recurse-submodules https://github.com/mikewolfe/ChIPseq_pipeline/
```

And you can change into the newly cloned `ChIPseq_pipeline` directory and test your installation with:
```
snakemake --use-conda --cores 10
```

Or if using a cluster with job management software you can run this
with an environment-specific
[profile](https://snakemake.readthedocs.io/en/v5.1.4/executable.html#profiles).
For example:
```
snakemake --use-conda --cores 10 --profile htcondor
```


This will run the entire pipeline using the provided test data consisting of small example fastqs heavily downsampled from real ChIP data. 

The first time you run the pipeline it will need to create dedicated `conda` environments for each module which will take some time. Afterwards, it will run quickly. For more information on using `conda` with `snakemake` including how to set things up to run offline check out the documentation [here](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management).

If everything runs smoothly you can then clean up and remove the results from the test data using:
```
snakemake clean_all --cores 1
```

# Running the pipeline

This pipeline uses `snakemake` to manage the workflow and familiarity with `snakemake` will help with getting the most out of the pipeline. Fortunately, `snakemake` has an excellent tutorial that can be found [here](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) if you are unfamiliar with this workflow management system.

## Input

This pipeline takes as input a [Portable Encaspulated Project (PEP)](http://pep.databio.org/en/latest/) which is essentially a `.csv` of samples and metadata together with a `.yaml` file allowing for extensions to the existing metadata.

The following required fields are needed for this pipeline
- `sample_name` - a unique identifier for each sample
- `filenameR1` - the base file name for read1 of a set of paired fastq files
- `filenameR2` - the base file name for read2 of a set of paired fastq files
- `file_path` - the path to where the files for a given sample live
- `input_sample` - which unique `sample_name` acts as the input for each extracted sample.
Input samples should leave this field blank.
- `genome` - what reference genome should this sample be aligned to

An example of a sample sheet is included at [pep/test_samples.csv](pep/test_samples.csv).

Additionally, the sample sheet can be augmented with a required `config.yaml` file. In the included test example this is used to replace the `file_path` field with a specific location. This example can be found at [pep/config.yaml](pep/config.yaml). 

The `pep/config.yaml` should be edited to point towards your `pep/samples.csv` file. If you create your `samples.csv` file with excel be sure to save it as `comma seperated values` not any of the other encodings for `.csv`.

## Configuration

The pipeline itself, including parameters controlling specific tools, is controlled by a `.yaml` file in [config/config.yaml](config/config.yaml). The included `config.yaml` has all possible options specified with comments describing what those options control.

## Rules

The pipeline is organized into modules each of which runs a specific task needed for ChIP-seq analysis.

- [workflow/rules/preprocessing.smk](workflow/rules/preprocessing.smk) includes rules for trimming raw reads for adapters and quality
- [workflow/rules/alignment.smk](workflow/rules/alignment.smk) includes rules for aligning samples to their reference genome
- [workflow/rules/coverage_and_norm.smk](workflow/rules/coverage_and_norm.smk) includes rules for calculating read coverage over the genome and performing within and between sample normalization
- [workflow/rules/peak_calling.smk](workflow/rules/peak_calling.smk) includes rules for calling ChIP-seq peaks
- [workflow/rules/quality_control.smk](workflow/rules/quality_control.smk) includes rules for performing summarizing quality control on the reads themselves and ChIP-seq specific quality control
- [workflow/rules/postprocessing.smk](workflow/rules/postprocessing.smk) includes rules for getting summaries of ChIP signals over specified regions or locations
- [workflow/rules/variant_calling.smk](workflow/rules/variant_calling.smk) includes rules for checking for mutations against the reference genome. Typically run on ChIP input samples. Will take awhile to run.

Each of these rules can be run individually using:
```
snakemake run_module_name --use-conda --cores 10
```

For example:
```
snakemake run_preprocessing --use-conda --cores 10
```

Additionally to remove the output of a given module run:
```
snakemake clean_module_name --use-conda --cores 1
```

For example:
```
snakemake clean_preprocessing --use-conda --cores 1
```

Many of later modules are dependent on the earlier modules and running a later module will run the required rules in an earlier module automatically.

# Issues with the pipeline

If you run into any issues with the pipeline and would like help please submit it to the Issues page. Please include your `config/config.yaml` file, your `pep/config.yaml` file, your `pep/samples.yaml` file, and the output from `snakemake` that includes your error. Please contact Christine Hustmyer hustmyer@wisc.edu for questions about analysis in Hustmyer et. al 2022 and mwolfe6@wisc.edu for overall questions about the ChIP-seq pipeline. 
