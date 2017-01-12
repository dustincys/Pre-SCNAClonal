# Pre-SCNAClonal

## Introduction

![Coupling GC bias](http://ww1.sinaimg.cn/large/61dccbaajw1fbcgaakjnfj20fa0b40v7.jpg "Coupling GC bias")

Tumor samples often present heterogeneity, containing not only multiple subpopulations of cancerous cells defined by distinct somatic mutations, but also the normal, non-cancerous cells.
While using the whole-genome sequencing data of tumor sample to reconstruct the subclonal composition determined by somatic copy number alternations (SCNAs), the absolute copy number and its subclonal population are both required to be estimated.
However, the raw reads in WGS data present coupling GC bias in the same SCNA genome segmentation of paired tumor and normal samples, which largely affect the absolute copy number estimation conducted by existing SCNV based subclonal inferring tools.
We provide Pre-SCNVClonal, a comprehensive package for both automatically and visually correcting coupling GC bias and baseline correction.
Pre-SCNVClonal could be strung together with the SCNV based subclonal inferring tool as a pipeline or run individually as needed.

## Requirement

- Python2.7

- python package:
 - `scipy`
 - `numpy`
 - `pickle`
 - `matplotlib`

## Manual for Visually correct GC bias

### Manually adjust the linear regression line

1. Click this button to pan axes with left mouse, zoom with right mouse

  ![Step1](http://ww4.sinaimg.cn/large/61dccbaajw1fbci9hckglj20hs0fata5.jpg "Step1")

2. Adjust the `alpha` and `area` bars to better visualize the stripes

  ![Step2](http://ww3.sinaimg.cn/large/61dccbaajw1fbci9hkr1lj20hs0fatbx.jpg "Step2")

3. Adjust the `Slope delta` and `Interception delta` bars to move linear
   regression line overlope one of the stripes

  ![Step3](http://ww4.sinaimg.cn/large/61dccbaajw1fbci9hqatbj20hs0fawhp.jpg "Step3")

4. Click `OK` button

### Manually eliminate the horizontal GC bias

1. Click this button off after the adjustment of `alpha` and `area` bars

  ![Step4](http://ww2.sinaimg.cn/large/61dccbaajw1fbci9hvl3qj20hs0faada.jpg "Step4")

2. Hold left mouse and move to pick out one stripe out, then Pre-SCNAClonal with
   redraw the linear regression line.

  ![Step5](http://ww1.sinaimg.cn/large/61dccbaajw1fbci9i06yij20hs0fajum.jpg "Step5")

3. Click `OK` button

## Usage of Pre-SCNAClonal

Generally, SCNAs based tumor composition reconstruction tools utilize the *read count*,
*B allele frequency (BAF)* information to estimate the abslute copy number and
tumor subclonal frequency.
According to the method of geting read count and BAF information, these
reconstruction tools could be grouped into two classes, 1) Get the read count
and BAF information totally from the BAM files, i.e. MixClone;
2) Get the read count and BAF position from heterozygous SNP position, then get the BAF
information from BAM files at the heterozygous SNP position. i.e. THetA;

### Example of MixClone like

#### usage

	usage: run.py MixClone [-h] [--pkl_path PKL_PATH]
			       [--max_copynumber MAX_COPYNUMBER]
			       [--subclone_num SUBCLONE_NUM]
			       [--baseline_thred_LOH BASELINE_THRED_LOH]
			       [--baseline_thred_APM BASELINE_THRED_APM]
			       [--pkl_flag PKL_FLAG] [--min_depth MIN_DEPTH]
			       [--min_base_qual MIN_BASE_QUAL]
			       [--min_map_qual MIN_MAP_QUAL]
			       [--process_num PROCESS_NUM]
			       [--gc_correction_method GC_CORRECTION_METHOD]
			       normal_bam tumor_bam reference_genome
			       input_filename_base segments_bed BICseq_bed_corrected

	positional arguments:
	  normal_bam            BAM file for normal sample.
	  tumor_bam             BAM file for tumor sample.
	  reference_genome      FASTA file for reference genome.
	  input_filename_base   Base name of the preprocessed input file to be
				created.
	  segments_bed          BED file for segments.
	  BICseq_bed_corrected  The name of corrected BICseq result file

	optional arguments:
	  -h, --help            show this help message and exit
	  --pkl_path PKL_PATH   Load the pkl path
	  --max_copynumber MAX_COPYNUMBER
				Set the maximum copy number
	  --subclone_num SUBCLONE_NUM
				Set the subclone number
	  --baseline_thred_LOH BASELINE_THRED_LOH
				The threshold of LOH sites fraction within each
				segment to define the segment is LOH, the range is
				[baseline_thred_LOH, 1]. Default is 0.16.
	  --baseline_thred_APM BASELINE_THRED_APM
				The threshold of average P and M SNP sites fraction
				within each segment to define the segment as baseline,
				the range is [baseline_thred_APM, 1]. Default is 0.2.
	  --pkl_flag PKL_FLAG   The pkl flag
	  --min_depth MIN_DEPTH
				Minimum reads depth required for both normal and tumor
				samples. Default is 20.
	  --min_base_qual MIN_BASE_QUAL
				Minimum base quality required. Default is 10.
	  --min_map_qual MIN_MAP_QUAL
				Minimum mapping quality required. Default is 10.
	  --process_num PROCESS_NUM
				Number of processes to launch for preprocessing.
				Default is 1.
	  --gc_correction_method GC_CORRECTION_METHOD
				The gc correction method, one of auto and visual

####Example

	python2.7 $PRESCNVCLONAL_PATH/run.py MixClone \
	    $normal_bam \
	    $tumor_bam \
	    $reference_genome \
	    $input_filename_base \
	    $segments_bed \
	    $BICseq_bed_corrected\
	    --pkl_path $pkl_path\
	    --max_copynumber $max_copynumber \
	    --subclone_num $subclone_num \
	    --baseline_thred_LOH $baseline_thred_LOH \
	    --baseline_thred_APM $baseline_thred_APM \
	    --min_depth $min_depth  \
	    --min_base_qual $min_base_qual \
	    --min_map_qual $min_map_qual \
	    --process_num $process_num \
	    --gc_correction_method auto \

### Example of THetA like

#### Usage

	usage: run.py THetA [-h] [--pkl_path PKL_PATH] [--pkl_flag PKL_FLAG]
			    [--max_copynumber MAX_COPYNUMBER]
			    [--subclone_num SUBCLONE_NUM]
			    [--sampleNumber SAMPLENUMBER]
			    [--gc_correction_method GC_CORRECTION_METHOD]
			    [--baseline_thred_LOH BASELINE_THRED_LOH]
			    [--baseline_thred_APM BASELINE_THRED_APM]
			    [--baseline_thred_SNPDENSITY BASELINE_THRED_SNPDENSITY]
			    [--gamma GAMMA] [--process_num PROCESS_NUM]
			    input_filename_base BICseq_bed BICseq_bed_corrected
			    tumor_SNP normal_SNP

	positional arguments:
	  input_filename_base   Base name of the preprocessed input file to be
				created.
	  BICseq_bed            BICseq result file
	  BICseq_bed_corrected  The name of corrected BICseq result file
	  tumor_SNP             tumor snp file, generated by THetA script
	  normal_SNP            normal snp file. generated by THetA script

	optional arguments:
	  -h, --help            show this help message and exit
	  --pkl_path PKL_PATH   Load the pkl path
	  --pkl_flag PKL_FLAG   The pkl flag
	  --max_copynumber MAX_COPYNUMBER
				Set the maximum copy number
	  --subclone_num SUBCLONE_NUM
				Set the subclone number
	  --sampleNumber SAMPLENUMBER
				Set the sample number for visual plot
	  --gc_correction_method GC_CORRECTION_METHOD
				The gc correction method, one of auto and visual
	  --baseline_thred_LOH BASELINE_THRED_LOH
				The threshold of LOH sites fraction within each
				segment to define the segment is LOH, the range is
				[baseline_thred_LOH, 1]. Default is 0.3.
	  --baseline_thred_APM BASELINE_THRED_APM
				The threshold of average P and M SNP sites fraction
				within each segment to define the segment as baseline,
				The range is [baseline_thred_APM, 1]. Default is 0.01.
	  --baseline_thred_SNPDENSITY BASELINE_THRED_SNPDENSITY
				The threshold of density of SNP sites in segment, The
				range is [baseline_thred_APM, 1]. Default is 0.01.
	  --gamma GAMMA         gamma parameter for determine heterozygous snp loci
	  --process_num PROCESS_NUM
				Number of processes to launch for preprocessing.
				Default is 1.


#### Example

	python2.7 ./run.py THetA \
	    $segments_bed \
	    $BICseq_bed_corrected\
	    $tumor_SNP_counts \
	    $nomal_SNP_counts \
	    --gc_correction_method [visual|auto]\
	    --baseline_thred_LOH 0.16 \
	    --baseline_thred_APM 0.6


