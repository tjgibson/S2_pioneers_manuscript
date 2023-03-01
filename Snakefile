##########################################################################################
##########################################################################################
###                                                                                    ###
###                                 Setup                                              ###
###                                                                                    ###
##########################################################################################
##########################################################################################

from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("6.0")

##### specify config file #####
configfile: 'config/config.yaml'

##########################################################################################
##########################################################################################
###                                                                                    ###
###                          ChIP-seq analysis                                         ###
###                                                                                    ###
##########################################################################################
##########################################################################################

module ChIPseq:
    snakefile: "https://github.com/tjgibson/NGS-workflow-chipseq/raw/main/workflow/Snakefile"
#     snakefile: "../../NGS_workflows/NGS-workflow-chipseq/workflow/Snakefile"
	config: config["ChIPseq"]
	prefix: "ChIPseq"
use rule * from ChIPseq as ChIPseq_*

# additional rules for processing ChIP-seq peaks
# rule regions_to_exclude:
# 	input:
# 		"ChIPseq/results/peaks/filtered/S2-WT_aZld_IP.narrowPeak",
# 		"ChIPseq/results/peaks/filtered/S2-WT_aGrh_IP.narrowPeak",
# 		"ChIPseq/results/peaks/filtered/S2-WT_aTwi.narrowPeak",
# 		
# 	output:
# 		"resources/zld_regions_to_exclude.bed",
# 		"resources/grh_regions_to_exclude.bed",
# 		"resources/twi_regions_to_exclude.bed",
# 	params:
# 		"config/zld_artifacts.bed",
# 		"config/grh_artifacts.bed",
# 		"config/twi_artifacts.bed",
# 	shell:
# 		"""
# 		cat {input[0]} {params[0]} | cut -f 1,2,3 > {output[0]}
# 		cat {input[1]} {params[1]} | cut -f 1,2,3 > {output[1]}
# 		cat {input[2]} {params[2]} | cut -f 1,2,3 > {output[2]}
# 		"""
rule regions_to_exclude:
	input:
		"ChIPseq/results/peaks/filtered/S2-WT_aZld_IP.narrowPeak",
		"ChIPseq/results/peaks/filtered/S2-WT_aGrh_IP.narrowPeak",
		"ChIPseq/results/peaks/filtered/S2-WT_aTwi.narrowPeak",
	output:
		"resources/zld_regions_to_exclude.bed",
		"resources/grh_regions_to_exclude.bed",
		"resources/twi_regions_to_exclude.bed",
	params:
		"config/zld_artifacts.bed",
		"config/grh_artifacts.bed",
		"config/twi_artifacts.bed",
	shell:
		"""
		cat {params[0]} | cut -f 1,2,3 > {output[0]}
		cat {params[1]} | cut -f 1,2,3 > {output[1]}
		cat {params[2]} | cut -f 1,2,3 > {output[2]}
		"""


rule process_ChIP_peaks:
	input:
		"ChIPseq/results/peaks/filtered/{sample}.narrowPeak",
		rules.regions_to_exclude.output
	output:
		"ChIPseq/results/peaks/final/{sample}.narrowPeak",
		"ChIPseq/results/peaks/final/{sample}.bed"
	conda:
		"workflow/envs/r_process_bed.yaml"
	params:
		extend_width=201,
		r_source= "workflow/scripts/utils.R",
	script:
		"workflow/scripts/process_peaks.R"

rule annotate_ChIP_peaks:
	input:
		"ChIPseq/results/peaks/final/{sample}.bed"
	output:
		"ChIPseq/results/peaks/final/{sample}_peak_annotations.tsv"
# 	conda:
# 		"workflow/envs/r_ChIPseeker.yaml"
	params:
		promoter_upstream= -500,
		promoter_downstream= 100
	script:
		"workflow/scripts/annotate_peaks_ChIPseeker.R"
		
# process previously published ChIP-seq data for chromatin proteins
module published_ChIPseq:
    snakefile: "https://github.com/tjgibson/NGS-workflow-chipseq/raw/main/workflow/Snakefile"
# 	snakefile: "../../NGS_workflows/NGS-workflow-chipseq/workflow/Snakefile"    	
	config: config["published_ChIPseq"]
	prefix: "published_ChIPseq"
use rule * from published_ChIPseq as published_ChIPseq_*

# additional rules for analysis of published data 
rule get_embryo_peak_seqs:
	input:
		"published_ChIPseq/results/peaks/merged/narrow/{sample}_peaks.narrowPeak"
	output:
		temp("published_ChIPseq/results/peaks/peak_sequences/{sample}.fasta")
	params:
		r_source= "workflow/scripts/utils.R",
	script:
		"workflow/scripts/get_embryo_peak_sequences.R"

rule get_embryo_motifs:
	input:
		"published_ChIPseq/results/peaks/peak_sequences/{sample}.fasta"
	output:
		multiext("published_ChIPseq/results/motifs/{sample}/meme", ".html", ".txt", ".xml")
	conda:
		"workflow/envs/meme.yaml"
	params:
		mode="zoops",
		n_motifs=10,
	shell:
		 "meme {input} -oc published_ChIPseq/results/motifs/{wildcards.sample} -revcomp -mod {params[0]} -nmotifs {params[1]} -dna"


def get_meme_motif_index(wildcards):
	motif_index = {
	"embryo-nc14_aZld": 3,
	"embryo-15-16H_aGrh": 3,
	"embryo-1-3H_aTwi": 5
	}
	
	return motif_index[wildcards.sample]
	


rule get_embryo_PWMs:
	input:
		meme_file="published_ChIPseq/results/motifs/{sample}/meme.txt",
	output:
		"published_ChIPseq/results/motifs/{sample}/{sample}_PWM.rds"
	params:
		meme_motif_number= get_meme_motif_index,
	script:
		"workflow/scripts/get_motif_PWM.R"
		

##########################################################################################
##########################################################################################
###                                                                                    ###
###                           RNA-seq analysis                                         ###
###                                                                                    ###
##########################################################################################
##########################################################################################

module RNAseq:
    snakefile:
        "https://github.com/tjgibson/NGS-workflow-RNAseq/raw/main/workflow/Snakefile"
#         "../NGS_workflows/NGS-workflow-RNAseq/workflow/Snakefile"
		config: config["RNAseq"]
		prefix: "RNAseq"
use rule * from RNAseq as RNAseq_*

# additional rules for processing RNA-seq data
# - Call differential expression, incorporating controls as covariates
# - Exclude any problem regions as needed (Zld, Grh, Twi, MtnA etc.)
# - Annotate genes with ChIP data

def get_ChIP_peaks_annotated(wildcards):
	motif_index = {
	"S2-Zld_RNAseq": "ChIPseq/results/peaks/final/S2-Zld_aZld_IP_peak_annotations.tsv",
	"S2-Grh_RNAseq": "ChIPseq/results/peaks/final/S2-Grh_aGrh_IP_peak_annotations.tsv",
	"S2-Twi_RNAseq": "ChIPseq/results/peaks/final/S2-Twi_aTwi_IP_peak_annotations.tsv",
	}
	
	return motif_index[wildcards.experiment]
	

rule annotate_RNAseq_results:
	input:
		"RNAseq/results/DEseq2/{experiment}_{contrast}_results.tsv",
		get_ChIP_peaks_annotated,
	output:
		"RNAseq/results/DEseq2/{experiment}_{contrast}_results_annotated.tsv"
# 	conda:
# 		"workflow/envs/r_ChIPseeker.yaml"
	script:
		"workflow/scripts/annotate_RNAseq_results.R"


##########################################################################################
##########################################################################################
###                                                                                    ###
###                          ATAC-seq analysis                                         ###
###                                                                                    ###
##########################################################################################
##########################################################################################

module ATACseq:
	snakefile: "https://github.com/tjgibson/NGS-workflow-ATACseq/raw/main/workflow/Snakefile"
# 	snakefile: "../../NGS_workflows/NGS-workflow-ATACseq/workflow/Snakefile"
	config: config["ATACseq"]
	prefix: "ATACseq"
use rule * from ATACseq as ATACseq_*

# additional rules for ATAC-seq data
# - Exclude any problem regions as needed (Zld, Grh, Twi, MtnA etc.)
# - Annotate peaks with ChIP data

rule filter_ATAC_results:
	input:
		"ATACseq/results/DEseq2/{experiment}_{contrast}_results.tsv"
	output:
		"ATACseq/results/DEseq2_results_filtered/{experiment}_{contrast}_results.tsv",
	conda:
		"workflow/envs/r_process_bed.yaml"
	script:
		"workflow/scripts/filter_ATAC_results.R"

rule annotate_ATAC_peaks:
	input:
		"ATACseq/results/peaks/final/{sample}.bed"
	output:
		"ATACseq/results/peaks/final/{sample}_peak_annotations.tsv"
# 	conda:
# 		"workflow/envs/r_ChIPseeker.yaml"
	params:
		promoter_upstream= -500,
		promoter_downstream= 100
	script:
		"workflow/scripts/annotate_peaks_ChIPseeker.R"

##########################################################################################
##########################################################################################
###                                                                                    ###
###                           histone CUT&RUN analysis                                 ###
###                                                                                    ###
##########################################################################################
##########################################################################################

# process CUT&RUN data
module histone_CUTandRUN:
	snakefile: "https://github.com/tjgibson/NGS-workflow-CUTandRUN/raw/master/workflow/Snakefile"
# 	snakefile: "../NGS_workflows/NGS-workflow-CUTandRUN/workflow/Snakefile"
	config: config["histone_CUTandRUN"]
	prefix: "histone_CUTandRUN"
use rule * from histone_CUTandRUN as histone_CUTandRUN_*

##########################################################################################
##########################################################################################
###                                                                                    ###
###                             TF CUT&RUN analysis                                    ###
###                                                                                    ###
##########################################################################################
##########################################################################################

# process CUT&RUN data
module TF_CUTandRUN:
	snakefile: "https://github.com/tjgibson/NGS-workflow-CUTandRUN/raw/master/workflow/Snakefile"
# 	snakefile: "../NGS_workflows/NGS-workflow-CUTandRUN/workflow/Snakefile"
	config: config["TF_CUTandRUN"]
	prefix: "TF_CUTandRUN"
use rule * from TF_CUTandRUN as TF_CUTandRUN_*



##########################################################################################
##########################################################################################
###                                                                                    ###
###                          integrated analysis                                       ###
###                                                                                    ###
##########################################################################################
##########################################################################################

# Annotate ChIP classes
# - Use output of above modules to annotate ChIP peaks as class I, II, or III
# - Annotate with additional information 
#     - Diff expression of proximal gene, FC, padj
#     - Diff accessibility, FC, padj
#     - Promoter/distal
#     - etc.

def get_ChIP_class_input(wildcards):
	return {
	"ChIP_peaks":config["annotate_ChIP_classes"][wildcards.factor]["ChIP_peaks"],
	"WT_ATAC_peaks":config["annotate_ChIP_classes"][wildcards.factor]["WT_ATAC_peaks"],
	"ATAC_results":config["annotate_ChIP_classes"][wildcards.factor]["ATAC_results"],
	"RNAseq_results":config["annotate_ChIP_classes"][wildcards.factor]["RNAseq_results"],
	"ChIP_feature_annotation":config["annotate_ChIP_classes"][wildcards.factor]["ChIP_feature_annotation"],
	"motif_instances":config["annotate_ChIP_classes"][wildcards.factor]["motif_instances"]
	}
rule annotate_ChIP_classes:
	input:
		unpack(get_ChIP_class_input),
		gtf_genome_annotation=rules.RNAseq_get_genome_annotation.output
	output:
		"results/ChIP_peak_classes/{factor}_ChIP_classes.tsv"
	params:
		r_source= "workflow/scripts/utils.R",
	script:
		"workflow/scripts/annotate_ChIP_classes.R"

rule split_ChIP_classes:
	input:
		"results/ChIP_peak_classes/{factor}_ChIP_classes.tsv"
	output:
		"results/ChIP_peak_classes/{factor}_class_I.bed",
		"results/ChIP_peak_classes/{factor}_class_II.bed",
		"results/ChIP_peak_classes/{factor}_class_III.bed",
		"results/ChIP_peak_classes/{factor}_class_II_and_III.bed",
		"results/ChIP_peak_classes/{factor}_class_I.fasta",
		"results/ChIP_peak_classes/{factor}_class_II.fasta",
		"results/ChIP_peak_classes/{factor}_class_III.fasta",
		"results/ChIP_peak_classes/{factor}_class_II_and_III.fasta"
	script:
		"workflow/scripts/split_ChIP_classes.R"

rule define_nonspecific_sites:
	input:
		"results/ChIP_peak_classes/zld_class_I.bed",
		"results/ChIP_peak_classes/grh_class_I.bed",
		"results/ChIP_peak_classes/twi_class_I.bed",
	output:
		"results/ChIP_peak_classes/nonspecific_sites/nonspecific_sites.fasta",
		"results/ChIP_peak_classes/nonspecific_sites/nonspecific_sites.bed",
	params:
		r_source= "workflow/scripts/utils.R",
	script:
		"workflow/scripts/get_nonspecific_peaks.R"

rule filter_ChIP_classes:
	input:
		"results/ChIP_peak_classes/{factor}_ChIP_classes.tsv",
		"results/ChIP_peak_classes/nonspecific_sites/nonspecific_sites.bed",
	output:
		"results/ChIP_peak_classes_filtered/{factor}_ChIP_classes.tsv",
		"results/ChIP_peak_classes_filtered/{factor}_ChIP_classes.bed",
	params:
		r_source= "workflow/scripts/utils.R",
	conda:
		"workflow/envs/r_process_bed.yaml"
	script:
		"workflow/scripts/filter_ChIP_classes.R"

rule split_ChIP_classes_filtered:
	input:
		"results/ChIP_peak_classes_filtered/{factor}_ChIP_classes.tsv"
	output:
		"results/ChIP_peak_classes_filtered/{factor}_class_I.bed",
		"results/ChIP_peak_classes_filtered/{factor}_class_II.bed",
		"results/ChIP_peak_classes_filtered/{factor}_class_III.bed",
		"results/ChIP_peak_classes_filtered/{factor}_class_II_and_III.bed",
		"results/ChIP_peak_classes_filtered/{factor}_class_I.fasta",
		"results/ChIP_peak_classes_filtered/{factor}_class_II.fasta",
		"results/ChIP_peak_classes_filtered/{factor}_class_III.fasta",
		"results/ChIP_peak_classes_filtered/{factor}_class_II_and_III.fasta"
	script:
		"workflow/scripts/split_ChIP_classes.R"

def get_motif_pwm(wildcards):
	motif_pwm = {
	"zld":"published_ChIPseq/results/motifs/embryo-nc14_aZld/embryo-nc14_aZld_PWM.rds",
	"grh":"published_ChIPseq/results/motifs/embryo-15-16H_aGrh/embryo-15-16H_aGrh_PWM.rds",
	"twi":"published_ChIPseq/results/motifs/embryo-1-3H_aTwi/embryo-1-3H_aTwi_PWM.rds"
	}
	
	return motif_pwm[wildcards.factor]
	

rule get_motif_instances:
	input:
		PWM=get_motif_pwm
	output:
		bed="results/motif_instances/{factor}_motifs.bed",
		tsv="results/motif_instances/{factor}_motifs.tsv",
		bw__motif_score="results/motif_instances/{factor}_motifs_score.bw",
		bw_motif_pvalue="results/motif_instances/{factor}_motifs_pvalue.bw",
	params:
		threshold=1e-3
	script:
		"workflow/scripts/get_PWM_instances.R"

rule annotate_motif_classes:
	input:
		all_peaks = lambda wildcards: config["annotate_motif_classes"][wildcards.factor]["ChIP_peaks"],
		class_I_peaks = lambda wildcards: config["annotate_motif_classes"][wildcards.factor]["class_I_peaks"],
		class_II_peaks = lambda wildcards: config["annotate_motif_classes"][wildcards.factor]["class_II_peaks"],
		class_III_peaks = lambda wildcards: config["annotate_motif_classes"][wildcards.factor]["class_III_peaks"],
		motifs = lambda wildcards: config["annotate_motif_classes"][wildcards.factor]["motif_instances"],
		other_tissue_peaks = lambda wildcards: config["annotate_motif_classes"][wildcards.factor]["other_tissue"],
		H3K27me3_peaks = "published_ChIPseq/results/peaks/individual/broad/GSE151983_S2_aH3K27me3_IP_peaks.broadPeak",
		H3K9me3_peaks = "published_ChIPseq/results/peaks/filtered/GSE160855_aH3K9me3.broadPeak",
		keep_chroms = "config/keep_chroms.txt",
		nonspecific_sites = "results/ChIP_peak_classes/nonspecific_sites/nonspecific_sites.bed",
	params:
		bichrom_window_size=config["bichrom_params"]["window_size"],
	output:
		"results/ChIP_tissue_classes/{factor}_motif_classes.tsv",
	script:
		"workflow/scripts/get_motif_classes.R"

rule annotate_tissue_classes:
	input:
		all_peaks = lambda wildcards: config["annotate_motif_classes"][wildcards.factor]["ChIP_peaks"],
		class_I_peaks = lambda wildcards: config["annotate_motif_classes"][wildcards.factor]["class_I_peaks"],
		class_II_peaks = lambda wildcards: config["annotate_motif_classes"][wildcards.factor]["class_II_peaks"],
		class_III_peaks = lambda wildcards: config["annotate_motif_classes"][wildcards.factor]["class_III_peaks"],
		motifs = lambda wildcards: config["annotate_motif_classes"][wildcards.factor]["motif_instances"],
		other_tissue_peaks = lambda wildcards: config["annotate_motif_classes"][wildcards.factor]["other_tissue"],
		H3K27me3_peaks = "published_ChIPseq/results/peaks/individual/broad/GSE151983_S2_aH3K27me3_IP_peaks.broadPeak",
		H3K9me3_peaks = "published_ChIPseq/results/peaks/filtered/GSE160855_aH3K9me3.broadPeak",
		keep_chroms = "config/keep_chroms.txt",
		nonspecific_sites = "results/ChIP_peak_classes/nonspecific_sites/nonspecific_sites.bed",
	params:
		bichrom_window_size=config["bichrom_params"]["window_size"],
	output:
		"results/ChIP_tissue_classes/{factor}_tissue_classes.tsv",
	script:
		"workflow/scripts/get_tissue_classes.R"


rule annotate_titration_classes:
	input:
		"ChIPseq/results/peaks/filtered/S2-Zld-500uM_aZld_IP.narrowPeak",
		"ATACseq/results/peaks/merged_by_sample/S2-WT_1000uM_summits.bed",
		"ATACseq/results/DEseq2/S2-Zld-titration_ATACseq_S2-Zld-500uM-vs-0uM_results.tsv",
		"ChIPseq/results/peaks/filtered/S2-Zld-1000uM_aZld_IP.narrowPeak",
		"ATACseq/results/peaks/merged_by_sample/S2-WT_1000uM_summits.bed",
		"ATACseq/results/DEseq2/S2-Zld-titration_ATACseq_S2-Zld-1000uM-vs-0uM_results.tsv",
		"ChIPseq/results/peaks/filtered/S2-Zld-1500uM_aZld_IP.narrowPeak",
		"ATACseq/results/peaks/merged_by_sample/S2-WT_1000uM_summits.bed",
		"ATACseq/results/DEseq2/S2-Zld-titration_ATACseq_S2-Zld-1500uM-vs-0uM_results.tsv",
		"ChIPseq/results/peaks/filtered/S2-Grh-25uM_aGrh_IP.narrowPeak",
		"ATACseq/results/peaks/merged_by_sample/FL_ATAC_S2-WT_100uM_summits.bed",
		"ATACseq/results/DEseq2/S2-Grh-titration_ATACseq_S2-Grh-25uM-vs-0uM_results.tsv",
		"ChIPseq/results/peaks/filtered/S2-Grh-100uM_aGrh_IP.narrowPeak",
		"ATACseq/results/peaks/merged_by_sample/FL_ATAC_S2-WT_100uM_summits.bed",
		"ATACseq/results/DEseq2/S2-Grh-titration_ATACseq_S2-Grh-100uM-vs-0uM_results.tsv",
		"ChIPseq/results/peaks/filtered/S2-Grh-400uM_aGrh_IP.narrowPeak",
		"ATACseq/results/peaks/merged_by_sample/FL_ATAC_S2-WT_100uM_summits.bed",
		"ATACseq/results/DEseq2/S2-Grh-titration_ATACseq_S2-Grh-400uM-vs-0uM_results.tsv",
		"ChIPseq/results/peaks/filtered/S2-HA-Twi-10uM_aHA_IP.narrowPeak",
		"ATACseq/results/peaks/merged_by_sample/Twi_ATAC_S2-WT_40uM_summits.bed",
		"ATACseq/results/DEseq2/S2-HA-Twi-titration_ATACseq_S2-HA-Twi-10uM-vs-0uM_results.tsv",
		"ChIPseq/results/peaks/filtered/S2-HA-Twi-40uM_aHA_IP.narrowPeak",
		"ATACseq/results/peaks/merged_by_sample/Twi_ATAC_S2-WT_40uM_summits.bed",
		"ATACseq/results/DEseq2/S2-HA-Twi-titration_ATACseq_S2-HA-Twi-40uM-vs-0uM_results.tsv",
		"ChIPseq/results/peaks/filtered/S2-HA-Twi-160uM_aHA_IP.narrowPeak",
		"ATACseq/results/peaks/merged_by_sample/Twi_ATAC_S2-WT_40uM_summits.bed",
		"ATACseq/results/DEseq2/S2-HA-Twi-titration_ATACseq_S2-HA-Twi-160uM-vs-0uM_results.tsv",
	output:
		"results/ChIP_titration_classes/zld_titration_classes.tsv",
		"results/ChIP_titration_classes/grh_titration_classes.tsv",
		"results/ChIP_titration_classes/twi_titration_classes.tsv",
	script:
		"workflow/scripts/define_titration_classes.R"

rule generate_fig5_RPKM_tables:
	input:
		"results/ChIP_tissue_classes/zld_tissue_classes.tsv",
		"results/ChIP_tissue_classes/grh_tissue_classes.tsv",
		"results/ChIP_tissue_classes/twi_tissue_classes.tsv",
		rules.ChIPseq_all.input,
	output:
		"results/ChIP_tissue_classes/zld_classes_titration_ChIP_rpkm.tsv",
		"results/ChIP_tissue_classes/zld_classes_titration_ATAC_rpkm.tsv",
		"results/ChIP_tissue_classes/grh_classes_titration_ChIP_rpkm.tsv",
		"results/ChIP_tissue_classes/grh_classes_titration_ATAC_rpkm.tsv",
		"results/ChIP_tissue_classes/twi_classes_titration_ChIP_rpkm.tsv",
		"results/ChIP_tissue_classes/twi_classes_titration_ATAC_rpkm.tsv",
	script:
		"workflow/scripts/generate_fig5_RPKM_tables.R"

##########################################################################################
##########################################################################################
###                                                                                    ###
###                        train deep learning models                                  ###
###                                                                                    ###
##########################################################################################
##########################################################################################
# run bichrom
include: "bichrom/workflow/rules/bichrom.smk"

# run BPnet
include: "BPnet/workflow/rules/BPnet.smk"



##########################################################################################
##########################################################################################
###                                                                                    ###
###                        generate figures for MS                                     ###
###                                                                                    ###
##########################################################################################
##########################################################################################
rule figure_1:
	input:
		Zld_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw",
		Zld_WT_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/S2-WT_1000uM_small.bw",
		Zld_Zld_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Zld_1000uM_small.bw",
		Zld_WT_RNAseq_bw = "RNAseq/results/bigwigs/zscore_normalized/merged/S2-WT_1000uM.bw",
		Zld_Zld_RNAseq_bw = "RNAseq/results/bigwigs/zscore_normalized/merged/S2-Zld_1000mM.bw",
		Grh_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw",
		Grh_WT_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/FL_ATAC_S2-WT_100uM_small.bw",
		Grh_Grh_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Grh_100uM_small.bw",
		Grh_WT_RNAseq_bw = "RNAseq/results/bigwigs/zscore_normalized/merged/S2-WT_100uM.bw",
		Grh_Grh_RNAseq_bw = "RNAseq/results/bigwigs/zscore_normalized/merged/S2-Grh_100mM.bw",
		zld_ChIP_classes = "results/ChIP_peak_classes/zld_ChIP_classes.tsv",
		grh_ChIP_classes = "results/ChIP_peak_classes/grh_ChIP_classes.tsv"
	output:
		"manuscript/figures/fig1.pdf"
	script:
		"workflow/scripts/fig1.R"

rule figure_2:
	input:
		Twi_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw",
		Twi_WT_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/Twi_ATAC_S2-WT_40uM_small.bw",
		Twi_Twi_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Twi_40uM_small.bw",
		Zld_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw",
		Grh_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw",
		twi_ChIP_classes = "results/ChIP_peak_classes/twi_ChIP_classes.tsv",
		zld_ChIP_classes = "results/ChIP_peak_classes/zld_ChIP_classes.tsv",
		grh_ChIP_classes = "results/ChIP_peak_classes/grh_ChIP_classes.tsv",
		twi_RNAseq_results = "RNAseq/results/DEseq2/S2-Twi_RNAseq_S2-Twi-vs-S2-WT_results.tsv",
		zld_RNAseq_results = "RNAseq/results/DEseq2/S2-Zld_RNAseq_S2-Zld-vs-S2-WT_results.tsv",
		grh_RNAseq_results = "RNAseq/results/DEseq2/S2-Grh_RNAseq_S2-Grh-vs-S2-WT_results.tsv"
	output:
		"manuscript/figures/fig2.pdf"
	script:
		"workflow/scripts/fig2.R"

rule figure_3:
	input:
  		zld_class_I_bed_fn = "results/ChIP_peak_classes/zld_class_I.bed",
  		grh_class_I_bed_fn = "results/ChIP_peak_classes/grh_class_I.bed",
  		twi_class_I_bed_fn = "results/ChIP_peak_classes/twi_class_I.bed",

		ns_sites_bed_fn = "results/ChIP_peak_classes/nonspecific_sites/nonspecific_sites.bed",

		S2_Zld_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw",
		S2_Grh_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw",
		S2_Twi_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw",
		H3K27ac_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K27ac.bw",
		Nej_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE64464_aNej.bw",
		H3K4me1_bw =  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me1.bw",
		H3K4me3_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me3.bw",
		H2AV_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_H2Av_IP.bw",

		zld_ChIP_classes_fn = "results/ChIP_peak_classes/zld_ChIP_classes.tsv",
		grh_ChIP_classes_fn = "results/ChIP_peak_classes/grh_ChIP_classes.tsv",
		twi_ChIP_classes_fn = "results/ChIP_peak_classes/twi_ChIP_classes.tsv",

		zld_motifs_fn = "results/motif_instances/zld_motifs.tsv",
		grh_motifs_fn = "results/motif_instances/grh_motifs.tsv",
		twi_motifs_fn = "results/motif_instances/twi_motifs.tsv",

		zld_WT_atac_fn = "ATACseq/results/peaks/merged_by_sample/S2-WT_1000uM_summits.bed",
		grh_WT_atac_fn = "ATACseq/results/peaks/merged_by_sample/FL_ATAC_S2-WT_100uM_summits.bed",
		twi_WT_atac_fn = "ATACseq/results/peaks/merged_by_sample/Twi_ATAC_S2-WT_40uM_summits.bed",
	output:
		"manuscript/figures/fig3.pdf"
	script:
		"workflow/scripts/fig3.R"

rule figure_4:
	input:
		zld_tissue_occupancy = "results/ChIP_tissue_classes/zld_tissue_classes.tsv",
		grh_tissue_occupancy = "results/ChIP_tissue_classes/grh_tissue_classes.tsv",
		twi_tissue_occupancy = "results/ChIP_tissue_classes/twi_tissue_classes.tsv",
		S2_ZLD_ChIP_bw =   "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw",
		embryo_Zld_ChIP_bw =  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-nc14_aZld.bw",
		brain_Zld_ChIP_bw =  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/brain_aZld.bw",
		S2_H3K27me3_bw =  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE151983_S2_aH3K27me3_IP.bw",
		S2_H3K9me3_bw =  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE160855_aH3K9me3.bw",
		S2_Grh_ChIP_bw =   "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw",
		embryo_Grh_ChIP_bw =  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-15-16H_aGrh.bw",
		wing_disc_Grh_ChIP_bw =  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/wing-disc_aGrh.bw",
		S2_Twi_ChIP_bw =  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw",
		embryo_Twi_ChIP_bw =  "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-1-3H_aTwi.bw",
		S2_Zld_ChIP_DMSO_bw =  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_DMSO_aZld.bw",
		S2_Zld_ChIP_taz_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_taz_aZld.bw",
		S2_Zld_H3K27me3_DMSO_bw = "histone_CUTandRUN/results/bigwigs/spikeIn_normalized/merged/S2-Zld_DMSO_aH3K27me3_large.bw",
		S2_Zld_H3K27me3_taz_bw = "histone_CUTandRUN/results/bigwigs/spikeIn_normalized/merged/S2-Zld_Taz_aH3K27me3_large.bw",
		S2_Grh_ChIP_DMSO_bw =  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_DMSO_aGrh.bw",
		S2_Grh_ChIP_taz_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_taz_aGrh.bw",
		S2_Grh_H3K27me3_DMSO_bw = "histone_CUTandRUN/results/bigwigs/spikeIn_normalized/merged/S2-Grh_DMSO_aH3K27me3_large.bw",
		S2_Grh_H3K27me3_taz_bw = "histone_CUTandRUN/results/bigwigs/spikeIn_normalized/merged/S2-Grh_Taz_aH3K27me3_large.bw",
		S2_Twi_ChIP_DMSO_bw =  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-HA-Twi_DMSO_aHA.bw",
		S2_Twi_ChIP_taz_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-HA-Twi_taz_aHA.bw",
		S2_Twi_H3K27me3_DMSO_bw = "histone_CUTandRUN/results/bigwigs/spikeIn_normalized/merged/S2-HA-Twi_DMSO_aH3K27me3_large.bw",
		S2_Twi_H3K27me3_taz_bw = "histone_CUTandRUN/results/bigwigs/spikeIn_normalized/merged/S2-HA-Twi_Taz_aH3K27me3_large.bw",
	output:
		"manuscript/figures/fig4.pdf"
	script:
		"workflow/scripts/fig4.R"

rule figure_5:
	input:
		zld_tissue_occupancy = "results/ChIP_tissue_classes/zld_tissue_classes.tsv",
		grh_tissue_occupancy = "results/ChIP_tissue_classes/grh_tissue_classes.tsv",
		twi_tissue_occupancy = "results/ChIP_tissue_classes/twi_tissue_classes.tsv",
		zld_titration_classes_fn = "results/ChIP_titration_classes/zld_titration_classes.tsv",
		grh_titration_classes_fn = "results/ChIP_titration_classes/grh_titration_classes.tsv",
		twi_titration_classes_fn = "results/ChIP_titration_classes/twi_titration_classes.tsv",
		zld_titration_ChIP_rpkm_fn = "results/ChIP_tissue_classes/zld_classes_titration_ChIP_rpkm.tsv",
		zld_titration_ATAC_rpkm_fn = "results/ChIP_tissue_classes/zld_classes_titration_ATAC_rpkm.tsv",
		grh_titration_ChIP_rpkm_fn = "results/ChIP_tissue_classes/grh_classes_titration_ChIP_rpkm.tsv",
		grh_titration_ATAC_rpkm_fn = "results/ChIP_tissue_classes/grh_classes_titration_ATAC_rpkm.tsv",
		twi_titration_ChIP_rpkm_fn = "results/ChIP_tissue_classes/twi_classes_titration_ChIP_rpkm.tsv",
		twi_titration_ATAC_rpkm_fn = "results/ChIP_tissue_classes/twi_classes_titration_ATAC_rpkm.tsv",

	output:
		"manuscript/figures/fig5.pdf"
	script:
		"workflow/scripts/fig5.R"

rule figure_6:
	input:
  		Zld_FL_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw",
		Zld_DBD_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld-DBD_aZld_IP.bw",
		Zld_WT_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/S2-WT_1000uM_small.bw",
		Zld_Zld_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Zld_1000uM_small.bw",

		Grh_FL_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw",
		Grh_DBD_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh-DBD_aGrh_IP.bw",
		Grh_WT_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/FL_ATAC_S2-WT_100uM_small.bw",
		Grh_Grh_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Grh_100uM_small.bw",

		zld_ChIP_classes = "results/ChIP_peak_classes/zld_ChIP_classes.tsv",
		grh_ChIP_classes = "results/ChIP_peak_classes/grh_ChIP_classes.tsv",

		zld_FL_atac_results = "ATACseq/results/DEseq2_results_filtered/S2-Zld_ATACseq_S2-Zld-FL-vs-S2-WT_results.tsv",
		zld_DBD_atac_results = "ATACseq/results/DEseq2_results_filtered/S2-Zld-DBD_ATACseq_S2-Zld-DBD-vs-S2-WT_results.tsv",
		grh_FL_atac_results = "ATACseq/results/DEseq2_results_filtered/S2-Grh_ATACseq_S2-Grh-FL-vs-S2-WT_results.tsv",
		grh_DBD_atac_results = "ATACseq/results/DEseq2_results_filtered/S2-Grh-DBD_ATACseq_S2-Grh-DBD-vs-S2-WT_results.tsv",

	output:
		"manuscript/figures/fig6.pdf"
	script:
		"workflow/scripts/fig6.R"

rule extended_data_fig_1:
	input:
  		RPKM_table_fn = "RNAseq/results/count_tables/S2-Grh_RNAseq_RPKM.tsv",
		zld_blot_image = "data/immunoblot_raw_images/2018-10-17_Zld_induction/2018-1018-144408_pub.tif",
		grh_blot_image = "data/immunoblot_raw_images/2018-11-08_Grh_induction/2018-1108-132834_pub.tif",
		S2_Zld_aZld_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw",
		S2_WT_aZld_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-WT_aZld_IP.bw",
		S2_Zld_IgG_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2_Zld_aIgG_IP.bw",
		zld_ChIP_classes_fn = "results/ChIP_peak_classes/zld_ChIP_classes.tsv",
		S2_Grh_aGrh_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw",
		S2_WT_aGrh_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-WT_aGrh_IP.bw",
		S2_Grh_IgG_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2_Grh_aIgG_IP.bw",
		grh_ChIP_classes_fn = "results/ChIP_peak_classes/grh_ChIP_classes.tsv",
		S2_Twi_aTwi_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw",
		S2_WT_aTwi_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-WT_aTwi.bw",
		S2_Twi_IgG_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2_Zld_aIgG_IP.bw",
		twi_ChIP_classes_fn = "results/ChIP_peak_classes/twi_ChIP_classes.tsv",
	output:
		"manuscript/figures/extended_data_fig1.pdf"
	script:
		"workflow/scripts/extended_data_fig1.R"

rule extended_data_fig_2:
	input:
		S2_Zld_CR_0H_bw = "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_0H_small.bw",
		S2_Zld_CR_4H_bw = "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_4H_small.bw",
		S2_Zld_CR_12H_bw = "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_12H_small.bw",
		S2_Zld_CR_24H_bw = "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_24H_small.bw",
		S2_Zld_CR_48H_bw = "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_48H_small.bw",
		S2_Grh_CR_0H_bw = "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_0H_small.bw",
		S2_Grh_CR_4H_bw = "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_4H_small.bw",
		S2_Grh_CR_12H_bw = "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_12H_small.bw",
		S2_Grh_CR_24H_bw = "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_24H_small.bw",
		S2_Grh_CR_48H_bw = "TF_CUTandRUN/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_48H_small.bw",
		zld_ChIP_classes_fn = "results/ChIP_peak_classes/zld_ChIP_classes.tsv",
		grh_ChIP_classes_fn = "results/ChIP_peak_classes/grh_ChIP_classes.tsv",
		zld_blot_image = "data/immunoblot_raw_images/2018-11-27_timecourse/zld_timecourse.tif",
		grh_blot_image = "data/immunoblot_raw_images/2018-11-27_timecourse/grh_timecourse.tif",
	output:
		"manuscript/figures/extended_data_fig2.pdf"
	script:
		"workflow/scripts/extended_data_fig2.R"


rule extended_data_fig_3:
	input:
  		zld_RNAseq_results_fn = "RNAseq/results/DEseq2/S2-Zld_RNAseq_S2-Zld-vs-S2-WT_results_annotated.tsv",
		grh_RNAseq_results_fn  = "RNAseq/results/DEseq2/S2-Grh_RNAseq_S2-Grh-vs-S2-WT_results_annotated.tsv",
		zld_ATAC_results_fn = "ATACseq/results/DEseq2_results_filtered/S2-Zld_ATACseq_S2-Zld-FL-vs-S2-WT_results.tsv",
		grh_ATAC_results_fn = "ATACseq/results/DEseq2_results_filtered/S2-Grh_ATACseq_S2-Grh-FL-vs-S2-WT_results.tsv",
		zld_ChIP_peaks_fn = "ChIPseq/results/peaks/final/S2-Zld_aZld_IP.narrowPeak",
		grh_ChIP_peaks_fn = "ChIPseq/results/peaks/final/S2-Grh_aGrh_IP.narrowPeak",
	output:
		"manuscript/figures/extended_data_fig3.pdf"
	script:
		"workflow/scripts/extended_data_fig3.R"

rule extended_data_fig_4:
	input:
  		twi_RNAseq_results_fn = "RNAseq/results/DEseq2/S2-Twi_RNAseq_S2-Twi-vs-S2-WT_results_annotated.tsv",
		twi_ATAC_results_fn = "ATACseq/results/DEseq2_results_filtered/S2-Twi_ATACseq_S2-Twi-vs-S2-WT-40uM_results.tsv",
		twi_ChIP_peaks_fn = "ChIPseq/results/peaks/final/S2-Twi_aTwi_IP.narrowPeak",
		twi_blot_image = "data/immunoblot_raw_images/2021-9-13_Twi-vs-embryo/Twi-vs-embryos.tif",
	output:
		"manuscript/figures/extended_data_fig4.pdf"
	script:
		"workflow/scripts/extended_data_fig4.R"
		
rule extended_data_fig_5:
	input:
		H3K27ac_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K27ac.bw",
		Nej_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE64464_aNej.bw",
		H3K4me1_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me1.bw",
		H3K4me3_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE85191_aH3K4me3.bw",
		H4K16ac_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE37865_aH4K16ac.bw",
		H2AV_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_H2Av_IP.bw",
		Rpb3_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_Rpb3_IP.bw",
		PolII_pS2_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_PolII_phosphoSer2_IP.bw",
		H3K36me3_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_H3K36me3_IP.bw",
		SSRP1_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_SSRP1_IP.bw",
		Spt16_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_Spt16_IP.bw",
		H3_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE129236_H3_IP.bw",
		H1_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE127227_aH1.bw",
		Pho_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE84502_aPho_IP.bw",
		Ez_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_Ez_IP.bw",
		Pc_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE24521_Pc_IP.bw",
		Psc_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE24521_Psc_IP.bw",
		Ph_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE24521_Ph_IP.bw",
		dRing_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_dRing_IP.bw",
		H3K27me3_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE151983_S2_aH3K27me3_IP.bw",
		H3K9me3_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE160855_aH3K9me3.bw",
		HP1a_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE56101_aHP1a.bw",
		M1BP_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_M1BP_IP.bw",
		GAF_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE101557_S2_GAGA_IP.bw",
		BEAF32_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE52887_aBEAF32_IP.bw",
		CTCF_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE41354_aCTCF.bw",
		CP190_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE41354_aCP190_IP.bw",
		Su_hw_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE41354_aSu-Hw_IP.bw",
		Mod_mdg4_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE41354_aMod-mdg4.bw",
		zld_chip_classes = "results/ChIP_peak_classes/zld_ChIP_classes.tsv",
		grh_chip_classes = "results/ChIP_peak_classes/grh_ChIP_classes.tsv",
		twi_chip_classes = "results/ChIP_peak_classes/twi_ChIP_classes.tsv",
		Zld_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw",
		Zld_WT_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/S2-WT_1000uM_small.bw",
		Zld_Zld_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Zld_1000uM_small.bw",
		Grh_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw",
		Grh_WT_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/FL_ATAC_S2-WT_100uM_small.bw",
		Grh_Grh_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Grh_100uM_small.bw",
		Twi_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw",
		Twi_WT_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/Twi_ATAC_S2-WT_40uM_small.bw",
		Twi_Twi_ATAC_bw = "ATACseq/results/bigwigs/zscore_normalized/merged/S2-Twi_40uM_small.bw",
		H3K27me3_ChIP_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/individual/GSE151983_S2_aH3K27me3_IP.bw",

	output:
		"manuscript/figures/extended_data_fig5.pdf"
	script:
		"workflow/scripts/extended_data_fig5.R"
		
rule extended_data_fig_6:
	input:
		S2_Zld_ChIP_bw =   "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw",
		S2_Grh_ChIP_bw =   "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw",
		S2_Twi_ChIP_bw =   "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw",
		MNase_seq_bw =  "data/published_datasets/Chereji_2019/GSE128689_RAW/GSM3934479_S2_exp2_seq1_40min.bw",
		zld_ChIP_classes_fn = "results/ChIP_peak_classes/zld_ChIP_classes.tsv",
		grh_ChIP_classes_fn = "results/ChIP_peak_classes/grh_ChIP_classes.tsv",
		twi_ChIP_classes_fn = "results/ChIP_peak_classes/twi_ChIP_classes.tsv",
		zld_motif_instances_fn = "results/motif_instances/zld_motifs.bed",
		grh_motif_instances_fn = "results/motif_instances/grh_motifs.bed",
		twi_motif_instances_fn = "results/motif_instances/twi_motifs.bed",
	output:
		"manuscript/figures/extended_data_fig6.pdf"
	script:
		"workflow/scripts/extended_data_fig6.R"
 		
rule extended_data_fig_7:
	input:
		S2_Zld_ChIP_fn = "ChIPseq/results/peaks/final/S2-Zld_aZld_IP.bed",
		nc14_Zld_ChIP_fn = "published_ChIPseq/results/peaks/individual/narrow/embryo-nc14_aZld_summits.bed",
		brain_Zld_ChIP_fn =  "../S2_pioneers/results/bed_files/peaks/brain_aZld.bed",
		Zld_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw",
		Zld_nc14_ChIP_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/individual/embryo-nc14_aZld.bw",
		S2_Grh_ChIP_fn = "ChIPseq/results/peaks/final/S2-Grh_aGrh_IP.bed",
		embryo_Grh_ChIP_fn = "published_ChIPseq/results/peaks/individual/narrow/embryo-15-16H_aGrh_summits.bed",
		wing_Grh_ChIP_fn = "published_ChIPseq/results/peaks/individual/narrow/wing-disc_aGrh_summits.bed",
		Grh_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw",
		Grh_embryo_ChIP_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-15-16H_aGrh.bw",
		S2_twi_ChIP_fn = "ChIPseq/results/peaks/final/S2-Twi_aTwi_IP.bed",
		embryo_twi_ChIP_fn = "published_ChIPseq/results/peaks/merged/narrow/embryo-1-3H_aTwi_summits.bed",
		Twi_ChIP_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Twi_aTwi_IP.bw",
		Twi_embryo_ChIP_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/embryo-1-3H_aTwi.bw",
		H3K27me3_ChIP_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/individual/GSE151983_S2_aH3K27me3_IP.bw",
		H3K9me3_ChIP_bw = "published_ChIPseq/results/bigwigs/zscore_normalized/merged/GSE160855_aH3K9me3.bw",
	output:
		"manuscript/figures/extended_data_fig7.pdf"
	script:
		"workflow/scripts/extended_data_fig7.R"
		
rule extended_data_fig_8:
	input:
		CR_spikeIn_counts_fn ="histone_CUTandRUN/results/scaling_factors/epiCypher_barcode_counts.tsv",
		taz_blot_aH3K27me3_image = "data/immunoblot_raw_images/2021-06-29_taz/anti-H3K27me3_3.tif",
		taz_blot_aTub_image = "data/immunoblot_raw_images/2021-06-29_taz/anti-tubulin_2.tif",
	output:
		"manuscript/figures/extended_data_fig8.pdf"
	script:
		"workflow/scripts/extended_data_fig8.R"
		
rule extended_data_fig_9:
	input:
		zld_titration_blot_aZld = "data/immunoblot_raw_images/2022-7-12_titration/Zld-titration_aZld.tif",
		grh_titration_blot_aGrh = "data/immunoblot_raw_images/2022-7-12_titration/Grh-titration_aGrh_short.tif",
		grh_titration_blot_aTub = "data/immunoblot_raw_images/2022-7-12_titration/Grh-titration_aTub.tif",
		twi_titration_blot_aTwi = "data/immunoblot_raw_images/2022-7-12_titration/Twi-titration_aTwi.tif",
		twi_titration_blot_aTub = "data/immunoblot_raw_images/2022-7-12_titration/Twi-titration_aTub.tif",
		zld_tissue_classes_fn = "results/ChIP_tissue_classes/zld_tissue_classes.tsv",
		grh_tissue_classes_fn = "results/ChIP_tissue_classes/grh_tissue_classes.tsv",
		twi_tissue_classes_fn = "results/ChIP_tissue_classes/twi_tissue_classes.tsv",
	output:
		"manuscript/figures/extended_data_fig9.pdf"
	script:
		"workflow/scripts/extended_data_fig9.R"
		
rule extended_data_fig_10:
	input:
		zld_DBD_blot = "data/immunoblot_raw_images/2021-10-14_Zld-DBD/Zld-Grh-DBDs.tif",
		grh_DBD_blot = "data/immunoblot_raw_images/2021-08-19_Grh-DBD/Grh-DBD.tif",
		zld_FL_ChIP_classes_fn = "results/ChIP_peak_classes/zld_ChIP_classes.tsv",
		zld_DBD_peaks_fn = "ChIPseq/results/peaks/filtered/S2-Zld-DBD_aZld_IP.narrowPeak",
		grh_FL_ChIP_classes_fn = "results/ChIP_peak_classes/grh_ChIP_classes.tsv",
		grh_DBD_peaks_fn = "ChIPseq/results/peaks/filtered/S2-Grh-DBD_aGrh_IP.narrowPeak",
		zld_FL_bw =  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld_aZld_IP.bw",
		zld_DBD_bw =  "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Zld-DBD_aZld_IP.bw",
		grh_FL_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh_aGrh_IP.bw",
		grh_DBD_bw = "ChIPseq/results/bigwigs/zscore_normalized/merged/S2-Grh-DBD_aGrh_IP.bw",	
		S2_Zld_FL_aZld_GFP = "data/DBD_immunostaining_raw_images/20211027 ZLD-FL_GFP_RGB_470X 520M (GFP).tif",
		S2_Zld_FL_aZld_DAPI = "data/DBD_immunostaining_raw_images/20211027 ZLD-FL_DAPI_RGB_395X 445M (DAPI).tif",
		S2_Zld_FL_aZld_BF = "data/DBD_immunostaining_raw_images/20211027 ZLD-FL_BF_RGB_BF.tif",
		S2_Zld_DBD_aZld_GFP = "data/DBD_immunostaining_raw_images/20211027 ZLD-DBD_GFP_RGB_470X 520M (GFP).tif",
		S2_Zld_DBD_aZld_DAPI = "data/DBD_immunostaining_raw_images/20211027 ZLD-DBD_DAPI_RGB_395X 445M (DAPI).tif",
		S2_Zld_DBD_aZld_BF = "data/DBD_immunostaining_raw_images/20211027 ZLD-DBD_BF_RGB_BF.tif",
		S2_WT_aZld_GFP = "data/DBD_immunostaining_raw_images/20211027 WT anti-ZLD_GFP_RGB_470X 520M (GFP).tif",
		S2_WT_aZld_DAPI = "data/DBD_immunostaining_raw_images/20211027 WT anti-ZLD_DAPI_RGB_395X 445M (DAPI).tif",
		S2_WT_aZld_BF = "data/DBD_immunostaining_raw_images/20211027 WT anti-ZLD_BF_RGB_BF.tif",
		S2_Grh_FL_aGrh_GFP = "data/DBD_immunostaining_raw_images/20211027 Grh-FL_GFP_RGB_470X 520M (GFP).tif",
		S2_Grh_FL_aGrh_DAPI = "data/DBD_immunostaining_raw_images/20211027 Grh-FL_DAPI_RGB_395X 445M (DAPI).tif",
		S2_Grh_FL_aGrh_BF = "data/DBD_immunostaining_raw_images/20211027 Grh-FL_BF_RGB_BF.tif",
		S2_Grh_DBD_aGrh_GFP = "data/DBD_immunostaining_raw_images/20211027 GRH-DBD_GFP_RGB_470X 520M (GFP).tif",
		S2_Grh_DBD_aGrh_DAPI = "data/DBD_immunostaining_raw_images/20211027 GRH-DBD_DAPI_RGB_395X 445M (DAPI).tif",
		S2_Grh_DBD_aGrh_BF = "data/DBD_immunostaining_raw_images/20211027 GRH-DBD_BF_RGB_BF.tif",
		S2_WT_aGrh_GFP = "data/DBD_immunostaining_raw_images/20211027 WT anti-GRH_GFP_RGB_470X 520M (GFP).tif",
		S2_WT_aGrh_DAPI =  "data/DBD_immunostaining_raw_images/20211027 WT anti-GRH_DAPI_RGB_395X 445M (DAPI).tif",
		S2_WT_aGrh_BF = "data/DBD_immunostaining_raw_images/20211027 WT anti-GRH_BF_RGB_BF.tif",
	output:
		"manuscript/figures/extended_data_fig10.pdf"
	script:
		"workflow/scripts/extended_data_fig10.R"

##########################################################################################
##########################################################################################
###                                                                                    ###
###                               target rules                                         ###
###                                                                                    ###
##########################################################################################
##########################################################################################
rule all:
    input:
    	rules.ChIPseq_all.input,
        rules.published_ChIPseq_all.input,
        rules.RNAseq_all.input,
        rules.ATACseq_all.input,
        rules.histone_CUTandRUN_all.input,
        rules.TF_CUTandRUN_all.input,
        "manuscript/figures/fig1.pdf",
        "manuscript/figures/fig2.pdf",
        "manuscript/figures/fig3.pdf",
        "manuscript/figures/fig4.pdf",
        "manuscript/figures/fig5.pdf",
        "manuscript/figures/fig6.pdf",
        "manuscript/figures/extended_data_fig1.pdf",
        "manuscript/figures/extended_data_fig2.pdf",
        "manuscript/figures/extended_data_fig3.pdf",
        "manuscript/figures/extended_data_fig4.pdf",
        "manuscript/figures/extended_data_fig5.pdf",
        "manuscript/figures/extended_data_fig6.pdf",
        "manuscript/figures/extended_data_fig7.pdf",
        "manuscript/figures/extended_data_fig8.pdf",
        "manuscript/figures/extended_data_fig9.pdf",
        "manuscript/figures/extended_data_fig10.pdf",
    default_target: True