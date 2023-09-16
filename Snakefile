import os
import re

configfile: "config.yaml"
gatk = config["gatk_path"]
ref_path0 = config["reference_panel_path"]
def fa_dict(path0: str) -> str:
	dir_path = os.path.dirname(path0)
	file_name = os.path.splitext(os.path.basename(path0))[0]
	output_file = os.path.join(dir_path, file_name + ".dict")
	return output_file

ref_dict_path = fa_dict(ref_path0)

rule end:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.filtered.vcf.gz", # trigger the rule HaplotypeCaller, may change later
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam.table.pdf", # trigger the rule AnalyzeCovariates
		f"{ref_dict_path}" # trigger the rule dict_index

# Prepare environment
rule index_reference:
	input:
		f"{config['reference_panel_path']}"
	output:
		f"{config['reference_panel_path']}.bwt",
		f"{config['reference_panel_path']}.pac",
		f"{config['reference_panel_path']}.ann",
		f"{config['reference_panel_path']}.sa",
		f"{config['reference_panel_path']}.amb"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"bwa index {config['reference_panel_path']}" # bwa index will create the index files with the fasta file
		
rule fastqc: # output is the html and zip files
	input:
		f"{config['R1_path']}",
		f"{config['R2_path']}"
	output:
		f"{config['output_dir']}/fastqc/{config['sample_name']}_fastqc.html"
		f"{config['output_dir']}/fastqc/{config['sample_name']}_fastqc.zip"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"fastqc {input} -o {config['output_dir']}/fastqc/"

rule align_readers:
	input:
		f"{config['reference_panel_path']}.bwt",
		f"{config['reference_panel_path']}.pac",
		f"{config['reference_panel_path']}.ann",
		f"{config['reference_panel_path']}.sa",
		f"{config['reference_panel_path']}.amb",
		f"{config['reference_panel_path']}",
		f"{config['R1_path']}",
		f"{config['R2_path']}"
	output:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.bam"
	conda:
		"./first_step_mamba.yml"
	threads: 8
	shell: #f"bwa mem -t {threads} {input.reference} {input.R1_path} {input.R2_path} | samtools sort > {output}" # bwa mem will align the reads to the reference panel
		f"""bwa mem {config['reference_panel_path']} {config['R1_path']} {config['R2_path']} | \\
		samtools sort > {config['output_dir']}/alignments/{config['sample_name']}.bwa.bam""" # bwa mem will align the reads to the reference panel

rule mark_duplicates:
	input:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.bam"
	output:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.bam",
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.metrics"
	shell:
		f"""{gatk} MarkDuplicates \\
			I={config['output_dir']}/alignments/{config['sample_name']}.bwa.bam \\
			O={config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.bam \\
			M={config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.metrics""" # sambamba markdup will mark the duplicates in the bam file

rule AddOrReplaceReadGroups: # Provide information for BQSR
	input:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.bam"
	output:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam"
	shell:
		f"""{gatk} AddOrReplaceReadGroups \\
			I={config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.bam \\
			O={config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam \\
			RGID={config['sample_name']} \\
			RGLB=lib1 \\
			RGPL={config['platform']} \\
			RGPU={config['lane']} \\
			RGSM={config['sample_name']} \\
			CREATE_INDEX=True
"""

			#--RGPU -PU: Read-Group platform unit (eg. run barcode)
			#--RGLB -LB: Read-Group library
			#--RGPL -PL: Read-Group platform (eg. illumina, solid)
			#--RGSM -SM: Read-Group sample name
			#--RGID -ID: Read-Group ID
			#--RGCN -CN: Read-Group sequencing center name
			#--RGDS -DS: Read-Group description
			#--RGDT -DT: Read-Group run date
			#--RGPI -PI: Read-Group predicted insert size
			#--RGPG -PG: Read-Group program group
			#--RGPM -PM: Read-Group platform model
			#--RGKS -KS: Read-Group key sequence
			#--RGFO -FO: Read-Group flow order

# PU = Platform Unit
# The PU holds three types of information, 
# the {FLOWCELL_BARCODE}.{LANE}.{SAMPLE_BARCODE}. 
# The {FLOWCELL_BARCODE} refers to the unique identifier for a particular flow cell. 
# The {LANE} indicates the lane of the flow cell and the {SAMPLE_BARCODE} is a sample/library-specific identifier. 
# Although the PU is not required by GATK but takes precedence over ID for base recalibration if it is present. 
# In the example shown earlier, two read group fields, ID and PU, appropriately differentiate flow cell lane, marked by .2, a factor that contributes to batch effects.

rule dict_index:
	input:
		f"{config['reference_panel_path']}"
	output:
		f"{ref_dict_path}"
	shell:
		f"{gatk} CreateSequenceDictionary REFERENCE={config['reference_panel_path']} OUTPUT={ref_dict_path}"

rule index_bam: 
	input:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam"
	output:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bai"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"samtools index {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam" # gatk build bam index will create the index file for the bam file

rule IndexFeatureFile: # index the known sites(dbsnp, vcf file, gz)
	input:
		f"{config['known_sites']}"
	output:
		f"{config['known_sites']}.tbi"
	shell:
		f"""{gatk} IndexFeatureFile \\
			-I {config['known_sites']}"""

rule BaseRecalibrator:
	input: 
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bai",
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam",
		f"{config['reference_panel_path']}",
		f"{config['known_sites']}",
		f"{config['known_sites']}.tbi"
	output: 
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.table"
	shell:
		f"""{gatk} BaseRecalibratorSpark \\
			-R {config['reference_panel_path']} \\
			-I {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam \\
			-O {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.table \\
			--known-sites {config['known_sites']}""" 

rule ApplyBQSR:
	input:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.table"
	output:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam"
	shell:
		f"""{gatk} ApplyBQSR -R {config['reference_panel_path']} \\
			-I {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam \\
			-O {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam \\
			--bqsr-recal-file {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.table"""

rule index_bam2:
	input:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam"
	output:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam.bai"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"samtools index {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam" # gatk build bam index will create the index file for the bam file

rule fasta_faidx:
	input:
		f"{config['reference_panel_path']}"
	output:
		f"{config['reference_panel_path']}.fai"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"samtools faidx {config['reference_panel_path']}"
	# samtools faidx will create the index file for the fasta file

rule BaseRecalibrator2:
	input:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam.bai",
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam",
		f"{config['reference_panel_path']}.fai",
		f"{config['reference_panel_path']}",
		f"{config['known_sites']}",
		f"{config['known_sites']}.tbi"
	output:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam.table"
	shell:
		f"""{gatk} BaseRecalibrator \\
			-R {config['reference_panel_path']} \\
			-I {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam \\
			-O {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam.table \\
			--known-sites {config['known_sites']}"""

rule AnalyzeCovariates:
	input:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam.table"
	output:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam.table.pdf"
	conda:
		"./gatk_R_plot.yml"
	shell:
		f"""{gatk} AnalyzeCovariates \\
			-before {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.table \\
			-after {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam.table \\
			-plots {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam.table.pdf"""


rule HaplotypeCaller:
	input:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam",
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam.bai",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.g.vcf.gz"
	shell:
		f"""{gatk} HaplotypeCallerSpark \\
			-R {config['reference_panel_path']} \\
			-I {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam \\
			-O {config['output_dir']}/variants/{config['sample_name']}.g.vcf.gz \\
			-ERC GVCF""" 
# Notice: HaplotypeCaller GVCF 
rule GenotypeGVCFs:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.g.vcf.gz",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.vcf.gz"
	shell:
		f"""{gatk} GenotypeGVCFs \\
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/{config['sample_name']}.g.vcf.gz \\
			-O {config['output_dir']}/variants/{config['sample_name']}.vcf.gz"""

rule VariantFiltration:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.vcf.gz",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.filtered.vcf.gz"
	shell:
		f"""{gatk} VariantFiltration \\
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/{config['sample_name']}.vcf.gz \\
			-O {config['output_dir']}/variants/{config['sample_name']}.filtered.vcf.gz \\
			--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \\
			--filter-name "my_snp_filter" """
	# QD<2.0、FS>60.0、SOR>3.0、MQ<40.0、MQRankSum < -12.5 和 ReadPosRankSum < -8.0

			# --filter-expression 'QUAL<30.0' --filter-name 'LOW_QUAL' \\
			# --filter-expression 'QD<2.0' --filter-name 'LOW_QD' \\
			# --filter-expression 'FS>60.0' --filter-name 'HIGH_FS' \\
			# --filter-expression 'MQ<40.0' --filter-name 'LOW_MQ' \\
			# --filter-expression 'MQRankSum<-12.5' --filter-name 'LOW_MQRS' \\
			# --filter-expression 'ReadPosRankSum<-8.0' --filter-name 'LOW_RPRS' \\
			# --filter-expression 'SOR>3.0' --filter-name 'HIGH_SOR \\

### Variant Quality Score Recalibration (VQSR)? ###