import os

configfile: "config.yaml"
gatk = config["gatk_path"]
ref_path0 = config["reference_panel_path"]
def fa_dict(path0: str) -> str:
	dir_path = os.path.dirname(path0)
	file_name = os.path.splitext(os.path.basename(path0))[0]
	output_file = os.path.join(dir_path, file_name + ".dict")
	return output_file

ref_dict_path = fa_dict(ref_path0)

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
		f"../results/fastqc/{config['sample_name']}_fastqc.html"
		f"../results/fastqc/{config['sample_name']}_fastqc.zip"
	conda:
		"./first_step_mamba.yml"
	shell:
		f"{config['path_to_mamba_env']}/bin/fastqc {input} -o ../results/fastqc/"

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
		f"../results/alignments/{config['sample_name']}.bwa.bam"
	conda:
		"./first_step_mamba.yml"
	threads: 8
	shell: #f"bwa mem -t {threads} {input.reference} {input.R1_path} {input.R2_path} | samtools sort > {output}" # bwa mem will align the reads to the reference panel
		f"{config['path_to_mamba_env']}/bin/bwa mem {config['reference_panel_path']} {config['R1_path']} {config['R2_path']} | {config['path_to_mamba_env']}/bin/samtools sort > ../results/alignments/{config['sample_name']}.bwa.bam" # bwa mem will align the reads to the reference panel

rule mark_duplicates:
	input:
		f"../results/alignments/{config['sample_name']}.bwa.bam"
	output:
		f"../results/alignments/{config['sample_name']}.bwa.markdup.bam",
		f"../results/alignments/{config['sample_name']}.bwa.markdup.metrics"
	shell:
		f"{gatk} MarkDuplicates I=../results/alignments/{config['sample_name']}.bwa.bam O=../results/alignments/{config['sample_name']}.bwa.markdup.bam M=../results/alignments/{config['sample_name']}.bwa.markdup.metrics" # sambamba markdup will mark the duplicates in the bam file

rule ADD_READ_GROUPS:
	input:
		f"../results/alignments/{config['sample_name']}.bwa.markdup.bam"
	output:
		f"../results/alignments/{config['sample_name']}.bwa.markdup.rg.bam"
	shell:
		f"{gatk} AddOrReplaceReadGroups I=../results/alignments/{config['sample_name']}.bwa.markdup.bam O=../results/alignments/{config['sample_name']}.bwa.markdup.rg.bam RGID={config['sample_name']} RGLB={config['sample_name']} RGPL=illumina RGPU=unit1 RGSM={config['sample_name']}" # gatk add read groups will add the read groups to the bam file
	
rule dict_index:
	input:
		f"{config['reference_panel_path']}"
	output:
		f"{ref_dict_path}"
	shell:
		f"{gatk} CreateSequenceDictionary REFERENCE={config['reference_panel_path']} OUTPUT={ref_dict_path}"

rule index_bam:
	input:
		f"../results/alignments/{config['sample_name']}.bwa.markdup.rg.bam"
	output:
		f"../results/alignments/{config['sample_name']}.bwa.markdup.rg.bai"
	shell:
		f"{config['path_to_mamba_env']}/bin/samtools index ../results/alignments/{config['sample_name']}.bwa.markdup.rg.bam" # gatk build bam index will create the index file for the bam file

rule variant_calling:
	input:
		f"../results/alignments/{config['sample_name']}.bwa.markdup.rg.bam",
		f"../results/alignments/{config['sample_name']}.bwa.markdup.rg.bai",
		f"{config['reference_panel_path']}",
		f"{ref_dict_path}"
	output:
		f"../results/variants/{config['sample_name']}.g.vcf.gz"
	shell:
		f"{gatk} HaplotypeCaller -R {config['reference_panel_path']} -I ../results/alignments/{config['sample_name']}.bwa.markdup.rg.bam -O ../results/variants/{config['sample_name']}.g.vcf.gz -ERC GVCF" # gatk haplotype caller will call the variants in the bam file

rule :
	input:
		f"../results/variants/{config['sample_name']}.g.vcf.gz"
	output:
		"../results/variants/.sentinel"
	shell:
		"echo 'emm' > ../results/variants/.sentinel"
