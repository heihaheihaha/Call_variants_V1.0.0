import os
import re
import sys

configfile: "config.yaml"
gatk = config["gatk_path"]
ref_path0 = config["reference_panel_path"]
def fa_dict(path0: str) -> str:
	dir_path = os.path.dirname(path0)
	file_name = os.path.splitext(os.path.basename(path0))[0]
	output_file = os.path.join(dir_path, file_name + ".dict")
	return output_file

ref_dict_path = fa_dict(ref_path0)

db = ""
for i in config["category"].split("-"):
    db += f"-{i} Y "

if ('ASD' in {config['panel_type']}):
	getVar2_flag = 1
else:
	getVar2_flag = 0

rule end:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.annovar", # trigger the rule HaplotypeCaller, may change later
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam.table.pdf", # trigger the rule AnalyzeCovariates
		f"{config['output_dir']}/variants/{config['sample_name']}.stat.xlsx",
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
	threads: 16
	shell: #f"bwa mem -t {threads} {input.reference} {input.R1_path} {input.R2_path} | samtools sort > {output}" # bwa mem will align the reads to the reference panel
		"bwa mem -t {threads}" + f" {config['reference_panel_path']} {config['R1_path']} {config['R2_path']} | samtools sort > {config['output_dir']}/alignments/{config['sample_name']}.bwa.bam" # bwa mem will align the reads to the reference panel

rule AddOrReplaceReadGroups: # Provide information for BQSR
	input:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.bam"
	output:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.rg.bam"
	shell:
		f"""{gatk} AddOrReplaceReadGroups \\
			I={config['output_dir']}/alignments/{config['sample_name']}.bwa.bam \\
			O={config['output_dir']}/alignments/{config['sample_name']}.bwa.rg.bam \\
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

rule mark_duplicates:
	input:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.rg.bam"
	output:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam",
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.metrics"
	shell:
		f"""{gatk} MarkDuplicatesSpark \\
			-I {config['output_dir']}/alignments/{config['sample_name']}.bwa.rg.bam \\
			-O {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam \\
			-M {config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.metrics""" # sambamba markdup will mark the duplicates in the bam file


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

rule BaseRecalibrator:
	input: 
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bai",
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam",
		f"{config['reference_panel_path']}",
		f"{config['reference_panel_path']}.fai",
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
		f"""{gatk} ApplyBQSRSpark -R {config['reference_panel_path']} \\
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

rule BaseRecalibrator2:
	input:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam.bai",
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam",
		f"{config['reference_panel_path']}",
		f"{config['known_sites']}",
		f"{config['known_sites']}.tbi"
	output:
		f"{config['output_dir']}/alignments/{config['sample_name']}.bwa.markdup.rg.bam.bqsr.bam.table"
	shell:
		f"""{gatk} BaseRecalibratorSpark \\
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

rule selectvariants_snp:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.vcf.gz",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.snp.vcf.gz"
	shell:
		f"""{gatk} SelectVariants \\
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/{config['sample_name']}.vcf.gz \\
			-select-type-to-include SNP \\
			-O {config['output_dir']}/variants/{config['sample_name']}.snp.vcf.gz"""

rule VariantFiltration:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.snp.vcf.gz",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.snp.passed.vcf.gz"
	shell:
		f"""{gatk} VariantFiltration \\
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/{config['sample_name']}.snp.vcf.gz \\
			-O {config['output_dir']}/variants/{config['sample_name']}.snp.passed.vcf.gz \\
			--filter-expression 'QUAL<30.0' --filter-name 'LOW_QUAL' \\
			--filter-expression 'QD<2.0' --filter-name 'LOW_QD' \\
			--filter-expression 'FS>60.0' --filter-name 'HIGH_FS' \\
			--filter-expression 'MQ<40.0' --filter-name 'LOW_MQ' \\
			--filter-expression 'MQRankSum<-12.5' --filter-name 'LOW_MQRS' \\
			--filter-expression 'ReadPosRankSum<-8.0' --filter-name 'LOW_RPRS' \\
			--filter-expression 'SOR>3.0' --filter-name 'HIGH_SOR' 
			"""
			# --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \\

### Variant Quality Score Recalibration (VQSR)? ###

rule perl_filter:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.snp.passed.vcf.gz"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.snp.filtered.vcf.gz"
	shell:
		r"""perl -ne 'chomp;if($_=~/^#/ || $_ =~ /PASS/){{print "$_\n"}}'""" + f" {config['output_dir']}/variants/{config['sample_name']}.snp.passed.vcf.gz" + ">" + f"{config['output_dir']}/variants/{config['sample_name']}.snp.filtered.vcf.gz"

rule FilterVar:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.snp.filtered.vcf.gz"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.snp.annovar"
	shell:
		f"""perl ./FilterVar.pl \\
		-in {config['output_dir']}/variants/{config['sample_name']}.snp.filtered.vcf.gz \\
		-out {config['output_dir']}/variants/{config['sample_name']}.snp.annovar \\
		-minhet 0.20 --wildsample -qual 30 -dp 20 -adp 10 -gq 30 --rough"""
	
# filter indel
rule selectvariants_indel:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.vcf.gz",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.indel.vcf.gz"
	shell:
		f"""{gatk} SelectVariants \\
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/{config['sample_name']}.vcf.gz \\
			-select-type-to-include INDEL \\
			-O {config['output_dir']}/variants/{config['sample_name']}.indel.vcf.gz"""

rule VariantFiltration_indel:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.indel.vcf.gz",
		f"{config['reference_panel_path']}"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.indel.passed.vcf.gz"
	shell:
		f"""{gatk} VariantFiltration \\
			-R {config['reference_panel_path']} \\
			-V {config['output_dir']}/variants/{config['sample_name']}.indel.vcf.gz \\
			-O {config['output_dir']}/variants/{config['sample_name']}.indel.passed.vcf.gz \\
			--filter-expression 'QUAL<30.0' --filter-name 'LOW_QUAL' \\
			--filter-expression 'QD<2.0' --filter-name 'LOW_QD' \\
			--filter-expression 'FS>200.0' --filter-name 'HIGH_FS' \\
			--filter-expression 'ReadPosRankSum<-20.0' --filter-name 'LOW_RPRS' \\
			--filter-expression 'SOR>10.0' --filter-name 'HIGH_SOR' 
			"""

rule perl_filter_indel:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.indel.passed.vcf.gz"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.indel.filtered.vcf.gz"
	shell:
		r"""perl -ne 'chomp;if($_=~/^#/ || $_ =~ /PASS/){{print "$_\n"}}'""" + f" {config['output_dir']}/variants/{config['sample_name']}.indel.passed.vcf.gz" + ">" + f"{config['output_dir']}/variants/{config['sample_name']}.indel.filtered.vcf.gz"

rule FilterVar_indel:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.indel.filtered.vcf.gz"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.indel.annovar"
	shell:
		f"""perl ./FilterVar.pl \\
		-in {config['output_dir']}/variants/{config['sample_name']}.indel.filtered.vcf.gz \\
		-out {config['output_dir']}/variants/{config['sample_name']}.indel.annovar \\
		-minhet 0.30 --wildsample -qual 30 -dp 20 -adp 10 -gq 30 --rough"""

rule CombineVariants:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.snp.annovar",
		f"{config['output_dir']}/variants/{config['sample_name']}.indel.annovar"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.annovar"
	shell:
		f"""cat {config['output_dir']}/variants/{config['sample_name']}.snp.annovar {config['output_dir']}/variants/{config['sample_name']}.indel.annovar > {config['output_dir']}/variants/{config['sample_name']}.annovar"""

rule ExtremeVar:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.annovar",
		f"{config['ExtremeVar_PATH']}"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.initial.extreme.xls",
		f"{config['output_dir']}/variants/{config['sample_name']}.initial.hg38_multianno.txt"
	shell:
		# "$ExtremeVar -in $dir[2]/$id2name{$i}.annovar -out $id2name{$i}.initial -Psoft $opt{Psoft} -maf $opt{MAFS} -reference $opt{ref} --extreme --extreme_all --remove -database $db_path  $db\n",
		f"""perl {config['ExtremeVar_PATH']} \\
		-in {config['output_dir']}/variants/{config['sample_name']}.annovar \\
		-out {config['output_dir']}/variants/{config['sample_name']}.initial \\
		-Psoft {config['Psoft']} \\
		-maf {config['MAFS']} \\
		-reference {config['reference_panel_name']} \\
		--extreme \\
		--extreme_all \\
		--remove \\
		-database {config['database_path']}  {db}
		"""

#"$ExtremeVar2 -out $id2name{$i}.initial --extreme --extreme_all -MPsoft $opt{MPsoft2} -reference $opt{ref}"
rule ExtremeVar2:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.initial.extreme.xls",
		f"{config['output_dir']}/variants/{config['sample_name']}.initial.hg38_multianno.txt",
		f"{config['ExtremeVar2_PATH']}"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.initial.tmp.xls"
	shell:
		f"""perl {config['ExtremeVar2_PATH']} \\
		-out {config['output_dir']}/variants/{config['sample_name']}.initial \\
		--extreme \\
		--extreme_all \\
		-MPsoft {config['MPsoft2']} \\
		-reference {config['reference_panel_name']}
		"""

rule awk_select:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.initial.tmp.xls"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.initial.xls"
	shell:
		r"""awk -F \"\\t\" 'BEGIN{{IGNORECASE=1}} NR==1 {{print $0}} NR>1 {{if($1~/Y/ || $105~/pathogenic/ || $105~/drug_response/ || ($107~/DM/ && $105 !~/benign/)) print $0}}' """+f"""{config['output_dir']}/variants/{config['sample_name']}.initial.tmp.xls"""+" | cut -f 1-10,12-26,30-33,37,41,45,51-53,55-60,64,69,79,85,88,93-94,97,100-111,115 > "+f"{config['output_dir']}/variants/{config['sample_name']}.initial.xls"

rule progress:
	input:
		f"{config['output_dir']}/variants/{config['sample_name']}.initial.xls",
		f"{config['progress_PATH']}", # perl script
		f"{config['TRANSCRIPT']}",
		f"{config['GENEMODEL']}",
		f"{config['HPO']}"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.extreme.xls"
	shell:
	# $process $id2name{$i}.initial.xls $transcript $gm $hpo > $id2name{$i}.extreme.xls
		f"perl {config['progress_PATH']} {config['output_dir']}/variants/{config['sample_name']}.initial.xls {config['TRANSCRIPT']} {config['GENEMODEL']} {config['HPO']}>{config['output_dir']}/variants/{config['sample_name']}.extreme.xls"

rule ACMG:
	input:
		f"{config['ACMG_perl_script_PATH']}",
		f"{config['output_dir']}/variants/{config['sample_name']}.extreme.xls"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.prefinal.xls"
	shell:
		f"perl {config['ACMG_perl_script_PATH']} -in {config['output_dir']}/variants/{config['sample_name']}.extreme.xls -out {config['output_dir']}/variants/{config['sample_name']}.prefinal.xls"

rule getVar2:
	input:
		f"{config['getVar2_script_PATH']}",
		f"{config['output_dir']}/variants/{config['sample_name']}.prefinal.xls"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.clinvar_HGMD.xls",
		f"{config['output_dir']}/variants/{config['sample_name']}.loose.xls",
		f"{config['output_dir']}/variants/{config['sample_name']}.strict.xls",
		f"{config['output_dir']}/variants/{config['sample_name']}.final.xls"
	shell:
		f"""if [{getVar2_flag}]
		then
			{config['getVar2_script_PATH']} -in {config['output_dir']}/variants/{config['sample_name']}.prefinal.xls -list {config['panel_PATH']} -acmglist {config['acmg78']} -o1 {config['output_dir']}/variants/{config['sample_name']}.clinvar_HGMD.xls -o2 {config['output_dir']}/variants/{config['sample_name']}.loose.xls -o3 {config['output_dir']}/variants/{config['sample_name']}.strict.xls -o4 {config['output_dir']}/variants/{config['sample_name']}.final.xls
		else
			perl {config['getVar2_script_PATH']} -in {config['output_dir']}/variants/{config['sample_name']}.prefinal.xls -acmglist {config['acmg78']} -o1 {config['output_dir']}/variants/{config['sample_name']}.clinvar_HGMD.xls -o2 {config['output_dir']}/variants/{config['sample_name']}.loose.xls -o3 {config['output_dir']}/variants/{config['sample_name']}.strict.xls -o4 {config['output_dir']}/variants/{config['sample_name']}.final.xls
		fi"""

rule get_result:
	input:
		f"{config['get_result_perl_script_PATH']}",
		f"{config['output_dir']}/variants/{config['sample_name']}.clinvar_HGMD.xls",
		f"{config['output_dir']}/variants/{config['sample_name']}.loose.xls",
		f"{config['output_dir']}/variants/{config['sample_name']}.strict.xls",
		f"{config['output_dir']}/variants/{config['sample_name']}.final.xls"
	output:
		f"{config['output_dir']}/variants/{config['sample_name']}.stat.xls"
	shell:
		f"""perl {config['get_result_perl_script_PATH']} \\
		{config['output_dir']}/variants/{config['sample_name']}.strict.xls \\
		{config['output_dir']}/variants/{config['sample_name']}.clinvar_HGMD.xls \\
		{config['output_dir']}/variants/{config['sample_name']}.loose.xls \\
		{config['output_dir']}/variants/{config['sample_name']}.final.xls \\
		{config['output_dir']}/variants/{config['sample_name']}.stat.xls"""



# rule get_result:
# input:
# 	f"{config['get_result_perl_script_PATH']}",
# 	f"{config['output_dir']}/variants/{config['sample_name']}.strict.xls",
# 	f"{config['output_dir']}/variants/{config['sample_name']}.clinvar_HGMD.xls"
# output:

# shell:
# 	f"perl {config['get_result_perl_script_PATH']} {config['output_dir']}/variants/{config['sample_name']}.strict.xls {config['output_dir']}/variants/{config['sample_name']}.clinvar_HGMD.xls "
