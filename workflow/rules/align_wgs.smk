rule merge_bam_files: 
	input:
		script="workflow/scripts/merge_bam_files.sh",
		wgs_dir=lambda wildcards: directory(config["wgs_dir"])
	output:
		temp("alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.tmp.bam")
	params:
		lambda wildcards: int(config["time_wgs_analysis"])
	threads:
		lambda wildcards: int(config["threads"])
	conda:
		"../envs/samtools.yaml"
	log:
		"logs/{sample}/{sample}.merge_bam_files.log"
	benchmark:
		"benchmarks/{sample}/{sample}.merge_bam_files.tsv"
	shell:
		"bash {input.script} --input_dir {input.wgs_dir} --output_bam {output} --threads {threads} --time_wgs_analysis {params} 2> {log}"

#rule samtools_sort_wgs:
#	input:
#		"alignments/{sample}/wgs/{sample}.minimap2.wgs.bam"
#	output:
#		"alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam"
#	threads:
#		lambda wildcards: int(config["threads"])
#	conda:
#		"../envs/samtools.yaml"
#	log:
#		"logs/{sample}/{sample}.samtools_sort_wgs.log"
#	benchmark:
#		"benchmarks/{sample}/{sample}.samtools_sort_wgs.tsv"
#	shell:
#		"samtools sort -@ {threads} -o {output} {input} 2> {log}"

rule filter_bam_reads:
	input:
		script="workflow/scripts/filter_bam_reads.py",
		bam="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.tmp.bam"
	output:
		bam="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam"
	conda:
		"../envs/filter_bam_reads.yaml"
	log:
		"logs/{sample}/{sample}.filter_bam_reads.log"
	benchmark:
		"benchmarks/{sample}/{sample}.filter_bam_reads.tsv"
	shell:
		"python {input.script} --input {input.bam} --output {output.bam} 2> {log}"


rule samtools_index_wgs:
	input:
		"alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam"
	output:
		"alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam.bai"
	threads: 
		lambda wildcards: int(config["threads"])
	conda: 
		"../envs/samtools.yaml"
	log:
		"logs/{sample}/{sample}.samtools_index_wgs.log"
	benchmark:
		"benchmarks/{sample}/{sample}.samtools_index_wgs.tsv"
	shell:
		"samtools index -@ {threads} {input} 2> {log}"


#rule minimap2_wgs:
#	input:
#		mmi="resources/GRCh38_full_analysis_set_plus_decoy_hla.fa",
#		fq="data/{sample}/wgs/{sample}.wgs.fastq.gz"
#	output:
#		"alignments/{sample}/wgs/{sample}.minimap2.{genome}.srt.bam"
#	threads: 
#		lambda wildcards: int(config["threads"])
#	benchmark:
#		"benchmarks/{sample}/{sample}.minimap2_align_{genome}.txt"
#	log:
#		"logs/{sample}/{sample}.minimap2_align_{genome}.log"
#	conda: 
#		"../envs/minimap2.yaml"
#	shell:
#		"minimap2 --MD -ax lr:hq -t {threads} {input.mmi} {input.fq} 2>{log}| samtools sort -@ {threads} -o {output} -  2>>{log}"



