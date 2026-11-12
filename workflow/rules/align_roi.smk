rule merge_fastq_files:	
	input:
		fusion_dir=lambda wildcards: directory(config["fusion_dir"]),
		roi_dir=lambda wildcards: directory(config["roi_dir"])
	output:
		"data/{sample}/targeted/roi/{sample}.targeted.fastq.gz"
	params:
		lambda wildcards: int(config["time_roi_analysis"])
	threads:
		1
	log:
		"logs/{sample}/{sample}.merge_fastq_files.log"
	benchmark:
		"benchmarks/{sample}/{sample}.merge_fastq_files.tsv"
	shell:
		"sleep $((3600 * {params})) && cat {input.fusion_dir}/*.fastq.gz {input.roi_dir}/*.fastq.gz > {output} 2>{log}"

rule minimap2_roi:
	input:
		"resources/GRCh38_full_analysis_set_plus_decoy_hla.masked.fa.mmi",
		"data/{sample}/targeted/roi/{sample}.targeted.fastq.gz"
	output:
		"alignments/{sample}/targeted/roi/{sample}.minimap2.targeted.srt.bam"
	threads:
		lambda wildcards: round(int(config["threads"]) / 2)
	log:
		"logs/{sample}/{sample}.minimap2_align_targeted.log"
	benchmark:
		"benchmarks/{sample}/{sample}.minimap2_align_targeted.tsv"
	conda: 
		"../envs/minimap2.yaml"
	shell:
		"minimap2 -ax splice -t {threads} {input} 2>{log}| samtools sort -@ {threads} -O bam -o {output} 2>>{log}"

rule ngmlr_roi:
	input:
		ref="resources/GRCh38_full_analysis_set_plus_decoy_hla.masked.fa",
		fq="data/{sample}/targeted/roi/{sample}.targeted.fastq.gz"
	output:
		"alignments/{sample}/targeted/roi/{sample}.ngmlr.targeted.srt.bam"
	threads:
		lambda wildcards: round(int(config["threads"]) / 2 - 0.5)
	log:
		"logs/{sample}/{sample}.ngmlr_samtools_alignment.log"
	benchmark:
		"benchmarks/{sample}/{sample}.ngmlr_align.tsv"
	conda: 
		"../envs/ngmlr.yaml"
	shell:
		"zcat {input.fq} | ngmlr -r {input.ref} -x ont -t {threads} 2>{log} | samtools sort -@ {threads} -o {output} - 2>>{log}"

rule samtools_index_minimap2_roi:
	input:
		"alignments/{sample}/targeted/roi/{sample}.minimap2.targeted.srt.bam"
	output:
		"alignments/{sample}/targeted/roi/{sample}.minimap2.targeted.srt.bam.bai"    
	threads:
		lambda wildcards: round(int(config["threads"]) / 2) 
	log:
		"logs/{sample}/{sample}.samtools_index_minimap2.log"
	benchmark:
		"benchmarks/{sample}/{sample}.samtools_index_minimap2.tsv" 
	conda: 
		"../envs/samtools.yaml"
	shell:
		"samtools index -@ {threads} {input}"

rule samtools_index_ngmlr_roi:
	input:
		"alignments/{sample}/targeted/roi/{sample}.ngmlr.targeted.srt.bam"
	output:
		"alignments/{sample}/targeted/roi/{sample}.ngmlr.targeted.srt.bam.bai"
	threads:
		lambda wildcards: round(int(config["threads"]) / 2 - 0.5)  
	log:
		"logs/{sample}/{sample}.samtools_index_ngmlr.log"
	benchmark:
		"benchmarks/{sample}/{sample}.samtools_index_ngmlr.tsv"	
	conda: 
		"../envs/samtools.yaml"
	shell:
		"samtools index -@ {threads} {input}"

rule samtools_view_flt3:
	input:
		bam="alignments/{sample}/targeted/roi/{sample}.ngmlr.targeted.srt.bam",
		bai="alignments/{sample}/targeted/roi/{sample}.ngmlr.targeted.srt.bam.bai"
	output:
		bam="alignments/{sample}/targeted/roi/{sample}.ngmlr.flt3.bam"
	threads: 
		lambda wildcards: round(int(config["threads"]) / 2 - 0.5)	
	log:
		"logs/{sample}/{sample}.ngmlr.samtools_view_flt3.log"
	benchmark:
		"benchmarks/{sample}/{sample}.ngmlr.samtools_view_flt3.tsv"
	conda:
		"../envs/samtools.yaml"
	shell:
		"samtools view -@ {threads} -O BAM -o {output} {input.bam} chr13:28033883-28034223"

rule samtools_index_flt3:
	input:
		"alignments/{sample}/targeted/roi/{sample}.ngmlr.flt3.bam"
	output:
		"alignments/{sample}/targeted/roi/{sample}.ngmlr.flt3.bam.bai"
	threads: 
		lambda wildcards: round(int(config["threads"]) / 2 - 0.5)
	log:
		"logs/{sample}/{sample}.samtools_index_flt3.log"
	benchmark:
		"benchmarks/{sample}/{sample}.ngmlr.samtools_index_flt3.tsv"
	conda:
		"../envs/samtools.yaml"
	shell:
		"samtools index -@ {threads} {input}"

