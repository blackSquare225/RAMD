rule count_500kbp_reads:
	input:
		script="workflow/scripts/count_fastq_files.sh",
		fusion_dir=lambda wildcards: directory(config["fusion_dir"])
	output:
		"data/{sample}/targeted/fusion/500kbp/{sample}.fusion.500kbp.fastq.gz",
	threads: 
		1
	log:
		"logs/{sample}/{sample}.count_500kbp_reads.log"
	benchmark:
		"benchmarks/{sample}/{sample}.count_500kbp_reads.tsv"
	shell:	
		"bash {input.script} --input_dir {input.fusion_dir} --output_fq {output} --n_files 125 2> {log}"

rule count_1Mbp_reads:
	input:
		script="workflow/scripts/count_fastq_files.sh",
		fusion_dir=lambda wildcards: directory(config["fusion_dir"])
	output:
		"data/{sample}/targeted/fusion/1Mbp/{sample}.fusion.1Mbp.fastq.gz"
	threads: 
		1
	log:
		"logs/{sample}/{sample}.count_1Mbp_reads.log"
	benchmark:
		"benchmarks/{sample}/{sample}.count_1Mbp_reads.tsv"
	shell:	
		"bash {input.script} --input_dir {input.fusion_dir} --output_fq {output} --n_files 250 2> {log}"

rule jaffal_500kbp:
	input:
		jaffal_bpipe="resources/JAFFA-version-2.3/tools/bin/bpipe",
		jaffal_groovy="resources/JAFFA-version-2.3/JAFFAL.groovy",
		fq="data/{sample}/targeted/fusion/500kbp/{sample}.fusion.500kbp.fastq.gz"
	output:
		"results/{sample}/targeted/fusion/500kbp/jaffa_results.csv"
	params:
		outdir="results/{sample}/targeted/fusion/500kbp",
	threads: 
		1
	conda:
		"../envs/jaffa.yaml"
	log:
		"logs/{sample}/{sample}.jaffal_500kbp.log"
	benchmark:
		"benchmarks/{sample}/{sample}.jaffal_500kbp.tsv"
	shell:
		"""
		mkdir -p {params.outdir} && \
		cd {params.outdir} && \
		../../../../../{input.jaffal_bpipe} \
		../../../../../{input.jaffal_groovy} \
		../../../../../{input.fq}
		"""

rule jaffal_1Mbp:
	input:
		jaffal_bpipe="resources/JAFFA-version-2.3/tools/bin/bpipe",
		jaffal_groovy="resources/JAFFA-version-2.3/JAFFAL.groovy",
		fq="data/{sample}/targeted/fusion/1Mbp/{sample}.fusion.1Mbp.fastq.gz"
	output:
		"results/{sample}/targeted/fusion/1Mbp/jaffa_results.csv"
	params:
		outdir="results/{sample}/targeted/fusion/1Mbp",
	threads:
		1
	conda:
		"../envs/jaffa.yaml"
	log:
		"logs/{sample}/{sample}.jaffal_1Mbp.log"
	benchmark:
		"benchmarks/{sample}/{sample}.jaffal_1Mbp.tsv"
	shell:
		"""
		mkdir -p {params.outdir} && \
		cd {params.outdir} && \
		../../../../../{input.jaffal_bpipe} \
		../../../../../{input.jaffal_groovy} \
		../../../../../{input.fq} 
		"""