rule sniffles_flt3:
	input:
		bam="alignments/{sample}/targeted/roi/{sample}.ngmlr.flt3.bam",
		bai="alignments/{sample}/targeted/roi/{sample}.ngmlr.flt3.bam.bai"
	output:
		"results/{sample}/targeted/roi/{sample}.ngmlr.sniffles.flt3.vcf"
	threads:
		lambda wildcards: round(int(config["threads"]) / 2 - 0.5)
	log:
		"logs/{sample}/{sample}.ngmlr.sniffles_flt3.log"
	benchmark:
		"benchmarks/{sample}/{sample}.ngmlr.sniffles_flt3.benchmark.tsv"
	conda:
		"../envs/sniffles.yaml"
	shell:
		"sniffles -m {input.bam} -v {output} -s 2 -t {threads} -l 5 --report-seq > {log} 2>&1" ## --genotype

rule varscan_snvs:
	input:
		ref="resources/GRCh38_full_analysis_set_plus_decoy_hla.masked.fa",
		bam="alignments/{sample}/targeted/roi/{sample}.minimap2.targeted.srt.bam",
		bai="alignments/{sample}/targeted/roi/{sample}.minimap2.targeted.srt.bam.bai",
		intervals="resources/ramd_intervals.bed",
	output:
		vcf="results/{sample}/targeted/roi/{sample}.varscan.snvs.vcf.gz",
		tbi="results/{sample}/targeted/roi/{sample}.varscan.snvs.vcf.gz.tbi"
	threads:
		1
	log:
		"logs/{sample}/{sample}.varscan_snvs.log"
	benchmark:
		"benchmarks/{sample}/{sample}.varscan_snvs.tsv"
	conda:
		"../envs/varscan.yaml"
	shell:
		"""
		samtools mpileup -d 1000000 -l {input.intervals} -f {input.ref} {input.bam} | \
		varscan mpileup2snp --output-vcf --min-avg-qual 0 --p-value 0.1 --min-var-freq 0.01 | \
		bgzip -c > {output.vcf} 2>{log} && \
		tabix {output.vcf} 2>>{log}
		"""

rule varscan_indels:
	input:
		ref="resources/GRCh38_full_analysis_set_plus_decoy_hla.masked.fa",
		bam="alignments/{sample}/targeted/roi/{sample}.minimap2.targeted.srt.bam",
		bai="alignments/{sample}/targeted/roi/{sample}.minimap2.targeted.srt.bam.bai",
		intervals="resources/ramd_intervals.bed",
	output:
		vcf="results/{sample}/targeted/roi/{sample}.varscan.indels.vcf.gz",
		tbi="results/{sample}/targeted/roi/{sample}.varscan.indels.vcf.gz.tbi"
	threads:
		1
	log:
		"logs/{sample}/{sample}.varscan_indels.log"
	benchmark:
		"benchmarks/{sample}/{sample}.varscan_indels.tsv"
	conda:
		"../envs/varscan.yaml"
	shell:
		"""
		samtools mpileup -d 1000000 -l {input.intervals} -f {input.ref} {input.bam} | \
		varscan mpileup2indel --output-vcf --min-avg-qual 0 --p-value 0.1 --min-var-freq 0.01 | \
		bgzip -c > {output.vcf} 2>{log} && \
		tabix -p vcf {output.vcf} 2>>{log}
		"""

rule varscan_concat:
	input:
		vcf_snv="results/{sample}/targeted/roi/{sample}.varscan.snvs.vcf.gz",
		vcf_indel="results/{sample}/targeted/roi/{sample}.varscan.indels.vcf.gz",
		tbi_snv="results/{sample}/targeted/roi/{sample}.varscan.snvs.vcf.gz.tbi",
		tbi_indel="results/{sample}/targeted/roi/{sample}.varscan.indels.vcf.gz.tbi",
		ref="resources/GRCh38_full_analysis_set_plus_decoy_hla.masked.fa"
	output:
		temp("results/{sample}/targeted/roi/{sample}.varscan.concat.tmp.vcf.gz")
	threads:
		1
	log:
		"logs/{sample}/{sample}.varscan_concat.log"
	benchmark:
		"benchmarks/{sample}/{sample}.varscan_concat.tsv"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"bcftools concat -a -Ov {input.vcf_snv} {input.vcf_indel} 2> {log} | "
		"bcftools norm -Oz -m - -f {input.ref} -Oz -o {output} 2>> {log}"

rule varscan_reheader:
	input:
		"results/{sample}/targeted/roi/{sample}.varscan.concat.tmp.vcf.gz"
	output:
		vcf="results/{sample}/targeted/roi/{sample}.varscan.vcf.gz",
		tbi="results/{sample}/targeted/roi/{sample}.varscan.vcf.gz.tbi",
		txt=temp("results/{sample}/targeted/roi/{sample}.rh.txt")
	threads:
		1
	log:
		"logs/{sample}/{sample}.varscan_reheader.log"
	benchmark:
		"benchmarks/{sample}/{sample}.varscan_reheader.tsv"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"echo {wildcards.sample} > {output.txt} && "
		"bcftools reheader -s {output.txt} -o {output.vcf} {input} 2> {log} && "
		"tabix -p vcf {output.vcf}"
