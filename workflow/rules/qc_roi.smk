rule mosdepth_minimap2_roi:
	input:
		bam="alignments/{sample}/targeted/roi/{sample}.minimap2.targeted.srt.bam",
		bai="alignments/{sample}/targeted/roi/{sample}.minimap2.targeted.srt.bam.bai",
		intervals="resources/ramd_intervals.bed"
	output:
		"qc/{sample}/targeted/roi/{sample}.minimap2_targeted.mosdepth.global.dist.txt",
		"qc/{sample}/targeted/roi/{sample}.minimap2_targeted.mosdepth.region.dist.txt",
		"qc/{sample}/targeted/roi/{sample}.minimap2_targeted.mosdepth.summary.txt",
		"qc/{sample}/targeted/roi/{sample}.minimap2_targeted.per-base.bed.gz.csi",
		"qc/{sample}/targeted/roi/{sample}.minimap2_targeted.per-base.bed.gz",
		"qc/{sample}/targeted/roi/{sample}.minimap2_targeted.regions.bed.gz.csi",
		"qc/{sample}/targeted/roi/{sample}.minimap2_targeted.regions.bed.gz"
	threads: 
		lambda wildcards: round(int(config["threads"]) / 2)
	log:
		"logs/{sample}/{sample}.mosdepth_minimap2_targeted.log"
	benchmark:
		"benchmarks/{sample}/{sample}.mosdepth_minimap2_targeted.tsv"
	conda: 
		"../envs/mosdepth.yaml"
	params:
		prefix="qc/{sample}/targeted/roi/{sample}.minimap2_targeted",
	shell:
		"mosdepth -t {threads} --by {input.intervals} {params.prefix} {input.bam} 2>{log}"

rule coverage_stats_roi:
	input:
		bed="qc/{sample}/targeted/roi/{sample}.minimap2_targeted.regions.bed.gz",
		script="workflow/scripts/coverage_stats.R"
	output:
		plot="qc/{sample}/targeted/roi/{sample}_cov_plot.pdf",
		tsv="qc/{sample}/targeted/roi/{sample}_cov_stats.tsv",
		csv="qc/{sample}/targeted/roi/{sample}_cov_stats.csv"
	threads: 
		1
	log:
		"logs/{sample}/{sample}.coverage_stats.log"
	benchmark:
		"benchmarks/{sample}/{sample}.coverage_stats.tsv"
	conda:
		"../envs/coverage_stats_roi.yaml"
	shell:
		"Rscript {input.script} {input.bed} {output.plot} {output.tsv} {output.csv}"
