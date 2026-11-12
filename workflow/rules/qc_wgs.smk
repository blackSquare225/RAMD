rule NanoPlot:
	input:
		bam="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam",
		bai="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam.bai"
	output:
		"qc/{sample}/wgs/{sample}_NanoStats.txt"
	log:
		"logs/{sample}/{sample}.NanoPlot.log"
	threads:
		config["threads"]
	conda:
		"../envs/nanoplot.yaml"
	params:
		outdir="qc/{sample}/wgs/",
		prefix="{sample}_"
	shell:
		"NanoPlot -t {threads} --bam {input.bam} --N50 --prefix {params.prefix} -o {params.outdir}"

rule mosdepth_minimap2_wgs:
	input:
		bam="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam",
		bai="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam.bai"
	output:
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.mosdepth.global.dist.txt",
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.mosdepth.region.dist.txt",
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.mosdepth.summary.txt",
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.per-base.bed.gz.csi",
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.per-base.bed.gz",
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.regions.bed.gz.csi",
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.regions.bed.gz"
	threads: 
		lambda wildcards: round(int(config["threads"]) / 2)
	log:
		"logs/{sample}/{sample}.mosdepth_minimap2_wgs.log"
	benchmark:
		"benchmarks/{sample}/{sample}.mosdepth_minimap2_wgs.tsv"
	conda: 
		"../envs/mosdepth.yaml"
	params:
		prefix="qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp",
	shell:
		"mosdepth --fast-mode -t {threads} --by 1000 -Q20 {params.prefix} {input.bam} 2>{log}"

rule multiqc:
	input:
		"qc/{sample}/wgs/{sample}_NanoStats.txt",
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.mosdepth.global.dist.txt",
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.mosdepth.region.dist.txt",
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.mosdepth.summary.txt",
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.per-base.bed.gz.csi",
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.per-base.bed.gz",
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.regions.bed.gz.csi",
		"qc/{sample}/wgs/{sample}.minimap2_wgs_1000bp.regions.bed.gz"
	output:
		"qc/{sample}/wgs/{sample}.multiqc.html"
	params:
		"qc/{sample}/wgs/"
	threads:
		1
	conda:
		"../envs/multiqc.yaml"
	log:
		"log/{sample}/multiqc.log"
	shell:
		"multiqc --force --filename {output} {params} 2> {log}"