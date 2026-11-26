rule readCounter_50kbp:
	input:
		bam="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam",
		bai="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam.bai"
	output:
		"results/{sample}/wgs/{sample}.50kbp.wig"
	threads:
		1
	log:
		"logs/{sample}/{sample}.readCounter_50kbp.log"
	benchmark:
		"benchmarks/{sample}/{sample}.readCounter_50kbp.tsv"
	conda:
		"../envs/hmmcopy.yaml"
	shell:
		"""
		readCounter \
			--window 50000 \
			--quality 20 \
			--chromosome $(samtools view -H {input.bam} | tail -n+2 | head -24 | cut -f2 | cut -f2 -d':' | tr '\n' ',' | sed 's/,$//') \
			{input.bam} > {output}
		"""

rule ichorCNA_50kbp:
	input:
		ichorCNA="resources/ichorCNA/scripts/runIchorCNA.R",
		wig="results/{sample}/wgs/{sample}.50kbp.wig",
		gcWig="resources/ichorCNA/inst/extdata/gc_hg38_50kb.wig",
		mapWig="resources/ichorCNA/inst/extdata/map_hg38_50kb.wig",
		centromere="resources/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt",
	output:
		pdf="results/{sample}/wgs/{sample}/{sample}_tpdf.pdf",
		seg="results/{sample}/wgs/{sample}.seg"
	params:
		"results/{sample}/wgs/"
	threads:
		1
	conda:
		"../envs/ichorCNA.yaml"
	shell:
		"""
        Rscript {input.ichorCNA} \
            --genomeBuild hg38 \
            --libdir resources/ichorCNA/ \
            --id {wildcards.sample} \
            --WIG {input.wig} \
            --ploidy "c(2,3)" \
            --normal "c(0.3,0.4,0.5,0.6,0.7,0.8,0.9)" \
            --maxCN 5 \
            --gcWig {input.gcWig} \
            --mapWig {input.mapWig} \
            --centromere {input.centromere} \
            --includeHOMD False \
            --chrs "paste0('chr', c(1:22, 'X', 'Y'))" \
            --chrTrain "paste0('chr', c(1:22))" \
            --estimateNormal True \
            --estimatePloidy True \
            --estimateScPrevalence True \
            --scStates "c(1,3)" \
            --outDir {params} 
        """

rule QDNAseq_50kbp:
	input:
		Rscript="workflow/scripts/QDNAseq.R",
		bam="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam",
		bai="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam.bai"
	output:
		bed="results/{sample}/wgs/{sample}.QDNAseq.50kbp.bed",
		seg="results/{sample}/wgs/{sample}.QDNAseq.50kbp.seg"
	params:
		"results/{sample}/wgs/"
	threads: 1
	log:
		"logs/{sample}/{sample}.minimap2.QDNAseq_50kbp.log"
	benchmark:
		"benchmarks/{sample}/{sample}.minimap2.QDNAseq_50kbp.tsv"
	#conda:
	#	"../envs/QDNAseq.yaml"
	shell:
		"Rscript {input.Rscript} 50 {input.bam} {params}"


rule parse_QDNAseq:
	input:
		Rscript="workflow/scripts/parse_QDNAseq.R",
		seg="results/{sample}/wgs/{sample}.QDNAseq.50kbp.seg",
	output:
		"results/{sample}/wgs/{sample}.QDNAseq.{window}.filt.bed"
	threads:
		1
	log:
		"logs/{sample}/{sample}.parse_QDNAseq.{window}.log"
	benchmark:
		"benchmarks/{sample}/{sample}.parse_QDNAseq.{window}.tsv"
	shell:
		"""
		Rscript {input.Rscript} {input.seg} {wildcards.sample} {output} 2> {log}
		"""

rule parse_ichorCNA:
	input:
		"results/{sample}/wgs/{sample}.seg"
	output:
		"results/{sample}/wgs/{sample}.ichorCNA.50kbp.filt.bed"
	threads:
		1
	shell:
		"""
		awk -v sample={wildcards.sample} 'BEGIN {{OFS=FS="\t"}} NR>1 && $5!="NEUT" {{print "chr"$2, $3, $4, $5, $6, $7, $8, "ichorCNA", sample}}' {input} > {output}
		"""

rule annot_ichorCNA:
	input:
		AnnotSV="resources/AnnotSV/bin/AnnotSV",
		bed="results/{sample}/wgs/{sample}.ichorCNA.50kbp.filt.bed"
	output:
		bed="results/{sample}/wgs/{sample}.ichorCNA.50kbp.annotSV.tsv"
	params:
		"results/{sample}/wgs/"
	threads:
		1
	log:
		"logs/{sample}/{sample}.annot_ichorCNA.50kbp.log"
	benchmark:
		"benchmarks/{sample}/{sample}.annot_ichorCNA.50kbp.tsv"
	shell:
		"""
		{input.AnnotSV} \
			-annotationMode full \
			-svtBEDcol 4 \
			-samplesidBEDcol 9 \
			-SVinputFile {input.bed} \
			-genomeBuild GRCh38 \
			-outputFile {output.bed} \
			-outputDir {params} \
			-overlap 60 && \
		sed -i 's/user#1/CN/;s/user#2/Bins/;s/user#3/Log2Ratio/;s/user#4/Tool/' {output.bed}
		"""


rule create_final_report_cnv:
	input:
		Rscript="workflow/scripts/create_final_report_cnv.R",
		QDNAseq_50kbp="results/{sample}/wgs/{sample}.QDNAseq.50kbp.annotSV.tsv",
		QDNAseq_100kbp="results/{sample}/wgs/{sample}.QDNAseq.100kbp.annotSV.tsv", 
		ichorCNA="results/{sample}/wgs/{sample}.ichorCNA.50kbp.annotSV.tsv",
	output:
		"results/{sample}/wgs/{sample}.report.cnv.xlsx"
	threads:
		1
	log:
		"logs/{sample}/{sample}.create_final_report_cnv.log"
	shell:
		"Rscript {input.Rscript} {input.QDNAseq_50kbp} {input.QDNAseq_100kbp} {input.ichorCNA} {output} 2> {log}"