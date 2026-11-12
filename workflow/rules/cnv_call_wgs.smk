rule QDNAseq_30kbp:
	input:
		Rscript="workflow/scripts/QDNAseq.R",
		bam="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam",
		bai="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam.bai"
	output:
		bed="results/{sample}/wgs/{sample}.QDNAseq.30kbp.bed",
		seg="results/{sample}/wgs/{sample}.QDNAseq.30kbp.seg"
	params:
		"results/{sample}/wgs/"
	threads: 
		1
	log:
		"logs/{sample}/{sample}.QDNAseq_30kbp.log"
	benchmark:
		"benchmarks/{sample}/{sample}.QDNAseq_30kbp.tsv"
	#conda:
	#	"../envs/QDNAseq.yaml"
	shell:
		"Rscript {input.Rscript} 30 {input.bam} {params}"

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


rule QDNAseq_100kbp:
	input:
		Rscript="workflow/scripts/QDNAseq.R",
		bam="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam",
		bai="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam.bai"
	output:
		bed="results/{sample}/wgs/{sample}.QDNAseq.100kbp.bed",
		seg="results/{sample}/wgs/{sample}.QDNAseq.100kbp.seg"
	params:
		"results/{sample}/wgs/"
	threads: 1
	log:
		"logs/{sample}/{sample}.minimap2.QDNAseq_100kbp.log"
	benchmark:
		"benchmarks/{sample}/{sample}.minimap2.QDNAseq_100kbp.tsv"
	#conda:
	#	"../envs/QDNAseq.yaml"
	shell:
		"Rscript {input.Rscript} 100 {input.bam} {params}"

rule QDNAseq_500kbp:
	input:
		Rscript="workflow/scripts/QDNAseq.R",
		bam="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam",
		bai="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam.bai"
	output:
		bed="results/{sample}/wgs/{sample}.QDNAseq.500kbp.bed",
		seg="results/{sample}/wgs/{sample}.QDNAseq.500kbp.seg"
	params:
		"results/{sample}/wgs/"
	threads: 1
	log:
		"logs/{sample}/{sample}.minimap2.QDNAseq_500kbp.log"
	benchmark:
		"benchmarks/{sample}/{sample}.minimap2.QDNAseq_500kbp.tsv"
	#conda:
	#	"../envs/QDNAseq.yaml"
	shell:
		"Rscript {input.Rscript} 500 {input.bam} {params}"

rule QDNAseq_1000kbp:
	input:
		Rscript="workflow/scripts/QDNAseq.R",
		bam="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam",
		bai="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam.bai"
	output:
		bed="results/{sample}/wgs/{sample}.QDNAseq.1000kbp.bed",
		seg="results/{sample}/wgs/{sample}.QDNAseq.1000kbp.seg"
	params:
		"results/{sample}/wgs/"
	threads: 1
	log:
		"logs/{sample}/{sample}.minimap2.QDNAseq_1000kbp.log"
	benchmark:
		"benchmarks/{sample}/{sample}.minimap2.QDNAseq_1000kbp.tsv"
	#conda:
	#	"../envs/QDNAseq.yaml"
	shell:
		"Rscript {input.Rscript} 1000 {input.bam} {params}"


rule parse_QDNAseq:
	input:
		Rscript="workflow/scripts/parse_QDNAseq.R",
		seg="results/{sample}/wgs/{sample}.QDNAseq.{window}.seg",
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

rule annot_QDNAseq:
	input:
		AnnotSV="resources/AnnotSV/bin/AnnotSV",
		bed="results/{sample}/wgs/{sample}.QDNAseq.{window}.filt.bed"
	output:
		bed="results/{sample}/wgs/{sample}.QDNAseq.{window}.annotSV.tsv"
	params:
		"results/{sample}/wgs/"
	threads:
		1
	log:
		"logs/{sample}/{sample}.annot_QDNAseq.{window}.log"
	benchmark:
		"benchmarks/{sample}/{sample}.annot_QDNAseq.{window}.tsv"
	shell:
		"""
		{input.AnnotSV} \
			-annotationMode full \
			-svtBEDcol 4 \
			-samplesidBEDcol 7 \
			-SVinputFile {input.bed} \
			-genomeBuild GRCh38 \
			-outputFile {output.bed} \
			-outputDir {params} \
			-overlap 60 && \
		sed -i 's/user#1/Bins/;s/user#2/Log2Ratio/;s/user#3/SampleID/' {output.bed}
		"""


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



#rule spectre_RemoveNs:
#	input:
#		ref="resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"
#	output:
#		"resources/metadata_file.mdr",
#	params:
#		outdir="resources/",
#		outfile="metadata_file.mdr"
#	threads: 
#		config["threads"]
#	conda:
#		"../envs/spectre.yaml"
#	log:
#		"logs/spectre_RemoveNs.log"
#	benchmark:
#		"benchmarks/spectre_RemoveNs.tsv"
#	shell:
#		"""
#		spectre removeNs \
#			--reference {input.ref} \
#			--output-dir {params.outdir} \
#			--output-file {params.outfile} \
#			--bin-size 1000 2> {log}
#		"""

#rule snf2json:
#	input:
#		"results/{sample}/{sample}.minimap2.{genome}.sniffles2.snf"
#	output:
#		"results/{sample}/{sample}.minimap2.{genome}.sniffles2.snfj.gz"
#	threads:
#		config["threads"]
#	log:
#		"logs/{sample}/{sample}.minimap2.{genome}.snf2json.log"
#	benchmark:
#		"benchmarks/{sample}/{sample}.minimap2.{genome}.snf2json.tsv"
#	conda:
#		"../envs/snf2json.yaml"
#	shell:
#		"snf2json {input} {output}"


rule spectre_CNVCaller:
	input:
		cov="qc/{sample}/{sample}.{aligner}.{genome}.regions.bed.gz",
		ref="resources/GRCh38_full_analysis_set_plus_decoy_hla",
		meta="resources/metadata_removeNs.mdr",
		blacklist="resources/grch38_blacklist_spectre.bed",
		snfj="results/{sample}/{sample}.{aligner}.{genome}.sniffles2.snfj.gz"
	output:
		"results/{sample}/{sample}.{aligner}.{genome}_spectre.vcf.gz.tbi"
	params:
		outdir="results/{sample}",
		name="{sample}.{aligner}.{genome}_spectre"
	threads: 
		config["threads"]
	log:
		"logs/{sample}/{sample}.{aligner}.{genome}.spectre_CNVCaller.log"
	benchmark:
		"benchmarks/{sample}/{sample}.{aligner}.{genome}.spectre_CNVCaller.tsv"
	conda:
		"../envs/spectre.yaml"
	shell:
		"""
		spectre CNVCaller \
		  --coverage {input.cov} \
		  --snfj {input.snfj} \
		  --sample-id {params.name} \
		  --output-dir {params.outdir} \
		  --reference {input.ref} \
		  --metadata {input.meta} \
		  --blacklist {input.blacklist} \
		  --threads {threads} 2>{log}
		""" 

#rule spectre_CNVCaller_cancerMode:
#	input:
#		cov="qc/{sample}/{sample}.{aligner}.{genome}.regions.bed.gz",
#		ref="resources/GRCh38_full_analysis_set_plus_decoy_hla",
#		meta="resources/metadata_removeNs.mdr",
#		blacklist="resources/grch38_blacklist_spectre.bed",
#		snfj="results/{sample}/{sample}.{aligner}.{genome}.sniffles2.snfj.gz"
#	output:
#		"results/{sample}/{sample}.{aligner}.{genome}_spectre_cancer.vcf.gz.tbi"
#	params:
#		outdir="results/{sample}",
#		name="{sample}.{aligner}.{genome}_spectre_cancer"
#	threads: 
#		config["threads"]
#	log:
#		"logs/{sample}/{sample}.{aligner}.{genome}.spectre_CNVCaller_cancerMode.log"
#	benchmark:
#		"benchmarks/{sample}/{sample}.{aligner}.{genome}.spectre_CNVCaller_cancerMode.tsv"
#	conda:
#		"../envs/spectre.yaml"
#	shell:
#		"""
#		spectre CNVCaller \
#		  --coverage {input.cov} \
#		  --snfj {input.snfj} \
#		  --sample-id {params.name} \
#		  --output-dir {params.outdir} \
#		  --reference {input.ref} \
#		  --metadata {input.meta} \
#		  --blacklist {input.blacklist} \
#		  --threads {threads}
#		  --cancer 2>{log}
#		""" 



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