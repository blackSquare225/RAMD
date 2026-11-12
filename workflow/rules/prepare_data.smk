rule download_gencode_annotation:
	output:
		"resources/gencode.v49.annotation.gtf"
	threads:
		1
	log:
		"logs/download_gencode_annotation.log"
	shell:
		"""
		wget --output-document=resources/gencode.v49.annotation.gtf.gz https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz 2> {log} && \
		bgzip -d resources/gencode.v49.annotation.gtf.gz		
		"""

rule parse_gtf_into_bed:
	input:
		"workflow/scripts/parse_gtf_into_bed.R",
		"resources/gencode.v49.annotation.gtf"
	output:
		"resources/gencode.v49.annotation.bed"
	threads:
		1
	conda:
		"../envs/parse_gtf_into_bed.yaml"
	log:
		"logs/parse_gtf_into_bed.log"
	shell:
		"Rscript {input} {output} 2> {log}"


#rule create_varscan_header:
#	input:
#		"results/testSample/testSample.varscan.vep.vcf.gz",
#	output:
#		"resources/header.varscan.tsv"
#	threads:
#		1
#	log:
#		"logs/testSample/testSample.create_varscan_header.log"
#	benchmark:
#		"benchmarks/testSample/testSample.create_varscan_header.tsv"
#	shell:
#		"""
#		zgrep -m1 "##INFO=<ID=CSQ" {input} | \
#		sed -n 's/.*Format: \\(.*\\)/\\1/p' | \
#		sed 's/>$//' | \
#		sed 's/"//g' | \
#		tr '|' '\\t' | \
#		awk 'OFS=FS="\\t"''{{print "CHROM", "POS", "REF", "ALT", $0, "GT", "RD", "AD", "FREQ", "RBQ", "ABQ", "PVAL", "SAMPLE"}}' > {output} 2>{log}      
#		"""

#rule create_sniffles_header:
#	input:
#		"results/testSample/testSample.ngmlr.sniffles.flt3.vep.vcf.gz",
#	output:
#		"resources/header_sniffles.tsv"
#	threads:
#		1
#	log:
#		"logs/testSample/testSample.create_sniffles_header.log"
#	benchmark:
#		"benchmarks/testSample/testSample.create_header.tsv"
#	shell:
#		"""
#		zgrep -m1 "##INFO=<ID=CSQ" {input} | \
#		sed -n 's/.*Format: \\(.*\\)/\\1/p' | \
#		sed 's/>$//' | \
#		sed 's/"//g' | \
#		tr '|' '\\t' | \
#		awk 'OFS=FS="\\t"''{{print "CHROM", "POS", "REF", "ALT", "PRECISE", "SVTYPE", "SVLEN", $0, "GT", "DR", "DV", "AF"}}' > {output} 2>{log}      
#		"""