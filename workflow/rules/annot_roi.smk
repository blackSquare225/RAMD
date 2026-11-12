rule sniffles_vep:
	input:
		ref="resources/GRCh38_full_analysis_set_plus_decoy_hla.masked.fa",
		vcf="results/{sample}/targeted/roi/{sample}.ngmlr.sniffles.flt3.vcf"
	output:
		calls="results/{sample}/targeted/roi/{sample}.ngmlr.sniffles.flt3.vep.vcf.gz",
		html="results/{sample}/targeted/roi/{sample}.ngmlr.sniffles.flt3.html"
	threads: 
		lambda wildcards: round(int(config["threads"]) / 2 - 0.5)
	log:
		"logs/{sample}/{sample}.ngmlr.sniffles.log"
	benchmark:
		"benchmarks/{sample}/{sample}.sniffles_vep.tsv"
	conda:
		"../envs/vep.yaml"
	shell:
		"vep -i {input.vcf} -o {output.calls} --compress_output bgzip --fork {threads} --offline --species homo_sapiens --stats_file {output.html} --assembly GRCh38 --cache --cache_version 113 --merged --dir_cache resources/vep/cache --fasta {input.ref} --format vcf --symbol --no_intergenic --merged --polyphen b --sift b --cache --pick --pick_allele --vcf 2>{log}"

rule sniffles_vcf_to_tsv:
	input:
		vcf="results/{sample}/targeted/roi/{sample}.ngmlr.sniffles.flt3.vep.vcf.gz",
		header="resources/header_sniffles.tsv"
	output:
		"results/{sample}/targeted/roi/{sample}.ngmlr.sniffles.vep.tsv"
	threads: 
		1
	log:
		"logs/{sample}/{sample}.sniffles_vcf_to_tsv.log"
	benchmark:
		"benchmarks/{sample}/{sample}.sniffles_vcf_to_tsv.tsv"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		cat {input.header} <(bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%REF\t%ALT\t%PRECISE\t%SVTYPE\t%SVLEN\t%CSQ\t[%GT\t%DR\t%DV]\t%AF\n" -d -A tab) > {output} 2>{log}
		"""

rule varscan_vep:
	input:
		vcf="results/{sample}/targeted/roi/{sample}.varscan.vcf.gz",
		ref="resources/GRCh38_full_analysis_set_plus_decoy_hla.masked.fa",
	output:
		calls="results/{sample}/targeted/roi/{sample}.varscan.vep.vcf.gz",
		html="results/{sample}/targeted/roi/{sample}.varscan.vep.html"
	threads: 
		lambda wildcards: round(int(config["threads"]) / 2)
	log:
		"logs/{sample}/{sample}.varscan_vep.log"
	benchmark:
		"benchmarks/{sample}/{sample}.varscan_vep.tsv"
	conda:
		"../envs/vep.yaml"
	shell:
		"vep -i {input.vcf} -o {output.calls} --compress_output bgzip --fork {threads} --everything --offline --species homo_sapiens --stats_file {output.html} --assembly GRCh38 --cache --cache_version 113 --merged --dir_cache resources/vep/cache --fasta {input.ref} --format vcf --symbol --no_intergenic --merged --polyphen b --sift b --cache --pick --pick_allele --vcf --plugin dbNSFP,resources/dbNSFP4.5a_grch38.gz,SIFT_converted_rankscore,SIFT_pred,SIFT4G_converted_rankscore,SIFT4G_pred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,Polyphen2_HVAR_score,Polyphen2_HVAR_pred,LRT_score,LRT_pred,MutationTaster_score,MutationTaster_pred,MutationAssessor_rankscore,MutationAssessor_pred,FATHMM_converted_rankscore,FATHMM_pred,PROVEAN_converted_rankscore,PROVEAN_pred,VEST4_rankscore,MetaSVM_rankscore,MetaSVM_pred,MetaLR_rankscore,MetaLR_pred,MetaRNN_rankscore,MetaRNN_pred,M-CAP_score,M-CAP_pred,DEOGEN2_score,DEOGEN2_pred,ClinPred_score,ClinPred_pred,LIST-S2_score,LIST-S2_pred,CADD_raw_rankscore,CADD_phred,fathmm-MKL_coding_score,fathmm-MKL_coding_pred,fathmm-XF_coding_score,fathmm-XF_coding_pred,GERP++_RS_rankscore,ESM1b_rankscore,ESM1b_pred 2>{log}"

rule varscan_vcf_to_tsv:
	input:
		vcf="results/{sample}/targeted/roi/{sample}.varscan.vep.vcf.gz",
		header="resources/header.varscan.tsv"
	output:
		"results/{sample}/targeted/roi/{sample}.varscan.vep.tsv"
	threads: 
		1
	log:
		"logs/{sample}/{sample}.varscan_vcf_to_tsv.log"
	conda:
		"../envs/bcftools.yaml"
	benchmark:
		"benchmarks/{sample}/{sample}.varscan_vcf_to_tsv.tsv"
	shell:
		"""
		cat {input.header} <(bcftools +split-vep {input.vcf} -f "%CHROM\t%POS\t%REF\t%ALT\t%CSQ\t[%GT\t%RD\t%AD\t%FREQ\t%RBQ\t%ABQ\t%PVAL\t%SAMPLE]\n" -d -A tab) > {output} 2>{log}
		"""
