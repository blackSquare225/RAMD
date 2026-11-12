rule score_variants_roi:
	input:
		script="workflow/scripts/score_variants.R",
		tsv="results/{sample}/targeted/roi/{sample}.varscan.vep.tsv",
	output:
		tsv=temp("results/{sample}/targeted/roi/{sample}.varscan.vep.scored.tmp.tsv")
	threads: 
		1
	log:
		"logs/{sample}/{sample}.score_variants.log"
	benchmark:
		"benchmarks/{sample}/{sample}.score_variants.tsv"
	conda:
		"../envs/score_variants.yaml" 
	shell:
		"""
		Rscript {input.script} {input.tsv} {output.tsv} 2> {log}
		"""

rule fill_with_NA_roi:
	input:
		"results/{sample}/targeted/roi/{sample}.varscan.vep.scored.tmp.tsv"
	output:
		"results/{sample}/targeted/roi/{sample}.varscan.vep.scored.tsv"	
	threads: 
		1
	log:
		"logs/{sample}/{sample}.fill_with_NA.log"
	benchmark:
		"benchmarks/{sample}/{sample}.fill_with_NA.tsv"
	shell:
		"""
		awk 'BEGIN{{FS=OFS="\t"}} {{for(i=1;i<=NF;i++) if($i==".") $i="NA"}}1' {input} > {output} 2> {log}
		"""

rule filter_variants_roi:
	input:
		script="workflow/scripts/filter_variants.R",
		tsv="results/{sample}/targeted/roi/{sample}.varscan.vep.scored.tsv",
		cosmic="resources/COSMIC_db.tsv",
		fp="resources/FP_db.tsv",
		benign="resources/BENIGN_db.tsv",
		vus="resources/VUS_db.tsv",
		diagnostic="resources/DIAGNOSTIC_db.tsv",
		public_data="resources/PUBLIC_AML_DATA.tsv"
	output:
		tsv="results/{sample}/targeted/roi/{sample}.varscan.vep.filtered.tsv"
	threads: 
		1
	log:
		"logs/{sample}/{sample}.filter_variants.log"
	benchmark:
		"benchmarks/{sample}/{sample}.filter_variants.tsv"
	conda:
		"../envs/filter_variants.yaml"
	shell:
		"""
		Rscript {input.script} {input.tsv} {input.cosmic} {input.fp} {input.benign} {input.vus} {input.diagnostic} {input.public_data} {output.tsv} 2> {log}
		"""

rule create_final_report_targeted: 
	input:
		script="workflow/scripts/create_final_report_targeted.R",
		fusion_500kbp="results/{sample}/targeted/fusion/500kbp/jaffa_results.csv",
		fusion_1Mbp="results/{sample}/targeted/fusion/1Mbp/jaffa_results.csv",
		snpIndel="results/{sample}/targeted/roi/{sample}.varscan.vep.filtered.tsv",
		flt3="results/{sample}/targeted/roi/{sample}.ngmlr.sniffles.vep.tsv",
	output:
		"results/{sample}/targeted/{sample}.report.targeted.xlsx"
	threads: 1
	log:
		"logs/{sample}/{sample}.create_final_report_targeted.log"
	benchmark:
		"benchmarks/{sample}/{sample}.create_final_report_targeted.tsv"
	conda:	
		"../envs/create_final_report.yaml"
	shell:
		"""
		Rscript {input.script} {input.fusion_500kbp} {input.fusion_1Mbp} {input.snpIndel} {input.flt3} {output} 2> {log}
		"""