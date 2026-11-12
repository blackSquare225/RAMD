rule sniffles2_wgs:
	input:
		bam="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam",
		bai="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam.bai",
		ref="resources/GRCh38_full_analysis_set_plus_decoy_hla.fa"
	output:
		vcf="results/{sample}/wgs/{sample}.sniffles2.vcf.gz",
		snf="results/{sample}/wgs/{sample}.sniffles2.snf"
	params:
		"results/{sample}/wgs/{sample}.sniffles2.vcf"
	threads: 
		lambda wildcards: int(config["threads"])
	log:
		"logs/{sample}/{sample}.sniffles2_wgs.log"
	benchmark:
		"benchmarks/{sample}/{sample}.sniffles2_wgs.tsv"
	conda:
		"../envs/sniffles2.7.yaml"
	shell:
		"""
		sniffles \
			--input {input.bam} \
			--vcf {params} \
			--snf {output.snf} \
			--output-rnames \
			--reference {input.ref} \
			--minsupport 1 \
			--minsvlen 50 \
			--allow-overwrite \
			--qc-bnd-filter-strand False \
			--mapq 20 \
			--threads {threads} 2>{log} && \
		bgzip {params}  
		""" 

rule sniffles2_parse_svtype:
	input:
		vcf="results/{sample}/wgs/{sample}.sniffles2.vcf.gz",
		excl="resources/GRCh38_full_analysis_set_plus_decoy_hla.exclude.txt",
	output:
		vcf_nd="results/{sample}/wgs/{sample}.sniffles2.nd.vcf",
		vcf_bnd="results/{sample}/wgs/{sample}.sniffles2.BND.vcf",
		vcf_inv="results/{sample}/wgs/{sample}.sniffles2.INV.vcf"
	threads:
		1
	log:
		"logs/{sample}/{sample}.sniffles2_parse_svtype.log"
	benchmark:
		"benchmarks/{sample}/{sample}.sniffles2_parse_svtype.tsv"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		zgrep -v -F -f {input.excl} {input.vcf} > {output.vcf_nd} 2> {log} && \
		bcftools view -i'INFO/SVTYPE == "BND"' -Ov -o {output.vcf_bnd} {output.vcf_nd} 2>> {log} && \
		bcftools view -i'INFO/SVTYPE == "INV" && INFO/SVLEN > 500000' -Ov -o {output.vcf_inv} {output.vcf_nd} 2>> {log}
		"""

rule sniffles2_svafotate:
	input:
		vcf="results/{sample}/wgs/{sample}.sniffles2.{svtype}.vcf",
		bed="resources/SVAFotate_core_SV_popAFs.GRCh38.v4.1.bed.gz",
	output:
		vcf="results/{sample}/wgs/{sample}.sniffles2.{svtype}.svafotate.vcf",
	threads:
		config["threads"]
	conda:
		"../envs/svafotate.yaml"
	log:
		"logs/{sample}/{sample}.sniffles2_svafotate.{svtype}.log"
	shell:
		"""
		svafotate annotate --minf 0.8 --ann best pops -c 0 -l 1000000 --cpu {threads} --bed {input.bed} --vcf {input.vcf} --out {output} 2> {log}
		"""

rule sniffles2_remove_pop_vars:
	input:
		"results/{sample}/wgs/{sample}.sniffles2.{svtype}.svafotate.vcf"
	output:
		"results/{sample}/wgs/{sample}.sniffles2.{svtype}.filtPopVars.vcf"
	threads:
		1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}/{sample}.sniffles2_remove_pop_vars.{svtype}.log"
	shell:
		"bcftools view -i'Max_AF <= 0.01' -Ov -o {output} {input} 2> {log}"

rule sniffles2_query_info_field:
	input:
		vcf_bnd="results/{sample}/wgs/{sample}.sniffles2.BND.filtPopVars.vcf",
		vcf_inv="results/{sample}/wgs/{sample}.sniffles2.INV.filtPopVars.vcf"
	output:
		tsv_bnd="results/{sample}/wgs/{sample}.sniffles2.BND.filtPopVars.INFO.tsv",
		tsv_inv="results/{sample}/wgs/{sample}.sniffles2.INV.filtPopVars.INFO.tsv"
	threads:
		1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}/{sample}.sniffles2_query_info_field.log"
	shell:
		"""
		cat <(echo -e "name\tQUAL\tPRECISE\tSVTYPE\tRNAMES\tCOVERAGE\tSTRAND\tSTDEV_POS\tVAF\tDR\tDV\tGT") \
		<(bcftools query -f'%ID\t%QUAL\t%PRECISE\t%SVTYPE\t%RNAMES\t%COVERAGE\t%STRAND\t%STDEV_POS\t%VAF\t[%DR\t%DV\t%GT]' {input.vcf_bnd}) > {output.tsv_bnd} 2> {log} && \
		cat <(echo -e "name\tQUAL\tPRECISE\tSVLEN\tSVTYPE\tRNAMES\tCOVERAGE\tSTDEV_LEN\tSTDEV_POS\tVAF\tDR\tDV\tGT") \
		<(bcftools query -f'%ID\t%QUAL\t%PRECISE\t%SVLEN\t%SVTYPE\t%RNAMES\t%COVERAGE\t%STDEV_LEN\t%STDEV_POS\t%VAF\t[%DR\t%DV\t%GT]' {input.vcf_inv}) > {output.tsv_inv} 2>> {log}
		"""

rule sniffles2_vcf_to_bedpe:
	input:
		"results/{sample}/wgs/{sample}.sniffles2.{svtype}.filtPopVars.vcf"
	output:
		"results/{sample}/wgs/{sample}.sniffles2.{svtype}.filtPopVars.bedpe"
	threads:
		1
	conda:
		"../envs/svtools.yaml"
	log:
		"logs/{sample}/{sample}.sniffles2_vcf_to_bedpe.{svtype}.log"
	shell:
		"svtools vcftobedpe --input {input} --output {output} 2> {log}"


rule sniffles2_bedpe_to_brk_bnd:
	input:
		bedpe_bnd="results/{sample}/wgs/{sample}.sniffles2.BND.filtPopVars.bedpe",
	output:
		brk1_bnd="results/{sample}/wgs/{sample}.sniffles2.BND.filtPopVars.brk1.bed",
		brk2_bnd="results/{sample}/wgs/{sample}.sniffles2.BND.filtPopVars.brk2.bed"
	threads:
		1
	log:
		"logs/{sample}/{sample}.sniffles2_bedpe_to_brk_bnd.log"
	shell:
		r"""
		if [ ! -s {input.bedpe_bnd} ]; then
			echo "Input is empty or missing → touching output"
			touch {output.brk1_bnd} {output.brk2_bnd}
		else
			echo "Input is not empty → processing..."
			grep -v '^#' {input.bedpe_bnd} | awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $2+2, $7}}' > {output.brk1_bnd} 2> {log} && \
			grep -v '^#' {input.bedpe_bnd} | awk 'BEGIN{{OFS="\t"}} {{print $4, $5, $5+2, $7}}' > {output.brk2_bnd} 2>> {log}        
		fi
		"""

rule sniffles2_bedpe_to_brk_inv:
	input:
		bedpe_inv="results/{sample}/wgs/{sample}.sniffles2.INV.filtPopVars.bedpe"
	output:
		brk1_inv="results/{sample}/wgs/{sample}.sniffles2.INV.filtPopVars.brk1.bed",
		brk2_inv="results/{sample}/wgs/{sample}.sniffles2.INV.filtPopVars.brk2.bed"
	threads:
		1
	log:
		"logs/{sample}/{sample}.sniffles2_bedpe_to_brk_inv.log"
	shell:
		r"""
		if [ ! -s {input.bedpe_inv} ]; then
			echo "Input is empty or missing → touching output"
			touch {output.brk1_inv} {output.brk2_inv}
		else
			echo "Input is not empty → processing..."
			grep -v '^#' {input.bedpe_inv} | awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $2+2, $7}}' > {output.brk1_inv} 2> {log} && \
			grep -v '^#' {input.bedpe_inv} | awk 'BEGIN{{OFS="\t"}} {{print $4, $5, $5+2, $7}}' > {output.brk2_inv} 2>> {log}
		fi
		"""

rule sniffles2_annotate_brk:
	input:
		"workflow/scripts/annotate_brk.R",
		"resources/gencode.v49.annotation.bed",
		"results/{sample}/wgs/{sample}.sniffles2.{svtype}.filtPopVars.{brk}.bed"
	output:
		"results/{sample}/wgs/{sample}.sniffles2.{svtype}.filtPopVars.{brk}.annot.bed"
	threads:
		1
	conda:
		"../envs/annotate_brk.yaml"
	log:
		"logs/{sample}/{sample}.sniffles2_annotate_{svtype}.{brk}.log"
	shell:
		"Rscript {input} {output} 2> {log}"

rule sniffles2_merge_annot_brk: 
	input:
		"workflow/scripts/merge_annot_brk.R",
		"results/{sample}/wgs/{sample}.sniffles2.{svtype}.filtPopVars.brk1.annot.bed",
		"results/{sample}/wgs/{sample}.sniffles2.{svtype}.filtPopVars.brk2.annot.bed",
		"results/{sample}/wgs/{sample}.sniffles2.{svtype}.filtPopVars.INFO.tsv"
	output:
		"results/{sample}/wgs/{sample}.sniffles2.{svtype}.filtPopVars.merged.annot.tsv"
	threads:
		1
	conda:
		"../envs/merge_annot_brk.yaml"
	log:
		"logs/{sample}/{sample}.sniffles2_merge_annot_brk.{svtype}.log"
	shell:
		"Rscript {input} {output} 2> {log}"

rule chr3_rearrangment:
	input:
		script="workflow/scripts/gata2_mecom.py",
		bam="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam",
		bed="resources/hg38.ncbiRefSeq.sorted.bed"
	output:
		"results/{sample}/wgs/{sample}.chr3rearrangments.annot.tsv"
	threads: 
		1 
	conda:
		"../envs/gata2_mecom.yaml"
	log:
		"logs/{sample}/{sample}.chr3_rearrangment.log"
	benchmark:
		"benchmarks/{sample}/{sample}.chr3_rearrangment.tsv"
	shell:
		"python {input.script} -b {input.bam} -a {input.bed} -r chr3:169000000-169088157 -o {output}"

rule create_final_report_sv:
	input:
		Rscript="workflow/scripts/create_final_report_sv.R",
		vcf_inv="results/{sample}/wgs/{sample}.sniffles2.INV.filtPopVars.merged.annot.tsv",
		vcf_bnd="results/{sample}/wgs/{sample}.sniffles2.BND.filtPopVars.merged.annot.tsv",
		chr3_rearr="results/{sample}/wgs/{sample}.chr3rearrangments.annot.tsv"
	output:
		"results/{sample}/wgs/{sample}.report.sv.xlsx"
	threads:
		1
	log:
		"logs/{sample}/{sample}.create_final_report_sv.log"
	shell:
		"Rscript {input.Rscript} {input.vcf_inv} {input.vcf_bnd} {input.chr3_rearr} {output} 2> {log}"