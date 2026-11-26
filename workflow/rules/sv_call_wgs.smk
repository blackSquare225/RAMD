rule severus_to:
	input:
		bam="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam",
		bai="alignments/{sample}/wgs/{sample}.minimap2.wgs.srt.bam.bai",
		vntr="resources/human_GRCh38_no_alt_analysis_set.trf.bed",
		pon="resources/PoN/PoN_1000G_hg38.tsv.gz"
	output:
		"results/{sample}/wgs/{sample}/severus_out/somatic_SVs/severus_somatic.vcf"
	params:
		"results/wgs/{sample}/severus_out/"
	threads:
		config["threads"]
	conda:
		"../envs/severus.yaml"
	log:
		"logs/{sample}/severus_out.log"
	shell:
		"""
		severus \
			--target-bam {input.bam} \
			--out-dir {params} \
			--target-sample {wildcards.sample} \
			--threads {threads} \
			--vntr-bed {input.vntr} \
			--PON {input.pon} \
			--min-support 3 \
			--min-sv-size 500 \
			2> {log}
		"""

rule severus_parse_svtype:
	input:
		vcf="results/{sample}/wgs/{sample}/severus_out/somatic_SVs/severus_somatic.vcf",
		excl="resources/GRCh38_full_analysis_set_plus_decoy_hla.exclude.txt",
	output:
		vcf_nd="results/{sample}/wgs/{sample}.severus.nd.vcf",
		vcf_bnd="results/{sample}/wgs/{sample}/{sample}.severus.BND.vcf",
		vcf_inv="results/{sample}/wgs/{sample}/{sample}.severus.INV.vcf",
		vcf_dup="results/{sample}/wgs/{sample}/{sample}.severus.DUP.vcf"
	threads:
		1
	log:
		"logs/{sample}/{sample}.severus_parse_svtype.log"
	benchmark:
		"benchmarks/{sample}/{sample}.severus_parse_svtype.tsv"
	conda:
		"../envs/bcftools.yaml"
	shell:
		"""
		grep -v -F -f {input.excl} {input.vcf} > {output.vcf_nd} 2> {log} && \
		bcftools view -i 'QUAL >= 5' {output.vcf_nd} | bcftools view -i'INFO/SVTYPE == "BND"' -Ov -o {output.vcf_bnd} - 2>> {log} && \
		bcftools view -i 'QUAL >= 5' {output.vcf_nd} | bcftools view -i'INFO/SVTYPE == "INV"' -Ov -o {output.vcf_inv} - 2>> {log} && \
		bcftools view -i 'QUAL >= 5' {output.vcf_nd} | bcftools view -i'INFO/SVTYPE == "DUP"' -Ov -o {output.vcf_dup} - 2>> {log}
		"""

rule severus_query_info_field:
	input:
		vcf_bnd="results/{sample}/wgs/{sample}.severus.BND.vcf",
		vcf_inv="results/{sample}/wgs/{sample}.severus.INV.vcf",
		vcf_dup="results/{sample}/wgs/{sample}.severus.DUP.vcf"
	output:
		tsv_bnd="results/{sample}/wgs/{sample}.severus.BND.INFO.tsv",
		tsv_inv="results/{sample}/wgs/{sample}.severus.INV.INFO.tsv",
		tsv_dup="results/{sample}/wgs/{sample}.severus.DUP.INFO.tsv"
	threads:
		1
	conda:
		"../envs/bcftools.yaml"
	log:
		"logs/{sample}/{sample}.severus_query_info_field.log"
	shell:
		"""
		cat <(echo -e "name\tQUAL\tSVTYPE\tSTRANDS\tDR\tDV\tGT") \
			<(bcftools query -f'%ID\t%QUAL\t%SVTYPE\t%STRANDS\t[%DR\t%DV\t%GT]\n' {input.vcf_bnd}) > {output.tsv_bnd} 2> {log} && \
		cat <(echo -e "name\tQUAL\tSVTYPE\tSTRANDS\tDR\tDV\tGT") \
			<(bcftools query -f'%ID\t%QUAL\t%SVTYPE\t%STRANDS\t[%DR\t%DV\t%GT]\n' {input.vcf_inv}) > {output.tsv_inv} 2>> {log} && \
		cat <(echo -e "name\tQUAL\tSVTYPE\tSTRANDS\tDR\tDV\tGT") \
			<(bcftools query -f'%ID\t%QUAL\t%SVTYPE\t%STRANDS\t[%DR\t%DV\t%GT]\n' {input.vcf_dup}) > {output.tsv_dup} 2>> {log}
		"""

rule severus_vcf_to_bedpe:
	input:
		"results/{sample}/wgs/{sample}.severus.{SVTYPE}.vcf"
	output:
		"results/{sample}/wgs/{sample}.severus.{SVTYPE}.bedpe"
	threads:
		1
	conda:
		"../envs/svtools.yaml"
	log:
		"logs/{sample}/{sample}.severus_vcf_to_bedpe.{SVTYPE}.log"
	shell:
		"svtools vcftobedpe --input {input} --output {output} 2> {log}"

rule severus_bedpe_to_brk_bnd:
	input:
		"results/{sample}/wgs/{sample}.severus.{SVTYPE}.bedpe",
	output:
		brk1="results/{sample}/wgs/{sample}.severus.{SVTYPE}.brk1.bed",
		brk2="results/{sample}/wgs/{sample}.severus.{SVTYPE}.brk2.bed"
	threads:
		1
	log:
		"logs/{sample}/{sample}.severus_bedpe_to_brk_bnd.{SVTYPE}.log"
	shell:
		r"""
		if [ ! -s {input} ]; then
			echo "Input is empty or missing → touching output"
			touch {output.brk1} {output.brk2}
		else
			echo "Input is not empty → processing..."
			grep -v '^#' {input} | awk 'BEGIN{{OFS="\t"}} {{print $1, $2, $2+2, $7}}' > {output.brk1} 2> {log} && \
			grep -v '^#' {input} | awk 'BEGIN{{OFS="\t"}} {{print $4, $5, $5+2, $7}}' > {output.brk2} 2>> {log}        
		fi
		"""

rule severus_annotate_brk:
	input:
		"workflow/scripts/annotate_brk.R",
		"results/svcallers_comparison/{sample}/{sample}.severus.{SVTYPE}.{brk}.bed",
		"resources/gencode_v49.bed",
	output:
		"results/{sample}/wgs/{sample}.severus.{SVTYPE}.{brk}.annot.bed"
	threads:
		1
	conda:
		"../envs/annotate_brk.yaml"
	log:
		"logs/{sample}/{sample}.cutesv_annotate_brk.{SVTYPE}.{brk}.log"
	shell:
		"Rscript {input} {output} 2> {log}"


rule severus_merge_annot_brk: 
	input:
		"workflow/scripts/merge_annot_brk.R",
		"results/{sample}/wgs/{sample}.severus.{SVTYPE}.brk1.annot.bed",
		"results/{sample}/wgs/{sample}.severus.{SVTYPE}.brk2.annot.bed",
		"results/{sample}/wgs/{sample}.severus.{SVTYPE}.INFO.tsv"
	output:
		"results/{sample}/wgs/{sample}.severus.{SVTYPE}.merged.annot.tsv"
	threads:
		1
	conda:
		"../envs/merge_annot_brk.yaml"
	log:
		"logs/{sample}/{sample}.severus_merge_annot_brk.{SVTYPE}.log"
	shell:
		"Rscript {input} {output} 2> {log}"

rule severus_flag_bnd_known:
	input:
		"workflow/scripts/SV_known_flag.R",
		"results/{sample}/wgs/{sample}.severus.BND.merged.annot.tsv",
		config["BND_TARGET"],
		config["ATLAS_TARGET"],
		config["GENEXUS_FUSION"],
		config["BND_FALSE_POSITIVES"]
	output:
		"results/{sample}/wgs/{sample}.severus.BND.merged.annot.flagged.tsv"
	conda: 
		"../envs/fusion_analysis.yaml"
	log:
		"logs/{sample}/severus_flag_bnd_known.log"
	shell:
		"Rscript {input} {output} 2>{log}"

rule severus_flag_inv_known:
	input:
		"workflow/scripts/INV_DUP_known_flag.R",
		"results/{sample}/wgs/{sample}.severus.INV.merged.annot.tsv",
		"resources/genexus_INV_DUP.bed",
	output:
		"results/{sample}/wgs/{sample}.severus.INV.merged.annot.flagged.tsv"
	conda: 
		"../envs/fusion_analysis.yaml"
	log:
		"logs/{sample}/severus_flag_inv_known.log"
	shell:
		"Rscript {input} {output} 2>{log}"

rule severus_flag_dup_known:
	input:
		"workflow/scripts/INV_DUP_known_flag.R",
		"results/{sample}/wgs/{sample}.severus.DUP.merged.annot.tsv",
		"resources/genexus_INV_DUP.bed",
	output:
		"results/{sample}/wgs/{sample}.severus.DUP.merged.annot.flagged.tsv"
	conda: 
		"../envs/fusion_analysis.yaml"
	log:
		"logs/{sample}/severus_flag_dup_known.log"
	shell:
		"Rscript {input} {output} 2>{log}"

rule severus_create_final_report_sv:
	input:
		"workflow/scripts/create_final_report_sv.R",
		"results/{sample}/wgs/{sample}.severus.INV.merged.annot.flagged.tsv",
		"results/{sample}/wgs/{sample}.severus.BND.merged.annot.flagged.tsv",
		"results/{sample}/wgs/{sample}.severus.DUP.merged.annot.flagged.tsv",
	output:
		"results/svcallers_comparison/{sample}/{sample}.severus.SV_final_report.xlsx"
	conda:
		"../envs/fusion_analysis.yaml"
	log:
		"logs/{sample}/severus_create_final_report_sv.log"
	shell:
		"Rscript {input} {output} 2> {log}"


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