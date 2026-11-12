#!/bin/bash
set -euxo pipefail

# Function to print usage information
usage() {
	echo "Usage: bash $0 --input_vcf <input_vcf/> --output_file <merged_bam>"
	echo ""
	echo "Options:"
	echo "   --AnnotSV               Path to annotSV executable file"
	echo "   --input_vcf_inv         Input vcf for INV from sniffles2.7"
	echo "   --input_vcf_bnd         Input vcf for BND from sniffles2.7"
	echo "   --tmp_dir               Path to tmp files"
	echo "   --sample_name           Name of the sample"
	echo "   --output_inv_tsv        Output file tsv for INV"
   	echo "   --output_bnd_tsv        Output file tsv for BND"
    exit 1
}

#Check if no arguments were provided
if [ "$#" -eq 0 ]; then
    echo "Error: No arguments provided."
    usage
fi

# Parse command-line arguments
input_vcf_inv=""
input_vcf_bnd=""
AnnotSV=""
tmp_dir=""
sample=""

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --input_vcf_inv)
            input_vcf_inv=$2
            shift 2
            ;;
        --input_vcf_bnd)
            input_vcf_bnd="$2"
            shift 2
            ;;
        --AnnotSV)
        	AnnotSV="$2"
        	shift 2
        	;;
        --tmp_dir)
            tmp_dir="$2"
            shift 2
            ;;
        --sample_name)
            sample="$2"
            shift 2
            ;;
        --output_inv_tsv)
            output_inv_tsv="$2"
            shift 2
            ;;
        --output_bnd_tsv)
            output_bnd_tsv="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Error: Invalid argument '$1'."
            usage
            ;;
    esac
done

if [[ -f "$output_inv_tsv" ]]; then 

	rm "$output_inv_tsv"

fi

######################
######## INV #########
######################
svtools vcftobedpe --input "$input_vcf_inv" --output "$tmp_dir"/"$sample".sniffles2.INV.filtPopVars.bedpe 

bcftools query -f'%INFO/STDEV_LEN\t%INFO/STDEV_POS\t%INFO/COVERAGE\t[%DR\t%DV]\t%INFO/VAF' "$input_vcf_inv" > "$tmp_dir"/"$sample".sniffles2.INV.infoField.tsv

grep -v '^#' "$tmp_dir"/"$sample".sniffles2.INV.filtPopVars.bedpe | \
awk '{OFS="\t"; print $1,$2-1001,$2+1000,$2,$4,$5-1001,$5+1000,$5,$7,$11}' > "$tmp_dir"/"$sample".sniffles2.INV.filtPopVars.tmp.bedpe 

rm "$tmp_dir"/"$sample".sniffles2.INV.filtPopVars.bedpe

header=false
cat "$tmp_dir"/"$sample".sniffles2.INV.filtPopVars.tmp.bedpe | while read line; do 

	echo "$line"
	echo "$line" | cut -f1,2,3,4,9,10 > "$tmp_dir"/"$sample".sniffles2.INV1.filtPopVars.tmp.bed
	echo "$line" | cut -f5,6,7,8,9,10  > "$tmp_dir"/"$sample".sniffles2.INV2.filtPopVars.tmp.bed

	"$AnnotSV" \
		-annotationMode full \
		-svtBEDcol 6 \
		-SVinputFile "$tmp_dir"/"$sample".sniffles2.INV1.filtPopVars.tmp.bed \
		-genomeBuild GRCh38 \
		-outputFile "$tmp_dir"/"$sample".sniffles2.INV1.filtPopVars.tmp.annotSV.tsv \
		-outputDir "$tmp_dir" \
		-overlap 60 
	
	"$AnnotSV" \
		-annotationMode full \
		-svtBEDcol 6 \
		-SVinputFile "$tmp_dir"/"$sample".sniffles2.INV2.filtPopVars.tmp.bed \
		-genomeBuild GRCh38 \
		-outputFile "$tmp_dir"/"$sample".sniffles2.INV2.filtPopVars.tmp.annotSV.tsv \
		-outputDir "$tmp_dir" \
		-overlap 60 

	rm "$tmp_dir"/"$sample".sniffles2.INV1.filtPopVars.tmp.bed "$tmp_dir"/"$sample".sniffles2.INV2.filtPopVars.tmp.bed
	
	if ! "$header"; then 

		paste "$tmp_dir"/"$sample".sniffles2.INV1.filtPopVars.tmp.annotSV.tsv \
			"$tmp_dir"/"$sample".sniffles2.INV2.filtPopVars.tmp.annotSV.tsv >> "$tmp_dir"/"$sample".sniffles2.INV1-2.filtPopVars.tmp.annotSV.tsv

		header=true		
	else

		paste <(tail -n+2 "$tmp_dir"/"$sample".sniffles2.INV1.filtPopVars.tmp.annotSV.tsv) \
			<(tail -n+2 "$tmp_dir"/"$sample".sniffles2.INV2.filtPopVars.tmp.annotSV.tsv) >> "$tmp_dir"/"$sample".sniffles2.INV1-2.filtPopVars.tmp.annotSV.tsv
	fi

	rm "$tmp_dir"/"$sample".sniffles2.INV1.filtPopVars.tmp.annotSV.tsv "$tmp_dir"/"$sample".sniffles2.INV2.filtPopVars.tmp.annotSV.tsv

done

echo $(head -1 "$tmp_dir"/"$sample".sniffles2.INV1-2.filtPopVars.tmp.annotSV.tsv | \
	sed 's/user#1/Annotate_Length/g;s/user#2/Real_bkr/g') SVLEN STDEV_LEN STDEV_POS COVERAGE DR DV VAF | \
	tr ' ' '\t' > "$tmp_dir"/"$sample".header.INV.annotSV.tsv

cat "$tmp_dir"/"$sample".header.INV.annotSV.tsv \
	<(paste <(tail -n+2 "$tmp_dir"/"$sample".sniffles2.INV1-2.filtPopVars.tmp.annotSV.tsv) "$tmp_dir"/"$sample".sniffles2.INV.infoField.tsv) > "$output_inv_tsv"

rm "$tmp_dir"/"$sample".sniffles2.INV.filtPopVars.tmp.bedpe \
	"$tmp_dir"/"$sample".sniffles2.INV1-2.filtPopVars.tmp.annotSV.tsv \
	"$tmp_dir"/"$sample".header.INV.annotSV.tsv


##############################
############ BND #############
##############################
svtools vcftobedpe --input "$input_vcf_bnd" --output "$tmp_dir"/"$sample".sniffles2.BND.filtPopVars.bedpe

bcftools query -f'%INFO/STDEV_POS\t%INFO/COVERAGE\t[%DR\t%DV]\t%INFO/VAF' "$input_vcf_inv" > "$tmp_dir"/"$sample".sniffles2.BND.infoField.tsv

grep -v '^#' "$tmp_dir"/"$sample".sniffles2.BND.filtPopVars.bedpe  | \
awk '{OFS="\t"; print $1,$2-1001,$2+1000,$2,$4,$5-1001,$5+1000,$5,$7,$11}' > "$tmp_dir"/"$sample".sniffles2.BND.filtPopVars.tmp.bedpe 

rm "$tmp_dir"/"$sample".sniffles2.BND.filtPopVars.bedpe

header=false
cat "$tmp_dir"/"$sample".sniffles2.BND.filtPopVars.tmp.bedpe | while read line; do 

	echo "$line"
	echo "$line" | cut -f1,2,3,4,9,10 > "$tmp_dir"/"$sample".sniffles2.BND1.filtPopVars.tmp.bed
	echo "$line" | cut -f5,6,7,8,9,10  > "$tmp_dir"/"$sample".sniffles2.BND2.filtPopVars.tmp.bed

	"$AnnotSV" \
		-annotationMode full \
		-svtBEDcol 6 \
		-SVinputFile "$tmp_dir"/"$sample".sniffles2.BND1.filtPopVars.tmp.bed \
		-genomeBuild GRCh38 \
		-outputFile "$tmp_dir"/"$sample".sniffles2.BND1.filtPopVars.tmp.annotSV.tsv \
		-outputDir "$tmp_dir" \
		-overlap 60 
	
	"$AnnotSV" \
		-annotationMode full \
		-svtBEDcol 6 \
		-SVinputFile "$tmp_dir"/"$sample".sniffles2.BND2.filtPopVars.tmp.bed \
		-genomeBuild GRCh38 \
		-outputFile "$tmp_dir"/"$sample".sniffles2.BND2.filtPopVars.tmp.annotSV.tsv \
		-outputDir "$tmp_dir" \
		-overlap 60 

	rm "$tmp_dir"/"$sample".sniffles2.BND1.filtPopVars.tmp.bed "$tmp_dir"/"$sample".sniffles2.BND2.filtPopVars.tmp.bed
	
	if ! "$header"; then 

		paste "$tmp_dir"/"$sample".sniffles2.BND1.filtPopVars.tmp.annotSV.tsv \
			"$tmp_dir"/"$sample".sniffles2.BND2.filtPopVars.tmp.annotSV.tsv >> "$tmp_dir"/"$sample".sniffles2.BND1-2.filtPopVars.tmp.annotSV.tsv

		header=true		
	else

		paste <(tail -n+2 "$tmp_dir"/"$sample".sniffles2.BND1.filtPopVars.tmp.annotSV.tsv) \
			<(tail -n+2 "$tmp_dir"/"$sample".sniffles2.BND2.filtPopVars.tmp.annotSV.tsv) >> "$tmp_dir"/"$sample".sniffles2.BND1-2.filtPopVars.tmp.annotSV.tsv
	fi

	#rm "$tmp_dir"/"$sample".sniffles2.BND1.filtPopVars.tmp.annotSV.tsv "$tmp_dir"/"$sample".sniffles2.BND2.filtPopVars.tmp.annotSV.tsv

done

echo $(head -1 "$tmp_dir"/"$sample".sniffles2.BND1-2.filtPopVars.tmp.annotSV.tsv | \
	sed 's/user#1/Annotate_Length/g;s/user#2/Real_bkr/g') SVLEN STDEV_LEN STDEV_POS COVERAGE DR DV VAF | \
	tr ' ' '\t' > "$tmp_dir"/"$sample".header.BND.annotSV.tsv

cat "$tmp_dir"/"$sample".header.BND.annotSV.tsv \
	<(paste <(tail -n+2 "$tmp_dir"/"$sample".sniffles2.BND1-2.filtPopVars.tmp.annotSV.tsv) "$tmp_dir"/"$sample".sniffles2.BND.infoField.tsv) > "$output_bnd_tsv"

rm "$tmp_dir"/"$sample".sniffles2.BND.filtPopVars.tmp.bedpe \
	"$tmp_dir"/"$sample".sniffles2.BND1-2.filtPopVars.tmp.annotSV.tsv \
	"$tmp_dir"/"$sample".header.BND.annotSV.tsv
