import pysam 
import re 
import pybedtools
import argparse
#import pandas as pd

def create_bed_windows(chrom, pos, right_window=100000, left_window=100000):
	"""
	Creates a BED window based on position and window size.

	Parameters:
		chrom (str): Chromosome name.
		pos (int): Position on the chromosome (1-based).
		window_size (int): Size of the window.

	Returns:
		tuple: Chromosome, start, and end positions of the window.
	"""
	
	start = max(0, int(pos) - int(left_window))
	end = int(pos) + int(right_window) - 1 # because bed are 0-based (chrom, pos-1, pos)
	return chrom, start, end


def find_and_annotate_chimeric_reads(bam_file, region, annotation_bed, output_tsv):
	"""
	Extract information and Annotates primary and supplementary alignments for each read.

	Parameters:
		bam_file (str): Path to the input BAM file.
		region (str): Genomic region in the format 'chr:start-end'.
		gene_bed_file (str): Path to the gene annotations BED file.

	Returns:
		list_of_tuples: Gene annotation for each  .
	"""
	chrom, positions = region.split(':')
	start, end = map(int, positions.split('-'))
	
	# Gene BED as BedTool object
	gene_bed = pybedtools.BedTool(annotation_bed)
	
	with open(output_tsv, "w") as f:

	# write the header 
		f.write("prim_chrom\tprim_start\tprim_end\tprim_pos\tprim_strand\tprim_mapq\tprim_mapped_length\tprim_genes\tsupp_chrom\tsupp_start\tsupp_end\tsupp_pos\tsupp_strand\tsupp_mapq\tsupp_mapped_length\tsupp_genes\tread_id\tread_quality\n")

		with pysam.AlignmentFile(bam_file, "rb") as bam:
			for read in bam.fetch(chrom, start, end):
				if read.has_tag("SA"):
					
					# Extract primary alignment
					read_id = read.query_name
					#read_quality = sum(read.query_qualities) / len(read.query_qualities)

					prim_chrom = chrom 
					prim_pos = int(read.reference_start)
					prim_strand = "+" if read.is_forward else "-"
					prim_mapq = read.mapping_quality
					prim_mapped_length = sum(length for op, length in read.cigartuples if op in {0, 2, 3})
					
					prim_chrom, prim_start, prim_end = create_bed_windows(chrom, prim_pos)

					# Create primary BED line and annotate
					primary_bed = pybedtools.BedTool(f"{prim_chrom}\t{prim_start}\t{prim_end}\t{prim_pos}\t{prim_strand}\t{prim_mapq}\t{prim_mapped_length}", from_string=True)
					primary_annotated = primary_bed.intersect(gene_bed, wa=True, wb=True)
					
					# Collapse gene annotations
					prim_genes = ";".join(set(line.split('\t')[-1].strip() for line in str(primary_annotated).splitlines()))
					
					# Extract supplementary alignments
					sa_tag = read.get_tag("SA")
					supplementary_alignments = sa_tag.split(';')[:-1]
					
					for supp_alignment in supplementary_alignments:
						supp_fields = supp_alignment.split(',')
						supp_chrom = supp_fields[0]
						supp_pos = int(supp_fields[1])
						supp_strand = supp_fields[2]  # Strand of the supplementary alignment
						supp_mapq = int(supp_fields[4])  # Mapping quality of the supplementary alignment
						supp_cigar = re.findall(r'(\d+)([A-Za-z=])', supp_fields[3]) # create a list of tuples from CIGAR string of the supplementary alignment
						supp_mapped_length = sum(int(length) for length, op in supp_cigar if op in "M=X") # Match or gap operations					
						
						supp_chrom, supp_start, supp_end = create_bed_windows(supp_chrom, supp_pos)

						# Create supplementary BED line and annotate
						supplementary_bed = pybedtools.BedTool(f"{supp_chrom}\t{supp_start}\t{supp_end}\t{supp_pos}\t{supp_strand}\t{supp_mapq}\t{supp_mapped_length}", from_string=True)
						supplementary_annotated = supplementary_bed.intersect(gene_bed, wa=True, wb=True)
						
						# Retain only gene annotations
						supp_genes = ";".join(set(line.split('\t')[-1].strip() for line in str(supplementary_annotated).splitlines()))
					
						annotatated_read = (prim_chrom, prim_start, prim_end, prim_pos, prim_strand, prim_mapq, prim_mapped_length, prim_genes, 
						supp_chrom, supp_start, supp_end, supp_pos, supp_strand, supp_mapq, supp_mapped_length, supp_genes, read_id)

						f.write("\t".join(map(str, annotatated_read)) + "\n")

def main():
	# Set up argument parser
	parser = argparse.ArgumentParser(description="Find supplementary alignments in a specified genomic region.")
	parser.add_argument("-b", "--bam_file", required=True, type= str, help="Path to the input BAM file.")
	parser.add_argument("-r", "--region", required=True, type= str, help="Genomic region in the format 'chr:start-end'.")
	parser.add_argument("-a", "--annotation_bed", required=True, type=str, help="bed file for genomic coordinates annotation")
	parser.add_argument("-o", "--output_tsv", required=True, type= str, help="Path to the output TSV file.")

	# Parse arguments
	args = parser.parse_args()

	# Call the function to find and annotate chimeric alignments 
	find_and_annotate_chimeric_reads(args.bam_file, args.region, args.annotation_bed, args.output_tsv)

	print(f"Annotated chimeric reads succesfully saved to {args.output_tsv}")

if __name__ == "__main__":
	main()