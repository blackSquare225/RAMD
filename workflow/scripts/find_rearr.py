import pysam 
import re 
import pybedtools
import argparse

parser = argparse.ArgumentParser(description="Fetch reads in a specific region and in a list of bed file coordinates at the same time")
parser.add_argument("-b", "--bam_file", required=True, type=str,help="Input BAM file")
parser.add_argument("-r", "--region", required=True, type=str, help="Define the region as chr:stat-end")
parser.add_argument("-l", "--bed_file", required=True, type=str, help="Define bed file with genomic coordinates")
parser.add_argument("-o", "--output_file", required=True, type=str, help="Define output bed file")
args = parser.parse_args()

bam_file=args.bam_file
rara_coordinates=args.region
#rara_coordinates="chr17:40309180-40357643"
bed_file=args.bed_file
output_bed=args.output_file

chrom_rara, positions_rara = rara_coordinates.split(':')
start_rara, end_rara = map(int, positions_rara.split('-'))

# write a pseudocode 
# 1. check for all primary alignments in reads that maps on RARA coordinates 
rara_reads = []
with pysam.AlignmentFile(bam_file, "rb") as bam:

	for read in bam.fetch(chrom_rara, start_rara, end_rara):
		
		if read.is_unmapped:
			continue

		if read.is_secondary:
			continue

		read_name = read.query_name
		prim_chrom = read.reference_name
		prim_start= int(read.reference_start) - 1 
		prim_mapped_length = sum(int(length) for op, length in read.cigartuples if op in {0,7,8,2,3})
		prim_end = prim_start + prim_mapped_length  
		
		if read.has_tag("SA"):
			sa_tag = read.get_tag("SA")
			supplementary_alignments = sa_tag.split(';')[:-1]

			for entry in supplementary_alignments:

				supp_chrom, supp_pos, supp_strand, supp_cigar, supp_mapq, supp_nm = entry.split(",")

				supp_cigar = re.findall(r'(\d+)([MIDNSHP=XB])', supp_cigar)
				supp_mapped_length = sum(int(length) for length, op in supp_cigar if op in "M=XDN")

				supp_start = int(supp_pos) -1
				supp_end = supp_start + supp_mapped_length

				rara_reads.append((supp_chrom, supp_start, supp_end, prim_chrom, prim_start, prim_end, read_name))


# import bed file 
rara_bed = pybedtools.BedTool(bed_file)
rara_reads = pybedtools.BedTool(rara_reads)

reads_int_bed = rara_reads.intersect(rara_bed)
reads_int_bed.saveas(output_bed, trackline='#All reads that map in RARA and in at least one of putative partner genes')



