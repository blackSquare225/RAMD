import pysam
import re
from intervaltree import IntervalTree
import csv
import argparse

parser = argparse.ArgumentParser(
                    description='Detect rearrangment in target gene'
                    )

parser.add_argument('-i', '--input_bam', required=True, help="Input BAM file")         
parser.add_argument('-r', '--region', required=True, help="Region specified as chr:start-end over which find rearrangment")
parser.add_argument('-o', '--output_tsv', required=True, help="output TSV file with breakpoint coordinates")      

args = parser.parse_args()

bam_file = args.input_bam
target_coord = args.region
output_tsv = args.output_tsv

# Parse target coordinates
chrom, coords = target_coord.split(":")
start, end = map(int, coords.split("-"))

# Build interval tree for target gene
target_tree = IntervalTree()
target_tree.addi(start, end)

# Define helper functions
def overlaps_target(chromosome, s, e):
    return chromosome == chrom and bool(target_tree.overlap(s, e))
    
def mapped_length(cigartuples):
    return sum(length for op, length in cigartuples if op in {0, 2, 3, 7, 8})

def mapped_length_from_cigar(cigarstring):
    cigartuples = re.findall(r"(\d+)([MIDNSHP=X])", cigarstring)
    return sum(int(length) for length, op in cigartuples if op in "M=XDN")


# Scan BAM file
candidate_reads = []

with pysam.AlignmentFile(bam_file, "rb") as bam:
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue

        read_name = read.query_name
        read_len_full = read.query_length 
        primary_mapped_len = mapped_length(read.cigartuples)
        prim_chrom = read.reference_name
        prim_start = read.reference_start
        prim_end = prim_start + mapped_length(read.cigartuples)

        target_in_primary = overlaps_target(prim_chrom, prim_start, prim_end)
        sa_records = []

        # Collect SA entries
        if read.has_tag("SA"):
            sa_tag = read.get_tag("SA")
            sa_entries = sa_tag.strip(";").split(";")
            for entry in sa_entries:
                sa_chrom, sa_pos, sa_strand, sa_cigar, sa_mapq, sa_nm = entry.split(",")
                sa_start = int(sa_pos) - 1
                sa_end = sa_start + mapped_length_from_cigar(sa_cigar)
                sa_records.append({
                    "chrom": sa_chrom,
                    "start": sa_start,
                    "end": sa_end,
                    "strand": sa_strand,
                    "cigar": sa_cigar,
                    "mapped_length": mapped_length_from_cigar(sa_cigar)
                })

        # -----------------------------
        # Case 1: Target in primary → save SA partners
        # -----------------------------
        if target_in_primary and sa_records:
            for sa in sa_records:
                candidate_reads.append({
                    "read": read_name,
                    "read_length_full": read_len_full,
                    "primary_mapped_length": primary_mapped_len,
                    "supplementary_mapped_length": sa["mapped_length"],
                    
                    "target_type": "primary",
                    "target_chrom": prim_chrom,
                    "target_start": prim_start,
                    "target_end": prim_end,
                    
                    "partner_type": "supplementary",
                    "partner_chrom": sa["chrom"],
                    "partner_start": sa["start"],
                    "partner_end": sa["end"],
                    "partner_strand": sa["strand"],
                    "partner_cigar": sa["cigar"]
                })

        # -----------------------------
        # Case 2: Target only in SA → primary is partner
        # -----------------------------
        elif not target_in_primary:
            for sa in sa_records:
                if overlaps_target(sa["chrom"], sa["start"], sa["end"]):
                    candidate_reads.append({
                        "read": read_name,
                        "read_length_full": read_len_full,
                        "primary_mapped_length": primary_mapped_len,
                        "supplementary_mapped_length": sa["mapped_length"],
                        
                        "target_type": "supplementary",
                        "target_chrom": sa["chrom"],
                        "target_start": sa["start"],
                        "target_end": sa["end"],
                        
                        "partner_type": "primary",
                        "partner_chrom": prim_chrom,
                        "partner_start": prim_start,
                        "partner_end": prim_end,
                        "partner_strand": sa["strand"],
                        "partner_cigar": sa["cigar"]
                    })


# Write TSV file
header_tsv1 = [
    "read",
    "read_length_full",
    "primary_mapped_length",
    "supplementary_mapped_length",
    "target_type", "target_chrom", "target_start", "target_end", "target_strand",
    "partner_type", "partner_chrom", "partner_start", "partner_end",
    "partner_strand", "partner_cigar"
]

with open(output_tsv, "w", newline="") as tsvfile:
    writer = csv.DictWriter(tsvfile, fieldnames=header, delimiter="\t")
    writer.writeheader()
    for row in candidate_reads:
        writer.writerow(row)

print(f"TSV written to: {output_tsv}")
