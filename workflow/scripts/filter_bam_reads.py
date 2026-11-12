#!/usr/bin/env python3

import pysam
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Filter a BAM file by read length using pysam."
    )
    parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file")
    parser.add_argument("-l", "--length", type=int, default=2000,
                        help="Minimum read length to keep (default: 2000 bp)")
    
    args = parser.parse_args()

    inbam = pysam.AlignmentFile(args.input, "rb")
    outbam = pysam.AlignmentFile(args.output, "wb", template=inbam)

    min_len = args.length

    for read in inbam:
        if read.query_length and read.query_length >= min_len:
            outbam.write(read)

    inbam.close()
    outbam.close()

    print(f"Filtering complete. Output written to: {args.output}")

if __name__ == "__main__":
    main()
