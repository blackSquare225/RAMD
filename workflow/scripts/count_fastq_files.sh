#!/bin/bash

# Function to print usage information
usage() {
    echo "Usage: bash $0 [--input_dir <fusion_dir> --output_fq <fastq.gz> --n_files <n_of_files_to_consider>]"
    echo ""
    echo "Options:"
    echo "  --input_dir   Directory containing the fusion fastq files (each file has 4000 reads by default)"
    echo "  --output_fq   Merged FASTQ output from the --n_files input files"
    echo "  --n_files     Threshold of number of files to merge"
    echo "  -h,--help     Display this help message."
    echo "Examples:"
    echo "bash $0 "
    echo "bash $0 "
    exit 1
}

# Check if no arguments were provided
if [ "$#" -eq 0 ]; then
    echo "Error: No arguments provided."
    usage
fi

# Parse command-line arguments
fusion_dir=""
output_fq=""
n_files=""

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --input_dir)
            fusion_dir=$2
            shift 2
            ;;
        --output_fq)
            output_fq="$2"
            shift 2
            ;;
        --n_files)
            n_files="$2"
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

# Wait until the directory exists
while [ ! -d "$fusion_dir" ]; do
    echo "Waiting for directory to be created..."
    sleep 1
done

# Optional: wait until the directory is not empty
while [ -z "$(ls -A "$fusion_dir")" ]; do
    echo "Directory $fusion_dir exists but is empty. Waiting for files..."
    sleep 1
done

echo "Directory $fusion_dir exist and it is not empty"

#if [ ! -d "$(dirname "$output_fq1")" ]; then 
#    echo "Creating output dir for 500kbp fusion"
#    mkdir -p "$(dirname "$output_fq1")"
#fi

#if [ ! -d "$(dirname "$output_fq2")" ]; then 
#    echo "Creating output dir for 1Mbp fusion"
#    mkdir -p "$(dirname "$output_fq2")"
#fi

# Loop into directory and process new files
while true; do
    
    # Count FASTQ files
    count=$(find "$fusion_dir" -maxdepth 1 -name '*.fastq.gz' | wc -l)
    echo "Current file count: $count"

    # If --n_files threshold has been reached
    if (( count >= n_files )); then
        
            echo "Finally reached (($n_files*4000)) reads, merging them into $output_fq..."
            
            files_to_merge=$(ls "$fusion_dir"/*.fastq.gz | head -n "$n_files")
            cat $files_to_merge > "$output_fq"
            
            if [[ $? -eq 0 ]]; then
        
                echo "Successfully created: $output_fq"

                break
        
            else
        
                echo "Error merging files!"
                exit 1
                
            fi
            
    fi

    sleep 30  

done




