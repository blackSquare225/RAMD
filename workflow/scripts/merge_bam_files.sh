#!/bin/bash 

# Function to print usage information
usage() {
	echo "Usage: bash $0 --input_dir <input_dir/> --output_file <merged_bam>"
	echo ""
	echo "Options:"
	echo "   --input_dir             Input dir in which find bam files to merge"
	echo "   --output_bam            Merged output bam file"
	echo "   --threads               Number of threads to use for samtools merge"
    echo "   --time_wgs_analysis     time to wait for starting bam merging"
    exit 1
}

# Check if no arguments were provided
if [ "$#" -eq 0 ]; then
    echo "Error: No arguments provided."
    usage
fi

# Parse command-line arguments
wgs_dir=""
output_bam=""
threads=1
time_wgs_analysis=0

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --input_dir)
            wgs_dir=$2
            shift 2
            ;;
        --output_bam)
            output_bam="$2"
            shift 2
            ;;
        --threads)
        	threads="$2"
        	shift 2
        	;;
        --time_wgs_analysis)
            time_wgs_analysis="$2"
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

# Waiting for the bam_pass directory and bam files
bam_pass_dir=($(find "$wgs_dir" -maxdepth 3 -name bam_pass -type d))
echo "$bam_pass_dir"

while true; do

    all_good=true

    for dir in "${bam_pass_dir[@]}"; do
        
        if [[ ! -d "$dir" ]]; then
        
            echo "Directory does not exist: $dir"
        
            all_good=false
        
        elif [[ -z "$(ls -A "$dir")" ]]; then
        
            echo "Directory is empty: $dir"
        
            all_good=false
        
        else
        
            echo "Directory ready: $dir"
        
        fi
    
    done

    if [[ "$all_good" == true ]]; then
        
        echo "All bam_pass directories exist and are not empty. Proceeding..."
        
        break
    
    fi

    echo "Not all directories are ready. Checking again in 5 seconds..."
    
    sleep 5

done


# Wait specified time before merging bam files 
echo "Analysis will start in $((time_wgs_analysis)) hours"

sleep $((3600 * time_wgs_analysis)) 

echo "ANALYSIS STARTED"

echo "Merging bam files"
ulimit -n 65535
samtools merge -@ "$threads" -o "$output_bam" $(find "$wgs_dir" -maxdepth 4 -name "*.bam" -type f)

