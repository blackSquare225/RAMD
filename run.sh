#!/bin/bash

# Function to print usage information
usage() {
    echo "Usage: bash $0 [--targeted --outdir <ont_file_location> --experiment_name <experiment_name> --roi_barcode <name_barcode_roi> --fusion_barcode <name_barcode_fusion> --time_roi_analysis <time_to_wait> --threads <num_of_threads>]"
    echo "Usage: bash $0 [--wgs --outdir <ont_file_location> --experiment_name <experiment_name> --time_wgs_analysis <time_to_wait> --n_wgs_run <num_of_run> --threads <num_of_threads>]"
    echo ""
    echo "Analysis mode:"
    echo "  --targeted            Analyse targeted data (requires --experiment_name --roi_dir and --fusion_dir)"
    echo "  --wgs                 Analyse wgs data (requires --experiment_name)"
    echo ""
    echo "Options:"
    echo "  --outdir              Path to output file location for sequencing run. (default /var/lib/minknow/data/RAMD)"
    echo "  --experiment_name     Name of the experiment in MinKNOW (sampleID must be named as <experiment_name>_targeted and <experiment_name>_wgs for --targeted and --wgs mode, respectively)"
    echo "  --roi_barcode         Name of the barcode used for roi multiplexing (only for --targeted mode)"
    echo "  --fusion_barcode      Name of the barcode used for fusion multiplexing (only for --targeted mode)"
    echo "  --time_roi_analysis   Time (integer) in hours to wait for starting roi analysis (only for --targeted mode). (default 3)"
    echo "  --time_wgs_analysis   Time (integer) in hours to wait for starting wgs analysis (only for --wgs mode). (default 50)"
    echo "  --n_wgs_run           Number (integer) of sequencing run (only for --wgs mode). (default 1)"
    echo "  --threads             Number of threads to use. (default: 3)"
    echo ""
    echo "  -h,--help             Display this help message."
    echo ""
    echo "Examples:"
    echo "bash $0 --targeted --experiment_name sample01 --roi_barcode barcode01 --fusion_barcode barcode02 --threads 10"
    echo "bash $0 --wgs --experiment_name sample01 --threads 10"
    exit 1
}

# Ensure that a conda environment is activated by checking CONDA_PREFIX
if [ -z "$CONDA_PREFIX" ]; then
    echo "Error: No Conda environment is activated. Please activate the snakemake environment with <conda activate snakemake>." >&2
    exit 1
fi

# Extract the name of the activated environment
ACTIVE_ENV=$(basename "$CONDA_PREFIX")

# Check that the activated environment is "snakemake"
if [ "$ACTIVE_ENV" != "snakemake" ]; then
    echo "Error: The activated environment is '$ACTIVE_ENV'. Please activate the snakemake environment with <conda activate snakemake>." >&2
    exit 1
fi

# Check if no arguments were provided
if [ "$#" -eq 0 ]; then
    echo "Error: No arguments provided."
    usage
fi

# Parse command-line arguments
targeted=false
wgs=false
outdir="/var/lib/minknow/data/RAMD"
experiment_name=""
roi_barcode=""
fusion_barcode=""
time_roi_analysis=3
time_wgs_analysis=50
n_wgs_run=1
threads=3

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --targeted)
            targeted=true
            shift
            ;;
        --wgs)
            wgs=true
            shift
            ;;
        --outdir)
            outdir="$2"
            shift 2
            ;;     
        --experiment_name)
            project_name="$2"
            shift 2
            ;;
        --targeted_dir)
            targeted_dir="$2"
            shift 2
            ;;
        --roi_barcode)
            roi_barcode="$2"
            shift 2
            ;; 
        --fusion_barcode)
            fusion_barcode="$2"
            shift 2
            ;;      
        --time_roi_analysis)
            time_roi_analysis="$2"
            shift 2
            ;;      
        --time_wgs_analysis)
            time_wgs_analysis="$2"
            shift 2
            ;;              
        --n_wgs_run)
            n_wgs_run="$2"
            shift 2
            ;;
        --threads)
            threads="$2"
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

# Ensure mutually exclusive options are enforced
if [ "$targeted" = true ] && [ "$wgs" = true ]; then
    echo "Error: You cannot select both targeted and wgs analysis."
    usage
fi

# Ensure at least one option is selected
if [ "$targeted" != true ] && [ "$wgs" != true ]; then
    echo "Error: You must select at least one analysis mode: targeted or wgs."
    usage
fi

# Check if output ont file location dir exist 
if [ ! -d "$outdir" ]; then
    echo "Directory does not exist: $outdir"
    usage
fi

# Create the config directory if it doesn't exist
#if [ ! -d "config/" ]; then 
#    mkdir config
#fi

## Handling targeted Mode
if [ "$targeted" = true ]; then

    # Ensure required arguments for targeted analysis
    if [ -z "$project_name" ] || [ -z "$roi_barcode" ] || [ -z "$fusion_barcode" ]; then
        echo "Error: targeted mode requires a project name (--project_name) and both roi (--roi_barcode) and fusion (--fusion_barcode) barcode names."
        usage
    fi
    
    if ! [[ "$threads" =~ ^[0-9]+$ ]] || ! [[ "$time_roi_analysis" =~ ^[0-9]+$ ]]; then
        echo "Error: targeted mode requires an integer number for both threads and time to wait for starting analysis"
        usage    
    fi

    # Generate valid barcodes 
    valid_barcodes=()
    for i in $(seq -w 1 96); do
        valid_barcodes+=("barcode$i")
    done

    # Flags to track if barcodes are valid
    roi_barcode_valid=false
    fusion_barcode_valid=false

    # Check if barcodes are valid 
    for b in "${valid_barcodes[@]}"; do
        if [[ "$roi_barcode" == "$b" ]]; then
            echo "Using ROI barcode: $b"
            roi_barcode_valid=true
        fi
        if [[ "$fusion_barcode" == "$b" ]]; then
            echo "Using Fusion barcode: $b"
            fusion_barcode_valid=true
        fi
    done

    # Exit if any barcode is invalid
    if [[ "$roi_barcode_valid" != true ]]; then
        echo "Error: roi barcode '$roi_barcode' is not valid."
        exit 1
    fi

    if [[ "$fusion_barcode_valid" != true ]]; then
        echo "Error: fusion barcode '$fusion_barcode' is not valid."
        exit 1
    fi

    # Step 1: Wait for the project directory 
    echo "Waiting for project directory to be created inside $outdir..."
    while true; do

        project_dir=$(find "$outdir" -mindepth 1 -maxdepth 1 -type d -name "$project_name")
        
        if [[ -n "$project_dir" && -d "$project_dir"  ]]; then
            
            echo "Project directory successfully created: $project_dir"
        
            break
        fi
        
        sleep 1     
    done    

    # Step 2: Wait for the project targeted directory inside the project directory
    echo "Waiting for project_targeted directory to be created inside $project_dir..."

    timeout=$((2 * 60))  
    elapsed=0

    while true; do
        
        targeted_dir=$(find "$project_dir" -mindepth 1 -maxdepth 1 -type d -name "${project_name}_targeted")


        if [[ -n "$targeted_dir" && -d "$targeted_dir" ]]; then
            echo "project_targeted directory successfully created: $targeted_dir"
            break
        fi

        sleep 1
        ((elapsed++))

        if (( elapsed >= timeout )); then
            echo "Error: project_targeted directory '${project_name}_targeted' was not created inside $project_dir within 30 minutes, maybe you have provided a wrong name."
            usage
        fi
    done

    # Step 3: Wait for the flowcell directory inside the flowcell directory
    echo "Waiting for a flowcell-run directory to be created inside $targeted_dir..."
    while true; do

        flowcell_dir=$(find "$targeted_dir" -mindepth 1 -maxdepth 1 -type d)

        if [[ -n "$flowcell_dir" && -d "$flowcell_dir"  ]]; then
            
            echo "Flowcell directory successfully created: $flowcell_dir"
        
            break
        fi
        
        sleep 1     
    done
 
    # Step 4: Wait for the fastq_pass directory inside the flowcell directory
    fastq_pass_dir="$flowcell_dir/fastq_pass"

    echo "Waiting for "$fastq_pass_dir" to be created..."
    while [ ! -d "$fastq_pass_dir" ]; do
    
        sleep 1  
    
    done

    echo "$fastq_pass_dir directory detected!"

    # Step 5: Wait for both roi_barcode and fusion_barcode directories inside fastq_pass
    fusion_dir="$fastq_pass_dir/$fusion_barcode"
    roi_dir="$fastq_pass_dir/$roi_barcode"
    
    echo "Waiting for both $fusion_dir and $roi_dir to be created..."
    while [[ ! -d "$fusion_dir" || ! -d "$roi_dir" ]]; do
        sleep 1  
    done
    
    echo "Both $fusion_dir and $roi_dir directories detected! Proceeding with analysis..."
    
    cat <<'EOF'
    ██████   █████  ███    ███ ██████
    ██   ██ ██   ██ ████  ████ ██   ██
    ██████  ███████ ██ ████ ██ ██   ██
    ██   ██ ██   ██ ██  ██  ██ ██   ██
    ██   ██ ██   ██ ██      ██ ██████
    EOF
    echo "Running targeted analysis"          
    snakemake targeted --config project_name="$project_name" roi_dir="$roi_dir" fusion_dir="$fusion_dir" time_roi_analysis="$time_roi_analysis" threads="$threads" --cores "$threads" --use-conda --rerun-incomplete


elif [ "$wgs" = true ]; then
    
    # Ensure required arguments for wgs analysis
    if [ -z "$project_name" ]; then
        echo "Error: wgs mode requires a project name (--project_name)"
        usage
    fi
     
    # Project directory must exist, because it was created in previous targeted analysis, check it
    project_dir="$outdir/$project_name"
    if [ ! -d "$project_dir" ]; then 
        echo "Project dir does not exist or you have provived a wrong project name" 
        usage
    fi

    if ! [[ "$threads" =~ ^[0-9]+$ ]] || ! [[ "$time_wgs_analysis" =~ ^[0-9]+$ ]] || ! [[ "$n_wgs_run" =~ ^[0-9]+$ ]]; then
        echo "Error: wgs mode requires an integer number for threads, time to wait and number of run"
        usage    
    fi

    # Step 1: Check for the creation of project_wgs directory 
    echo "Waiting for project_wgs directory to be created inside $project_dir..."

    timeout=$((2 * 60))  
    elapsed=0

    while true; do
        
        wgs_dir=$(find "$project_dir" -mindepth 1 -maxdepth 1 -type d -name "${project_name}_wgs")

        if [[ -n "$wgs_dir" && -d "$wgs_dir" ]]; then
            echo "project_wgs directory successfully created: $wgs_dir"
            break
        fi

        sleep 1
        ((elapsed++))

        if (( elapsed >= timeout )); then
            echo "Error: project_wgs directory '${project_name}_wgs' was not created inside $project_dir within 2 minutes, maybe you have provided a wrong name."
            usage
        fi
    done

    # Step 2: Wait for the flowcell directory inside the flowcell directory
    echo "Waiting for a flowcell-run directory to be created inside $wgs_dir..."
    while true; do

        flowcell_dir=($(find "$wgs_dir" -mindepth 1 -maxdepth 1 -type d))

        if [[ -n "$flowcell_dir" ]]; then

            if (( ${#flowcell_dir[@]} == n_wgs_run )); then 
            
                echo "successfully created ${#flowcell_dir[@]} flowcells directory"
        
                break
            
            fi
        
        fi
        
        sleep 1     
    done

    echo "Running wgs analysis"
    snakemake wgs --config project_name="$project_name" wgs_dir="$wgs_dir" time_wgs_analysis="$time_wgs_analysis" threads="$threads" --cores "$threads" --use-conda --rerun-incomplete

fi
