#!/bin/bash
set -e

# 1. GWAS data (user input)
CONFIG=$1
SCRIPT_DIR=`yq .script.path "${CONFIG}"`
OUTPUT=`yq .input.output ${CONFIG}`
REPORT=/opt/workflow/workdir/`yq .input.report_dir ${CONFIG}`


copy_directory() {
    local src_dir=$1
    local dest_dir=$2

    if [ ! -d "$src_dir" ]; then
        echo "Error: Source directory does not exist: $src_dir"
        return 1
    fi

    if [ ! -d "$dest_dir" ]; then
        mkdir -p "$dest_dir"
        if [ $? -ne 0 ]; then
            echo "Error: Failed to create destination directory: $dest_dir"
            return 1
        fi
    fi

    cp -a "$src_dir/." "$dest_dir"
    if [ $? -eq 0 ]; then
        echo "Copy from $src_dir to $dest_dir successful."
    else
        echo "Error: Copy failed."
        return 1
    fi
}

# GWAS
copy_directory $OUTPUT/GWAS/manhattan_plot $REPORT/GWAS/manhattan_plot
copy_directory $OUTPUT/COJO/summary $REPORT/COJO/summary
copy_directory $OUTPUT/Clumping/summary $REPORT/Clumping/summary
copy_directory $OUTPUT/MAGMA/summary $REPORT/MAGMA/summary
copy_directory $OUTPUT/mBATcombo/summary $REPORT/mBATcombo/summary

# V2G
copy_directory $OUTPUT/V2G/plot $REPORT/V2G/plot
copy_directory $OUTPUT/Wu_adj/summary $REPORT/Wu_adj/summary

# L2G
copy_directory $OUTPUT/MAGIC/summary $REPORT/MAGIC/summary
copy_directory $OUTPUT/L2G/summary $REPORT/L2G/summary

# Network
copy_directory $OUTPUT/PoPS/PoPS_score $REPORT/PoPS/PoPS_score
copy_directory $OUTPUT/DEPICT/output $REPORT/DEPICT/output
copy_directory $OUTPUT/RWR_PPR/summary $REPORT/RWR_PPR/summary

# GAMMA
# copy_directory $OUTPUT/GAMMA/text $REPORT/GAMMA/text
copy_directory $OUTPUT/GAMMA/score $REPORT/GAMMA/score
copy_directory $OUTPUT/GAMMA/plot $REPORT/GAMMA/plot