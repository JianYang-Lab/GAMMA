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
copy_directory $OUTPUT/SMR/SMR_Portal $REPORT/plot
rm -rf $REPORT/plot/*.database