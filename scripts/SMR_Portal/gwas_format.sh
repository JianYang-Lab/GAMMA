#!/bin/bash
set -e

CONFIG=$1

INPUT=`yq .input.gwas_raw "${CONFIG}"`
OUTPUT=`yq .input.gwas "${CONFIG}"`

echo INPUT:$INPUT 
echo OUTPUT:$OUTPUT

/opt/workflow/bin/gwas_format -i "$INPUT" -o "$OUTPUT" -s 1
