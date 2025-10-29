#!/bin/bash
set -e

CONFIG=$1

INPUT=`yq .input.gwas_raw "${CONFIG}"`
OUTPUT=`yq .input.gwas "${CONFIG}"`

echo INPUT:$INPUT 
echo OUTPUT:$OUTPUT

# find the target column
if [[ $INPUT == *.gz ]]; then
    headstr=$(zcat "$INPUT" | head -n 1)
else
    headstr=$(head -n 1 "$INPUT")
fi

expected_header=("snp" "a1" "a2" "freq|af1|af|frq" "b|beta" "se" "p" "n")
header=$(head -n 1 "$INPUT" | tr '[:upper:]' '[:lower:]')
column_indices=()
IFS=$'\t' read -ra header_fields <<< "$headstr"
for col in "${expected_header[@]}"; do
    regex="^\s*($col)\s*$"
    index=-1
    for i in "${!header_fields[@]}"; do
        if [[ "${header_fields[i],,}" =~ $regex ]]; then
            index=$i
            break
        fi
    done

    if [[ $index == -1 ]]; then
        echo "ERROR: can't find column $col"
        exit 1
    fi

    column_indices+=("$index")
done

# Print the standard output.
awk_param=""
for column in "${column_indices[@]}"; do
    awk_param+="\$"$((${column} + 1))"\"\\t\""
done
awk_param=$(echo "$awk_param" | sed 's/"\\t"$//') 

if [[ $INPUT == *.gz ]]; then
    zcat "$INPUT" | awk "{ print ${awk_param} }" > "$OUTPUT"
else
    awk -v header="$gwas_header" "{ print ${awk_param} }" "$INPUT" > "$OUTPUT"
fi
sed -i "1s/.*/SNP	A1	A2	freq	b	se	P	N/" "$OUTPUT"

echo "input file passed verification"
