#!/bin/bash

CONFIG_template=""
trait_name="trait"
GWAS_DATA=""
WORK_DIR="/storage/yangjianLab/guoyazhou/GAMMA_git"
OUTPUT="/storage/yangjianLab/guoyazhou/GAMMA_git_data"
SCRIPT_DIR="/storage/yangjianLab/guoyazhou/GAMMA_github/gamma-script/scripts"
report_dir=""
qtl_list=""
maf=""
peqtl_smr=""
peqtl_heidi=""
mr_enable=""
MR_snp_min=""
xQTL=""

parse_options() {
    while [[ $# -gt 0 ]]; do
		local value="${1#*=}"
		if [[ "$value" == "None" ]] || [[ -z "$value" ]]; then
			shift
			continue
		fi
        case "$1" in

            --CONFIG_template=*)
                CONFIG_template="${1#*=}"
                shift
                ;;

            --trait_name=*)
                trait_name="${1#*=}"
                shift
                ;;

            --GWAS_DATA=*)
                GWAS_DATA="${1#*=}"
                shift
                ;;

            --WORK_DIR=*)
                WORK_DIR="${1#*=}"
                shift
                ;;

            --OUTPUT=*)
                OUTPUT="${1#*=}"
                shift
                ;;

            --SCRIPT_DIR=*)
                SCRIPT_DIR="${1#*=}"
                shift
                ;;

            --report_dir=*)
                report_dir="${1#*=}"
                shift
                ;;

            --qtl_list=*)
                qtl_list="${1#*=}"
                shift
                ;;

            --xQTL=*)
                xQTL="${1#*=}"
                shift
                ;;

            --maf=*)
                maf="${1#*=}"
                shift
                ;;

            --peqtl_smr=*)
                peqtl_smr="${1#*=}"
                shift
                ;;

            --peqtl_heidi=*)
                peqtl_heidi="${1#*=}"
                shift
                ;;

            --mr_enable=*)
                mr_enable="${1#*=}"
                shift
                ;;

            --MR_snp_min=*)
                MR_snp_min="${1#*=}"
                shift
                ;;

            -*)
                echo "Error: Unknown option $1" >&2
                exit 1
                ;;

            *)
                echo "Error: Unknown argument $1" >&2
                exit 1
                ;;
        esac
    done
}

parse_options "$@"

mkdir -p ${OUTPUT}/GWAS/COJO_format
COJO_file=${OUTPUT}/GWAS/COJO_format/${trait_name}.txt

# CONFIG_template=`yq .yaml.template "${CONFIG_template}"`
mkdir -p ${WORK_DIR}/yaml_file 
CONFIG=${WORK_DIR}/yaml_file/${trait_name}.yaml 
cp ${CONFIG_template} ${CONFIG}

yq -i ".input.trait = \"$trait_name\"" "$CONFIG"
yq -i ".input.gwas_raw = \"/opt/workflow/workdir/$GWAS_DATA\"" "$CONFIG"
yq -i ".input.gwas = \"$COJO_file\"" "$CONFIG"
yq -i ".input.output = \"$OUTPUT\"" "$CONFIG"

yq -i ".script.work_path = \"$WORK_DIR\"" "$CONFIG"
yq -i ".script.path = \"$SCRIPT_DIR\"" "$CONFIG"

yq -i ".yaml.config = \"$CONFIG\"" "$CONFIG"
yq -i ".input.report_dir = \"$report_dir\"" "$CONFIG"

if [ -n "${maf}" ]; then
	yq -i ".smr.maf = $maf" "$CONFIG"
fi
if [ -n "${peqtl_smr}" ]; then
	yq -i ".smr.peqtl_smr = $peqtl_smr" "$CONFIG"
fi
if [ -n "${peqtl_heidi}" ]; then
	yq -i ".smr.peqtl_heidi = $peqtl_heidi" "$CONFIG"
fi
if [ -n "${MR_snp_min}" ]; then
	yq -i ".mr.MR_snp_min = $MR_snp_min" "$CONFIG"
fi
if [ -n "${mr_enable}" ]; then
	yq -i ".mr.enable = $mr_enable" "$CONFIG"
fi



# make new QTL_list
QTL_list_old=${SCRIPT_DIR}/L2G/MAGIC_data/QTL_data/QTL_list.txt
QTL_list_new=${WORK_DIR}/yaml_file/QTL_list.txt
> $QTL_list_new
IFS=',' read -r -a elements <<< "$qtl_list"
for element in "${elements[@]}"; do
	grep -E "^${element}" "$QTL_list_old" >> "$QTL_list_new"
done
yq -i ".magic.QTL_list = \"$QTL_list_new\"" "$CONFIG"

# make new QTL_clumped_list
QTL_clumped_list_old=${SCRIPT_DIR}/L2G/MR/QTL_clumped/QTL_clumped_list.txt
QTL_clumped_list_new=${WORK_DIR}/yaml_file/QTL_clumped_list.txt
> $QTL_clumped_list_new
yq -i ".mr.QTL_clumped_list = \"$QTL_clumped_list_new\"" "$CONFIG"
if [ "${mr_enable}" == "True" ]; then
	IFS=',' read -r -a elements <<< "$qtl_list"
	for element in "${elements[@]}"; do
		grep -E "^${element}" "$QTL_clumped_list_old" >> "$QTL_clumped_list_new"
	done
fi

if [ -n "${xQTL}" ]; then
	TGZ_FILE=/opt/workflow/workdir/${xQTL}
	TGZ_DIR="${WORK_DIR}/SMR/xQTL"
	mkdir -p "$TGZ_DIR"
	
	tar -xzvf "$TGZ_FILE" -C "$TGZ_DIR"
    find "$TGZ_DIR" -type f -name "*.gz" | while read -r file; do
        gzip -d "$file"
    done
	BESD_FILE=$(find "$TGZ_DIR" -type f -name "*.besd" | head -n 1)
	
	if [ -z "$BESD_FILE" ]; then
		echo "No .besd file found in the extracted directory."
		exit 1
	fi
	
	xQTL_PREFIX="${BESD_FILE%.*}"  
	echo -e "xQTL\t${xQTL_PREFIX}\tFALSE" >> "$QTL_list_new"
fi
