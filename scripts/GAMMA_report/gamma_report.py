import argparse
import os
import re
import json
import time
import requests
import openai
from openai import AzureOpenAI
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--disease", help = "Input disease name")
parser.add_argument("-y", "--gamma_path", help = "Gamma task output file path (Task ID path)")
parser.add_argument("-g1", "--gene1", help = "Report generation target gene1, if no gene picked, genes with highest GAMMA will be represented", default=None)
parser.add_argument("-g2", "--gene2", help = "Report generation target gene2, if no gene picked, genes with highest GAMMA will be represented", default=None)
parser.add_argument("-g3", "--gene3", help = "Report generation target gene3, if no gene picked, genes with highest GAMMA will be represented", default=None)
parser.add_argument("-o", "--output_path", help = "Report outpur path")
parser.add_argument("-k", "--key_path", help = "Azure key path, first row key, second row endpoint")

args = parser.parse_args()

disease = args.disease
gamma_path = args.gamma_path
genes = [args.gene1, args.gene2, args.gene3]
report_out = args.output_path
key_path = args.key_path

api_version = '2023-05-15'
api_type = 'azure'
api_model = 'gpt-4'

##### key and endpoint #####
with open(key_path) as file:
    lines = [line.rstrip() for line in file]
    file.close()
api_key = lines[0]
api_base = lines[1]
##### key and endpoint #####

# gets the API Key from environment variable AZURE_OPENAI_API_KEY
client = AzureOpenAI(
    # https://learn.microsoft.com/en-us/azure/ai-services/openai/reference#rest-api-versioning
    api_version=api_version,
    # https://learn.microsoft.com/en-us/azure/cognitive-services/openai/how-to/create-resource?pivots=web-portal#create-a-resource
    azure_endpoint=api_base,
    api_key=api_key
)

def gene_summary(gene, disease):
    message_gene = messages + [{
        "role": "user", 
        "content": f"Please summarize the current studied relationship of gene {gene} and disease {disease} in one sentence start with {gene}."
    }]
    response_gene = client.chat.completions.create(
        model=api_model,
        messages=message_gene,
        temperature=0,
    )
    return response_gene.choices[0].message.content

print('Start report generation.')
### Backgrounds (GPT)
messages=[
    {"role": "system", "content": "You are a competent professor in human population genetics area, good at literature reading and summarizing, expert in analyzing data in genetics field."},
]
messages.append({"role":"user", "content": f"Write the background of {disease} studies in a scientific writing way."})
response = client.chat.completions.create(
    model=api_model,
    messages=messages,
    temperature=0,
)
messages.append({"role":response.choices[0].message.role, 
                 "content":response.choices[0].message.content})
background = response.choices[0].message.content
print('GPT done.')

### GWAS
gwas_intro = "GWAS have revealed important biological insights into complex diseases, offering the promise of the identification of novel drug targets and opportunities for treatments. Retrospective analyses have estimated that drugs with targets supported by human genetic evidence are twice as likely to be approved for clinical use."
fig_1_mht = f"Fig. 1: Manhattan plot displaying {disease} GWAS P values. The x-axis is the chromosome position, and the y-axis is the -log10(P-values), representing the association between each variant and {disease}. The red dotted line represents the genome-wide significance level (P=5e-8)."

gwas = pd.read_csv(os.path.join(gamma_path, 'GWAS/COJO_format/trait.txt'), sep='\t')
gwas_n = gwas['N'].max()

### clumping
clumping = pd.read_csv(os.path.join(gamma_path, 'Clumping/summary/trait.locus'), sep='\t')
clumping = clumping.sort_values(['chr', 'start']).reset_index(drop=True)
n_loci = len(clumping)
clumping_intro = f"To discern quasi-independent genetic variants within the GWAS loci, we performed Clumping analysis. This methodology leverages both the summary statistics from our {gwas_n} people GWAS analysis and a reference panel consisting of genotype data from 10,000 individuals, randomly selected from the UK Biobank, to assess the linkage disequilibrium (LD) patterns among SNPs. The Clumping procedure, typically executed through software tools such as PLINK, identifies and retains a 'lead' SNP with the lowest p-value within a locus while clumping together other SNPs in LD with the lead SNP. SNPs are clumped if they are within a specified physical distance and have an LD (r²) above a certain threshold with the lead SNP, ensuring that the retained SNPs are representative of independent association signals. Using this approach, we isolated {n_loci} distinct loci (Table 1). With these independently associated loci identified, we will advance to more targeted gene prioritization analyses, concentrating on the genetic variants that have emerged from the Clumping process."

### Gene prioritization
###################################################
prio_methods = "Copy the NAFLD methods here!"     #
###################################################

gamma = pd.read_csv(os.path.join(gamma_path, 'GAMMA/score/trait_GAMMA.summary'), sep='\t')
gamma = gamma.sort_values('GAMMA', ascending=False).reset_index(drop=True)
gamma_locus = gamma.loc[gamma.Lead_SNP.isin(clumping.Lead_SNP)]
gamma_locus = gamma_locus.loc[gamma_locus['GAMMA']>=1]

n_genes_per_locus = gamma_locus.groupby('Lead_SNP').size()
max_gamma_per_locus = gamma_locus.groupby('Lead_SNP')['GAMMA'].max()
min_gamma_per_locus = gamma_locus.groupby('Lead_SNP')['GAMMA'].min()

gene_max_per_locus = gamma_locus.loc[gamma_locus.groupby('Lead_SNP')['GAMMA'].idxmax()].sort_values('GAMMA', ascending=False).reset_index(drop=True)

n_locus_with_gene = sum([i >= 1 for i in n_genes_per_locus])
n_genes = sum(n_genes_per_locus)
max_gamma = max(max_gamma_per_locus)
min_max_gamma = min(max_gamma_per_locus)
min_gamma = min(min_gamma_per_locus)
n_max_3 = sum([i>=3 for i in max_gamma_per_locus])
n_max_6 = sum([i>=6 for i in max_gamma_per_locus])

if n_locus_with_gene > 10:
    n_top_gamma = 10
    n_locus_over = True
else:
    n_top_gamma = n_locus_with_gene
    n_locus_over = False

top_gamma_gene_list = list(gene_max_per_locus.gene_name)
complementary_gene_list = gamma_locus.loc[gamma_locus['GAMMA'] >= 3].gene_name
complementary_gene_list = [i for i in complementary_gene_list if i not in top_gamma_gene_list]
if len(complementary_gene_list) > 5:
    n_complementary = 5
else:
    n_complementary = len(complementary_gene_list)
complementary_gene_list = complementary_gene_list[:n_complementary]

gamma_summary_1 = f"In total, we prioritized {n_genes} genes across {n_loci} loci, with gamma scores ranging from {min_gamma} to {max_gamma}. The evidence supporting each gene is summarized in Fig.2. "
gamma_summary_2 = f"Of the {n_loci} GWAS loci, {n_max_3} prioritized genes have GAMMA score larger than 3. "

if n_top_gamma >= 2:
    if n_locus_over: 
        gamma_summary_3 = "The genes with the highest score in each locus – for example, " + ', '.join(top_gamma_gene_list[:n_top_gamma]) + ', and so on, ' 
    else:
        gamma_summary_3 = "The genes with the highest score in each locus – namely, " + ', '.join(top_gamma_gene_list[:n_top_gamma-1]) + ', and ' + top_gamma_gene_list[n_top_gamma-1] + ', '
    
    gamma_summary_3 = gamma_summary_3 + f"are more likely to be the causal genes for {disease}. However, this does not eliminate the possibility of other prioritized genes having the potential of being therapeutic targets for {disease}. We suggest that genes such as "
    
    if n_complementary >= 2:
        gamma_summary_3 = gamma_summary_3 + ', '.join(complementary_gene_list[:n_complementary-1]) + ', and ' + top_gamma_gene_list[n_top_gamma]
    elif n_complementary == 1:
        gamma_summary_3 = gamma_summary_3 + complementary_gene_list[0]
    
    gamma_summary_3 = gamma_summary_3 + ", which exhibit higher gamma scores (>=3), can also be considered as potential therapeutic targets."

gamma_summary = gamma_summary_1 + gamma_summary_2 + gamma_summary_3
fig_2_gamma = f"Fig. 2: Gene prioritization for {disease}. Summary of GAMMA scores for the prioritized genes underlying the {n_top_gamma} GWAS loci with highest GAMMA scores. The genes were prioritized at three different confidence levels, including GAMMA_V2G, GAMMA_xQTL, GAMMA_DEPICT + GAMMA_POPS, corresponding to variant-, locus-, and network-level evidence, respectively. The leftmost squares represent the position of the GWAS loci. The loci and potential risk genes are ranked based on their GAMMA scores, summarized by GAMMA_V2G, GAMMA_xQTL, GAMMA_DEPICT, and GAMMA_POPS, and colored in red. The GAMMA_V2G category assesses whether the conditionally independent variant associated with the prioritized genes has supporting evidence from functional annotations. These annotations include proximity to the transcription start site (TSS) or exon of the target gene, involvement in enhancer-promoter regulation, or participation in PCHi-C. The GAMMA_xQTL category shows the evidence from integrative analysis of GWAS and xQTL through SMR, COLOC, FUSION, GSMR, and OPERA analysis. Different types of evidence are color-coded according to their respective categories, with deeper colors indicating higher scores or more significant signals as defined in the legend."

### For the picked genes
gene_summary_all = []
genes = [i for i in genes if i in gamma.gene_name.values]

if len(genes)>0:
    pick_gene_name = genes
else:
    pick_gene_name = gene_max_per_locus['gene_name'][:3]
pick_gene_info = [gamma.loc[gamma['gene_name']==gene].iloc[0] for gene in pick_gene_name]
pick_gene_gamma = [gamma['GAMMA'][gamma['gene_name']==gene].iloc[0] for gene in pick_gene_name]

if len(genes)==0:
    pick_gene_intro = f"Here we highlighted three genes, {pick_gene_name[0]}, {pick_gene_name[1]}, and {pick_gene_name[2]} with highest corresponding GAMMA scores among all loci, namely {pick_gene_gamma[0]}, {pick_gene_gamma[1]}, and {pick_gene_gamma[2]}. "
else:
    pick_gene_intro = f"Here we highlighted {' and '.join(pick_gene_name)}, with corresponding GAMMA scores {' and '.join([str(i) for i in pick_gene_gamma])}. "

exon_annotation = pd.read_csv(os.path.join(gamma_path, 'V2G/summary/trait_Exon.annot.summary'), sep='\t', header=None)

print('Start gene report generation.')
# Generate paragraph for each gene prioritized
for pick_loc_id in range(len(pick_gene_name)):
    gene = pick_gene_name[pick_loc_id]
    gene_info = pick_gene_info[pick_loc_id]

    print(gene)
    # GPT one sentence summary
    gene_summary_gpt = gene_summary(gene, disease)

    # Number of SNPs test (Only one for clumping, kept anyway)
    gene_locus_clumping = clumping.loc[clumping['Lead_SNP']==gene_info['Lead_SNP']].sort_values('P')
    gene_locus_min_p_value = gene_locus_clumping['P'].min()

##############################################
    # Karma and Clumping mapping problem!    #
##############################################
    if len(gene_locus_clumping)==0:
        gene_summary_snps = f" In our analysis, {gene} is not in any clumping locus."
        gene_final = gene_summary_gpt + gene_summary_snps
        gene_summary_all.append(gene_final)
        continue
    
    gene_summary_snps = f" In our analysis, {len(gene_locus_clumping)} variant was prioritized in this locus through fine-mapping, and the lead variant was {gene_locus_clumping['Lead_SNP'].iloc[0]} with the most significant p-value of {gene_locus_min_p_value}. "
    # Min p value test
    if gene_locus_min_p_value > 5e-9:
        gene_summary_snps = gene_summary_snps + f"The lead SNP {gene_locus_clumping['Lead_SNP'].iloc[0]} did not exhibit a highly significant signal but merely passed the suggestive significance threshold. "

    # Get credible set
    cs_path = os.path.join(gamma_path, 'Wu_adj/Result')
    credible_sets = os.listdir(cs_path)

    snp_credible_set_files = [[j for j in credible_sets if (j.startswith(f'trait_{i}') and j.endswith('.CS'))][0] for i in gene_locus_clumping['Lead_SNP']]
    snp_credible_set = [pd.read_csv(os.path.join(cs_path, i), sep='\t') for i in snp_credible_set_files]

    # Exon test
    gene_exon_snp_list = exon_annotation.loc[exon_annotation[10]==gene][3]
    gene_exon_overlap = any(i in list(gene_locus_clumping['Lead_SNP']) for i in list(gene_exon_snp_list))
    if gene_info['Exon'] == 1:
        if gene_exon_overlap:
            gene_summary_exon = f"The lead variant, {' and '.join(gene_locus_clumping['Lead_SNP'])}, was located in the exon of the protein-coding gene {gene}. We also "
        else:
            gene_summary_exon = "We "
        
        credible_set_in_exon = [gene_locus_clumping['Lead_SNP'].iloc[i] for i,cs_snps in enumerate(snp_credible_set) if any(snp in list(cs_snps.SNP) for snp in list(gene_exon_snp_list))]
        credible_set_in_exon = [i for i in credible_set_in_exon if i not in gene_locus_clumping['Lead_SNP'].values]
        
        gene_summary_exon = gene_summary_exon + "prioritized a credible set of SNPs of each lead SNP through fine-mapping. "
        
        if len(credible_set_in_exon) > 0:
            gene_summary_exon = gene_summary_exon + f"Among these variants, some variants in the credible set of {' and '.join(credible_set_in_exon)} were located in the exons of {gene}. "
        else:
            gene_summary_exon = gene_summary_exon + f"No variants in the credible set other than {gene_locus_clumping['Lead_SNP'].iloc[0]} located in the exons of {gene}. "
    else:
        gene_summary_exon = f"The lead variant and the prioritized credible set of SNPs through fine-mapping are not located in the exon of {gene}. "

    # Where are the SNPs
    gen_locus_snp_credible_set = []
    for snp_idx in range(len(gene_locus_clumping['Lead_SNP'])):
        snp_idx_cset = snp_credible_set[snp_idx]
        snp_idx_cset['before_start'] = snp_idx_cset.POS < gene_info.start
        snp_idx_cset['after_end'] = snp_idx_cset.POS > gene_info.end
        if gene_info.strand == '-':
            snp_idx_cset.rename(columns={"before_start": "after_end", "after_end": "before_start"}, inplace=True)
        snp_idx_cset['within'] = (~snp_idx_cset['before_start']) & (~snp_idx_cset['after_end'])
        snp_idx_cset['exon'] = snp_idx_cset['SNP'].isin(gene_exon_snp_list)
        snp_idx_cset['lead_snp'] = gene_locus_clumping['Lead_SNP'].iloc[snp_idx]
        
        gen_locus_snp_credible_set.append(snp_idx_cset)

    gen_locus_snp_credible_set = pd.concat(gen_locus_snp_credible_set)

    # Closest TSS
    CTSS = pd.read_csv(os.path.join(gamma_path, 'V2G/summary/trait_closestTSS.annot.summary'), sep='\t', header=None)
    CTSS = CTSS.loc[CTSS[5]==gene]

    if gene_info['ClosestTSS'] > 0:
        CTSS = pd.read_csv(os.path.join(gamma_path, 'V2G/summary/trait_closestTSS.annot.summary'), sep='\t', header=None)
        CTSS = CTSS.loc[CTSS[5]==gene]
        gene_CTSS = gen_locus_snp_credible_set.loc[gen_locus_snp_credible_set['SNP'].isin(CTSS[3])]
        gene_summary_CTSS = f"The transciption start site (TSS) of {gene} is the closest TSS for {len(gene_CTSS)} variants in the credible set."
    else:
        gene_summary_CTSS = ""

    # ABC
    if gene_info['ABC'] > 0:
        ABC_data = pd.read_csv(os.path.join(gamma_path, 'V2G/summary/trait_ABC.annot.summary'), sep='\t', header=None)
        ABC_data = ABC_data.loc[ABC_data[8]==gene]
        gene_ABC = gen_locus_snp_credible_set.loc[gen_locus_snp_credible_set['SNP'].isin(ABC_data[3])]
        gene_summary_ABC = "The activity-by-contact (ABC) model suggested that genetic variants reside in "
        ABC_pos = []
        if sum(gene_ABC['before_start'] > 0):
            ABC_pos.append("5'")
        if sum(gene_ABC['after_end'] > 0):
            ABC_pos.append("3'")
        if sum(gene_ABC['exon'] > 0):
            ABC_pos.append("exon")
        else:
            if sum(gene_ABC['within'] > 0):
                ABC_pos.append("UTR or intron")
        gene_summary_ABC = gene_summary_ABC + '/'.join(ABC_pos) + f" regions that interact with the promoter of {gene}. "
    else:
        gene_summary_ABC = ""

    # EpiMap / RoadMap
    epi_road = []
    if gene_info['EpiMap'] > 0:
        epi_road.append('EpiMap')
    if gene_info['RoadMap'] > 0:
        epi_road.append('RoadMap')

    if len(epi_road) > 0:
        gene_summary_er = "Data from " + ' and '.join(epi_road) + f" demonstrated that these genetic variants in the regulatory regions exhibit co-accessibility with the promoter regions of {gene}. "
    else:
        gene_summary_er = ""

    # PCHi-C
    if gene_info['PCHiC'] > 0:
        gene_summary_pchic = f"PCHi-C data further confirmed that these genetic variants displayed interactions with the promoter of {gene}. "
    else:
        gene_summary_pchic = ""

    # Transistion
    if (gene_summary_ABC=='') & (gene_summary_er=='') & (gene_summary_pchic==''):
        gene_summary_ABC = "The variant to gene prioritization methods we used, namely, ABC, EpiMap, RoadMap, and PCHi-C did not find an interaction between the fine-mapped variants and the gene promoter. "
    if gene_info['GAMMA_xQTL'] > 0:
        gene_summary_transition = "In addition to the variant-level evidence, we also found locus-level evidence by integrating GWAS with xQTL data. "
    else:
        gene_summary_transition = f"In addition to the variant-level methods, we also tested locus-level evidence by integrating GWAS with xQTL data using SMR, COLOC, and OPERA, but found no evidence that {gene} is colocalized with {disease}. "

    # xQTL
    detailed_gamma = pd.read_csv(os.path.join(gamma_path, 'GAMMA/feature/trait_GAMMA.feature'), sep='\t')
    xQTL = detailed_gamma.loc[detailed_gamma.gene_name == gene]
    xqtl_columns = xQTL.columns

    # SMR
    smr_sets = ['e', 'p', 's']
    smr_most_sig = [[i] for i in smr_sets]

    for i in range(3):
        quan_trait = smr_sets[i]

        xqtl_psmr = xQTL[[i for i in xqtl_columns if i.startswith(f'p_SMR_{quan_trait}QTL_')]].melt()
        xqtl_pheidi = xQTL[[i for i in xqtl_columns if i.startswith(f'p_HEIDI_SMR_{quan_trait}QTL_')]].melt()

        smr_test_n = sum(~xqtl_psmr.value.isna())
        
        if gene_info[f'{quan_trait}SMR'] < 0.05 / max(smr_test_n, 1):
            min_psmr_entry = xqtl_psmr.loc[xqtl_pheidi.value > 0.05].sort_values('value')
            if len(min_psmr_entry) > 0:
                min_psmr_entry = min_psmr_entry.iloc[0]
                smr_tissue = re.sub(f'p_SMR_{quan_trait}QTL_', '', min_psmr_entry.variable)
                smr_p_value = min_psmr_entry.value
                smr_most_sig[i].append(smr_tissue)
                smr_most_sig[i].append(smr_p_value)

    gene_summary_smr = ""

    eps_dict = {'e':'gene expression', 'p':'protein expression', 's':'splicing isoform expression'}

    passed_smr = []
    for quan_trait_entry in smr_most_sig:
        if len(quan_trait_entry) == 3:
            passed_smr.append(f'{quan_trait_entry[0]}QTLs')

    if len(passed_smr) > 0:
        gene_summary_smr = gene_summary_smr + f"We discovered that these {disease} associated variants are {' and '.join(passed_smr)} for {gene}. "

    for quan_trait_entry in smr_most_sig:
        if len(quan_trait_entry) == 3:
            gene_summary_smr = gene_summary_smr + f"SMR analysis in {quan_trait_entry[1]} {quan_trait_entry[0]}QTL shows most significant {eps_dict[quan_trait_entry[0]]} asscociation results among all tissues tested. "

    # COLOC
    coloc_sets = ['e', 'p', 's']
    coloc_most_sig = [[i] for i in coloc_sets]

    for i in range(3):
        quan_trait = coloc_sets[i]

        xqtl_pp4 = xQTL[[i for i in xqtl_columns if i.startswith(f'pp4_COLOC_{quan_trait}QTL')]].melt()
        
        if gene_info[f'{quan_trait}COLOC'] > 0.8:
            max_pp4_entry = xqtl_pp4.sort_values('value', ascending=False).iloc[0]
            coloc_tissue = re.sub(f'(pp4_COLOC_{quan_trait}QTL_)(.*)', '\\2', max_pp4_entry.variable)
            coloc_pp4 = max_pp4_entry.value
            coloc_most_sig[i].append(coloc_tissue)
            coloc_most_sig[i].append(coloc_pp4)

    gene_summary_coloc = ""

    passed_coloc = []
    for quan_trait_entry in coloc_most_sig:
        if len(quan_trait_entry) == 3:
            passed_coloc.append(f'{quan_trait_entry[0]}QTLs')

    if len(passed_coloc) > 0:
        gene_summary_coloc = gene_summary_coloc + f"Colocalization results showed that {disease} associated variants are colocalzied with {' and '.join(passed_coloc)} for {gene}. "

    for quan_trait_entry in coloc_most_sig:
        if len(quan_trait_entry) == 3:
            gene_summary_coloc = gene_summary_coloc + f"We found {quan_trait_entry[1]} {quan_trait_entry[0]}QTL shows highest post-probability of colocalization among all tissues tested. "

    # GSMR / FUSION
    xqtl_gsmr = xQTL['p_GSMR'].iloc[0]
    xqtl_fusion = xQTL[[i for i in xqtl_columns if i.startswith(f'p_FUSION')]].melt()
    fusion_test_n = sum(~xqtl_fusion.value.isna())
    gene_summary_gsmr_fusion = ""

    gsmr_fusion = []
    if xqtl_gsmr < 0.05:
        gsmr_fusion.append('GSMR')

    if gene_info[f'FUSION'] < 0.05 / max(fusion_test_n, 1):
        gsmr_fusion.append('FUSION')
        
    if len(gsmr_fusion) > 0:
        gene_summary_gsmr_fusion = gene_summary_gsmr_fusion + f"{' and '.join(gsmr_fusion)} analysis also showed a significant result. "

    # Conclusion
    gene_summary_conclusion = ""
    if gene_info['GAMMA'] >= 3:
        gene_summary_conclusion = gene_summary_conclusion + f"All the evidence above suggests that {gene} may be the target gene underlying this locus, supported {gene} as one promising therapeutic target for {disease}. "

    # Gene Add up
    gene_final = gene_summary_gpt + gene_summary_snps + gene_summary_exon + gene_summary_CTSS + \
            gene_summary_ABC + gene_summary_er + gene_summary_pchic + \
            gene_summary_transition + \
            gene_summary_smr + gene_summary_coloc + gene_summary_gsmr_fusion + \
            gene_summary_conclusion

    gene_summary_all.append(gene_final)

### Out put
out_list = [background, gwas_intro, fig_1_mht, clumping_intro, prio_methods, gamma_summary, fig_2_gamma, pick_gene_intro] + gene_summary_all

with open(report_out, 'w') as f:
    for item in out_list:
        f.write("%s\n\n" % item)