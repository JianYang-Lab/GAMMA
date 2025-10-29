args <- commandArgs(trailingOnly = TRUE)
trait_name=args[1]
GWAS_DATA=args[2]
OUTPUT=args[3]
reference_all_bim_file=args[4]

suppressMessages({
library(data.table)
library(dplyr)
})

GWAS=fread(GWAS_DATA, header=T,col.names=c("SNP","A1","A2","AF1","BETA","SE","P","N"))
GWAS$P=as.numeric(GWAS$P)


GRCh38_reference=fread(reference_all_bim_file, col.names=c("CHR","SNP","MOL","POS","A1","A2"))
GWAS$CHR=NA
GWAS$POS=NA
index=match(GWAS$SNP, GRCh38_reference$SNP, nomatch=0)
GWAS$CHR[which(index!=0)]=GRCh38_reference$CHR[index]
GWAS$POS[which(index!=0)]=GRCh38_reference$POS[index]
GWAS=GWAS[which(index!=0),]

GWAS=GWAS[,c("SNP","CHR","POS","P")]
GWAS=GWAS[which(GWAS$P < 0.05),]

GWAS$P[which(GWAS$P <= 1e-300)]=2.225e-308
write.table(GWAS, paste0(OUTPUT,"/GWAS/manhattan_plot/",trait_name,"_manhattan_plot.txt"),sep = "\t", quote = F, na = "NA", row.names = F, col.names = T)

# GAMMA summary ---------------

GAMMA=fread(paste0(OUTPUT, "/GAMMA/score/",trait_name,"_GAMMA.summary"))
maxGAMMA <- GAMMA %>%
  filter(!is.na(GWAS_LOCUS)) %>%  # Filtering non-NA GWAS_LOCUS values right away
  group_by(GWAS_LOCUS) %>%
  filter(GAMMA == max(GAMMA, na.rm = TRUE)) %>%  # Filtering for max GAMMA within each group
  arrange(desc(`Highest Status Reached Value`)) %>%  # Ensuring correct column reference
  slice(1) %>%  # Selecting the first entry post-arrangement
  ungroup()

maxGAMMA$P=NA
index=match(maxGAMMA$Lead_SNP, GWAS$SNP, nomatch=0)
maxGAMMA$P[which(index!=0)]=GWAS$P[index]
maxGAMMA=maxGAMMA[,c("GWAS_LOCUS","gene_name","Lead_SNP","P")]

fwrite(maxGAMMA,paste0(OUTPUT,"/GWAS/manhattan_plot/",trait_name,"_manhattan_maxGAMMA.txt"),sep = "\t", quote = F, na = "NA", row.names = F, col.names = T)

closestGene = GAMMA %>%
  filter(!is.na(GWAS_LOCUS)) %>%  # Filtering non-NA GWAS_LOCUS values right away
  group_by(GWAS_LOCUS) %>%
  filter(DistanceTSS == min(DistanceTSS, na.rm=TRUE))  %>%
  arrange(desc(`Highest Status Reached Value`)) %>% 
  slice(1) %>% 
  ungroup()

closestGene$P=NA
index=match(closestGene$Lead_SNP, GWAS$SNP, nomatch=0)
closestGene$P[which(index!=0)]=GWAS$P[index]
closestGene=closestGene[,c("GWAS_LOCUS","gene_name","Lead_SNP","P")]

fwrite(closestGene,paste0(OUTPUT,"/GWAS/manhattan_plot/",trait_name,"_manhattan_closestGene.txt"),sep = "\t", quote = F, na = "NA", row.names = F, col.names = T)