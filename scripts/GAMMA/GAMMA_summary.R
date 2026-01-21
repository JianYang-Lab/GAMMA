args <- commandArgs(trailingOnly = TRUE)
trait_name=args[1]
OUTPUT=args[2]
Pharmaprojects_data_file=args[3]
gamma_gene_file=args[4]
R_functions_file=args[5]
MeSH_id=args[6]

print(trait_name)
print(OUTPUT)
print(Pharmaprojects_data_file)
print(gamma_gene_file)
print(R_functions_file)
print(MeSH_id)

source(R_functions_file)
suppressMessages({
library(data.table);
library(stringr);
library(dplyr)
})

### ********************************************************************************************
# genecode summary
## ********************************************************************************************
# 0. loading gene annotation -----------------------------------------------------

result=fread(gamma_gene_file)
result$GWAS_LOCUS=NA;result$Lead_SNP=NA;result$Lead_SNP_BP=NA;
# 1. COJO locus -----------------------------------------------------
# print("cojo test-----")
COJO_locus=fread(paste0(OUTPUT, "/Clumping/summary/",trait_name,".locus"), header=T)
for(j in 1:nrow(COJO_locus)){
  chr=as.numeric(gsub("chr","",COJO_locus$chr[j]))
  start=COJO_locus$start[j]
  end=COJO_locus$end[j]
  locus=COJO_locus$GWAS_LOCUS[j]
  lead_snp=COJO_locus$Lead_SNP[j]
  lead_snp_bp=COJO_locus$Lead_SNP_BP[j]

  index=which(result$chr==chr & result$start<=end & result$start>=start & result$end<=end & result$end>=start)
  result$GWAS_LOCUS[index]=locus
  result$Lead_SNP[index]=lead_snp
  result$Lead_SNP_BP[index]=lead_snp_bp
}

# print("test2----------")
## ********************************************************************************************
# Pharmaprojects_2022 Nelson Paper
## ********************************************************************************************

Pharmaprojects_data=fread(Pharmaprojects_data_file)
index=which(Pharmaprojects_data$MeSH_id == MeSH_id)
clinical_data = Pharmaprojects_data[index,]

result$`Highest Status Reached Value`=NA
index=match(result$Gene_ID, clinical_data$`Target Entrez Gene ID`, nomatch=0)
result$`Highest Status Reached Value`[which(index!=0)] = clinical_data$`Highest Status Reached Value`[index]


result <- result %>%
  mutate(
    preclinical = if_else(is.na(`Highest Status Reached Value`), 0, as.integer(`Highest Status Reached Value` >= 0)),
    phase1 = if_else(is.na(`Highest Status Reached Value`), 0, as.integer(`Highest Status Reached Value` >= 1)),
    phase2 = if_else(is.na(`Highest Status Reached Value`), 0, as.integer(`Highest Status Reached Value` >= 2)),
    phase3 = if_else(is.na(`Highest Status Reached Value`), 0, as.integer(`Highest Status Reached Value` >= 3)),
    approved = if_else(is.na(`Highest Status Reached Value`), 0, as.integer(`Highest Status Reached Value` == 4))
  )

print("----------")
## ********************************************************************************************
# Load Gene-based test such as MAGMA mBAT ... method name represent p values
## ********************************************************************************************
MAGMA_file=paste0(OUTPUT, "/MAGMA/summary/",trait_name,"_chrALL.genes.out")
result$z_MAGMA=result$MAGMA=NA
result$GAMMA_MAGMA=0

if(file.exists(MAGMA_file)){
data=fread(MAGMA_file,head=T,stringsAsFactors=F,data.table=F)

index=match(result$gene_id,  data$GENE, nomatch=0)
result$z_MAGMA[which(index!=0)]=data$ZSTAT[index]
result$MAGMA[which(index!=0)]=data$P[index]

result$GAMMA_MAGMA=ifelse(result$MAGMA < 0.05/length(which(!is.na(result$MAGMA))), 1, 0)
}



print("----------")
mBATcombo_file=paste0(OUTPUT, "/mBATcombo/summary/",trait_name,".gene.assoc.mbat")

result$mBATcombo=NA
result$GAMMA_mBATcombo=0
if(file.exists(mBATcombo_file)){
data=fread(mBATcombo_file,head=T,stringsAsFactors=F,data.table=F)

index=match(result$gene_id,  data$Gene, nomatch=0)
result$mBATcombo[which(index!=0)]=data$P_mBATcombo[index]

result$GAMMA_mBATcombo=ifelse(result$mBATcombo < 0.05/length(which(!is.na(result$mBATcombo))), 1, 0)
}

result$GAMMA_GWAS=0
result$GAMMA_GWAS=apply(result[,c("GAMMA_MAGMA","GAMMA_mBATcombo")],1,function(x) sum(x,na.rm=T))

print("----------")
## ********************************************************************************************
# GAMMA summary
## ********************************************************************************************

GAMMA_V2G=fread(paste0(OUTPUT,"/V2G/score/",trait_name,"_GAMMA_V2G.summary"))
GAMMA_V2G$GAMMA_ClosestTSS=GAMMA_V2G$ClosestTSS
GAMMA_V2G$GAMMA_Exon=GAMMA_V2G$Exon
GAMMA_V2G$GAMMA_ABC=GAMMA_V2G$ABC
GAMMA_V2G$GAMMA_EpiMap=GAMMA_V2G$EpiMap
GAMMA_V2G$GAMMA_RoadMap=GAMMA_V2G$RoadMap
GAMMA_V2G$GAMMA_PCHiC=GAMMA_V2G$PCHiC


GAMMA_xQTL=fread(paste0(OUTPUT,"/L2G/score/",trait_name,"_GAMMA_xQTL.summary"))
GAMMA_Network=fread(paste0(OUTPUT,"/Network/score/",trait_name,"_GAMMA_Network.summary"))

GAMMA_V2G_feature=fread(paste0(OUTPUT,"/V2G/feature/",trait_name,"_GAMMA_V2G.feature"))
GAMMA_L2G_feature=fread(paste0(OUTPUT,"/L2G/feature/",trait_name,"_GAMMA_xQTL.feature"))
GAMMA_Network_feature=fread(paste0(OUTPUT,"/Network/feature/",trait_name,"_GAMMA_Network.feature"))


# GAMMA <- merge(result, GAMMA_V2G[,-c("Gene_ID","gene_id","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP")], by = "gene_name", all = TRUE)
# GAMMA <- merge(GAMMA, GAMMA_xQTL[,-c("Gene_ID","gene_id","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP")], by = "gene_name", all = TRUE)
# GAMMA <- merge(GAMMA, GAMMA_Network[,-c("Gene_ID","gene_id","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP")], by = "gene_name", all = TRUE)
GAMMA <- cbind(result, GAMMA_V2G[,-c("Gene_ID","gene_id","gene_name","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP")])
GAMMA <- cbind(GAMMA, GAMMA_xQTL[,-c("Gene_ID","gene_id","gene_name","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP")])
GAMMA <- cbind(GAMMA, GAMMA_Network[,-c("Gene_ID","gene_id","gene_name","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP")])




GAMMA$GAMMA=apply(GAMMA[,c("GAMMA_GWAS","GAMMA_V2G","GAMMA_xQTL","GAMMA_Network")],1,sum)


GAMMA=GAMMA[,c("Gene_ID","gene_id","gene_name","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP","Highest Status Reached Value", "preclinical", "phase1", "phase2","phase3", "approved",
					# GAMMA summary
					"GAMMA","GAMMA_GWAS","GAMMA_V2G","GAMMA_xQTL","GAMMA_Network",
					# GAMMA GWAS information ---- gene-based test
					"MAGMA","mBATcombo",
					# GAMMA V2G information
					"DistanceTSS","ClosestTSS","Exon","ABC","EpiMap","RoadMap","PCHiC",
					# GAMMA L2G (xQTL) inforamtion
					"eSMR","sSMR","pSMR","eCOLOC","sCOLOC","pCOLOC","FUSION","GSMR", "MAGIC","eMAGIC", "sMAGIC", "pMAGIC", "mMAGIC", "hMAGIC", "caMAGIC",
					# GAMMA Network information
					"PoPS","DEPICT","PPR","RWR",
					#### GAMMA xQTL corresponding QTL name
					"eSMR_QTL_name","sSMR_QTL_name","pSMR_QTL_name", "eCOLOC_QTL_name", "sCOLOC_QTL_name", "pCOLOC_QTL_name","FUSION_QTL_name","eMAGIC_QTL_name", "sMAGIC_QTL_name", "pMAGIC_QTL_name", "mMAGIC_QTL_name", "hMAGIC_QTL_name", "caMAGIC_QTL_name",
					"eSMR_p_HEIDI","sSMR_p_HEIDI","pSMR_p_HEIDI","eSMR_probeID","sSMR_probeID","pSMR_probeID","eCOLOC_probeID","sCOLOC_probeID","pCOLOC_probeID","FUSION_probeID","eMAGIC_probeID", "sMAGIC_probeID", "pMAGIC_probeID", "mMAGIC_probeID", "hMAGIC_probeID", "caMAGIC_probeID",
					#### GAMMA GWAS information
					"GAMMA_MAGMA","GAMMA_mBATcombo",
					#### GAMMA V2G information	
					"GAMMA_ClosestTSS","GAMMA_Exon","GAMMA_ABC","GAMMA_EpiMap","GAMMA_RoadMap","GAMMA_PCHiC",
					#### GAMMA L2G (xQTL) inforamtion
					"GAMMA_eSMR","GAMMA_sSMR","GAMMA_pSMR","GAMMA_eCOLOC","GAMMA_sCOLOC","GAMMA_pCOLOC","GAMMA_FUSION","GAMMA_GSMR", "GAMMA_MAGIC","GAMMA_eMAGIC", "GAMMA_sMAGIC", "GAMMA_pMAGIC", "GAMMA_mMAGIC", "GAMMA_hMAGIC", "GAMMA_caMAGIC",
					#### GAMMA Network information
					"GAMMA_PoPS","GAMMA_DEPICT","GAMMA_RWR","GAMMA_PPR")]


# GAMMA summary
GAMMA_feature=cbind(GAMMA[,c("Gene_ID","gene_id","gene_name","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP","Highest Status Reached Value", "preclinical", "phase1", "phase2","phase3", "approved",
					"GAMMA","GAMMA_GWAS","GAMMA_V2G","GAMMA_xQTL","GAMMA_Network")],
					# GAMMA GWAS information 
					result[,c("MAGMA","z_MAGMA", "mBATcombo")],
					# GAMMA V2G information
					GAMMA[,c("DistanceTSS","ClosestTSS","Exon","ABC","EpiMap","RoadMap","PCHiC")],
					# GAMMA L2G (xQTL) inforamtion
					GAMMA_L2G_feature[,-c("Gene_ID","gene_id","gene_name","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP","GAMMA_xQTL")],
					# GAMMA Network feature
					GAMMA[,c("PoPS","DEPICT","PPR","RWR")])

# GAMMA_feature=cbind(result,
					# GAMMA_V2G_feature[,-c("Gene_ID","gene_id","gene_name","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP")],
					# GAMMA_L2G_feature[,-c("Gene_ID","gene_id","gene_name","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP")],
					# GAMMA_Network_feature[,-c("Gene_ID","gene_id","gene_name","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP")])

write.table(GAMMA,paste0(OUTPUT,"/GAMMA/score/",trait_name,"_GAMMA.summary"),row=F,col=T,quo=F,sep="\t")
write.table(GAMMA_feature,paste0(OUTPUT,"/GAMMA/feature/",trait_name,"_GAMMA.feature"),row=F,col=T,quo=F,sep="\t")


index=which(!is.na(GAMMA$GWAS_LOCUS))
GAMMA_plot=GAMMA[index,]
write.table(GAMMA_plot,paste0(OUTPUT,"/GAMMA/plot/",trait_name,"_GAMMA_plot.summary"),row=F,col=T,quo=F,sep="\t")
