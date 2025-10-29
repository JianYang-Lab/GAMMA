args <- commandArgs(trailingOnly = TRUE)
trait_name=args[1]
OUTPUT=args[2]
gamma_gene_file=args[3]
R_functions_file=args[4]

source(R_functions_file)
suppressMessages({
library(data.table)
library(stringr)
library(dplyr)
})

### ******************************************************************************************** 
# genecode summary
## ********************************************************************************************
# 0. loading gene annotation -----------------------------------------------------
result=fread(gamma_gene_file)
result$GWAS_LOCUS=NA;result$Lead_SNP=NA;result$Lead_SNP_BP=NA;

# 1. COJO locus -----------------------------------------------------

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

## ********************************************************************************************
# GAMMA Network summary
## ********************************************************************************************


# *****************************
# 1. loading PoPS_v2 results
# *****************************

PoPS_file = paste0(OUTPUT, "/PoPS/PoPS_score/", trait_name, "_PoPS.preds")
result$PoPS=NA;
result$GAMMA_PoPS=0;

if(file.exists(PoPS_file)){
data=fread(PoPS_file,head=T,stringsAsFactors=F,data.table=F)
index=match(result$gene_id,data$ENSGID,nomatch=0)
result$PoPS[which(index!=0)]=data$PoPS_Score[index]


maxcpp = data.frame(result %>%
  group_by(GWAS_LOCUS) %>%
  filter(if (all(is.na(PoPS))) {NA} else {PoPS == max(PoPS, na.rm = TRUE)}))

index=match(result$gene_id,maxcpp$gene_id,nomatch=0)
result$GAMMA_PoPS[which(index!=0)]=1

}



# *****************************
# 2. loading DEPICT results
# *****************************
DEPICT_file=paste0(OUTPUT, "/DEPICT/output/",trait_name,"_clean_geneprioritization.txt")
result$DEPICT=NA
result$GAMMA_DEPICT=0

if(file.exists(DEPICT_file)){
data=fread(DEPICT_file,head=T,stringsAsFactors=F,data.table=F)
index=match(result$gene_id,data[which(data[,"False discovery rate"]=="<=0.01"),"Ensembl gene ID"],nomatch=0)
result$GAMMA_DEPICT[which(index!=0)]=1

index=match(result$gene_id,data$`Ensembl gene ID`,nomatch=0)
result$DEPICT[which(index!=0)]=data$`Nominal P value`[index]
}


# *****************************
# 3. loading RWR and PPR results
# *****************************

result$RWR=result$PPR=NA
result$GAMMA_RWR=result$GAMMA_PPR=0

GAMMA_threshold="0"
RWR_PPR_file=paste0(OUTPUT,"/RWR_PPR/summary/",trait_name,"_",GAMMA_threshold,"_RWR_PPR.txt")


if(file.exists(RWR_PPR_file)){
	data=fread(RWR_PPR_file)
	data <- data %>%
		mutate(
			RWR_Score_Normalized = (RWR_Score - min(RWR_Score, na.rm = TRUE)) / (max(RWR_Score, na.rm = TRUE) - min(RWR_Score, na.rm = TRUE)),
      modified_RWR_Score_Normalized = (modified_RWR_Score - min(modified_RWR_Score, na.rm = TRUE)) / (max(modified_RWR_Score, na.rm = TRUE) - min(modified_RWR_Score, na.rm = TRUE)),
			PPR_Score_Normalized = (PPR_Score - min(PPR_Score, na.rm = TRUE)) / (max(PPR_Score, na.rm = TRUE) - min(PPR_Score, na.rm = TRUE))
			)
index=match(result$gene_id, data$gene_id, nomatch=0)
result$RWR[which(index!=0)]=data$modified_RWR_Score_Normalized[index]
result$PPR[which(index!=0)]=data$PPR_Score_Normalized[index]

maxcpp=data.frame(result %>% group_by(GWAS_LOCUS) %>% filter(RWR==max(RWR,na.rm = T)))
index=match(result$gene_id,maxcpp$gene_id,nomatch=0)
result$GAMMA_RWR[which(index!=0)]=1

maxcpp=data.frame(result %>% group_by(GWAS_LOCUS) %>% filter(PPR==max(PPR,na.rm = T)))
index=match(result$gene_id,maxcpp$gene_id,nomatch=0)
result$GAMMA_PPR[which(index!=0)]=1
}

## ********************************************************************************************
# Construct GAMMA results Network score
## ********************************************************************************************
result$GAMMA_Network=apply(result[,c("GAMMA_DEPICT","GAMMA_PoPS","GAMMA_RWR","GAMMA_PPR")],1,sum)


GAMMA_Network_name_list=c(## Gene information 
	"Gene_ID","gene_id","gene_name","chr","start","end","strand","GWAS_LOCUS","Lead_SNP","Lead_SNP_BP",
	## GAMMA Network information
	"GAMMA_Network", "GAMMA_PoPS","GAMMA_DEPICT","GAMMA_RWR","GAMMA_PPR",
	"PoPS","DEPICT","PPR","RWR")


result=result[,..GAMMA_Network_name_list]
feature=result

write.table(result,paste0(OUTPUT,"/Network/score/",trait_name,"_GAMMA_Network.summary"),row=F,col=T,quo=F,sep="\t")
write.table(feature,paste0(OUTPUT,"/Network/feature/",trait_name,"_GAMMA_Network.feature"),row=F,col=T,quo=F,sep="\t")


