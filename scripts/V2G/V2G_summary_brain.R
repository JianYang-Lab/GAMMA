args <- commandArgs(trailingOnly = TRUE)
trait_name=args[1]
OUTPUT=args[2]
gencode_file=args[3]
gamma_gene_file=args[4]
R_functions_file=args[5]
# locus_path=args[3]

source(R_functions_file)
suppressMessages({
library(data.table)
library(stringr)
library(dplyr)
library(GenomicRanges);
library(regioneR);
library(Repitools)
})

### ******************************************************************************************** 
# genecode summary
## ********************************************************************************************
# 0. loading gene annotation
# gencode_file="/storage/yangjianLab/qiting/data/annotation/gencode/gencode.v40.GRCh38.gene.annotation.bed"
gencode=fread(gencode_file,head=T,stringsAsFactors=F,data.table=F)
gencode$V1 <- sub("^chr", "", gencode$V1)
gencode$gene_id <- sub("\\..*", "", gencode$gene_id)

index=which(gencode$gene_type=="protein_coding" & gencode$V1!="M" & gencode$V1!="X" & gencode$V1!="Y")
gencode1=gencode[index,]
gencode1$gene_tss = ifelse(gencode1$V7 == "+", gencode1$V4, gencode1$V5)

result=fread(gamma_gene_file)
result$ClosestTSS=result$Exon=result$ABC=result$EpiMap=result$RoadMap=result$PCHiC=result$GAMMA_V2G=0;
result$GWAS_LOCUS=NA;result$Lead_SNP=NA;result$Lead_SNP_BP=NA;

# 1. COJO locus 
if(is.null(locus_path)){
  locus_path=paste0(OUTPUT, "/Clumping/summary/",trait_name,".locus")
}
COJO_locus=fread(locus_path, header=T)

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
# GAMMA V2G summary
## ********************************************************************************************
TSS=fread(paste0(OUTPUT,"/V2G/summary/",trait_name,"_closestTSS.annot.summary"),head=F,stringsAsFactors=F,data.table=F)
Exon=fread(paste0(OUTPUT,"/V2G/summary/",trait_name,"_Exon.annot.summary"),head=F,stringsAsFactors=F,data.table=F)
ABC=fread(paste0(OUTPUT,"/V2G/summary/",trait_name,"_ABC.annot.summary"),head=F,stringsAsFactors=F,data.table=F)
EpiMap=fread(paste0(OUTPUT,"/V2G/summary/",trait_name,"_EpiMap.annot.summary"),head=F,stringsAsFactors=F,data.table=F)
RoadMap=fread(paste0(OUTPUT,"/V2G/summary/",trait_name,"_RoadMap.annot.summary"),head=F,stringsAsFactors=F,data.table=F)
PCHiC=fread(paste0(OUTPUT,"/V2G/summary/",trait_name,"_PCHiC.annot.summary"),head=F,stringsAsFactors=F,data.table=F)



## ********************************************************************************************
# Brain tissue
## ********************************************************************************************
ABC_brain_related_cells <- c("Astrocyte 1",
                            "Astrocyte 2", 
                            "astrocyte-ENCODE", 
                            "BBB Endothelial", 
                            "Bipolar neuron from iPSC-ENCODE", 
                            "Ciliated", 
                            "Dendritic cell treated with Lipopolysaccharide", 
                            "Enteric Neuron", 
                            "GABA 1", 
                            "GABA 2", 
                            "Glutamatergic 1", 
                            "Glutamatergic 2", 
                            "Microglia", 
                            "Oligodendrocyte", 
                            "Oligo Precursor")
EpiMap_brain_related_cells = c("brain",
                              "neurosph",
                              "pns")
RoadMap_brain_related_cells = c("E007",
                                "E009",
                                "E010",
                                "E012",
                                "E053",
                                "E054",
                                "E067",
                                "E068",
                                "E069",
                                "E070",
                                "E071",
                                "E072",
                                "E073",
                                "E074",
                                "E081",
                                "E082")
PCHiC_brain_related_cells = c("FC",
                              "HCmerge",
                              "NPC")


ABC=ABC[which(ABC$V11 %in% ABC_brain_related_cells),]
EpiMap=EpiMap[which(EpiMap$V11 %in% EpiMap_brain_related_cells),]
RoadMap=RoadMap[which(RoadMap$V12 %in% RoadMap_brain_related_cells),]
PCHiC=PCHiC[which(PCHiC$V10 %in% PCHiC_brain_related_cells),]



## ********************************************************************************************
# Brain tissue
## ********************************************************************************************

index=match(TSS$V6,result$gene_name,nomatch=0)
result$ClosestTSS[index]=1
index=match(Exon$V9,result$gene_id,nomatch=0)
result$Exon[index]=1
feature=result

if(nrow(ABC)>0){
maxABC=data.frame(ABC %>% group_by(V4) %>% filter(V10==max(V10)))
index=match(maxABC$V9,result$gene_name,nomatch=0)
result$ABC[index]=1

maxABC_score=data.frame(maxABC %>% group_by(V9) %>%  filter(V10 == max(V10)) %>% distinct(V9, .keep_all = TRUE))
index=match(feature$gene_name, maxABC_score$V9,nomatch=0)
feature$ABC[which(index!=0)]=maxABC_score$V10[index]

rmdupABC=as.data.frame(ABC %>% group_by(V4,V9) %>% filter(V10==max(V10)))
rmdupABC$gene_tss=NA
index=match(rmdupABC$V9, gencode1$gene_name, nomatch=0)
rmdupABC$gene_tss[which(index!=0)]=gencode1$gene_tss[index]
rmdupABC=rmdupABC[which(index!=0),]
rmdupABC$method="ABC"

rmdupABC = rmdupABC %>% group_by(V4) %>% mutate(highlight_index = if_else(V10 == max(V10), 1, 0)) %>% ungroup() %>% data.frame()
}

if(nrow(EpiMap)>0){
maxEpiMap=data.frame(EpiMap %>% group_by(V4) %>% filter(V10==max(V10)))
index=match(maxEpiMap$V9,result$gene_id,nomatch=0)
result$EpiMap[index]=1

maxEpiMap_score=data.frame(maxEpiMap %>% group_by(V9) %>%  filter(V10 == max(V10)) %>% distinct(V9, .keep_all = TRUE))
index=match(feature$gene_id, maxEpiMap_score$V9,nomatch=0)
feature$EpiMap[which(index!=0)]=maxEpiMap_score$V10[index]

rmdupEpiMap=as.data.frame(EpiMap %>% group_by(V4, V9) %>% filter(V10==max(V10)))
rmdupEpiMap$gene_tss=NA
index=match(rmdupEpiMap$V9, gencode1$gene_id, nomatch=0)
rmdupEpiMap$V9[which(index!=0)]=gencode1$gene_name[index]
rmdupEpiMap$gene_tss[which(index!=0)]=gencode1$gene_tss[index]
rmdupEpiMap=rmdupEpiMap[which(index!=0),]
rmdupEpiMap$method="EpiMap"

rmdupEpiMap = rmdupEpiMap %>% group_by(V4) %>% mutate(highlight_index = if_else(V10 == max(V10), 1, 0)) %>% ungroup() %>% data.frame()
}

if(nrow(RoadMap)>0){
maxRoadMap=data.frame(RoadMap %>% group_by(V4) %>% filter(V10==max(V10)))
index=match(maxRoadMap$V9,result$gene_id,nomatch=0)
result$RoadMap[index]=1

maxRoadMap_score=data.frame(maxRoadMap %>% group_by(V9) %>%  filter(V10 == max(V10)) %>% distinct(V9, .keep_all = TRUE))
index=match(feature$gene_id, maxRoadMap_score$V9,nomatch=0)
feature$RoadMap[which(index!=0)]=maxRoadMap_score$V10[index]

rmdupRoadMap=as.data.frame(RoadMap %>% group_by(V4, V9) %>% filter(V10==max(V10)))
rmdupRoadMap$gene_tss=NA
index=match(rmdupRoadMap$V9, gencode1$gene_id, nomatch=0)
rmdupRoadMap$gene_tss[which(index!=0)]=gencode1$gene_tss[index]
rmdupRoadMap$V9[which(index!=0)]=gencode1$gene_name[index]
rmdupRoadMap$V11=rmdupRoadMap$V12
rmdupRoadMap=rmdupRoadMap[which(index!=0),] %>% select(-V12)
rmdupRoadMap$method="RoadMap"

rmdupRoadMap = rmdupRoadMap %>% group_by(V4) %>% mutate(highlight_index = if_else(V10 == max(V10), 1, 0)) %>% ungroup() %>% data.frame()
}

if(nrow(PCHiC)>0){
maxPCHiC=data.frame(PCHiC %>% group_by(V4) %>% filter(V11==max(V11)))
index=match(maxPCHiC$V9,result$gene_name,nomatch=0)
result$PCHiC[index]=1

maxPCHiC_score=data.frame(maxPCHiC %>% group_by(V9) %>%  filter(V11 == max(V11)) %>% distinct(V9, .keep_all = TRUE))
index=match(feature$gene_name, maxPCHiC_score$V9,nomatch=0)
feature$PCHiC[which(index!=0)]=maxPCHiC_score$V11[index]

rmdupPCHiC=as.data.frame(PCHiC %>% group_by(V4, V9) %>% filter(V11==max(V11)))
rmdupPCHiC$gene_tss=NA
index=match(rmdupPCHiC$V9, gencode1$gene_name, nomatch=0)
rmdupPCHiC$gene_tss[which(index!=0)]=gencode1$gene_tss[index]
rmdupPCHiC=rmdupPCHiC[which(index!=0),]
rmdupPCHiC_tmp=rmdupPCHiC$V10
rmdupPCHiC$V10=rmdupPCHiC$V11
rmdupPCHiC$V11=rmdupPCHiC_tmp
rmdupPCHiC$method="PCHiC"

rmdupPCHiC = rmdupPCHiC %>% group_by(V4) %>% mutate(highlight_index = if_else(V10 == max(V10), 1, 0)) %>% ungroup() %>% data.frame()
}

result$GAMMA_V2G=apply(result[,c("ClosestTSS","Exon","ABC", "EpiMap", "RoadMap","PCHiC")],1,sum)
feature$GAMMA_V2G=result$GAMMA_V2G

result=result[,c("gene_id","gene_name", "chr","start","end","strand",
		"GWAS_LOCUS", "Lead_SNP", "Lead_SNP_BP",
		"GAMMA_V2G", "ClosestTSS","Exon","ABC", "EpiMap", "RoadMap","PCHiC")]

feature=feature[,c("gene_id","gene_name", "chr","start","end","strand",
		"GWAS_LOCUS", "Lead_SNP", "Lead_SNP_BP",
		"GAMMA_V2G", "ClosestTSS","Exon","ABC", "EpiMap", "RoadMap","PCHiC")]

plot=rbind(rmdupABC, rmdupEpiMap, rmdupRoadMap, rmdupPCHiC)

write.table(result,paste0(OUTPUT,"/V2G/score/",trait_name,"_GAMMA_V2G_brain.summary"),row=F,col=T,quo=F,sep="\t")
write.table(feature,paste0(OUTPUT,"/V2G/feature/",trait_name,"_GAMMA_V2G_brain.feature"),row=F,col=T,quo=F,sep="\t")
write.table(plot,paste0(OUTPUT,"/V2G/plot/",trait_name,"_GAMMA_V2G_brain.plot"),row=F,col=T,quo=F,sep="\t")



